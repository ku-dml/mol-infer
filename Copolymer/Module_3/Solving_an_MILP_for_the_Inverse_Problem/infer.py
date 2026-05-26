from src.config import ConfigLoader
from src_milp.seed_graph import create_seed_graph
from src_inverter.inverter import create_inverter
from src_milp.assignment import create_assignment, check_assignment
from src.utils import _map_optional_list

import os, subprocess, sys, argparse
# import pulp_modified as pulp

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("config", help="config filename")
	parser.add_argument("output_prefix", help='output_prefix')

	args = parser.parse_args()

	return args

def main(args):
	cfg = ConfigLoader(args.config)
	cfg.exp_start()
	output_prefix = args.output_prefix

	### set up necessary information for seed graphs
	sg_list = dict()
	for (sg, sgl) in cfg.graphs.items():
		new_sg = create_seed_graph(sgl.type, sg)
		new_sg.read_seed_graph(sgl.filename)
		new_sg.read_fringe_trees(sgl.fringe_filename)
		new_sg.set_up_seed_graph()
		sg_list[sg] = new_sg

	### build C2 parts of variables and constraints
	MILP = pulp.LpProblem(name="mol_infer")
	for sg in cfg.graphs:
		new_sg = sg_list[sg]
		new_sg.set_up_variables()
		MILP = new_sg.add_constraints_C2(MILP)
	 
	### build C1 parts of variables and constraints
	inv_list = dict()
	for (inv, invl) in cfg.learning_models.items():
		new_inv = create_inverter(invl.type, inv)
		new_inv.set_up_inverter(invl)
		inv_list[inv] = new_inv

	### connect C1 and C2
	asgm_list = dict()
	output_list = dict()
	for asgml in cfg.assignments_list:
		inv = _map_optional_list(asgml.learning_model, inv_list)
		sg = _map_optional_list(asgml.seed_graph, sg_list)
		asgm = create_assignment(sg=sg, inv=inv, asgml=asgml)
		check_assignment(asgm)
		MILP = asgm.connect_C1_C2(MILP)
		asgm_list[asgml.ind] = asgm

	### solve and output
	MILP.writeLP(f"{output_prefix}.lp")
	cfg.init_end()
	print(f"Initializing Time: {cfg.init_end_time - cfg.start_time:.3f}")

	num_vars = len(MILP.variables())
	num_vars_ints = len([var for var in MILP.variables() if var.cat == pulp.LpInteger])
	num_vars_bins = len([v.name for v in MILP.variables() if (v.cat == pulp.LpBinary 
		or (v.cat == "Integer" and v.upBound != None and v.lowBound != None and round(v.upBound - v.lowBound) == 1))])
	num_vars_reals = len([v.name for v in MILP.variables() if v.cat == "Continuous"])
	num_constraints = len([c for c in MILP.constraints.items()])
	
	print("Number of variables:", num_vars)
	print(" - Integer :", num_vars_ints)
	print(" - Binary  :", num_vars_bins)
	print("Number of constraints:", num_constraints)

	# Only CPLEX now
	CPLEX = cfg.prepare_solver()
	MILP.solve(CPLEX)
	cfg.exp_end()

	if pulp.LpStatus[MILP.status] == "Optimal":
		output_status = "Feasible"
	else:
		output_status = pulp.LpStatus[MILP.status]
	print("Status:", output_status)

	if pulp.LpStatus[MILP.status] == "Optimal":
		for asgml in cfg.assignments_list:
			# y_MILP = asgm_list[asgml.ind].output.value()
			# y_scale = y_MILP * (inv_list[asgml.learning_model].y_max - inv_list[asgml.learning_model].y_min) + inv_list[asgml.learning_model].y_min
			y_MILP, y_scale = asgm_list[asgml.ind].get_output()

			print(f"Assignment: {asgml.ind}, MILP output: {y_MILP:.3f}, scaled: {y_scale:.3f}")
	print(f"Time: {cfg.end_time - cfg.start_time:.3f} s")

	if pulp.LpStatus[MILP.status] == "Optimal":
	
		# Output SDF file
		for sg in cfg.graphs:
			new_sg = sg_list[sg]
			sdf_filename_prefix = f"{output_prefix}_{new_sg.name}"
			new_sg.print_sdf(sdf_filename_prefix)

		# Output file of partition which will be used in graph generation
		for sg in cfg.graphs:
			new_sg = sg_list[sg]
			partition_filename = f"{output_prefix}_{new_sg.name}_partition.txt"
			new_sg.print_partition_file(partition_filename)

		## Check the calculated descriptors
		if cfg._debug:
			for asgml in cfg.assignments_list:
				if isinstance(asgml.seed_graph, str):
					sg = sg_list[asgml.seed_graph]
					asgm = asgm_list[asgml.ind]
					sdf_filename_prefix = f"{output_prefix}_{sg.name}"
					test_prefix = f"{output_prefix}_ASMG{asgml.ind}_test_tmp"
					y_insp, y_insp_scale = asgm.inspection(sdf_filename_prefix, test_prefix)
					print(f"Assignment: {asgml.ind}. Inspection output: {y_insp:.3f}, scaled: {y_insp_scale:.3f}")
				else: # List[str]
					sgs = [sg_list[_sg] for _sg in asgml.seed_graph]
					asgm = asgm_list[asgml.ind]
					sdf_filename_prefix = f"{output_prefix}"
					# sdf_filename_prefix = [f"{output_prefix}_{sg.name}" for sg in sgs]
					test_prefix = f"{output_prefix}_ASMG{asgml.ind}_test_tmp"
					y_insp, y_insp_scale = asgm.inspection(sdf_filename_prefix, test_prefix)
					print(f"Assignment: {asgml.ind}. Inspection output: {y_insp:.3f}, scaled: {y_insp_scale:.3f}")
			
if __name__ == '__main__':
	main(get_args())



