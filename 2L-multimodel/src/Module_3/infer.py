from ANN.infer_graph_2L_fc import ANN_add_vars_constraints_to_MILP
from LR.infer_2LMM_LLR import LR_add_vars_constraints_to_MILP
from twolayered_MILP_2LMM_L import print_sdf_file, print_gstar_file
from create_base_MILP import create_base_MILP

import sys, time
import yaml
import pulp

####################################################
# IMPORTANT:
# Please specify the path of the cplex solver here
CPLEX_PATH= \
"/Applications/CPLEX_Studio2211/cplex/bin/x86-64_osx/cplex"
# CPLEX_PATH= \
# "/opt/cplex_12.10/cplex/bin/x86-64_linux/cplex"

CPLEX_MSG = False
CPLEX_TIMELIMIT = 0
solver_type = 1 # change 1 to 2 if use CBC solver
std_eps = 1e-5
fv_gen_name = "./2LMM_v018/FV_2LMM_V018"  # file of fv generator used to generate fv to check
####################################################

def main(argv):
    start = time.time()

    ########## preparation ##########
    # Read config file
    try:
        with open("./config.yaml") as file:
            config = yaml.safe_load(file)
    except Exception as e:
        print("Error: reading the config file")
        print(e)
        sys.exit(1)

    if config["output_prefix"] == None or \
        config["instance_file"] == None or \
        config["fringe_tree_file"] == None or \
        config["normalized_dataset_filename"] == None or \
        config["input_data"] == None:
        print("Error: Please specify (output_prefix), (instance_file), (fringe_tree_file), (normalized_dataset_filename), (input_data) in config.yaml")
        sys.exit(1)

    output_prefix = config["output_prefix"]
    instance_file = config["instance_file"]
    fringe_tree_file = config["fringe_tree_file"]
    normalized_dataset_filename = config["normalized_dataset_filename"]

    print("output_prefix:", output_prefix)
    print("instance_file:", instance_file)
    print("fringe_tree_file:", fringe_tree_file)
    print("normalized_dataset_filename:", normalized_dataset_filename)
    
    ## create MILP
    MILP = pulp.LpProblem("MultiModel")
    base_var = ()

    base_var = create_base_MILP(MILP, instance_file, fringe_tree_file, normalized_dataset_filename)

    MILP_results = [None] * len(config["input_data"])
    for (index, item) in enumerate(config["input_data"]):
        match item["model"]:
            case "ANN":
                print("\tmodel: ANN")
                y = ANN_add_vars_constraints_to_MILP(config, index, MILP, base_var)
                MILP_results[index] = (y)
            case "LR":
                print("\tmodel: LR")
                y, y_min, y_max = LR_add_vars_constraints_to_MILP(config, index, MILP, base_var)
                MILP_results[index] = (y, y_min, y_max)
            case _:
                print("\tmodel {} not supported".format(item["model"]))
                sys.exit(1)
    
    ########## solve and output ##########
    # Output all MILP variables and constraints
    MILP.writeLP(output_prefix + ".lp")
    
    init_end = time.time()
    print("Initializing Time:", "{:.3f}".format(init_end - start))
    
    num_vars = len(MILP.variables())
    num_ints = len([var for var in MILP.variables() if var.cat == pulp.LpInteger])
    bins = [v.name for v in MILP.variables()
                                if (v.cat == pulp.LpBinary or
                                (v.cat == "Integer" and
                                v.upBound and v.lowBound != None and
                                round(v.upBound - v.lowBound) == 1))]
    
    reals = [v.name for v in MILP.variables() if v.cat == "Continuous"]
    
    num_constraints = len(
                [c for c in MILP.constraints.items()])
    
    print("Number of variables:", num_vars)
    print(" - Integer :", num_ints)
    print(" - Binary  :", len(bins))
    print("Number of constraints:", num_constraints, end = "\n\n", flush = True)

    # Solve MILP
    if solver_type == 1:
        if CPLEX_TIMELIMIT > 0:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG,
                               timeLimit = CPLEX_TIMELIMIT)
        else:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG)
        # print("Start Solving Using CPLEX...")
        MILP.solve(CPLEX)
        solve_end = time.time()
    else:
        # print("Start Solving Using Coin-OR...")
        MILP.solve()
        solve_end = time.time()
    # layer_num = len(weights)
    
    # Output result of solving MILP to command-line
    # print("# -------- solve status --------")
    if pulp.LpStatus[MILP.status] == "Optimal":
        output_status = "Feasible"
    else:
        output_status = pulp.LpStatus[MILP.status]
    print("Status:", output_status)
    
    # y_star = (tree_vars['y'].value())
    if pulp.LpStatus[MILP.status] == "Optimal":
        for (index, item) in enumerate(config["input_data"]):
            match item["model"]:
                case "ANN":
                    y = MILP_results[index]
                    y_star = y.value()
                    print("index", index, " y*:", "{:.3f}".format(y_star))
                case "LR":
                    y, y_min, y_max = MILP_results[index]
                    y_star = y.value()
                    y_star_scaled = y_star * (y_max - y_min) + y_min
                    print("index", index, " y*:", "{:.3f}".format(y_star_scaled))
    
    print("Solving Time:", "{:.3f}".format(solve_end - init_end))

    # ############################################
    # # The following block of code is used to print out the value of feature vector #
    # # and the value of all variables used in MILP in to a file "test.txt" #
    # ############################################
    # #
    with open(output_prefix + "_test_all.txt", "w") as f:
        # f.write("######### Feature Vector ############\n")
        # for i in range(1, num_fv):
        #     f.write(stringoutput[i] + str(y[(1, i)].value()) + "\n")
        # f.write("\n\n\n")
        # f.write("#####################################\n")
    
        for var in MILP.variables():
            if type(var) is dict:
                for x in var:
                    if var[x].value() is None:
                        pass
                    elif abs(var[x].value()) > 0.0000001:
                        f.write(f"{var[x].name}: {var[x].value()}\n")
            elif var.value() is None:
                pass
            elif abs(var.value()) > 0.0000001:
                f.write(f"{var.name}: {var.value()}\n")
    
        f.write("\n")

    # #############################################
    
    if pulp.LpStatus[MILP.status] == "Optimal":

        strF, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, \
        n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G, mass_n, \
        t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, \
        head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F, \
        v_T, v_F, alpha_C, alpha_T, alpha_F, \
        beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, \
        e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        t_T, E_C, ch_LB, ch_UB, \
        set_F_v, set_F_E = base_var
    
        # Output SDF file
        outputfilename = output_prefix + ".sdf"
    
        index_C, index_T, graph_adj, graph_ind = print_sdf_file(
            t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, Lambda_int, Lambda_ex,
            head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F,
            n_G, v_T, v_F, alpha_C, alpha_T, alpha_F,
            beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F,
            e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F,
            outputfilename
        )

        # Output file of partition which will be used in graph generation
        outputfilename = output_prefix + "_partition.txt"
    
        print_gstar_file(
            graph_ind, chi_T,
            t_C, t_T, index_C, index_T, graph_adj, E_C, ch_LB, ch_UB,
            I_ge_one, I_ge_two, I_zero_one, I_equal_one,
            set_F_v, set_F_E,
            outputfilename
        )

        # ## Check the calculated descriptors
        # outputfilename = output_prefix + ".sdf"
        # test_prefix = output_prefix + "_test_tmp"
        # subprocess.run([fv_gen_name, f"{prop}.sdf", f"{prop}_test", outputfilename, test_prefix],
        #                stdout=subprocess.DEVNULL)
        # # os.system(f"{fv_gen_name} {prop}.sdf {prop}_test {outputfilename} ttt")
        # y_predict = lr_inverter.inspection(f"{test_prefix}_desc_norm.csv", normalized_dataset_filename, LR_filename, x_hat, std_eps)

        # print("Inspection value:  ",  y_predict[0])
        # print("Inspection value (scaled):  ", y_predict[0] * (y_max - y_min) + y_min)

if __name__ == '__main__':
    main(sys.argv)
