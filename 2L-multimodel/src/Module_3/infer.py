"""
multi-model 
"""
import sys
import time
import subprocess
import yaml
import pulp

from ANN.infer_graph_2L_fc import ANN_add_vars_constraints_to_MILP
from LR.infer_2LMM_LLR import LR_add_vars_constraints_to_MILP
from twolayered_MILP_2LMM_L import print_sdf_file, print_gstar_file
from create_base_MILP import create_base_MILP


import ANN.ann_inverter as ann_inverter
import LR.lr_inverter as lr_inverter

####################################################
# IMPORTANT:
# Please specify the path of the cplex solver here
CPLEX_PATH = "/Applications/CPLEX_Studio2211/cplex/bin/x86-64_osx/cplex"
# CPLEX_PATH= \
# "/opt/cplex_12.10/cplex/bin/x86-64_linux/cplex"

CPLEX_MSG = False
CPLEX_TIMELIMIT = 0
SOLVER_CASE = 1  # change 1 to 2 if use CBC solver
STD_EPS = 1e-5
# file of fv generator used to generate fv to check
FV_GEN_NAME = "sample_instance/2LMM_v019/FV_2LMM_V019"
####################################################


def print_help():
    """Print help message."""
    print("Usage: python3 infer.py [config.yaml]")
    print("Example: python3 infer.py config.yaml")


def config_get(config, key):
    """Get value from config file."""
    value = config.get(key)
    if value is None:
        print(f"Error: Please specify ({key}) in config.yaml")
        sys.exit(1)
    print(f"{key}:", value)
    return value


def main(argv):
    """Main function."""
    start = time.time()

    ########## preparation ##########
    # Read config file
    if len(argv) != 2:
        print_help()
        sys.exit(1)

    config_filename = argv[1]
    try:
        with open(f"{config_filename}", "r", encoding="utf-8") as file:
            config = yaml.safe_load(file)
    except Exception as ex:
        print(f"Error: {ex}")
        sys.exit(1)

    output_prefix = config_get(config, "output_prefix")
    instance_file = config_get(config, "instance_file")
    fringe_tree_file = config_get(config, "fringe_tree_file")
    original_dataset_filename = config_get(config, "original_dataset_filename")

    ## Create MILP
    MILP = pulp.LpProblem("MultiModel")
    base_var = ()

    base_var = create_base_MILP(
        MILP, instance_file, fringe_tree_file, original_dataset_filename
    )

    milp_results = [()] * len(config["input_data"])
    for index, item in enumerate(config["input_data"]):
        match item["model"]:
            case "ANN":
                milp_results[index] = ANN_add_vars_constraints_to_MILP(
                    config, index, MILP, base_var
                )
            case "LR":
                milp_results[index] = LR_add_vars_constraints_to_MILP(
                    config, index, MILP, base_var
                )
            case _:
                print(f"\tmodel {item['model']} not supported")
                sys.exit(1)

    ########## solve and output ##########
    # Output all MILP variables and constraints
    MILP.writeLP(output_prefix + ".lp")

    init_end = time.time()
    print("Initializing Time:", f"{init_end - start:.3f}")

    num_vars = len(MILP.variables())
    num_ints = len([var for var in MILP.variables() if var.cat == pulp.LpInteger])
    bins = [v.name for v in MILP.variables()
                                if (v.cat == pulp.LpBinary or
                                (v.cat == "Integer" and
                                v.upBound and v.lowBound is not None and
                                round(v.upBound - v.lowBound) == 1))]

    num_constraints = len(MILP.constraints.items())

    print("Number of variables:", num_vars)
    print(" - Integer :", num_ints)
    print(" - Binary  :", len(bins))
    print("Number of constraints:", num_constraints, end="\n\n", flush=True)

    # Solve MILP
    if SOLVER_CASE == 1:
        if CPLEX_TIMELIMIT > 0:
            cplex = pulp.CPLEX(
                path=CPLEX_PATH, msg=CPLEX_MSG, timeLimit=CPLEX_TIMELIMIT
            )
        else:
            cplex = pulp.CPLEX(path=CPLEX_PATH, msg=CPLEX_MSG)
        MILP.solve(cplex)
        solve_end = time.time()
    else:
        MILP.solve()
        solve_end = time.time()

    # Output result of solving MILP to command-line
    if pulp.LpStatus[MILP.status] == "Optimal":
        output_status = "Feasible"
    else:
        output_status = pulp.LpStatus[MILP.status]
    print("Status:", output_status)

    if pulp.LpStatus[MILP.status] == "Optimal":
        for index, item in enumerate(config["input_data"]):
            match item["model"]:
                case "ANN":
                    y, y_min, y_max, _, _, _ = milp_results[index]
                    y_star = y.value()
                    y_star_scaled = y_star * (y_max - y_min) + y_min
                    print("index", index, " y*:", f"{y_star_scaled:.3f}")
                case "LR":
                    y, y_min, y_max, _, _, _ = milp_results[index]
                    y_star = y.value()
                    y_star_scaled = y_star * (y_max - y_min) + y_min
                    print("index", index, " y*:", f"{y_star_scaled:.3f}")

    print("Solving Time:", f"{solve_end - init_end:.3f}")

    ############################################
    # The following block of code is used to print out the value of feature vector #
    # and the value of all variables used in MILP in to a file "test.txt" #
    ############################################
    with open(output_prefix + "_test_all.txt", "w") as file:
        for var in MILP.variables():
            if isinstance(var, dict):
                for x in var:
                    if var[x].value() is None:
                        pass
                    elif abs(var[x].value()) > 0.0000001:
                        file.write(f"{var[x].name}: {var[x].value()}\n")
            elif var.value() is None:
                pass
            elif abs(var.value()) > 0.0000001:
                file.write(f"{var.name}: {var.value()}\n")

        file.write("\n")

    ##############################################

    if pulp.LpStatus[MILP.status] == "Optimal":
        strF, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, \
        n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G, mass_n, \
        t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, \
        head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F, \
        v_T, v_F, alpha_C, alpha_T, alpha_F, \
        beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, \
        e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        E_C, ch_LB, ch_UB, \
        set_F_v, set_F_E = base_var

        # Output SDF file
        outputfilename = output_prefix + ".sdf"
    
        index_C, index_T, graph_adj, graph_ind = print_sdf_file(
            t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, Lambda_int, Lambda_ex,
            head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F,
            n_G, v_T, v_F, alpha_C, alpha_T, alpha_F,
            beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F,
            e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F,
            outputfilename,
        )
        
        # Output file of partition
        outputfilename = output_prefix + "_partition.txt"

        print_gstar_file(
            graph_ind, chi_T,
            t_C, t_T, index_C, index_T, graph_adj, E_C, ch_LB, ch_UB,
            I_ge_one, I_ge_two, I_zero_one, I_equal_one,
            set_F_v, set_F_E,
            outputfilename,
        )

        ## Check the calculated descriptors
        outputfilename = output_prefix + ".sdf"
        test_prefix = output_prefix + "_test_tmp"
        for index, item in enumerate(config["input_data"]):
            prop = config["input_data"][index]["prefix"]
            subprocess.run([FV_GEN_NAME, f"{prop}.sdf", f"{prop}_test", outputfilename, test_prefix],
                stdout=subprocess.DEVNULL , check=False)
            match item["model"]:
                case "ANN":
                    y, y_min, y_max, x_hat, ann_training_data_norm_filename, ann = milp_results[index]
                    y_predict, _ = ann_inverter.inspection(f"{test_prefix}_desc_norm.csv", ann_training_data_norm_filename, ann, x_hat, STD_EPS)
                case "LR":
                    y, y_min, y_max, x_hat, lr_filename, normalized_dataset_filename = milp_results[index]
                    y_star = y.value()
                    y_star_scaled = y_star * (y_max - y_min) + y_min
                    y_predict = lr_inverter.inspection(f"{test_prefix}_desc_norm.csv", normalized_dataset_filename, lr_filename, x_hat, STD_EPS)
                    y_predict = y_predict[0]
            print("Inspection value:  ",  y_predict)
            print("Inspection value (scaled):  ", y_predict * (y_max - y_min) + y_min)

if __name__ == "__main__":
    main(sys.argv)
