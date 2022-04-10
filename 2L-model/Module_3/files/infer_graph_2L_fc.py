from MILP_2L_fc import *
from read_instance_2L_fc import *
import ann_inverter

import subprocess

####################################################
# IMPORTANT:
# Please specify the path of the cplex solver here
CPLEX_PATH= \
"/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/cplex"
# "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex"
# CPLEX_PATH= \
# "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"
# CPLEX_PATH= \
# "/opt/cplex_12.10/cplex/bin/x86-64_linux/cplex"
# CPLEX_PATH= \
# "/usr/local/cplex/bin/cplex"

# for bash script
if len(sys.argv) == 8:
    CPLEX_PATH = sys.argv[7]

CPLEX_MSG = False
CPLEX_TIMELIMIT = 3600
####################################################

def main(argv):

    if len(argv) < 7:
        sys.exit()

    start = time.time()

    set_Lambda = prepare_CG_element_info()

    prop = argv[1]
    target_value = float(argv[2])
    instance_file = argv[3]
    fringe_tree_file = argv[4]
    outputfileprefix = argv[5]
    solver_type = int(argv[6])

    # prop = "Kow3"
    # target_value = -100
    # instance_file = "instance_a.txt"
    # fringe_tree_file = "ins_a_fringe.txt"  # fringe trees that will be used in MILP
    # outputfileprefix = "Kow3"

#     try:
#         prop, target_value, instance_file, outputfileprefix =\
#             sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
#     except:
#         sys.stderr.write('''
# usage: {} (ANNfile_prefix) (target_value) (instance_file) (outputfile_prefix)
#
#   - For ANNfile, you need *_desc.csv, *_biases.txt and *_weights.txt
#
# '''.format(sys.argv[0]))
#         exit(1)

    #input file of bias
    ann_bias_filename = "{}_biases.txt".format(prop)
    #input file of weight
    ann_weights_filename = "{}_weights.txt".format(prop)  
     #input file of fv
    ann_training_data_filename = "{}_desc.csv".format(prop)
    ann_training_data_norm_filename = "{}_desc_norm.csv".format(prop)

    fv_fringe_tree_filename = "{}_fringe.txt".format(prop)    # all fringe trees used in learning

    training_data = ann_inverter.read_training_data(
        ann_training_data_norm_filename)
    des = ann_inverter.read_fv_descriptors(ann_training_data_norm_filename)
    des = des[1:] # Drop the "CID" descriptor
    weights, biases = ann_inverter.read_trained_ANN(ann_weights_filename,
                                                    ann_bias_filename)
    ann = ann_inverter.ANN(weights, biases)
    milp_ann_constants = ann_inverter.initialize_constants(
        ann, training_data)
    ann_a, ann_b, ann_b_hat, ann_c, ann_z = milp_ann_constants

    # solver_type = 1 # change 1 to 2 if use CBC solver
    
    V_C, E_C, \
    E_ge_two, E_ge_one, E_zero_one, E_equal_one, \
    I_ge_two, I_ge_one, I_zero_one, I_equal_one, \
    ell_LB, ell_UB, n_LB_int, n_UB_int, \
    n_LB, n_star, rho, \
    ch_LB, ch_UB, bl_LB, bl_UB, \
    Lambda, Lambda_dg_int, Gamma_int_ac, Gamma_int, \
    Lambda_star, na_LB, na_UB, Lambda_int, \
    na_LB_int, na_UB_int, ns_LB_int, ns_UB_int, \
    ac_LB_int, ac_UB_int, ec_LB_int, ec_UB_int, \
    bd2_LB, bd2_UB, bd3_LB, bd3_UB, dg_LB, dg_UB = read_seed_graph(instance_file)

    set_F, Lambda_ex, strF = prepare_fringe_trees(fringe_tree_file)

    set_Lambda = prepare_Lambda(set_Lambda, Lambda)

    F_Lambda = prepare_F_Lambda(Lambda)

    set_Lambda_dg, Lambda_dg, epsilon, epsilon_dg, Code_Lambda_dg, \
    Code_Lambda, Code_Lambda_int, Code_Lambda_ex, \
    MAX_CODE, MAX_CODE_dg, MAX_CODE_int, MAX_CODE_ex, \
    Code_Gamma_ec_int, Code_Gamma_ac_int, val, mass = \
        prepare_Lambda_dg(set_Lambda, Lambda, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_ac)

    Gamma_int_less, Gamma_int_equal, Gamma_int_great, \
    Gamma_int_ac_less, Gamma_int_ac_equal, Gamma_int_ac_great, \
    Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, \
    Gamma_tilde_ac_F, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, \
    Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, \
    Gamma_tilde_ec_F, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF = \
        prepare_Gamma_ac(Gamma_int, Gamma_int_ac, Code_Lambda, Code_Lambda_dg)

    Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F = \
        prepare_Lambda_tilde(Lambda_dg_int)

    t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, c_F, m_C, \
    head_C, tail_C, tail_F, E_C_plus, E_C_minus, \
    I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, \
    I_zero_one_plus, I_zero_one_minus, I_equal_one_plus, I_equal_one_minus, delta_i = \
        prepare_dataset_for_scheme_graph(n_star, rho, V_C, E_C, E_ge_two, E_ge_one,
        n_LB_int, n_UB_int, bl_LB, bl_UB, ch_LB, ch_UB, ell_UB, I_ge_two, I_ge_one,
        I_equal_one, I_zero_one, bd2_LB, bd3_LB, val, Lambda_star)

    Code_F, n_psi, deg_r, beta_r, atom_r, ht, F_C, F_T, F_F, \
    n_C, n_T, n_F, F_Cp, F_Tp, F_Fp, set_F_E, set_F_v, n_H, \
    na_alpha_ex, alpha_r = prepare_fringe_tree(set_F, F_Lambda, 
        V_C, t_T, t_F, val, Lambda_int, Lambda_ex, Code_Lambda)

    # print(Lambda_ex)
    # print(ht)
    # for psi in (set_F + F_Lambda):
    #     print(str(Code_F[psi]) + " " + str(n_psi[Code_F[psi]]) + " " + \
    #         str(ht[Code_F[psi]]) + " " + str(alpha_r[Code_F[psi]]) + " " + \
    #         str(deg_r[Code_F[psi]]) + " " + str(beta_r[Code_F[psi]]))

    # print(n_H)
    # print(na_alpha_ex)

    MILP = pulp.LpProblem(name="MILP_2layered_fc")

    # Prepare variables of Appendix A.1 Selecting Core-vertices and Core-edges
    e_C, v_T, e_T, chi_T, clr_T, delta_chi_T, deg_tilde_C_plus, \
    deg_tilde_C_minus = prepare_variables_selecting_core(
        t_C, k_C_tilde, k_C, t_T, m_C, ell_LB, ell_UB, bl_LB, bl_UB)

    # Add constraints of Appendix A.1 Selecting Core-vertices and Core-edges
    MILP = add_constraints_selecting_core(
        MILP,
        t_C, k_C_tilde, k_C, t_T, m_C, n_LB_int, n_UB_int, ell_LB, ell_UB, bl_LB, bl_UB,
        I_equal_one, I_equal_one_minus, I_equal_one_plus, I_ge_two, I_ge_one,
        I_ge_one_minus, I_ge_one_plus, I_zero_one_minus, I_zero_one_plus,
        e_C, v_T, e_T, delta_chi_T,
        chi_T, clr_T, deg_tilde_C_plus, deg_tilde_C_minus)
    
    # Prepare variables of Appendix A.2 Constraints for Including Internal Vertices and Edges
    bl_G, v_F, e_F, chi_F, clr_F, delta_chi_F, bl, n_G_int = prepare_variables_internal_vertices_and_edges(
        t_F, t_T, t_C, t_C_tilde, c_F, V_C, E_ge_one, E_ge_two, I_ge_one, I_ge_two, 
        bl_LB, bl_UB, E_C, n_LB_int, n_UB_int)

    # Add constraints of Appendix A.2 Constraints for Including Internal Vertices and Edges
    MILP = add_constraints_internal_vertices_and_edges(
        MILP,
        t_T, t_F, t_C, t_C_tilde, c_F, tail_F, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C,
        delta_chi_F, chi_T, e_F, v_F, v_T, chi_F, clr_F, bl, bl_G, n_G_int)

    # Prepare variables of Appendix A.3 Constraints for Including Fringe-trees
    n_G, h_T, h_C, h_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        deg_ex_C, deg_ex_T, deg_ex_F, sigma = \
        prepare_variables_fringe_trees(n_LB, n_star, rho, ch_LB, ch_UB, t_T, t_C, t_F, 
            n_T, n_C, n_F, delta_i, I_ge_two, I_ge_one, V_C, E_ge_one, E_ge_two, 
            v_T, v_F, F_C, F_T, F_F, Code_F, ht, F_Lambda)

    # Add constraints of Appendix A.3 Constraints for Including Fringe-trees
    MILP = add_constraints_fringe_trees(
        MILP,
        t_T, t_F, t_C, t_C_tilde, c_F, I_ge_one, I_ge_two, rho, n_star, 
        n_T, n_C, n_F, ch_LB, ch_UB, E_C, F_C, F_T, F_F, F_Cp, F_Tp, F_Fp,
        F_Lambda, Code_F, val, n_psi, deg_r, ht,
        delta_chi_F, delta_chi_T, e_F, v_F, v_T, h_T, h_C, h_F,
        sigma, delta_fr_C, delta_fr_T, delta_fr_F,
        chi_T, chi_F, clr_F, n_G, deg_ex_C, deg_ex_T, deg_ex_F)

    # Prepare variables of A.4 Descriptor for the Number of Specified Degree
    deg_C, deg_T, deg_F, deg_CT, deg_TC, \
        delta_dg_C, delta_dg_T, delta_dg_F, dg, \
        deg_int_C, deg_int_T, deg_int_F, \
        delta_int_dg_C, delta_int_dg_T, delta_int_dg_F, \
        dg_int = prepare_variables_degree(t_C, t_T, t_F, n_C, n_T, n_F, 
            delta_i, n_star, dg_LB, dg_UB, rho)

    # Add constraints of A.4 Descriptor for the Number of Specified Degree
    MILP = add_constraints_degree(
        MILP,
        t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, delta_i,
        I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, rho,
        F_C, F_T, F_F, F_Cp, F_Tp, F_Fp, Code_F, deg_r,
        delta_dg_C, delta_dg_T, delta_dg_F, delta_int_dg_C, delta_int_dg_T, delta_int_dg_F,
        e_T, e_F, v_F, v_T, delta_chi_T, delta_chi_F, delta_fr_C,
        deg_C, deg_T, deg_F, deg_CT, deg_TC, dg, deg_int_C, deg_int_T, deg_int_F,
        dg_int, deg_tilde_C_minus, deg_tilde_C_plus, deg_ex_C, deg_ex_T, deg_ex_F)

    # Prepare variables of A.5 Assigning Multiplicity
    beta_T, beta_F, beta_C, beta_plus, beta_minus, beta_in, beta_ex_C, beta_ex_T, beta_ex_F, \
        delta_beta_T, delta_beta_F, delta_beta_C, \
        delta_beta_plus, delta_beta_minus, delta_beta_in, bd_int, \
        bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF = \
            prepare_variables_multiplicity(t_C, t_T, t_F, n_C, n_T, n_F, c_F, delta_i,
                I_equal_one, I_zero_one, I_ge_one, I_ge_two, E_C,
                bd2_LB, bd2_UB, bd3_LB, bd3_UB, n_UB_int, n_star, rho)
            
    # Add constraints of A.5 Assigning Multiplicity
    MILP = add_constraints_multiplicity(
        MILP,
        t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, c_F, delta_i,
        head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two,
        E_C, bd2_LB, bd2_UB, bd3_LB, bd3_UB, rho, F_C, F_T, F_F,
        Code_F, beta_r, e_C, e_T, e_F, v_F, v_T,
        delta_chi_T, delta_chi_F, delta_beta_C, delta_beta_T, delta_beta_F,
        delta_beta_plus, delta_beta_minus, delta_beta_in,
        delta_fr_C, delta_fr_T, delta_fr_F,
        beta_C, beta_T, beta_F, beta_plus, beta_minus, beta_in,
        bd_int, bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF,
        beta_ex_C, beta_ex_T, beta_ex_F)

    # Prepare variables of A.6 Assigning Chemical Elements and Valence Condition
    beta_CT, beta_TC, beta_CF, beta_TF, alpha_C, alpha_T, alpha_F, MASS, \
        delta_alpha_C, delta_alpha_T, delta_alpha_F, na, na_int, na_C, na_T, na_F, \
        na_ex_C, na_ex_T, na_ex_F, delta_hyd_C, delta_hyd_T, delta_hyd_F, hydg, \
        na_ex = \
        prepare_variables_chemical_elements(t_C, t_T, t_F, n_C, n_T, n_F, delta_i,
            n_star, Lambda, epsilon, na_LB, na_UB, rho, na_LB_int, na_UB_int, 
            Lambda_int, Lambda_ex, MAX_CODE, MAX_CODE_int, MAX_CODE_ex,
            Code_Lambda_int, Code_Lambda_ex, F_C, F_T, F_F, Code_F, alpha_r, na_alpha_ex,
            n_H)

    # Add constraints of A.6 Assigning Chemical Elements and Valence Condition
    MILP = add_constraints_chemical_elements(
        MILP,
        t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, c_F, delta_i, E_C,
        I_equal_one, I_zero_one, I_ge_one, I_ge_one_plus, I_ge_one_minus,
        I_ge_two, I_ge_two_plus, I_ge_two_minus, val, mass, 
        Lambda, Code_Lambda, Lambda_star, epsilon, rho, na_LB_int, na_UB_int,
        Lambda_int, Lambda_ex, MAX_CODE, MAX_CODE_int, MAX_CODE_ex,
        Code_Lambda_int, Code_Lambda_ex, F_C, F_T, F_F, F_Lambda, Code_F, alpha_r, na_alpha_ex, n_H,
        v_T, v_F, e_T, e_F, delta_chi_T, delta_chi_F, chi_T, chi_F, delta_alpha_C, delta_alpha_T, delta_alpha_F,
        delta_fr_C, delta_fr_T, delta_fr_F, delta_hyd_C, delta_hyd_T, delta_hyd_F,
        deg_C, deg_T, deg_F, beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF,
        beta_plus, beta_minus, beta_in, alpha_C, alpha_T, alpha_F, MASS, na, na_int,
        na_C, na_T, na_F, na_ex_C, na_ex_T, na_ex_F, beta_ex_C, beta_ex_T, beta_ex_F, na_ex, bd_int, hydg)

    # Prepare variables of A.7 Constraints for Bounds on the Number of Bonds
    bd_T = prepare_variables_number_of_bounds(k_C, t_T, bd_T)

    # Add constraints of A.7 Constraints for Bounds on the Number of Bonds
    MILP = add_constraints_number_of_bounds(
        MILP,
        t_T, k_C, E_C, I_equal_one, I_zero_one, bd2_LB, bd2_UB, bd3_LB, bd3_UB,
        chi_T, delta_beta_C, delta_beta_T, delta_beta_plus, delta_beta_minus, bd_T
    )

    # Prepare variables of A.8 Descriptors for the Number of Adjacency-configuration
    ac_int, ac_C, ac_T, ac_F, ac_CT, ac_TC, ac_CF, ac_TF, \
        delta_ac_C, delta_ac_T, delta_ac_F, \
        alpha_CT, alpha_TC, alpha_CF, alpha_TF, \
        Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, \
        Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus, \
        Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus, Delta_ac_TF_minus, \
        delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF = \
            prepare_variables_adjacency_configuration(
                t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F, n_T, n_C, n_F,
                delta_i, rho, Gamma_int_ac, Gamma_tilde_ac_C, Gamma_tilde_ac_T,
                Gamma_tilde_ac_F, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, Gamma_tilde_ac_CF,
                Gamma_tilde_ac_TF, ac_LB_int, ac_UB_int, Lambda_int, Lambda_ex, MAX_CODE_int, MAX_CODE_ex)
            
    # Add constraints of A.8 Descriptors for the Number of Adjacency-configuration
    MILP = add_constraints_adjacency_configuration(
        MILP,
        t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, k_C, k_C_tilde, m_C, c_F, tail_C, head_C,
        Gamma_int_ac, Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F, Gamma_tilde_ac_CT,
        Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, Gamma_int_ac_less, Gamma_int_ac_equal, 
        ac_LB_int, ac_UB_int, Lambda_int, Lambda_ex,
        Code_Lambda_int, Code_Lambda_ex, MAX_CODE_int, MAX_CODE_ex, rho, delta_i,
        delta_alpha_C, delta_alpha_T, delta_alpha_F, delta_beta_C, delta_beta_T, delta_beta_F,
        delta_beta_plus, delta_beta_minus, delta_beta_in, delta_chi_T, delta_chi_F, chi_T, chi_F,
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,
        e_C, e_T, e_F, v_T, v_F, Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus,
        Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, Delta_ac_CT_plus, Delta_ac_CT_minus,
        Delta_ac_TC_plus, Delta_ac_TC_minus, Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus,
        Delta_ac_TF_minus, alpha_C, alpha_T, alpha_F, alpha_CT, alpha_TC, alpha_CF, alpha_TF, beta_C,
        beta_T, beta_F, beta_plus, beta_minus, beta_in, ac_C, ac_T, ac_F, ac_CT, ac_TC, ac_CF, ac_TF,
        ac_int)

    # Prepare variables of A.9 Descriptor for the Number of Chemical Symbols
    ns_int, delta_ns_C, delta_ns_T, delta_ns_F = \
        prepare_variables_chemical_symbols(
            t_C, t_T, t_F, n_C, n_T, n_F, delta_i,
            n_star, n_LB_int, n_UB_int, Lambda_dg_int, epsilon_dg
        )

    # # Add constraints of A.9 Descriptor for the Number of Chemical Symbols
    MILP = add_constraints_chemical_symbols(
        MILP,
        t_C, t_T, t_F, n_C, n_T, n_F,
        Lambda_dg_int,
        Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F,
        Code_Lambda, Code_Lambda_int, Code_Lambda_ex, epsilon_dg, delta_i, rho,
        delta_ns_C, delta_ns_T, delta_ns_F,
        alpha_C, alpha_T, alpha_F, deg_C, deg_T, deg_F, ns_int)

    # Prepare variables of A.10 Descriptor for the Number of Edge-configurations
    ec_int, ec_C, ec_T, ec_F, ec_CT, ec_TC, ec_CF, ec_TF, \
        delta_ec_C, delta_ec_T, delta_ec_F, \
        delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF, \
        deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF, \
        Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, \
        Delta_ec_F_plus, Delta_ec_F_minus, \
        Delta_ec_CT_plus, Delta_ec_CT_minus, Delta_ec_TC_plus, Delta_ec_TC_minus, \
        Delta_ec_CF_plus, Delta_ec_CF_minus, Delta_ec_TF_plus, Delta_ec_TF_minus = \
            prepare_variables_edge_configuration(
                t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F,
                n_T, n_C, n_F, delta_i, rho, Gamma_int, 
                Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_F, 
                Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF,
                ec_LB_int, ec_UB_int)
    
    # Add constraints of A.10 Descriptor for the Number of Edge-configurations
    MILP = add_constraints_edge_configuration(
        MILP,
        t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, k_C, k_C_tilde, m_C, c_F,
        tail_C, head_C, Gamma_int, Gamma_int_less, Gamma_int_equal,
        Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_F, Gamma_tilde_ec_CT,
        Gamma_tilde_ec_TC, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF,
        Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F, Gamma_tilde_ac_CT,
        Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF,
        rho, delta_i, Code_Gamma_ac_int, Code_Gamma_ec_int,
        delta_dg_C, delta_dg_T, delta_dg_F, delta_chi_T, delta_chi_F, chi_T, chi_F,
        delta_beta_C, delta_beta_T, delta_beta_F, delta_beta_plus, delta_beta_minus, delta_beta_in,
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,
        delta_ec_C, delta_ec_T, delta_ec_F, delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF,
        e_C, e_T, e_F, v_T, v_F, deg_C, deg_T, deg_F,
        deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF,
        ec_C, ec_T, ec_F, ec_CT, ec_TC, ec_CF, ec_TF,
        ec_int,Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, 
        Delta_ec_F_plus, Delta_ec_F_minus,Delta_ec_CT_plus, Delta_ec_CT_minus, 
        Delta_ec_TC_plus, Delta_ec_TC_minus, Delta_ec_CF_plus, Delta_ec_CF_minus, 
        Delta_ec_TF_plus, Delta_ec_TF_minus)

    # Prepare variables of A.11 Descriptor for the Number of Fringe-configurations
    fc = prepare_variables_fringe_configuration(t_C, t_T, t_F, set_F, Code_F)

    # Add constraints of A.11 Descriptor for the Number of Edge-configurations
    MILP = add_constraints_fringe_configuration(
        MILP,
        t_C, t_T, t_F, set_F, Code_F, delta_fr_C, delta_fr_T, delta_fr_F, fc)

    fv_fringe_tree, index_set = read_fringe_tree(fv_fringe_tree_filename, strF)
    
    descriptors, num_fv, mass_ind, max_dcp, min_dcp, avg_dcp, sd_dcp, forbidden_node = prepare_fv(
        ann_training_data_filename, Lambda_int, Lambda_ex, Gamma_int, fv_fringe_tree,
        index_set, n_G, n_G_int, MASS, dg, dg_int, hydg, bd_int, na_int, na_ex, ec_int, fc
    )

    milp_ann_vars = ann_inverter.initialize_lp_variables(
        ann, ann_a, ann_b, forbidden_node)
    # ann_descriptor_variables = ann_inverter.get_input_layer_variables(
    #     ann, milp_ann_vars, des, forbidden_node)

    eps = 0.02
    MILP = ann_inverter.build_MILP_ReLU(
        MILP,
        ann,
        milp_ann_vars,
        milp_ann_constants,
        target_value,
        eps,
        forbidden_node
    )

    _, y, _ = milp_ann_vars

    MILP, mass_n = add_constraints_mass_n(MILP, n_LB, n_star, n_G, MASS)

    x_hat, x_tilde = prepare_variables_nor_std_fv(num_fv)
    MILP = add_constraints_nor_std_fv(
        MILP,
        num_fv, mass_ind, descriptors, mass_n,
        max_dcp, min_dcp, avg_dcp, sd_dcp,
        x_hat, x_tilde,
        0.01,
        forbidden_node
    )

    MILP = add_constraints_ANN(MILP,
                        x_hat,
                        num_fv,
                        y, mass_ind,
                        mass_n,
                        forbidden_node)

    # Output all MILP variables and constraints
    MILP.writeLP("./MILP_2L_fc.lp")


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
    print("Number of constraints:", num_constraints)


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

    layer_num = len(weights)

    # Output result of solving MILP to command-line
    # print("# -------- solve status --------")
    if pulp.LpStatus[MILP.status] == "Optimal":
        output_status = "Feasible"
    else:
        output_status = pulp.LpStatus[MILP.status]
        
    y_star = (y[(layer_num + 1, 1)].value())
    print("Status:", output_status)
    print("Solving Time:", "{:.3f}".format(solve_end - init_end))
    print("MILP y*:", "{:.3f}".format(y_star))
  
    # ############################################
    # # The following block of code is used to print out the value of feature vector #
    # # and the value of all variables used in MILP in to a file "test.txt" #
    # ############################################
    # #
    with open(outputfileprefix + "_test_all.txt", "w") as f:
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
                    elif var[x].value() > 0:
                        f.write(f"{var[x].name}: {var[x].value()}\n")
            elif var.value() is None:
                pass
            elif var.value() > 0:
                f.write(f"{var.name}: {var.value()}\n")
    
        f.write("\n")
    #
    # #############################################

    if MILP.status == 1:

        # Output SDF file
        outputfilename = outputfileprefix + ".sdf"

        index_C, index_T, graph_adj = print_sdf_file(
            t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, Lambda_int, Lambda_ex,
            head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F,
            n_G, v_T, v_F, alpha_C, alpha_T, alpha_F,
            beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F,
            e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F,
            outputfilename
        )

        # Output file of partition which will be used in graph generation
        outputfilename = outputfileprefix + "_partition.txt"

        print_gstar_file(
            n_G, chi_T,
            t_C, t_T, index_C, index_T, graph_adj, E_C, ch_LB, ch_UB,
            I_ge_one, I_ge_two, I_zero_one, I_equal_one,
            set_F_v, set_F_E,
            outputfilename
        )
        
        # # ## Check the calculated descriptors
        # # for dd, var in ann_descriptor_variables.items():
        # #     print(dd, var.value())
        
        # ## Check the feature vector of the inferred graph
        # # call subprocess externally to calculate 
        # # the fv of the calculated graph
        # fv_result = subprocess.run(
        #     ["./fv", outputfileprefix + ".sdf"], 
        #         capture_output=True, text=True
        #     )  
        # fv_output = fv_result.stdout.split("\n")
        # fv_names = fv_output[0].split(",")[1:]
        # fv_values = [float(v) for v in fv_output[1].split(",")[1:]]
        
        # dv_dict = {d:v for d, v in zip(fv_names, fv_values)}
        
        # x_dagger = list()
        # x_star = list()
        # x_different = False
        # x_difference = dict()
            
        # for desc_name in des:
        #     if (desc_name in fv_names):
        #         x_dagger.append(dv_dict[desc_name])                
        #     else:
        #         x_dagger.append(0)
        #     x_star.append(ann_descriptor_variables[desc_name].value())
        #     if abs(x_dagger[-1] - x_star[-1]) > 0.001:
        #         x_different = True
        #         x_difference[desc_name] = (x_dagger[-1], x_star[-1])
                
        # y_dagger = ann.propagate(x_dagger)[0]
        # y_star = ann.propagate(x_star)[0]
        
        # print("ANN propagated y*:", "{:.3f}".format(y_star))
        # if x_different:
        #     print("Descriptor mismatch:")
        #     print("Descriptor: (x*, x)")
        #     for dname, dval in x_difference.items():
        #         print(f"{dname}:", dval)
        #     print("="*20)
        #     print("Inferred graph y:", "{:.3f}".format(y_dagger))
        #     print("Deviation:", 
        #           "{:.3f}".format(abs(target_value - y_dagger)*100/target_value), "%")
        # else:
        #     print("All descriptors match")
            
        
        # print("Descriptor,\tMILP,\tf(G)")
        # for dd in des:
        #     ys = ann_descriptor_variables[dd].value()
        #     if dd in fv_names:
        #         yd = dv_dict[dd]
        #     else:
        #         yd = 0
        #     print(f"{dd}\t{ys}\t{yd}")
        


if __name__ == '__main__':
    main(sys.argv)
