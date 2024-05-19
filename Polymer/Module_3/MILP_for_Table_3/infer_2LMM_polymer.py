from twolayered_MILP_2LMM_polymer import *
from read_instance_2layer_2LMM_polymer import *
import lr_inverter

import os, subprocess, sys


"""
MILP for 2LMM-LLR-polymer

Previous Varsion: milp_2LMM_polymer_0627_LR

add rank by IDO
"""

####################################################
# IMPORTANT:
# Please specify the path of the cplex solver here
CPLEX_PATH= \
"/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex"
# CPLEX_PATH= \
# "/opt/cplex_12.10/cplex/bin/x86-64_linux/cplex"
# CPLEX_PATH= \
# "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"

CPLEX_MSG = False
CPLEX_TIMELIMIT = 0
solver_type = 1 # change 1 to 2 if use CBC solver
std_eps = 1e-5
fv_gen_name = "./2LMM_polymer_v05.1/FV_2LMM_polymer_v051"  # file of fv generator used to generate fv to check
####################################################

def main(argv):

    # prop = "AmpD"
    # target_value_lb = 0
    # target_value_ub = 10
    # instance_file = "instance_a_polymer.txt"
    # fringe_tree_file = "ins_a_fringe_2LMM.txt"
    # output_prefix = "test_0717_a"
    ########## preparation ##########
    ## process arguments
    if len(argv) < 13:
        sys.stderr.write('usage: {} (prop_prefix_1)(prop_prefix_2)(prop_prefix_3)(target_value_lb_1)(target_value_ub_1)(target_value_lb_2)(target_value_ub_2)(target_value_lb_3)(target_value_ub_3)(instance_file)(FT_file)(output_prefix)\n\n'.format(argv[0]))
        sys.exit()
    prop1 = argv[1]
    prop2 = argv[2]
    prop3 = argv[3]
    target_value_lb1 = float(argv[4])
    target_value_ub1 = float(argv[5])
    target_value_lb2 = float(argv[6])
    target_value_ub2 = float(argv[7])
    target_value_lb3 = float(argv[8])
    target_value_ub3 = float(argv[9])
    instance_file = argv[10]
    fringe_tree_file = argv[11]
    output_prefix = argv[12]

    fq = None

    ### preprocess
    start = time.time()
    set_Lambda = prepare_CG_element_info()

    # file for linear regression
    LR_filename1 = "{}_linreg.txt".format(prop1)
    # file for original csv
    original_dataset_filename1 = "{}_desc.csv".format(prop1)
    # file for normalized csv
    normalized_dataset_filename1 = "{}_desc_norm.csv".format(prop1)
    # file for fringe trees
    fv_fringe_tree_filename1 = "{}_fringe.txt".format(prop1)  # all fringe trees used in learning
    # value file
    value_filename1 = "{}_values.txt".format(prop1)

    ### read dataset
    dataset1 = lr_inverter.read_training_data(normalized_dataset_filename1)
    des1 = lr_inverter.read_fv_descriptors(normalized_dataset_filename1)
    des1 = des1[1:]  # Drop the "CID" descriptor

     # file for linear regression
    LR_filename2 = "{}_linreg.txt".format(prop2)
    # file for original csv
    original_dataset_filename2 = "{}_desc.csv".format(prop2)
    # file for normalized csv
    normalized_dataset_filename2 = "{}_desc_norm.csv".format(prop2)
    # file for fringe trees
    fv_fringe_tree_filename2 = "{}_fringe.txt".format(prop2)  # all fringe trees used in learning
    # value file
    value_filename2 = "{}_values.txt".format(prop2)

    ### read dataset
    dataset2 = lr_inverter.read_training_data(normalized_dataset_filename2)
    des2 = lr_inverter.read_fv_descriptors(normalized_dataset_filename2)
    des2 = des2[1:]  # Drop the "CID" descriptor

     # file for linear regression
    LR_filename3 = "{}_linreg.txt".format(prop3)
    # file for original csv
    original_dataset_filename3 = "{}_desc.csv".format(prop3)
    # file for normalized csv
    normalized_dataset_filename3 = "{}_desc_norm.csv".format(prop3)
    # file for fringe trees
    fv_fringe_tree_filename3 = "{}_fringe.txt".format(prop3)  # all fringe trees used in learning
    # value file
    value_filename3 = "{}_values.txt".format(prop3)

    ### read dataset
    dataset3 = lr_inverter.read_training_data(normalized_dataset_filename3)
    des3 = lr_inverter.read_fv_descriptors(normalized_dataset_filename3)
    des3 = des3[1:]  # Drop the "CID" descriptor

    ### create MILP object
    MILP = pulp.LpProblem(name="MILP_2LMM_polymer")
    
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
    bd2_LB, bd2_UB, bd3_LB, bd3_UB, \
    dg_LB, dg_UB, \
    Gamma_ac_lnk, Gamma_lnk, ac_LB_lnk, ac_UB_lnk, \
    ec_LB_lnk, ec_UB_lnk, num_E_C_lnk, E_C_lnk, \
    n_LB_lnk, n_UB_lnk, \
    ac_LB_lf, ac_UB_lf, ac_LB_lf_common, ac_UB_lf_common, r_GC, \
    ns_LB_cnt, ns_UB_cnt = read_seed_graph(instance_file)

    Lambda, Lambda_dg_int, Lambda_int, Gamma_int, Gamma_int_ac, Gamma_lnk, Gamma_ac_lnk = \
        prepare_lambda_pre( 
            normalized_dataset_filename1,
            Lambda, Lambda_dg_int, Lambda_int, Gamma_int, Gamma_int_ac, Gamma_lnk, Gamma_ac_lnk
        )

    Lambda, Lambda_dg_int, Lambda_int, Gamma_int, Gamma_int_ac, Gamma_lnk, Gamma_ac_lnk = \
        prepare_lambda_pre( 
            normalized_dataset_filename2,
            Lambda, Lambda_dg_int, Lambda_int, Gamma_int, Gamma_int_ac, Gamma_lnk, Gamma_ac_lnk
        )
   
    Lambda, Lambda_dg_int, Lambda_int, Gamma_int, Gamma_int_ac, Gamma_lnk, Gamma_ac_lnk = \
        prepare_lambda_pre( 
            normalized_dataset_filename3,
            Lambda, Lambda_dg_int, Lambda_int, Gamma_int, Gamma_int_ac, Gamma_lnk, Gamma_ac_lnk
        )


    set_F, Lambda_ex, strF, fc_LB, fc_UB = prepare_fringe_trees(fringe_tree_file, Lambda)

    set_Lambda = prepare_Lambda(set_Lambda, Lambda)

    # F_Lambda = prepare_F_Lambda(Lambda)

    set_Lambda_dg, Lambda_dg, epsilon, epsilon_dg, Code_Lambda_dg, \
    Code_Lambda, Code_Lambda_int, Code_Lambda_ex, \
    MAX_CODE, MAX_CODE_dg, MAX_CODE_int, MAX_CODE_ex, \
    Code_Gamma_ec_int, Code_Gamma_ac_int, val, mass = \
        prepare_Lambda_dg(set_Lambda, Lambda, Lambda_int, Lambda_dg_int, Lambda_ex, Gamma_int, Gamma_int_ac)

    Gamma_int_less, Gamma_int_equal, Gamma_int_great, \
    Gamma_int_ac_less, Gamma_int_ac_equal, Gamma_int_ac_great, \
    Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, \
    Gamma_tilde_ac_F, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, \
    Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, \
    Gamma_tilde_ec_F, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF, \
    Gamma_lnk_less, Gamma_lnk_equal, Gamma_lnk_great, \
    Gamma_ac_lnk_less, Gamma_ac_lnk_equal, Gamma_ac_lnk_great, \
    Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great = \
        prepare_Gamma_ac(Gamma_int, Gamma_int_ac, Code_Lambda, Code_Lambda_dg, Gamma_lnk, Gamma_ac_lnk, ns_LB_cnt, ns_UB_cnt)

    Gamma_lf_ac, ac_LB_lf, ac_UB_lf = prepare_Gamma_lf_ac(set_Lambda,
        ac_LB_lf, ac_UB_lf, ac_LB_lf_common, ac_UB_lf_common
    )

    Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F = \
        prepare_Lambda_tilde(Lambda_dg_int)

    t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, c_F, m_C, \
    head_C, tail_C, tail_F, E_C_plus, E_C_minus, \
    I_lnk, \
    I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, \
    I_zero_one_plus, I_zero_one_minus, I_equal_one_plus, I_equal_one_minus, delta_i, n_lnk_equal_one = \
        prepare_dataset_for_scheme_graph(n_star, rho, V_C, E_C, E_C_lnk,  E_ge_two, E_ge_one, E_equal_one,
        n_LB_int, n_UB_int, bl_LB, bl_UB, ch_LB, ch_UB, ell_UB, I_ge_two, I_ge_one,
        I_equal_one, I_zero_one,bd2_LB, bd3_LB, val, Lambda_star)

    Code_F, n_psi_H, deg_r_H, deg_r_hyd, beta_r, atom_r, ht_H, F_C, F_T, F_F, \
    n_C, n_T, n_F, F_Cp, F_Tp, F_Fp, set_F_E, set_F_v, \
    na_alpha_ex, alpha_r, deg_fr, v_ion, ac_psi_lf = prepare_fringe_tree(set_F,
        V_C, t_T, t_F, val, Lambda_int, Lambda_ex, Code_Lambda_int, Gamma_lf_ac)

    # print(Lambda_ex)
    # print(ht_H)
    # for psi in set_F:
    #     print(str(Code_F[psi]) + " " + str(n_psi_H[Code_F[psi]]) + " " + \
    #         str(ht_H[Code_F[psi]]) + " " + str(alpha_r[Code_F[psi]]) + " " + \
    #         str(deg_r_H[Code_F[psi]]) + " " + str(deg_r_hyd[Code_F[psi]]) + " " + str(beta_r[Code_F[psi]]))
    #
    # # print(n_H)
    # print(na_alpha_ex)

    # MILP = pulp.LpProblem(name="MILP_2layered_fc")

    # Prepare variables of Appendix A.1 Selecting Core-vertices and Core-edges
    e_C, v_T, e_T, chi_T, clr_T, delta_chi_T, deg_tilde_C_plus, \
    deg_tilde_C_minus, n_lnk, rank_G = prepare_variables_selecting_core(
        t_C, k_C_tilde, k_C, t_T, m_C, n_LB_int, n_UB_int, n_LB_lnk, n_UB_lnk, ell_LB, ell_UB, bl_LB, bl_UB)
    
    # Add constraints of Appendix A.1 Selecting Core-vertices and Core-edges
    MILP = add_constraints_selecting_core(
        MILP,
        t_C, k_C_tilde, k_C, t_T, m_C, n_LB_int, n_UB_int, ell_LB, ell_UB, bl_LB, bl_UB,
        I_equal_one, I_equal_one_minus, I_equal_one_plus, I_ge_two, I_ge_one,
        I_ge_one_minus, I_ge_one_plus, I_zero_one, I_zero_one_minus, I_zero_one_plus, n_lnk_equal_one, I_lnk, r_GC,
        e_C, v_T, e_T, delta_chi_T,
        chi_T, clr_T, deg_tilde_C_plus, deg_tilde_C_minus, n_lnk, rank_G)
    
    # Prepare variables of Appendix A.2 Constraints for Including Internal Vertices and Edges
    n_G_int, v_F, e_F, chi_F, clr_F, delta_chi_F, bl = prepare_variables_internal_vertices_and_edges(
        t_F, t_T, t_C, t_C_tilde, c_F, n_LB_int, n_UB_int, V_C, E_ge_one, E_ge_two, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C)
    
    # Add constraints of Appendix A.2 Constraints for Including Internal Vertices and Edges
    MILP = add_constraints_internal_vertices_and_edges(
        MILP,
        t_T, t_F, t_C, t_C_tilde, c_F, tail_F, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C,
        delta_chi_F, chi_T, e_F, v_F, v_T, chi_F, clr_F, bl, n_G_int)
    
    # Prepare variables of Appendix A.3 Constraints for Including Fringe-trees
    n_G, h_T, h_C, h_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        deg_ex_C, deg_ex_T, deg_ex_F, hyddeg_C, hyddeg_T, hyddeg_F, eledeg_C, eledeg_T, eledeg_F, sigma, ac_lf = \
        prepare_variables_fringe_trees(n_LB, n_star, rho, ch_LB, ch_UB, t_T, t_C, t_F,
            n_T, n_C, n_F, delta_i, I_ge_two, I_ge_one, V_C, E_ge_one, E_ge_two,
            v_T, v_F, F_C, F_T, F_F, Code_F, Gamma_lf_ac, ac_LB_lf, ac_UB_lf)
    
    # Add constraints of Appendix A.3 Constraints for Including Fringe-trees
    MILP = add_constraints_fringe_trees(
        MILP,
        t_T, t_F, t_C, t_C_tilde, c_F, I_ge_one, I_ge_two, rho, n_star,
        n_T, n_C, n_F, ch_LB, ch_UB, E_C, F_C, F_T, F_F, F_Cp, F_Tp, F_Fp,
        Code_F, val, n_psi_H, deg_r_H, deg_r_hyd, ht_H, v_ion, ac_psi_lf, Gamma_lf_ac,
        delta_chi_F, delta_chi_T, e_F, v_F, v_T, h_T, h_C, h_F,
        sigma, delta_fr_C, delta_fr_T, delta_fr_F,
        chi_T, chi_F, clr_F, n_G, deg_ex_C, deg_ex_T, deg_ex_F, 
        hyddeg_C, hyddeg_T, hyddeg_F, eledeg_C, eledeg_T, eledeg_F, ac_lf)
    
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
        F_C, F_T, F_F, F_Cp, F_Tp, F_Fp, Code_F, deg_fr,
        delta_dg_C, delta_dg_T, delta_dg_F, delta_int_dg_C, delta_int_dg_T, delta_int_dg_F,
        e_T, e_F, v_F, v_T, delta_chi_T, delta_chi_F, delta_fr_C, delta_fr_T, delta_fr_F,
        deg_C, deg_T, deg_F, deg_CT, deg_TC, dg, deg_int_C, deg_int_T, deg_int_F,
        dg_int, deg_tilde_C_minus, deg_tilde_C_plus, deg_ex_C, deg_ex_T, deg_ex_F, hyddeg_C, hyddeg_T, hyddeg_F)
    
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
        na_ex_C, na_ex_T, na_ex_F, \
        na_ex = \
        prepare_variables_chemical_elements(t_C, t_T, t_F, n_C, n_T, n_F, delta_i,
            n_star, Lambda, epsilon, na_LB, na_UB, rho, na_LB_int, na_UB_int,
            Lambda_int, Lambda_ex, MAX_CODE, MAX_CODE_int, MAX_CODE_ex,
            Code_Lambda_int, Code_Lambda_ex, F_C, F_T, F_F, Code_F, alpha_r, na_alpha_ex
        )
    
    # Add constraints of A.6 Assigning Chemical Elements and Valence Condition
    MILP = add_constraints_chemical_elements(
        MILP,
        t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, c_F, delta_i, E_C,
        I_equal_one, I_zero_one, I_ge_one, I_ge_one_plus, I_ge_one_minus,
        I_ge_two, I_ge_two_plus, I_ge_two_minus, val, mass,
        Lambda, Code_Lambda, Lambda_star, epsilon, rho, na_LB_int, na_UB_int,
        Lambda_int, Lambda_ex, MAX_CODE, MAX_CODE_int, MAX_CODE_ex,
        Code_Lambda_int, Code_Lambda_ex, F_C, F_T, F_F, Code_F, alpha_r, na_alpha_ex,
        v_T, v_F, e_T, e_F, delta_chi_T, delta_chi_F, chi_T, chi_F, delta_alpha_C, delta_alpha_T, delta_alpha_F,
        delta_fr_C, delta_fr_T, delta_fr_F, 
        deg_C, deg_T, deg_F, beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF,
        beta_plus, beta_minus, beta_in, alpha_C, alpha_T, alpha_F, MASS, na, na_int,
        na_C, na_T, na_F, na_ex_C, na_ex_T, na_ex_F, beta_ex_C, beta_ex_T, beta_ex_F, na_ex, bd_int)
    
    # Prepare variables of A.7 Constraints for Bounds on the Number of Bonds
    bd_T = prepare_variables_number_of_bounds(k_C, t_T, bd_T)
    
    # Add constraints of A.7 Constraints for Bounds on the Number of Bonds
    MILP = add_constraints_number_of_bounds(
        MILP,
        t_T, k_C, E_C, I_equal_one, I_zero_one, bd2_LB, bd2_UB, bd3_LB, bd3_UB,
        chi_T, delta_beta_C, delta_beta_T, delta_beta_plus, delta_beta_minus, bd_T
    )
    
    # Prepare variables of A.8 Descriptors for the Number of Adjacency-configuration  # IDO
    ac_int, ac_C, ac_T, ac_F, ac_CT, ac_TC, ac_CF, ac_TF, \
        delta_ac_C, delta_ac_T, delta_ac_F, \
        alpha_CT, alpha_TC, alpha_CF, alpha_TF, \
        Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, \
        Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus, \
        Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus, Delta_ac_TF_minus, \
        delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,  \
        ac_lnk, ac_C_lnk, ac_T_lnk, ac_CT_lnk, ac_TC_lnk, delta_ac_T_lnk = \
            prepare_variables_adjacency_configuration(
                t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F, n_T, n_C, n_F,
                delta_i, rho, Gamma_int_ac, Gamma_tilde_ac_C, Gamma_tilde_ac_T,
                Gamma_tilde_ac_F, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, Gamma_tilde_ac_CF,
                Gamma_tilde_ac_TF, ac_LB_int, ac_UB_int, Lambda_int, Lambda_ex, MAX_CODE_int, MAX_CODE_ex, 
                Gamma_ac_lnk, ac_LB_lnk, ac_UB_lnk)
    
    # Add constraints of A.8 Descriptors for the Number of Adjacency-configuration IDO
    MILP = add_constraints_adjacency_configuration(
        MILP,
        t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, k_C, k_C_tilde, m_C, c_F, tail_C, head_C,
        Gamma_int_ac, Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F, Gamma_tilde_ac_CT,
        Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, Gamma_int_ac_less, Gamma_int_ac_equal,
        ac_LB_int, ac_UB_int, Lambda_int, Lambda_ex,
        Code_Lambda_int, Code_Lambda_ex, MAX_CODE_int, MAX_CODE_ex, rho, delta_i,
        Gamma_ac_lnk, Gamma_ac_lnk_less, Gamma_ac_lnk_equal,
        delta_alpha_C, delta_alpha_T, delta_alpha_F, delta_beta_C, delta_beta_T, delta_beta_F,
        delta_beta_plus, delta_beta_minus, delta_beta_in, delta_chi_T, delta_chi_F, chi_T, chi_F,
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,
        e_C, e_T, e_F, v_T, v_F, Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus,
        Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, Delta_ac_CT_plus, Delta_ac_CT_minus,
        Delta_ac_TC_plus, Delta_ac_TC_minus, Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus,
        Delta_ac_TF_minus,
        delta_ac_T_lnk,
        clr_T,
        alpha_C, alpha_T, alpha_F, alpha_CT, alpha_TC, alpha_CF, alpha_TF, beta_C,
        beta_T, beta_F, beta_plus, beta_minus, beta_in, ac_C, ac_T, ac_F, ac_CT, ac_TC, ac_CF, ac_TF,
        ac_int,
        ac_lnk, ac_C_lnk, ac_T_lnk, ac_CT_lnk, ac_TC_lnk, I_lnk, n_lnk)
    
    # Prepare variables of A.9 Descriptor for the Number of Chemical Symbols
    ns_int, delta_ns_C, delta_ns_T, delta_ns_F = \
        prepare_variables_chemical_symbols(
            t_C, t_T, t_F, n_C, n_T, n_F, delta_i,
            n_star, n_LB_int, n_UB_int, ns_LB_int, ns_UB_int, Lambda_dg_int, epsilon_dg
        )
    
    # Add constraints of A.9 Descriptor for the Number of Chemical Symbols
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
        Delta_ec_CF_plus, Delta_ec_CF_minus, Delta_ec_TF_plus, Delta_ec_TF_minus, \
        ec_lnk, ec_C_lnk, ec_T_lnk, ec_CT_lnk, ec_TC_lnk, delta_ec_T_lnk, delta_cnt = \
            prepare_variables_edge_configuration(
                t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F,
                n_T, n_C, n_F, delta_i, rho, Gamma_int,
                Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_F,
                Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF,
                ec_LB_int, ec_UB_int, ec_LB_lnk, ec_UB_lnk, Gamma_lnk, Lambda_dg_int,
                Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great)
    
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
        Gamma_lnk, Gamma_lnk_less, Gamma_lnk_equal, Lambda_dg_int, ns_LB_cnt, ns_UB_cnt, I_lnk, 
        Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great,
        delta_dg_C, delta_dg_T, delta_dg_F, delta_chi_T, delta_chi_F, chi_T, chi_F,
        delta_beta_C, delta_beta_T, delta_beta_F, delta_beta_plus, delta_beta_minus, delta_beta_in,
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,
        delta_ec_C, delta_ec_T, delta_ec_F, delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF,
        e_C, e_T, e_F, v_T, v_F,
        delta_ec_T_lnk, delta_cnt,
        clr_T, deg_C, deg_T, deg_F,
        deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF,
        ec_C, ec_T, ec_F, ec_CT, ec_TC, ec_CF, ec_TF,
        ec_int,Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus,
        Delta_ec_F_plus, Delta_ec_F_minus,Delta_ec_CT_plus, Delta_ec_CT_minus,
        Delta_ec_TC_plus, Delta_ec_TC_minus, Delta_ec_CF_plus, Delta_ec_CF_minus,
        Delta_ec_TF_plus, Delta_ec_TF_minus,
        ec_lnk, ec_T_lnk, ec_C_lnk, ec_CT_lnk, ec_TC_lnk, n_lnk)
    
    # Prepare variables of A.11 Descriptor for the Number of Fringe-configurations
    fc = prepare_variables_fringe_configuration(t_C, t_T, t_F, set_F, Code_F, fc_LB, fc_UB)
    
    # Add constraints of A.11 Descriptor for the Number of Edge-configurations
    MILP = add_constraints_fringe_configuration(
        MILP,
        t_C, t_T, t_F, set_F, Code_F, delta_fr_C, delta_fr_T, delta_fr_F, fc)
    
    fv_fringe_tree1, index_set1 = read_fringe_tree(fv_fringe_tree_filename1, strF)

    fv_fringe_tree2, index_set2 = read_fringe_tree(fv_fringe_tree_filename2, strF)

    fv_fringe_tree3, index_set3 = read_fringe_tree(fv_fringe_tree_filename3, strF)
    
    descriptors1, num_fv1, mass_ind1, max_dcp1, min_dcp1, avg_dcp1, sd_dcp1, forbidden_node1, I_integer1, I_nonneg1 = prepare_fv(
        original_dataset_filename1, Lambda_int, Lambda_ex, Lambda_dg_int, Gamma_int, Gamma_int_less, Gamma_int_equal,
        Gamma_lnk, Gamma_lnk_less, Gamma_lnk_equal, Gamma_lf_ac, fv_fringe_tree1,
        index_set1, n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, fq, ac_lf, n_lnk, ec_lnk, rank_G, 
        delta_cnt, Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great
    )

    descriptors2, num_fv2, mass_ind2, max_dcp2, min_dcp2, avg_dcp2, sd_dcp2, forbidden_node2, I_integer2, I_nonneg2 = prepare_fv(
        original_dataset_filename2, Lambda_int, Lambda_ex, Lambda_dg_int, Gamma_int, Gamma_int_less, Gamma_int_equal,
        Gamma_lnk, Gamma_lnk_less, Gamma_lnk_equal, Gamma_lf_ac, fv_fringe_tree2,
        index_set2, n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, fq, ac_lf, n_lnk, ec_lnk, rank_G, 
        delta_cnt, Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great
    )

    descriptors3, num_fv3, mass_ind3, max_dcp3, min_dcp3, avg_dcp3, sd_dcp3, forbidden_node3, I_integer3, I_nonneg3 = prepare_fv(
        original_dataset_filename3, Lambda_int, Lambda_ex, Lambda_dg_int, Gamma_int, Gamma_int_less, Gamma_int_equal,
        Gamma_lnk, Gamma_lnk_less, Gamma_lnk_equal, Gamma_lf_ac, fv_fringe_tree3,
        index_set3, n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, fq, ac_lf, n_lnk, ec_lnk, rank_G, 
        delta_cnt, Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great
    )

    # A12
    MILP, mass_n = add_constraints_mass_n(MILP, n_LB, n_star, n_G, na_LB, na_UB, na_ex, MASS)

    x_hat1, x_tilde1 = prepare_variables_nor_std_fv(num_fv1, 1)
    MILP = add_constraints_nor_std_fv(
        MILP,
        num_fv1, mass_ind1, descriptors1, mass_n,
        max_dcp1, min_dcp1, avg_dcp1, sd_dcp1,
        x_hat1, x_tilde1,
        std_eps,
        forbidden_node1,
        1
    )

    y_min1, y_max1 = get_value(value_filename1)

    x_hat2, x_tilde2 = prepare_variables_nor_std_fv(num_fv2, 2)
    MILP = add_constraints_nor_std_fv(
        MILP,
        num_fv2, mass_ind2, descriptors2, mass_n,
        max_dcp2, min_dcp2, avg_dcp2, sd_dcp2,
        x_hat2, x_tilde2,
        std_eps,
        forbidden_node2,
        2
    )

    y_min2, y_max2 = get_value(value_filename2)

    x_hat3, x_tilde3 = prepare_variables_nor_std_fv(num_fv3, 3)
    MILP = add_constraints_nor_std_fv(
        MILP,
        num_fv3, mass_ind3, descriptors3, mass_n,
        max_dcp3, min_dcp3, avg_dcp3, sd_dcp3,
        x_hat3, x_tilde3,
        std_eps,
        forbidden_node3,
        3
    )

    y_min3, y_max3 = get_value(value_filename3)

    ########## Inverse problem: LR part ##########
    ### read LR file to obtain constants
    fp1 = open(LR_filename1)
    LR1 = lr_inverter.read_LR(fp1)
        
    ### prepare variables (tree part)
    weight_var1, y1 = LR1.build_weight_var(I_integer1, I_nonneg1, 1)
    
    ### add constraints (tree part)
    target_value_lb1 = (target_value_lb1 - y_min1) / (y_max1 - y_min1)
    target_value_ub1 = (target_value_ub1 - y_min1) / (y_max1 - y_min1)
    
    LR1.build_constraints(MILP, target_value_lb1, target_value_ub1, 1)
    fp1.close()

    MILP = add_constraints__LR(MILP, x_hat1, num_fv1, LR1, mass_ind1, mass_n, forbidden_node1, "1")

    fp2 = open(LR_filename2)
    LR2 = lr_inverter.read_LR(fp2)
        
    ### prepare variables (tree part)
    weight_var2, y2 = LR2.build_weight_var(I_integer2, I_nonneg2, 2)
    
    ### add constraints (tree part)
    target_value_lb2 = (target_value_lb2 - y_min2) / (y_max2 - y_min2)
    target_value_ub2 = (target_value_ub2 - y_min2) / (y_max2 - y_min2)
    
    LR2.build_constraints(MILP, target_value_lb2, target_value_ub2, 2)
    fp2.close()

    MILP = add_constraints__LR(MILP, x_hat2, num_fv2, LR2, mass_ind2, mass_n, forbidden_node2, "2")

    fp3 = open(LR_filename3)
    LR3 = lr_inverter.read_LR(fp3)
        
    ### prepare variables (tree part)
    weight_var3, y3 = LR3.build_weight_var(I_integer3, I_nonneg3, 3)
    
    ### add constraints (tree part)
    target_value_lb3 = (target_value_lb3 - y_min3) / (y_max3 - y_min3)
    target_value_ub3 = (target_value_ub3 - y_min3) / (y_max3 - y_min3)
    
    LR3.build_constraints(MILP, target_value_lb3, target_value_ub3, 3)
    fp3.close()

    MILP = add_constraints__LR(MILP, x_hat3, num_fv3, LR3, mass_ind3, mass_n, forbidden_node3, "3")

    ########## solve and output ##########
    # Output all MILP variables and constraints
    # MILP.writeLP("./MILP_2LMH_DT_" + output_prefix + ".lp")
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
    print("Solving Time:", "{:.3f}".format(solve_end - init_end))
    
    # layer_num = len(weights)
    
    # Output result of solving MILP to command-line
    # print("# -------- solve status --------")
    if pulp.LpStatus[MILP.status] == "Optimal":
        output_status = "Feasible"
    else:
        output_status = pulp.LpStatus[MILP.status]
    
    y_star1 = y1.value()
    y_star_scaled1 = y_star1 * (y_max1 - y_min1) + y_min1
    y_star2 = y2.value()
    y_star_scaled2 = y_star2 * (y_max2 - y_min2) + y_min2
    y_star3 = y3.value()
    y_star_scaled3 = y_star3 * (y_max3 - y_min3) + y_min3
    print("Status:", output_status)
    print("Solving Time:", "{:.3f}".format(solve_end - init_end))
    # y_star = y.value()
    # y_star_scaled = y_star * (y_max - y_min) + y_min
    print("MILP y*1:", "{:.3f}".format(y_star1))
    print("MILP y*1 (scaled):", "{:.3f}".format(y_star_scaled1))
    print("MILP y*2:", "{:.3f}".format(y_star2))
    print("MILP y*2 (scaled):", "{:.3f}".format(y_star_scaled2))
    print("MILP y*3:", "{:.3f}".format(y_star3))
    print("MILP y*3 (scaled):", "{:.3f}".format(y_star_scaled3))

    print(f"n_G : {n_G.value()}, n_G_int : {n_G_int.value()}")
    
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
    #
    # #############################################
    
    if MILP.status == 1:
    
        # Output SDF file
        outputfilename = output_prefix + ".sdf"
    
        index_C, index_T, graph_adj, graph_ind = print_sdf_file(
            t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, Lambda_int, Lambda_ex,
            head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, I_lnk, set_F, Code_F, E_C,
            n_G, v_T, v_F, alpha_C, alpha_T, alpha_F,
            beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F,
            e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, delta_cnt, fq,
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

        ###############################
        ######## test ######
        # outputfilename = output_prefix + ".sdf"
        # test_prefix = output_prefix + "_test_tmp"
        # subprocess.run([fv_gen_name, f"{prop}.sdf", f"{prop}_test", outputfilename, test_prefix],
        #                stdout=subprocess.DEVNULL)
        # # os.system(f"{fv_gen_name} {prop}.sdf {prop}_test {outputfilename} ttt")
        # y_predict = ann_inverter.inspection(f"{test_prefix}_desc_norm.csv", ann_training_data_norm_filename, ann, x_hat, std_eps)

        # print("Inspection value:  ",  y_predict)
        # print("Inspection value (scaled):  ", y_predict * (y_max - y_min) + y_min)
        # outputfilename = output_prefix + ".sdf"
        # test_prefix = output_prefix + "_test_tmp1"
        # subprocess.run([fv_gen_name, f"{prop1}.sdf", f"{prop1}_test", outputfilename, test_prefix],
        #                stdout=subprocess.DEVNULL)
        # # os.system(f"{fv_gen_name} {prop}.sdf {prop}_test {outputfilename} ttt")
        # y_predict1 = lr_inverter.inspection(f"{test_prefix}_desc_norm.csv", normalized_dataset_filename1, LR_filename1, x_hat1, std_eps)

        # print("Inspection value:  ",  y_predict1[0])
        # print("Inspection value (scaled):  ", y_predict1[0] * (y_max1 - y_min1) + y_min1)

        # test_prefix = output_prefix + "_test_tmp2"
        # subprocess.run([fv_gen_name, f"{prop2}.sdf", f"{prop2}_test", outputfilename, test_prefix],
        #                stdout=subprocess.DEVNULL)
        # # os.system(f"{fv_gen_name} {prop}.sdf {prop}_test {outputfilename} ttt")
        # y_predict2 = lr_inverter.inspection(f"{test_prefix}_desc_norm.csv", normalized_dataset_filename2, LR_filename2, x_hat2, std_eps)

        # print("Inspection value:  ",  y_predict2[0])
        # print("Inspection value (scaled):  ", y_predict2[0] * (y_max2 - y_min2) + y_min2)

        # test_prefix = output_prefix + "_test_tmp3"
        # subprocess.run([fv_gen_name, f"{prop3}.sdf", f"{prop3}_test", outputfilename, test_prefix],
        #                stdout=subprocess.DEVNULL)
        # # os.system(f"{fv_gen_name} {prop}.sdf {prop}_test {outputfilename} ttt")
        # y_predict3 = lr_inverter.inspection(f"{test_prefix}_desc_norm.csv", normalized_dataset_filename3, LR_filename3, x_hat3, std_eps)

        # print("Inspection value:  ",  y_predict3[0])
        # print("Inspection value (scaled):  ", y_predict3[0] * (y_max3 - y_min3) + y_min3)

        ###############################

        #check if solutions satisfy original constrains   #cao
        # check(
        #     Lambda, Lambda_int, Lambda_dg,
        #     Gamma_int, Gamma_lnk, Gamma_int_ac, Gamma_ac_lnk, set_F,
        #     na, na_int, ns_int, ec_int, ac_int, ec_lnk, ac_lnk, fc, n_lnk,
        #     na_LB, na_UB, na_LB_int, na_UB_int, ns_LB_int, ns_UB_int,
        #     ec_LB_int, ec_UB_int,  ac_LB_int, ac_UB_int,
        #     ec_LB_lnk, ec_UB_lnk, ac_LB_lnk, ac_UB_lnk,
        #     fc_LB, fc_UB, n_LB_lnk,  n_UB_lnk, Code_F
        # )
        
        # ## Check the calculated descriptors
        # for dd, var in ann_descriptor_variables.items():
        #     print(dd, var.value())
        
        ## Check the feature vector of the inferred graph
        # call subprocess externally to calculate 
        # the fv of the calculated graph
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
