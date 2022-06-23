from cyclic_graphs_MILP_ec_id import *
from read_instance_BH_cyclic import *
import ann_inverter
from NeosClient import solveNeosMILP

import subprocess

####################################################
# IMPORTANT:
# Please specify the path of the cplex solver here
CPLEX_PATH= \
"/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux/cplex"
# "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/cplex"
# "/opt/cplex_12.10/cplex/bin/x86-64_linux/cplex"
# "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex"
# "/opt/cplex_12.10/cplex/bin/x86-64_linux/cplex"

CPLEX_MSG = False
CPLEX_TIMELIMIT = 0
####################################################

def main(argv):
    start = time.time()

    set_Lambda = prepare_CG_element_info()

    try:
        prop, target_value, instance_file, outputfileprefix, solver_type =\
            sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5])
    except:
        sys.stderr.write('''
usage: {} (ANNfile_prefix) (target_value) (chemical_specification_file) (outputfile_prefix) (solver_type)

  - For ANNfile, you need *_desc.csv, *_biases.txt and *_weights.txt
  - Solver type: 1 - CPLEX, 2 - CBC

'''.format(sys.argv[0]))
        exit(1)

    # for bash script
    if len(sys.argv) >= 7:
        if solver_type == 1:
            CPLEX_PATH = sys.argv[6]
        elif solver_type == 3:  # NEOS solver
            USER_EMAIL = sys.argv[6]

    #input file of bias
    ann_bias_filename = "{}_biases.txt".format(prop) 
    #input file of weight
    ann_weights_filename = "{}_weights.txt".format(prop)  
     #input file of fv
    ann_training_data_filename = "{}_desc.csv".format(prop)

    training_data = ann_inverter.read_training_data(
        ann_training_data_filename)
    des = ann_inverter.read_fv_descriptors(ann_training_data_filename)
    des = des[1:] # Drop the "CID" descriptor
    weights, biases = ann_inverter.read_trained_ANN(ann_weights_filename,
                                                    ann_bias_filename)
    ann = ann_inverter.ANN(weights, biases)
    milp_ann_constants = ann_inverter.initialize_constants(
        ann, training_data)
    ann_a, ann_b, ann_b_hat, ann_c, ann_z = milp_ann_constants
    milp_ann_vars = ann_inverter.initialize_lp_variables(
        ann, ann_a, ann_b)
    ann_descriptor_variables = ann_inverter.get_input_layer_variables(
        ann, milp_ann_vars, des)
    
    V_C, E_C, E_ge_two, E_ge_one, E_zero_one, E_equal_one, \
    I_ge_two, I_ge_one, I_zero_one, I_equal_one, \
    ell_LB, ell_UB, cs_LB, cs_UB, bl_LB, bl_UB, ch_LB, ch_UB, \
    Lambda, Lambda_dg_co, Lambda_dg_nc, Gamma_co, Gamma_in, Gamma_ex, \
    ns_LB, ns_UB, ns_LB_co, ns_UB_co, ns_LB_nc, ns_UB_nc, Lambda_star, \
    bd2_LB, bd2_UB, bd3_LB, bd3_UB, dg4_UB_nc, \
    n_LB, n_star, rho, r_Gc, dg_LB, dg_UB, \
    na_LB, na_UB, na_LB_co, na_UB_co, na_LB_nc, na_UB_nc, \
    Lambda_co, Lambda_nc, ec_LB_co, ec_UB_co, ec_LB_in, ec_UB_in, ec_LB_ex, ec_UB_ex, \
    Gamma_co_ac, Gamma_in_ac, Gamma_ex_ac, \
    ac_LB_co, ac_UB_co, ac_LB_in, ac_UB_in, ac_LB_ex, ac_UB_ex = \
        read_seed_graph(instance_file)

    d_max = 4
    if dg4_UB_nc == 0:
        d_max = 3

    set_Lambda = prepare_Lambda(set_Lambda, Lambda)

    Code_Lambda, MAX_CODE, \
    set_Lambda_dg, Lambda_dg, Code_Lambda_dg, MAX_CODE_dg, \
    val, mass, epsilon, epsilon_dg, \
    Code_Lambda_co, Code_Lambda_nc, MAX_CODE_co, MAX_CODE_nc, \
    Code_Gamma_ec_co, Code_Gamma_ec_in, Code_Gamma_ec_ex, Code_Gamma_ac_co, Code_Gamma_ac_in, Code_Gamma_ac_ex = \
        prepare_Lambda_dg(set_Lambda, Lambda, Lambda_co, Lambda_nc, 
            Gamma_co, Gamma_in, Gamma_ex, Gamma_co_ac, Gamma_in_ac, Gamma_ex_ac)


    t_C, t_T, t_F, k_C, k_C_tilde, c_F, m_C, head_C, tail_C, tail_F, E_C_plus, E_C_minus, \
    n_C, n_T, n_F, delta_i, Cld_T, Cld_C, Cld_F, prt_T, prt_C, prt_F, P_prc_T, P_prc_C, P_prc_F, Dsn_T, Dsn_C, Dsn_F, \
    j_T, j_C, j_F, t_C_tilde, I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, I_zero_one_plus, I_zero_one_minus, \
    I_equal_one_plus, I_equal_one_minus, m_UB_co, m_UB_nc, m_UB, \
    tie_breaker_T, tie_breaker_C, tie_breaker_F = prepare_dataset_for_scheme_graph(
        n_star, rho, V_C, E_C, E_ge_two, E_ge_one, cs_LB, cs_UB, bl_LB, bl_UB, ch_LB, ch_UB, ell_UB,
        I_ge_two, I_ge_one, I_equal_one, I_zero_one, bd2_LB, bd3_LB,
        val, Lambda_star, d_max, r_Gc)


    Gamma_co_less, Gamma_co_equal, Gamma_co_great, \
        Gamma_co_ac_less, Gamma_co_ac_equal, Gamma_co_ac_great, \
        Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, \
        Gamma_tilde_ac_F, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, Gamma_tilde_ac_C_ex, Gamma_tilde_ac_T_ex, Gamma_tilde_ac_F_ex, \
        Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, \
        Gamma_tilde_ec_F, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF, Gamma_tilde_ec_C_ex, Gamma_tilde_ec_T_ex, Gamma_tilde_ec_F_ex = \
        prepare_Gamma_ac(
            Gamma_co, Gamma_in, Gamma_ex, Gamma_co_ac, Gamma_in_ac, Gamma_ex_ac, Code_Lambda, Code_Lambda_dg)

    Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F, Lambda_tilde_dg_C_nc, Lambda_tilde_dg_T_nc, Lambda_tilde_dg_F_nc = \
        prepare_Lambda_tilde(Lambda_dg_co, Lambda_dg_nc)

    MILP = pulp.LpProblem(name="MILP_cyclic_Branch")

    # Prepare variables of Appendix A.1 Selecting Core-vertices and Core-edges
    e_C, v_T, e_T, chi_T, clr_T, delta_chi_T, \
    deg_tilde_C_plus, deg_tilde_C_minus, cs = prepare_variables_selecting_core(
        t_C, t_T, k_C, m_C, ell_LB, ell_UB, bl_LB, bl_UB, cs_LB, cs_UB)

    # Add constraints of Appendix A.1 Selecting Core-vertices and Core-edges
    MILP = add_constraints_selecting_core(
        MILP,
        t_C, t_T, k_C,
        I_equal_one, I_equal_one_minus, I_equal_one_plus,
        I_ge_two, I_ge_one, I_ge_one_minus, I_ge_one_plus,
        I_zero_one_minus,
        I_zero_one_plus,
        e_C, v_T, e_T, delta_chi_T,
        chi_T, clr_T, deg_tilde_C_plus, deg_tilde_C_minus, cs)
    
    # Prepare variables of Appendix A.2 Constraints for Including Internal Vertices and Edges
    bl_G, v_F, e_F, chi_F, clr_F, delta_chi_F, \
    sigma, sigma_C, sigma_T, sigma_F, bl = prepare_variables_internal_vertices_and_edges(
        t_F, t_T, t_C, t_C_tilde, c_F, V_C, E_ge_one, E_ge_two, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C)

    # Add constraints of Appendix A.2 Constraints for Including Internal Vertices and Edges
    MILP = add_constraints_internal_vertices_and_edges(
        MILP,
        t_T, t_F, t_C, t_C_tilde, c_F, tail_F, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C,
        delta_chi_F, delta_chi_T, e_F, v_F, v_T, sigma, sigma_C, sigma_T, sigma_F,
        chi_F, clr_F, bl, bl_G)

    # Prepare variables of Appendix A.3 Constraints for Including Fringe-trees
    n_G, ch_G, v_T, v_C, v_F, h_T, h_C, h_F, sigma = \
        prepare_variables_fringe_trees(n_LB, n_star, rho, ch_UB, t_T, t_C, t_F, n_T, n_C, n_F,
            delta_i, I_ge_two, I_ge_one, V_C, E_ge_one, E_ge_two, v_T, v_F, sigma)

    # Add constraints of Appendix A.3 Constraints for Including Fringe-trees
    MILP = add_constraints_fringe_trees(
        MILP,
        t_T, t_F, t_C, t_C_tilde, c_F, I_ge_one, I_ge_two, rho, n_star,
        j_T, j_F, j_C, n_T, n_C, n_F, Cld_T, Cld_C, Cld_F, P_prc_T, P_prc_C, P_prc_F, delta_i, ch_LB, ch_UB, E_C,
        delta_chi_F, delta_chi_T, e_F, v_F, v_C, v_T, h_T, h_C, h_F, sigma, sigma_C, sigma_T, sigma_F,
        chi_F, clr_F, n_G, ch_G
    )

    # Prepare variables of A.4 Descriptor for the Number of Specified Degree
    deg_C, deg_T, deg_F, delta_dg_C, delta_dg_T, delta_dg_F, dg, dg_co, \
    dg_C, dg_T, dg_nc, dg_in, dg_ex, dg_Cp, dg_Tp, dg_Fp = \
        prepare_variables_degree(t_C, t_T, t_F, n_C, n_T, n_F, delta_i, n_star, cs_LB, cs_UB, dg_LB, dg_UB, rho)

    # Add constraints of A.4 Descriptor for the Number of Specified Degree
    MILP = add_constraints_degree(
        MILP,
        t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, Cld_T, Cld_C, Cld_F, Dsn_T, Dsn_C, Dsn_F, delta_i, dg4_UB_nc,
        I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, rho,
        delta_dg_C, delta_dg_T, delta_dg_F, e_T, e_F, v_F, v_C, v_T, delta_chi_T, delta_chi_F,
        deg_C, deg_T, deg_F, dg, dg_co, dg_C, dg_T, dg_nc, dg_in, dg_ex, dg_Cp, dg_Tp, dg_Fp,
        deg_tilde_C_minus, deg_tilde_C_plus
    )

    # Prepare variables of A.5 Assigning Multiplicity
    beta_T, beta_F, beta_C, beta_plus, beta_minus, beta_in, delta_beta_T, delta_beta_F, delta_beta_C, \
        delta_beta_plus, delta_beta_minus, delta_beta_in, bd, bd_co, bd_in, bd_ex, \
        bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF, bd_Cp, bd_Tp, bd_Fp = \
            prepare_variables_multiplicity(t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, c_F, delta_i,
                I_equal_one, I_zero_one, I_ge_one, I_ge_two, E_C,
                bd2_LB, bd2_UB, bd3_LB, bd3_UB, m_UB, m_UB_co, m_UB_nc, n_star, rho)
            
    # Add constraints of A.5 Assigning Multiplicity
    MILP = add_constraints_multiplicity(
        MILP,
        t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, c_F, delta_i, head_C, tail_C, Dsn_C, Dsn_T, Dsn_F,
        I_equal_one, I_zero_one, I_ge_one, I_ge_two, E_C, bd2_LB, bd2_UB, bd3_LB, bd3_UB, rho,
        e_C, e_T, e_F, v_F, v_C, v_T,
        delta_chi_T, delta_chi_F, delta_beta_C, delta_beta_T, delta_beta_F, delta_beta_plus, delta_beta_minus, delta_beta_in,
        beta_C, beta_T, beta_F, beta_plus, beta_minus, beta_in,
        bd, bd_co, bd_in, bd_ex, bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF, bd_Cp, bd_Tp, bd_Fp
    )

    # Prepare variables of A.6 Assigning Chemical Elements and Valence Condition
    beta_CT, beta_TC, beta_CF, beta_TF, alpha_C, alpha_T, alpha_F, delta_alpha_C, delta_alpha_T, delta_alpha_F, MASS, n_H, \
        na, na_co, na_C, na_T, na_nc, na_in, na_ex, na_Cp, na_Tp, na_Fp = \
        prepare_variables_chemical_elements(t_C, t_T, t_F, n_C, n_T, n_F, delta_i, n_star, 
            MAX_CODE, MAX_CODE_co, MAX_CODE_nc,
            Lambda, Lambda_co, Lambda_nc, epsilon,
            na_LB, na_UB, na_LB_co, na_UB_co, na_LB_nc, na_UB_nc, rho)

    # Add constraints of A.6 Assigning Chemical Elements and Valence Condition
    MILP = add_constraints_chemical_elements(
        MILP,
        t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, c_F, delta_i, E_C,
        I_equal_one, I_zero_one, I_ge_one, I_ge_one_plus, I_ge_one_minus, I_ge_two, I_ge_two_plus, I_ge_two_minus,
        Cld_C, Cld_T, Cld_F, Dsn_C, Dsn_T, Dsn_F, val, mass,
        Lambda, Lambda_co, Lambda_nc, Code_Lambda, Code_Lambda_co, Code_Lambda_nc, Lambda_star, epsilon, rho,
        e_T, e_F, delta_chi_T, delta_chi_F, delta_alpha_C, delta_alpha_T, delta_alpha_F,
        deg_C, deg_T, deg_F, beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF,
        beta_plus, beta_minus, beta_in,
        alpha_C, alpha_T, alpha_F, MASS, n_H,
        na, na_co, na_C, na_T, na_nc, na_in, na_ex, na_Cp, na_Tp, na_Fp, bd_co, bd_in, bd_ex
    )

    # Prepare variables of A.7 Constraints for Bounds on the Number of Bonds
    bd_T = prepare_variables_number_of_bounds(k_C, t_T, bd_T)

    # Add constraints of A.7 Constraints for Bounds on the Number of Bonds
    MILP = add_constraints_number_of_bounds(
        MILP,
        t_T, k_C, E_C, I_equal_one, I_zero_one, bd2_LB, bd2_UB, bd3_LB, bd3_UB,
        delta_chi_T, delta_beta_C, delta_beta_T, delta_beta_plus, delta_beta_minus, bd_T
    )

    # Prepare variables of A.8 Descriptors for the Number of Adjacency-configuration
    ac_C, ac_T, ac_F, ac_Cp, ac_Tp, ac_Fp, ac_CT, ac_TC, ac_CF, ac_TF, \
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_Cp, delta_ac_Tp, delta_ac_Fp, \
        ac_co, ac_in, ac_ex, ac_nc, \
        delta_alpha_CT, delta_alpha_TC, delta_alpha_CF, delta_alpha_TF, \
        alpha_CT, alpha_TC, alpha_CF, alpha_TF, \
        Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, \
        Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus, \
        Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus, Delta_ac_TF_minus, \
        delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF = \
            prepare_variables_adjacency_configuration(
                t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F, n_T, n_C, n_F,
                delta_i, Dsn_C, Dsn_T, Dsn_F, rho, m_UB,
                Gamma_co_ac, Gamma_in_ac, Gamma_ex_ac,
                Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F,
                Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF,
                Gamma_tilde_ac_C_ex, Gamma_tilde_ac_T_ex, Gamma_tilde_ac_F_ex,
                ac_LB_co, ac_UB_co, ac_LB_in, ac_UB_in, ac_LB_ex, ac_UB_ex,
                Lambda_co, Lambda_nc, MAX_CODE_co, MAX_CODE_nc
            )
            
    # Add constraints of A.8 Descriptors for the Number of Adjacency-configuration
    MILP = add_constraints_adjacency_configuration(
        MILP,
        t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, k_C, k_C_tilde, m_C, c_F, tail_C, head_C,
        Dsn_C, Dsn_T, Dsn_F, prt_C, prt_T, prt_F,
        Gamma_co_ac, Gamma_in_ac, Gamma_ex_ac, Gamma_co_ac_less, Gamma_co_ac_equal,
        Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F,
        Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF,
        Gamma_tilde_ac_C_ex, Gamma_tilde_ac_T_ex, Gamma_tilde_ac_F_ex,
        rho, delta_i, Lambda_co, Lambda_nc,
        Code_Lambda_co, Code_Lambda_nc, MAX_CODE_co, MAX_CODE_nc,
        delta_alpha_C, delta_alpha_T, delta_alpha_F, delta_beta_C, delta_beta_T, delta_beta_F,
        delta_beta_plus, delta_beta_minus, delta_beta_in, delta_chi_T, delta_chi_F,
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_Cp, delta_ac_Tp, delta_ac_Fp,
        delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,
        e_C, e_T, e_F, v_C, v_T, v_F,
        delta_alpha_CT, delta_alpha_TC, delta_alpha_CF, delta_alpha_TF,
        Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus,
        Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus,
        Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus, Delta_ac_TF_minus,
        alpha_C, alpha_T, alpha_F, alpha_CT, alpha_TC, alpha_CF, alpha_TF,
        beta_C, beta_T, beta_F, beta_plus, beta_minus, beta_in,
        ac_C, ac_T, ac_F, ac_Cp, ac_Tp, ac_Fp,
        ac_CT, ac_TC, ac_CF, ac_TF, ac_co, ac_in, ac_ex, ac_nc
    )

    # Prepare variables of A.9 Descriptor for the Number of Chemical Symbols
    ns_co, ns_nc, delta_ns_C, delta_ns_T, delta_ns_F = \
        prepare_variables_chemical_symbols(
            t_C, t_T, t_F, n_C, n_T, n_F, delta_i,
            n_star, cs_LB, cs_UB, Lambda_dg_co, Lambda_dg_nc, epsilon_dg
        )

    # Add constraints of A.9 Descriptor for the Number of Chemical Symbols
    MILP = add_constraints_chemical_symbols(
        MILP,
        t_C, t_T, t_F, n_C, n_T, n_F,
        Lambda_dg_co, Lambda_dg_nc,
        Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F,
        Lambda_tilde_dg_C_nc, Lambda_tilde_dg_T_nc, Lambda_tilde_dg_F_nc,
        Code_Lambda, Code_Lambda_co, Code_Lambda_nc, epsilon_dg, delta_i,
        delta_ns_C, delta_ns_T, delta_ns_F,
        alpha_C, alpha_T, alpha_F, deg_C, deg_T, deg_F, ns_co, ns_nc
    )

    # Prepare variables of A.10 Descriptor for the Number of Edge-configurations
    ec_C, ec_T, ec_F, ec_Cp, ec_Tp, ec_Fp, ec_CT, ec_TC, ec_CF, ec_TF, \
        delta_ec_C, delta_ec_T, delta_ec_F, delta_ec_Cp, delta_ec_Tp, delta_ec_Fp, \
        delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF, ec_co, ec_in, ec_ex, ec_nc, \
        delta_dg_CT, delta_dg_TC, delta_dg_CF, delta_dg_TF, \
        deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF, \
        Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, Delta_ec_F_plus, Delta_ec_F_minus, \
        Delta_ec_CT_plus, Delta_ec_CT_minus, Delta_ec_TC_plus, Delta_ec_TC_minus, \
        Delta_ec_CF_plus, Delta_ec_CF_minus, Delta_ec_TF_plus, Delta_ec_TF_minus = \
            prepare_variables_edge_configuration(
                t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F,
                n_T, n_C, n_F, delta_i, Dsn_C, Dsn_T, Dsn_F, rho, m_UB,
                Gamma_co, Gamma_in, Gamma_ex,
                Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_F, Gamma_tilde_ec_CT,
                Gamma_tilde_ec_TC, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF,
                Gamma_tilde_ec_C_ex, Gamma_tilde_ec_T_ex, Gamma_tilde_ec_F_ex,
                ec_LB_co, ec_UB_co, ec_LB_in, ec_UB_in, ec_LB_ex, ec_UB_ex
            )

    # Add constraints of A.10 Descriptor for the Number of Edge-configurations
    MILP = add_constraints_edge_configuration(
        MILP,
        t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, k_C, k_C_tilde, m_C, c_F,
        tail_C, head_C, Dsn_C, Dsn_T, Dsn_F, prt_C, prt_T, prt_F,
        Gamma_co, Gamma_in, Gamma_ex, Gamma_co_less, Gamma_co_equal,
        Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_F, Gamma_tilde_ec_CT,
        Gamma_tilde_ec_TC, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF,
        Gamma_tilde_ec_C_ex, Gamma_tilde_ec_T_ex, Gamma_tilde_ec_F_ex,
        Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F, Gamma_tilde_ac_CT,
        Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF,
        Gamma_tilde_ac_C_ex, Gamma_tilde_ac_T_ex, Gamma_tilde_ac_F_ex,
        rho, delta_i,
        Code_Gamma_ac_co, Code_Gamma_ac_in, Code_Gamma_ac_ex,
        Code_Gamma_ec_co, Code_Gamma_ec_in, Code_Gamma_ec_ex,
        delta_dg_C, delta_dg_T, delta_dg_F, delta_chi_T, delta_chi_F,
        delta_beta_C, delta_beta_T, delta_beta_F, delta_beta_plus, delta_beta_minus, delta_beta_in,
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_Cp, delta_ac_Tp, delta_ac_Fp,
        delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,
        delta_ec_C, delta_ec_T, delta_ec_F, delta_ec_Cp, delta_ec_Tp, delta_ec_Fp,
        delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF,
        e_C, e_T, e_F, v_C, v_T, v_F, delta_dg_CT, delta_dg_TC, delta_dg_CF, delta_dg_TF,
        deg_C, deg_T, deg_F,
        deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF,
        ec_C, ec_T, ec_F, ec_Cp, ec_Tp, ec_Fp, ec_CT, ec_TC, ec_CF, ec_TF,
        ec_co, ec_in, ec_ex, ec_nc,
        Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, Delta_ec_F_plus, Delta_ec_F_minus,
        Delta_ec_CT_plus, Delta_ec_CT_minus, Delta_ec_TC_plus, Delta_ec_TC_minus,
        Delta_ec_CF_plus, Delta_ec_CF_minus, Delta_ec_TF_plus, Delta_ec_TF_minus
    )

    # Add constraints of fringe tree tie breaker
    # MILP = add_constraints_tie_breaker(
    #     MILP,
    #     t_C, t_T, t_F, tie_breaker_C, tie_breaker_T, tie_breaker_F, delta_i,
    #     alpha_C, alpha_T, alpha_F, beta_C, beta_T, beta_F
    # )

    descriptors, num_fv, desc_index_v, desc_index_e, mass_ind = prepare_fv(
        ann_training_data_filename,
        Lambda_dg_co, Lambda_dg_nc, Gamma_co, Gamma_in, Gamma_ex,
        n_G, cs, ch_G, bl_G, MASS, bd_co, bd_in, bd_ex, n_H,
        dg_co, dg_nc, ns_co, ns_nc, ec_co, ec_in, ec_ex
    )

    eps = 0.02
    MILP = ann_inverter.build_MILP_ReLU(
        MILP,
        ann,
        milp_ann_vars,
        milp_ann_constants,
        target_value,
        eps
    )

    _, y, _ = milp_ann_vars

    MILP, mass_n = add_constraints_mass_n(MILP, n_LB, n_star, n_G, MASS)

    MILP = add_constraints_ANN(MILP,
                        descriptors,
                        num_fv,
                        y, mass_ind,
                        mass_n)

    # Output all MILP variables and constraints
    # MILP.writeLP("./MILP_cyclic_Branch.lp")

    init_end = time.time()
    print("Initializing Time:", "{:.3f}".format(init_end - start),
            file=sys.stderr, flush = True)

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
    
    # print("Number of variables:", num_vars)
    # print(" - Integer :", num_ints)
    # print(" - Binary  :", len(bins))
    # print("Number of constraints:", num_constraints)


    # Solve MILP
    if solver_type == 1:
        if CPLEX_TIMELIMIT > 0:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH, 
                               msg = CPLEX_MSG,
                               timeLimit = CPLEX_TIMELIMIT)
        else:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG)
        print("Start Solving Using CPLEX...",
                file = sys.stderr, flush = True)        
        MILP.solve(CPLEX)
        solve_end = time.time()
    elif solver_type == 2:
        print("Start Solving Using Coin-OR...",
                file = sys.stderr, flush = True)
        MILP.solve()
        solve_end = time.time()
    elif solver_type == 3: # online NEOS solver
        print("Start solving using NEOS...")
        lp_filename = outputfileprefix+'.lp'
        MILP.writeLP(lp_filename)
        variables_dict = solveNeosMILP(USER_EMAIL, lp_filename)
        x = variables_dict['x']
        y = variables_dict['y']
        n_G = variables_dict['n_G']
        v_C = variables_dict['v_C']
        v_T = variables_dict['v_T']
        v_F = variables_dict['v_F']
        alpha_C = variables_dict['alpha_C']
        alpha_T = variables_dict['alpha_T']
        alpha_F = variables_dict['alpha_F']
        beta_C = variables_dict['beta_C']
        beta_T = variables_dict['beta_T']
        beta_F = variables_dict['beta_F']
        beta_CT = variables_dict['beta_CT']
        beta_TC = variables_dict['beta_TC']
        beta_CF = variables_dict['beta_CF']
        beta_TF = variables_dict['beta_TF']
        delta_chi_T = variables_dict['delta_chi_T']
        delta_chi_F = variables_dict['delta_chi_F']
        e_T = variables_dict['e_T']
        e_F = variables_dict['e_F']
        chi_T = variables_dict['chi_T']
        solve_end = time.time()
    else:
        raise Exception("Unknown solver type: " + solver_type)

    layer_num = len(weights)

    # Output result of solving MILP to command-line
    # print("# -------- solve status --------")
    if pulp.LpStatus[MILP.status] == "Optimal" or solver_type == 3:
        output_status = "Feasible"
        y_star = (y[(layer_num + 1, 1)].value())        
        print("Status:", output_status)
        print("MILP y*:", "{:.3f}".format(y_star))
    else:
        output_status = pulp.LpStatus[MILP.status]        
        print("Status:", output_status)
    
      
  
    # ############################################
    # # The following block of code is used to print out the value of feature vector #
    # # and the value of all variables used in MILP in to a file 
    # # "{outputfileprefix}_test_all.txt" where outputfileprefix is 
    # # an input parameter
    # ############################################
    # #
    # with open(outputfileprefix + "_test_all.txt", "w", newline='\n') as f:
    #     # f.write("######### Feature Vector ############\n")
    #     # for i in range(1, num_fv):
    #     #     f.write(stringoutput[i] + str(y[(1, i)].value()) + "\n")
    #     # f.write("\n\n\n")
    #     # f.write("#####################################\n")
    
    #     for var in MILP.variables():
    #         if type(var) is dict:
    #             for x in var:
    #                 if var[x].value() is None:
    #                     pass
    #                 elif var[x].value() > 0:
    #                     f.write(f"{var[x].name}: {var[x].value()}\n")
    #         elif var.value() is None:
    #             pass
    #         elif var.value() > 0:
    #             f.write(f"{var.name}: {var.value()}\n")
    
    #     f.write("\n")
    #
    # #############################################

    if MILP.status == 1 or solver_type == 3:

        # Output SDF file
        outputfilename = outputfileprefix + ".sdf"

        index_C, index_T, graph_adj = print_sdf_file(
            t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, delta_i, Lambda, Lambda_co, Lambda_nc,
            prt_T, prt_C, prt_F, head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two,
            n_G, v_C, v_T, v_F, alpha_C, alpha_T, alpha_F,
            beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, delta_chi_T, delta_chi_F,
            e_T, e_F,
            outputfilename
        )

        # Output file of partition which will be used in graph generation
        outputfilename = outputfileprefix + "_partition.txt"

        print_gstar_file(
            n_G, chi_T,
            t_C, t_T, index_C, index_T, graph_adj, E_C, ch_LB, ch_UB,
            I_ge_one, I_ge_two, I_zero_one, I_equal_one,
            outputfilename
        )
        
        # ## Check the calculated descriptors
        # for dd, var in ann_descriptor_variables.items():
        #     print(dd, var.value())
        
        #x_star = list()
        
        #for desc_name in des:
        #    x_star.append(ann_descriptor_variables[desc_name].value())
            
        #y_star = ann.propagate(x_star)[0]
        
        #print("ANN propagated y*:", "{:.3f}".format(y_star))
        
    print("Solving Time:", "{:.3f}".format(solve_end - init_end))
        


if __name__ == '__main__':
    main(sys.argv)
