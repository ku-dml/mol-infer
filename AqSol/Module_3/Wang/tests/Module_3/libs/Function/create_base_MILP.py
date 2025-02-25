from .twolayered_MILP_2LMM_L import *
from .read_instance_2layer_2LMM_L import *

def create_base_MILP(MILP, instance_file, fringe_tree_file):
    set_Lambda = prepare_CG_element_info()
    
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
    bd2_LB, bd2_UB, bd3_LB, bd3_UB, dg_LB, dg_UB, \
    ac_LB_lf, ac_UB_lf, ac_LB_lf_common, ac_UB_lf_common, r_GC = read_seed_graph(instance_file)

    set_F, Lambda_ex, strF, fc_LB, fc_UB = prepare_fringe_trees(fringe_tree_file, Lambda)

    set_Lambda = prepare_Lambda(set_Lambda, Lambda)

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
    Gamma_tilde_ec_F, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF = \
        prepare_Gamma_ac(Gamma_int, Gamma_int_ac, Code_Lambda, Code_Lambda_dg)
    
    Gamma_lf_ac, ac_LB_lf, ac_UB_lf = prepare_Gamma_lf_ac(set_Lambda,
        ac_LB_lf, ac_UB_lf, ac_LB_lf_common, ac_UB_lf_common
    )

    Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F = \
        prepare_Lambda_tilde(Lambda_dg_int)

    t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, c_F, m_C, \
    head_C, tail_C, tail_F, E_C_plus, E_C_minus, \
    I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, \
    I_zero_one_plus, I_zero_one_minus, I_equal_one_plus, I_equal_one_minus, delta_i = \
        prepare_dataset_for_scheme_graph(n_star, rho, V_C, E_C, E_ge_two, E_ge_one,
        n_LB_int, n_UB_int, bl_LB, bl_UB, ch_LB, ch_UB, ell_UB, I_ge_two, I_ge_one,
        I_equal_one, I_zero_one, bd2_LB, bd3_LB, val, Lambda_star)

    Code_F, n_psi_H, deg_r_H, deg_r_hyd, beta_r, atom_r, ht_H, F_C, F_T, F_F, \
    n_C, n_T, n_F, F_Cp, F_Tp, F_Fp, set_F_E, set_F_v, \
    na_alpha_ex, alpha_r, deg_fr, v_ion, ac_psi_lf = prepare_fringe_tree(set_F,
        V_C, t_T, t_F, val, Lambda_int, Lambda_ex, Code_Lambda_int, Gamma_lf_ac)

    # Prepare variables of Appendix A.1 Selecting Core-vertices and Core-edges
    e_C, v_T, e_T, chi_T, clr_T, delta_chi_T, deg_tilde_C_plus, \
    deg_tilde_C_minus, n_G_int, rank_G = prepare_variables_selecting_core(
        t_C, k_C_tilde, k_C, t_T, m_C, n_LB_int, n_UB_int, ell_LB, ell_UB, bl_LB, bl_UB)
    
    # Add constraints of Appendix A.1 Selecting Core-vertices and Core-edges
    MILP = add_constraints_selecting_core(
        MILP,
        t_C, k_C_tilde, k_C, t_T, m_C, n_LB_int, n_UB_int, ell_LB, ell_UB, bl_LB, bl_UB,
        I_equal_one, I_equal_one_minus, I_equal_one_plus, I_ge_two, I_ge_one,
        I_ge_one_minus, I_ge_one_plus, I_zero_one, I_zero_one_minus, I_zero_one_plus, r_GC,
        e_C, v_T, e_T, delta_chi_T,
        chi_T, clr_T, deg_tilde_C_plus, deg_tilde_C_minus, n_G_int, rank_G)
    
    # Prepare variables of Appendix A.2 Constraints for Including Internal Vertices and Edges
    v_F, e_F, chi_F, clr_F, delta_chi_F, bl = prepare_variables_internal_vertices_and_edges(
        t_F, t_T, t_C, t_C_tilde, c_F, V_C, E_ge_one, E_ge_two, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C)
    
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
        na_C, na_T, na_F, na_ex_C, na_ex_T, na_ex_F, beta_ex_C, beta_ex_T, beta_ex_F, na_ex, bd_int,
        eledeg_C, eledeg_T, eledeg_F)
    
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
    fc = prepare_variables_fringe_configuration(t_C, t_T, t_F, set_F, Code_F, fc_LB, fc_UB)
    
    # Add constraints of A.11 Descriptor for the Number of Edge-configurations
    MILP = add_constraints_fringe_configuration(
        MILP,
        t_C, t_T, t_F, set_F, Code_F, delta_fr_C, delta_fr_T, delta_fr_F, fc)

    # A12
    MILP, mass_n = add_constraints_mass_n(MILP, n_LB, n_star, n_G, na_UB, na_ex, MASS)
    
    return strF, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac,\
        n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G, mass_n, \
        t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, \
        head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F, \
        v_T, v_F, alpha_C, alpha_T, alpha_F, \
        beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, \
        e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        E_C, ch_LB, ch_UB, \
        set_F_v, set_F_E