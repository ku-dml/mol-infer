"""
Copyright (C) 2020
by Discrete Mathematics Lab, 
   Department of Applied Mathematics and Physics, 
   Graduate School of Informatics,
   Kyoto University

Licensed under the MIT license, see License.txt
"""

"""
infer_acyclic_graphs.py

This file implements functions that given 
files that contain the weight and bias values of
a trained ANN, solve the inverse problem to
get a vector of resources for constructing 
an acyclic graph
"""

import os
from acyclic_graphs_MILP import *

####################################################
# IMPORTANT:
# Please specify the path of the cplex solver here
CPLEX_PATH= \
"/opt/ibm/cplex_12.10/cplex/bin/x86-64_linux/cplex"
#####################################################

def main(argv):
    start = time.time()

    target_value = float(argv[1])   # target value
    n_star = int(argv[2])
    dia_star = int(argv[3])
    k_star = int(argv[4])           # only can be 2
    d_max = int(argv[5])            # only can be 3 or 4
    bl_star = int(argv[6])
    bh_star = int(argv[7])
    solver_type = int(argv[8])      # 1 for CPLEX, 2 for COIN-OR

    prop = argv[9]                  # name of the chemical property

    # for bash script
    if len(sys.argv) >= 11:
        CPLEX_PATH = sys.argv[10]

    #input file of bias
    ann_bias_filename = "{}_biases.txt".format(prop) 
    #input file of weight
    ann_weights_filename = "{}_weights.txt".format(prop)  
     #input file of fv
    ann_training_data_filename = "{}_desc.csv".format(prop)
    
    #AD min and max filename used in AD
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ad_min_filename = os.path.join(script_dir, "AD_min.txt")
    ad_max_filename = os.path.join(script_dir, "AD_max.txt")
    
    training_data = ann_inverter.read_training_data(
        ann_training_data_filename)
    des = ann_inverter.read_fv_descriptors(ann_training_data_filename)
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
    
    MILP = pulp.LpProblem(name="MILP_Acyclic_Branch")
    
    dia_star, d_max, k_star, bh_star, t_star, s_star,\
    c_star, root_num, bl_star, E_B, E1_plus, E1_minus, tail, head,\
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,\
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S, \
    P_S_prc, P_T_prc, Bcset = prepare_dataset_for_scheme_graph(
        n_star, d_max, dia_star, k_star, bl_star, bh_star)
    
    MAX_BOND = 3
    MAX_VAL = 4
    
    set_Lambda = prepare_CG_element_info()
    set_Lambda = prepare_CG_from_fv(set_Lambda, ann_training_data_filename)
    Lambda = [ele.symbol for ele in set_Lambda]
    Extra_Lambda = Lambda[:]
    Extra_Lambda.append("e")
    Code_Lambda = {ele: i + 1 for i, ele in enumerate(Lambda)}
    Code_Lambda["e"] = 0
    MAX_CODE = len(Lambda)
    
    val = {ele.symbol: ele.valence for ele in set_Lambda}
    mass = {ele.symbol: ele.mass for ele in set_Lambda}
    Gamma_equal = {(atom, atom, k) for atom in Lambda 
        for k in range(1, MAX_BOND + 1) if k <= val[atom] - 1}
    Gamma_less = {(atom1, atom2, k) for atom1 in Lambda 
        for atom2 in Lambda for k in range(1, MAX_BOND + 1)
        if mass[atom1] < mass[atom2] and k <= min(val[atom1], val[atom2])
                  and k <= max(val[atom1], val[atom2]) - 1}
    Gamma_great = {(atom2, atom1, k) for (atom1, atom2, k) in Gamma_less}
    Gamma = Gamma_equal | Gamma_less | Gamma_great
    Gamma_zero = {(atom1, atom2, 0) for atom1 in Extra_Lambda 
        for atom2 in Extra_Lambda}
    Gamma_ALL = Gamma | Gamma_zero

    # # Nc constants 200618
    # a_1 = {(atom, 1) for atom in Lambda}
    # a_2 = {(atom, 2) for atom in Lambda if val[atom] >= 2}
    # a_3 = {(atom, 3) for atom in Lambda if val[atom] >= 3}
    # a_1_1 = {(atom, 1, 1) for atom in Lambda if val[atom] >= 2}
    # a_1_2 = {(atom, 1, 2) for atom in Lambda if val[atom] >= 3}
    # a_1_3 = {(atom, 1, 3) for atom in Lambda if val[atom] >= 4}
    # a_2_2 = {(atom, 2, 2) for atom in Lambda if val[atom] >= 4}
    # a_1_1_1 = {(atom, 1, 1, 1) for atom in Lambda if val[atom] >= 3}
    # a_1_1_2 = {(atom, 1, 1, 2) for atom in Lambda if val[atom] >= 4}
    # a_1_1_1_1 = {(atom, 1, 1, 1, 1) for atom in Lambda if val[atom] >= 4}

    # Ncset = a_1 | a_2 | a_3 | a_1_1 | a_1_2 | a_1_3 | a_2_2 |\
    #         a_1_1_1 | a_1_1_2 | a_1_1_1_1

    # Ncset_e = Ncset.copy()
    # Ncset_e.add('eps')

    # ce_nc = {nu: Code_Lambda[nu[0]] for nu in Ncset}
    # ce_nc["eps"] = 0

    # dg_nc = {nu: len(nu) - 1 for nu in Ncset}
    # dg_nc["eps"] = 0

    # ml = {nu: sum(nu[1:]) for nu in Ncset}
    # ml["eps"] = 0
    
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus = \
        prepare_variables_4_2(
            n_star, dia_star, k_star, d_max, bl_star, bh_star, 
            s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
            n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
            V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S, 
            P_S_prc, P_T_prc, Bcset,
            MAX_BOND, MAX_VAL, MAX_CODE,
            Lambda, Extra_Lambda, Code_Lambda, val, mass, 
            Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_2(MILP, 
        n_star, dia_star, k_star, d_max, bl_star, bh_star, 
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass, 
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus
        )
    tilde, u, v, e = prepare_variables_4_3(
        n_star, dia_star, k_star, d_max, bl_star, bh_star, 
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass, 
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_3(MILP,
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus,
        tilde, u, v, e
        )
    alpha_tilde = prepare_variables_4_4(
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    beta_tilde, beta_hat = prepare_variables_4_5(
        n_star, dia_star, k_star, d_max, bl_star, bh_star, 
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass, 
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_5(MILP,
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus,
        tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat
        )
    delta_alpha, delta_beta_tilde, delta_beta_hat = \
        prepare_variables_4_6(
            n_star, dia_star, k_star, d_max, bl_star, bh_star, 
            s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
            n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
            V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
            P_S_prc, P_T_prc, Bcset,
            MAX_BOND, MAX_VAL, MAX_CODE,
            Lambda, Extra_Lambda, Code_Lambda, val, mass, 
            Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_6(MILP,
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus,
        tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
        delta_alpha, delta_beta_tilde, delta_beta_hat
        )
    ce_in, ce_ex, ce, Mass, b_in, b_ex, b, n_H = \
        prepare_variables_4_7(
            n_star, dia_star, k_star, d_max, bl_star, bh_star, 
            s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
            n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
            V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
            P_S_prc, P_T_prc, Bcset,
            MAX_BOND, MAX_VAL, MAX_CODE,
            Lambda, Extra_Lambda, Code_Lambda, val, mass, 
            Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_7(MILP,
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus,
        tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
        delta_alpha, delta_beta_tilde, delta_beta_hat, ce_in, ce_ex, ce,
        Mass, b_in, b_ex, b, n_H
        )
    deg, delta_deg, dg_in, dg_ex = prepare_variables_4_8(
        n_star, dia_star, k_star, d_max, bl_star, bh_star, 
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass, 
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_8(MILP,
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus,
        tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
        delta_alpha, delta_beta_tilde, delta_beta_hat, ce_in, ce_ex, ce,
        Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex
        )
    delta_tau, delta_tau_hat, ac_in, ac_ex, ac = \
        prepare_variables_4_9(
            n_star, dia_star, k_star, d_max, bl_star, bh_star, 
            s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
            n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
            V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
            P_S_prc, P_T_prc, Bcset,
            MAX_BOND, MAX_VAL, MAX_CODE,
            Lambda, Extra_Lambda, Code_Lambda, val, mass, 
            Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_9(MILP,
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus,
        tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
        delta_alpha, delta_beta_tilde, delta_beta_hat, ce_in, ce_ex, ce,
        Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex,
        delta_tau, delta_tau_hat, ac_in, ac_ex, ac, n_H
        )
    bc, bc_in, bc_ex, delta_dc, delta_dc_hat = prepare_variables_4_10(
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        )
    MILP = add_constraints_4_10(MILP,
        n_star, dia_star, k_star, d_max, bl_star, bh_star,
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass,
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus,
        tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
        delta_alpha, delta_beta_tilde, delta_beta_hat, ce_in, ce_ex, ce,
        Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex,
        delta_tau, delta_tau_hat, ac_in, ac_ex, ac, n_H,
        bc, bc_in, bc_ex, delta_dc, delta_dc_hat
        )

    # nc, delta_nc = prepare_variables_B_8(
    #     n_star, dia_star, k_star, d_max, bl_star, bh_star,
    #     s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    #     n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    #     V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    #     P_S_prc, P_T_prc, Dcset, 
    #     MAX_BOND, MAX_VAL, MAX_CODE,
    #     Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    #     Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    #     a_3, a_1_3, Ncset, Ncset_e, ce_nc, dg_nc, ml
    #     )
    # MILP = add_constraints_B_8(MILP,
    #     n_star, dia_star, k_star, d_max, bl_star, bh_star,
    #     s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    #     n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    #     V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    #     P_S_prc, P_T_prc, Dcset,
    #     MAX_BOND, MAX_VAL, MAX_CODE,
    #     Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    #     Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    #     a_3, a_1_3, Ncset, Ncset_e, ce_nc, dg_nc, ml,
    #     a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    #     tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    #     delta_alpha, delta_beta_tilde, delta_beta_hat,
    #     ce_bt, ce_ft, ce, Mass, b_bt, b_ft, b, deg, delta_deg, dg,
    #     delta_tau, delta_tau_hat, ac_bt, ac_ft, ac, n_H,
    #     dc, delta_dd, delta_dd_hat,
    #     nc, delta_nc
    #     )
    
    
    descriptors, stringoutput, num_fv, ceset, acset, \
    desc_index_v, desc_index_e = prepare_fv(ann_training_data_filename,
        n_star, dia_star, k_star, d_max, bl_star, bh_star, 
        s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
        n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
        V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
        P_S_prc, P_T_prc, Bcset,
        MAX_BOND, MAX_VAL, MAX_CODE,
        Lambda, Extra_Lambda, Code_Lambda, val, mass, 
        Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
        a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
        tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
        delta_alpha, delta_beta_tilde, delta_beta_hat, ce_in, ce_ex, ce, 
        Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex,
        delta_tau, delta_tau_hat, ac_in, ac_ex, ac, n_H, 
        bc, bc_in, bc_ex, delta_dc, delta_dc_hat 
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
    
    MILP = add_constraints_ANN(MILP,
                            descriptors,
                            num_fv,
                            y, n_star)
    
    MILP = add_constraints_AD(MILP,
                           n_star,
                           num_fv, ceset, acset,
                           ce_in, ce_ex, ac_in, ac_ex, y,
                           Lambda,Gamma_equal, Gamma_less,
                           desc_index_v, desc_index_e,
                           ann_training_data_filename,
                           ad_min_filename,
                           ad_max_filename)
    
    # MILP.writeLP("./MILP_Acyclic_Branch.lp")
    # print(len(weights))
    layer_num = len(weights)
    
    if solver_type == 1:
        CPLEX = pulp.CPLEX(path=CPLEX_PATH, msg=None)
        # print("Start Solving Using CPLEX...")
        init_end = time.time()
        MILP.solve(CPLEX)
        solve_end = time.time()
    else:
        # print("Start Solving Using Coin-OR...")
        init_end = time.time()
        MILP.solve()
        solve_end = time.time()
    
    
    print("Status:" + pulp.LpStatus[MILP.status] + "\n")
    print("Initializing Time: " + str(init_end - start) + "\n")
    print("Solving Time: " + str(solve_end - init_end) + "\n")
    print("Value: " + str(y[(layer_num + 1, 1)].value()) + "\n")
    
    # strtemp = ""
    # if MILP.status == -1:
    #     strtemp = "Infeasible, "
    # else:
    #     strtemp = "Feasible, "
    # strtemp = strtemp + str(solve_end - init_end) + ", "
    # if MILP.status != -1:
    #     strtemp = strtemp + str(y[(layer_num + 1, 1)].value())
    
    # print(strtemp)

    outputfileprefix = "{}_tv{}_n{}_dia{}_k{}_dmax{}_bl{}_bh{}_solver{}".format(prop, target_value, n_star, dia_star, k_star, d_max, bl_star, bh_star, solver_type)

    # for bash script
    if len(sys.argv) >= 12:
        outputfileprefix = sys.argv[11]
    
    # ############################################
    # # The following block of code is used to print out the value of feature vector #
    # # and the value of all variables used in MILP in to a file "test.txt" #
    # ############################################
    # #
    # with open(outputfileprefix + "_test_all.txt", "w") as f:
    #     f.write("######### Feature Vector ############\n")
    #     for i in range(1, num_fv):
    #         f.write(stringoutput[i] + str(y[(1, i)].value()) + "\n")
    #     f.write("\n\n\n")
    #     f.write("#####################################\n")
    #
    #     for var in MILP.variables():
    #         if type(var) is dict:
    #             for x in var:
    #                 if var[x].value() is None:
    #                     pass
    #                 elif var[x].value() > 0:
    #                     f.write(var[x].name + ": " + str(var[x].value()) +
    #                          "\n")
    #         else:
    #             f.write(var.name + ": " + str(var.value()) + "\n")
    #
    #     f.write("\n")

    # #############################################

    if MILP.status != -1:

        outputfilename = outputfileprefix + ".txt"

        print_gstar_file(
            dia_star, d_max,
            Bcset,
            set_Lambda,
            Gamma_equal, Gamma_less,
            ce_in, ce_ex, ac_in, ac_ex, 
            bc_in, bc_ex, dg_in, dg_ex,
            outputfilename)

        outputfilename = outputfileprefix + ".sdf"

        print_sdf_file(
            n_star, E_B, tail, head,
            s_star, t_star, n_S_tree, n_T_tree,
            prt_S, prt_T,
            Lambda,
            u, v, alpha_tilde, beta_tilde, beta_hat,
            outputfilename)


if __name__ == '__main__':
    main(sys.argv)
