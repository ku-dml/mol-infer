import ANN.ann_inverter as ann_inverter
from twolayered_MILP_2LMM_L import *

def ANN_add_vars_constraints_to_MILP(config, index, base_MILP, base_var):

    try:
        prop = config["input_data"][index]["prefix"]
        target_value = config["input_data"][index]["target_value"]
    except:
        print("\tError: Please specify (prefix), (target_value) in config.yaml")
        sys.exit(1)

    print("\tprefix:", prop)
    print("\ttarget_value:", target_value, "\n")

    #input file of bias
    ann_bias_filename = "{}_biases.txt".format(prop)
    #input file of weight
    ann_weights_filename = "{}_weights.txt".format(prop)  
    #input file of fv
    ann_training_data_filename = "{}_desc.csv".format(prop)
    ann_training_data_norm_filename = "{}_desc_norm.csv".format(prop)
    # all fringe trees used in learning
    fv_fringe_tree_filename = "{}_fringe.txt".format(prop)    
    
    strF, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, \
    n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G, mass_n, \
        t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, \
        head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F, \
        v_T, v_F, alpha_C, alpha_T, alpha_F, \
        beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, \
        e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        t_T, E_C, ch_LB, ch_UB, \
        set_F_v, set_F_E = base_var
    
    ########## Inverse problem: ANN part ##########

    fv_fringe_tree, index_set = read_fringe_tree(fv_fringe_tree_filename, strF)

    descriptors, num_fv, mass_ind, max_dcp, min_dcp, avg_dcp, sd_dcp, forbidden_node, I_integer, I_nonneg = prepare_fv(
        ann_training_data_filename, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, fv_fringe_tree,
        index_set, n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G
    )

    tmp_MILP = pulp.LpProblem(name="ANN")

    std_eps = 0.01
    x_hat, x_tilde = prepare_variables_nor_std_fv(num_fv)
    tmp_MILP = add_constraints_nor_std_fv(
        tmp_MILP,
        num_fv, mass_ind, descriptors, mass_n,
        max_dcp, min_dcp, avg_dcp, sd_dcp,
        x_hat, x_tilde,
        std_eps,
        forbidden_node
    )

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

    milp_ann_vars = ann_inverter.initialize_lp_variables(
        ann, ann_a, ann_b, forbidden_node)

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

    MILP = add_constraints_ANN(MILP,
                        x_hat,
                        num_fv,
                        y, mass_ind,
                        mass_n,
                        forbidden_node)
    
    ########## add variables and constraints to base MILP ##########
    for var in MILP.variables():
        var.name = "ANN_{}_{}".format(index, var.name)
        pass

    for cons in MILP.constraints.values():
        cons.name = "ANN_{}_{}".format(index, cons.name)
        base_MILP += cons

    layer_num = len(weights)
    return y[(layer_num + 1, 1)]
