from Module_3.libs.ANN import ann_inverter
from Module_3.libs.Function.twolayered_MILP_2LMM_L import *
from Module_3.libs.Function.read_instance_2layer_2LMM_L import *
from Module_3.libs.Function.read_instance_2layer_2LMM_L import *
from Module_3.libs import pulp_modified as pulp

def ANN_add_vars_constraints_to_MILP(
    ann_biases_filename: str,
    ann_weights_filename: str,
    ann_training_data_filename: str,
    ann_training_data_norm_filename: str,
    fv_fringe_tree_filename: str,
    values_filename: str,
    target_value_lb: float,
    target_value_ub: float,
    MILP: pulp.LpProblem,
    base_var: tuple,
    index: int,
):
    """Add variables and constraints to MILP
    
    Arguments:
        ann_biases_filename {str} -- *_biases.txt
        ann_weights_filename {str} -- *_weights.txt
        ann_training_data_filename {str} -- *_desc.csv
        ann_training_data_norm_filename {str} -- *_desc_norm_selected.csv
        fv_fringe_tree_filename {str} -- *_fringe.txt (used in learning)
        values_filename {str} -- *_values.txt
        target_value_lb {float} -- Lower bound of target value
        target_value_ub {float} -- Upper bound of target value
    """
    strF, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, \
        n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G, mass_n, \
        t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, \
        head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F, \
        v_T, v_F, alpha_C, alpha_T, alpha_F, \
        beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, \
        e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        E_C, ch_LB, ch_UB, \
        set_F_v, set_F_E = base_var
    
    ########## Inverse problem: ANN part ##########

    fv_fringe_tree, index_set = read_fringe_tree(fv_fringe_tree_filename, strF)

    descriptors, num_fv, mass_ind, max_dcp, min_dcp, avg_dcp, sd_dcp, forbidden_node, I_integer, I_nonneg, fv_list = prepare_fv(
        ann_training_data_norm_filename, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, fv_fringe_tree,
        index_set, n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G
    )

    max_dcp, min_dcp, avg_dcp, sd_dcp = prepare_max_min(ann_training_data_filename)

    std_eps = 1e-5
    x_hat, x_tilde = prepare_variables_nor_std_fv(num_fv, prop=index)
    MILP = add_constraints_nor_std_fv(
        MILP,
        num_fv, mass_ind, descriptors, fv_list, mass_n,
        max_dcp, min_dcp, avg_dcp, sd_dcp,
        x_hat, x_tilde,
        std_eps,
        forbidden_node,
        prop=index
    )

    training_data = ann_inverter.read_training_data(
        ann_training_data_norm_filename)
    des = ann_inverter.read_fv_descriptors(ann_training_data_norm_filename)
    des = des[1:] # Drop the "CID" descriptor
    weights, biases = ann_inverter.read_trained_ANN(ann_weights_filename,
                                                    ann_biases_filename)
    ann = ann_inverter.ANN(weights, biases)
    milp_ann_constants = ann_inverter.initialize_constants(
        ann, training_data)
    ann_a, ann_b, ann_b_hat, ann_c, ann_z = milp_ann_constants

    milp_ann_vars = ann_inverter.initialize_lp_variables(
        ann, ann_a, ann_b, forbidden_node, property_name=index)

    y_min, y_max = get_value(values_filename)
    target_value_lb = (target_value_lb - y_min) / (y_max - y_min)
    target_value_ub = (target_value_ub - y_min) / (y_max - y_min)

    MILP = ann_inverter.build_MILP_ReLU(
        MILP,
        ann,
        milp_ann_vars,
        milp_ann_constants,
        target_value_lb,
        target_value_ub,
        forbidden_node,
        property_name=index
    )
    _, y, _ = milp_ann_vars

    MILP = add_constraints_ANN(MILP,
                        x_hat,
                        num_fv,
                        y, mass_ind,
                        mass_n,
                        forbidden_node,
                        prop=index)
    
    layer_num = len(weights)
    return y[(layer_num + 1, 1)], y_min, y_max, x_hat, ann_training_data_norm_filename, ann
