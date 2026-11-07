from Module_3.libs.RF import rf_inverter
from Module_3.libs.Function.read_instance_2layer_2LMM_L import *
from Module_3.libs.Function.twolayered_MILP_2LMM_L import *

def RF_add_vars_constraints_to_MILP(config, index, MILP, base_var):
    """Add variables and constraints to MILP"""
    ########## preparation ##########
    prop = config.get_with_index(index, "prefix")
    target_value_lb = config.get_with_index(index, "target_value_lower_bound")
    target_value_ub = config.get_with_index(index, "target_value_upper_bound")
    ########## preparation ##########

    # file for original csv
    original_dataset_filename = f"{prop}_desc.csv"
    # file for normalized csv
    normalized_dataset_filename = f"{prop}_desc_norm_selected.csv"
    # file for fringe trees
    fv_fringe_tree_filename = f"{prop}_fringe.txt"  # all fringe trees used in learning
    # value file
    value_filename = f"{prop}_values.txt"
    
    strF, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, \
        n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G, mass_n, \
        t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, \
        head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, set_F, Code_F, \
        v_T, v_F, alpha_C, alpha_T, alpha_F, \
        beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, \
        e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        E_C, ch_LB, ch_UB, \
        set_F_v, set_F_E = base_var
        
    
    fv_fringe_tree, index_set = read_fringe_tree(fv_fringe_tree_filename, strF)

    descriptors, num_fv, mass_ind, max_dcp, min_dcp, avg_dcp, sd_dcp, forbidden_node, I_integer, I_nonneg, fv_list = prepare_fv(
        normalized_dataset_filename, Lambda_int, Lambda_ex, Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_lf_ac, fv_fringe_tree,
        index_set, n_G, n_G_int, MASS, dg, dg_int, bd_int, na_int, na_ex, ec_int, fc, ac_lf, rank_G
    )

    max_dcp, min_dcp, avg_dcp, sd_dcp = prepare_max_min(original_dataset_filename)

    x_hat, x_tilde = prepare_variables_nor_std_fv(num_fv, prop=index)
    std_eps = 1e-5
    MILP = add_constraints_nor_std_fv(
        MILP,
        num_fv, mass_ind, descriptors, fv_list, mass_n,
        max_dcp, min_dcp, avg_dcp, sd_dcp,
        x_hat, x_tilde,
        std_eps,
        forbidden_node,
        prop=index
    )

    y_min, y_max = get_value(value_filename)
    target_value_lb = (target_value_lb - y_min) / (y_max - y_min)
    target_value_ub = (target_value_ub - y_min) / (y_max - y_min)
    
    ########## Inverse problem: RF part ##########
    rf_filename = f"{prop}_RF.txt"
    rf_model = rf_inverter.read_rf(rf_filename)
    rf_inv = rf_inverter.RFRegInv(n_estimators=len(rf_model.dt_list), property_name=index)
    rf_inv.build_var(rf_model, property_name=index)
    x_hat_dict = {i: x_hat[i + 1] for i in range(num_fv - 1)}
    ##### RF is not made for quad desc
    MILP = rf_inv.build_constraints(MILP, target_value_lb, target_value_ub, x_hat_dict, rf_model, property_name=index)

    y = rf_inv.ensembled_val
    return y, y_min, y_max, x_hat, fv_list, rf_model, num_fv