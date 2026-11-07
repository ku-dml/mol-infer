from . import ann_inverter
from ..Function.twolayered_MILP_2LMM_L import *
from ..Function.read_instance_2layer_2LMM_L import *
from ..Function.create_base_MILP import create_base_MILP


def get_value(filename, property_names):
    """
    Read minimum and maximum values for the requested properties.

    The properties are selected by their column names, so they may appear
    in any order in the value file.
    """
    num_properties = len(property_names)

    y_min = [0.0] * num_properties
    y_max = [0.0] * num_properties
    data_found = False

    with open(filename, "r") as file:
        header = [
            column.strip()
            for column in file.readline().strip().split(",")
        ]

        missing_properties = [
            property_name
            for property_name in property_names
            if property_name not in header
        ]

        if missing_properties:
            raise ValueError(
                "Properties not found in value file: "
                + ", ".join(missing_properties)
            )

        property_indices = [
            header.index(property_name)
            for property_name in property_names
        ]

        for line in file:
            values = [
                value.strip()
                for value in line.strip().split(",")
            ]

            if not values or len(values) <= max(property_indices):
                continue

            try:
                y_values = [
                    float(values[column_index])
                    for column_index in property_indices
                ]
            except ValueError:
                continue

            if not data_found:
                y_min = y_values.copy()
                y_max = y_values.copy()
                data_found = True
            else:
                for property_index in range(num_properties):
                    if y_values[property_index] < y_min[property_index]:
                        y_min[property_index] = y_values[property_index]

                    if y_values[property_index] > y_max[property_index]:
                        y_max[property_index] = y_values[property_index]

    if not data_found:
        raise ValueError(
            f"No valid property values were found in {filename}."
        )

    print("\nFinal values:")

    for property_index, property_name in enumerate(property_names):
        print(
            f"  {property_name}: "
            f"y_min = {y_min[property_index]}, "
            f"y_max = {y_max[property_index]}"
        )

    return y_min, y_max


def ANN_add_vars_constraints_to_MILP(config, index, MILP, base_var):
    """Add variables and constraints to MILP."""

    ########## preparation ##########
    prop = config.get_with_index(index, "prefix")

    target_value_lb = config.get_with_index(
        index,
        "target_value_lower_bound"
    )
    target_value_ub = config.get_with_index(
        index,
        "target_value_upper_bound"
    )

    property_names = config.input_data[index].get(
        "property_names"
    )

    if property_names is None:
        raise ValueError(
            "property_names must be provided for the ANN model "
            "in the configuration file."
        )

    if isinstance(property_names, str):
        property_names = [
            property_name.strip()
            for property_name in property_names.split(",")
            if property_name.strip()
        ]
    else:
        property_names = list(property_names)

    num_properties = len(property_names)

    if num_properties == 0:
        raise ValueError(
            "At least one ANN property must be provided."
        )

    if len(target_value_lb) != num_properties:
        raise ValueError(
            "Number of target lower bounds does not match "
            f"the number of properties: "
            f"{len(target_value_lb)} != {num_properties}."
        )

    if len(target_value_ub) != num_properties:
        raise ValueError(
            "Number of target upper bounds does not match "
            f"the number of properties: "
            f"{len(target_value_ub)} != {num_properties}."
        )

    # Input file of bias
    ann_bias_filename = f"{prop}_biases.txt"

    # Input file of weight
    ann_weights_filename = f"{prop}_weights.txt"

    # Input file of FV
    ann_training_data_filename = f"{prop}_desc.csv"
    ann_training_data_norm_filename = (
        f"{prop}_desc_norm_selected.csv"
    )

    # All fringe trees used in learning
    fv_fringe_tree_filename = f"{prop}_fringe.txt"

    # Value file
    value_filename = f"{prop}_values.txt"

    (
        strF,
        Lambda_int,
        Lambda_ex,
        Gamma_int,
        Gamma_int_less,
        Gamma_int_equal,
        Gamma_lf_ac,
        n_G,
        n_G_int,
        MASS,
        dg,
        dg_int,
        bd_int,
        na_int,
        na_ex,
        ec_int,
        fc,
        ac_lf,
        rank_G,
        mass_n,
        t_C,
        t_T,
        t_F,
        t_C_tilde,
        n_C,
        n_T,
        n_F,
        c_F,
        Lambda,
        head_C,
        tail_C,
        I_equal_one,
        I_zero_one,
        I_ge_one,
        I_ge_two,
        set_F,
        Code_F,
        v_T,
        v_F,
        alpha_C,
        alpha_T,
        alpha_F,
        beta_C,
        beta_T,
        beta_F,
        beta_CT,
        beta_TC,
        beta_CF,
        beta_TF,
        chi_T,
        chi_F,
        e_T,
        e_F,
        delta_fr_C,
        delta_fr_T,
        delta_fr_F,
        E_C,
        ch_LB,
        ch_UB,
        set_F_v,
        set_F_E,
    ) = base_var

    ########## Inverse problem: ANN part ##########

    fv_fringe_tree, index_set = read_fringe_tree(
        fv_fringe_tree_filename,
        strF
    )

    (
        descriptors,
        num_fv,
        mass_ind,
        max_dcp,
        min_dcp,
        avg_dcp,
        sd_dcp,
        forbidden_node,
        I_integer,
        I_nonneg,
        fv_list,
    ) = prepare_fv(
        ann_training_data_norm_filename,
        Lambda_int,
        Lambda_ex,
        Gamma_int,
        Gamma_int_less,
        Gamma_int_equal,
        Gamma_lf_ac,
        fv_fringe_tree,
        index_set,
        n_G,
        n_G_int,
        MASS,
        dg,
        dg_int,
        bd_int,
        na_int,
        na_ex,
        ec_int,
        fc,
        ac_lf,
        rank_G
    )

    max_dcp, min_dcp, avg_dcp, sd_dcp = prepare_max_min(
        ann_training_data_filename
    )

    std_eps = 1e-5

    x_hat, x_tilde = prepare_variables_nor_std_fv(
        num_fv,
        prop=index
    )

    MILP = add_constraints_nor_std_fv(
        MILP,
        num_fv,
        mass_ind,
        descriptors,
        fv_list,
        mass_n,
        max_dcp,
        min_dcp,
        avg_dcp,
        sd_dcp,
        x_hat,
        x_tilde,
        std_eps,
        forbidden_node,
        prop=index
    )

    training_data = ann_inverter.read_training_data(
        ann_training_data_norm_filename
    )

    des = ann_inverter.read_fv_descriptors(
        ann_training_data_norm_filename
    )
    des = des[1:]  # Drop the "CID" descriptor

    weights, biases = ann_inverter.read_trained_ANN(
        ann_weights_filename,
        ann_bias_filename
    )

    ann = ann_inverter.ANN(weights, biases)

    milp_ann_constants = ann_inverter.initialize_constants(
        ann,
        training_data
    )

    (
        ann_a,
        ann_b,
        ann_b_hat,
        ann_c,
        ann_z,
    ) = milp_ann_constants

    milp_ann_vars = ann_inverter.initialize_lp_variables(
        ann,
        ann_a,
        ann_b,
        forbidden_node,
        property_name=index
    )

    y_min, y_max = get_value(
        value_filename,
        property_names
    )

    target_value_lb_norm = []
    target_value_ub_norm = []

    for property_index in range(num_properties):
        property_range = (
            y_max[property_index]
            - y_min[property_index]
        )

        if property_range == 0:
            raise ValueError(
                "The property range is zero for "
                f"{property_names[property_index]}."
            )

        normalized_lower_bound = (
            (
                target_value_lb[property_index]
                - y_min[property_index]
            )
            / property_range
        )

        normalized_upper_bound = (
            (
                target_value_ub[property_index]
                - y_min[property_index]
            )
            / property_range
        )

        target_value_lb_norm.append(
            normalized_lower_bound
        )
        target_value_ub_norm.append(
            normalized_upper_bound
        )

    target_value_lb = target_value_lb_norm
    target_value_ub = target_value_ub_norm

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

    MILP = add_constraints_ANN(
        MILP,
        x_hat,
        num_fv,
        y,
        mass_ind,
        mass_n,
        forbidden_node,
        prop=index
    )

    layer_num = len(weights)

    num_ann_outputs = len(biases[-1])

    if num_ann_outputs != num_properties:
        raise ValueError(
            "ANN output count does not match property count: "
            f"{num_ann_outputs} ANN outputs, "
            f"{num_properties} properties."
        )

    y_outputs = [
        y[(layer_num + 1, output_index)]
        for output_index in range(1, num_ann_outputs + 1)
    ]

    return (
        y_outputs,
        y_min,
        y_max,
        x_hat,
        ann_training_data_norm_filename,
        ann
    )
