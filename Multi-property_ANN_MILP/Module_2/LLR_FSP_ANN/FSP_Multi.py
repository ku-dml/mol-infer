# -*- coding: utf-8 -*-
"""
Created on Sun May 18 03:03:10 2025

@author: Raveena Rai

Forward Stepwise Descriptor Selection for Multiple Properties

This script performs greedy forward stepwise descriptor selection using
5-fold cross-validation. At each step, every unselected descriptor is
evaluated together with the descriptors selected in the previous steps.

For each candidate descriptor:
    1. Train and evaluate the selected regression model using 5-fold CV.
    2. Calculate R2 separately for each target property.
    3. Combine the property-specific R2 values using equal weights.
    4. Calculate the median weighted R2 across the five folds.
    5. Select the descriptor with the highest median weighted R2.

Supported regression types:
    multi : Multiple linear regression using LinearRegression
    lasso : Lasso regression using Lasso

Usage:
    python FSP.py <lasso_csv_name> <target_file> <num_descriptors> \
<property1> [<property2> ...] <regression_type>

The Lasso CSV is read automatically from the ``log`` folder in the current
working directory.

"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold


N_SPLITS = 5
RANDOM_STATE = 0
LASSO_ALPHA = 0.0001
LASSO_MAX_ITER = 100000


def read_data_file(file_name):
    """
    Read a comma-separated or whitespace-separated data file.
    """
    data = pd.read_csv(file_name)

    if len(data.columns) == 1:
        whitespace_data = pd.read_csv(
            file_name,
            sep=r"\s+",
            engine="python",
        )

        if len(whitespace_data.columns) > 1:
            data = whitespace_data

    data.columns = data.columns.str.strip()
    return data


def create_model(regression_type, alpha=LASSO_ALPHA):
    """
    Create the selected regression model.
    """
    if regression_type == "multi":
        return LinearRegression()

    if regression_type == "lasso":
        return Lasso(
            alpha=alpha,
            max_iter=LASSO_MAX_ITER,
        )

    raise ValueError(
        "Invalid regression type. Use 'multi' or 'lasso'."
    )


def forward_stepwise_selection(
    X,
    y,
    num_descriptors,
    regression_type,
    weights,
    alpha=LASSO_ALPHA,
):
    """
    Perform greedy forward stepwise selection using 5-fold CV.
    """
    n_features = X.shape[1]

    if num_descriptors > n_features:
        raise ValueError(
            f"Requested descriptors: {num_descriptors}; "
            f"available descriptors: {n_features}."
        )

    selected_features = []

    kf = KFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=RANDOM_STATE,
    )

    for step in range(num_descriptors):
        best_feature = None
        best_median_score = -np.inf

        for feature_index in range(n_features):
            if feature_index in selected_features:
                continue

            current_features = selected_features + [feature_index]
            fold_weighted_scores = []

            for train_index, test_index in kf.split(X):
                X_train = X[train_index][:, current_features]
                X_test = X[test_index][:, current_features]
                y_train = y[train_index]
                y_test = y[test_index]

                model = create_model(
                    regression_type=regression_type,
                    alpha=alpha,
                )
                model.fit(X_train, y_train)
                y_pred = model.predict(X_test)

                if y_pred.ndim == 1:
                    y_pred = y_pred.reshape(-1, 1)

                if y_test.ndim == 1:
                    y_test = y_test.reshape(-1, 1)

                property_r2_scores = [
                    r2_score(
                        y_test[:, property_index],
                        y_pred[:, property_index],
                    )
                    for property_index in range(y.shape[1])
                ]

                weighted_r2 = sum(
                    weight * score
                    for weight, score in zip(
                        weights,
                        property_r2_scores,
                    )
                )
                fold_weighted_scores.append(weighted_r2)

            median_weighted_r2 = float(
                np.median(fold_weighted_scores)
            )

            if median_weighted_r2 > best_median_score:
                best_feature = feature_index
                best_median_score = median_weighted_r2

        if best_feature is None:
            raise RuntimeError(
                f"No descriptor could be selected at step {step + 1}."
            )

        selected_features.append(best_feature)

        print(
            f"Step {step + 1}/{num_descriptors}: "
            f"selected descriptor index {best_feature}, "
            f"median weighted CV R2 = {best_median_score:.6f}"
        )

    return selected_features


def main(args):
    """
    Read command-line inputs, load the Lasso-generated CSV automatically
    from the log folder, run FSP, and save the reduced descriptor file.
    """
    if len(args) < 7:
        print(
            "Usage:\n"
            "python FSP.py <lasso_csv_name> <target_file> "
            "<num_descriptors> <property1> "
            "[<property2> ...] <regression_type>\n\n"
            "Example:\n"
            "python FSP.py LLR_ANN_alpha_r2_desc_norm.csv "
            "norm.txt 20 alpha r2 multi"
        )
        sys.exit(1)

    # The first command-line argument is only the Lasso CSV filename.
    # It is read automatically from the existing log folder.
    lasso_file_name = args[1].strip()
    target_file = Path(args[2]).expanduser()

    try:
        num_descriptors = int(args[3])
    except ValueError as error:
        raise ValueError(
            "num_descriptors must be an integer."
        ) from error

    # The final argument is the regression type. All arguments between
    # num_descriptors and regression_type are target-property names.
    properties = [
        property_name.strip()
        for property_name in args[4:-1]
    ]
    regression_type = args[-1].strip().lower()

    # LLR_Multi creates its output inside log under the directory from
    # which the command is run. FSP reads the named CSV from there.
    log_directory = Path.cwd() / "log"
    lasso_feature_file = log_directory / Path(lasso_file_name).name

    if num_descriptors <= 0:
        raise ValueError(
            "num_descriptors must be greater than zero."
        )

    if regression_type not in {"multi", "lasso"}:
        raise ValueError(
            "The final argument must be 'multi' or 'lasso'."
        )

    if not properties:
        raise ValueError(
            "At least one target-property name is required."
        )

    if not log_directory.is_dir():
        raise FileNotFoundError(
            f"The log folder was not found: {log_directory}"
        )

    if not lasso_feature_file.is_file():
        raise FileNotFoundError(
            "Lasso feature file was not found in the log folder: "
            f"{lasso_feature_file}"
        )

    if not target_file.is_file():
        raise FileNotFoundError(
            f"Target file not found: {target_file}"
        )

    feature_df = read_data_file(lasso_feature_file)
    target_df = read_data_file(target_file)

    if "CID" not in feature_df.columns:
        raise ValueError(
            "The Lasso feature file must contain a 'CID' column."
        )

    required_target_columns = ["CID"] + properties
    missing_columns = [
        column
        for column in required_target_columns
        if column not in target_df.columns
    ]

    if missing_columns:
        raise ValueError(
            "The following columns were not found in the "
            f"target file: {missing_columns}"
        )

    target_df_selected = target_df[
        required_target_columns
    ].copy()

    target_df_selected.replace(
        ["x", "X", ""],
        pd.NA,
        inplace=True,
    )

    for property_name in properties:
        target_df_selected[property_name] = pd.to_numeric(
            target_df_selected[property_name],
            errors="coerce",
        )

    target_df_selected.dropna(
        subset=properties,
        inplace=True,
    )

    feature_df["CID"] = (
        feature_df["CID"]
        .astype(str)
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )

    target_df_selected["CID"] = (
        target_df_selected["CID"]
        .astype(str)
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )

    feature_df.drop_duplicates(
        subset="CID",
        keep="first",
        inplace=True,
    )

    target_df_selected.drop_duplicates(
        subset="CID",
        keep="first",
        inplace=True,
    )

    descriptor_columns = [
        column
        for column in feature_df.columns
        if column != "CID"
    ]

    if not descriptor_columns:
        raise ValueError(
            "No descriptor columns were found in the Lasso feature file."
        )

    for descriptor_name in descriptor_columns:
        feature_df[descriptor_name] = pd.to_numeric(
            feature_df[descriptor_name],
            errors="coerce",
        )

    merged_df = feature_df.merge(
        target_df_selected,
        on="CID",
        how="inner",
    )

    merged_df.replace(
        [np.inf, -np.inf],
        np.nan,
        inplace=True,
    )

    merged_df.dropna(
        subset=descriptor_columns + properties,
        inplace=True,
    )

    if len(merged_df) < N_SPLITS:
        raise ValueError(
            f"At least {N_SPLITS} matched molecules are required "
            f"for {N_SPLITS}-fold cross-validation."
        )

    X = merged_df[
        descriptor_columns
    ].to_numpy(dtype=float)

    y = merged_df[
        properties
    ].to_numpy(dtype=float)

    if y.ndim == 1:
        y = y.reshape(-1, 1)

    weights = [
        1.0 / len(properties)
    ] * len(properties)

    selected_features = forward_stepwise_selection(
        X=X,
        y=y,
        num_descriptors=num_descriptors,
        regression_type=regression_type,
        weights=weights,
        alpha=LASSO_ALPHA,
    )

    selected_descriptor_names = [
        descriptor_columns[index]
        for index in selected_features
    ]

    property_label = "_".join(properties)

    feature_output_file = (
        log_directory
        / f"LLR_FSP_{property_label}_"
          f"{num_descriptors}_desc_norm.csv"
    )

    merged_df[
        ["CID"] + selected_descriptor_names
    ].to_csv(
        feature_output_file,
        index=False,
    )

    print(f"Lasso feature file read from: {lasso_feature_file}")
    print(f"Saved reduced features to: {feature_output_file}")


if __name__ == "__main__":
    main(sys.argv)
