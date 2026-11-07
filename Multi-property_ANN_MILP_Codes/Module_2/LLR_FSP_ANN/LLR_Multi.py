# -*- coding: utf-8 -*-
"""
Multi-property Lasso descriptor selection.

The script runs Lasso separately for every requested target property and then
combines the selected descriptor indices using either union or intersection.

Command format:
    python LLR_ANN_Lasso_updated.py \
        <feature_csv> <target_csv> <number_of_properties> \
        <property_1> ... <property_n> <union|intersection>


"""

from __future__ import annotations

import re
import sys
import gc
import ctypes
from pathlib import Path

import numpy as np
import pandas as pd

from src import FUN, LR


# ============================================================
# MEMORY RELEASE HELPER
# ============================================================

try:
    _libc = ctypes.CDLL("libc.so.6")
except OSError:
    _libc = None


def release_unused_memory():
    """Release unused Python/Linux heap memory without changing calculations."""
    gc.collect()
    if _libc is not None:
        try:
            _libc.malloc_trim(0)
        except Exception:
            pass


# ============================================================
# LASSO SETTINGS
# ============================================================

lambda_list_lin = [
    0.000001,
    0.00001,
    0.0001,
    0.001,
    0.01,
    0.05,
    0.1,
]

lambda_list_quad = [
    0.000001,
    0.00001,
    0.0001,
    0.001,
]

STEP_lin = 20
STEP_quad = 5

RANDOM_SEED_BASE = 1000
MAX_NUM = 1000000


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def safe_filename_part(value: str) -> str:
    """Convert a property name into a safe filename component."""
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip())
    cleaned = cleaned.strip("_")
    return cleaned or "property"


# ============================================================
# FIND BEST LASSO LAMBDA
# ============================================================

def LLR_ANN_Lasso_pre(
    x,
    y,
    log_filename,
    quad=True,
    metric="R2",
    Sbar=FUN.Status_Bar(),
):
    """Run the existing Lasso pre-selection routine and return best lambda."""
    if quad:
        lambda_list = lambda_list_quad[:]
        step = STEP_quad
    else:
        lambda_list = lambda_list_lin[:]
        step = STEP_lin

    best_lambda, _, _ = LR.Lasso_pre(
        x,
        y,
        str(log_filename),
        lambda_list,
        step,
        metric=metric,
        Sbar=Sbar,
    )

    return best_lambda


# ============================================================
# SELECT NONZERO LASSO FEATURES
# ============================================================

def shrink_dataset(
    x,
    y,
    best_lambda,
    fv_list,
):
    """Select descriptors having nonzero Lasso coefficients."""
    nonzero_linear = 0
    nonzero_quad = 0

    _, arr = LR.learn_Lasso_eval(
        x,
        y,
        best_lambda,
    )

    # Preserve the original fallback behaviour.
    if len(arr) == 0:
        print(
            "Warning: No nonzero descriptor was selected. "
            "All descriptors will be retained."
        )
        arr = list(range(x.shape[1]))

    arr = [int(index) for index in arr]
    xsel = x[:, arr]

    for index in arr:
        if "*" in str(fv_list[index]):
            nonzero_quad += 1
        else:
            nonzero_linear += 1

    return (
        xsel,
        arr,
        nonzero_linear,
        nonzero_quad,
    )


# ============================================================
# MAIN FUNCTION
# ============================================================

def main(args):
    feature_file = args[1]
    target_file = args[2]
    number_of_properties = int(args[3])
    properties = args[4:4 + number_of_properties]
    selection_type = args[4 + number_of_properties].lower()

    # The output directory is based on where the command is executed,
    # not on the location of the Python script or input files.
    log_directory = Path.cwd() / "log"
    log_directory.mkdir(parents=True, exist_ok=True)

    property_filename_parts = [
        safe_filename_part(prop)
        for prop in properties
    ]
    property_string = "_".join(property_filename_parts)

    print("=" * 70)
    print("MULTI-PROPERTY LASSO FEATURE SELECTION")
    print("=" * 70)
    print("Feature file:", feature_file)
    print("Target file:", target_file)
    print("Number of properties:", number_of_properties)
    print("Properties:", properties)
    print("Feature combination type:", selection_type)
    print("=" * 70)

    # --------------------------------------------------------
    # LOAD FEATURE AND TARGET DATA
    # --------------------------------------------------------

    feature_df = pd.read_csv(feature_file)
    target_df = pd.read_csv(target_file)

    feature_df.columns = feature_df.columns.astype(str).str.strip()
    target_df.columns = target_df.columns.astype(str).str.strip()

    if "CID" not in feature_df.columns:
        raise ValueError("'CID' column was not found in the feature file.")

    if "CID" not in target_df.columns:
        raise ValueError("'CID' column was not found in the target file.")

    missing_properties = [
        prop
        for prop in properties
        if prop not in target_df.columns
    ]

    if missing_properties:
        raise ValueError(
            "The following requested properties were not found in the target "
            f"file: {', '.join(missing_properties)}"
        )

    # --------------------------------------------------------
    # SELECT REQUESTED TARGET PROPERTIES
    # --------------------------------------------------------

    selected_columns = ["CID"] + properties
    target_df1 = target_df[selected_columns].copy()

    target_df1.replace(
        ["x", "X", ""],
        pd.NA,
        inplace=True,
    )

    for prop in properties:
        target_df1[prop] = pd.to_numeric(
            target_df1[prop],
            errors="coerce",
        )

    target_df1 = target_df1.dropna(subset=properties)

    # --------------------------------------------------------
    # MATCH FEATURE AND TARGET FILES USING CID
    # --------------------------------------------------------

    feature_df = feature_df.copy()

    feature_df["CID"] = (
        feature_df["CID"]
        .astype(str)
        .str.strip()
    )

    target_df1["CID"] = (
        target_df1["CID"]
        .astype(str)
        .str.strip()
    )

    feature_df_filtered = (
        feature_df
        .drop_duplicates(subset="CID", keep="first")
    )

    target_df1_filtered = (
        target_df1
        .drop_duplicates(subset="CID", keep="first")
    )

    merged_df = (
        feature_df_filtered
        .merge(
            target_df1_filtered,
            on="CID",
            how="inner",
        )
        .sort_values("CID")
        .reset_index(drop=True)
    )

    if merged_df.empty:
        raise ValueError(
            "No matching CID values were found between the feature and "
            "target files."
        )

    print("\nLength of merged DataFrame:", len(merged_df))

    # --------------------------------------------------------
    # EXTRACT DESCRIPTOR MATRIX
    # --------------------------------------------------------

    descriptor_columns = [
        column
        for column in feature_df_filtered.columns
        if column != "CID"
    ]

    if not descriptor_columns:
        raise ValueError("The feature file contains no descriptor columns.")

    # Make descriptor conversion explicit so invalid cells are reported.
    descriptor_data = merged_df[descriptor_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )

    invalid_descriptor_columns = descriptor_data.columns[
        descriptor_data.isna().any()
    ].tolist()

    if invalid_descriptor_columns:
        raise ValueError(
            "The following descriptor columns contain missing or non-numeric "
            f"values after CID matching: {', '.join(invalid_descriptor_columns)}"
        )

    X = descriptor_data.to_numpy(dtype=float)
    fv_list = descriptor_columns

    # Preserve only the values needed later.
    cid_column = merged_df["CID"].reset_index(drop=True).copy()
    target_arrays = {
        prop: merged_df[prop].to_numpy(dtype=float, copy=True)
        for prop in properties
    }

    # Release duplicate pandas objects. Numerical values in X/targets are unchanged.
    del descriptor_data
    del feature_df
    del target_df
    del target_df1
    del feature_df_filtered
    del target_df1_filtered
    del merged_df
    release_unused_memory()

    print("Original number of descriptors:", X.shape[1])

    # --------------------------------------------------------
    # RUN LASSO SEPARATELY FOR EACH PROPERTY
    # --------------------------------------------------------

    selected_feature_indices: set[int] | None = None

    for property_number, prop in enumerate(properties, start=1):
        print("\n" + "-" * 70)
        print(
            f"Property {property_number}/{number_of_properties}: {prop}"
        )
        print("-" * 70)

        y = target_arrays[prop]

        property_filename = safe_filename_part(prop)
        log_filename = (
            log_directory
            / f"LLR_{property_filename}_lasso_log.txt"
        )

        # Remove an old log with the same name so a new run does not append
        # results to data generated by a previous run.
        if log_filename.exists():
            log_filename.unlink()

        best_lambda = LLR_ANN_Lasso_pre(
            X,
            y,
            log_filename,
            quad=False,
            metric="R2",
        )

        print(f"Best lambda for {prop}: {best_lambda}")

        (
            selected_x_temp,
            selected_indices,
            nonzero_linear,
            nonzero_quad,
        ) = shrink_dataset(
            X,
            y,
            best_lambda,
            fv_list,
        )

        current_feature_set = set(selected_indices)

        if selected_feature_indices is None:
            selected_feature_indices = current_feature_set
        elif selection_type == "union":
            selected_feature_indices.update(current_feature_set)
        else:
            selected_feature_indices.intersection_update(
                current_feature_set
            )

        print(
            #f"Selected features for {prop} "
            f"(count {len(selected_indices)}):"
        )
        #print(selected_indices)
        print("Number of linear descriptors:", nonzero_linear)
        #print("Number of quadratic descriptors:", nonzero_quad)
        #print("Lasso log saved as:", log_filename)

        # These temporary objects are no longer needed after updating the
        # union/intersection feature set.
        del selected_x_temp
        del selected_indices
        del current_feature_set
        del y
        release_unused_memory()

    if selected_feature_indices is None:
        raise RuntimeError("Lasso did not return a descriptor set.")

    selected_feature_indices_sorted = sorted(selected_feature_indices)

    # --------------------------------------------------------
    # CREATE FINAL SELECTED FEATURE SET
    # --------------------------------------------------------

    print("\n" + "=" * 70)
    print(
        f"FINAL {selection_type.upper()} OF SELECTED DESCRIPTORS"
    )
    print("=" * 70)
    print(
        f"Total {selection_type} descriptors:",
        len(selected_feature_indices_sorted),
    )
    print(
        f"{selection_type.capitalize()} descriptor indices:"
    )
    print(selected_feature_indices_sorted)

    selected_feature_names = [
        descriptor_columns[index]
        for index in selected_feature_indices_sorted
    ]

    #print(f"\nSelected {selection_type} descriptor names:")
    #for descriptor_name in selected_feature_names:
        #print(descriptor_name)

    # --------------------------------------------------------
    # CREATE THE ONE COMBINED SELECTED-DESCRIPTOR FILE
    # --------------------------------------------------------

    selected_features_df = pd.DataFrame(
        X[:, selected_feature_indices_sorted],
        columns=selected_feature_names,
    )

    final_df = pd.concat(
        [cid_column, selected_features_df],
        axis=1,
    )

    output_filename = (
        log_directory
        / f"LLR_{property_string}_desc_norm.csv"
    )

    final_df.to_csv(
        output_filename,
        index=False,
    )

    # print(
        # f"\nSaved selected {selection_type} descriptors to:"
    # )
    # print(output_filename)
    print(f"Lasso-selected descriptor file: {output_filename.name}")
    

    # Release remaining large arrays after the output file is written.
    del selected_features_df
    del final_df
    del target_arrays
    del X
    release_unused_memory()

    print("\nProcessing completed successfully.")
    


# ============================================================
# RUN SCRIPT
# ============================================================

if __name__ == "__main__":
    try:
        main(sys.argv)
    except Exception as error:
        print(f"\nError: {error}", file=sys.stderr)
        sys.exit(1)
