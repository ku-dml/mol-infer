# -*- coding: utf-8 -*-
"""
Generalized Optuna-based ANN architecture search + 10x5 CV.
Supports any number of target properties (minimum 1) and a list of
descriptor counts applied sequentially to the same feature file.

The original experimental settings are preserved while the output dimension is generalized:
- Global seed = 0
- TensorFlow/Keras deterministic settings
- TPESampler(seed=0)
- No external train/test split; CV runs on the complete matched dataset
- KFold(n_splits=5, shuffle=True, random_state=1)
- 10 Optuna trials
- 10 x 5-fold final CV
- ReLU hidden layers and linear output
- Adam learning rate = 0.001
- MSE loss
- rho = 0.98 and maxitr = 10000
- Equal-weight weighted-sum Optuna score for any number of targets
- Combined statistics = average of per-target median values (Eq. (3))
"""

import os
import sys
import random
import warnings
import gc
import shutil
import numpy as np

# ===== GLOBAL REPRODUCIBILITY SETTINGS =====
SEED = 0
os.environ["PYTHONHASHSEED"] = str(SEED)
os.environ["TF_DETERMINISTIC_OPS"] = "1"
os.environ["TF_CUDNN_DETERMINISTIC"] = "1"
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

warnings.filterwarnings("ignore", category=UserWarning)

random.seed(SEED)
np.random.seed(SEED)

import optuna
import tensorflow as tf
import keras
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam

# TensorFlow / Keras seeds
# These are intentionally kept identical to the original two-property code.
tf.random.set_seed(SEED)
keras.utils.set_random_seed(SEED)
tf.get_logger().setLevel("FATAL")
optuna.logging.set_verbosity(optuna.logging.WARNING)

# ===== OUTPUT DIRECTORIES =====
# The ANN_Result folder and all four subfolders must already exist.
# This script does not create any folder automatically.
BASE_OUTPUT_DIR = os.path.join(os.getcwd(), "ANN_Result")
ARCH_DIR = os.path.join(BASE_OUTPUT_DIR, "ANN_pre_result")
CV_DIR = os.path.join(BASE_OUTPUT_DIR, "ANN_Eval")
RESULT_DIR = os.path.join(BASE_OUTPUT_DIR, "Result")
WB_DIR = os.path.join(BASE_OUTPUT_DIR, "weight_biases")


# ==========================================================
#                    SINGLE PROGRESS BAR
# ==========================================================
def show_progress(completed_steps, total_steps):
    """Display one progress bar on the same terminal line."""
    bar_length = 50
    progress = completed_steps / total_steps
    filled_length = int(bar_length * progress)
    bar = "=" * filled_length + " " * (bar_length - filled_length)

    print(
        f"\r[{bar}] {progress * 100:6.2f}%",
        end="",
        flush=True,
    )

    if completed_steps >= total_steps:
        print()


# ==========================================================
#                       MODEL
# ==========================================================
def build_ann(architecture, input_shape, num_targets, alpha=1e-5):
    """
    architecture: list of hidden units, e.g. [101, 111, ...]
    input_shape: (n_features,)
    num_targets: number of output properties (minimum 1)

    The unused alpha argument is retained to preserve the original function style.
    """
    model = Sequential()
    model.add(
        Dense(
            units=architecture[0],
            activation="relu",
            input_shape=input_shape,
        )
    )
    for units in architecture[1:]:
        model.add(Dense(units=units, activation="relu"))

    model.add(Dense(units=num_targets, activation="linear"))
    return model


def early_stop_reached(r2_scores_train, rho):
    """
    Preserve the exact original two-target stopping expression when two
    properties are supplied; use the same condition across all targets for
    all multi-property tasks.
    """
    if len(r2_scores_train) == 2:
        return (
            r2_scores_train[0] >= rho
            and r2_scores_train[1] >= rho
        )

    return all(score >= rho for score in r2_scores_train)


def calculate_weighted_score(median_test_scores):
    """
    Equal-weight weighted-sum formula.

    For two targets, this deliberately uses the exact original expression:
        0.5 * score_1 + 0.5 * score_2

    For three or more targets, the same weighted-sum formula is extended using weights
    that sum to one. This is not implemented with np.mean().
    """
    num_targets = len(median_test_scores)

    if num_targets == 1:
        return float(median_test_scores[0])

    if num_targets == 2:
        w1 = 0.5
        w2 = 0.5
        return (
            w1 * median_test_scores[0]
            + w2 * median_test_scores[1]
        )

    weights = np.full(num_targets, 1.0 / num_targets, dtype=float)
    weighted_score = sum(
        weights[i] * median_test_scores[i]
        for i in range(num_targets)
    )
    return float(weighted_score)


# ==========================================================
#                     OPTUNA OBJECTIVE
# ==========================================================
def objective(
    trial,
    X_all,
    y_norm_all,
    kf,
    properties,
    progress_status,
    total_progress_steps,
    rho=0.98,
    maxitr=10000,
    batch_size=None,
    base_tag="",
):
    """
    Optuna objective:
      - uses normalized targets
      - performs 5-fold CV
      - stops a fold when every target reaches R² >= rho
      - returns the equal-weight weighted median testing R²
    """
    if batch_size is None:
        batch_size = len(X_all)

    num_targets = y_norm_all.shape[1]

    # Original Optuna search space and parameter names
    n_layers = trial.suggest_int("n_layers", 1, 8)
    architecture = []
    for i in range(n_layers):
        n_neurons = trial.suggest_int(f"n_units_layer_{i}", 16, 125)
        architecture.append(n_neurons)

    fold_train_r2_scores = []
    fold_test_r2_scores = []

    for train_index, test_index in kf.split(X_all):
        X_train_fold, X_test_fold = X_all[train_index], X_all[test_index]
        y_train_fold_norm = y_norm_all[train_index]
        y_test_fold_norm = y_norm_all[test_index]

        # Release Keras/TensorFlow state from the previous fold before
        # constructing a new model. This prevents model/graph accumulation.
        tf.keras.backend.clear_session()

        model = build_ann(
            architecture,
            input_shape=(X_all.shape[1],),
            num_targets=num_targets,
        )
        model.compile(
            loss="mean_squared_error",
            optimizer=Adam(learning_rate=0.001),
        )

        for t in range(10, maxitr + 1, 10):
            model.fit(
                X_train_fold,
                y_train_fold_norm,
                epochs=1,
                batch_size=batch_size,
                verbose=0,
            )

            y_train_pred_norm = model(X_train_fold, training=False).numpy()
            r2_scores_train = [
                r2_score(
                    y_train_fold_norm[:, i],
                    y_train_pred_norm[:, i],
                )
                for i in range(num_targets)
            ]

            if early_stop_reached(r2_scores_train, rho):
                break

        y_test_pred_norm = model(X_test_fold, training=False).numpy()
        r2_scores_test = [
            r2_score(
                y_test_fold_norm[:, i],
                y_test_pred_norm[:, i],
            )
            for i in range(num_targets)
        ]

        fold_train_r2_scores.append(r2_scores_train)
        fold_test_r2_scores.append(r2_scores_test)

        progress_status["completed"] += 1
        show_progress(
            progress_status["completed"],
            total_progress_steps,
        )

        # Explicitly release fold-specific TensorFlow objects.
        del model
        del y_train_pred_norm
        del y_test_pred_norm
        tf.keras.backend.clear_session()
        gc.collect()

    median_train_scores = np.median(
        np.asarray(fold_train_r2_scores),
        axis=0,
    )
    median_test_scores = np.median(
        np.asarray(fold_test_r2_scores),
        axis=0,
    )

    weighted_median_test_score = calculate_weighted_score(
        median_test_scores
    )

    # ===== ARCHITECTURE LOG FILE =====
    if base_tag:
        arch_filename = os.path.join(
            ARCH_DIR,
            f"architecture_result-MAE_{base_tag}.txt",
        )
    else:
        arch_filename = os.path.join(
            ARCH_DIR,
            "architecture_result-MAE_default.txt",
        )

    with open(arch_filename, "a", encoding="utf-8") as file:
        file.write(f"Architecture: {architecture}\n")
        for i, property_name in enumerate(properties):
            file.write(
                f"  Median Training R² ({property_name}): "
                f"{median_train_scores[i]:.6f}\n"
            )
        for i, property_name in enumerate(properties):
            file.write(
                f"  Median Testing R² ({property_name}): "
                f"{median_test_scores[i]:.6f}\n"
            )
        file.write(
            "  Weighted Median Testing R²: "
            f"{weighted_median_test_score:.6f}\n\n"
        )

    return weighted_median_test_score


# ==========================================================
#            MAIN 10x5 CV WITH MAE/NMAE
# ==========================================================
def learn_ANN_eval(
    X_all,
    y_norm_all,
    y_orig_all,
    properties,
    num_descriptors,
    base_tag,
    y_min,
    y_range,
    rho=0.98,
    maxitr=10000,
    batch_size=None,
):
    """
    Runs Optuna architecture selection followed by 10 x 5-fold CV
    on the complete matched dataset. No external train/test split is used.
    """
    if batch_size is None:
        batch_size = len(X_all)

    num_targets = y_norm_all.shape[1]
    kf = KFold(n_splits=5, shuffle=True, random_state=1)

    # 10 Optuna trials x 5 folds + 10 repetitions x 5 folds = 100
    total_progress_steps = 100
    progress_status = {"completed": 0}
    show_progress(0, total_progress_steps)

    study = optuna.create_study(
        sampler=optuna.samplers.TPESampler(seed=0),
        direction="maximize",
    )
    study.optimize(
        lambda trial: objective(
            trial,
            X_all,
            y_norm_all,
            kf,
            properties,
            progress_status,
            total_progress_steps,
            rho,
            maxitr,
            batch_size,
            base_tag,
        ),
        n_trials=10,
    )

    best_trial = study.best_trial
    best_architecture = [
        best_trial.params[f"n_units_layer_{i}"]
        for i in range(best_trial.params["n_layers"])
    ]

    final_train_r2_scores = []
    final_test_r2_scores = []
    cv_fold_records = []

    fold_train_mae = []
    fold_test_mae = []
    fold_train_nmae = []
    fold_test_nmae = []

    for iteration in range(10):
        for fold, (train_index, test_index) in enumerate(kf.split(X_all)):
            X_tr_fold = X_all[train_index]
            X_test_fold = X_all[test_index]
            y_tr_fold_norm = y_norm_all[train_index]
            y_test_fold_norm = y_norm_all[test_index]
            y_tr_fold_orig = y_orig_all[train_index]
            y_test_fold_orig = y_orig_all[test_index]

            # Release Keras/TensorFlow state from the previous fold before
            # constructing a new model. This prevents model/graph accumulation.
            tf.keras.backend.clear_session()

            model = build_ann(
                best_architecture,
                input_shape=(X_all.shape[1],),
                num_targets=num_targets,
            )
            model.compile(
                loss="mean_squared_error",
                optimizer=Adam(learning_rate=0.001),
            )

            for t in range(10, maxitr + 1, 10):
                model.fit(
                    X_tr_fold,
                    y_tr_fold_norm,
                    epochs=1,
                    batch_size=batch_size,
                    verbose=0,
                )

                y_tr_pred_norm = model(X_tr_fold, training=False).numpy()
                train_r2 = [
                    r2_score(
                        y_tr_fold_norm[:, d],
                        y_tr_pred_norm[:, d],
                    )
                    for d in range(num_targets)
                ]

                if early_stop_reached(train_r2, rho):
                    break

            y_test_pred_norm = model(X_test_fold, training=False).numpy()
            test_r2 = [
                r2_score(
                    y_test_fold_norm[:, d],
                    y_test_pred_norm[:, d],
                )
                for d in range(num_targets)
            ]

            final_train_r2_scores.append(train_r2)
            final_test_r2_scores.append(test_r2)
            # Combined R² for this fold: average across all target properties.
            combined_train_r2 = float(np.mean(np.asarray(train_r2)))
            combined_test_r2 = float(np.mean(np.asarray(test_r2)))

            cv_fold_records.append(
                {
                    "iteration": iteration + 1,
                    "fold": fold + 1,
                    "train_r2": train_r2,
                    "test_r2": test_r2,
                    "combined_train_r2": combined_train_r2,
                    "combined_test_r2": combined_test_r2,
                }
            )

            # Back-transform predictions and compute MAE/NMAE
            y_tr_pred_orig = y_tr_pred_norm * y_range + y_min
            y_test_pred_orig = y_test_pred_norm * y_range + y_min

            mae_train = np.mean(
                np.abs(y_tr_pred_orig - y_tr_fold_orig),
                axis=0,
            )
            mae_test = np.mean(
                np.abs(y_test_pred_orig - y_test_fold_orig),
                axis=0,
            )

            nmae_train = mae_train / y_range
            nmae_test = mae_test / y_range

            fold_train_mae.append(mae_train)
            fold_test_mae.append(mae_test)
            fold_train_nmae.append(nmae_train)
            fold_test_nmae.append(nmae_test)

            write_weights_biases(
                model,
                iteration + 1,
                fold + 1,
                properties,
                num_descriptors,
                base_tag,
            )

            progress_status["completed"] += 1
            show_progress(
                progress_status["completed"],
                total_progress_steps,
            )

            # Explicitly release fold-specific TensorFlow objects and
            # prediction arrays before the next model is created.
            del model
            del y_tr_pred_norm
            del y_test_pred_norm
            del y_tr_pred_orig
            del y_test_pred_orig
            tf.keras.backend.clear_session()
            gc.collect()

    # ---------- R² statistics ----------
    final_train_r2_arr = np.asarray(final_train_r2_scores)
    final_test_r2_arr = np.asarray(final_test_r2_scores)

    median_train_r2_per_dimension = np.median(
        final_train_r2_arr,
        axis=0,
    )
    median_test_r2_per_dimension = np.median(
        final_test_r2_arr,
        axis=0,
    )

    # Combined R² (Eq. (3)): average of the per-target median R² values
    median_train_r2_combined = np.mean(median_train_r2_per_dimension)
    median_test_r2_combined = np.mean(median_test_r2_per_dimension)

    # ---------- MEDIAN-CLOSEST REAL FOLD (same L1 technique as previous code) ----------
    # No retraining and no change to any R2/MAE/NMAE calculation.
    # All 50 original final-CV weights/biases remain saved in WB_DIR.
    # This only identifies the real fold closest to the per-target median
    # Testing R2 vector and copies that fold's existing files to RESULT_DIR.
    median_testing_vector = np.asarray(
        median_test_r2_per_dimension,
        dtype=float,
    )

    median_closest_idx = None
    median_closest_l1 = float("inf")

    for idx, r2_vec in enumerate(final_test_r2_scores):
        r2_vec = np.asarray(r2_vec, dtype=float)
        l1_score = float(
            np.sum(np.abs(r2_vec - median_testing_vector))
        )
        if l1_score < median_closest_l1:
            median_closest_l1 = l1_score
            median_closest_idx = idx

    median_closest_iteration = median_closest_idx // 5 + 1
    median_closest_fold = median_closest_idx % 5 + 1
    median_closest_test_r2 = final_test_r2_scores[median_closest_idx]

    source_weights_file = os.path.join(
        WB_DIR,
        f"qm9_{median_closest_iteration}_{median_closest_fold}_"
        f"{base_tag}_weights.txt",
    )
    source_biases_file = os.path.join(
        WB_DIR,
        f"qm9_{median_closest_iteration}_{median_closest_fold}_"
        f"{base_tag}_biases.txt",
    )

    median_weights_file = os.path.join(
        RESULT_DIR,
        f"qm9_median_closest_{base_tag}_weights.txt",
    )
    median_biases_file = os.path.join(
        RESULT_DIR,
        f"qm9_median_closest_{base_tag}_biases.txt",
    )

    shutil.copyfile(source_weights_file, median_weights_file)
    shutil.copyfile(source_biases_file, median_biases_file)

    median_info_file = os.path.join(
        RESULT_DIR,
        f"qm9_median_closest_{base_tag}_info.txt",
    )
    with open(median_info_file, "w", encoding="utf-8") as file:
        file.write(
            f"Median Testing R2 per target: "
            f"{median_test_r2_per_dimension}\n"
        )
        file.write(
            f"Selected Iteration: {median_closest_iteration}\n"
        )
        file.write(
            f"Selected Fold: {median_closest_fold}\n"
        )
        file.write(
            f"Selected Fold Testing R2 per target: "
            f"{median_closest_test_r2}\n"
        )
        file.write(
            f"L1 closeness score: {median_closest_l1}\n"
        )
        file.write(
            f"Source weights: {source_weights_file}\n"
        )
        file.write(
            f"Source biases: {source_biases_file}\n"
        )
        file.write(
            f"Copied median-closest weights: {median_weights_file}\n"
        )
        file.write(
            f"Copied median-closest biases: {median_biases_file}\n"
        )

    print("\n===== MEDIAN-CLOSEST FINAL-CV FOLD (L1) =====")
    print(
        "Median Testing R2 per target: "
        f"{median_test_r2_per_dimension}"
    )
    print(
        f"Selected Iteration, Fold: "
        f"{median_closest_iteration}, {median_closest_fold}"
    )
    print(
        "Selected Fold Testing R2 per target: "
        f"{median_closest_test_r2}"
    )
    print(f"L1 closeness score: {median_closest_l1}")
    print(f"Median-closest weights saved to: {median_weights_file}")
    print(f"Median-closest biases saved to: {median_biases_file}")

    # ---------- MAE / NMAE statistics ----------
    train_mae_arr = np.vstack(fold_train_mae)
    test_mae_arr = np.vstack(fold_test_mae)
    train_nmae_arr = np.vstack(fold_train_nmae)
    test_nmae_arr = np.vstack(fold_test_nmae)

    median_train_mae = np.median(train_mae_arr, axis=0)
    median_test_mae = np.median(test_mae_arr, axis=0)
    median_train_nmae = np.median(train_nmae_arr, axis=0)
    median_test_nmae = np.median(test_nmae_arr, axis=0)

    # Combined MAE/NMAE: average of the per-target median values
    median_train_mae_combined = np.mean(median_train_mae)
    median_test_mae_combined = np.mean(median_test_mae)
    median_train_nmae_combined = np.mean(median_train_nmae)
    median_test_nmae_combined = np.mean(median_test_nmae)

    print(
        f"Best architecture: {best_architecture} "
        f"with Median Testing R2 score: {best_trial.value}"
    )
    print(
        "\n================ FINAL MEDIAN RESULTS "
        "(CV on full dataset) ================"
    )
    print(f"Targets: {', '.join(properties)}")
    print(
        "Median Train R2 per target (norm): "
        f"{median_train_r2_per_dimension}"
    )
    print(
        "Median Testing R2 per target (norm): "
        f"{median_test_r2_per_dimension}"
    )
    print(
        "Median combined Train R2 (norm, average): "
        f"{median_train_r2_combined}"
    )
    print(
        "Median combined Testing R2 (norm, average): "
        f"{median_test_r2_combined}"
    )
    print("------------------------------------------------------")
    print(f"Median Train MAE  (orig units): {median_train_mae}")
    print(f"Median Testing MAE  (orig units): {median_test_mae}")
    print(
        "Median combined Train MAE (orig units): "
        f"{median_train_mae_combined}"
    )
    print(
        "Median combined Testing MAE (orig units): "
        f"{median_test_mae_combined}"
    )
    print(f"Median Train NMAE: {median_train_nmae}")
    print(f"Median Testing NMAE: {median_test_nmae}")
    print(
        "Median combined Train NMAE: "
        f"{median_train_nmae_combined}"
    )
    print(
        "Median combined Testing NMAE: "
        f"{median_test_nmae_combined}"
    )
    print("=====================================================================\n")

    cv_filename = os.path.join(
        CV_DIR,
        f"cross_validation_{base_tag}.txt",
    )
    with open(cv_filename, "w", encoding="utf-8") as file:
        file.write(f"Targets: {', '.join(properties)}\n\n")

        # Detailed 10 x 5-fold R² values belong in ANN_Eval.
        for record in cv_fold_records:
            file.write(
                f"Iteration {record['iteration']}, "
                f"Fold {record['fold']}, "
                f"Training R2: {record['train_r2']}, "
                f"Testing R2: {record['test_r2']}, "
                f"Combined Training R2: {record['combined_train_r2']}, "
                f"Combined Testing R2: {record['combined_test_r2']}\n"
            )

        file.write("\n")
        file.write(
            "Median Training R2 score: "
            f"{median_train_r2_per_dimension}\n"
        )
        file.write(
            "Median Testing R2 score: "
            f"{median_test_r2_per_dimension}\n"
        )
        file.write(
            "Combined Median Training R2 score: "
            f"{median_train_r2_combined}\n"
        )
        file.write(
            "Combined Median Testing R2 score: "
            f"{median_test_r2_combined}\n"
        )

    return {
        "best_architecture": best_architecture,
        "median_train_r2_per_dim": median_train_r2_per_dimension,
        "median_test_r2_per_dim": median_test_r2_per_dimension,
        "median_train_r2_combined": median_train_r2_combined,
        "median_test_r2_combined": median_test_r2_combined,
        "median_train_mae": median_train_mae,
        "median_test_mae": median_test_mae,
        "median_train_mae_combined": median_train_mae_combined,
        "median_test_mae_combined": median_test_mae_combined,
        "median_train_nmae": median_train_nmae,
        "median_test_nmae": median_test_nmae,
        "median_train_nmae_combined": median_train_nmae_combined,
        "median_test_nmae_combined": median_test_nmae_combined,
        "median_closest_iteration": median_closest_iteration,
        "median_closest_fold": median_closest_fold,
        "median_closest_test_r2": median_closest_test_r2,
        "median_closest_l1": median_closest_l1,
        "median_closest_weights_file": median_weights_file,
        "median_closest_biases_file": median_biases_file,
        "median_closest_info_file": median_info_file,
    }


# ==========================================================
#                 WEIGHTS / BIASES WRITER
# ==========================================================
def write_weights_biases(
    model,
    iteration,
    fold,
    properties,
    num_descriptors,
    base_tag,
):
    """Save each trained model's weights and biases."""
    if base_tag:
        tag = base_tag
    else:
        tag = f"{'_'.join(properties)}_{num_descriptors}desc"

    weights_filename = os.path.join(
        WB_DIR,
        f"qm9_{iteration}_{fold}_{tag}_weights.txt",
    )
    biases_filename = os.path.join(
        WB_DIR,
        f"qm9_{iteration}_{fold}_{tag}_biases.txt",
    )

    weights_and_biases = model.get_weights()
    weights = weights_and_biases[::2]
    biases = weights_and_biases[1::2]
    num_features = weights[0].shape[0]
    num_targets = len(properties)

    with open(weights_filename, "w", encoding="utf-8") as file:
        file.write(str(num_features) + " ")
        for layer in model.layers:
            if (
                hasattr(layer, "units")
                and layer.activation.__name__ == "relu"
            ):
                file.write(str(layer.units) + " ")

        # Dynamic output-neuron count: equal to the number of properties
        file.write(f"{num_targets}\n")

        for item in weights:
            for i in range(item.shape[0]):
                for j in range(item.shape[1]):
                    if abs(item[i][j]) > 10 ** (-308):
                        file.write(str(item[i][j]) + " ")
                    else:
                        file.write("0 ")
                file.write("\n")

    with open(biases_filename, "w", encoding="utf-8") as file:
        for item in biases:
            for i in range(item.shape[0]):
                file.write(str(item[i]) + "\n")


# ==========================================================
#                       MAIN
# ==========================================================
def split_descriptor_counts_and_properties(arguments):
    """
    Split command-line values into:
      1. one or more consecutive integer descriptor counts;
      2. one or more target-property names.

    Example:
        50 45 35 25 20 15 alpha r2

    becomes:
        descriptor_counts = [50, 45, 35, 25, 20, 15]
        properties = ["alpha", "r2"]
    """
    descriptor_counts = []
    property_start_index = None

    for index, value in enumerate(arguments):
        try:
            descriptor_count = int(value)
        except ValueError:
            property_start_index = index
            break

        if descriptor_count < 1:
            raise ValueError(
                "Every descriptor count must be a positive integer."
            )

        descriptor_counts.append(descriptor_count)

    if not descriptor_counts:
        raise ValueError(
            "Provide at least one descriptor count before the property names."
        )

    if property_start_index is None:
        raise ValueError(
            "Provide at least one target property after the descriptor counts."
        )

    properties = arguments[property_start_index:]

    # Once property names begin, another integer is ambiguous and is therefore
    # rejected instead of silently treating it as a property.
    numeric_properties = [
        value for value in properties
        if value.lstrip("+-").isdigit()
    ]
    if numeric_properties:
        raise ValueError(
            "All descriptor counts must appear consecutively before the "
            "property names. Unexpected numeric value(s) after a property: "
            + ", ".join(numeric_properties)
        )

    # Remove repeated descriptor counts while preserving the supplied order.
    unique_descriptor_counts = []
    seen_counts = set()

    for count in descriptor_counts:
        if count not in seen_counts:
            seen_counts.add(count)
            unique_descriptor_counts.append(count)

    return unique_descriptor_counts, properties


def reset_reproducibility_state():
    """
    Reset random states before each descriptor-count run.

    This makes a run inside the descriptor list behave like a separate
    execution of the original single-descriptor-count command.
    """
    tf.keras.backend.clear_session()
    random.seed(SEED)
    np.random.seed(SEED)
    tf.random.set_seed(SEED)
    keras.utils.set_random_seed(SEED)
    gc.collect()


def main(args):
    """
    Command format:
      python ANN_general_multi_descriptor_list.py \
          feature.csv normalized.txt original.txt \
          descriptor_count1 [descriptor_count2 ...] \
          property1 [property2 ...]

    The same feature file is used for every descriptor count. For each count,
    the script selects the first N descriptor columns after CID, runs the full
    Optuna + 10x5 CV workflow, and saves results with N in the output filename.
    """
    if len(args) < 6:
        print(
            "Usage:\n"
            "python ANN_general_multi_descriptor_list.py "
            "desc.csv norm.txt unnormalized.txt "
            "descriptor_count1 [descriptor_count2 ...] "
            "property1 [property2 ...]\n\n"
            "Examples:\n"
            "Single descriptor count (old format remains valid):\n"
            "python ANN_general_multi_descriptor_list.py "
            "desc.csv norm.txt original.txt 20 alpha r2\n\n"
            "Descriptor-count list:\n"
            "python ANN_general_multi_descriptor_list.py "
            "desc.csv norm.txt original.txt "
            "50 45 35 25 20 15 alpha r2\n\n"
            "Multiple target properties:\n"
            "python ANN_general_multi_descriptor_list.py "
            "desc.csv norm.txt original.txt "
            "50 45 35 25 20 15 alpha r2 zpve cv gap lumo"
        )
        sys.exit(1)

    feature_file = args[1]
    norm_target_file = args[2]
    unscaled_target_file = args[3]

    descriptor_counts, properties = (
        split_descriptor_counts_and_properties(args[4:])
    )

    num_targets = len(properties)

    if num_targets < 1:
        raise ValueError("Provide at least 1 property name.")

    if len(set(properties)) != num_targets:
        raise ValueError("Property names must not be repeated.")

    print(f"Number of target properties: {num_targets}")
    print(f"Target properties: {', '.join(properties)}")
    print(
        "Descriptor counts to run: "
        + ", ".join(str(value) for value in descriptor_counts)
    )

    # ----- Load feature vectors only once -----
    feature_df = pd.read_csv(feature_file)
    feature_df.columns = feature_df.columns.str.strip()

    if "CID" not in feature_df.columns:
        raise ValueError("Feature file must contain a 'CID' column.")

    feature_df["CID"] = feature_df["CID"].astype(str)

    # ----- Load normalized targets only once -----
    target_norm_df = pd.read_csv(norm_target_file)
    target_norm_df.columns = target_norm_df.columns.str.strip()

    if "CID" not in target_norm_df.columns:
        raise ValueError(
            "Normalized target file must contain a 'CID' column."
        )

    target_norm_df["CID"] = target_norm_df["CID"].astype(str)

    missing_norm = [
        column for column in properties
        if column not in target_norm_df.columns
    ]
    if missing_norm:
        raise ValueError(
            "Columns missing from normalized target file: "
            + ", ".join(missing_norm)
        )

    target_norm_df = target_norm_df[["CID"] + properties]
    target_norm_df = target_norm_df.dropna(subset=properties)
    target_norm_df = target_norm_df.drop_duplicates(
        subset=["CID"],
        keep="first",
    )

    # ----- Load original/unscaled targets only once -----
    target_orig_df = pd.read_csv(unscaled_target_file)
    target_orig_df.columns = target_orig_df.columns.str.strip()

    if "CID" not in target_orig_df.columns:
        raise ValueError(
            "Unscaled target file must contain a 'CID' column."
        )

    target_orig_df["CID"] = target_orig_df["CID"].astype(str)

    missing_orig = [
        column for column in properties
        if column not in target_orig_df.columns
    ]
    if missing_orig:
        raise ValueError(
            "Columns missing from unscaled target file: "
            + ", ".join(missing_orig)
        )

    target_orig_df = target_orig_df[["CID"] + properties]
    target_orig_df = target_orig_df.dropna(subset=properties)
    target_orig_df = target_orig_df.drop_duplicates(
        subset=["CID"],
        keep="first",
    )

    # ----- Intersection of CIDs across all three files -----
    intersection = (
        set(feature_df["CID"])
        & set(target_norm_df["CID"])
        & set(target_orig_df["CID"])
    )

    if not intersection:
        raise ValueError(
            "No common CID values were found across the three input files."
        )

    feature_df_filtered = feature_df[
        feature_df["CID"].isin(intersection)
    ].copy()
    target_norm_filtered = target_norm_df[
        target_norm_df["CID"].isin(intersection)
    ].copy()
    target_orig_filtered = target_orig_df[
        target_orig_df["CID"].isin(intersection)
    ].copy()

    # Stable CID order.
    feature_df_filtered = feature_df_filtered.sort_values("CID")

    merged_norm = feature_df_filtered.merge(
        target_norm_filtered,
        on="CID",
        how="left",
    )

    original_column_names = {
        property_name: property_name + "_orig"
        for property_name in properties
    }
    target_orig_filtered = target_orig_filtered.rename(
        columns=original_column_names
    )

    merged_orig = feature_df_filtered.merge(
        target_orig_filtered,
        on="CID",
        how="left",
    )

    y_norm_all = merged_norm[properties].values

    y_orig_columns = [
        property_name + "_orig"
        for property_name in properties
    ]
    y_orig_all = merged_orig[y_orig_columns].values

    # ----- Global min/range for original targets -----
    y_min = y_orig_all.min(axis=0)
    y_max = y_orig_all.max(axis=0)
    y_range = y_max - y_min
    y_range[y_range == 0] = 1.0

    available_descriptors = feature_df_filtered.shape[1] - 1
    excessive_counts = [
        count for count in descriptor_counts
        if count > available_descriptors
    ]

    if excessive_counts:
        raise ValueError(
            f"The feature file contains {available_descriptors} descriptor "
            "columns, but the following requested count(s) are too large: "
            + ", ".join(str(value) for value in excessive_counts)
        )

    completed_runs = []
    failed_runs = []

    # ======================================================
    # Run the complete ANN workflow for each descriptor count
    # ======================================================
    for run_index, num_descriptors in enumerate(
        descriptor_counts,
        start=1,
    ):
        base_tag = (
            f"{'_'.join(properties)}_{num_descriptors}desc"
        )

        print("\n" + "=" * 74)
        print(
            f"DESCRIPTOR RUN {run_index}/{len(descriptor_counts)}"
        )
        print(f"Number of descriptors: {num_descriptors}")
        print(f"Targets: {', '.join(properties)}")
        print(f"Output tag: {base_tag}")
        print("=" * 74)

        descriptor_start_idx = 1  # CID is in column 0
        descriptor_end_idx = (
            descriptor_start_idx + num_descriptors
        )

        X = merged_norm.iloc[
            :,
            descriptor_start_idx:descriptor_end_idx,
        ].values

        # Reset seeds so each list item is equivalent to a fresh standalone run.
        reset_reproducibility_state()

        try:
            results = learn_ANN_eval(
                X,
                y_norm_all,
                y_orig_all,
                properties,
                num_descriptors,
                base_tag,
                y_min,
                y_range,
                rho=0.98,
                maxitr=10000,
                batch_size=None,
            )

            write_results_to_files(
                results,
                properties,
                base_tag,
            )

            completed_runs.append(num_descriptors)

            print(
                f"\nCompleted descriptor count: {num_descriptors}"
            )

        except KeyboardInterrupt:
            print(
                "\nProcessing stopped by the user during descriptor "
                f"count {num_descriptors}."
            )
            raise

        except Exception as error:
            failed_runs.append((num_descriptors, str(error)))
            print(
                f"\nERROR for {num_descriptors} descriptors: {error}"
            )
            print(
                "The script will continue with the next descriptor count."
            )

        finally:
            tf.keras.backend.clear_session()
            gc.collect()

    print("\n" + "=" * 74)
    print("ALL DESCRIPTOR-COUNT RUNS FINISHED")
    print("=" * 74)

    if completed_runs:
        print(
            "Completed counts: "
            + ", ".join(str(value) for value in completed_runs)
        )
    else:
        print("Completed counts: none")

    if failed_runs:
        print("Failed counts:")
        for count, error in failed_runs:
            print(f"  {count}: {error}")
    else:
        print("Failed counts: none")


def write_results_to_files(results, properties, base_tag):
    """Write a short summary of all key metrics."""
    if base_tag:
        results_filename = os.path.join(
            RESULT_DIR,
            f"results_R2_MAE_{base_tag}.txt",
        )
    else:
        results_filename = os.path.join(
            RESULT_DIR,
            "results_R2_MAE.txt",
        )

    with open(results_filename, "w", encoding="utf-8") as file:
        file.write(f"Targets: {', '.join(properties)}\n")
        file.write(
            f"Best Architecture: {results['best_architecture']}\n\n"
        )
        file.write("=== R2 (normalized target space) ===\n")
        file.write(
            "Median Train R2 per target: "
            f"{results['median_train_r2_per_dim']}\n"
        )
        file.write(
            "Median Testing R2 per target: "
            f"{results['median_test_r2_per_dim']}\n"
        )
        file.write(
            "Median combined Train R2 (average): "
            f"{results['median_train_r2_combined']}\n"
        )
        file.write(
            "Median combined Testing R2 (average): "
            f"{results['median_test_r2_combined']}\n\n"
        )
        file.write("=== MAE / NMAE (original units) ===\n")
        file.write(
            "Median Train MAE  (orig units): "
            f"{results['median_train_mae']}\n"
        )
        file.write(
            "Median Testing MAE (orig units): "
            f"{results['median_test_mae']}\n"
        )
        file.write(
            "Median combined Train MAE (orig units): "
            f"{results['median_train_mae_combined']}\n"
        )
        file.write(
            "Median combined Testing MAE (orig units): "
            f"{results['median_test_mae_combined']}\n\n"
        )
        file.write(
            f"Median Train NMAE: {results['median_train_nmae']}\n"
        )
        file.write(
            f"Median Testing NMAE: {results['median_test_nmae']}\n"
        )
        file.write(
            "Median combined Train NMAE: "
            f"{results['median_train_nmae_combined']}\n"
        )
        file.write(
            "Median combined Testing NMAE: "
            f"{results['median_test_nmae_combined']}\n"
        )

        file.write("\n=== Median-closest final-CV model (L1) ===\n")
        file.write(
            "Selection: minimum L1 distance between each fold's "
            "per-target Testing R2 vector and the per-target "
            "median Testing R2 vector.\n"
        )
        file.write(
            f"Selected Iteration: "
            f"{results['median_closest_iteration']}\n"
        )
        file.write(
            f"Selected Fold: {results['median_closest_fold']}\n"
        )
        file.write(
            f"Selected Fold Testing R2 per target: "
            f"{results['median_closest_test_r2']}\n"
        )
        file.write(
            f"L1 closeness score: "
            f"{results['median_closest_l1']}\n"
        )
        file.write(
            f"Median-closest weights: "
            f"{results['median_closest_weights_file']}\n"
        )
        file.write(
            f"Median-closest biases: "
            f"{results['median_closest_biases_file']}\n"
        )


# ==========================================================
#                   SCRIPT ENTRY POINT
# ==========================================================
if __name__ == "__main__":
    main(sys.argv)