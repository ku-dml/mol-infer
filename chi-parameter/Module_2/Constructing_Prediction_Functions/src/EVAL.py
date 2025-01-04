import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model

from src.ELR_actfun import *
from src import LR, ALR, ANN, RF, DT

import pulp 

import time,sys,copy,itertools,math,warnings

warnings.simplefilter('ignore')

####################################################

Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 10000
ZERO_TOL = 0.000001

STEP = 20

MAX_NUM = 1000000

####################################################

def eval_separate_plane(
    x_D1, y_D1_org, fv_list_D1, log_filename_D1, params_D1, output_prefix_D1,
    x_D2, y_D2_org, fv_list_D2, log_filename_D2, params_D2, output_prefix_D2
):

    R2_train = []
    R2_test = [] 

    R2_train_D1 = []
    R2_test_D1 = []
    R2_train_D2 = []
    R2_test_D2 = []

    sum_time_D1 = 0
    sum_time_D2 = 0

    x_sel_D1 = x_D1
    x_sel_D2 = x_D2

    # re-normalization
    min_y_D1 = np.amin(y_D1_org)
    max_y_D1 = np.amax(y_D1_org)
    y_D1 = (y_D1_org - min_y_D1) / (max_y_D1 - min_y_D1)

    min_y_D2 = np.amin(y_D2_org)
    max_y_D2 = np.amax(y_D2_org)
    y_D2 = (y_D2_org - min_y_D2) / (max_y_D2 - min_y_D2)

    with open(log_filename_D1, 'a') as f:
        f.write(f"EVALUATION:\n")
    with open(log_filename_D2, 'a') as f:
        f.write(f"EVALUATION:\n")

    # if params_D1["way"] == "ALR":
    #     epsilon_result_D1, alpha_result_D1, selected_D1, I_minus_D1, I_plus_D1 = ALR.learn_ALR_pre_var1(x_sel_D1, y_D1, a=params_D1["lambda"])
    # if params_D2["way"] == "ALR":
    #     epsilon_result_D2, alpha_result_D2, selected_D2, I_minus_D2, I_plus_D2 = ALR.learn_ALR_pre_var1(x_sel_D2, y_D2, a=params_D2["lambda"])

    x = np.concatenate((x_D1[:, 0], x_D2[:, 0]), axis=0)

    seed_success = 0
    split_seed = 1

    while seed_success < Times_eva:
        kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)

        try:
            ff = 0
            for train, test in kf.split(x):
                ff += 1
                train_D1 = [_x for _x in train if _x < len(x_D1)]
                test_D1 = [_x for _x in test if _x < len(x_D1)]
                train_D2 = [_x - len(x_D1) for _x in train if _x >= len(x_D1)]
                test_D2 = [_x - len(x_D1) for _x in test if _x >= len(x_D1)]

                x_train_D1 = x_sel_D1[train_D1]
                y_train_D1 = y_D1[train_D1]
                x_test_D1 = x_sel_D1[test_D1]
                y_test_D1 = y_D1[test_D1]
                x_train_D2 = x_sel_D2[train_D2]
                y_train_D2 = y_D2[train_D2]
                x_test_D2 = x_sel_D2[test_D2]
                y_test_D2 = y_D2[test_D2]

                y_train = np.concatenate((y_D1_org[train_D1], y_D2_org[train_D2]), axis=None)
                y_test = np.concatenate((y_D1_org[test_D1], y_D2_org[test_D2]), axis=None)

                time_start = time.time()
                if params_D1["way"] == "Lasso":
                    output_filename = output_prefix_D1 + f"{split_seed}_{ff}_linreg.txt"
                    lr, y_predict_train_D1, y_predict_test_D1, _ = \
                        LR.learn_LLR_eval(x_train_D1, y_train_D1, x_test_D1, y_test_D1, a=params_D1["lambda"])
                    with open(output_filename, 'w') as f:
                        f.write(f"{x_sel_D1.shape[1]}\n")
                        for _w in lr.coef_:
                            f.write(f"{_w} ")
                        f.write(f"\n{lr.intercept_}\n")
                elif params_D1["way"] == "BSP":
                    output_filename = output_prefix_D1 + f"{split_seed}_{ff}_linreg.txt"
                    lr, y_predict_train_D1, y_predict_test_D1 = \
                        LR.learn_MLR_eval(x_train_D1, y_train_D1, x_test_D1, y_test_D1)
                    with open(output_filename, 'w') as f:
                        f.write(f"{x_sel_D1.shape[1]}\n")
                        for _w in lr.coef_:
                            f.write(f"{_w} ")
                        f.write(f"\n{lr.intercept_}\n")
                elif params_D1["way"] == "ANN":
                    output_filename = output_prefix_D1 + f"{split_seed}_{ff}"
                    ann, y_predict_train_D1, y_predict_test_D1 = \
                        ANN.learn_ANN_eval(x_train_D1, y_train_D1, x_test_D1, y_test_D1, 
                            arch=params_D1["arch"], rho=params_D1["rho"], maxitr=params_D1["maxitr"])
                    ANN.write_weights_biases(ann[0], output_filename)
                elif params_D1["way"] == "RF":
                    params_D1_for_RF = {k: v for k, v in params_D1.items() if k != "way" and k != "selected_fv_filename"}
                    _, _, y_predict_train_D1, y_predict_test_D1 = \
                        RF.learn_random_forest(x_train_D1, y_train_D1, x_test_D1, y_test_D1, params=params_D1_for_RF)
                # elif params_D1["way"] == "ALR":
                #     _, y_predict_train_D1, y_predict_test_D1, _, _ = \
                #         ALR.learn_ALR_eval(x_train_D1, y_train_D1, x_test_D1, y_test_D1, fv_list_D1, params_D1["lambda"], 
                #             epsilon_result_D1, alpha_result_D1, selected_D1, I_minus_D1, I_plus_D1)
                elif params_D1["way"] == "DT":
                    params_D1_for_DT = {k: v for k, v in params_D1.items() if k != "way" and k != "selected_fv_filename"}
                    _, _, y_predict_train_D1, y_predict_test_D1 = \
                        DT.learn_decision_tree(x_train_D1, y_train_D1, x_test_D1, y_test_D1, params=params_D1_for_DT)
                time_end = time.time()
                time_D1 = time_end - time_start
                sum_time_D1 += time_D1

                time_start = time.time()
                if params_D2["way"] == "Lasso":
                    output_filename = output_prefix_D2 + f"{split_seed}_{ff}_linreg.txt"
                    lr, y_predict_train_D2, y_predict_test_D2, _ = \
                        LR.learn_LLR_eval(x_train_D2, y_train_D2, x_test_D2, y_test_D2, a=params_D2["lambda"])
                    with open(output_filename, 'w') as f:
                        f.write(f"{x_sel_D2.shape[1]}\n")
                        for _w in lr.coef_:
                            f.write(f"{_w} ")
                        f.write(f"\n{lr.intercept_}\n")
                elif params_D2["way"] == "BSP":
                    output_filename = output_prefix_D2 + f"{split_seed}_{ff}_linreg.txt"
                    lr, y_predict_train_D2, y_predict_test_D2 = \
                        LR.learn_MLR_eval(x_train_D2, y_train_D2, x_test_D2, y_test_D2)
                    with open(output_filename, 'w') as f:
                        f.write(f"{x_sel_D2.shape[1]}\n")
                        for _w in lr.coef_:
                            f.write(f"{_w} ")
                        f.write(f"\n{lr.intercept_}\n")
                elif params_D2["way"] == "ANN":
                    output_filename = output_prefix_D2 + f"{split_seed}_{ff}"
                    ann, y_predict_train_D2, y_predict_test_D2 = \
                        ANN.learn_ANN_eval(x_train_D2, y_train_D2, x_test_D2, y_test_D2, 
                            arch=params_D2["arch"], rho=params_D2["rho"], maxitr=params_D2["maxitr"])
                    ANN.write_weights_biases(ann[0], output_filename)
                elif params_D2["way"] == "RF":
                    params_D2_for_RF = {k: v for k, v in params_D2.items() if k != "way" and k != "selected_fv_filename"}
                    _, _, y_predict_train_D2, y_predict_test_D2 = \
                        RF.learn_random_forest(x_train_D2, y_train_D2, x_test_D2, y_test_D2, params=params_D2_for_RF)
                # elif params_D2["way"] == "ALR":
                #     _, y_predict_train_D2, y_predict_test_D2, _, _ = \
                #         ALR.learn_ALR_eval(x_train_D2, y_train_D2, x_test_D2, y_test_D2, fv_list_D2, params_D2["lambda"], 
                #             epsilon_result_D2, alpha_result_D2, selected_D2, I_minus_D2, I_plus_D2)
                elif params_D2["way"] == "DT":
                    params_D2_for_DT = {k: v for k, v in params_D2.items() if k != "way" and k != "selected_fv_filename"}
                    _, _, y_predict_train_D2, y_predict_test_D2 = \
                        DT.learn_decision_tree(x_train_D2, y_train_D2, x_test_D2, y_test_D2, params=params_D2_for_DT)
                time_end = time.time()
                time_D2 = time_end - time_start
                sum_time_D2 += time_D2

                y_predict_train_D1_org = y_predict_train_D1 * (max_y_D1 - min_y_D1) + min_y_D1
                y_predict_test_D1_org = y_predict_test_D1 * (max_y_D1 - min_y_D1) + min_y_D1
                y_predict_train_D2_org = y_predict_train_D2 * (max_y_D2 - min_y_D2) + min_y_D2
                y_predict_test_D2_org = y_predict_test_D2 * (max_y_D2 - min_y_D2) + min_y_D2
                
                y_predict_train = np.concatenate((y_predict_train_D1_org, y_predict_train_D2_org), axis=None)
                y_predict_test = np.concatenate((y_predict_test_D1_org, y_predict_test_D2_org), axis=None)

                r2train = r2_score(y_train, y_predict_train)
                r2test = r2_score(y_test, y_predict_test)

                time_tmp = time_D1 + time_D2

                if np.isnan(r2train):
                    r2train = 1.0
                if np.isnan(r2test):
                    r2test = 1.0

                with open(log_filename_D1, 'a') as f:
                    f.write(f"{split_seed}\t{ff}\t{r2train}\t{r2test}\t{time_tmp}\t{r2_score(y_train_D1, y_predict_train_D1)}\t{r2_score(y_test_D1, y_predict_test_D1)}\n")
                with open(log_filename_D2, 'a') as f:
                    f.write(f"{split_seed}\t{ff}\t{r2train}\t{r2test}\t{time_tmp}\t{r2_score(y_train_D2, y_predict_train_D2)}\t{r2_score(y_test_D2, y_predict_test_D2)}\n")
                # print(split_seed, r2train, r2test, r2train_mlr, r2test_mlr)

                R2_train.append(r2train)
                R2_test.append(r2test)

                r2train_D1 = r2_score(y_train_D1, y_predict_train_D1)
                r2test_D1 = r2_score(y_test_D1, y_predict_test_D1)
                r2train_D2 = r2_score(y_train_D2, y_predict_train_D2)
                r2test_D2 = r2_score(y_test_D2, y_predict_test_D2)

                if np.isnan(r2train_D1):
                    r2train_D1 = 1.0
                if np.isnan(r2test_D1):
                    r2test_D1 = 1.0
                if np.isnan(r2train_D2):
                    r2train_D2 = 1.0
                if np.isnan(r2test_D2):
                    r2test_D2 = 1.0

                R2_train_D1.append(r2train_D1)
                R2_test_D1.append(r2test_D1)
                R2_train_D2.append(r2train_D2)
                R2_test_D2.append(r2test_D2)

                # selected_linear.append(nonzero_linear)
                # selected_quadratic.append(nonzero_quadratic)

            seed_success += 1
        except:
            pass

        split_seed += 1

    return R2_train, R2_test, R2_train_D1, R2_test_D1, R2_train_D2, R2_test_D2, sum_time_D1, sum_time_D2

def eval_single(x, y, fv_list, log_filename, params):
    R2_train = []
    R2_test = [] 

    sum_time = 0

    x_sel = x

    # re-normalization
    min_y = np.amin(y)
    max_y = np.amax(y)
    y = (y - min_y) / (max_y - min_y)

    with open(log_filename, 'a') as f:
        f.write(f"EVALUATION:\n")
    
    if params["way"] == "ALR":
        epsilon_result, alpha_result, selected, I_minus, I_plus = ALR.learn_ALR_pre_var1(x_sel, y, a=params["lambda"])
   
    seed_success = 0
    split_seed = 1

    while seed_success < Times_eva:
        kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)

        try:
            ff = 0
            for train, test in kf.split(x):
                ff += 1
                x_train = x_sel[train]
                y_train = y[train]
                x_test = x_sel[test]
                y_test = y[test]

                time_start = time.time()
                if params["way"] == "Lasso":
                    _, y_predict_train, y_predict_test, _ = \
                        LR.learn_LLR_eval(x_train, y_train, x_test, y_test, a=params["lambda"])
                elif params["way"] == "BSP":
                    _, y_predict_train, y_predict_test = \
                        LR.learn_MLR_eval(x_train, y_train, x_test, y_test)
                elif params["way"] == "ANN":
                    _, y_predict_train, y_predict_test = \
                        ANN.learn_ANN_eval(x_train, y_train, x_test, y_test, 
                            arch=params["arch"], rho=params["rho"], maxitr=params["maxitr"])
                elif params["way"] == "RF":
                    params_for_RF = {k: v for k, v in params.items() if k != "way" and k != "selected_fv_filename"}
                    _, _, y_predict_train, y_predict_test = \
                        RF.learn_random_forest(x_train, y_train, x_test, y_test, params=params_for_RF)
                elif params["way"] == "ALR":
                    # _, y_predict_train, y_predict_test, _, _ = \
                    #     ALR.learn_ALR_eval(x_train, y_train, x_test, y_test, fv_list, params["lambda"], 
                    #         epsilon_result, alpha_result, selected, I_minus, I_plus)
                    _, r2train, r2test, _ = \
                        ALR.learn_ALR_pre(x_train, y_train, x_test, y_test, params["lambda"], 
                            epsilon_result, alpha_result, selected, I_minus, I_plus)
                elif params["way"] == "DT":
                    params_for_DT = {k: v for k, v in params.items() if k != "way" and k != "selected_fv_filename"}
                    _, _, y_predict_train, y_predict_test = \
                        DT.learn_decision_tree(x_train, y_train, x_test, y_test, params=params_for_DT)
                time_end = time.time()
                time_D1 = time_end - time_start
                sum_time += time_D1

                if params["way"] != "ALR":
                    r2train = r2_score(y_train, y_predict_train)
                    r2test = r2_score(y_test, y_predict_test)

                time_tmp = time_D1 

                if np.isnan(r2train):
                    r2train = 1.0
                if np.isnan(r2test):
                    r2test = 1.0

                with open(log_filename, 'a') as f:
                    f.write(f"{split_seed}\t{ff}\t{r2train}\t{r2test}\t{time_tmp}\n")
                    # f.write(f"{len(train_D1)}\t{len(train_D2)}\t{len(test_D1)}\t{len(test_D2)}\n")
                    # f.write(f"{r2_score(y_train, y_predict_train)}\t{r2_score(y_test, y_predict_test)}\n")
                # print(split_seed, r2train, r2test, r2train_mlr, r2test_mlr)

                R2_train.append(r2train)
                R2_test.append(r2test)

                # R2_train_D1.append(r2_score(y_train_D1, y_predict_train_D1))
                # R2_test_D1.append(r2_score(y_test_D1, y_predict_test_D1))
                # R2_train_D2.append(r2_score(y_train_D2, y_predict_train_D2))
                # R2_test_D2.append(r2_score(y_test_D2, y_predict_test_D2))

                # selected_linear.append(nonzero_linear)
                # selected_quadratic.append(nonzero_quadratic)

            seed_success += 1
        except:
            pass

        split_seed += 1

    return R2_train, R2_test, sum_time


