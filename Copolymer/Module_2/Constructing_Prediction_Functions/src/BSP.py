import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model

# from ELR_actfun import *
from src import LR, ALR

import pulp 

import time,sys,copy,itertools,math,warnings

warnings.simplefilter('ignore')

####################################################

Times_pre_1 = 5
Times_pre_2 = 5
Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 10000
ZERO_TOL = 0.000001

lambda_list = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1]
STEP = 20

MAX_NUM = 1000000

####################################################


def BSP(x, y, arr_D, best_D, best_train, best_median, best_min, log_filename):

    TsR = {i: [] for i in range(len(arr_D))}
    TrR = {i: [] for i in range(len(arr_D))}
    MED = {i: 0 for i in range(len(arr_D))}
    MIN = {i: 0 for i in range(len(arr_D))}
    TRA = {i: 0 for i in range(len(arr_D))}

    tmp_D = None
    tmp_median = None
    tmp_min = None

    with open(log_filename, 'a') as f:
        f.write(f"BSP, len={len(arr_D)}\n")  
        f.write(f"now_size\tto_remove\tsplit_seed\tr2train\tr2test\ttime\n")

    for i in range(len(arr_D)):
        x_sel = x[:, [arr_D[j] for j in range(len(arr_D)) if j != i]]

        for split_seed in range(1, Times_pre_2+1):
            # 10000 is added because this is for preliminary experiments
            kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

            for train, test in kf.split(x_sel):
                time_start = time.time()


                # change this to MLR
                # epsilon_result, alpha_result, selected, I_minus, I_plus = learn_ALR_pre_var1(x_sel, y, a=0)
                # _, r2train, r2test, ttt = learn_ALR_pre(x_sel[train], y[train], x_sel[test], y[test], 0,
                #     epsilon_result, alpha_result, selected, I_minus, I_plus)

                _, r2train, r2test = LR.learn_MLR(x_sel[train], y[train], x_sel[test], y[test])

                time_end = time.time()
                TsR[i].append(r2test)
                TrR[i].append(r2train)

                with open(log_filename, 'a') as f:
                    f.write(f"{len(arr_D)}\t{i}\t{split_seed}\t{r2train}\t{r2test}\t{time_end - time_start}\n")
                # print(lmd, split_seed, r2train, r2test, ttt)

        with open(log_filename, 'a') as f:
            f.write(f"BSP, SUMMARY: ")
            f.write(f"{len(arr_D)}\t{i}\t{np.median(np.array(TrR[i]))}\t{np.median(np.array(TsR[i]))}\n")

        MED[i] = math.floor(np.median(TsR[i])*10**4)/(10**4)
        MIN[i] = math.floor(np.amin(TsR[i])*10**4)/(10**4)
        TRA[i] = math.floor(np.median(TrR[i])*10**4)/(10**4)

        if tmp_D is None:
            tmp_D = [arr_D[j] for j in range(len(arr_D)) if j != i]
            tmp_median = MED[i]
            tmp_min = MIN[i]
            tmp_train = TRA[i]
        elif MED[i] > tmp_median:
            tmp_D = [arr_D[j] for j in range(len(arr_D)) if j != i]
            tmp_median = MED[i]
            tmp_min = MIN[i]
            tmp_train = TRA[i]
        elif MED[i] == tmp_median and MIN[i] > tmp_min:
            tmp_D = [arr_D[j] for j in range(len(arr_D)) if j != i]
            tmp_median = MED[i]
            tmp_min = MIN[i]
            tmp_train = TRA[i]

    if tmp_median > best_median:
        best_D = tmp_D[:]
        best_median = tmp_median
        best_min = tmp_min
        best_train = tmp_train
    elif tmp_median == best_median and tmp_min > best_min:
        best_D = tmp_D[:]
        best_median = tmp_median
        best_min = tmp_min
        best_train = tmp_train

    return tmp_D, best_D, best_train, best_median, best_min, tmp_median, tmp_min

def BSP_main(x, y, best_train, best_median, best_min, log_filename):
    x_size = x.shape[1]
    best_D = [j for j in range(x_size)]
    arr_D = best_D[:]

    with open(log_filename, 'a') as f:
        f.write(f"BSP_MAIN\n")  

    for _ in range(x_size - 1):
        best_median_tmp = best_median
        best_min_tmp = best_min
        best_train_tmp = best_train
        
        time_start = time.time()

        arr_D, best_D, best_train, best_median, best_min, tmp_median, tmp_min = BSP(x, y, arr_D, best_D, best_train, best_median, best_min, log_filename)
        time_end = time.time()
        
        with open(log_filename, 'a') as f:
            f.write(f"BSP_MAIN, SUMMARY: ")
            f.write(f"{len(arr_D) + 1}\t{arr_D}\t{best_train_tmp}\t{best_median_tmp}\t{best_min_tmp}\t{tmp_median}\t{tmp_min}\t{time_end - time_start}\n")

        # early halt
        # if tmp_median < best_median_tmp or (tmp_median == best_median_tmp and tmp_min <= best_min_tmp):
        #     break

    return best_D, best_train, best_median, best_min

def LLR_ALR_pre(x, y_org, fv_list, log_filename):
    y = np.array(y_org, copy=True)

    # re-normalization
    min_y = np.amin(y)
    max_y = np.amax(y)
    y = (y - min_y) / (max_y - min_y)

    with open(log_filename, 'a') as f:
        f.write(f"LAMBDA_1\n")
        f.write(f"lmd\tnonzero_linear\tnonzero_quad\tsplit_seed\tr2train\tr2test\ttime\n")

    TsR = {lmd: [] for lmd in lambda_list}
    TrR = {lmd: [] for lmd in lambda_list}

    for lmd in lambda_list:
        _, selected_fv = LR.learn_Lasso_eval(x, y, lmd)
        x_sel = x[:, selected_fv]

        if len(selected_fv) == 0:
            TsR[lmd] = [-1000000000]
            TrR[lmd] = [-1000000000]
            continue

        nonzero_linear = 0
        nonzero_quadratic = 0
        for i in selected_fv:
            if "*" in fv_list[i]:
                nonzero_quadratic += 1
            else:
                nonzero_linear += 1

        for split_seed in range(1, Times_pre_1+1):
            # 10000 is added because this is for preliminary experiments
            kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

            for train, test in kf.split(x_sel):
                time_start = time.time()

                # change this to MLR
                # epsilon_result, alpha_result, selected, I_minus, I_plus = learn_ALR_pre_var1(x_sel, y, a=0)
                # _, r2train, r2test, ttt = learn_ALR_pre(x_sel[train], y[train], x_sel[test], y[test], 0,
                #     epsilon_result, alpha_result, selected, I_minus, I_plus)

                _, r2train, r2test = LR.learn_MLR(x_sel[train], y[train], x_sel[test], y[test])

                time_end = time.time()
                TsR[lmd].append(r2test)
                TrR[lmd].append(r2train)


                with open(log_filename, 'a') as f:
                    f.write(f"{lmd}\t{nonzero_linear}\t{nonzero_quadratic}\t{split_seed}\t{r2train}\t{r2test}\t{time_end - time_start}\n")
                # print(lmd, split_seed, r2train, r2test, ttt)

        with open(log_filename, 'a') as f:
            f.write("SUMMARY:")
            f.write(f"{lmd}\t{nonzero_linear}\t{nonzero_quadratic}\t{np.median(np.array(TsR[lmd]))}\n")

    TsR_avg = {lmd: np.median(np.array(TsR[lmd])) for lmd in TsR}

    best_lambda = max(TsR_avg, key=TsR_avg.get)

    lambda_list_2 = list()

    if best_lambda == lambda_list[0]:
        d1 = 0
        d2 = lambda_list[0]
        d3 = lambda_list[1]

        d = (d3 - d2) / STEP
        lambda_list_2 = [d2 + d * i for i in range(STEP + 1)]

    elif best_lambda == lambda_list[-1]:
        d1 = lambda_list[-2]
        d2 = lambda_list[-2]
        d3 = lambda_list[-1]

        d = (d3 - d2) / STEP
        lambda_list_2 = [d2 + d * i for i in range(STEP + 1)]
    else:
        idx = [i for (i, lmd) in enumerate(lambda_list) if lmd == best_lambda]
        d1 = lambda_list[idx[0] - 1]
        d2 = lambda_list[idx[0]]
        d3 = lambda_list[idx[0] + 1]

        d = (d2 - d1)/ (STEP / 2)
        lambda_list_2 = [d1 + d * i for i in range(int(STEP / 2))]
        d = (d3 - d2)/ (STEP / 2)
        lambda_list_2_1 = [d2 + d * i for i in range(int(STEP / 2) + 1)]
        lambda_list_2.extend(lambda_list_2_1)

    for lmd in lambda_list_2:
        if lmd not in TsR:
            TsR[lmd] = []
            TrR[lmd] = []

    best_D = None
    best_median = None
    best_min = None

    with open(log_filename, 'a') as f:
        f.write(f"LAMBDA_2\n")
        f.write(f"lmd\tnonzero_linear\tnonzero_quad\tsplit_seed\tr2train\tr2test\ttime\n")

    for lmd in lambda_list_2:
        # if len(TsR[lmd]) == Times_pre * Fold:
        #     continue
        if len(TsR[lmd]) != 0:
            TsR[lmd] = []
            TrR[lmd] = []

        # limit size of |D_i| for BSP
        # THRESHOLD = math.floor(200 + 100 * (1000 / (len(y) + 250)))
        THRESHOLD = math.floor(150 + 10 * (1000 / (len(y) + 200)))

        _, selected_fv = LR.learn_Lasso_eval(x, y, lmd, THRESHOLD)
        x_sel = x[:, selected_fv]

        if len(selected_fv) == 0:
            TsR[lmd] = [-1000000000]
            TrR[lmd] = [-1000000000]
            continue

        nonzero_linear_before = 0
        nonzero_quadratic_before = 0
        for i in selected_fv:
            if "*" in fv_list[i]:
                nonzero_quadratic_before += 1
            else:
                nonzero_linear_before += 1

        for split_seed in range(1, Times_pre_2+1):
            # 10000 is added because this is for preliminary experiments
            kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

            for train, test in kf.split(x_sel):
                time_start = time.time()

                # change this to MLR
                # epsilon_result, alpha_result, selected, I_minus, I_plus = learn_ALR_pre_var1(x_sel, y, a=0)
                # _, r2train, r2test, ttt = learn_ALR_pre(x_sel[train], y[train], x_sel[test], y[test], 0,
                #     epsilon_result, alpha_result, selected, I_minus, I_plus)

                _, r2train, r2test = LR.learn_MLR(x_sel[train], y[train], x_sel[test], y[test])

                time_end = time.time()
                TsR[lmd].append(r2test)
                TrR[lmd].append(r2train)

                with open(log_filename, 'a') as f:
                    f.write(f"{lmd}\t{nonzero_linear_before}\t{nonzero_quadratic_before}\t{split_seed}\t{r2train}\t{r2test}\t{time_end - time_start}\n")
                # print(lmd, split_seed, r2train, r2test, ttt)

        with open(log_filename, 'a') as f:
            f.write(f"BEFORE BSP, THRESHOLD={THRESHOLD} {lmd}\t{nonzero_linear_before}\t{nonzero_quadratic_before}\t{np.median(np.array(TsR[lmd]))}\n")

        tmp_median = math.floor(np.median(TsR[lmd])*10**4)/(10**4)
        tmp_min = math.floor(np.amin(TsR[lmd])*10**4)/(10**4)
        tmp_train = math.floor(np.median(TrR[lmd])*10**4)/(10**4)

        tmp_D, tmp_train, tmp_median, tmp_min = BSP_main(x_sel, y, tmp_train, tmp_median, tmp_min, log_filename)
        # tmp_D = [i for i in range(len(selected_fv))]

        with open(log_filename, 'a') as f:
            f.write(f"AFTER BSP, {lmd}\t{selected_fv}\t{tmp_D}\n")

        nonzero_linear = 0
        nonzero_quadratic = 0
        for i in tmp_D:
            if "*" in fv_list[selected_fv[i]]:
                nonzero_quadratic += 1
            else:
                nonzero_linear += 1

        with open(log_filename, 'a') as f:
            f.write(f"AFTER BSP, {lmd}\t{nonzero_linear}\t{nonzero_quadratic}\t{tmp_median}\n")

        if best_median is None:
            best_D = [selected_fv[i] for i in tmp_D]
            best_median = tmp_median
            best_min = tmp_min  
            best_linear_before = nonzero_linear_before
            best_quadratic_before = nonzero_quadratic_before   
            best_train = tmp_train       
        elif tmp_median > best_median:
            best_D = [selected_fv[i] for i in tmp_D]
            best_median = tmp_median
            best_min = tmp_min
            best_linear_before = nonzero_linear_before
            best_quadratic_before = nonzero_quadratic_before
            best_train = tmp_train
        elif tmp_median == best_median and tmp_min > best_min:
            best_D = [selected_fv[i] for i in tmp_D]
            best_median = tmp_median
            best_min = tmp_min
            best_linear_before = nonzero_linear_before
            best_quadratic_before = nonzero_quadratic_before
            best_train = tmp_train

    return best_D, best_train, best_median, best_linear_before, best_quadratic_before


