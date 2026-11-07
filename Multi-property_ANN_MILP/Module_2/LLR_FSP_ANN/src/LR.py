import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model

import time,sys,copy,itertools,math,warnings

from src import CLSF, FUN

Times_pre = 5
Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 10000
ZERO_TOL = 0.000001

lambda_list = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1]
STEP = 20

warnings.simplefilter('ignore')

####################################################
def learn_MLR(x_train, y_train, x_test, y_test, metric='R2'):
    try:
        mlr = LinearRegression()
        mlr.fit(x_train, y_train)
    except:
        mlr = Lasso(alpha=0, max_iter=10**5)
        mlr.fit(x_train, y_train)
    # nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])
    # selected_fv = []

    # for (i, w) in enumerate(lasso.coef_):
    #     if abs(w) >= ZERO_TOL:
    #         # print(fv_list[i])
    #         selected_fv.append(i)

    if metric == 'R2':
        r2_train = mlr.score(x_train, y_train)
        r2_test = mlr.score(x_test, y_test)
    elif metric == 'BACC':
        r2_train, r2_test = CLSF.calc_bacc(mlr, x_train, y_train, x_test, y_test)
    elif metric == 'AUCROC':
        r2_train, r2_test = CLSF.calc_aucroc(mlr, x_train, y_train, x_test, y_test)

    return (mlr, r2_train, r2_test)

def learn_MLR_eval(x_train, y_train, x_test, y_test):
    try:
        mlr = LinearRegression()
        mlr.fit(x_train, y_train)
    except:
        mlr = Lasso(alpha=0, max_iter=10**5)
        mlr.fit(x_train, y_train)
    # nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])
    # selected_fv = []

    # for (i, w) in enumerate(lasso.coef_):
    #     if abs(w) >= ZERO_TOL:
    #         # print(fv_list[i])
    #         selected_fv.append(i)
    return (mlr, mlr.predict(x_train), mlr.predict(x_test))

def learn_Lasso_pre(x_train, y_train, x_test, y_test, a=1.0, metric='R2'):
    lasso = Lasso(alpha=a, max_iter=10**5)
    lasso.fit(x_train, y_train)
    if x_test.any() != None and y_test.any() != None:
        if metric == 'R2':
            r2train = lasso.score(x_train, y_train)
            r2test = lasso.score(x_test,y_test)
        elif metric == 'BACC':
            r2train, r2test = CLSF.calc_bacc(lasso, x_train, y_train, x_test, y_test)
        elif metric == 'AUCROC':
            r2train, r2test = CLSF.calc_aucroc(lasso, x_train, y_train, x_test, y_test)
    else:
        r2train = None
        r2test = None
    nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])
    return (lasso, nonzero, r2train, r2test)

def learn_Lasso_eval(x_train, y_train, a=1.0, threshold=None):
    lasso = Lasso(alpha=a, max_iter=10**5)
    lasso.fit(x_train, y_train)
    # nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])
    selected_fv = []

    for (i, w) in enumerate(lasso.coef_):
        if abs(w) >= ZERO_TOL:
            # print(fv_list[i])
            selected_fv.append(i)

    if threshold is not None:
        weight_array = np.array([abs(w) if abs(w) >= ZERO_TOL else 0 for w in lasso.coef_])
        if threshold > len(selected_fv):
            threshold = len(selected_fv)
        if threshold != 0:
            selected_fv = np.argpartition(weight_array, -threshold)[-threshold:] # select highest weights
        else:
            selected_fv = []

    # print(len(selected_fv), threshold)

    return (lasso, selected_fv)

def learn_LLR_eval(x_train, y_train, x_test, y_test, a=1.0):
    lasso = Lasso(alpha=a, max_iter=10**5)
    lasso.fit(x_train, y_train)
    # nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])
    selected_fv = []

    for (i, w) in enumerate(lasso.coef_):
        if abs(w) >= ZERO_TOL:
            # print(fv_list[i])
            selected_fv.append(i)

    # print(len(selected_fv), threshold)

    return lasso, lasso.predict(x_train), lasso.predict(x_test), selected_fv

def Lasso_pre(x, y, log_filename, lambda_list=lambda_list, STEP=STEP, metric="R2", Sbar=FUN.Status_Bar()):

    TsR = {lmd: [] for lmd in lambda_list}
    TrR = {lmd: [] for lmd in lambda_list}

    Sbar.set_total_times((len(lambda_list) + STEP - 2) * Times_pre * Fold)
    Sbar.output()

    for lmd in lambda_list:
        for split_seed in range(1, Times_pre+1):
            # 10000 is added because this is for preliminary experiments
            kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

            for train, test in kf.split(x):
                time_start = time.time()
                _, nonzero, r2train, r2test = learn_Lasso_pre(x[train], y[train], x[test], y[test], a=lmd, metric=metric)
                time_end = time.time()
                TsR[lmd].append(r2test)
                TrR[lmd].append(r2train)
                with open(log_filename, 'a') as f:
                    f.write(f"{lmd}, {split_seed}, {r2train}, {r2test}, {nonzero}, {time_end - time_start}\n")
                    f.close()  

                Sbar.append(time_end - time_start)

            with open(log_filename, 'a') as f:
                f.write(f"{lmd}, {np.mean(np.array(TsR[lmd]))}\n")
                f.close()  

            Sbar.output()

    TsR_avg = {lmd: np.mean(np.array(TsR[lmd])) for lmd in TsR}

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

    true_lambda_2 = 0
    for lmd in lambda_list_2:
        if lmd not in TsR:
            TsR[lmd] = []
            TrR[lmd] = []
            true_lambda_2 += 1
    Sbar.sub_total_times((STEP - 2 - true_lambda_2) * Times_pre * Fold)
    Sbar.output()

    for lmd in lambda_list_2:
        for split_seed in range(1, Times_pre+1):
            # 10000 is added because this is for preliminary experiments
            kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

            if len(TsR[lmd]) == Times_pre * Fold:
                continue
            for train, test in kf.split(x):
                time_start = time.time()
                _, nonzero, r2train, r2test = learn_Lasso_pre(x[train], y[train], x[test], y[test], a=lmd, metric=metric)
                time_end = time.time()
                TsR[lmd].append(r2test)
                TrR[lmd].append(r2train)
                with open(log_filename, 'a') as f:
                    f.write(f"{lmd}, {split_seed}, {r2train}, {r2test}, {nonzero}, {time_end - time_start}\n")
                    f.close()     

                Sbar.append(time_end - time_start)

            with open(log_filename, 'a') as f:
                f.write(f"{lmd}, {np.mean(np.array(TsR[lmd]))}\n")
                f.close()

            Sbar.output()

    TsR_avg = {lmd: np.median(np.array(TsR[lmd])) for lmd in TsR}
    TrR_avg = {lmd: np.median(np.array(TrR[lmd])) for lmd in TrR}

    best_lambda = max(TsR_avg, key=TsR_avg.get) 

    return best_lambda, TrR_avg[best_lambda], TsR_avg[best_lambda]

def write_file(x_sel_length, lr, output_filename):
    with open(output_filename, 'w') as f:
        f.write(f"{x_sel_length}\n")
        for _w in lr.coef_:
            f.write(f"{_w} ")
        f.write(f"\n{lr.intercept_}\n")
