import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model

from src.ELR_actfun import *

import pulp 

import time,sys,copy,itertools,math,warnings

warnings.simplefilter('ignore')

####################################################
# IMPORTANT:
# Please specify the path of the cplex solver here
# CPLEX_PATH= \
# "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/cplex"
CPLEX_PATH= \
"/opt/cplex_12.10/cplex/bin/x86-64_linux/cplex"
# CPLEX_PATH= \
# "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"

CPLEX_MSG = False
CPLEX_TIMELIMIT = 0
solver_type = 1 # change 1 to 2 if use CBC solver
std_eps = 1e-5

# list of values of lambda used in Lasso (preliminary)
lambda_list = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1]
# lambda_list = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10, 25, 50, 100]

list_act_fun = [act_fun_0, act_fun_x2, act_fun_x2_alt]
# list_act_fun = [act_fun_0, act_fun_1, act_fun_2, act_fun_3, act_fun_4, act_fun_x2, act_fun_x2_alt, act_fun_half_p25, act_fun_half_p75, act_fun_half_q2_p25, act_fun_half_q2_p75]# test 4

num_act_fun = len(list_act_fun)

Times_pre = 1
Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 10000
ZERO_TOL = 0.000001

STEP = 20

MAX_NUM = 1000000

####################################################
def prepare_variables_ALR(
    K, train
):
    w = {(q, i): pulp.LpVariable(f"w({q},{i})", 0, cat=pulp.LpContinuous) for q in range(num_act_fun) for i in range(K)}
    b = pulp.LpVariable(f"b", cat=pulp.LpContinuous)
    alpha = {i: pulp.LpVariable(f"alpha({i})", 0, cat=pulp.LpContinuous) for i in range(num_act_fun)}
    # alpha = {i: pulp.LpVariable(f"alpha({i})", 0, cat=pulp.LpContinuous) for i in range(1, 3)}

    b_bar = pulp.LpVariable(f"b_bar", 0, cat=pulp.LpContinuous)
    Delta = {j: pulp.LpVariable(f"Delta({j})", 0, cat=pulp.LpContinuous) for j in range(len(train))}

    return w, b, alpha, b_bar, Delta

def build_constraints_ALR(
    MILP,
    x_train, y_train, K,
    llambda, I_minus, I_plus,
    w, b, alpha, b_bar, Delta
):

    MILP += pulp.lpSum(alpha[j] for j in range(num_act_fun)) == 1, "Lassp_var1"

    for i in range(len(x_train)):
        MILP += Delta[i] >= act_fun_all(y_train[i], list_act_fun, [alpha[q] for q in range(num_act_fun)]) - \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(q, j)] for q in range(num_act_fun)]) 
                    for j in I_plus]) + \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(q, j)] for q in range(num_act_fun)])
                    for j in I_minus]) - b, f"Lasso_var_2_{i}_1"
        MILP += -Delta[i] <= act_fun_all(y_train[i], list_act_fun, [alpha[q] for q in range(num_act_fun)]) - \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(q, j)] for q in range(num_act_fun)])
                    for j in I_plus]) + \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(q, j)] for q in range(num_act_fun)])
                    for j in I_minus]) - b, f"Lasso_var_2_{i}_2"

    MILP += b_bar >= b, f"Lasso_var2_3_1"
    MILP += b_bar >= -b, f"Lasso_var2_3_2"

    # Target function
    MILP += 1 / (2 * len(x_train)) * pulp.lpSum([Delta[i] for i in range(len(x_train))]) + \
            llambda * (pulp.lpSum([w[(q, j)] for q in range(num_act_fun) for j in range(K)]) + b_bar)

    return MILP

def prepare_variables_ALR_var1(
    K, train
):
    w = {(q, i): pulp.LpVariable(f"w({q},{i})", 0, cat=pulp.LpContinuous) for q in range(1) for i in range(K)}
    b = pulp.LpVariable(f"b", cat=pulp.LpContinuous)
    # alpha = {i: pulp.LpVariable(f"alpha({i})", 0, cat=pulp.LpContinuous) for i in range(num_act_fun)}
    # alpha = {i: pulp.LpVariable(f"alpha({i})", 0, cat=pulp.LpContinuous) for i in range(1, 3)}

    b_bar = pulp.LpVariable(f"b_bar", 0, cat=pulp.LpContinuous)
    Delta = {j: pulp.LpVariable(f"Delta({j})", 0, cat=pulp.LpContinuous) for j in range(len(train))}

    return w, b, b_bar, Delta

def build_constraints_ALR_var1(
    MILP,
    x_train, y_train, K,
    llambda, epsilon, I_minus, I_plus, selected,
    w, b, alpha, b_bar, Delta
):
    for i in range(len(x_train)):
        MILP += Delta[i] >= act_fun_all(y_train[i], list_act_fun, [alpha[q] for q in range(num_act_fun)]) - \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(0, j)] * epsilon[(q, j)] for q in range(num_act_fun)]) 
                    for j in I_plus if selected[j]]) + \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(0, j)] * epsilon[(q, j)] for q in range(num_act_fun)])
                    for j in I_minus if selected[j]]) - b, f"Lasso_var_2_{i}_1"
        MILP += -Delta[i] <= act_fun_all(y_train[i], list_act_fun, [alpha[q] for q in range(num_act_fun)]) - \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(0, j)] * epsilon[(q, j)] for q in range(num_act_fun)])
                    for j in I_plus if selected[j]]) + \
                pulp.lpSum([act_fun_all(x_train[i][j], list_act_fun, [w[(0, j)] * epsilon[(q, j)] for q in range(num_act_fun)])
                    for j in I_minus if selected[j]]) - b, f"Lasso_var_2_{i}_2"

    MILP += b_bar >= b, f"Lasso_var2_3_1"
    MILP += b_bar >= -b, f"Lasso_var2_3_2"

    # Target function
    MILP += 1 / (2 * len(x_train)) * pulp.lpSum([Delta[i] for i in range(len(x_train))]) + \
            llambda * (pulp.lpSum([w[(0, j)] * epsilon[(q, j)] for q in range(num_act_fun) for j in range(K)]) + b_bar)

    return MILP

def learn_ALR_pre_var1(x_train, y_train, a):
    numdata = x_train.shape[0]
    numfeature = x_train.shape[1]

    I_minus = list()
    I_plus = list()
    for j in range(numfeature):
        corr = np.corrcoef(x_train[:, j], y_train)
        if corr[0][1] > 0:
            I_plus.append(j)
        if corr[0][1] < 0:
            I_minus.append(j)

    #########################################################
    start_time = time.time()

    MILP = pulp.LpProblem("MILP_ALR", pulp.LpMinimize)
    w, b, alpha, b_bar, Delta = prepare_variables_ALR(numfeature, x_train)
    MILP = build_constraints_ALR(
        MILP, 
        x_train, y_train, numfeature, 
        a, I_minus, I_plus, 
        w, b, alpha, b_bar, Delta
    )

    # Solve MILP
    if solver_type == 1:
        if CPLEX_TIMELIMIT > 0:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG,
                               timeLimit = CPLEX_TIMELIMIT)
        else:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG)
        # print("Start Solving Using CPLEX...")
        MILP.solve(CPLEX)
        solve_end = time.time()
    else:
        # print("Start Solving Using Coin-OR...")
        MILP.solve()
        solve_end = time.time()

    # for _q in range(num_act_fun):
    #     for j in range(numfeature):
    #         print(f"w({_q}, {j}) = {w[(_q, j)].value()}")

    #########################################################
    epsilon_sum = dict()
    for j in range(numfeature):
        try:
            epsilon_sum[j] = sum(w[(_q, j)].value() for _q in range(num_act_fun))
        except:
            epsilon_sum[j] = 0
    epsilon_result = {(q, j): w[(q, j)].value() / epsilon_sum[j]
                            if abs(epsilon_sum[j]) > ZERO_TOL and w[(q, j)].value() is not None else 0 
                        for q in range(num_act_fun) for j in range(numfeature)}
    selected = {j: True if abs(epsilon_sum[j]) > ZERO_TOL else False for j in range(numfeature)}
    alpha_result = {i: alpha[i].value() for i in range(num_act_fun)}
    

    return epsilon_result, alpha_result, selected, I_minus, I_plus



def learn_ALR_pre(x_train, y_train, x_test, y_test, a, epsilon_result, alpha_result, selected, I_minus, I_plus):
    numdata = x_train.shape[0]
    numfeature = x_train.shape[1]

    #########################################################
    start_time = time.time()

    MILP = pulp.LpProblem("MILP_ALR", pulp.LpMinimize)
    w, b, b_bar, Delta = prepare_variables_ALR_var1(numfeature, x_train)
    MILP = build_constraints_ALR_var1(
        MILP, 
        x_train, y_train, numfeature, 
        a, epsilon_result, I_minus, I_plus, selected,
        w, b, alpha_result, b_bar, Delta
    )

    # Solve MILP
    if solver_type == 1:
        if CPLEX_TIMELIMIT > 0:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG,
                               timeLimit = CPLEX_TIMELIMIT)
        else:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG)
        # print("Start Solving Using CPLEX...")
        MILP.solve(CPLEX)
        solve_end = time.time()
    else:
        # print("Start Solving Using Coin-OR...")
        MILP.solve()
        solve_end = time.time()

    #########################################################

    reg_Lasso = linear_model.LinearRegression()
    reg_Lasso.fit(x_train, y_train)
    reg_Lasso.intercept_ = b.value()
    # for j in range(numfeature):
    #     print(j, numfeature, selected[j], w[(0, j)].value())
    reg_Lasso.coef_ = np.array([w[(0, j)].value() if j in I_plus and selected[j] and w[(0, j)].value() is not None
        else -w[(0, j)].value() if j in I_minus and selected[j] and w[(0, j)].value() is not None else 0
        for j in range(numfeature)])

    x_train = np.array([[act_fun_all(_x[j], list_act_fun, [epsilon_result[(q, j)] for q in range(num_act_fun)])
                for j in range(numfeature)] for _x in x_train])
    x_test = np.array([[act_fun_all(_x[j], list_act_fun, [epsilon_result[(q, j)] for q in range(num_act_fun)])
                for j in range(numfeature)] for _x in x_test])
    y_train = np.array([act_fun_all(_y, list_act_fun, [alpha_result[q] for q in range(num_act_fun)]) for _y in y_train])
    y_test = np.array([act_fun_all(_y, list_act_fun, [alpha_result[q] for q in range(num_act_fun)]) for _y in y_test])

    r2train = reg_Lasso.score(x_train, y_train)
    r2test = reg_Lasso.score(x_test, y_test)

    return reg_Lasso, r2train, r2test, solve_end - start_time



def calc_L1(pred, y):
    
    y_avg = np.mean(y)
    l1_sum = sum([abs(_y - y_avg) for _y in y])
    l1_u = sum([abs(_y - _p) for (_y, _p) in zip(y, pred)])
    
    return 1 - l1_u / l1_sum

def learn_ALR_eval(x_train, y_train, x_test, y_test, fv_list, a, epsilon_result, alpha_result, selected, I_minus, I_plus):
    numdata = x_train.shape[0]
    numfeature = x_train.shape[1]

    #########################################################
    start_time = time.time()

    MILP = pulp.LpProblem("MILP_ALR", pulp.LpMinimize)
    w, b, b_bar, Delta = prepare_variables_ALR_var1(numfeature, x_train)
    MILP = build_constraints_ALR_var1(
        MILP, 
        x_train, y_train, numfeature, 
        a, epsilon_result, I_minus, I_plus, selected,
        w, b, alpha_result, b_bar, Delta
    )

    # Solve MILP
    if solver_type == 1:
        if CPLEX_TIMELIMIT > 0:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG,
                               timeLimit = CPLEX_TIMELIMIT)
        else:
            CPLEX = pulp.CPLEX(path = CPLEX_PATH,
                               msg = CPLEX_MSG)
        # print("Start Solving Using CPLEX...")
        MILP.solve(CPLEX)
        solve_end = time.time()
    else:
        # print("Start Solving Using Coin-OR...")
        MILP.solve()
        solve_end = time.time()

    #########################################################

    reg_Lasso = linear_model.LinearRegression()
    reg_Lasso.fit(x_train, y_train)
    reg_Lasso.intercept_ = b.value()
    reg_Lasso.coef_ = np.array([w[(0, j)].value() if j in I_plus and selected[j] and w[(0, j)].value() is not None
        else -w[(0, j)].value() if j in I_minus and selected[j] and w[(0, j)].value() is not None else 0
        for j in range(numfeature)])

    x_train = np.array([[act_fun_all(_x[j], list_act_fun, [epsilon_result[(q, j)] for q in range(num_act_fun)])
                for j in range(numfeature)] for _x in x_train])
    x_test = np.array([[act_fun_all(_x[j], list_act_fun, [epsilon_result[(q, j)] for q in range(num_act_fun)])
                for j in range(numfeature)] for _x in x_test])
    y_train = np.array([act_fun_all(_y, list_act_fun, [alpha_result[q] for q in range(num_act_fun)]) for _y in y_train])
    y_test = np.array([act_fun_all(_y, list_act_fun, [alpha_result[q] for q in range(num_act_fun)]) for _y in y_test])

    r1train = calc_L1(reg_Lasso.predict(x_train), y_train)
    r1test = calc_L1(reg_Lasso.predict(x_test), y_test)
    r2train = reg_Lasso.score(x_train, y_train)
    r2test = reg_Lasso.score(x_test, y_test)

    nonzero_linear = 0
    nonzero_quadratic = 0
    for (i, w) in enumerate(reg_Lasso.coef_):
        if abs(w) >= ZERO_TOL:
            # print(fv_list[i])
            # selected_fv.append(fv_list[i])
            if "*" in fv_list[i]:
                nonzero_quadratic += 1
            else:
                nonzero_linear += 1

    return reg_Lasso, reg_Lasso.predict(x_train), reg_Lasso.predict(x_test), nonzero_linear, nonzero_quadratic

def ALR_pre(x, y, log_filename):

    TsR = {lmd: [] for lmd in lambda_list}
    TrR = {lmd: [] for lmd in lambda_list}

    for lmd in lambda_list:
        # _, selected_fv = learn_Lasso_eval(x, y, lmd)
        # x_sel = x[:, selected_fv]

        # if len(selected_fv) == 0:
        #     TsR[lmd] = [-1000000000]
        #     continue

        for split_seed in range(1, Times_pre+1):
            # 10000 is added because this is for preliminary experiments
            kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

            for train, test in kf.split(x):
                epsilon_result, alpha_result, selected, I_minus, I_plus = learn_ALR_pre_var1(x, y, a=lmd)
                _, r2train, r2test, ttt = learn_ALR_pre(x[train], y[train], x[test], y[test], lmd,
                    epsilon_result, alpha_result, selected, I_minus, I_plus)
                TsR[lmd].append(r2test)
                TrR[lmd].append(r2train)

                with open(log_filename, 'a') as f:
                    f.write(f"{lmd}\t{split_seed}\t{r2train}\t{r2test}\t{ttt}\n")
                # print(lmd, split_seed, r2train, r2test, ttt)

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

    for lmd in lambda_list_2:
        if len(TsR[lmd]) == Times_pre * Fold:
            continue

        # _, selected_fv = learn_Lasso_eval(x, y, lmd)
        # x_sel = x[:, selected_fv]

        # if len(selected_fv) == 0:
        #     TsR[lmd] = [-1000000000]
        #     continue

        for split_seed in range(1, Times_pre+1):
            # 10000 is added because this is for preliminary experiments
            kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

            for train, test in kf.split(x):
                epsilon_result, alpha_result, selected, I_minus, I_plus = learn_ALR_pre_var1(x, y, a=lmd)
                _, r2train, r2test, ttt = learn_ALR_pre(x[train], y[train], x[test], y[test], lmd,
                    epsilon_result, alpha_result, selected, I_minus, I_plus)
                TsR[lmd].append(r2test)
                TrR[lmd].append(r2train)

                with open(log_filename, 'a') as f:
                    f.write(f"{lmd}\t{split_seed}\t{r2train}\t{r2test}\t{ttt}\n")
                # print(lmd, split_seed, r2train, r2test, ttt)

    TsR_avg = {lmd: np.median(np.array(TsR[lmd])) for lmd in TsR}
    # print(TsR_avg)

    best_lambda = max(TsR_avg, key=TsR_avg.get)

    return best_lambda, np.median(np.array(TrR[best_lambda])), np.median(np.array(TsR[best_lambda]))
