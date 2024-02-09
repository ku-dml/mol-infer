#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############
#   -- v4, 10/21
############

import numpy as np
import pandas as pd
import pulp
import sys, time, warnings

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

MAX_NUM = 1000000

# minimun dataset size to run for separated subset
MIN_SIZE_SUBSET = 5

# minimun dataset size to run
MIN_SIZE = 200

# max percentage of instances inside intersenction
PERCENTAGE_INSEC = 0.1

####################################################
def read_dataset(data_csv, value_txt):
    
    ### read files ###
    # read the csv and the observed values
    x = pd.read_csv(data_csv, index_col=0)
    value = pd.read_csv(value_txt)

    ### prepare data set ###
    # prepare CIDs
    CIDs = np.array(x.index)
    # # prepare target, train, test arrays
    # # target = np.array(value['a'])
    # # construct dictionary: CID to feature vector
    # fv_dict = {}
    # for cid,row in zip(CIDs, x.values):
    #     fv_dict[cid] = row
    # # construct dictionary: CID to target value
    target_dict = {}
    for cid, val in zip(np.array(value['CID']), np.array(value['a'])):
        target_dict[cid] = val
    # # check CIDs: target_values_filename should contain all CIDs that appear in descriptors_filename
    # for cid in CIDs:
    #     if cid not in target_dict:
    #         sys.stderr.write('error: misses the target value of CID {}\n'.format(cid))
    #         exit(1)
    
    y = np.array([target_dict[cid] for cid in CIDs])

    return x, y, CIDs

def y_f(y, y_obs):
    ans = abs(y - y_obs) ** 2
    # ans = abs(y - y_obs) ** (0.1)
    return ans

def solve_MILP(x, y, y_obs):
    MILP = pulp.LpProblem("find_hyperplane", pulp.LpMinimize)

    K = len(x[0])
    N = len(x)

    w = {i: pulp.LpVariable(f"w({i})", cat=pulp.LpContinuous) for i in range(K)}
    b = pulp.LpVariable("b", cat=pulp.LpContinuous)
    delta = {i: pulp.LpVariable(f"delta({i})", 0, cat=pulp.LpContinuous) for i in range(N)}

    MILP += pulp.lpSum([delta[i] for i in range(N)]), "target"

    for i, (_x, _y) in enumerate(zip(x, y)):
        if _y <= y_obs:
            MILP += delta[i] >= pulp.lpSum(w[j] * _x[j] for j in range(K)) - b + y_f(_y, y_obs)
        else:
            MILP += delta[i] >= -(pulp.lpSum(w[j] * _x[j] for j in range(K)) - b) + y_f(_y, y_obs)
        
        if _y == 0.0:
            MILP += pulp.lpSum(w[j] * _x[j] for j in range(K)) - b <= 0

        if _y == 1.0:
            MILP += pulp.lpSum(w[j] * _x[j] for j in range(K)) - b >= 0

    # MILP.writeLP("test.lp")

    # Solve MILP
    solve_begin = time.time()
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

    w_value = {i: w[i].value() if w[i].value() is not None else 0 for i in range(K)}
    b_value = b.value()

    return w_value, b_value, solve_end - solve_begin

def check(w, b, x):
    K = len(x)
    return sum(w[i] * x[i] for i in range(K)) - b

def split_D1_D2(x, y, w, b):
    D_1_ind = list()
    D_2_ind = list()

    for i, (_x, _y) in enumerate(zip(x, y)):
        if check(w, b, _x) <= 0:
            D_1_ind.append(i)
        else:
            D_2_ind.append(i)

    return x[D_1_ind], y[D_1_ind], x[D_2_ind], y[D_2_ind]

def evaluation_exp(x, y, w, b, lambda_D1, lambda_D2, fv_list):

    time_start = time.time()

    R2_train = []
    R2_test = [] 
    selected_linear = []
    selected_quadratic = []

    for split_seed in range(1, Times_eva+1):
        kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)

        for train, test in kf.split(x):
            x_train = x[train]
            y_train = y[train]
            x_test = x[test]
            y_test = y[test]

            x_train_D1, y_train_D1, x_train_D2, y_train_D2 = split_D1_D2(x_train, y_train, w, b)
            x_test_D1, y_test_D1, x_test_D2, y_test_D2 = split_D1_D2(x_test, y_test, w, b)

            var_y_train = np.var(y_train) * len(y_train)
            var_y_test = np.var(y_test) * len(y_test)

            Lasso_D1, selected_fv_D1, nonzero_linear_D1, nonzero_quadratic_D1 = learn_Lasso_eval(x_train_D1, y_train_D1, fv_list, lambda_D1)
            Lasso_D2, selected_fv_D2, nonzero_linear_D2, nonzero_quadratic_D2 = learn_Lasso_eval(x_train_D2, y_train_D2, fv_list, lambda_D2)

            y_train_D1_predict = Lasso_D1.predict(x_train_D1)
            y_test_D1_predict = Lasso_D1.predict(x_test_D1)
            y_train_D2_predict = Lasso_D2.predict(x_train_D2)
            y_test_D2_predict = Lasso_D2.predict(x_test_D2)

            err_y_train = np.sum(np.square(y_train_D1 - y_train_D1_predict)) + np.sum(np.square(y_train_D2 - y_train_D2_predict))  
            err_y_test = np.sum(np.square(y_test_D1 - y_test_D1_predict)) + np.sum(np.square(y_test_D2 - y_test_D2_predict))  

            r2train = 1 - err_y_train / var_y_train 
            r2test = 1 - err_y_test / var_y_test          

            R2_train.append(r2train)
            R2_test.append(r2test)

            selected_fv = set(selected_fv_D1) | (set(selected_fv_D2))

            nonzero_linear = sum([1 for _fv in selected_fv if '*' not in _fv])
            nonzero_quadratic = sum([1 for _fv in selected_fv if '*' in _fv])

            selected_linear.append(nonzero_linear)
            selected_quadratic.append(nonzero_quadratic)

    time_end = time.time()

    return R2_train, R2_test, selected_linear, selected_quadratic, time_end - time_start

def main(argv):

    # list of values of theta
    theta_list = [0.001, 0.003, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
    # theta_list = [0.4]

    OUTPUT = False

    try:
        x_org, y, CIDs = read_dataset(argv[1], argv[2])
        try:
            theta_list = [float(argv[3])]
            OUTPUT = True
        except:
            pass
    except:
        sys.stderr.write("usage: {} (input_data.csv)(input_values.txt)\n\n".format(sys.argv[0]))
        exit(1)

    x = x_org.values

    # y_obs = np.median(y)

    # print(f"property\t#size\t#descs\ttheta\ta_min_1\ta_max_1\ta_min_2\ta_max_2\t#intersection\ttime_LP\tseparated\t|D_1|\t|D_2|\tlambda for D_1\ttime for D_1, pre\tlambda for D_2\ttime for D_2, pre\tMedian of train R^2\tMin of train R^2\tMax of train R^2\tMedian of test R^2\tMin of test R^2\tMax of test R^2\t#1-D selected desc\t#2-D selected desc\ttime for eval")

    min_y = np.amin(y)
    max_y = np.amax(y)
    y = (y - min_y) / (max_y - min_y)

    sp_loc = argv[1].find('_desc_norm')
    sp_loc2 = argv[1].rfind('/')
    if sp_loc2 == -1:
        prop_name = argv[1][:sp_loc]
    else:
        prop_name = argv[1][sp_loc2 + 1:sp_loc]

    data_size = len(y)

    filename = argv[1]

    # if len(y) < MIN_SIZE:
    #     output_str = f"{prop_name}\t{data_size}\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t"
    #     # print(output_str)
    #     return

    for theta in theta_list:

        try:
            w, b, solve_time = solve_MILP(x, y, theta)
        except:
            print(f"theta={theta}, error when solving LP")
            continue
        # print(w, b)
        # print(f"eps: {eps}")

        a_min_1 = MAX_NUM
        a_min_2 = MAX_NUM
        a_max_1 = -MAX_NUM
        a_max_2 = -MAX_NUM

        D_1_ind = list()   # indices of D_1 in D
        D_2_ind = list()   # indices of D_2 in D

        for i, (_x, _y) in enumerate(zip(x, y)):
            if check(w, b, _x) <= 0:
                D_1_ind.append(i)
                if _y < a_min_1:
                    a_min_1 = _y
                if _y > a_max_1:
                    a_max_1 = _y
            else:
                D_2_ind.append(i)
                if _y < a_min_2:
                    a_min_2 = _y
                if _y > a_max_2:
                    a_max_2 = _y

        num_insec = len([_y for _y in y if _y < a_max_1 and _y > a_min_2])

        output_str = f"{prop_name}\t{data_size}\t{x.shape[1]}\t{theta}\t{a_min_1}\t{a_max_1}\t{a_min_2}\t{a_max_2}\t{len(D_1_ind)}\t{len(D_2_ind)}\t{num_insec}\t{num_insec/data_size}\t{solve_time}\t"

        if a_max_1 < a_min_2:
            output_str += "Yes\t"
        else:
            output_str += "No\t"

        # output_str += f"{len(D_1_ind)}\t{len(D_2_ind)}\t"

        if OUTPUT:
            x_new_D1 = x_org.loc[CIDs[D_1_ind]]
            x_new_D2 = x_org.loc[CIDs[D_2_ind]]
            y_D1 = y[D_1_ind]
            y_D2 = y[D_2_ind]
            y_new_D1 = pd.DataFrame(data={'CID': CIDs[D_1_ind], 'a': y_D1}).set_index('CID')
            y_new_D2 = pd.DataFrame(data={'CID': CIDs[D_2_ind], 'a': y_D2}).set_index('CID')
            # y_new_D1 = pd.DataFrame(index={'CID': CIDs[D_1_ind]}, data={'a': y_D1})
            # y_new_D1 = pd.DataFrame(index={'CID': CIDs[D_2_ind]}, data={'a': y_D2})

            output_filename_x_D1 = filename[:sp_loc] + f"_theta{theta}_D1" + filename[sp_loc:]
            output_filename_x_D2 = filename[:sp_loc] + f"_theta{theta}_D2" + filename[sp_loc:]
            output_filename_y_D1 = filename[:sp_loc] + f"_theta{theta}_D1_values.txt"
            output_filename_y_D2 = filename[:sp_loc] + f"_theta{theta}_D2_values.txt"

            x_new_D1.to_csv(output_filename_x_D1, sep=',')
            x_new_D2.to_csv(output_filename_x_D2, sep=',')
            y_new_D1.to_csv(output_filename_y_D1, sep=',')
            y_new_D2.to_csv(output_filename_y_D2, sep=',')


        # if len(D_1_ind) < MIN_SIZE_SUBSET or len(D_2_ind) < MIN_SIZE_SUBSET or num_insec > PERCENTAGE_INSEC * data_size:
        #     output_str += "---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t"
        # else:

        #     x_D1 = x[D_1_ind]
        #     y_D1 = y[D_1_ind]
        #     x_D2 = x[D_2_ind]
        #     y_D2 = y[D_2_ind]

        #     try:
        #         lambda_D1, time_D1 = Lasso_pre(x_D1, y_D1)
        #         lambda_D2, time_D2 = Lasso_pre(x_D2, y_D2)

        #         try:
        #             R2_train, R2_test, nonzero_linear, nonzero_quadratic, time_eval = evaluation_exp(x, y, w, b, lambda_D1, lambda_D2, fv_list)

        #             output_str += f"{lambda_D1}\t{time_D1}\t"
        #             output_str += f"{lambda_D2}\t{time_D2}\t"
        #             output_str += f"{np.median(R2_train)}\t{np.amin(R2_train)}\t{np.amax(R2_train)}\t"
        #             output_str += f"{np.median(R2_test)}\t{np.amin(R2_test)}\t{np.amax(R2_test)}\t"
        #             output_str += f"{np.mean(nonzero_linear)}\t{np.mean(nonzero_quadratic)}\t{time_eval}\t"
        #         except:
        #             # evaluation exp fails because of small size of D_1 or D_2
        #             output_str += f"{lambda_D1}\t{time_D1}\t{lambda_D2}\t{time_D2}\t---\t---\t---\t---\t---\t---\t---\t---\t"
        #     except:
        #         # preliminary exp fails because of small size of D_1 or D_2
        #         output_str += "---\t---\t---\t---\t---\t---\t---\t---\t---\t---\t"

        print(output_str)
# 

if __name__ == "__main__":
    main(sys.argv)
    # main((0,"./data_regression_var0/At_large_var0_desc_norm.csv","./data_regression_var0/At_large_norm_values.txt"))
    # main((0,"GP-KRI-ADs_org_quadratic_desc_norm.csv","GP-KRI-ADs_BP_norm_values.txt"))



