#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.neural_network import MLPRegressor,MLPClassifier
from sklearn.feature_selection import SelectKBest, f_regression,f_classif
from sklearn.metrics import balanced_accuracy_score
from sklearn.linear_model import Lasso
from CLSF import calc_bacc,calc_aucroc
import sys, time, itertools, csv, math

ZERO_TOL = 0.000001
Times = 10
Fold = 5
metric = "BACC"
quad = False
MaxItr = 10000
def learn_Lasso(metric, x_train, y_train, x_test, y_test, Weight, a=1.0):
    lasso = Lasso(alpha=a, max_iter=10 ** 5)
    if metric == "r2":
        lasso.fit(x_train, y_train)
    else:
        lasso.fit(x_train, y_train, sample_weight=Weight)

    if metric == 'r2':
        r2_train = lasso.score(x_train, y_train)
        r2_test = lasso.score(x_test, y_test)
    elif metric == 'AUCROC':
        r2_train, r2_test = calc_aucroc(lasso, x_train, y_train, x_test, y_test)
    elif metric == 'BACC':
        r2_train, r2_test = calc_bacc(lasso, x_train, y_train, x_test, y_test)

    nonzero = len([w for w in lasso.coef_ if abs(w) >= ZERO_TOL])

    return lasso, nonzero, r2_train, r2_test

def adjust_lambda(Lambda,x,y):
    testR2_score_for_each_lambda = {}
    logLmd = math.log10(abs(Lambda))
    power_of_ten = math.floor(logLmd)
    for lmda_tmp in np.arange(start=Lambda - 10 ** power_of_ten + 10 ** (power_of_ten - 1),
                          stop=Lambda + 10 ** power_of_ten + 10 ** (power_of_ten - 1), step=10 ** (power_of_ten - 1)):
        lmda = float('%.2g' % lmda_tmp)
        Ten_test_r2 = []
        print("Lambda\t{}".format(lmda))
        for split_seed in range(1, Times + 1): #10
            kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)
            fold = 0
            Ts = []
            for train, test in kf.split(x):
                fold += 1
                c1 = len([cls for cls in y[train] if cls == 1])
                c0 = len(y[train]) - c1
                W = []
                for cls in y[train]:
                    if cls == 1:
                        W.append(1.0 / c1)
                    else:
                        W.append(1.0 / c0)
                _, _, _, r2test = learn_Lasso(metric,x[train], y[train], x[test], y[test], Weight=W, a=lmda)
                Ts.append(r2test)
            Ten_test_r2.extend(Ts) # 保存了10轮5个数字，一共50个数
        testR2_score_for_each_lambda[lmda] = np.median(Ten_test_r2) # 保存一个lambda的对应的50个值的median，一共保存20个数字
    lambda_with_highest_score = max(testR2_score_for_each_lambda, key=testR2_score_for_each_lambda.get)
    return lambda_with_highest_score, power_of_ten

def get_best_lmd(initial_lmd,x,y):

    selected_fv = []
    lambda_with_highest_score, power_of_ten = adjust_lambda(initial_lmd, x, y)
    while True:
        logLmd2 = math.log10(abs(lambda_with_highest_score))  #
        power_of_ten2 = math.floor(logLmd2)  # 计算获得的最佳的lambda的科学计数法的位数。
        if power_of_ten != power_of_ten2:
            lambda_with_highest_score, power_of_ten = adjust_lambda(lambda_with_highest_score, x, y)
        elif power_of_ten == power_of_ten2:
            lasso = Lasso(alpha = lambda_with_highest_score, max_iter=10 ** 5)
            lasso.fit(x,y)

            for (i, w) in enumerate(lasso.coef_):
                if abs(w) >= ZERO_TOL:
                    selected_fv.append(i)
            break

    return selected_fv

def read_dataset(data_csv, value_txt):
    ### read files ###
    # read the csv and the observed values
    fv = pd.read_csv(data_csv)
    value = pd.read_csv(value_txt)

    ### prepare data set ###
    # prepare CIDs
    CIDs = np.array(fv['CID'])
    # prepare target, train, test arrays
    target = np.array(value['a'])
    # construct dictionary: CID to feature vector
    fv_dict = {}
    for cid, row in zip(CIDs, fv.values[:, 1:]):  # [From:to(row),From:to(column)]
        fv_dict[cid] = row
    # construct dictionary: CID to target value
    target_dict = {}
    for cid, val in zip(np.array(value['CID']), np.array(value['a'])):
        target_dict[cid] = val
    # check CIDs: target_values_filename should contain all CIDs that appear in descriptors_filename
    for cid in CIDs:
        if cid not in target_dict:
            sys.stderr.write('error: {} misses the target value of CID {}\n'.format(value_txt, cid))
            exit(1)
    # construct x and y so that the CIDs are ordered in ascending order
    CIDs.sort()
    x = np.array([fv_dict[cid] for cid in CIDs])
    y = np.array([target_dict[cid] for cid in CIDs])
    return x, y

def get_basic_architectures():
    return ((10,), (30,20))

def get_architectures(K):
    A = []
    F = []
    for i in range(1,4):
        n = max(int(i*K/4), 10)
        if not n in F:
            A.append((n,))
            F.append(n)
    # the numbers of hidden nodes are same over layers
    for n in F:
        for H in range(2,6):
            A.append((n,)*H)
    # the numbers of hidden nodes decrease over layers
    for n in F:
        a = [n,]
        for H in range(2,5):
            num = max(int(a[-1]*3/4), 3)
            a.append(num)
            A.append(tuple(a))
    return A

def learn_ANN_pre(metric, x_train, y_train, x_test, y_test, arch):

    best_score_test = None
    score_train_best = None
    # R = []

    if metric == "r2":
        est = MLPRegressor(activation='relu', solver='adam',
                           alpha=0.00001, hidden_layer_sizes=arch,
                           random_state=1, early_stopping=False)
    else:
        est = MLPClassifier(activation='relu', solver='adam',
                            alpha=0.00008, hidden_layer_sizes=arch,
                            random_state=1, early_stopping=False)
    est.warm_start = False

    for t in range(10, MaxItr + 1, 50):
        est.max_iter = t
        est.fit(x_train, y_train)
        est.warm_start = True

        if metric == 'r2':
            score_train = est.score(x_train, y_train)
            score_test = est.score(x_test, y_test)
        else:
            score_train = balanced_accuracy_score(y_true = y_train, y_pred = est.predict(x_train))
            score_test = balanced_accuracy_score(y_true = y_test, y_pred = est.predict(x_test))
        # R.append((t, score_train, score_test))

    # for ele in R:
        if best_score_test is None:
            best_score_test = score_test
            score_train_best = score_train
        elif best_score_test < score_test:
            best_score_test = score_test
            score_train_best = score_train

    return score_train_best, best_score_test

def get_basic_architectures():
    return ((10,), (30,20))

def kbest(x, y, k, metric):
    if metric == "r2":
        selector = SelectKBest(score_func=f_regression, k=k)
    else:
        selector = SelectKBest(score_func=f_classif, k=k)
    selector.fit(x, y)
    return selector.get_support(indices=True)


def tuning_1(x, y):
    best_test_score = None
    best_params = None
    fv_size = x.shape[1]  #
    k_cand = [math.floor(fv_size * i / 10) for i in range(5, 11)]  #
    Arch = get_architectures(len(x[0]))
    Arch += get_basic_architectures()

    for k in k_cand:  #

        arr = kbest(x, y, k, metric)
        x_selected = x[:, arr]

        for arch in Arch:

            R2_test = []

            for split_seed in range(1, Times + 1):
                kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)
                for train, test in kf.split(x_selected):
                    score_train, score_test = learn_ANN_pre(metric, x_selected[train], y[train], x_selected[test],
                                                            y[test], arch)
                    R2_test.append(score_test)

            if best_test_score is None:
                best_test_score = np.median(R2_test)
                best_params = arch
                x_final = x_selected

            elif np.median(R2_test) > best_test_score:
                best_test_score = np.median(R2_test)
                best_params = arch
                x_final = x_selected

    return best_params, x_final

def tuning(x_selected, y):

    best_test_score = None
    best_params = None
    # fv_size = x_selected.shape[1]

    Arch = get_architectures(len(x_selected[0]))
    Arch += get_basic_architectures()

    for arch in Arch:

        R2_test = []

        for split_seed in range(1, Times + 1):
            kf = KFold(n_splits=Fold, shuffle=True, random_state = split_seed)
            for train, test in kf.split(x_selected):
                score_train, score_test = learn_ANN_pre(metric, x_selected[train], y[train], x_selected[test], y[test], arch)
                R2_test.append(score_test)

        if best_test_score is None:
            best_test_score = np.median(R2_test)
            best_params = arch

        elif np.median(R2_test) > best_test_score:
            best_test_score = np.median(R2_test)
            best_params = arch

    return best_params

def train_ANN(x_selected,y,params):

    excel = [[] for _ in range(0,19)]
    Ten_test_r2 = []
    Ten_train_r2 = []
    Ten_time = []

    excel[17].append("params")
    excel[17].append(params)
    excel[17].append(metric)
    excel[18] = ["", "", "Round1", "Round2", "Round3", "Round4", "Round5", "Average", "Median", "Min", "Max"]

    # arr = kbest(x, y, params[0], metric)
    # x_selected = x[:, arr]

    for split_seed in range(1, Times + 1):
        kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)

        Tr, Tr_xls = [], ["", "Train"]
        Ts, Ts_xls = [], ["", "Test"]
        Tim, Tim_xls = [], ["", "Time(s)"]

        Tr_xls[0] = f"5-CV({split_seed})"

        for train, test in kf.split(x_selected):
            start_time = time.time()
            score_train, score_test = learn_ANN_pre(metric, x_selected[train],y[train], x_selected[test], y[test], params)
            comp_time = time.time() - start_time
            Tr.append(score_train)
            Ts.append(score_test)
            Tim.append(comp_time)

        for ele in Tr:
            Tr_xls.append("{:.6f}".format(ele))
        Tr_xls.append("{:.6f}".format(np.mean(Tr)))
        Tr_xls.append("{:.6f}".format(np.median(Tr)))
        Tr_xls.append("{:.6f}".format(np.min(Tr)))
        Tr_xls.append("{:.6f}".format(np.max(Tr)))
        excel.append(Tr_xls)

        for ele in Ts:
            Ts_xls.append("{:.6f}".format(ele))
        Ts_xls.append("{:.6f}".format(np.mean(Ts)))
        Ts_xls.append("{:.6f}".format(np.median(Ts)))
        Ts_xls.append("{:.6f}".format(np.min(Ts)))
        Ts_xls.append("{:.6f}".format(np.max(Ts)))
        excel.append(Ts_xls)

        for ele in Tim:
            Tim_xls.append("{:.6f}".format(ele))
        Tim_xls.append("{:.6f}".format(np.mean(Tim)))
        Tim_xls.append("{:.6f}".format(np.median(Tim)))
        Tim_xls.append("{:.6f}".format(np.min(Tim)))
        Tim_xls.append("{:.6f}".format(np.max(Tim)))
        excel.append(Tim_xls)

        print("{}\tTrain".format(split_seed), end="")
        for v in Tr:
            print("\t{:.6f}".format(v), end="")
        print()
        print(" \tTest", end="")
        for v in Ts:
            print("\t{:.6f}".format(v), end="")
        print()
        print(" \tTime", end="")
        for v in Tim:
            print("\t{:.6f}".format(v), end="")
        print()

        Ten_test_r2.extend(Ts)
        Ten_train_r2.extend(Tr)
        Ten_time.extend(Tim)

    scoreMax = '{:.6f}'.format(np.max(Ten_test_r2))
    scoreMin = '{:.6f}'.format(np.min(Ten_test_r2))
    median_score ='{:.6f}'.format(np.median(Ten_test_r2))

    excel.append(["", "", "", "", "", "", "Train R^2 over 50 trials", "{:.6f}".format(np.mean(Ten_train_r2)),
                  "{:.6f}".format(np.median(Ten_train_r2)), "{:.6f}".format(np.min(Ten_train_r2)),
                  "{:.6f}".format(np.max(Ten_train_r2))])
    excel.append(["", "", "", "", "", "", "Test R^2 over 50 trials", "{:.6f}".format(np.mean(Ten_test_r2)),
                  "{:.6f}".format(np.median(Ten_test_r2)), "{:.6f}".format(np.min(Ten_test_r2)),
                  "{:.6f}".format(np.max(Ten_test_r2))])
    excel.append(["", "", "", "", "", "", "Time(s) over 50 trials", "{:.6f}".format(np.mean(Ten_time)),
                  "{:.6f}".format(np.median(Ten_time)), "{:.6f}".format(np.min(Ten_time)),
                  "{:.6f}".format(np.max(Ten_time))])

    return median_score, scoreMax, scoreMin, excel

def output_result(file_path,output_path):
    result_list = []
    result_list.append(["dataset","Atom-Valence","#of feasible molecules","# of descriptors","# of select descriptors","Architecture","Median of test R^2 over 50 trials ","Min of test R^2 over 50 trials ","Max of test R^2 over 50 trials"])
    with open(file_path, 'r', newline='', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        file_list = [row for row in csv_reader]
        file.close()
    i = 0
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        for data,value in file_list:
            name = data.split('/')[-2]
            print(data)
            x, y = read_dataset(data, value)
            # selected_fv_index = get_best_lmd(1,x,y)
            # if len(selected_fv_index) == 0:
            #     x_selected = x
            # else:
            #     x_selected = x[:, selected_fv_index]
            # best_params = tuning(x_selected,y)
            best_params, x_selected = tuning_1(x, y)
            sub_result = []
            sub_result.append(name)
            sub_result.append("")
            sub_result.append(len(x_selected))
            sub_result.append(len(x[0]))
            sub_result.append(len(x_selected[0]))
            p = list(best_params)
            p.append(1)
            p.insert(0, len(x_selected[0]))
            sub_result.append(p)
            median_score, scoreMax, scoreMin, excel = train_ANN(x_selected, y, best_params)
            excel[10].append("Number of molecules")
            excel[10].append(len(x))
            excel[10].append("Number of descriptor")
            excel[10].append(len(x[0]))
            excel = pd.DataFrame(excel)
            excel.to_excel(writer, sheet_name=f"{name}_{i}", index=False, header=False)
            i += 1
            sub_result.append(median_score)
            sub_result.append(scoreMin)
            sub_result.append(scoreMax)
            result_list.append(sub_result)
        result_list = pd.DataFrame(result_list)
        result_list.to_excel(writer, sheet_name="summary", index=False, header=False)

output_result("/Users/sbw/Desktop/ML/ML_data_path/Classification_norm_pos.csv","")