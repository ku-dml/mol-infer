#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.linear_model import Lasso
import sys, time, csv

Times = 10
Fold = 5
ZERO_TOL = 0.000001
metric = "r2"
################################################
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
    return CIDs, x, y
################################################
def learn_Lasso(metric, x_train, y_train, x_test, y_test, a=1.0):
    lasso = Lasso(alpha=a, max_iter=10 ** 5)
    lasso.fit(x_train, y_train)

    if metric == 'r2':
        r2_train = lasso.score(x_train, y_train)
        r2_test = lasso.score(x_test, y_test)
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
        # print("Lambda\t{}".format(lmda))
        for split_seed in range(1, Times + 1): #10
            kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)
            fold = 0
            Ts = []
            for train, test in kf.split(x):#5
                fold += 1
                _, _, _, r2test = learn_Lasso(metric,x[train], y[train], x[test], y[test], a=lmda)
                Ts.append(r2test)
            Ten_test_r2.extend(Ts) 
        testR2_score_for_each_lambda[lmda] = np.median(Ten_test_r2)
    lambda_with_highest_score = max(testR2_score_for_each_lambda, key=testR2_score_for_each_lambda.get)
    return lambda_with_highest_score, power_of_ten

def train_Lasso(best_lambda,x,y):
    Ten_train_r2 = []
    Ten_test_r2 = []
    Ten_NonZ = []
    Ten_time = []
    excel = [[] for _ in range(0,19)]
    excel[17].append("Best_lambda:")
    excel[17].append(best_lambda)
    excel[17].append("Evaluation metric:")
    excel[17].append(metric)
    excel[18] = ["", "", "Round1", "Round2", "Round3", "Round4", "Round5", "Average", "Median", "Min", "Max"]
    for split_seed in range(1, Times + 1):

        kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)
        fold = 0
        Tr, Tr_xls = [], ["", "Train"]
        Ts, Ts_xls = [], ["", "Test"]
        Tim, Tim_xls = [], ["", "Time(s)"]
        NonZ, NonZ_xls = [], ["", "#Nonzero weight"]
        Tr_xls[0] = f"5-CV({split_seed})"
        for train, test in kf.split(x):
            fold += 1
            start_time = time.time()
            _, nonzero, r2train, r2test = learn_Lasso(metric, x[train], y[train], x[test], y[test], a=best_lambda)
            comp_time = time.time() - start_time
            Tr.append(r2train)
            Ts.append(r2test)
            Tim.append(comp_time)
            NonZ.append(nonzero)

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

        for ele in NonZ:
            NonZ_xls.append("{:.6f}".format(ele))
        NonZ_xls.append("{:.6f}".format(np.mean(NonZ)))
        NonZ_xls.append("{:.6f}".format(np.median(NonZ)))
        NonZ_xls.append("{:.6f}".format(np.min(NonZ)))
        NonZ_xls.append("{:.6f}".format(np.max(NonZ)))
        excel.append(NonZ_xls)

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
        print(" \tNonzero", end="")
        for v in NonZ:
            print("\t{}".format(v), end="")
        print()
        Ten_train_r2.extend(Tr)
        Ten_test_r2.extend(Ts)
        Ten_NonZ.extend(NonZ)
        Ten_time.extend(Tim)
    excel.append(["", "", "", "", "", "", "Train R^2 over 50 trials", "{:.6f}".format(np.mean(Ten_train_r2)),
                  "{:.6f}".format(np.median(Ten_train_r2)), "{:.6f}".format(np.min(Ten_train_r2)),
                  "{:.6f}".format(np.max(Ten_train_r2))])
    excel.append(["", "", "", "", "", "", "Test R^2 over 50 trials", "{:.6f}".format(np.mean(Ten_test_r2)),
                  "{:.6f}".format(np.median(Ten_test_r2)), "{:.6f}".format(np.min(Ten_test_r2)),
                  "{:.6f}".format(np.max(Ten_test_r2))])
    excel.append(["", "", "", "", "", "", "Time(s) over 50 trials", "{:.6f}".format(np.mean(Ten_time)),
                  "{:.6f}".format(np.median(Ten_time)), "{:.6f}".format(np.min(Ten_time)),
                  "{:.6f}".format(np.max(Ten_time))])
    excel.append(["", "", "", "", "", "", "#nonzero weights over 50 trials", "{:.6f}".format(np.mean(Ten_NonZ)),
                  "{:.6f}".format(np.median(Ten_NonZ)), "{:.6f}".format(np.min(Ten_NonZ)),
                  "{:.6f}".format(np.max(Ten_NonZ))])
    return Ten_test_r2, Ten_NonZ, Ten_train_r2, Ten_time, excel


def main(argv):
    fv_path = argv[1]
    output = argv[2]
    initial_lambda = float(argv[3])

    result = list()
    result.append(
        ["Dataset","#of feasible molecules","#of descriptors","Best lambda", "Median of test R^2 over 50 trials ","Min of test R^2 over 50 trials ",
         "Max of test R^2 over 50 trials","Avg of # of nonzero weights "])

    with open(fv_path, 'r', newline='', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        file_list = [row for row in csv_reader]
        file.close()

    i = 1

    with pd.ExcelWriter(output, engine='openpyxl') as writer:

        for data, value in file_list:
            name = data.split('/')[-2]
            sub_row = list()
            sub_row.append(name)
            try:
                CIDs, x, y = read_dataset(data, value)
                Initial_lambda = float(initial_lambda)
            except:
                sys.stderr.write("usage:(dataset_path_file.csv)(output)(lambda)\n\n")
                exit(1)
            lambda_with_highest_score, power_of_ten = adjust_lambda(Initial_lambda, x, y)
            while True:
                logLmd2 = math.log10(abs(lambda_with_highest_score))  #
                power_of_ten2 = math.floor(logLmd2)  
                if power_of_ten != power_of_ten2:
                    lambda_with_highest_score, power_of_ten = adjust_lambda(lambda_with_highest_score, x, y)
                elif power_of_ten == power_of_ten2:
                    All_test_r2, All_NonZ, All_train_r2, All_time, excel = train_Lasso(lambda_with_highest_score, x, y)
                    # excel[10].append("")
                    # excel[10].append(len(x))
                    # excel[10].append("")
                    # excel[10].append(len(x[0]))
                    excel = pd.DataFrame(excel)
                    excel.to_excel(writer, sheet_name=f"{name}_{i}", index=False,header=False)
                    i+=1
                    test_max = '{:.6f}'.format(np.max(All_test_r2))
                    test_min = '{:.6f}'.format(np.min(All_test_r2))
                    test_median = '{:.6f}'.format(np.median(All_test_r2))
                    AVGNonZ = '{:.2f}'.format(np.average(All_NonZ))
                    sub_row.append(len(x))
                    sub_row.append(len(x[0]))
                    sub_row.append(lambda_with_highest_score)
                    sub_row.append(test_median)
                    sub_row.append(test_min)
                    sub_row.append(test_max)
                    sub_row.append(AVGNonZ)
                    break
            result.append(sub_row)
        result = pd.DataFrame(result)
        result.to_excel(writer, sheet_name="summary", index=False, header=False)

if __name__ == '__main__':
    # main("/Sample_instance/Sample_dataset_path.txt", "/Sample_instance/Lasso_reg.xlsx", 1)
    main(sys.argv)