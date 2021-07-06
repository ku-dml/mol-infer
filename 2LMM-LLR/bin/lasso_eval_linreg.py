#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.linear_model import Lasso
import copy, sys, time

Times = 10
Fold = 5
ZERO_TOL = 0.000001

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
    for cid,row in zip(CIDs, fv.values[:,1:]):
        fv_dict[cid] = row
    # construct dictionary: CID to target value
    target_dict = {}
    for cid, val in zip(np.array(value['CID']), np.array(value['a'])):
        target_dict[cid] = val
    # check CIDs: target_values_filename should contain all CIDs that appear in descriptors_filename
    for cid in CIDs:
        if cid not in target_dict:
            sys.stderr.write('error: {} misses the target value of CID {}\n'.format(target_values_filename, cid))
            exit(1)
    # construct x and y so that the CIDs are ordered in ascending order
    CIDs.sort()
    x = np.array([fv_dict[cid] for cid in CIDs])
    y = np.array([target_dict[cid] for cid in CIDs])
    return (CIDs,x,y)

    
################################################
def learn_Lasso(x_train, y_train, x_test, y_test, a=1.0):
    lasso = Lasso(alpha=a, max_iter=10**5)
    lasso.fit(x_train, y_train)
    r2train = lasso.score(x_train,y_train)
    r2test = lasso.score(x_test,y_test)
    nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])
    return (lasso, nonzero, r2train, r2test)

################################################

try:
    CIDs, x, y = read_dataset(sys.argv[1], sys.argv[2])
    linreg = sys.argv[3]
    lmd = float(sys.argv[4])
except:
    sys.stderr.write("usage: {} (input_data.csv)(input_values.txt)(output_linreg.txt)(lambda)\n\n".format(sys.argv[0]))
    exit(1)

### experiment ###

print("Lambda\t{}".format(lmd))
f = open(sys.argv[1])
arr = f.readline().split(',')
K = len(arr)-1
print("NumDesc\t{}".format(K))
f.close()

Conf = []
    
for split_seed in range(1, Times+1):
    kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)

    fold = 0
    Tr = []
    Ts = []
    Tim = []
    NonZ = []
    for train, test in kf.split(x):
        fold += 1
        start_time = time.time()
        lasso, nonzero, r2train, r2test = learn_Lasso(x[train], y[train], x[test], y[test], a=lmd)
        comp_time = time.time() - start_time
        Tr.append(r2train)
        Ts.append(r2test)
        Tim.append(comp_time)
        NonZ.append(nonzero)
        #if best_conf == None or r2test > best_conf[1]:
        #    best_lasso = copy.deepcopy(lasso)
        #    best_conf = (r2train, r2test, comp_time, nonzero, split_seed, fold)
        Conf.append((r2train, r2test, comp_time, nonzero, split_seed, fold, copy.deepcopy(lasso)))
        
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
    
### output to XX_linreg.txt ###

Conf.sort(key=lambda a: -a[1])
conf = Conf[int(len(Conf)/2)]

fp = open(linreg, "w")
fp.write( "# data file : {}\n".format(sys.argv[1]))
fp.write( "# value file: {}\n".format(sys.argv[2]))
fp.write(f"# R^2 train:            {conf[0]}\n")
fp.write(f"# R^2 test:             {conf[1]}\n")
fp.write(f"# comp time:            {conf[2]}\n")
fp.write(f"# nonzeros:             {conf[3]}\n")
fp.write(f"# split_seed in KFold:  {conf[4]}\n")
fp.write(f"# fold:                 {conf[5]}\n")
fp.write(f"{K}\n")
flag = True
lasso = conf[-1]
for v in lasso.coef_:
    if flag == False:
        fp.write(" ")
    flag = False
    fp.write("{}".format(v))
fp.write("\n{}\n".format(lasso.intercept_))
fp.close()
