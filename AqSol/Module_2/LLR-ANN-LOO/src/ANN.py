import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model

from src import LR

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

lambda_list_lin = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1]
lambda_list_quad = [0.000001, 0.00001, 0.0001, 0.001]

Times_pre = 1
Times_eva = 10
Fold = 5
MaxItr = 10000
RANDOM_SEED_BASE = 1000

STEP_lin = 20
STEP_quad = 5
MAX_NUM = 1000000

# for ANN
V = [x/100 for x in range(-100,101)]
ITR_WEIGHT = 1.5

####################################################
def get_basic_architectures():
    return ((10,), (30,20))

def get_architectures(K):
    A = []
    F = []
    for i in range(1,3):
        n = max(int(i*K/3), 10)
        if not n in F:
            A.append((n,))
            F.append(n)
    # the numbers of hidden nodes are same over layers
    for n in F:
        for H in range(2,3):
            A.append((n,)*H)
    # the numbers of hidden nodes decrease over layers
    for n in F:
        a = [n,]
        for H in range(2,3):
            num = max(int(a[-1]*2/3), 3)
            a.append(num)
            A.append(tuple(a))
    return A


def LLR_ANN_Lasso_pre(
    x, y, log_filename, quad=True
):
    if quad:
        lambda_list = lambda_list_quad[:]
        STEP = STEP_quad
    else:
        lambda_list = lambda_list_lin[:]
        STEP = STEP_lin

    best_lambda, _, _ = LR.Lasso_pre(x, y, log_filename, lambda_list, STEP)

    return best_lambda

# select K descriptors according to KBest (F-value)
def shrink_dataset(x, y, best_lambda, fv_list):
    nonzero_linear = 0
    nonzero_quad = 0

    _, arr = LR.learn_Lasso_eval(x, y, best_lambda)

    if len(arr) == 0:
        arr = list(range(x.shape[1]))

    xsel = x[:,arr]

    for i in arr:
        if '*' in fv_list[i]:
            nonzero_quad += 1
        else:
            nonzero_linear += 1

    return xsel, arr, nonzero_linear, nonzero_quad

def learn_ANN_pre(x_train, y_train, x_test, y_test, arch, maxitr):

    R = []
    reg = MLPRegressor(activation='relu', solver='adam',
                       alpha=1e-5, hidden_layer_sizes=arch,
                       random_state=1, early_stopping=False)
    reg.warm_start = False
    start_time = time.time()
    # learn ANN, but stop the learning at itr=t in order to record stats
    for t in range(10, MaxItr+1, 10):
        reg.max_iter = t
        reg.fit(x_train, y_train)
        reg.warm_start = True                

        # calculate the prediction score (R^2)
        r2train = reg.score(x_train,y_train)
        # r2test = reg.score(x_test,y_test)
        r2test = reg.predict(x_test)[0]
        R.append((t,r2train, r2test))

    return R

def LLR_ANN_pre(x, y, fv_list, log_filename, quad=True):
    # dictionaly sanmed Best record best score for each seed and ANN architectures
    Best = {}
    output_str = ''

    best_lambda = LLR_ANN_Lasso_pre(x, y, log_filename, quad)
    xsel, arr, nonzero_linear, nonzero_quad = shrink_dataset(x, y, best_lambda, fv_list)
        
    Arch = get_architectures(len(xsel[0]))
    Arch += get_basic_architectures()

    n = x.shape[0]
    n_list = list(range(n))

    for arch in Arch:

        Res = {v: [] for v in V}

        # MaxTrR2 = None
        # for split_seed in range(1, Times_pre+1):
        #     kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)
        #     fold = 0
        #     for train, test in kf.split(xsel):
        #         fold += 1
        #         start_time = time.time()
        #         R = learn_ANN_pre(xsel[train], y[train], xsel[test], y[test], arch, MaxItr)
        #         comp_time = time.time() - start_time
        #         Res_per_round = dict()
        #         for v in V:
        #             Res_per_round[v] = None
        #         for itr,r2train,r2test in R:
        #             if MaxTrR2 == None or r2train > MaxTrR2:
        #                 MaxTrR2 = r2train
        #             for v in V:
        #                 if Res_per_round[v]==None and r2train>=v:
        #                     Res_per_round[v] = (itr,r2train,r2test)
        #         for v in V:
        #             Res[v].append(Res_per_round[v])

        y_predict_V = {v: [] for v in V}
        MaxTrR2 = None

        for i in range(n):
            start_time = time.time()

            train = n_list[:]
            train.remove(i)
            test = [i]
            
            R = learn_ANN_pre(xsel[train], y[train], xsel[test], y[test], arch, MaxItr)
            comp_time = time.time() - start_time
            Res_per_round = dict()
            for v in V:
                Res_per_round[v] = None
            for itr,r2train,r2test in R:
                if MaxTrR2 == None or r2train > MaxTrR2:
                    MaxTrR2 = r2train
                for v in V:
                    if Res_per_round[v]==None and r2train>=v:
                        Res_per_round[v] = (itr,r2train,r2test)
            for v in V:
                Res[v].append(Res_per_round[v])
                if Res_per_round[v] is not None:
                    y_predict_V[v].append(Res_per_round[v][2])
                else:
                    y_predict_V[v].append(None)

        Output = []
        MaxTsR2 = None
        for v in V:
            if not None in Res[v]:
                itr_mean = np.mean([a[0] for a in Res[v]])
                TrR2_mean = np.mean([a[1] for a in Res[v]])
                # TsR2_mean = np.mean([a[2] for a in Res[v]])
                TsR2_mean = r2_score(y, y_predict_V[v])
                Output.append( (v, itr_mean, TrR2_mean, TsR2_mean) )
                if MaxTsR2 == None or TsR2_mean > MaxTsR2:
                    MaxTsR2 = TsR2_mean

        archstr = "_".join(list(map(str,arch)))

        # --- header --- #
        output_str += f"#K:\t{len(x[0])}\n"
        output_str += f"#reduced_K:\t{len(xsel[0])}\n"
        output_str += f"#arch:\t{archstr}\n"

        # --- best result is output in advance --- #
        for T in Output:
            if T[3] >= MaxTsR2:
                output_str += f"#best:\t{T[0]}\t{T[1]}\t{T[2]}\t{T[3]}\n"
                Best[(best_lambda, arch, T[0], int(T[1] * ITR_WEIGHT))] = (T[3],arr, nonzero_linear, nonzero_quad, T[2])
                break
        # --- all results --- #
        for T in Output:
            output_str += f"{T[0]}\t{T[1]}\t{T[2]}\t{T[3]}"
            if T[3] >= MaxTsR2:
                output_str += "\t*"
            output_str += "\n"
        # --- MaxTrR2 --- #
        if len(Output)==0:
            output_str += f"#warning:\tMaxTrR2={MaxTrR2}\n"

        with open(log_filename, 'a') as f:
            f.write(f"{output_str}\n")
            f.close()

    if len(Best) == 0:
        best = None
    else:
        sorted_Best = sorted(Best.items(), key=lambda x:x[1][0], reverse=True)
        best = next(iter(sorted_Best))
    
    return best, best[1][4], best[1][0]

def learn_ANN_eval(x_train, y_train, x_test, y_test, arch, rho, maxitr):

    reg = MLPRegressor(activation='relu', solver='adam',
                       alpha=1e-5, hidden_layer_sizes=arch,
                       random_state=1, early_stopping=False)
    reg.warm_start = False
    start_time = time.time()
    # learn ANN, but stop the learning at itr=t in order to record stats
    for t in range(10, maxitr+1, 10):
        reg.max_iter = t
        reg.fit(x_train, y_train)
        reg.warm_start = True                
        # calculate the prediction score (R^2)
        r2train = reg.score(x_train,y_train)
        # r2test = reg.score(x_test,y_test)
        # R = (reg, t, r2train, r2test, time.time()-start_time)
        if r2train >= rho:
            break
            
    return reg, reg.predict(x_train), reg.predict(x_test)

# def ANN_evaluation_exp(x, y, arch, rho, maxitr):

#     R2_train = []
#     R2_test = [] 

#     for split_seed in range(1, Times_eva+1):
#         kf = KFold(n_splits=Fold, shuffle=True, random_state=split_seed)

#         for train, test in kf.split(selected_x):
#             x_train = x[train]
#             y_train = y[train]
#             x_test = x[test]
#             y_test = y[test]

#             _, maxitr, r2train, r2test, time = learn_ANN_eval(x_train, y_train, x_test, y_test, arch, rho, maxitr)  

#             R2_train.append(r2train)
#             R2_test.append(r2test)

#     return R2_train, R2_test

def write_weights_biases(reg, filename):
    """ 
    This function will write to disk 2 files, called
        "filename_weights.txt" and
        "filename_biases.txt"
    containing the weights and biases 
    of a trained artificial neural network reg given as an argument
    """
    # initialize separate filenames
    weights_filename = filename + "_weights.txt"
    biases_filename = filename + "_biases.txt"
    
    # Get the weights and biases from the trained MLP regressor
    weights = reg.coefs_
    biases = reg.intercepts_
    num_features = weights[0].shape[0]
    
    # Write the weights to file weights_filename
    with open(weights_filename, 'w') as f:
        f.write(str(num_features) + ' ')
        for i in range(len(reg.hidden_layer_sizes)):
            f.write(str(reg.hidden_layer_sizes[i]) + ' ')
        f.write('1\n')
        for item in weights:
            for i in range(item.shape[0]):
                for j in range(item.shape[1]):
                    if abs(item[i][j]) > 10**(-308):
                        f.write(str(item[i][j]) + ' ')
                    else:
                        f.write('0 ')
                f.write('\n')
    
    # Write the biases to a file biases_filename
    with open(biases_filename, 'w') as f:
        for item in biases:
            for i in range(item.shape[0]):
                f.write(str(item[i])+ '\n')
