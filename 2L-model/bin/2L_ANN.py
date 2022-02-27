#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  ANN learner for 2L (two-layered) model 

  Copyright 2021
  Discrete Mathematics Laboratory at Kyoto University (KU-DML)
  Released Under the MIT License
  https://opensource.org/licenses/mit-license.php
"""

############### import ###############
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
import time,sys,copy,itertools

import warnings
warnings.simplefilter('ignore') # ignore warnings

############### global variables ###############
CV = 5                # number of folds in cross-validation
ANN_SEED = (1,2)    # random seed set for MLPRegressor
SPLIT_SEED = (1,)  # random seed set for KFold
################################################

R_key = ["R2train", "R2test", "R2all",\
         "MAEtrain", "MAEtest", "MAEall",\
         "time", "reg"]
Usage = '''
usage: {} (input.csv)(input.txt)(output)(maxitr)(architecture)
                  
  - input.csv    ... CSV file generated in Module 1.
  - input.txt    ... CSV file that contains pairs of CID and an observed value.
  - output       ... The filestem of output files. 
  - maxitr       ... Num of iterations in neural network learning
  - architecture ... Num of nodes in hydden layer(s)

The program will generate output_{{biases,weights}}.txt.
These two files and the three files output in Module 1 are required in Module 3. 
'''.format(sys.argv[0])


############### write ANN data ###############
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
    

############################################################
def train_ANN(descriptors_filename, target_values_filename, architecture,\
              ANN_seed, split_seed, T):

    ########## preprocess ##########
    ### read files ###
    # read the training and target data
    fv = pd.read_csv(descriptors_filename) 
    value = pd.read_csv(target_values_filename)      

    ### prepare training set ###
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
    # obtain numbers of examples and features
    numdata = x.shape[0]
    numfeature = x.shape[1]

    ### prepare learning ###
    # initialize an ANN - MLP regressor
    reg = MLPRegressor(activation='relu', solver='adam',
                       alpha=1e-5, hidden_layer_sizes=architecture,
                       random_state=ANN_seed, early_stopping=False)
    # initalize array that stores the result
    R = {} # R[<key>][fold][t]
    for key in R_key:
        R[key] = []
        for fold in range(CV):
            R[key].append(dict())
    # separate the data randomly for cross-validation
    kf = KFold(n_splits=CV, shuffle=True, random_state=split_seed)
    fold = -1
    ### start learning experiments ###
    for train, test in kf.split(x):
        fold += 1
        x_train, x_test, y_train, y_test = x[train], x[test], y[train], y[test]
        reg.warm_start = False
        start = time.time()
        print("\n\n### (ANN_seed, split_seed)=({},{}), fold={}/{} ###".format(ANN_seed, split_seed, fold+1, CV))
        print("# t\ttrain\ttest\ttime")
        # learn ANN, but stop the learning at itr=t in order to record stats
        for t in T:
            reg.max_iter = t
            reg.fit(x_train, y_train)
            reg.warm_start = True                
            # obtain the prediction to compute MAE
            pred = reg.predict(x)
            pred_train = reg.predict(x_train)
            pred_test = reg.predict(x_test)
            # calculate the prediction score (R^2)
            R["R2train"][fold][t] = reg.score(x_train,y_train)
            R["R2test"][fold][t] = reg.score(x_test,y_test)
            R["R2all"][fold][t] = reg.score(x,y)
            # calculate MAE 
            R["MAEtrain"][fold][t] = mean_absolute_error(y_train,pred_train)
            R["MAEtest"][fold][t] = mean_absolute_error(y_test,pred_test)
            R["MAEall"][fold][t] = mean_absolute_error(y,pred) 
            # store time and ref
            R["time"][fold][t] = time.time() - start
            R["reg"][fold][t] = copy.deepcopy(reg)
            print("{}\t{:.4f}\t{:.4f}\t{:.4f}".format(t, R["R2train"][fold][t], R["R2test"][fold][t], R["time"][fold][t]))
  
    return R

def get_best_cv_and_t(Result, key, T):
    best_r = best_fold = best_t = best_ts = best_trial = None
    for r in Result.keys():
        R = Result[r]
        for t in T:
            x = [R[key][fold][t] for fold in range(CV)]
            avg = np.mean(x)
            if best_ts == None or avg > best_ts:
                best_r = r
                best_t = t
                best_ts = avg
                best_trial = max(x)
                best_fold = x.index(best_trial)
            
    return (best_r, best_fold, best_t, best_ts, best_trial)

############################################################
def main(argv):
    if len(argv)<6:
        print(Usage)
        sys.exit()

    ##### Parse the command line arguments #####
    descriptors_filename = argv[1]
    target_values_filename = argv[2]
    output_filename = argv[3]
    maxitr = int(argv[4])
    architecture = tuple(int(a) for a in argv[5:])
    
    # initialize recording step
    T = [t for t in range(100, min(1000, maxitr), 100)]
    T += [t for t in range(1000, maxitr+1, 1000)]
    
    ##### Perform 5-fold validation with the given training data #####
    Result = {}
    SEED = list(itertools.product(ANN_SEED, SPLIT_SEED))
    for (ANN_seed, split_seed) in SEED:
        Result[(ANN_seed, split_seed)] = train_ANN(descriptors_filename,
                                                   target_values_filename,
                                                   architecture,
                                                   ANN_seed, split_seed,
                                                   T)
    (best_r, best_fold, best_t, best_ts, best_trial) = get_best_cv_and_t(Result, "R2test", T)

    ANN_seed, split_seed = best_r
    tr = Result[(ANN_seed,split_seed)]["R2train"][best_fold][best_t]
    tim = Result[(ANN_seed,split_seed)]["time"][best_fold][best_t]
    print("\n\n========== COMPLETE ==========\n")
    print("Best R^2 value among all cross-validations:\t{:.6f}\n".format(best_ts))
    print("The performance of the selected neural network:")
    print("\tR^2 value for training set:\t{:.6f}".format(tr))
    print("\tR^2 value for test set:\t{:.6f}\n".format(best_trial))
    #print("# TIME:\t{}".format(tim))
    #print("# ANN_SEED:\t{}".format(ANN_seed))
    #print("# SPLIT_SEED:\t{}".format(split_seed))
    #print("# FOLD:\t{}".format(best_fold))
    #print("# ITR:\t{}".format(best_t))
    
    best_regressor = Result[(ANN_seed,split_seed)]["reg"][best_fold][best_t]
    
    # Write the weights and biases of the regressor to files
    write_weights_biases(best_regressor, output_filename)
    
    print ("The data for the selected neural network is stored in")
    print ("{}_biases.txt and {}_weights.txt.".format(output_filename, output_filename))
    print ("They are required in Module 3.\n")
    
main(sys.argv)
