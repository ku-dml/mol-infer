#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Copyright (C) 2020
by Discrete Mathematics Lab, 
   Department of Applied Mathematics and Physics, 
   Graduate School of Informatics,
   Kyoto University

Licensed under the MIT license, see License.txt
"""

"""
scikit_chemgraph_learning.py

This file implements functions that given a file with
descriptor values of chemical compounds and a file with target values,
performs 5-fold learning using an artificial neural network (MLP regressor)
and stores the weights and biases of the training iteration that achieved
the highest R^2 test score over the 5 trials. 

Command line arguments include
 - filename with a list of descriptors
 - filename with a list of target values
 - filename for the output weights and biases
 - network architecture, as a list of hidden layer sizes
 
The program will output on the terminal the R^2 and MAE scores, and
time taken for training for each 
trial in the 5-fold cross-validation, 
and the averages over the 5 trials at the end.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
import time
import sys

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
    

def train_ANN(descriptors_filename, target_values_filename, architecture):
    """
    Given filenames of a file containing a list of descriptors
    and a file containing target values, and a tuple of integers
    giving the number of nodes in hidden layers of an ANN,
    perform 5-fold cross-validation learning by using the 
    descriptors and target values with an ANN of the given architectire,
    and return the trained ANN (MLP regressor) that achieved the highest
    R^2 test score over the test data.
    """
    
    # read the training and target data
    fv = pd.read_csv(descriptors_filename) 
    value = pd.read_csv(target_values_filename) 
    
    cids = set(fv["CID"])
    
    # print(value)
    
    to_drop = list()
    for i, cc in enumerate(value["CID"]):
        # print (i, cc)
        if cc not in cids:
            to_drop.append(i)
    value.drop(to_drop, inplace=True)
    
    to_drop = list()
    cids = set(value["CID"])
    for i, cc in enumerate(fv["CID"]):
        if cc not in cids:
            to_drop.append(i)
    fv.drop(to_drop, inplace = True)
            
            
    # print("#"*40)
    # print(value)
    
    # prepare target, train, test arrays
    target = np.array(value['a'])
    
    x = fv.values[:,1:]
    y = target
    
    print('n range = [{},{}]'.format(fv['n'].min(),fv['n'].max()))
    print('a range = [{},{}]'.format(value['a'].min(),value['a'].max()))
    
    numdata = x.shape[0]
    numfeature = x.shape[1]
    print('#instances = {}, {}'.format(numdata, y.shape[0]))
    print('#features = {}'.format(numfeature))
    
    # Initialize an artificial neural network - MLP regressor
    reg = MLPRegressor(activation='relu', solver='adam',
                       alpha=1e-5, hidden_layer_sizes=architecture,
                       random_state=1, max_iter=100000000)
    
    score_R2 = np.array([])
    score_MAE = np.array([])
    time_train = np.array([])
    fold_info = np.array([])
    
    # Separate the data randomly for cross-validation
    kf = KFold(n_splits=5, shuffle=True, random_state=21)
    fold = 0
    for train, test in kf.split(x):
        fold += 1
        x_train, x_test, y_train, y_test = x[train], x[test], y[train], y[test]
        fold_info = np.append(fold_info, 
                              'D{}({},{})'.format(fold, 
                                        x_train.shape[0], x_test.shape[0]))
    
        start = time.time()
        reg.fit(x_train, y_train)
        end = time.time()
        timetemp = end-start
        
        # print('training time: {}'.format(timetemp))
        time_train = np.append(time_train, timetemp)
    
        pred = reg.predict(x)
        pred_train = reg.predict(x_train)
        pred_test = reg.predict(x_test)
           
        # calculate the prediction score (R^2)
        R2train = reg.score(x_train,y_train)
        R2test = reg.score(x_test,y_test)
        R2all = reg.score(x,y) 
        # print('R2 score train = {}'.format(R2train))
        # print('R2 score test = {}'.format(R2test))
        # print('R2 score all = {}'.format(R2all))
        temp = np.array([R2train, R2test, R2all]).reshape(1,3)
        score_R2 = np.append(score_R2, temp)
        
        # check the test R2 score and store the regressor with the highest one
        if (fold == 1):
            best_regressor = reg
            best_R2_score = R2test
        else:
            if (R2test > best_R2_score):
                best_regressor = reg
                best_R2_score = R2test
    
    score_R2 = score_R2.reshape(5, 3)
    fold_info = fold_info.reshape(5, 1)
    time_train = time_train[:,np.newaxis]
    to_print = np.hstack((fold_info, score_R2, time_train))
    for i in range(5):
        print(f"{architecture}; ", 
              '; '.join(str(j) for j in to_print[i].tolist()), 
              "; ",  flush = True)
    avg_time = np.mean(time_train)
    tot_time = np.sum(time_train)
    avg_testR2 = np.mean(score_R2, 0)
    print(f"{architecture}; ",
          "Average; {}; {}; {}".format("; ".join([str(aa) 
                                                  for aa in avg_testR2]), 
                                       avg_time, 
                                       tot_time),  flush = True)
    
    
    return best_regressor


def main(argv):
    if (len(argv) < 4):
        print("""
              Please supply at least 3 command line arguments:
                  - Descriptors as training data
                  - Target values as training data
                  - Filename for the output weights/biases files
                  
              The program will now terminate.""")
        sys.exit()
    # else:
    # Parse the command line arguments
    descriptors_filename = argv[1]
    print(descriptors_filename)
    target_values_filename = argv[2]
    output_filename = argv[3]
    architecture = tuple(int(a) for a in argv[4:])
    
    # Perform 5-fold validation with the given training data
    # and return the regressor that achieves highest R^2 test score
    
    print(architecture, time.strftime("%X"), flush = True)
    best_regressor = train_ANN(descriptors_filename,
                               target_values_filename,
                               architecture)
        
    # Write the weights and biases of the regressor to files
    write_weights_biases(best_regressor, 
                         output_filename)
    
main(sys.argv)
