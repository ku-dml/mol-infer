#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mol-infer_ANN.py

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

This is an improved version of scikit_chemgraph_learning.py
which is the first ANN learner in our mol-infer project. 
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
import time
import sys

MAX_ITER = 10000000000
EARLY_STOPPING = False

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
    with open(weights_filename, 'w', newline='\n') as f:
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
    with open(biases_filename, 'w', newline='\n') as f:
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

    The set of CIDs in descriptor_filename should be
    a subset of that in target_values_filename.
    """
    
    # read the training and target data
    fv = pd.read_csv(descriptors_filename) 
    print('{} contains {} vectors for {} (=CID+{}) features.'.format(descriptors_filename, fv.shape[0], fv.shape[1], fv.shape[1]-1))
    value = pd.read_csv(target_values_filename)       
    print('{} contains {} target values.'.format(target_values_filename, value.shape[0]))

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
        
    # CID check: target_values_filename should contain all cids that appear in descriptors_filename
    for cid in CIDs:
        if cid not in target_dict:
            sys.stderr.write('error: {} misses the target value of CID {}\n'.format(target_values_filename, cid))
            exit(1)
            
    # construct x and y so that the CIDs are ordered in ascending order
    CIDs.sort()
    x = np.array([fv_dict[cid] for cid in CIDs])
    y = np.array([target_dict[cid] for cid in CIDs])
    
    print('n range = [{},{}]'.format(fv['n'].min(),fv['n'].max()))
    print('a range = [{},{}]'.format(min(y), max(y)))
    
    numdata = x.shape[0]
    numfeature = x.shape[1]
    print('#instances = {}'.format(numdata))
    print('#features = {}'.format(numfeature))
    
    # Initialize an artificial neural network - MLP regressor
    reg = MLPRegressor(activation='relu', solver='adam',
                       alpha=1e-5, hidden_layer_sizes=architecture,
                       random_state=1, max_iter=MAX_ITER,
                       early_stopping=EARLY_STOPPING)
    
    score_R2 = np.array([])
    score_MAE = np.array([])
    time_train = np.array([])
    
    # Separate the data randomly for cross-validation
    kf = KFold(n_splits=5, shuffle=True, random_state=4)
    fold = 0
    for train, test in kf.split(x):
        fold += 1
        x_train, x_test, y_train, y_test = x[train], x[test], y[train], y[test]
        print('\nD{}: train {}, test {}'.format(fold, 
              x_train.shape[0], x_test.shape[0]))
    
        start = time.time()
        reg.fit(x_train, y_train)
        end = time.time()
        timetemp = end-start
        
        print('training time: {}'.format(timetemp))
        time_train = np.append(time_train, timetemp)
    
        pred = reg.predict(x)
        pred_train = reg.predict(x_train)
        pred_test = reg.predict(x_test)
           
        # calculate the prediction score (R^2)
        R2train = reg.score(x_train,y_train)
        R2test = reg.score(x_test,y_test)
        R2all = reg.score(x,y) 
        print('R2 score train = {}'.format(R2train))
        print('R2 score test = {}'.format(R2test))
        print('R2 score all = {}'.format(R2all))
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
    
        # calculate mean absolute error (MAE)
        MAEtrain = mean_absolute_error(y_train,pred_train)
        MAEtest = mean_absolute_error(y_test,pred_test)
        MAEall = mean_absolute_error(y,pred) 
        print('MAE score train = {}'.format(MAEtrain))
        print('MAE score test = {}'.format(MAEtest))
        print('MAE score all = {}'.format(MAEall))
        temp = np.array([MAEtrain, MAEtest, MAEall]).reshape(1,3)
        score_MAE = np.append(score_MAE, temp)
        
    score_R2 = score_R2.reshape(5,3)
    score_MAE = score_MAE.reshape(5,3)
    time_train = time_train[:,np.newaxis]
    to_print = np.hstack((score_R2, time_train))
    for i in range(5):
        print(' '.join(str(j) for j in to_print[i].tolist()))
    avg_time = np.mean(time_train)
    print('Average time = {}'.format(avg_time))
    avg_testR2 = np.mean(score_R2, 0)[1]
    print('Average R2 test score = {}'.format(avg_testR2))
    avg_testMAE = np.mean(score_MAE, 0)[1]
    print('Average MAE test score = {}'.format(avg_testMAE))
    
    return best_regressor


def main(argv):
    if (len(argv) < 5):
        print("""
usage: {} (TrSet_FeatureVector.csv)(TrSet_TargetValue.csv)(Output_Name)(NN_Architecture)
                  
  - TrSet_FeatureVector.csv ... The CSV file generated by our feature vector generator in Module 1.
  - TrSet_TargetValue.csv ... The CSV file that contains target values.
    = The first column must contain CIDs.
    = The second column must contain target values and be named \'a\'. 
    = All CIDs appearing in TrSet_FeatureVector.csv must be contained in TrSet_TargetValue.csv. 
  - OutputName_{{biases,weights}}.txt will be output. They contain data on the constructed ANN. 
  - NN_Architecture ... Number of nodes in hydden layer(s).
        """.format(sys.argv[0]))
        sys.exit()

    # Parse the command line arguments
    descriptors_filename = argv[1]
    target_values_filename = argv[2]
    output_filename = argv[3]
    architecture = tuple(int(a) for a in argv[4:])
    
    # Perform 5-fold validation with the given training data
    # and return the regressor that achieves highest R^2 test score
    best_regressor = train_ANN(descriptors_filename,
                               target_values_filename,
                               architecture)
    
    # Write the weights and biases of the regressor to files
    write_weights_biases(best_regressor, output_filename)
    
    
main(sys.argv)
    
    


        
      
