# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 10:26:17 2023

@author: Muniba
"""
"""
Forward Stepwise Selection with Leave one out Cross-Validation

This script performs forward stepwise selection to select a set of descriptors from a given feature vector file
by using either multi-linear regression or Lasso linear regression as the evaluation function. The feature
selection process is conducted with leave one out cross-validation, and the best descriptors are selected
based on the R2_test.
Note: For leave one out cross-validation R2_test is calculated on overall predicted values for test set

Usage: python forward_selection_cv.py <feature_file> <observed_file> <num_descriptors> [--regression {multi,lasso}]

"""
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Lasso
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut

def forward_leave_one_out_selection(X, y, num_descriptors, regression_type, alpha=0.0001):
    n_features = X.shape[1]
    #print(n_features)
    selected_features = []
    best_model = None
    best_median_score = float('-inf')
    R2_train = None
    R2_test = None
    
    for _ in range(num_descriptors):
        scores = []
        best_r2 = None
        best_fv = None
        
        for feature in range(n_features):
            if feature in selected_features:
                continue

            # Add the new feature to the selected features
            current_features = selected_features + [feature]

            # Perform leave-one-out cross-validation with linear regression or Lasso regression
            if regression_type == "multi":
                model = LinearRegression()
            elif regression_type == "lasso":
                model = Lasso(alpha=alpha)
            else:
                raise ValueError("Invalid regression type. Please specify 'multi' for Multi-Linear Regression or 'lasso' for Lasso Regression.")
            
            loo = LeaveOneOut()
            cv_scores = []
            overall_y_test_pred = []
            overall_y_test_values =[]
            r2_train_values =[]
            for train_idx, test_idx in loo.split(X):
                model.fit(X[train_idx][:, current_features], y[train_idx])
                y_test_pred = model.predict(X[test_idx][:, current_features])
                y_train_pred = model.predict(X[train_idx][:, current_features])
                r2_train = r2_score(y[train_idx], y_train_pred)
                overall_y_test_pred.append(y_test_pred)
                overall_y_test_values.append(y[test_idx])
                r2_train_values.append(r2_train)
                #cv_scores.append(r2_score(y[test_idx], y_pred))
            median_r2_train = np.median(r2_train_values)
            r2_test = r2_score(overall_y_test_values, overall_y_test_pred)
            #median_score = np.median(cv_scores)

            # Store the median R2 score
            scores.append(r2_test)
            if best_r2 is None or r2_test > best_r2:
                best_fv = feature
                best_r2 = r2_test
        #print(scores)    
        # Add the best feature to the selected features
        selected_features.append(best_fv)

        # Update the best model
        if regression_type == "multi":
            best_model = LinearRegression()
        elif regression_type == "lasso":
            best_model = Lasso(alpha=alpha)
        best_model.fit(X[:, selected_features], y)
    #print("R2_test:",best_r2)
    #print("R2_train:",median_r2_train)
    print(f"selected_descriptors:{num_descriptors}   R2_train:{median_r2_train}   R2_test:{best_r2} ")
    return best_model, selected_features

def main(sys):
   # feature_file = args[1]
   # observed_file = args[2]
   # num_descriptors = args[3]
   # regression_type = args[4]

    # Load feature vectors from CSV file
    feature_df = pd.read_csv(feature_file)

    # Load observed values from text file
    observed_df = pd.read_csv(observed_file)

    # Find CIDs in feature vectors not present in observed values
    missing_cids = set(feature_df['CID']) - set(observed_df['CID'])

    if missing_cids:
        print("The following CIDs are present in the feature vectors file but not in the observed values file:")
        print(missing_cids)

    # Merge feature vectors and observed values based on CID
    merged_df = feature_df.merge(observed_df, on='CID')

    # Extract feature vectors and observed values
    X = merged_df.iloc[:, 1:-1].values
    y = merged_df.iloc[:, -1].values

    best_model, selected_features = forward_leave_one_out_selection(X, y, num_descriptors, regression_type)

    #print(selected_features)
    # Save selected descriptors to a CSV file
    output_file = 'selected_descriptors.csv'
    descriptors_df = pd.DataFrame({'Selected Descriptors': selected_features})
    descriptors_df.to_csv(output_file, index=False)
    #print(descriptors_df)
    #print("file of selected desc")
    
    # Create a new feature vector file with only the values of the selected descriptors
    selected_features_indices = np.array(selected_features)
    #print(selected_features_indices)
    #print()
    selected_descriptors = feature_df.columns[1:][selected_features_indices]
    selected_values = merged_df[['CID'] + selected_descriptors.tolist()]

    selected_file_name = f"selected_feature_vectors_{num_descriptors}_{regression_type}.csv"
    selected_values.to_csv(selected_file_name, index=False)
    #print("new feature file generated")

    #print(f"{num_descriptors} {median_r2_train} {best_r2} ")
    #print("Selected features:", selected_features)
   # print("Coefficients:", best_model.coef_)
    #print("Intercept:", best_model.intercept_)


if __name__ == '__main__':
    if len(sys.argv) != 5:
      print("Usage: python forward_selection_cv.py <feature_file> <observed_file> <num_descriptors> [--regression {multi,lasso}]")
    else:
        feature_file = sys.argv[1]
        observed_file = sys.argv[2]
        num_descriptors = int(sys.argv[3])
        regression_type = sys.argv[4]
    #args = (0, "PROTAC_desc_norm.csv", "PROTAC_logS_norm.txt", 4, "multi")
    main(sys.argv)
