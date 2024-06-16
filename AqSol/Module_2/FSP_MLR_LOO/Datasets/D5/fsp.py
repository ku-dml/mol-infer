# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 02:53:37 2023

@author: Naveed

Forward Stepwise Selection with Cross-Validation

This script performs forward stepwise selection to select a set of descriptors from a given feature vector file
by using either multi-linear regression or Lasso linear regression as the evaluation function. The feature
selection process is conducted 1 time with 5-fold cross-validation, and the best descriptors are selected
based on the median R2 score.

Usage: python forward_selection_cv.py <feature_file> <observed_file> <num_descriptors> [--regression {multi,lasso}]

Arguments:
  feature_file       Path to the feature vector file in CSV format
  observed_file      Path to the observed values file in text format
  num_descriptors    Number of descriptors to select
  --regression       Regression type to use for evaluation (multi - Multi-linear Regression, lasso - Lasso Linear Regression)
 The default vlaue of alpha for lassor is 0.001
Example: python forward_selection_cv.py Klamt_D_desc_norm.csv klamt_D_logS_norm.txt 5 --regression lasso
"""
import argparse
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Lasso
from sklearn.metrics import r2_score
from sklearn.model_selection import cross_val_score, KFold

def forward_stepwise_selection(X, y, num_descriptors, regression_type, alpha=0.0001):
    n_features = X.shape[1]
    print(X.shape[1])
    
    selected_features = []
    best_model = None
    best_median_score = float('-inf')

    for _ in range(num_descriptors):
        
        scores = []
        best_r2 = None
        best_fv = None

        for feature in range(n_features):
            #print(feature)
            #print(range(n_features))
            #print(n_features)
            if feature in selected_features:
                #print(feature)
                continue

            # Add the new feature to the selected features
            current_features = selected_features + [feature]

            # Perform 5-fold cross-validation with linear regression or Lasso regression
            if regression_type == "multi":
                model = LinearRegression()
            elif regression_type == "lasso":
                model = Lasso(alpha=alpha)
            else:
                raise ValueError("Invalid regression type. Please specify 'multi' for Multi-Linear Regression or 'lasso' for Lasso Regression.")

            cv_scores = cross_val_score(model, X[:, current_features], y, cv=KFold(n_splits=2, shuffle=True, random_state=0))
            median_score = np.median(cv_scores)

            # Store the median R2 score
            scores.append(median_score)
            if best_r2 is None or median_score > best_r2:
                best_fv = feature
                best_r2 = median_score

        # Find the feature that maximizes the median R2 score
        print("SCORES:",scores)
        #print(len(scores))
        #best_feature = np.argmax(scores)
        #print(best_feature)
        #print(scores[24])
        #best_median_score = scores[best_feature]
        

        # Add the best feature to the selected features
        selected_features.append(best_fv)
        #print(selected_features)

        # Update the best model
        if regression_type == "multi":
            best_model = LinearRegression()
        elif regression_type == "lasso":
            best_model = Lasso(alpha=alpha)
        best_model.fit(X[:, selected_features], y)

    return best_model, selected_features

def main(args):
    feature_file = args[1]
    observed_file = args[2]
    num_descriptors = args[3]
    regression_type = args[4]

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

    best_model, selected_features = forward_stepwise_selection(X, y, num_descriptors, regression_type)

    print(selected_features)
    # Save selected descriptors to a CSV file
    output_file = 'selected_descriptors.csv'
    descriptors_df = pd.DataFrame({'Selected Descriptors': selected_features})
    descriptors_df.to_csv(output_file, index=False)
    #print(descriptors_df)
    print("file of selected desc")
    
    # Create a new feature vector file with only the values of the selected descriptors
    selected_features_indices = np.array(selected_features)
    print(selected_features_indices)
    print()
    selected_descriptors = feature_df.columns[1:][selected_features_indices]
    selected_values = merged_df[['CID'] + selected_descriptors.tolist()]

    selected_file_name = f"selected_feature_vectors_{num_descriptors}_{regression_type}.csv"
    selected_values.to_csv(selected_file_name, index=False)
    print("new feature file generated")

    #print("Selected features:", selected_features)
   # print("Coefficients:", best_model.coef_)
    #print("Intercept:", best_model.intercept_)


if __name__ == '__main__':
    args = (0, "Boobier2_desc_norm.csv","Boobier2_logS_norm.txt", 80, "multi")
    main(args)