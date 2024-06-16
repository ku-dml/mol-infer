
import sys
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def main(desc_data_filename, target_values_filename):
    # Load normalized desc data from CSV and target values from TXT
    desc_data = pd.read_csv(desc_data_filename)
    target_values = np.loadtxt(target_values_filename, skiprows=1, delimiter=',')
    cids = desc_data['CID']  # Assuming 'CID' is the column name for compound IDs

    num_descriptors = len(desc_data.columns) - 1  # Exclude the 'CID' column
    
    sp_loc = desc_data_filename.rfind('/')
    dataset_name = desc_data_filename[sp_loc+1:].split('_')[0]
    
    loo = LeaveOneOut()

    results = []  # List to store R-squared values and CIDs
    test_predictions = []  # List to store test set predictions
    test_targets = []  # List to store test set target values

    train_r2_values = []  # List to store train set R-squared values

    for train_index, test_index in loo.split(desc_data):
        X_train, X_test = desc_data.iloc[train_index], desc_data.iloc[test_index]
        y_train, y_test = target_values[train_index], target_values[test_index]

        model = LinearRegression()
        model.fit(X_train, y_train)

        y_train_pred = model.predict(X_train)
        train_r2 = r2_score(y_train, y_train_pred)
        train_r2_values.append(train_r2)

        y_test_pred = model.predict(X_test)

        # Append results to the list
        results.append({'CID': cids.iloc[test_index[0]], 'Train_R2': train_r2})

        # Save test set predictions and target values
        test_predictions.append(y_test_pred[0])
        test_targets.append(y_test[0])

    # Calculate average R-squared values
    max_train_r2 = np.max(train_r2_values)
    average_train_r2 = np.mean(train_r2_values)
    min_train_r2 = np.min(train_r2_values)

    print(f"{dataset_name}  {num_descriptors}  {average_train_r2:.4f}  {max_train_r2:.4f}  {min_train_r2:.4f}  {r2_score(test_targets, test_predictions):.4f}")

    # Create a DataFrame from the results list
    results_df = pd.DataFrame(results)

    # Ensure the log directory exists
    log_dir = 'log'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # Save results to a new CSV file in the log directory
    results_csv_path = os.path.join(log_dir, f"{dataset_name}_r2_results.csv")
    results_df.to_csv(results_csv_path, index=False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script_name.py desc_data.csv target_values.txt")
    else:
        desc_data_filename = sys.argv[1]
        target_values_filename = sys.argv[2]
        main(desc_data_filename, target_values_filename)
