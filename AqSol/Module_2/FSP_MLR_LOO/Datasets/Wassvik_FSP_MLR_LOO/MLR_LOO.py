
import sys
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
    average_train_r2 = np.mean(train_r2_values)

    print(f"Number of Descriptors: {num_descriptors} Average R-squared (Training Set): {average_train_r2:.4f} Overall R-squared for Test Set: {r2_score(test_targets, test_predictions):.4f}")

    # Create a DataFrame from the results list
    results_df = pd.DataFrame(results)

    # Save results to a new CSV file
    results_df.to_csv('r2_results.csv', index=False)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script_name.py desc_data.csv target_values.txt")
    else:
        desc_data_filename = sys.argv[1]
        target_values_filename = sys.argv[2]
        main(desc_data_filename, target_values_filename)
