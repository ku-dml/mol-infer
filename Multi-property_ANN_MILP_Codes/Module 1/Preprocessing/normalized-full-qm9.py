# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 12:23:28 2025

@author: Dell Vostro
"""

import os
import sys
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

def minmax_normalize_qm9(input_csv, output_txt):
    """
    Normalize all numeric QM9 property columns using Min–Max scaling.
    Save the result as a comma-separated TXT file.
    """

    print("Loading QM9 file:", input_csv)
    df = pd.read_csv(input_csv)

    # Identify numeric columns (properties)
    numeric_cols = df.select_dtypes(include=['int64', 'float64']).columns.tolist()

    print("Numeric columns found for scaling:")
    for col in numeric_cols:
        print("  -", col)

    # Apply Min–Max normalization
    scaler = MinMaxScaler()
    df[numeric_cols] = scaler.fit_transform(df[numeric_cols])

    # Save as TXT but comma-separated
    df.to_csv(output_txt, index=False, sep=",")

    print("\nNormalized TXT file saved as:", output_txt)


# -------------------------------------------------------
# Run the function
# -------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        print("Usage: python normalized-full-qm9.py input.csv [output.txt]")
        sys.exit(1)

    input_file = sys.argv[1]                                  # <-- input CSV
    if len(sys.argv) == 3:
        output_file = sys.argv[2]                             # <-- output TXT
    else:
        output_file = os.path.splitext(input_file)[0] + "_norm_values.txt"
    minmax_normalize_qm9(input_file, output_file)
