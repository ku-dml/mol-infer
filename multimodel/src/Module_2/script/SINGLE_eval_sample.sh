#!/bin/bash

echo 'property_D1	n_D1	K_D1	K_linear_D1	K_quad_D1	Way_D1	R2_train_D1			R2_test_D1 			Time_D1 Total_time'   


python SINGLE_eval.py    ./data_multi/Mp_small_var0_desc_norm.csv                   ./data_multi/Mp_small_norm_values.txt   -o ./output/Mp_small_LLR  -l 0.00019       
python SINGLE_eval.py    ./data_multi/Bp_small_var0_desc_norm.csv                   ./data_multi/Bp_small_norm_values.txt   -o ./output/Bp_small_RF  -rf ./log/RF_Bp_small_var0_desc_norm.csv max_features 1.0 min_samples_split 3 min_samples_leaf 1 max_depth 20 n_estimators 50      
python SINGLE_eval.py    ./data_multi/Fp_small_var0_desc_norm.csv                   ./data_multi/Fp_small_norm_values.txt   -o ./output/Fp_small_ANN  -ann ./log/LLR_ANN_Fp_small_var0_desc_norm.csv 0.89 1523 15 15 15 