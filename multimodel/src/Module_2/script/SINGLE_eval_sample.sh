#!/bin/bash

echo 'property_D1	n_D1	K_D1	K_linear_D1	K_quad_D1	Way_D1	R2_train_D1			R2_test_D1 			Time_D1 Total_time'   

# example for Bp

Module_1_output=../Module_1/sample_instance/output

python SINGLE_eval.py $Module_1_output/Bp_small_desc_norm.csv sample_instance/input/Bp_small_norm_values.txt -o ./output/Bp_small_ANN -ann ./log/LLR_ANN_Bp_small_desc_norm.csv 0.93 424 24 18 13

python SINGLE_eval.py $Module_1_output/Bp_small_desc_norm.csv sample_instance/input/Bp_small_norm_values.txt -o ./output/Bp_small_LLR -l 0.0001	

python SINGLE_eval.py $Module_1_output/Bp_small_desc_norm.csv sample_instance/input/Bp_small_norm_values.txt -o ./output/Bp_small_RF -rf ./log/RF_Bp_small_desc_norm.csv max_features 1.0 min_samples_split 3 min_samples_leaf 1 max_depth 12 n_estimators 100
