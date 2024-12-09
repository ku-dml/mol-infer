#!/bin/bash

echo 'property_D1	n_D1	K_D1	K_linear_D1	K_quad_D1	Way_D1	R2_train_D1			R2_test_D1 			Time_D1 Total_time'   

props=(Bp Fp Kow Lp Mp Sl)

Module_1_output=../Module_1/sample_instance/output

for prop in ${props[@]}; do
  # LR
  echo LR $prop start
  python3 SINGLE_eval.py   $Module_1_output/${prop}_small_desc_norm.csv sample_instance/input/${prop}_small_norm_values.txt  -o  sample_instance/output/${prop}_small_LLR  -l 0.00019

  # RF
  echo RF $prop start
  python3 SINGLE_eval.py   $Module_1_output/${prop}_small_desc_norm.csv sample_instance/input/${prop}_small_norm_values.txt  -o  sample_instance/output/${prop}_small_RF  -rf ./log/RF_${prop}_small_var0_desc_norm.csv max_features 1.0 min_samples_split 3 min_samples_leaf 1 max_depth 20 n_estimators 50

  # ANN
  echo ANN $prop start
  python3 SINGLE_eval.py   $Module_1_output/${prop}_small_desc_norm.csv sample_instance/input/${prop}_small_norm_values.txt  -o  sample_instance/output/${prop}_small_ANN  -ann ./log/LLR_ANN_${prop}_small_var0_desc_norm.csv 0.89 1523 15 15 15
done



   
					