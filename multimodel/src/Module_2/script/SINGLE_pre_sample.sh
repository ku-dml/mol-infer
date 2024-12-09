#!/bin/bash

echo 'property_D1	n_D1	K_D1	K_linear_D1	K_quad_D1	Way_D1	R2_train_D1	R2_test_D1 	Time_D1 		C_D1	C_D2'

props=(Bp Fp Kow Lp Mp Sl)

Module_1_output=../Module_1/sample_instance/output

for prop in ${props[@]}; do
  # LR
  echo LR $prop start
  python3 SINGLE_pre.py   $Module_1_output/${prop}_small_desc_norm.csv sample_instance/input/${prop}_small_norm_values.txt  -l
  # RF
  echo RF $prop start
  python3 SINGLE_pre.py   $Module_1_output/${prop}_small_desc_norm.csv sample_instance/input/${prop}_small_norm_values.txt  -rf
  # ANN
  echo ANN $prop start
  python3 SINGLE_pre.py   $Module_1_output/${prop}_small_desc_norm.csv sample_instance/input/${prop}_small_norm_values.txt  -ann
done



   
					