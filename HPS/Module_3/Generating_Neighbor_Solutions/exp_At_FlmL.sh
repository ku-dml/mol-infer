#!/bin/sh

## For experiments for GNS of At and FlmL

echo ' infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/At_c_370_380 ./p_max/p_max_delta_r_0.05.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/At_c_370_380 ./p_max/p_max_delta_r_0.05.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_0.05.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_0.05.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_0.05.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_0.05.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_0.05.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.05_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_0.05.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl

echo ' infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/At_c_370_380 ./p_max/p_max_delta_r_0.15.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/At_c_370_380 ./p_max/p_max_delta_r_0.15.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_0.15.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_0.15.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_0.15.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_0.15.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_0.15.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_0.15_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_0.15.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl

echo ' infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/At_c_370_380 ./p_max/p_max_delta_r_1e-1.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/At_c_370_380 ./p_max/p_max_delta_r_1e-1.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_1e-1.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_1e-1.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_1e-1.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_1e-1.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_1e-1.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-1_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_1e-1.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl

echo ' infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/At_c_370_380 ./p_max/p_max_delta_r_1e-2.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/At/At 370 380 ./MILP_result/At/instance_c_2LMM.txt ./MILP_result/At/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/At_c_370_380 ./p_max/p_max_delta_r_1e-2.txt -d1 ANN -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_1e-2.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_a_2LMM.txt ./MILP_result/FlmL/ins_a_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/FlmL_a_-0.55_-0.50 ./p_max/p_max_delta_r_1e-2.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_1e-2.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.45 -0.40 ./MILP_result/FlmL/instance_b4_2LMM.txt ./MILP_result/FlmL/ins_b4_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/FlmL_b4_-0.45_-0.40 ./p_max/p_max_delta_r_1e-2.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl
echo ' infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_1e-2.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl'
python infer_2LMM_SEP_GNS.py ./MILP_result/FlmL/FlmL -0.55 -0.50 ./MILP_result/FlmL/instance_c_2LMM.txt ./MILP_result/FlmL/ins_c_fringe_2LMM.txt ./MILP_result/sdf/delta_e-2_p6/FlmL_c_-0.55_-0.50 ./p_max/p_max_delta_r_1e-2.txt -d1 R-MLR -d2 R-MLR -tp 6 -t ./Mp/Mp ./Sl/Sl