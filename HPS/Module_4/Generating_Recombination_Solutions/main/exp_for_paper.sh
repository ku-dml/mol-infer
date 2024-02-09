#!/bin/sh

echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Next property At' 


echo '~~a'
			./generate_isomers.o ../../results/sdf/At_a_240_250.sdf 10 1000000 0 300 1000000 100 Results_At_a_240_250.sdf ../../results/sdf/At_a_240_250_partition.txt ../../results/At/ins_a_fringe_2LMM.txt  
echo '~~b1'                                                                                                                                                                          
			./generate_isomers.o ../../results/sdf/At_b1_240_250.sdf 10 1000000 0 300 1000000 100 Results_At_b1_240_250.sdf ../../results/sdf/At_b1_240_250_partition.txt ../../results/At/ins_b1_fringe_2LMM.txt  
echo '~~b2'                                                                                                                                                                          
			./generate_isomers.o ../../results/sdf/At_b2_240_250.sdf 10 1000000 0 300 1000000 100 Results_At_b2_240_250.sdf ../../results/sdf/At_b2_240_250_partition.txt ../../results/At/ins_b2_fringe_2LMM.txt  
echo '~~b3'                                                                                                                                                                          
			./generate_isomers.o ../../results/sdf/At_b3_190_200.sdf 10 1000000 0 300 1000000 100 Results_At_b3_190_200.sdf ../../results/sdf/At_b3_190_200_partition.txt ../../results/At/ins_b3_fringe_2LMM.txt  
echo '~~b4'                                                                                                                                                                          
			./generate_isomers.o ../../results/sdf/At_b4_300_310.sdf  10 1000000 0 300 1000000 100 Results_At_b4_300_310.sdf  ../../results/sdf/At_b4_300_310_partition.txt  ../../results/At/ins_b4_fringe_2LMM.txt
echo '~~c'                                                                                                                                                                           
			./generate_isomers.o ../../results/sdf/At_c_360_370.sdf  10 1000000 0 300 1000000 100 Results_At_c_360_370.sdf  ../../results/sdf/At_c_360_370_partition.txt  ../../results/At/ins_c_fringe_2LMM.txt
echo '~~d'                                                                                                                                                                           
			./generate_isomers.o ../../results/sdf/At_d_230_240.sdf  10 1000000 0 300 1000000 100 Results_At_d_230_240.sdf  ../../results/sdf/At_d_230_240_partition.txt  ../../results/At/ins_d_fringe_2LMM.txt



echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Next property FlmL' 


echo '~~a'
			./generate_isomers.o ../../results/sdf/FlmL_a_-0.55_-0.50.sdf 10 1000000 0 300 1000000 100 Results_FlmL_a_-0.55_-0.50.sdf ../../results/sdf/FlmL_a_-0.55_-0.50_partition.txt ../../results/FlmL/ins_a_fringe_2LMM.txt  
echo '~~b1'                                                                                                                                                                          
			./generate_isomers.o ../../results/sdf/FlmL_b1_-0.15_-0.10.sdf 10 1000000 0 300 1000000 100 Results_FlmL_b1_-0.15_-0.10.sdf ../../results/sdf/FlmL_b1_-0.15_-0.10_partition.txt ../../results/FlmL/ins_b1_fringe_2LMM.txt  
echo '~~b2'                                                                                                                                                                          
			./generate_isomers.o ../../results/sdf/FlmL_b2_-0.50_-0.45.sdf 10 1000000 0 300 1000000 100 Results_FlmL_b2_-0.50_-0.45.sdf ../../results/sdf/FlmL_b2_-0.50_-0.45_partition.txt ../../results/FlmL/ins_b2_fringe_2LMM.txt  
echo '~~b3'                                                                                                                                                                          
			./generate_isomers.o ../../results/sdf/FlmL_b3_-0.50_-0.45.sdf 10 1000000 0 300 1000000 100 Results_FlmL_b3_-0.50_-0.45.sdf ../../results/sdf/FlmL_b3_-0.50_-0.45_partition.txt ../../results/FlmL/ins_b3_fringe_2LMM.txt  
echo '~~b4'                                                                                                                                                                         
			./generate_isomers.o ../../results/sdf/FlmL_b4_-0.45_-0.40.sdf  10 1000000 0 300 1000000 100 Results_FlmL_b4_-0.45_-0.40.sdf  ../../results/sdf/FlmL_b4_-0.45_-0.40_partition.txt  ../../results/FlmL/ins_b4_fringe_2LMM.txt
echo '~~c'                                                                                                                                                                           
			./generate_isomers.o ../../results/sdf/FlmL_c_-0.55_-0.50.sdf  10 1000000 0 300 1000000 100 Results_FlmL_c_-0.55_-0.50.sdf  ../../results/sdf/FlmL_c_-0.55_-0.50_partition.txt  ../../results/FlmL/ins_c_fringe_2LMM.txt
echo '~~d'                                                                                                                                                                           
			./generate_isomers.o ../../results/sdf/FlmL_d_0.20_0.25.sdf  10 1000000 0 300 1000000 100 Results_FlmL_d_0.20_0.25.sdf  ../../results/sdf/FlmL_d_0.20_0.25_partition.txt  ../../results/FlmL/ins_d_fringe_2LMM.txt


