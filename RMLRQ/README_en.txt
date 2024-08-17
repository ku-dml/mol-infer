There are four modules of our model, where the Modules 1 and 4 are the same as the previous work. 
For details, please check Modules 1 and 4 at ku-dml/mol-infer/2LMM-LLR/.

This folder contains details related to Modules 2 and 3.

In Module 2, we generate the quadratic descriptors and use R-MLR to obtain prediction functions.  
This folder consists of two sub-folders:
	The folder "generate_quadratic_descritpors" contains the code for generating quadratic descriptors. 
	The folder "construct_prediction_functions" contains the code of using BSP procedure to select descriptors which 
	are used to construct prediction functions based on MLR.

In Module 3, we solve MILP to generate chemical graphs. 
There are two sub-folders of this folder:.
	The folder "MILP" contains the code for generating chemical graphs. 
	Please refer to Module 3 in ku-dml/mol-infer/2LMM-LLR for more details.
	The folder "GNS" contains the code for using grid-neighbor-search to generate more chemical graphs. 
	Please refer to ku-dml/mol-infer/Grid-neighbor-search for more details.
	
Each folder has its own readme, codes, and example instances. 

Reference: "Molecular Design Based on Integer Programming and Quadratic Descriptors in a Two-layered Model" by
Jianshen Zhu, Naveed Ahmed Azam, Shengjuan Cao, Ryota Ido, Kazuya Haraguchi, Liang Zhao, Hiroshi Nagamochi and Tatsuya Akutsu

