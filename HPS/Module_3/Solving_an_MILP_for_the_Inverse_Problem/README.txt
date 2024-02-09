This is the MILP code used to generate the results in Table 3&4.

The experiment for the paper can be redone by executing:
	python exp_At.py
	python exp_FlmL.py

Before running the code, it is necessary to specify the CPLEX location and compile the feature generator as below:
	g++ -o FV_2LMM_V019 ./2LMM_V019/fv_2LMM.cpp -O2 -Wall -std=c++20
