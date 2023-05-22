For conducting the experiment in the paper, please uncomment the part in infer_2LMM_L.py from line 352 and line 373.

Usage:
	python infer_2LMM_L_q.py (prefix for the property) (lower bound for predicted value		upper bound for predicted value) (instance file.txt) (fringe tree file.txt) (prefix for output)
Sample:
	python infer_2LMM_L_q.py ./instances/Bp/Bp 225 235 ./instances/Bp/instance_a_2LMM.txt 
		./instances/Bp/ins_a_fringe_2LMM.txt test_Bp_a_225_235

Please refer to Module 3 in ku-dml/mol-infer/2LMM-LLR for more details.
