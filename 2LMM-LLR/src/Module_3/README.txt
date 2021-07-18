For conducting the experiment in the paper, 
please uncomment the part in infer_2LMM_L.py from line 352 and line 373.

One sample command to run the code:
python infer_2LMM_L.py ./sample_instance/Hc 1900 1920 ./sample_instance/instance_b4_test_2LMM.txt 
		./sample_instance/ins_b4_test_fringe_2LMM.txt test

python infer_2LMM_L.py (prefix for the property) (lower bound for predicted value) 
		(upper bound for predicted value) (instance file) (fringe tree file) (prefix for output)
