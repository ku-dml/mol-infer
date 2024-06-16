Preliminary experiment:

Usage:
	python SINGLE_pre.py (.csv) (_values.txt) -l(lasso)/-ann(ANN)
Sample:
	python SINGLE_pre.py ./Datasets/D5_desc_norm.csv ./Datasets/D5_norm_values.txt  -ann


The output will be something like:
D5	91	118	118	0	ANN	0.9306658282010498	0.7347627803279586	4768.666820526123		-ann ./log/LLR_ANN_D5_desc_norm.csv 0.93 179 22 14 	

Please copy the last term (starting from -ann, or -l if using Lasso), and use the following command to run the evaluation experiment.

Usage:
	python SINGLE_eval.py (.csv) (_values.txt) (-..., the last term in the result of preliminary experiment)
Sample:
	python SINGLE_eval.py ./Datasets/D5_desc_norm.csv ./Datasets/D5_norm_values.txt -ann ./log/LLR_ANN_D5_desc_norm.csv 0.93 179 22 14

