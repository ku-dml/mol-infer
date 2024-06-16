Preliminary experiment:

Usage:
	python SINGLE_pre.py (.csv) (_values.txt) -l(lasso)/-ann(ANN)
Sample:
	python SINGLE_pre.py ./Datasets/Wang_desc_norm.csv ./Datasets/Wang_norm_values.txt  -ann


The output will be something like:
Wang	1414	405	405	0	ANN	0.9900755018624954	0.8400401975314331	6435.184849977493		-ann ./log/LLR_ANN_Wang_desc_norm.csv 0.99 7389 96 96 	

Please copy the last term (starting from -ann, or -l if using Lasso), and use the following command to run the evaluation experiment.

Usage:
	python SINGLE_eval.py (.csv) (_values.txt) (-..., the last term in the result of preliminary experiment)
Sample:
	python SINGLE_eval.py ./Datasets/Wang_desc_norm.csv ./Datasets/Wang_norm_values.txt -ann ./log/LLR_ANN_Wang_desc_norm.csv 0.99 7389 96 96

