# Preliminary experiment:

Usage:

	python SINGLE_pre.py (.csv) (_values.txt) -l(lasso)/-ann(ANN)/-rf(RF)
Sample:

	python SINGLE_pre.py ../Module_1/sample_instance/output/Bp_small_desc_norm.csv sample_instance/input/Bp_small_norm_values.txt -l


The output will be something like:

```
Bp_small	400	15	15	0	LASSO	0.549641087647406	0.41578203853379403	33.90013933181763 -l ./log/LLR_LASSO_Bp_small_desc_norm.csv 0.54 258 30 20
```

Please copy the last term (starting from -ann, or -l if using Lasso), and use the following command to run the evaluation experiment.

# Evaluation experiment:
Usage:
  
	python SINGLE_eval.py (.csv) (_values.txt) (-o (OUTPUT_PREFIX)) (-..., the last term in the result of preliminary experiment)
Sample:

	Bp_small	400	15	15	0	LASSO	0.549641087647406	0.41578203853379403	33.90013933181763		-l ./log/LLR_LASSO_Bp_small_desc_norm.csv 0.54 258 30 20

