Please check Github for the code to generate linear descriptors.

To run evaluation experiments (there is no preliminary experiment for MLR) with MLR:

Usage:
	python SINGLE_eval.py (.csv) (_values.txt) (-o (OUTPUT_PREFIX)) -mlr
Sample:
	python SINGLE_eval.py ./Datasets/D5_desc_norm.csv ./Datasets/D5_norm_values.txt -o ./output/At_small_ANN -mlr
