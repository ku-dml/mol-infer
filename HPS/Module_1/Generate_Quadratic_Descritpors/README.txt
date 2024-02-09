This is the code used to generate the feature vector file with quadratic terms used for the learning procedure of RMLR.

Usage:
	python generate_quadratic_descriptors.py (.csv) (_values.txt) (size of obtained quadratic descriptors)

Sample:
	Python generate_quadratic_descriptors.py ./instances/At_large_var0_desc_norm.csv ./instances/At_large_norm_values.txt 5000
The output file will be ./instances/At_large_var0_quadratic_h5000_desc_norm.csv
