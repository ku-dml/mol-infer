Preliminary experiment to select descriptors by using BSP:

Usage:
	python SINGLE_pre.py (.csv) (_values.txt) -rbsp
Sample:
	python SINGLE_pre.py ./instances/At_large_var0_quadratic_h5000_desc_norm.csv ./instances/At_large_norm_values.txt  -rbsp
Once the code finishes, one file containing the selected descriptors will be saved in ./log folder.




Construct prediction functions from the selected descriptors:

Usage:
	python SINGLE_eval.py (.csv) (_values.txt) -rbsp (./log/MLR_based_....csv (the file of the selected descriptors))
Sample:
	python SINGLE_eval.py ./instances/At_large_var0_quadratic_h5000_desc_norm.csv ./instances/At_large_norm_values.txt -rbsp 
	                      ./log/MLR_based_BSP_At_large_var0_quadratic_h5000_desc_norm.csv

