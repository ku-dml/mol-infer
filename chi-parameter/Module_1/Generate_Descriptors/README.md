## Generate Descriptors

This folder contains the code to generate the linear descriptors of each kinds of molecules, defined in Section 2.3, which will be then combined into one feature vector in [Combine Descriptors](/chi-parameter/Module_1/Combine_Descriptors) and used for the learning stage in [Module 2](/chi-parameter/Module_2).

Usage:

First compile the `c++` files as following on terminal:

```
g++ -o fv_2LMM_polymer fv_2LMM_polymer.cpp -O2 -Wall -std=c++20
```

This will generate an executable file named `fv_2LMM_polymer`,
and then generate the linear feature descriptor files by using the command:

```
./fv_2LMM_polymer DATASET.sdf OUTPUT
```

Here:
- DATASET.sdf: the input sdf file containing the **monomer-representation** of the polymer data set, and
- OUTPUT: the prefix for the output files.
  
When the generator finishes normally, it will generate the following files:
- OUTPUT_desc.csv: The csv file of the linear descriptors of the dataset **BEFORE** normalization.
- OUTPUT_desc_norm.csv: The csv file of the linear descriptors of the dataset **AFTER** normalization.
- OUTPUT_fringe.txt: The file containing the _fringe tree_ information of the dataset, which will be used in the stage of solving inverse problem (Module 3).

A sample usage:

```
./fv_2LMM_polymer ./sample_instance/sample1.sdf ./sample_instance/output
```

