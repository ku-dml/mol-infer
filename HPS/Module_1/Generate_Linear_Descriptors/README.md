## Generate Linear Descriptors

This folder contains the code to generate the linear descriptors defined in Section 3, paragraph "Feature Function", which will be used for the learning stage in [Module 2](/HPS/Module_2).

Usage:

First compile the `c++` files as following on terminal:

```
g++ -o fv_2LMM fv_2LMM.cpp -O2 -Wall -std=c++20
```

This will generate an executable file named `fv_2LMM`,
and then generate the linear feature desctipor files by using command:

```
./fv_2LMM DATASET.sdf OUTPUT
```

Here:
- DATASET.sdf: the input sdf file containing molecular information of the data set, and
- OUTPUT: the prefix for the output files.
  
When the generator finishes normally, it will generate the following files:
- OUTPUT_desc.csv: The csv file of the linear descriptors of the dataset **BEFORE** normalization.
- OUTPUT_desc_norm.csv: The csv file of the linear descriptors of the dataset **AFTER** normalization.
- OUTPUT_fringe.txt: The file containing the _fringe tree_ information of the dataset, which will be used in the stage of solving inverse problem (Module 3).

A sample usage:

```
./fv_2LMM ./sample_instance/sample1.sdf ./sample_instance/output
```

