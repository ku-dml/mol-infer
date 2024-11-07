<p align="center">
  <a href="/multimodel/README.md">English</a>
  ·
  <a href="/multimodel/README_jp.md">日本語</a>
</p>

# mol-infer/multimodel

This folder contains the codes and instances used to infer a chemical graph with multiple target values.

The folder structure is organized as follows:
1. Module 3
    - Solving an MILP for the inverse problem
    - Inferring a chemical graph with target values for multiple properties

# Module 3

Module 3 consists of the stage of solving an MILP formulation for the inverse problem to infer a chemical graph.

input:
- Instance_file
- Fringe_tree_file
- Prediction functions obtained from Module 2
- Target value ranges for multiple properties

output:
- Chemical graph with the predicted property values inside the target value ranges


to check the calculated descriptors
```bash
cd libs/2LMM_v019/
make FV_2LMM_V019
```

sample command to run the code
```bash
python3 -m Module_3 config/config.yaml
```


# config file
The config file is a yaml file containing the following parameters:
- instance_file: path to the instance file
- fringe_tree_file: path to the fringe tree file
- output_prefix: prefix for the output files
- input_data: a list of dictionaries containing the following parameters:
  - model: the model used for prediction (LR, RF, ANN)
  - prefix: the prefix for the input data files
  - target_value_lower_bound: the lower bound of the target value range
  - target_value_upper_bound: the upper bound of the target value range
