<p align="center">
  <a href="/multimodel/Module_3/README.md">English</a>
  ·
  <a href="/multimodel/Module_3/README_jp.md">日本語</a>
</p>

# Module 3

Module 3 consists of the stage of solving an MILP formulation for the inverse problem to infer a chemical graph.

input:
- Instance_file
- Fringe_tree_file
- Prediction functions obtained from Module 2
- Target value ranges for the properties

output:
- Chemical graph with the predicted property values inside the target value ranges


to check the calculated descriptors
```bash
make -C Module_3/libs/2LMM_v019/ FV_2LMM_V019
```

sample command
```bash
python3 -m Module_3 Module_3/config/config.yaml
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
