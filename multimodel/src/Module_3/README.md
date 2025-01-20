<p align="center">
  <a href="/multimodel/src/Module_3/README.md">English</a>
  ·
  <a href="/multimodel/src/Module_3/README_jp.md">日本語</a>
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
python3 -m Module_3 --config-name=config_sample
```


# config file
For the settings of the config file, refer to config/config_sample.yaml.
