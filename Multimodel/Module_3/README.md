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