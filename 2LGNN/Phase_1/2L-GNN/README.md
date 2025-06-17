## 2L-GNN

This folder contains the codes used to construct prediction functions, which is described in Section 3.1.
The code is particularly written for the QM9 dataset.

The parameter of the GNN model can be changed by directly modifying the value of the variable `best_params` in the file `GNN_main_qm9.py`.

### A Sample Usage

Usage:

```
python GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o OUTPUT -val PROP
```

Here:
- ./qm9_all.sdf, ./qm9_complete_values.csv: The files for QM9 dataset, which can be found in [instances_for_paper](/2LGNN/instances_for_paper).
- OUTPUT: The prefix of the output files which contains the information of weights about the prediction function and will be used in [Phase 2](/2LGNN/Phase_2)
- PROP: The name of the property in QM9 dataset, e.g., homo, lumo, mu.

The code will generate the following files when finishing normally:
- result_file.txt: A file containing the summary of the learning process;
- ./model/OUTPUT.pth: A file containing the information of weights about the prediction function;
- ./model/OUTPUT_config.txt: A file containing the information of parameters about the prediction function.