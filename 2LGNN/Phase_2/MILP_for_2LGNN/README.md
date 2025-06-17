## MILP_for_2L-GNN

This folder contains the code to generate **ONE** chemical graph by solving an MILP that is described in Section 3.2, "Phase 2: Inverse QSAR/QSPR Phase", once given:
- files describing the prediction functions constructed by 2L-GNN that is obtained from [Phase 1](/2LGNN/Phase_1);
- files describing the specification $\sigma$;
- two real numbers $y_l$ and $y_u$, specifying the range of target $\[ y_l, y_u \]$.

Also the installation of [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) are necessary.

> [!NOTE]
> The CPLEX path should be specified in 'infer_2LMM_GNN.py' before running the code.

Usage:

```
python infer_2LMM_GNN.py PRED_FUNC y_l y_u SPEC.txt FRINGE.txt OUTPUT
```

Here:
- PRED_FUNC: the prefix of the necessary files describing the information of the constructed prediction function. Detailed information will be expressed below;
- y_l, y_u: the two real numbers $y_l$ and $y_u$;
- SPEC.txt: the file describing the topological specification except fringe trees;
- FRINGE.txt: the file describing the available fringe trees in the topological specification, and;
- OUTPUT: the prefix of output files.

To describe the information of the constructed prediction function, the following files are necessary, which can be obtained in [Phase 1](/2LGNN/Phase_1/2L-GNN):
- PRED_FUNC.pth: file containing the weights and bias of the constructed 2L-GNN model, and;
- PRED_FUNC_config.txt: file containing the necessary information about the parameter of the 2L-GNN model.

When the code finishes normally, it will generate the following files:
- OUTPUT.sdf: (when MILP is feasible) the sdf file of the inferred chemical graph.
  
If MILP is considered as "Infeasible" by the solver, it means such a chemical graph does not exist under the given constraints/specifications.

A sample usage:

```
python infer_2LMM_GNN.py ./sample_instance/lumo -9.50 -9.00 ./sample_instance/instance_e_2LMM.txt /sample_instance/ins_e_fringe.txt lumo_-9.50_-9.00
```

