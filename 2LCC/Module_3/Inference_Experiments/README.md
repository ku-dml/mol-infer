## Inference Experiments

This folder contains the code to generate **ONE** chemical graph by solving an MILP that is described in Section 5.3, "Inference Experiments (Stage 4)", once given:
- a file describing the weights and bias of the prediction functions obtained from [Module 2](/2LCC/Module_2);
- files describing the specification $\sigma$;
- two real numbers $y_l$ and $y_u$, specifying the range of target $\[ y_l, y_u \]$.

Also the installation of [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) are necessary.

> [!NOTE]
> The CPLEX path should be specified in 'infer_graph_LR.py' before running the code.

Usage:

```
python infer_graph_LR.py DATASET SEED_TREE.txt FRINGE.txt OUTPUT OUTPUT_INSTANCE_NAME y_l y_u
```

Here:
- DATASET: the prefix of the necessary files describing the information of the data set used. Detailed information will be expressed below;
- SEED.txt: the file describing the specification (including the seed tree) except fringe trees;
- FRINGE.txt: the file describing the available fringe trees in the specification;
- OUTPUT: the name of the output SDF file;
- OUTPUT_INSTANCE_NAME: the name of the output chemical graph, and;
- y_l, y_u: the two real numbers $y_l$ and $y_u$.

To describe the information of the data set and the prediction functions, the following files are necessary:
- DATASET_desc.csv: descriptor file of the original data set **BEFORE NORMALIZATION**, obtained in [Module 1](/2LCC/Module_1);
- DATASET_fringe.txt: file containing all the fringe tree information of the original data set, obtained in [Module 1](/2LCC/Module_1);
- DATASET_values.txt: file containing the observed values of the original data set **BEFORE NORMALIZATION**;

Moreover, the following files containing the information of the prediction functions are necessary, which can be obtained in [Module 2/Constructing Prediction Functions](/2LCC/Module_2/Constructing_Prediction_Functions):
- DATASET_LR.txt: file containing the weights and bias of Lasso linear regression.

When the code finishes normally, it will generate the following files:
- OUTPUT: (when MILP is feasible) the sdf file of the inferred chemical graph.
  
If MILP is considered as "Infeasible" by the solver, it means such a chemical graph does not exist under the given constraints/specifications.

A sample usage:

```
python infer_graph_LR.py ./sample_instance/IhcLiq ./sample_instance/IhcLiq_3_a0.txt ./sample_instance/IhcLiq_fringe.txt ./IhcLiq_3_a0.sdf IhcLiq_3_a0 75.32 1956.1
```

