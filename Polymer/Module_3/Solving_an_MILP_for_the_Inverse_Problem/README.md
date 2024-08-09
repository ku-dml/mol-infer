## Solving an MILP for the Inverse Problem

This folder contains the code to generate **ONE** chemical graph by solving an MILP formulation for the inverse problem of Lasso linear regression, described in Section 4, paragraph "Stage 4", once given:
- a file describing the weights and bias of the prediction functions obtained from [Module 2](/Polymer/Module_2);
- files describing the topological specification $\sigma$;
- two real numbers $y_l$ and $y_u$, specifying the range of target $\[ y_l, y_u \]$.
  
Also the installation of [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) are necessary.
The CPLEX path should be specified in 'infer_2LMM_polymer.py' before running the code.

We refer [\[23\]](https://arxiv.org/abs/2107.02381) for more details of the MILP formulation and other things like topological specification. 

Usage:

```
python infer_2LMM_polymer.py DATASET y_l y_u SPEC.txt FRINGE.txt OUTPUT (fq_l fq_u)
```

Here:
- DATASET: the prefix of the necessary files describing the information of the data set used. Detailed information will be expressed below;
- y_l, y_u: the two real numbers $y_l$ and $y_u$;
- SPEC.txt: the file describing the topological specification except fringe trees;
- FRINGE.txt: the file describing the available fringe trees in the topological specification;
- OUTPUT: the prefix of output files;
- (fq_l, fq_u): lower and upper bounds of the descriptor 'fq' only used for the property Prm.

To describe the information of the data set and the prediction functions, the following files are necessary:
- DATASET_desc.csv: descriptor file of the original data set **BEFORE NORMALIZATION**, obtained in [Module 1](/Polymer/Module_1);
- DATASET_desc_norm.csv: descriptor file of the original data set **AFTER NORMALIZATION**, obtained in [Module 1](/Polymer/Module_1);
- DATASET_fringe.txt: file containing all the fringe tree information of the original data set, obtained in [Module 1](/Polymer/Module_1);
- DATASET_values.txt: file containing the observed values of the original data set **BEFORE NORMALIZATION**;

Moreover, the following files containing the information of the prediction functions are necessary, which can be obtained in [Module 2/Constructing Prediction Functions](/Polymer/Module_2/Constructing_Prediction_Functions):
- DATASET_linreg.txt: file containing the weights and bias of linear regression.

When the code finishes normally, it will generate the following files:
- OUTPUT.lp: the lp file describing the MILP formulation;
- OUTPUT.sdf: (when MILP is feasible) the sdf file of the inferred chemical graph, and;
- OUTPUT_partition.txt: (when MILP is feasible) the partition file of the inferred chemical graph, which will be used in [Module 4](/Polymer/Module_4).
  
If MILP is considered as "Infeasible" by the solver, it means such a chemical graph does not exist under the given constraints/specifications.


A sample usage:

```
python infer_2LMM_polymer.py ./sample_instance/AmpD 0.885 0.890 ./sample_instance/instance_a_polymer.txt ./sample_instance/ins_a_fringe_2LMM.txt ./sample_instance/At_0.885_0.890 
```
