## Inferring a Polymer with Target Values in Multiple Properties

This folder contains the code to generate **ONE** chemical graph for **THREE** properties at the same time by solving an MILP formulation for the inverse problem of Lasso linear regression, described in Section 4, paragraph "Inferring a Polymer with Target Values in Multiple Properties", once given:
- fils describing the weights and biases of the prediction functions obtained from [Module 2](/Polymer/Module_2);
- files describing the topological specification $\sigma$;
- two real numbers $y_l^i$ and $y_u^i$, specifying the range of target $\[ y_l^i, y_u^i \]$ for each property, i=1,2,3.
  
Also the installation of [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) are necessary.
The CPLEX path should be specified in 'infer_2LMM_polymer.py' before running the code.

We refer [\[23\]](https://arxiv.org/abs/2107.02381) for more details of the MILP formulation and other things like topological specification. 

Usage:

```
python infer_2LMM_polymer.py DATASET_1 DATASET_2 DATASET_3 y_l^1 y_u^1 y_l^2 y_u^2 y_l^3 y_u^3 SPEC.txt FRINGE.txt OUTPUT
```

Here:
- DATASET_i: the prefixes of the necessary files describing the information of the data set used for the three properties, i=1,2,3. Detailed information will be expressed below;
- y_l^i, y_u^i: the two real numbers $y_l^i$ and $y_u^i$ for each property, i=1,2,3;
- SPEC.txt: the file describing the topological specification except fringe trees;
- FRINGE.txt: the file describing the available fringe trees in the topological specification;
- OUTPUT: the prefix of output files.

To describe the information of the data set and the prediction functions, the following files are necessary:
- DATASET_i_desc.csv: descriptor files of the original data set **BEFORE NORMALIZATION**, obtained in [Module 1](/Polymer/Module_1) for each property, i=1,2,3;
- DATASET_i_desc_norm.csv: descriptor files of the original data set **AFTER NORMALIZATION**, obtained in [Module 1](/Polymer/Module_1) for each property, i=1,2,3;
- DATASET_i_fringe.txt: files containing all the fringe tree information of the original data set, obtained in [Module 1](/Polymer/Module_1) for each property, i=1,2,3;
- DATASET_i_values.txt: files containing the observed values of the original data set **BEFORE NORMALIZATION** for each property, i=1,2,3;

Moreover, the following files containing the information of the prediction functions are necessary, which can be obtained in [Module 2/Constructing Prediction Functions](/Polymer/Module_2/Constructing_Prediction_Functions):
- DATASET_i_linreg.txt: files containing the weights and bias of linear regression for each property, i=1,2,3.

When the code finishes normally, it will generate the following files:
- OUTPUT.lp: the lp file describing the MILP formulation;
- OUTPUT.sdf: (when MILP is feasible) the sdf file of the inferred chemical graph, and;
- OUTPUT_partition.txt: (when MILP is feasible) the partition file of the inferred chemical graph, which will be used in [Module 4](/Polymer/Module_4).
  
If MILP is considered as "Infeasible" by the solver, it means such a chemical graph does not exist under the given constraints/specifications.


A sample usage:

```
python infer_2LMM_polymer.py ./sample_instance/AmpD/AmpD ./sample_instance/HcLiq/HcLiq ./sample_instance/Tg/Tg 1.200 1.224 624.0 628.0 171.0 174.0 ./sample_instance/instance_P1_polymer.txt ./sample_instance/ins_P1_fringe_polymer.txt ./sample_instance/AmpD_HcLiq_Tg_1.200_1.224_624.0_628.0_171.0_174.0
```
