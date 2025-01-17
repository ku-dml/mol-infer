## Solving an MILP for the Inverse Problem

This folder contains the code to generate **ONE** chemical graph by solving an MILP formulation for the inverse problem of the machine learning method, described in Sections 4.2 and 4.3, particularly for the case of $\chi$-parameter, once given:
- files describing the weights/parameters of the prediction functions for the two subsets obtained from [Module 2](/chi-parameter/Module_2);
- files describing the topological specification $\sigma$;
- two real numbers $y_l$ and $y_u$, specifying the range of target $\[ y_l, y_u \]$.
  
Also the installation of [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) are necessary.
The CPLEX path should be specified in 'infer_2LMM_polymer_for_chi.py' before running the code.

We refer [\[22\]](https://arxiv.org/abs/2107.02381) for more details of the MILP formulation and other things like topological specification. 

Usage:

```
python infer_2LMM_polymer_for_chi.py DATASET y_l y_u SPEC.txt FRINGE.txt SOLVENT_fv.csv OUTPUT -ml (...) -T (...)
```

Here:
- DATASET: the prefix of the necessary files describing the information of the data set used. Detailed information will be expressed below;
- y_l, y_u: the two real numbers $y_l$ and $y_u$;
- SPEC.txt: the file describing the topological specification except fringe trees;
- FRINGE.txt: the file describing the available fringe trees in the topological specification;
- SOLVENT_fv.csv: the file containing the feature vector for the fixed solvent;
- OUTPUT: the prefix of output files;
- -ml: specifying the learning method (`ANN` or `LR`(for R-MLR) or `RF`), and;
- -T: the fixed temperature.

To describe the information of the data set, and the prediction functions, the following files are necessary:
- DATASET_desc.csv: linear descriptor file of the original data set **BEFORE NORMALIZATION**, obtained in [Module 1](/chi-paramter/Module_1);
- DATASET_desc_norm.csv: linear descriptor file of the original data set **AFTER NORMALIZATION**, obtained in [Module 1](/chi-paramter/Module_1);
- DATASET_fringe.txt: file containing all the fringe tree information of the original data set of the solute, obtained in [Module 1](/chi-paramter/Module_1), and;
- DATASET_values.txt: file containing the observed values of the original data set **BEFORE NORMALIZATION**.

Moreover, the following files containing the information of the prediction functions are necessary, which can be obtained in [Module 2/Constructing Prediction Functions](/chi-paramter/Module_2/Constructing_Prediction_Functions):
- DATASET_desc_norm_selected.csv: the file containing selected feature vectors during the feature selection stage of learning;
- DATASET_linreg.txt: (only when using R-MLR) file containing the weights and bias of linear regression of R-MLR;
- DATASET_weights.txt, DATASET_biases.txt: (only when using ANN) files containing the weights and biases of the neural networks, and;
- DATASET_rf.txt: (only when using RF) file containing the weights of random forest.


When the code finishes normally, it will generate the following files:
- OUTPUT.lp: the lp file describing the MILP formulation;
- OUTPUT.sdf: (when MILP is feasible) the sdf file of the inferred chemical graph, and;
- OUTPUT_partition.txt: (when MILP is feasible) the partition file of the inferred chemical graph, which will be used in [Module 4](/chi-parameter/Module_4).
  
If MILP is considered as "Infeasible" by the solver, it means such a chemical graph does not exist under the given constraints/specifications.


A sample usage:

```
python infer_2LMM_polymer_for_chi.py ./sample_instance/Chi_JOCTA 0.2 0.7 ./sample_instance/instance_c1_polymer.txt ./sample_instance/ins_c1_fringe.txt ./sample_instance/solvent_fv_Chi_JOCTA_PE.csv ./sample_instance/Chi_JOCTA_c1_0.2_0.7 -ml ANN -T 298.0
```
