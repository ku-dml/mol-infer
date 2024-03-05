## Solving an MILP for the Inverse Problem

This folder contains the code to generate **ONE** chemical graph by solving an MILP formulation for the inverse problem of the machine learning method, described in Section 6.2, paragraph "Solving an MILP for the Inverse Problem", once given:
- files describing the hyperplane used to split the data set to two subsets;
- files describing the weights/parameters of the prediction functions for the two subsets obtained from [Module 2](HPS/Module_2);
- files describing the topological specification $\sigma$;
- two real numbers $y_l$ and $y_u$, specifying the range of target $\[ y_l, y_u \]$.
  
Also the installation of [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) are necessary.
The CPLEX path should be specified in 'infer_2LMM_SEP.py' before running the code.

We refer [\[22\]](https://arxiv.org/abs/2107.02381) for more details of the MILP formulation and other things like topological specification. 

Usage:

```
python infer_2LMM_SEP.py DATASET y_l y_u SPEC.txt FRINGE.txt OUTPUT -d1 (...) -d2 (...)
```

Here:
- DATASET: the prefix of the necessary files describing the information of the data set used. Detailed information will be expressed below;
- y_l, y_u: the two real numbers $y_l$ and $y_u$;
- SPEC.txt: the file describing the topological specification except fringe trees;
- FRINGE.txt: the file describing the available fringe trees in the topological specification;
- OUTPUT: the prefix of output files;
- -d1 (...): specifying the learning method (`ANN` or `R-MLR`(RLR)) used for the first subset, and;
- -d2 (...): specifying the learning method (`ANN` or `R-MLR`(RLR)) used for the second subset.

To describe the information of the data set, separating hyperplane, and the prediction functions, the following files are necessary:
- DATASET_desc.csv: linear descriptor file of the original data set **BEFORE NORMALIZATION**, obtained in [Module 1](HPS/Module_1);
- DATASET_desc_norm.csv: linear descriptor file of the original data set **AFTER NORMALIZATION**, obtained in [Module 1](HPS/Module_1);
- DATASET_fringe.txt: file containing all the fringe tree information of the original data set, obtained in [Module 1](HPS/Module_1);
- DATASET_values.txt: file containing the observed values of the original data set **BEFORE NORMALIZATION**;
- DATASET_sep.txt: file describing the hyperplane used to split the data set, obtained in [Module 2/Splitting Data Sets via Hyperplanes](HPS/Module_2/Splitting_Data_Sets_via_Hyperplane);
- DATASET_D1_values.txt: file containing the observed values of the first subset, obtained in [Module 2/Splitting Data Sets via Hyperplanes](HPS/Module_2/Splitting_Data_Sets_via_Hyperplane), and
- DATASET_D2_values.txt: file containing the observed values of the second subset, obtained in [Module 2/Splitting Data Sets via Hyperplanes](HPS/Module_2/Splitting_Data_Sets_via_Hyperplane);

Moreover, the following files containing the information of the prediction functions are necessary, which can be obtained in [Module 2/Constructing Prediction Functions](HPS/Module_2/Constructing_Prediction_Functions):
- DATASET_D1_desc_norm_selected.csv: file containing the selected/reduced features for the first subset during the computation;
- DATASET_D2_desc_norm_selected.csv: file containing the selected/reduced features for the second subset during the computation;
- DATASET_Di_linreg.txt: (only when using RLR) file containing the weights and bias of linear regression of RLR, i=1,2, and;
- DATASET_Di_weights.txt, DATASET_Di_biases.txt: (only when using ANN) files containing the weights and biases of the neural networks, i=1,2.

When the code finishes normally, it will generate the following files:
- OUTPUT.lp: the lp file describing the MILP formulation;
- OUTPUT.sdf: (when MILP is feasible) the sdf file of the inferred chemical graph, and;
- OUTPUT_partition.txt: (when MILP is feasible) the partition file of the inferred chemical graph, which will be used in [Module 4](HPS/Module_4).
  
If MILP is considered as "Infeasible" by the solver, it means such a chemical graph does not exist under the given constraints/specifications.


A sample usage:

```
python infer_2LMM_SEP.py ./sample_instance/At 240 250 ./sample_instance/instance_a_2LMM.txt ./sample_instance/ins_a_fringe_2LMM.txt ./sample_instance/At_240_250 -d1 -ANN -d2 R-MLR
```
