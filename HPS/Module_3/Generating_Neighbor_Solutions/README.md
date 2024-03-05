## Generating Neighbor Solutions

This folder contains the code to generate the *neighbor solutions*, or *Grid Neighbor Search* (GNS), described in Section 6.2, paragraph "Generating Neighbor Solutions" by solving an MILP formulation for the inverse problem of the machine learning method systematically, once given:
- files describing the hyperplane used to split the data set to two subsets;
- files describing the weights/parameters of the prediction functions for the two subsets obtained from [Module 2](HPS/Module_2);
- files describing the topological specification $\sigma$;
- two real numbers $y_l$ and $y_u$, specifying the range of target $\[ y_l, y_u \]$;
- a set of additional linear constraints.
  
Also the installation [PuLP](https://coin-or.github.io/pulp/index.html) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) are necessary.
The CPLEX path should be specified in 'infer_2LMM_SEP.py' before running the code.
The default time limit for CPLEX to solve the augmented MILP is set to 300s.

We refer [\[24\]](https://www.computer.org/csdl/proceedings-article/bibm/2021/09669710/1A9VAbXVZJu) for the idea and [GNS package](Grid-neighbor-search) for more details.

Usage:

```
python infer_2LMM_SEP.py DATASET y_l y_u SPEC.txt FRINGE.txt OUTPUT GNS_PARAMETER.txt -d1 (...) -d2 (...) -t (...)
```

Here:
- DATASET: the prefix of the necessary files describing the information of the data set used. Detailed information will be expressed below;
- y_l, y_u: the two real numbers $y_l$ and $y_u$;
- SPEC.txt: the file describing the topological specification except fringe trees;
- FRINGE.txt: the file describing the available fringe trees in the topological specification;
- OUTPUT: the prefix of output files;
- GNS_PARAMETER.txt: the file describing the parameters related to GNS;
- -d1 (...): specifying the learning method (`ANN` or `R-MLR`(RLR)) used for the first subset;
- -d2 (...): specifying the learning method (`ANN` or `R-MLR`(RLR)) used for the second subset, and
- -t (...): the set of additional linear constraints (see [GNS](Grid-neighbor-search) for more details of this).

To describe the information of the data set, separating hyperplane, and the prediction functions, the following files are necessary:
- DATASET_desc.csv: linear descriptor file of the original data set **BEFORE NORMALIZATION**, obtained in [Module 1](HPS/Module_1);
- DATASET_desc_norm.csv: linear descriptor file of the original data set **AFTER NORMALIZATION**, obtained in [Module 1](HPS/Module_1);
- DATASET_fringe.txt: file containing all the fringe tree information of the original data set, obtained in [Module 1](HPS/Module_1);
- DATASET_values.txt: file containing the observed values of the original data set **BEFORE NORMALIZATION**;
- DATASET_sep.txt: file describing the hyperplane used to split the data set, obtained in [Module 2/Splitting Data Sets via Hyperplanes](HPS/Module_2/Splitting_Data_Sets_via_Hyperplane);
- DATASET_D1_values.txt: file containing the observed values of the first subset, obtained in [Module 2/Splitting Data Sets via Hyperplanes](HPS/Module_2/Splitting_Data_Sets_via_Hyperplane), and;
- DATASET_D2_values.txt: file containing the observed values of the second subset, obtained in [Module 2/Splitting Data Sets via Hyperplanes](HPS/Module_2/Splitting_Data_Sets_via_Hyperplane);

Moreover, the following files containing the information of the prediction functions are necessary, which can be obtained in [Module 2/Constructing Prediction Functions](HPS/Module_2/Constructing_Prediction_Functions):
- DATASET_D1_desc_norm_selected.csv: file containing the selected/reduced features for the first subset during the computation;
- DATASET_D2_desc_norm_selected.csv: file containing the selected/reduced features for the second subset during the computation;
- DATASET_Di_linreg.txt: (only when using RLR) file containing the weights and bias of linear regression of RLR, i=1,2, and;
- DATASET_Di_weights.txt, DATASET_Di_biases.txt: (only when using ANN) files containing the weights and biases of the neural networks, i=1,2.

When the code finishes normally, it will generate the following files:
- OUTPUT.lp: the lp file describing the MILP formulation;
- OUTPUT.sdf: (when MILP is feasible) the sdf file of the inferred chemical graph, and;
- OUTPUT_partition.txt: (when MILP is feasible) the partition file of the inferred chemical graph, which will be used in [Module 4](HPS/Module_4);

and for the augmented MILP corresponding to the i-th subspace, it will generate the following files:
- OUTPUT_i.lp: the lp file describing the MILP formulation corresponding to the i-th subspace;
- OUTPUT_i.sdf: (when MILP is feasible) the sdf file of the inferred chemical graph corresponding to the i-th subspace, and;
- OUTPUT_i_partition.txt: (when MILP is feasible) the partition file of the inferred chemical graph corresponding to the i-th subspace, which will be used in [Module 4](HPS/Module_4).
  
If MILP is considered as "Infeasible" by the solver, it means such a chemical graph does not exist under the given constraints/specifications. 

And finally:
- OUTPUT_GS_log.csv: a summary file containing all the results obtained for the given inputs.


A sample usage:

```
python infer_2LMM_SEP_GNS.py ./sample_instance/At 370 380 ./sample_instance/instance_c_2LMM.txt ./sample_instance/ins_c_fringe_2LMM.txt ./sample_instance/At_370_380 ./sample_instance/p_max_delta_r_1e-1.txt -d1 -ANN -d2 R-MLR -t ./sample_instance/constraints/Mp/Mp ./sample_instance/constraints/Sl/Sl
```
