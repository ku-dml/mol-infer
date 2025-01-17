## Constructing Prediction Functions

This folder contains the codes used to construct the prediction functions between the given descriptor file of the data set (`SINGLE_xxx`).

In both cases, we conduct the experiment in two steps, namely:
- the preliminary experiment (`xxx_pre`) to specify the hyperparameters, and
- the evaluation experiment (`xxx_eval`) to evaluate the learning performance by 10 times 5-cross validation, and construct the prediction function files that will be used later in the stage of solving inverse problem ([Module 3](/chi-parameter/Module_3)).

### Learning with the whole data set

Here we explain the usage of learning for the data set (`SINGLE_xxx`) by Lasso linear regression. 

#### Preliminary experiment:

```
python SINGLE_pre.py DATASET_desc_norm.csv DATASET_values.txt (-ann)
```

Here:
- DATASET_desc_norm.csv (DATASET_hK_desc_norm.csv): the linear(quadratic) descriptor files generated in [Module 1](/HPS/Module_1) of the whole given data set;
- DATASET_values.txt: the file containing observed value information of the data set, and;
- (-ann): representing the learning method artificial neural network. (`-rbsp` is used for R-MLR, and `-rf` is used for random forest)

#### Evaluation experiment:

```
python SINGLE_eval.py DATASET_desc_norm.csv DATASET_values.txt (-ann) (...)
```

Basically the inputs are the same as the ones used for the preliminary experiment, except:
- (-ann) (...): the learning method **AND** the corresponding parameters obtained in the preliminary experiment (see the example below).

#### A sample usage:

For preliminary experiment:

```
python SINGLE_pre.py ./sample_instance/Chi_Aoki_T1-1_desc_norm.csv ./sample_instance/Chi_Aoki_values.txt -ann
```

The output will be like:

```
Chi_Aoki_T-1	1190	218	218	0	ANN	0.9902814631598476	0.782947304131302	3361.1716463565826		-ann ./log/LLR_ANN_Chi_Aoki_T-1_var0_desc_norm.csv 0.99 9911 68 68 68 	
```

And for evaluation experiment:

```
python SINGLE_eval.py ./sample_instance/Chi_Aoki_T1-1_desc_norm.csv ./sample_instance/Chi_Aoki_values.txt -ann ./log/LLR_ANN_Chi_Aoki_T-1_var0_desc_norm.csv 0.99 9911 68 68 68 
```

Here `-ann ./log/LLR_ANN_Chi_Aoki_T-1_var0_desc_norm.csv 0.99 9911 68 68 68 ` at the end of the output is the selected hyperparameter obtained from the preliminary experiment. (The file `./log/LLR_ANN_Chi_Aoki_T-1_var0_desc_norm.csv` will be used in [Module 3](/chi-parameter/Module_3).)

When the code finishes normally, it will output the information about the learning performance like:

| Data set | \#instance | \#descriptors | \#linear descriptors | \#quadratic descriptors | learning method | metric | median of train R<sup>2</sup> | min of train R<sup>2</sup> | max of train R<sup>2</sup> | median of test R<sup>2</sup> | min of test R<sup>2</sup> | max of test R<sup>2</sup> | running time(sec) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Chi_Aoki_T-1 | 1190 | 91 | 91 | 0 | ANN | R2 | 0.9902780508167244 | 0.9900185173321279 | 0.9919624783939708 | 0.8009358575320049 | 0.6454732122572822 | 0.8856786111564742 | 229.5918996334076 |

Besides above, the code will also generate two summary files containing the information of the R<sup>2</sup> score for each cross-validation under the name `SINGLE_eval_DATASET_log.txt` in `./log_eval/` folder. One line contains the information of one cross-validation, and is like:

| i-th time | j-th cross-validation | train R<sup>2</sup> | test R<sup>2</sup> | running time(sec) | 
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 0.9907923878733901 | 0.7899040462933161 | 5.587889671325684 | 

Also the files about the constructed prediction functions in `./pred_func/` folder, which will be used in [Module 3](/chi-parameter/Module_3):
- DATASET_i_j_weight.txt: the file containing the weights of neural network constructed for the data set for the i-th time j-th cross-validation,
- DATASET_i_j_biases.txt: the file containing the biases of neural network constructed for the data set for the i-th time j-th cross-validation.

