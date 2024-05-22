## Constructing Prediction Functions

This folder contains the codes used to construct the prediction functions between the given descriptor file of the data set (`SINGLE_xxx`).

In both cases, we conduct the experiment in two steps, namely:
- the preliminary experiment (`xxx_pre`) to specify the hyperparameters, and
- the evaluation experiment (`xxx_eval`) to evaluate the learning performance by 10 times 5-cross validation, and construct the prediction function files that will be used later in the stage of solving inverse problem ([Module 3](Polymer/Module_3)).

### Learning with the whole data set

Here we explain the usage of learning for the data set (`SINGLE_xxx`) by Lasso linear regression. 

#### Preliminary experiment:

```
python SINGLE_pre.py DATASET_desc_norm.csv DATASET_values.txt -l
```

Here:
- DATASET_desc_norm.csv (DATASET_hK_desc_norm.csv): the linear(quadratic) descriptor files generated in [Module 1](/HPS/Module_1) of the whole given data set;
- DATASET_values.txt: the file containing observed value information of the data set, and;
- -l: representing the learning method Lasso linear regression.

#### Evaluation experiment:

```
python SINGLE_eval.py DATASET_desc_norm.csv DATASET_values.txt -l (...)
```

Basically the inputs are the same as the ones used for the preliminary experiment, except:
- -l (...): the learning method **AND** the corresponding parameters obtained in the preliminary experiment (see the example below).

#### A sample usage:

For preliminary experiment:

```
python SINGLE_pre.py ./sample_instance/At_large_var0_desc_norm.csv ./sample_instance/At_large_norm_values.txt -l
```

The output will be like:

```
At_large_var0	448	254	254	0	Lasso	0.5797154425805995	0.408987110958712	167.51508474349976		-l 0.00073
```

And for evaluation experiment:

```
python SINGLE_eval.py ./sample_instance/At_large_var0_desc_norm.csv ./sample_instance/At_large_norm_values.txt -l 0.00073
```

Here `-l 0.00073` at the end of the output is the selected hyperparameter obtained from the preliminary experiment.

When the code finishes normally, it will output the information about the learning performance like:

| Data set | \#instance | \#descriptors | \#linear descriptors | \#quadratic descriptors | learning method | median of train R<sup>2</sup> | min of train R<sup>2</sup> | max of train R<sup>2</sup> | median of test R<sup>2</sup> | min of test R<sup>2</sup> | max of test R<sup>2</sup> | running time(sec) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| At_large_var0 | 448 | 254 | 254 | 0 | Lasso | 0.5773279046164392 | 0.5366019574027185 | 0.6211274634517285 | 0.3911713431555095 | 0.038081164737786555 | 0.521242342190777 | 0.16965603828430176 |

Besides above, the code will also generate two summary files containing the information of the R<sup>2</sup> score for each cross-validation under the name `SINGLE_eval_DATASET_log.txt` in `./log_eval/` folder. One line contains the information of one cross-validation, and is like:

| i-th time | j-th cross-validation | train R<sup>2</sup> | test R<sup>2</sup> | running time(sec) | 
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | 0.8420611575174265 | 0.6587311490954544 | 0.08189535140991211 | 

Also the files about the constructed prediction functions in `./pred_func/` folder, which will be used in [Module 3](Polymer/Module_3):
- DATASET_i_j_linreg.txt: the file containing the weights and bias of Lasso linear regression constructed for the data set for the i-th time j-th cross-validation.

