## Constructing Prediction Functions

This folder contains the codes used to construct prediction functions that is described in Section 5.2, "ML Experiments (Stage 3)".

### Preliminary experiment:

```
python Lasso.py DATASET_desc_norm.csv DATASET_values.txt ALPHA
```

Here:
- DATASET_desc_norm.csv: the 2LCC descriptor files generated in [Module 1](/2LCC/Module_1);
- DATASET_values.txt: the file containing observed value information of the data set, and;
- ALPHA: ???.

(To Song-san: Please explain what ALPHA is and move what you described in the pdf file HERE!!!)
(Also describe what the output will be !!!)


### Evaluation experiment:

~~After obtaining the best parameters, please use "lasso_eval_linreg.py" to generate the prediction function.~~

(To Song-san: Also describe how to use the code "lasso_eval_linreg.py" like the above one !!!)


#### A sample usage:

For preliminary experiment:

```
python LASSO.py ./sample_instance/sample_desc_norm.csv ./sample_instance/sample_norm_values.txt 1
```
