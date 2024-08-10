## Constructing Prediction Functions

This folder contains the codes used to construct prediction functions that is described in Section 5.2, "ML Experiments (Stage 3)".

### Preliminary experiment:

The `Lasso.py` script will automatically search for the best parameter $\lambda$ in Lasso linear regression and output an `.xlsx` file containing detailed ML results. 

Usage:

```
python Lasso.py DATASET_fv_and_value.txt OUTPUT.xlsx LAMBDA
```

Here:
- DATASET_fv_and_value.txt: a file containing the path of the 2LCC descriptor file obtained in [Module 1](/2LCC/Module_1) and the path of the corresponding observed value file for one data set in each line, one line should be organized like `DATASET_desc_norm.csv, DATASET_norm_values.txt`;
- OUTPUT.xlsx: the filename of the output `.xlsx` file, and;
- LAMBDA: the initial parameter $\lambda$ in Lasso linear regression.

For the detailed progress of selecting the parameter $\lambda$, please check the file `how_to_tunning_lambda.pdf`.

<!--
About the strategy of tunning Lambda is as following:
1. Given initial $\lambda_0 = \hat{\lambda}_0 \times 10^{k_0}$, $\hat{\lambda}_0$ must be 1, I'll start from $\lambda_0 = 1$.

2. Search in range $[\lambda_0 - 10^{k_0}, \lambda_0 + 10^{k_0}]$, step $= 10^{k_0-1}$, find the $\lambda$ in above 20 values that makes $\text{TestR}^2$ are biggest as $\lambda_1$, represent $\lambda_1$ as $\lambda_1 = \hat{\lambda}_1 \times 10^{k_1}$.

3. Repeat process 2, until $k_n = k_{n-1}$, then $\lambda_n$ was selected as the best $\lambda$.

For example:

$\lambda_0 = 1 = 1 \times 10^0$, then searching in $\{0.1, 0.2, 0.3, \dots, 2\}$ suppose that $\Rightarrow 
\lambda_1 = 0.1 = 1 \times 10^{-1}$, then search in $\{0.01, 0.02, \dots, 0.2\}$, suppose that $\Rightarrow 
\lambda_2 = 0.03 = 3 \times 10^{-2}$, then search in $\{0.021, 0.022, \dots, 0.4\}$, suppose that $\Rightarrow 
\lambda_3 = 0.022 = 2.2 \times 10^{-2}$, $k_2 = k_3$, so take the $\lambda_3$ as best $\lambda$.
-->

#### A sample usage:

For preliminary experiment:

```
python LASSO.py ./sample_instance/sample_fv_value_path.txt ./sample_instance/sample_result.xlsx 1
```

### Generate prediction function:

After obtaining the best parameter $\lambda$ in the preliminary experiment, then use `lasso_eval_linreg.py` to generate the prediction function.

Usage:

```
python lasso_eval_linreg.py DATASET_desc_norm.csv DATASET_norm_values.txt OUTPUT.txt best_lambda
```

Here:
- DATASET_desc_norm.csv: the 2LCC feature vector file of the data set obtained in [Module 1](/2LCC/Module_1);
- DATASET_norm_values.txt: the corresponding observed value file of the data set;
- OUTPUT: the name of the output file which contains the information about the prediction function and will be used in [Module 3](/2LCC/Module_3), and;
- best_lambda: the value of the parameter $\lambda$ found in the preliminary experiment.

#### A sample usage:

For generating prediction function:

```
python lasso_eval_linreg.py ./sample_instance/sample_desc_norm.csv ./sample_instance/sample_norm_values.txt ./sample_instance/sample_LR.txt 0.001
```

