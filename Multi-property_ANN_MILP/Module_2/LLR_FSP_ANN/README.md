**LLR\_Multi experiment:**



Usage:

&#x09;python  LLR\_Multi.py (.csv)  (\_values.txt)  (number\_of\_properties)  (property\_1 ... property\_n)  (union/intersection)



Sample:

&#x09;python  LLR\_Multi.py  ./Datasets/qm9\_desc\_norm.csv  ./Datasets/qm9\_norm\_values.txt  2  alpha  r2  union





The Lasso procedure is applied separately to each specified property. The selected descriptor sets are then combined using either union or intersection.



The output will be something like:

&#x09;./log/LLR\_ANN\_alpha\_lasso\_log.txt

&#x09;./log/LLR\_ANN\_r2\_lasso\_log.txt

&#x09;./log/LLR\_ANN\_alpha\_r2\_desc\_norm.csv



The individual log files contain the Lasso results for each property. The combined descriptor file contains the shared descriptors selected for the specified property combination and is used as the input feature file for the FSP\_Multi experiment.



==============================================================================



**FSP\_Multi experiment:**



Usage:

&#x09;python  FSP\_Multi.py  (.csv)  (\_values.txt)  (number\_of\_descriptors)  (property\_1 \[property\_2 ...])  (multi/lasso)



Sample:

&#x09;python  FSPMulti.py  ./log/LLR\_ANN\_alpha\_r2\_desc\_norm.csv  ./Datasets/qm9\_norm\_values.txt  20  alpha  r2  multi





The number of target properties is determined automatically from the property names supplied in the command. The final argument specifies the regression model used during forward stepwise selection:

&#x09;multi = multiple linear regression

&#x09;lasso = Lasso regression



The output will be something like:

&#x09;	LLR\_FSP\_alpha\_r2\_20\_desc\_norm.csv



The file contains the CID column and the selected descriptor values. Please use the selected feature-vector file as the descriptor input for the ANN\_Multi\_Pre\_Eval experiment.



==============================================================================



**ANN\_Multi\_Pre\_Eval experiment:**



Usage:

&#x09;python  ANN\_Multi\_Pre\_Eval.py  (.csv)  (\_norm\_values.txt)  (\_original\_values.txt)  (number\_of\_descriptors\_list)  (property\_1 \[property\_2 ...])



Sample:

&#x09;python  ANN\_Multi\_Pre\_Eval.py  ./selected\_feature\_vectors\_alpha\_r2\_20\_multi.csv  ./Datasets/qm9\_norm\_values.txt  ./Datasets/qm9\_values.txt  20 15  alpha  r2





This code performs both the preliminary ANN architecture search and the final evaluation experiment in a single run. It supports one or more target properties.



The ANN\_Result folder and the following four subfolders must already exist:

&#x09;ANN\_Result/ANN\_pre\_result

&#x09;ANN\_Result/ANN\_Eval

&#x09;ANN\_Result/Result

&#x09;ANN\_Result/weight\_biases



The output will be saved as follows:

&#x09;ANN\_Result/ANN\_pre\_result

&#x09;	Optuna architecture-search results



&#x09;ANN\_Result/ANN\_Eval

&#x09;	Iteration-wise and fold-wise training, testing, and average R2 scores



&#x09;ANN\_Result/Result

&#x09;	Final R2, MAE, and NMAE summary



&#x09;ANN\_Result/weight\_biases

&#x09;	ANN weight and bias files for all 10 x 5-fold evaluation models



The code does not create any folder automatically.

