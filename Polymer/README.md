<p align="center">
  <a href="/Polymer/README.md">English</a>
  ·
  <a href="/Polymer/README_jp.md">日本語</a>
</p>

# mol-infer/Polymer

This folder contains the codes and instances used for the paper "A Method for Inferring Polymers Based on Linear Regression and Integer Programming".

- The paper is submitted to IEEE/ACM Transactions on Computational Biology and Bioinformatics.
   - The supplementary material is available as [Supplementary_Materials.pdf](Supplementary_Materials.pdf). 
- The preprint is available at https://arxiv.org/abs/2109.02628.
- All folder/file names in this repository are based on the table/figure numbers in the submitted paper.


The folder structure is organized as follows:
1. Module 1
   - Generate Descriptors (described in Section 3)
   - Preprocessors
1. Module 2
   - Constructing Prediction Functions (Section 4, Stage 3)
1. Module 3
   - Solving an MILP for the Inverse Problem (Section 4, Stage 4, Table 2)
   - Inferring a Polymer with Target Values in Multiple Properties (Section 4, Stage 4, Table 3)
1. Module 4
   - Generating Recombination Solutions (Section 4, Stage 5)
1. Instances for the paper
   - Original data sets used in the paper
   - Instance files used for Tables 2 and 3
  
Please refer each folder for detailed descriptions.

