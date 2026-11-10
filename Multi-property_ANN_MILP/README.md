This folder contains the codes and datasets used in the following paper:

Raveena Rai, Naveed Ahmed Azam, Jianshen Zhu, Kazuya Haraguchi, Liang Zhao, and Tatsuya Akutsu, "Mathematically Guaranteed Multi-Property Inverse Molecular Design via Shared Representation and Mixed-Integer Linear Programming."

The folder structure is organized as follows:

## Module 1

- **a.** Pre-processing: Extracting SDF files from the QM9 dataset using RDKit and eliminating compounds, as described in Section 2.1. The dataset used in the experiments is described in Section 3.
- **b.** Generating molecular descriptors, as described in Section 2. The experimental results based on these descriptors are given in Section 3.1.

## Module 2

Constructing the prediction function using LLR–FSP–ANN-Multi, as described in Section 2.1. The experimental results are given in Section 3.1, namely shared descriptor selection (Section 3.1.1), selected ANN architectures (Section 3.1.2), and performance of the selected ANN (Section 3.1.3).

## Module 3

Solving the MILP formulation for the inverse design problem with multiple desired properties, as described in Section 2.2. The experimental results are given in Section 3.2.

Please refer to each folder for detailed descriptions.
