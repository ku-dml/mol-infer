# Preprocessor
The directory contains two scripts that are used to preprocess your SDF file.

Note that [SDF files that are used in the paper](../../instances_for_paper) are already processed by these two scripts. You don't need to run the scripts for the SDF files in the above directory. 

The two scripts are: 

- **eliminate.py** (strongly recommended to run on your SDF before starting Module 1)
- **limit_atoms.py** (rather optional)

## eliminate.py
Our framework does not deal with molecules that are represented by special chemical graphs; e.g. a chemical graph with less than four carbons, a disconnected chemical graph, etc. **eliminate.py** builds an SDF that is the subset of feasible molecules in a given SDF. 

Let _A.sdf_ be your SDF file. 
To run **eliminate.py** to generate _B.sdf_, 
```
$ python eliminate.py A.sdf B.sdf
```
This _B.sdf_ can be used as an input of the [feature vector generator](../Generate_Linear_Descriptors). 

You can use **Example.sdf** in this directory for example. 
```
$ python eliminate.py Example.sdf Example_eli.sdf
```
Executing this, you will have **Example_eli.sdf**. You will also see that **Example_eli.sdf** contains 1000 molecules, where **Example.sdf** contains 1294 molecules. 

## limit_atoms.py
For a given SDF, you can extract molecules that do not contain any other atoms than specified ones, as we do in the experiments in Section 6.1 of the paper. If you want to extract molecules in B.sdf that do not contain any other atoms than C, N, O and H, run
```
$ python limit_atoms.py B.sdf C N O H
```
Then you will have an SDF file named _B_C_N_O_H.sdf_. Each molecule in this SDF consists of some of C, N, O and H and does not contain any other atoms. Of course this _B_C_N_O_H.sdf_ can be also used as an input of the [feature vector generator](../Generate_Linear_Descriptors). 

For example, 
```
$ python limit_atoms.py Example_eli.sdf C O N S Cl H
```
Then you will have **Example_eli_C_O_N_S_Cl_H.sdf** that contains 899 molecules, all of which consist of some of C, O, N, S, Cl and H. 

