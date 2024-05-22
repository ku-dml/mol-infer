# Preprocessors
The directory contains several scripts that are used to preprocess your SDF file.

Note that [SDF files that are used in the paper](../../instances_for_paper) are already processed by these scripts. You don't need to run the scripts for the SDF files in the above directory. 

The scripts are: 

- **polymer_converter_mol_to_sdf** (converting mol file for polymer to SDF file)
- **eliminate.py** (strongly recommended to run on your SDF before starting Module 1)
- **limit_atoms.py** (rather optional)
- **contract_e.cpp** (converting polymer to its monomer-representation)

## polymer_converter_mol_to_sdf.py
A file to convert the mol file representing polymers to the SDF format.

Suppose that the mol files are in the folder ./sample_instance/mol, run **polymer_converter_mol_to_sdf.py** to get a SDF file. We assume that the polymer data set is represented in the format like the ones in that folder.

```
$ python polymer_converter_mol_to_sdf.py ./sample_instance/mol ./sample_instance/sample_e*.sdf
```

Executing this, you will have _sample_e*.sdf_, which will be used for the following preprocessing procedure.


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
$ python eliminate.py ./sample_instance/sample_e*.sdf ./sample_instance/sample_eli_e*.sdf 
```
Executing this, you will have _sample_eli_e*.sdf_. You will also see that _sample_eli_e*.sdf_ contains 108 molecules, where _sample_e*.sdf_ contains 127 molecules. 

## limit_atoms.py
For a given SDF, you can extract molecules that do not contain any other atoms than specified ones, as we do in the experiments in Section 4 of the paper. If you want to extract molecules in B.sdf that do not contain any other atoms than C, N, O and H, run
```
$ python limit_atoms.py B.sdf C N O H
```
Then you will have an SDF file named _B_C_N_O_H.sdf_. Each molecule in this SDF consists of some of C, N, O and H and does not contain any other atoms. The order of atoms is not essential but affects the output file name. Of course this _B_C_N_O_H.sdf_ can be also used as an input of the [feature vector generator](../Generate_Descriptors) after converting into the monomer-representation. 

For example, 
```
$ python limit_atoms.py ./sample_instance/sample_eli_e*.sdf C O N H e*
```
Then you will have _sample_eli_C_O_N_H_e*.sdf_ that contains 86 molecules, all of which consist of some of C, O, N, H and e* (a dummy atom for polymer). 


##  contract_e.cpp
Finally, you need to convert the polymer to the **monomer-representation** described in Section 2.

First, compile the `c++` files as follows on the terminal:
```
g++ -o contract_e.o contract_e.cpp -O3 -std=c++11
```

This will generate an executable file named `contract_e.o`,
and then use the following command:
```
./contract_e.o ./sample_instance/sample_eli_C_O_N_H_e*.sdf ./sample_instance/sample_eli_C_O_N_H_contract.sdf
```

This will generate a **monomer-representation** file _sample_eli_C_O_N_H_contract.sdf_ for the polymers that can be directly used as an input of the [feature vector generator](../Generate_Descriptors).
