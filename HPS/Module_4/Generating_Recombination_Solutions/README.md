** UNDER CONSTRUCTION (2024.2.19) **

## Generating Recombination Solutions

This folder contains the code to generate chemical isomers from the output of MILP by using a dynamic programming algorithm, described in Section 6.2. We refer [\[22\]](https://arxiv.org/abs/2107.02381) and [2LMM-LLR](2LMM-LLR) for more details about the details.
 

Usage:

First compile the `c++` files as following on terminal:

```
g++ -o ./main/generate_isomers ./main/generate_isomers.cpp -O3 -std=c++11
```

This will generate an executable file `./main/generate_isomers`, and then generate the chemical isomers by using command:

```
./main/generate_isomers RES.sdf 4 100000 5 10 10000 2 OUTPUT.sdf RES_partition.txt FRINGE.txt
```

Here:
- RES.sdf: the generated chemical graph in sdf format obtained in [Module 3](HPS/Module_3);
- "4": time limit for each stage of program;
- "100000": an upper bound on the number of feature vectors stored in memory;
- "5": number of sample chemical graphs stored for each feature vector;
- "10": an upper bound on time for enumeration of paths;
- "10000": an upper bound on the number of total paths stored during the computation;
- "2": an upper bound on the number of output chemical isomers;
- OUTPUT.sdf: output file name;
- RES_partition.txt: file containing the partition information of RES.sdf obtained in [Module 3](HPS/Module_3), and
- FRINGE.txt: file containing the information of available fringe trees used as part of input for [Module 3](HPS/Module_3) to obtain RES.sdf.

When the code finishes normally, it will generate the following file:
- OUTPUT.sdf: the file containing generated chemical isomers to RES.sdf under the same topological configuration described in RES_partition.txt and FRINGE.txt.


A sample usage:

```
./main/generate_isomers ./sample_instance/sample.sdf 10 1000000 0 300 1000000 100 ./sample_instance/sample_output.sdf ./sample_instance/sample_partition.txt ./sample_instance/sample_fringe_tree.txt
```
