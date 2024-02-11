** UNDER CONSTRUCTION (2024.2.11) **

## Generating Recombination Solutions

This folder contains the code to generate chemical isomers from the output of MILP by using a dynamic programming algorithm.


Usage:

First compile the 'c++' files as following on terminal:

```
g++ -o ./main/generate_isomers ./main/generate_isomers.cpp -O3 -std=c++11
```

This will generate an executable file './main/generate_isomers', and then generate the chemical isomers by using command:

```
./main/generate_isomers 
```

Here:



A sample usage:

```
./main/generate_isomers ./sample_instance/sample.sdf 10 1000000 0 300 1000000 100 ./sample_instance/sample_output.sdf ./sample_instance/sample_partition.txt ./sample_instance/sample_fringe_tree.txt
```
