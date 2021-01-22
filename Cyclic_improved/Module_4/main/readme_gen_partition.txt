The code gen_partition.cpp is a code used to generate standard partition file given a SDF file. 
Here standard partition means, we consider each vertex in the 
polymer topology of a given graph as base vertex, where
a self-loop will be partitioned into two e-components.


Sample command to run the code is:

g++ -o gen_partition gen_partition.cpp -O3 -std=c++11
./gen_partition instance_1.sdf instance_1_partition.txt
   (File name of input graph)  (File name of output partition file)

A sample partition file is:
9                     // number of base vertices
21 # C                // the index of base vertex in SDF, here the atom is shown after #
0 0 0                 // first two numbers represent range of core height, for the third number, 0 means this v-component can be changed, 1 means this v-component is fixed. 
22 # C
0 0 0
23 # C
0 0 0
24 # C
0 1 1
26 # C
0 0 1
27 # C
0 0 1
30 # C
0 0 1
31 # C
0 0 1
1 # C
0 1 1
13                   // number of base edges
21 19 20 15 22 # C1C1N1C1C // the indices of base edges in SDF in a certain ordering that adjacent two indices represent an edge in SDF, here the path is shown after #
0 2 1                // first two numbers represent range of core height, for the third number, 0 means this e-component can be changed, 1 means this e-component is fixed. 
21 22 # C1C
0 0 0
21 24 # C1C
0 0 0
22 23 # C1C
0 0 0
23 13 14 27 # C1C1N1C
0 0 0
23 26 # C1C
0 0 0
24 26 # C1C
0 0 0
24 28 30 # C1C1C
0 0 0
26 27 # C2C
0 0 0
27 29 30 # C1C1C
0 0 0
30 10 11 12 31 # C2C1C2C1C
0 4 0
31 1 # C1C
0 0 1
1 3 4 32 9 8 31 # C1N1C1C2C1N1C
0 3 1
