To compute isomers of a given graph.
The output graphs will have exactly same ch as that of input graph.

an example of commands to run the program:

g++ -o main main.cpp -O3 -std=c++11
./main instance_1.sdf 3600 10000000 5 100 output.sdf instance_1.txt
   (File name of input graph)  (time limit) (vector size limit) (number of sample trees per vector) (number of output graphs) (file name of output SDF files) (file name of information of the partition)