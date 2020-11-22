# Manual for Module 4

*todo: full document in Markdown*

Last-updated: September 11th, 2020

See [Module4_manual_en.pdf](Module4_manual_en.pdf), or [Module4_manual_jp.pdf](Module4_manual_jp.pdf) for Japanese.

Compile
```
g++ -std=c++11 -Wall -O3 -o 2-branches/2-branches 2-branches/2-branches.cpp 
g++ -std=c++11 -Wall -O3 -o 3-branches/3-branches 3-branches/3-branches.cpp 
```

Here is an example of commands to run the program:

```
./2-branches 4logKow_tv9.0_n50_dia25_k2_dmax4_bl2_bh1_solver1.txt 3600 10000000 100  output.sdf  4logKow_tv9.0_n50_dia25_k2_dmax4_bl2_bh1_solver1.sdf
```
The command line arguments are the File name of input resources, the time limit, the vector size limit, the number of output graphs, the name of output SDF, and the file name of the sample sdf.
