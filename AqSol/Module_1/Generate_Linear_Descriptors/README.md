# Generate Linear Descriptors
This folder conatains code to generate the linear descriptors defined in section 3.1.1 which will be used for the learning stage in Module 2.
First compile the c++ files as following on terminal:

This will generate an executable file named fv_2LMM, and then generate the linear descriptors files by using following command:
```sh
./fv_2LMM input.sdf OUTPUT
```
where,
- input.sdf: the input sdf file (after elimmination), and 
- OUTPUT: the prefix for output files

This code will generate following files:
- **OUTPUT_desc.csv:** The csv file of the linear descriptors of the dataset before normalization.
- **OUTPUT_desc_norm.csv:** The csv file of the linear descriptors of the dataset after normalization.
- **OUTPUT_fringe.txt:** The file containing the fringe tree information of the dataset, which will be used in the stage of solving inverse problem (Module 3).

A sample usage:
```sh
./fv_2LMM ./Data/D5_eli.sdf ./Data/D5
```