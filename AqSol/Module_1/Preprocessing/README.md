# Preprocessing
This directory contains two preprocessing codes that are as follow:
- Ex_SDF_NS.py
- eli_v062.py

#### Ex_SDF_NS.py
All datasets collected are in .csv file  where containing SMILES, NAME or SLN of compounds along with their logS values.  Ex_SDF_NS.py generates an SDF for given compounds by using pubchempy library.
Let input.csv be your original data file. To extract sdf from pubchem database,  
```sh
$ python Ex_SDF_NS.py input.csv 
```
Executing this, you will have 
- input.sdf (required for eli_v062.py)
- input_values.csv (containing CID's and logS values for compounds avaibale on pubchem)
- input_norm_values.txt.(conatining CID's and normlaized logS values as this file will be required in Module 2)

#### eli_v062.py 
Our framework doesn't work with chemical compounds that have special chemical structures, like those with fewer than four carbon atoms or disconnected structures. The code eli_v062.py creates an SDF file containing only the valid chemical compounds from a given SDF file.
Let Input.sdf be your SDF file. To run eli_v062.py to generate input_eli.sdf,
```sh
$ python eli_v062.py Input.sdf
```
You can use an example instance D5.sdf in this directory.
```sh
$ python eli_v062.py D5.sdf
```
Executing this, you will have D5_eli.sdf. Also you can see that D5_eli.sdf contains 91 compounds, whereas D5.sdf contains 93 compounds.