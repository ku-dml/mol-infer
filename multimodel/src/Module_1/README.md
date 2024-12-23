<p align="center">
  <a href="/multimodel/src/Module_1/README.md">English</a>
  ·
  <a href="/multimodel/src/Module_1/README_jp.md">日本語</a>
</p>

# Module 1: Feature vector calculation

Module 1 consists of the stage of generating the feature vectors of a given data set (sdf format file).

# Preparation
## Compile fv_2LMM.cpp
```bash
cd src/Module_1
g++ -O2 -Wall -std=c++20 -o FV_2LMM_V019 fv_2LMM.cpp
```

## Data preparation
Delete unnecessary data with eliminate.py
```bash
python eliminate.py (data.sdf)
```

sample
```bash
python eliminate.py sample_instance/input/Bp_small.sdf
```

(optional)
Limit the atomic species with limit_atoms.py
```bash
python limit_atoms.py (data.sdf) (atomic species)
```

sample
```bash
python limit_atoms.py sample_instance/input/Bp_small.sdf C O N S Cl H
```

# Run the program
## Calculate feature vectors

```bash
./FV_2LMM_V019 (data.sdf) (output_prefix)
```

Sample script
```bash
./FV_2LMM_V019 sample_instance/input/Bp_small.sdf sample_instance/output/Bp_small 
```
