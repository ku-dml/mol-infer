## Combine Descriptors

This folder contains the code to combine descriptors generated in [Generate Descritpors](/chi-parameter/Module_1/Generate_Descriptors) as decribed in Section 3, which will be used for the learning stage in [Module 2](/chi-parameter/Module_2).

Here the code is designed particularly for the case of $\chi$-paramter. The temperature information is contained in the name of the molecule like "XXX_273.0" (the temperature is 273.0K).

Usage:

```
python combine_chi_desc_T.py SOLUTE.csv SOLVENT.csv OUTPUT -aT (k1) (k2) ...
```

Here:
- SOLUTE.csv: the descriptor file of the solutes;
- SOLVENT.csv: the descriptor file of the solvents;
- OUTPUT: the prefix for the output files, and
- -aT (k1) (k2) ...: add $T^{k1}$, $T^{k2}$ ... in the feature vector.
  
When the generator finishes normally, it will generate the following files:
- OUTPUT_desc.csv: The csv file of the combined descriptors of the data set **BEFORE** normalization.
- OUTPUT_desc_norm.csv: The csv file of the combined descriptors of the data set **AFTER** normalization.

A sample usage:

```
python combine_chi_desc_T.py ./sample_instance/Chi_Aoki_ch1_desc.csv ./sample_instance/Chi_Aoki_ch2_desc.csv -aT -1
```

