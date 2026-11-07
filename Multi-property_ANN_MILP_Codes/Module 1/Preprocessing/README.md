# Preprocessing



This directory contains the preprocessing codes used to prepare the dataset for descriptor generation and the learning stage.

The preprocessing is performed in the following order:



1. Select a random subset from the original CSV dataset.
2. Normalize the property values of the selected subset.
3. Generate the SDF records of the selected compounds from their SMILES strings using RDKit.
4. Eliminate unsupported or invalid chemical structures.
5. Select the required number of compounds from the eliminated SDF file.



The preprocessing scripts are:



* random\_sample.py
* normalized-full-qm9.py
* RDKit-SDF-extraction.py
* eli\_v062.py
* First\_n\_SDF\_selection.py

###### Random Subset Selection



The original dataset is provided as a CSV file containing SMILES strings and the corresponding molecular property values. The random\_sample.py script randomly selects the required number of compounds from the original CSV dataset.

Run:

python random\_sample.py input.csv number\_of\_samples



A sample usage:

python random\_sample.py qm9.csv 11000

Executing this command will generate:



qm9\_seed\_75\_11000.csv



###### Normalization of Property Values



The normalized-full-qm9.py script applies Min–Max scaling to every numeric column of the sampled CSV file, so each property value lies between 0 and 1. Non-numeric columns, such as the SMILES strings, are kept unchanged. The result is saved as a comma-separated TXT file.

Run:



python normalized-full-qm9.py input.csv \[output.txt]



The output file name is optional. If it is not given, the output is saved next to the input file with the same name followed by \_norm\_values.txt. Paths that contain spaces must be written in double quotes.

A sample usage:



python normalized-full-qm9.py qm9\_seed\_75\_11000.csv



Executing this command will generate:



qm9\_seed\_75\_11000\_norm\_values.txt



This file contains the normalized property values and is used as the normalized target file in the learning stage. The script prints the names of the columns it scaled, so check this list to confirm that only property columns were normalized.



###### RDKit SDF Extraction



The RDKit-SDF-extraction.py script reads the SMILES strings from the CID column of the sampled CSV file and uses RDKit to generate a 3D structure for each compound. It adds hydrogens, embeds the molecule with ETKDG (retrying with random coordinates and then a 2D-to-3D fallback if needed), optimizes it with UFF where possible, and writes all structures to an SDF file. The SMILES string is stored as the name of each SDF record. Invalid SMILES are skipped and reported.

The input CSV file and the output SDF file are given on the command line. Paths that contain spaces must be written in double quotes.

Run:



python RDKit-SDF-extraction.py input.csv output.sdf



A sample usage:



python RDKit-SDF-extraction.py qm9\_seed\_75\_11000.csv qm9\_seed\_75\_11000.sdf



Executing this command will generate:



qm9\_seed\_75\_11000.sdf



The script prints the total number of compounds and the number successfully written to the SDF file.



###### Elimination of Compounds



The framework does not support compounds containing certain special chemical structures, such as disconnected structures or compounds that do not satisfy the required structural conditions.

The eli\_v062.py script removes these compounds and generates an SDF file containing only the valid compounds.

Run:



python eli\_v062.py input.SDF



A sample usage:

python eli\_v062.py qm9\_seed\_75\_11000.sdf



Executing this command will generate:

qm9\_seed\_75\_11000\_eli.sdf



###### Selection of the Required Number of Compounds



The First\_n\_SDF\_selection.py script selects the first required number of compounds from the eliminated SDF file.

Run:



python First\_n\_SDF\_selection.py input\_eli.SDF number\_of\_compounds



A sample usage for selecting 500 compounds:



python First\_n\_SDF\_selection.py qm9\_seed\_75\_11000\_eli.sdf 500



Executing this command will generate:



qm9\_seed\_75\_11000\_eli\_first\_500.sdf



This file can be used in the next stage for generating molecular descriptors.

