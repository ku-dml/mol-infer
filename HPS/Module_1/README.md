# Module 1

Module 1 consists of the stage of generating the feature vectors (defined in Section 3, paragraph "Feature Function") of a given data set (sdf format file).

- Generate Linear Descriptors:
  
  Generate the _linear descriptors_ file of a given data set.
  
- Generate Quadratic Desctiptors:
  
  Generate the _quadratic descriptors_ file of a given data set that is used in the learning method RLR.

- Preprocessors:
  
  Scripts that are used for preprocessing _your_ SDF files.
  - **eliminate.py** builds an _SDF_ file that is the subset of feasible molecules in a given SDF file.
  - **limit_atoms.py** builds an _SDF_ file that is the subset of molecules containing only specified atoms in a given SDF file. 
  - Note that [SDF files that are used in the paper](../instances_for_paper) are already processed by these two scripts. You don't need to run the scripts for the SDF files in the above directory. 
