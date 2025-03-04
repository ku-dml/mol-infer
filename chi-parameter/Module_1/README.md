# Module 1

Module 1 consists of the stage of generating the feature vectors (defined in Sections 2 and 3) of a given polymer data set (sdf format file).

- Generate Descriptors:
  
  Generate the descriptor file of a given **polymer** data set that is described under the **monomer-representation** described in Section 2, paragraph "Polymers". Please refer the folder [Preprocessors](./Preprocessors) for how to generate this.

- Combine Descriptors:
  
  Generate the file containing the **combined** feature vector that is described in Section 3, particularly for the case of $\chi$-parameter.
  
- Preprocessors:
  - Two scripts that are used to preprocess _your_ SDF file;
  - Two scripts to converse the file representing **polymer** to the **monomer-representation** described in Section 2.3.2.
  - Note that [SDF files that are used in the paper](../instances_for_paper) are already processed by these scripts. You don't need to run the scripts for the SDF files in the above directory. 
