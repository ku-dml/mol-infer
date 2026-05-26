# Module 1

Module 1 generates feature vectors for a given polymer dataset provided as an SDF file.

The script used to create the mixed-vector representation is dataset-specific. For the datasets used in the paper, please refer to [`datasets.zip`](../instances_for_paper/datasets.zip).

- **Generate Descriptors**

  Generates the descriptor file for a given **polymer** dataset described in the **monomer representation** introduced in Section A.3.2. Please refer to the [`Preprocessors`](./Preprocessors) folder for instructions on how to generate this representation.

- **Generate Quadratic Descriptors**

  Generates the _quadratic descriptor_ file for a given dataset. These descriptors are used in the R-MLR learning method.

- **Preprocessors**

  The [`Preprocessors`](./Preprocessors) folder contains:
  
  - two scripts used to preprocess your SDF file;
  - two scripts used to convert a file representing a **polymer** into the **monomer representation** described in Section A.3.2.

  Note that the SDF files used in the paper and provided in [`instances_for_paper`](../instances_for_paper) have already been processed by these scripts. Therefore, you do not need to run these scripts for the SDF files in that directory.
  