---
title: "Readme for the Mol-Infer project"
date: "March 15, 2021"
author: "Discrete Mathematics Lab, Kyoto University"
---

<p align="center">
  <a href="/README.md">English</a>
  ·
  <a href="/README_jp.md">日本語</a>
</p>

# mol-infer: Molecular Infering

Mol-infer is a project developed by the Discrete Mathematics Lab at Kyoto Univerisity (ku-dml).
After many years of research on original graph algorithms for infering molecular structures,
we decided to open-source our programs for public use.
If you found it was useful in your work,
please consider citing our paper(s) as well as this GitHub repository.

## Overview of packages

We have uploaded our programs in the following packages, where each package is assigned one subfolder. 
All packages have a similar algorithmic structure and consist of four modules.
However, modules are **NOT** compatible between packages, since they use different algorithms.
You should think of different packages as different projects.


**Notice:** Please visit each package for details.
- Some packages may not be fully prepared. 
- In addition to structural assumptions, we may make other assumptions on input chemical graphs. For example, in [Cyclic](Cyclic/) and [Cyclic_improved](Cyclic_improved/), the graphs should be *2-lean* as well as cyclic.

### Package [Polymer](Polymer/) (May 2024)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:**
  - R. Ido, S. Cao, J. Zhu, N. A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu, A Method for Inferring Polymers Based on Linear Regression and Integer Programming, arXiv 2109.02628, 2021, https://arxiv.org/abs/2109.02628.

### Package [HPS](HPS/) (Mar 2024)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:**
  - J. Zhu, N. A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu, Molecular Design Based on Integer Programming and Splitting Data Sets by Hyperplanes, arXiv 2305.00801, 2023, https://arxiv.org/abs/2305.00801.

### Package [RMLRQ](RMLRQ/) (May 2023)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:**
  - J. Zhu, N. A. Azam, S. Cao, R. Ido, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu, Molecular Design Based on Integer Programming and Quadratic Descriptors in a Two-layered Model, arXiv 2209.13527, 2022, https://arxiv.org/abs/2209.13527.

### Package [Grid-neighbor-search](Grid-neighbor-search/) (Dec 2021)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:** 
  - N. A. Azam, J. Zhu, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu. Molecular Design Based on Artificial Neural Networks, Integer Programming and Grid Neighbor Search, *Proceedings of IEEE International Conference on Bioinformatics and Biomedicine (BIBM)*, 2021.

### Package [ALR](ALR/) (Sep 2021)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:** 
  - J. Zhu, K. Haraguchi, H. Nagamochi and T. Akutsu: Adjustive Linear Regression and Its Application to the Inverse QSAR, *Proceedings of the 15th International Joint Conference on Biomedical Engineering Systems and Technologies (BIOSTEC 2022) - Volume 3: BIOINFORMATICS*, 2021, https://doi.org/10.5220/0010853700003123.

### Package [2LMM-LLR](2LMM-LLR/) (Jul 2021)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:** 
  - J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Inverse QSAR Method Based on Linear Regression and Integer Programming, 2021, http://arxiv.org/abs/2107.02381.

### Package [2L-model](2L-model/) (Mar 2021)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:**
  - Y. Shi, J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Inverse QSAR Method Based on a Two-Layered Model and Integer Programming, *International Journal of Molecular Sciences*, **22**(6), 2021, https://doi.org/10.3390/ijms22062847. 

### Package [Cyclic_improved](Cyclic_improved/) (Jan 2021)
- **Input graphs:** Cyclic graphs.
- **Reference:**
  - J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Improved Integer Programming Formulation for Inferring Chemical Compounds with Prescribed Topological Structures, *Proceedings of IEA/AIE 2021 conference* (https://ieaaie2021.wordpress.com), 2021, accepted.

### Package [Cyclic](Cyclic/) (Nov 2020)
- **Input graphs:** Cyclic graphs.
- **References:**
  - J. Zhu, N.A. Azam, F. Zhang, A. Shurbevski, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: A Novel Method for Inferring of Chemical Compounds with Prescribed Topological Substructures Based on Integer Programming, 2020, submitted. 
  - T. Akutsu and H. Nagamochi: A novel method for inference of chemical compounds with prescribed topological substructures based on integer programming, 2020, https://arxiv.org/abs/2010.09203.

### Package [Acyclic](Acyclic/) (Sep 2020)
- **Input graphs:** Graphs with no cycle (i.e., tree structured graphs)
- **Reference:**
  - N.A. Azam, J. Zhu, Y. Sun, Y. Shi, A. Shurbevski, L. Zhao, H. Nagamochi and T. Akutsu, A Novel Method for Inference of Acyclic Chemical Compounds with Bounded Branch-height Based on Artificial Neural Networks and Integer Programming, 2020, https://arxiv.org/abs/2009.09646.


## Requirement

A standard C++ compiler and Python with some standard packages. See each package for detail please.

## Acknowledgement

This project is partially supported by JSPS Grant (KAKENHI) 18H04113 and 22H00532.
