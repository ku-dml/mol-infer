---
title: "Readme for the Mol-Infer project"
date: "Mar 14, 2021"
author: "Discrete Mathematics Lab, Kyoto University"
<STYLE TYPE="text/css">
h1{
	background-color: #d1f1ff;
}
h2{
	background-color: #ddf8ff;
}
h3{
	background-color: #eeffff;
}
body{
	font-size: 20pt;
}</STYLE>
---

# mol-infer: Molecular Infering

Mol-infer is a project developed by the Discrete Mathematics Lab at Kyoto Univerisity (ku-dml).
After many years research on original graph algorithms for infering molecular,
we decided to open-source our programs for public use.
If you found it was useful in your research,
please consider to cite our paper(s) as well as this GitHub repository.


## Overview of packages

We have uploaded our programs on the following packages, where each package is assigned one subfolder. All packages have similar algorithm structure and consist of four modules. However, their modules are **NOT** compatible between packages since they use different algorithms. So you should think any two packages are two different projects.


**Notice:** Please visit each package for detail.
- Some packages may not be fully prepared. 
- In addition to structural assumptions, we may make other assumptions on input chemical graphs. For example, in [Cyclic](Cyclic/) and [Cyclic_improved](Cyclic_improved/), the graphs should be *2-lean* as well as cyclic. 


### Package [2L-model](2L-model/) (Mar 2021)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:**
  - Y. Shi, J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Inverse QSAR Method Based on a Two-Layered Model and Integer Programming, *International Journal of Molecular Sciences*, **22**(6), 2021, https://doi.org/10.3390/ijms22062847. 

### Package [Cyclic_improved](Cyclic_improved/) (Jan 2021)
- **Input graphs:** Cyclic graphs.
- **Reference:**
  - J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Improved Integer Programming Formulation for Inferring Chemical Compounds with Prescribed Topological Structures, Proceedings of IEA/AIE 2021 conference, 2021, accepted.

### Package [Cyclic](Cyclic/) (Nov 2020)
- **Input graphs:** Cyclic graphs.
- **References:**
  - J. Zhu, N.A. Azam, F. Zhang, A. Shurbevski, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: A Novel Method for Inferring of Chemical Compounds with Prescribed Topological Substructures Based on Integer Programming, 2020, submitted. 
  - T. Akutsu and H. Nagamochi: A novel method for inference of chemical compounds with prescribed topological substructures based on integer programming, 2020, https://arxiv.org/abs/2010.09203.

### Package [Acyclic](Acyclic/) (Sep 2020)
- **Input graphs:** Graphs with no cycle (i.e., tree structured graphs)
- **References:** N.A. Azam, J. Zhu, Y. Sun, Y. Shi, A. Shurbevski, L. Zhao, H. Nagamochi and T. Akutsu, A Novel Method for Inference of Acyclic Chemical Compounds with Bounded Branch-height Based on Artificial Neural Networks and Integer Programming, 2020, https://arxiv.org/abs/2009.09646.


## Requirement

A standard C++ compiler and Python with some standard packages. See each package for detail please.

## Acknowledgement

This project is partially supported by JSPS Grant (KAKENHI) 18H04113.
