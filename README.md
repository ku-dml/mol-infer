---
title: "Readme for the Mol-Infer project"
date: "March 15, 2021"
author: "Discrete Mathematics Lab, Kyoto University"
---

Notice: at the point of writing, this repository hosts research work that is under progress.
Some parts are under development and will be updated in the future.
Please visit each package for detail.

# mol-infer: Molecular Infering

Mol-infer is a project developed by the Discrete Mathematics Lab at Kyoto Univerisity (ku-dml).
After many years of research on original graph algorithms for infering molecular structures,
we decided to open-source our programs for public use.
If you found it was useful in your research, please consider to cite our paper(s) as well as this GitHub repository.

> (Acyclic package) N.A. Azam, J. Zhu, Y. Sun, Y. Shi, A. Shurbevski, L. Zhao, H. Nagamochi and T. Akutsu, A Novel Method for Inference of Acyclic Chemical Compounds with Bounded Branch-height Based on Artificial Neural Networks and Integer Programming, 2020, https://arxiv.org/abs/2009.09646.

> (Cyclic package) T. Akutsu and H. Nagamochi, A novel method for inference of chemical compounds with prescribed topological substructures based on integer programming, 2020, https://arxiv.org/abs/2010.09203.

> (Cyclic_improved package) J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu, An Improved Integer Programming Formulation for Inferring Chemical Compounds with Prescribed Topological Structures, 2021, submitted.

> (2L-model package) Y. Shi, J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi, T. Akutsu, A Two-layered Model for Inferring Chemical Compounds with Integer Programming, 2021, submitted.

## Introduction

This project consists of four packages:
+ [Acyclic package for graphs with no cycle](Acyclic/);
+ [Cyclic package for graphs with cycle(s)](Cyclic/);
+ [Cyclic improved package for graphs with cycle(s)](Cyclic_improved/); and
+ [Two-layered model package](2L-model/).

All packages have a similar algorithmic structure, thus have modules with the same names.
However, their modules are NOT compatible, since they use different algorithms.
You should think of different packages as different projects.

## Requirement

A standard C++ compiler and Python with some standard packages. See each package for detail please.

## Acknowledgement

This project is partially supported by JSPS Grant (KAKENHI) 18H04113.
