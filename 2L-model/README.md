---
title: "Readme for the 2L-model Package of Project Mol-Infer"
date: "March 20, 2021"
author: "Discrete Mathematics Lab, Kyoto University"
---

<p align="center">
  <a href="/2L-model/README.md">English</a>
  ·
  <a href="/2L-model/README_jp.md">日本語</a>
</p>

**STATUS**

Todo:

+ English manual of Modules 1 and 2
+ Japanese manual of Module 3

We have carefully prepared this repository. If one finds bugs or mistakes, please contact us so that we can fix them
in the next release. Thank you.

---

## General info: mol-infer: Molecular Infering

Mol-infer is a project developed by the Discrete Mathematics Lab at Kyoto Univerisity (ku-dml). See [the top page](https://github.com/ku-dml/mol-infer) for more detail.

## Introduction of the 2L-model Package

This package consists of four modules.

+ Module 1 calculates descriptors. See [Module 1](Module_1/) for detail.
+ Module 2 constructs an *Artificial Neural Network* (ANN) that learns from known chemical compounds (given by FVs) and their properties. Thus this ANN can be used to infer the property of a given chemical compound. See [Module 2](Module_2/) for detail.
+ Module 3 implements a *Mixed-Integer Linear Programming* (MILP) that solves the inverse ANN problem.
[Module 3](Module_3/) for detail.
+ Module 4 generates graphs (partial enumeration). See [Module 4](Module_4/) for detail.

In order to understand how they deal with these tasks, one may need to read our [paper](https://doi.org/10.3390/ijms22062847). There is also an illustration of the flow in Japanese: see [JPEG file](illustration.jpg) or [PDF file](flow_jp.pdf).

## Compile and Usage

See the user's manual in each module please.

## How Input and Output Files Are Used
![how input and output files are used](illustration.jpg)

