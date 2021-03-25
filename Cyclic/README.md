---
title: "Readme for the Cyclic Package of Project Mol-Infer"
date: "Feb 15, 2020"
author: "Discrete Mathematics Lab, Kyoto University"
---

**STATUS**

We have carefully prepared this repository. If one finds bugs or mistakes, please contact us so that we can fix them
in the next release. Thank you.

---

## General info: mol-infer: Molecular Infering

Mol-infer is a project developed by the Discrete Mathematics Lab at Kyoto Univerisity (ku-dml). See [the top page](https://github.com/ku-dml/mol-infer) for more detail.

## Introduction of the Cyclic Package

This package consists of four modules.

+ Module 1 calculates descriptors. See [Module 1](Module_1/) for detail.
+ Module 2 constructs an *Artificial Neural Network* (ANN) that learns from known chemical compounds (given by FVs) and their properties. Thus this ANN can be used to infer the property of a given chemical compound. See [Module 2](Module_2/) for detail.
+ Module 3 implements a *Mixed-Integer Linear Programming* (MILP) that solves the inverse ANN problem.
[Module 3](Module_3/) for detail.
+ Module 4 generates graphs (partial enumeration). See [Module 4](Module_4/) for detail.

In order to understand how they deal with these tasks, one may need to read our [paper](https://arxiv.org/abs/2010.09203).

## Compile and Usage

For quick instructions on how to get started please check the [Graphical Instruction](cyclic_flow.pdf).
For more detailed information, please see the user's manual in each module.
