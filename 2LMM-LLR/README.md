---
title: "Readme for the 2LMM-LLR Package of Project Mol-Infer"
date: "July 4, 2021"
author: "Discrete Mathematics Lab, Kyoto University"
---

**STATUS**

We have carefully prepared this repository. If one finds bugs or mistakes, please contact us so that we can fix them
in the next release. Thank you.

---

## General info: mol-infer: Molecular Infering

Mol-infer is a project developed by the Discrete Mathematics Lab at Kyoto Univerisity (ku-dml). See [the top page](https://github.com/ku-dml/mol-infer) for more detail.

## Introduction of the 2LMM-LLR Package

This package consists of four modules.

+ Module 1 calculates descriptors. See [Module 1](Module_1/) for detail. 
+ Module 2 constructs a prediction function by using *Lasso Linear Regression* (LLR). See [Module 2](Module_2/) for detail.
+ Module 3 implements a *Mixed-Integer Linear Programming* (MILP) that solves the inverse problem.
[Module 3](Module_3/) for detail.
+ Module 4 generates graphs (partial enumeration). See [Module 4](Module_4/) for detail.

In order to understand how they deal with these tasks, one may need to read our [paper](https://arxiv.org/abs/2107.02381).

## Quickstart visual guide

Please check the image below on how data and files are used and passed through different modules
of the 2LMM-LLR Package.
For more details on usage and compiling, please see the user's manual in each module

![Data flow illustration](/2LMM-LLR/doc/2LMM-LLR_flow.PNG)
