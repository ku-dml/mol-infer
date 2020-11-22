---
title: "Readme for the Mol-Infer project"
date: "Nov 22, 2020"
author: "Discrete Mathematics Lab, Kyoto University"
---

Notice: at the point of writing, this repository has not been fully prepared. Please visit each package for detail.

# mol-infer: Molecular Infering

Mol-infer is a project developed by the Discrete Mathematics Lab at Kyoto Univerisity (ku-dml).
After many years research on original graph algorithms for infering molecular,
we decided to open-source our programs for public use.
If you found it was useful in your research, please consider to cite our following paper(s) as well as this GitHub repository.

> (Acyclic package) N.A. Azam, J. Zhu, Y. Sun, Y. Shi, A. Shurbevski, L. Zhao, H. Nagamochi and T. Akutsu, A Novel Method for Inference of Acyclic Chemical Compounds with Bounded Branch-height Based on Artificial Neural Networks and Integer Programming, 2020, https://arxiv.org/abs/2009.09646.

> (Cyclic package) TBA

## Introduction

This project consists of two packages: one for [Acyclic](acyclic graphs) and another for [Cycle](cyclic graphs). They has similar algorithm structure, thus may have modules with the same names. However, their modules are NOT compatible since they use different algorithms. So you should think these two packges are two different projects.

## Requirement

A standard C++ compiler and Python with some standard packages. See each package for detail please.

## Acknowledgement

This project is partially supported by JSPS Grant (KAKENHI) 18H04113.
