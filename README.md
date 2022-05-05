---
title: "Readme for the Mol-Infer project"
date: "March 15, 2021"
author: "Discrete Mathematics Lab, Kyoto University"
---


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

### Package [Grid-neighbor-search](Grid-neighbor-search/) (Dec 2021)
- **Input graphs:** Arbitrary graphs (i.e., both cyclic and acyclic graphs can be treated at the same time)
- **Reference:** 
  - N. A. Azam, J. Zhu, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu. Molecular Design Based on Artificial Neural Networks, Integer Programming and Grid Neighbor Search, *Proceedings of IEEE International Conference on Bioinformatics and Biomedicine (BIBM)*, 2021.


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

## クイックスタート

0. Python>=3.7とCPLEXソルバーが必要である。\
   Windowsユーザーには、Cygwinをインストール必要がある。\
   そして、Cygwinをインストールする時、`make`機能が有効であることが必要である。
1. `mol-infer_config.sh`を編集し、CPLEXのパスを設定する。
2. ターミナル（WindowsユーザーはCygwin Terminal）で：
   ```
   cd mol-infer
   make
   python -m venv python-venv
   (Windows):       python-venv\Scripts\pip install -r requirements.txt
   (Linux & MacOS): python-venv/bin/pip install -r requirements.txt
   mkdir experiments
   cd experiments
   bash ../mol-infer.sh
   ```

ターミナルに指示が表示されるはず。\
各パッケージには2つのオプションが含まれている：`train`と`infer`。
- `train`\
  分子の特徴ベクトルを計算し、ANNやLLRを学習させる（訓練する）。
- `infer`\
  `train`の出力を入力とする。\
  必要なプロパティを持つ分子を生成する。

計算例を実行するには、すべての入力を空欄にしてenterキーを押してください。デフォルトの値が使用されます。\
入力ファイルの準備や適切なパラメータの選択方法については、各パッケージのドキュメントを参照してください。

## Quick Start

0. Python>=3.7 and a CPLEX solver are required.\
   Fow Windows users: Cygwin installed with `make` enabled is required.
1. Edit `mol-infer_config.sh`, set CPLEX path.
2. In terminal (Use Cygwin Terminal if you are using Windows):
   ```
   cd mol-infer
   make
   python -m venv python-venv
   (Windows):       python-venv\Scripts\pip install -r requirements.txt
   (Linux & MacOS): python-venv/bin/pip install -r requirements.txt
   mkdir experiments
   cd experiments
   bash ../mol-infer.sh
   ```

Instructions should appear on your terminal.\
Each package includes two options: `train` and `infer`.
- `train`\
  Calculate descriptors for molecules and train an ANN or LLR.
- `infer`\
  Takes the output of `train` as input.\
  Generate molecules with desired properties.

To run examples, leave every input blank and press enter. Default values will be used.\
Please see the documents in each package for details about how to prepare input files and choose the proper parameters.

## Acknowledgement

This project is partially supported by JSPS Grant (KAKENHI) 18H04113.
