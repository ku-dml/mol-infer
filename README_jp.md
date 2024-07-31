---
title: "Mol-Inferプロジェクトについて"
date: "2024年7月9日"
author: "京都大学　情報学研究科　離散数理分野研究室"
---

<p align="center">
  <a href="/README.md">English</a>
  ·
  <a href="/README_jp.md">日本語</a>
</p>

# mol-infer: 化合物構造推定

Mol-inferプロジェクトは京都大学の離散数理分野研究室（ku-dml）が開発したプロジェクトです．
多年にわたる分子構造推定のためのグラフアルゴリズムの研究の後，我々はこれらのプログラムをオープンソース化しました．
このプロジェクトがあなたの研究に役立った場合は，私たちの論文とこのGitHubリポジトリを引用していただけると幸いです．

## パッケージの概要

プログラムはそれぞれのパッケージにサブフォルダを割り当て，アップロードしています．
全てのパッケージは似たよったアルゴリズム構造を持ち4つのモジュールで構成されていますが，
それぞれのモジュールは異なるアルゴリズムを使用しているため，**互換性はありません**．
パッケージごとに異なるプロジェクトとして考えてください．

**注意:** 詳細については，各パッケージをご覧ください．
- 何個かのパッケージは完全には準備されていない場合があります．
- 構造的な仮定に加えて，入力化学グラフについて他の仮定を行うことがあります．例えば，[Cyclic](Cyclic/)と[Cyclic_improved](Cyclic_improved/)では，グラフは*2-lean*である必要があります．

### Package [Polymer](Polymer/) (May 2024)
- **入力グラフ:** 任意のグラフ（i.e.，cyclicとacyclic両方のグラフを同時に扱うことができます）
- **参考文献:**
  - R. Ido, S. Cao, J. Zhu, N. A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu, A Method for Inferring Polymers Based on Linear Regression and Integer Programming, arXiv 2109.02628, 2021, https://arxiv.org/abs/2109.02628.

### Package [HPS](HPS/) (Mar 2024)
- **入力グラフ:** 任意のグラフ（i.e.，cyclicとacyclic両方のグラフを同時に扱うことができます）
- **参考文献:**
  - J. Zhu, N. A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu, Molecular Design Based on Integer Programming and Splitting Data Sets by Hyperplanes, arXiv 2305.00801, 2023, https://arxiv.org/abs/2305.00801.

### Package [RMLRQ](RMLRQ/) (May 2023)
- **入力グラフ:** 任意のグラフ（i.e.，cyclicとacyclic両方のグラフを同時に扱うことができます）
- **参考文献:**
  - J. Zhu, N. A. Azam, S. Cao, R. Ido, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu, Molecular Design Based on Integer Programming and Quadratic Descriptors in a Two-layered Model, arXiv 2209.13527, 2022, https://arxiv.org/abs/2209.13527.

### Package [Grid-neighbor-search](Grid-neighbor-search/) (Dec 2021)
- **入力グラフ:** 任意のグラフ（i.e.，cyclicとacyclic両方のグラフを同時に扱うことができます）
- **参考文献:**
  - N. A. Azam, J. Zhu, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu. Molecular Design Based on Artificial Neural Networks, Integer Programming and Grid Neighbor Search, *Proceedings of IEEE International Conference on Bioinformatics and Biomedicine (BIBM)*, 2021.

### Package [ALR](ALR/) (Sep 2021)
- **入力グラフ:** 任意のグラフ（i.e.，cyclicとacyclic両方のグラフを同時に扱うことができます）
- **参考文献:**
  - J. Zhu, K. Haraguchi, H. Nagamochi and T. Akutsu: Adjustive Linear Regression and Its Application to the Inverse QSAR, *Proceedings of the 15th International Joint Conference on Biomedical Engineering Systems and Technologies (BIOSTEC 2022) - Volume 3: BIOINFORMATICS*, 2021, https://doi.org/10.5220/0010853700003123.

### Package [2LMM-LLR](2LMM-LLR/) (Jul 2021)
- **入力グラフ:** 任意のグラフ（i.e.，cyclicとacyclic両方のグラフを同時に扱うことができます）
- **参考文献:** 
  - J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Inverse QSAR Method Based on Linear Regression and Integer Programming, 2021, http://arxiv.org/abs/2107.02381.

### Package [2L-model](2L-model/) (Mar 2021)
- **入力グラフ:** 任意のグラフ（i.e.，cyclicとacyclic両方のグラフを同時に扱うことができます）
- **参考文献:**
  - Y. Shi, J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Inverse QSAR Method Based on a Two-Layered Model and Integer Programming, *International Journal of Molecular Sciences*, **22**(6), 2021, https://doi.org/10.3390/ijms22062847. 

### Package [Cyclic_improved](Cyclic_improved/) (Jan 2021)
- **入力グラフ:** cyclicグラフ
- **参考文献:**
  - J. Zhu, N.A. Azam, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: An Improved Integer Programming Formulation for Inferring Chemical Compounds with Prescribed Topological Structures, *Proceedings of IEA/AIE 2021 conference* (https://ieaaie2021.wordpress.com), 2021, accepted.

### Package [Cyclic](Cyclic/) (Nov 2020)
- **入力グラフ:** cyclicグラフ
- **参考文献:**
  - J. Zhu, N.A. Azam, F. Zhang, A. Shurbevski, K. Haraguchi, L. Zhao, H. Nagamochi and T. Akutsu: A Novel Method for Inferring of Chemical Compounds with Prescribed Topological Substructures Based on Integer Programming, 2020, submitted. 
  - T. Akutsu and H. Nagamochi: A novel method for inference of chemical compounds with prescribed topological substructures based on integer programming, 2020, https://arxiv.org/abs/2010.09203.

### Package [Acyclic](Acyclic/) (Sep 2020)
- **入力グラフ:** cyclicを持たないグラフ（i.e.，木構造のグラフ）
- **参考文献:**
  - N.A. Azam, J. Zhu, Y. Sun, Y. Shi, A. Shurbevski, L. Zhao, H. Nagamochi and T. Akutsu, A Novel Method for Inference of Acyclic Chemical Compounds with Bounded Branch-height Based on Artificial Neural Networks and Integer Programming, 2020, https://arxiv.org/abs/2009.09646.


## 要件

標準のC++コンパイラとPythonといくつかの標準パッケージが必要です．詳細は各パッケージを参照してください．

## 謝辞

このプロジェクトは一部JSPS Grant (KAKENHI) 18H04113と22H00532に支援していただいてます．
