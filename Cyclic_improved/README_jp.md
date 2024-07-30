---
title: "Mol-Inferプロジェクトにおける改良型Cyclicパッケージの使い方"
date: " 2021年2月15日"
author: "京都大学　情報学研究科　離散数理分野研究室"
---

<p align="center">
  <a href="/Cyclic_improved/README_en.md">English</a>
  ·
  <a href="/Cyclic_improved/README_jp.md">日本語</a>
</p>

**お願い**

もしバグや間違いを見つけた場合は次回のリリースで修正するのでお知らせください．

---

### Molecular Inferingについて

Mol-inferプロジェクトとは，京都大学の離散数理分野研究室（ku-dml）が開発したプロジェクトです．詳細は[トップページ](https://github.com/ku-dml/mol-infer)を参照してください.

## Cyclicパッケージの紹介

このパッケージは4つのモジュールで構成されています．
This package consists of four modules.

+ モジュール１では記述子を導出します．詳細は[Module 1](Module_1/)を参照してください．
+ Module 1 calculates descriptors. See [Module 1](Module_1/) for detail.
+ モジュール２では*Artificial Neural Network*（ANN）を構築します．このANNは与えられた化合物の物性を推論するために使用されます．詳細は[Module 2](Module_2/)を参照してください．
+ モジュール３では，モジュール２で学習されたANNと目標物性から**グラフ記述子のベクトル**（すなわちFV）を推論する*Mixed-Integer Linear Programming*（MILP）を実装します．詳細は[Module 3](Module_3/)を参照してください．
+ モジュール４ではグラフを生成します（部分列挙）．詳細は[Module 4](Module_4/)を参照してください．

モジュールがこれらのタスクをどう処理するかを理解したい場合は，[論文](https://arxiv.org/abs/2010.09203)を参照してください．

## コンパイルと使い方

それぞれのモジュールのユーザーマニュアルを参照してください．
