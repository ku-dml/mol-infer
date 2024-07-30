---
title: "Mol-Inferプロジェクトについて"
date: "2021年9月22日"
author: "京都大学　情報学研究科　離散数理分野研究室"
---

<p align="center">
  <a href="/Acyclic/README_en.md">English</a>
  ·
  <a href="/Acyclic/README_jp.md">日本語</a>
</p>

注意: 現在（2021年9月22日）このリポジトリはほぼ完成しています．ドキュメントは確認中です．

# Molecular Inferingについて

Mol-inferプロジェクトは，京都大学の離散数理分野研究室（ku-dml）が開発したプロジェクトです．
多くの年月をかけて分子推論のためのオリジナルなグラフアルゴリズムに関する研究を行った後，
我々はプログラムをオープンソース化し，一般に公開することにしました．
もし研究に役立つと感じた場合は，以下の論文とこのGitHubリポジトリを引用していただけると幸いです．

> N.A. Azam, J. Zhu, Y. Sun, Y. Shi, A. Shurbevski, L. Zhao, H. Nagamochi and T. Akutsu, A Novel Method for Inference of Acyclic Chemical Compounds with Bounded Branch-height Based on Artificial Neural Networks and Integer Programming, 2020, https://arxiv.org/abs/2009.09646.

## はじめに

このプロジェクトは以下の4つのモジュールで構成されています．

![プロジェクトの概要](images/overview.png)

+ モジュール１では化合物の標準的な*Structure Data File*（SDF）で与えられた化合物に対して*特徴ベクトル*（FV）を計算します．ここで，FVは化合物からベクトルデータへのマッピングです．詳細は[Module 1](Module_1/)を参照してください．
+ モジュール２では既知の化合物（FVで与えられる）とその性質から学習する*Artificial Neural Network*（ANN）を構築します．このANNは，与えられた化合物の物性を推論するために使用できます．詳細は[Module 2](Module_2/)を参照してください．
+ モジュール３では，モジュール２で学習されたANNと目標物性から**グラフ記述子のベクトル**（すなわちFV）を推論する*Mixed-Integer Linear Programming*（MILP）を実装します．詳細は[Module 3](Module_3/)を参照してください．
+ モジュール４では，与えられたFVに基づいてFVを満たす非循環化学グラフを生成する2-branch構造特性に基づいた*部分列挙*を提供します（[論文](https://arxiv.org/abs/2009.09646)を参照してください）．
最初のプログラムは2-branch-number *2*のグラフを出力し，2番目のプログラムでは2-branch-number *4*のグラフを出力します．詳細は[Module 4](Module_4/)を参照してください．

モジュールがこれらのタスクをどう処理するかを理解したい場合は，[論文](https://arxiv.org/abs/2009.09646)を参照してください．

## 要件

モジュール1とモジュール4はC++で書かれています．
ISO C++ 2011 standardと互換性のあるコンパイラであれば動作するはずです．
Linux Mint 18, 19, 20でg++をテストしており，インストールされていない場合は次のコマンドでインストールできます．
```shell
$ sudo apt install g++
```

モジュール2とモジュール3はPythonで書かれています．
Python 3とライブラリとしてscikit-learn, PuLP, pandas, numpyが必要です．
Linux Mint 18/19/20では次のコマンドでインストールできます．
```
$ sudo apt install python3-sklearn python3-pulp python3-pandas python3-numpy
```

### その他OSについて

他のOSで他のプログラムを実行する方法についての情報を歓迎します．

## コンパイルと使い方

それぞれのモジュールのユーザーマニュアルを参照してください．

## 謝辞

このプロジェクトはJSPS Grant（KAKENHI）18H04113に一部サポートしていただいてます．
This project is partially supported by JSPS Grant (KAKENHI) 18H04113.
