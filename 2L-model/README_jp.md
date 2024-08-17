---
title: "Mol-Inferプロジェクトにおける2L-LLRパッケージの使い方"
date: "2021年3月20日"
author: "京都大学　情報学研究科　離散数理分野研究室"
---

<p align="center">
  <a href="/2L-model/README_en.md">English</a>
  ·
  <a href="/2L-model/README_jp.md">日本語</a>
</p>

**お願い**

Todo:
+ Module1とModule2の英語マニュアルを追加
+ Module3の日本語マニュアルを追加

もしバグや間違いを見つけた場合は次回のリリースで修正するのでお知らせください．


---

## Mol-Inferプロジェクトについて

Mol-inferとは，京都大学の離散数理分野研究室（ku-dml）が開発したプロジェクトです．詳細は[トップページ](https://github.com/ku-dml/mol-infer)を参照してください.

## 2L-modelパッケージの紹介

このパッケージは4つのモジュールで構成されています．

+ モジュール1では記述子を導出します．詳細は[Module 1](/Module_1/)を参照してください．
+ モジュール2では*Artificial Neural Network* (ANN)を用いて予測関数を構築し,与えられた化合物の物性値の予測に用いられます．詳細は[Module 2](/Module_2/)を参照してください．
+ モジュール3では逆問題を解く*Mixed-Integer Linear Programming*（MILP）を実装します．詳細は[Module 3](/Module_3/)を参照してください．
+ モジュール4ではグラフを生成します（部分列挙）．詳細は[Module 4](/Module_4/)を参照してください．

これらのモジュールがどうやってこれらのタスクを処理するかをより具体的に理解したい場合は，[論文](https://doi.org/10.3390/ijms22062847)を参照してください．

## コンパイル及び使い方について

それぞれのモジュールのユーザーマニュアルを参照してください．

## 入出力ファイルの使い方
![入出力ファイルの使い方](illustration.jpg)

