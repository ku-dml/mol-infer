---
title: "Mol-Inferプロジェクトにおける2LMM-LLRパッケージの使い方"
update date: "2021年7月4日"
author: "京都大学　情報学研究科　離散数理分野研究室"
---

<p align="center">
  <a href="/2LMM-LLR/README.md">English</a>
  ·
  <a href="/2LMM-LLR/README_jp.md">日本語</a>
</p>

**お願い**

もしバグや間違いを見つけた場合は次回のリリースで修正するのでお知らせください．

---

## Molecular Inferingについて
Mol-inferプロジェクトとは，京都大学の離散数理分野研究室（ku-dml）が開発したプロジェクトです．詳細は[トップページ](https://github.com/ku-dml/mol-infer)を参照してください.

## 2LMM-LLRパッケージの紹介

このパッケージは4つのモジュールで構成されています．

+ モジュール1では記述子を導出します．詳細は[Module 1](src/Module_1/)を参照してください．
+ モジュール2では*Lasso Linear Regression*（LLR）を用いて予測関数を構築します．詳細は[Module 2](src/Module_2/)を参照してください．
+ モジュール3では逆問題を解く*Mixed-Integer Linear Programming*（MILP）を実装します．詳細は[Module 3](src/Module_3/)を参照してください．
+ モジュール4ではグラフを生成します（部分列挙）．詳細は[Module 4](src/Module_4/)を参照してください．

このパッケージが具体的にどうやってこれらのタスクを処理するかを理解したい場合は，[論文](https://arxiv.org/abs/2107.02381)を参照してください．

## クイックスタートビジュアルガイド

使い方を示した[動画](https://www.youtube.com/watch?v=NfWlQbfs3qg)を公開しています．

2LMM-LLRパッケージ内でこれらのモジュールを通じてデータ及びファイルがどのように処理されているかを記してある画像を以下に表示します．
より細かい使い方やコンパイル方法については，各モジュールのユーザーマニュアルを参照してください．

![Data flow illustration](/2LMM-LLR/doc/2LMM-LLR_flow.PNG)
