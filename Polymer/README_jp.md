<p align="center">
  <a href="/Polymer/README_en.md">English</a>
  ·
  <a href="/Polymer/README_jp.md">日本語</a>
</p>

# mol-infer/Polymer

このフォルダーは，"A Method for Inferring Polymers Based on Linear Regression and Integer Programming"の論文で使用されたコードとインスタンスを含んでいます．

- この論文はIEEE/ACM Transactions on Computational Biology and Bioinformaticsに提出しました．
   - 付録は[Supplementary_Materials.pdf](Supplementary_Materials.pdf)を参照してください．
- プレプリントは https://arxiv.org/abs/2109.02628 で入手できます．
- このリポジトリ内の全てのフォルダ/ファイル名は提出された論文の表/図番号に基づいています．

フォルダー構造は以下のようになっています．
1. モジュール１
   - 記述子の生成（セクション3で説明）
   - 前処理
1. モジュール２
   - 予測関数の構築（セクション4, ステージ3）
1. モジュール３
   - 逆問題のためのMILPの計算（セクション4, ステージ4, 表2）
   - 複数の物性に対する目標値を持つポリマーの推論（セクション4, ステージ4, 表3）
1. モジュール４
   - 再組み合わせ解の生成（セクション4, ステージ5）
1. 論文用インスタンス
   - 論文で使用された元のデータセット
   - 表2と表3で使用されたインスタンスファイル

詳細については、各フォルダを参照してください。
