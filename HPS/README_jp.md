<p align="center">
  <a href="/HPS/README_en.md">English</a>
  ·
  <a href="/HPS/README_jp.md">日本語</a>
</p>

# mol-infer/HPS
このフォルダには，この論文で使用されたコードとデータセットが含まれています．

<p>Jianshen Zhu, Naveed Ahmed Azam, Kazuya Haraguchi, Liang Zhao, Hiroshi Nagamochi and Tatsuya Akutsu:<br>
  Molecular Design Based on Integer Programming and Splitting Data Sets by Hyperplanes.<br>
  <i>IEEE/ACM Transactions on Computational Biology and Bioinformatics</i>, 2024. DOI: 10.1109/TCBB.2024.3402675.
</p>

- https://ieeexplore.ieee.org/document/10534864/
- プレプリントは https://arxiv.org/abs/2305.00801 で入手できます．
- 全てのフォルダ/ファイル名はTCBB論文の表/図番号に基づいています．
- このディレクトリ内のファイルを使用する場合は，TCBB論文を引用してください．

-----
フォルダー構造は以下のようになっています．
1. モジュール１
   - 線形記述子の生成（セクション3で説明）
   - 二次記述子の生成（セクション3で説明）
   - 前処理
1. モジュール２
   - ハイパープレーンを使ったデータセットの分割（セクション5）
   - 予測関数の構築（セクション4）
1. モジュール３
   - 逆問題のためのMILPの計算（セクション6.2）
   - 近傍解の生成（セクション6.2）
1. モジュール４
   - 再組み合わせ解の生成（セクション6.2）
1. 論文用インスタンス
   - 論文で使用された元のデータセット
   - 表2, 3, 4で使用されたインスタンスファイル
  
詳細については、各フォルダを参照してください。

