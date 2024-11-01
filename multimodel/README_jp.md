<p align="center">
  <a href="/Multimodel/README.md">English</a>
  ·
  <a href="/Multimodel/README_jp.md">日本語</a>
</p>

# mol-infer/Multimodel

このフォルダーには, 複数の所望の物性値を持つ化合物を推論するコードとインスタンスが含まれています．

フォルダー構造は以下のようになっています．
1. モジュール 3
    - 逆問題のためのMILPの計算
    - 複数の物性に対する目標値を持つ化合物の推論

# Module 3
モジュール 3 は, 化合物の化学グラフを推論するための逆問題のためのMILP式を解くステージから構成されます．

入力:
- instance_file
- fringe_tree_file
- モジュール 2 から得られた予測関数
- 複数の物性の目標値範囲

出力:
- 予測された物性値が目標値の範囲内にある化学グラフ

計算された記述子を確認するには
```bash
cd libs/2LMM_v019/
make FV_2LMM_V019
```

コードを実行するサンプルコマンド
```bash
python3 -m Module_3 config/config.yaml
```

# config ファイル
config ファイルは以下のパラメータを含む yaml ファイルです:
- instance_file: インスタンスファイルへのパス
- fringe_tree_file: 外縁木ファイルへのパス
- output_prefix: 出力ファイルのプレフィックス
- input_data: 以下のパラメータを含む辞書のリスト
  - model: 予測に使用されるモデル (LR, RF, ANN)
  - prefix: 入力データファイルのプレフィックス
  - target_value_lower_bound: 目標値範囲の下限
  - target_value_upper_bound: 目標値範囲の上限
