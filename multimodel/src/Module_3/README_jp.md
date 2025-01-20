<p align="center">
  <a href="/multimodel/src/Module_3/README.md">English</a>
  ·
  <a href="/multimodel/src/Module_3/README_jp.md">日本語</a>
</p>

# Module 3
モジュール 3 は, 化合物の化学グラフを推論するための逆問題のためのMILP式を解くステージから構成されます．

入力:
- instance_file
- fringe_tree_file
- モジュール 2 から得られた予測関数
- 物性の目標値範囲

出力:
- 予測された物性値が目標値の範囲内にある化学グラフ

計算された記述子を確認するには
```bash
make -C Module_3/libs/2LMM_v019/ FV_2LMM_V019
```

サンプルコマンド
```bash
python -m Module_3 --config-name=config_sample
```

# config ファイル
config ファイルの設定は config/config_sample.yaml を参照してください．
