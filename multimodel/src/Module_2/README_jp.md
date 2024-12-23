<p align="center">
  <a href="/multimodel/src/Module_2/README.md">English</a>
  ·
  <a href="/multimodel/src/Module_2/README_jp.md">日本語</a>
</p>

# Module 2: 

モジュール2は、Lasso Linear Regression (LLR)、Random Forest (RF)、またはArtificial Neural Network (ANN)を使用して予測関数を構築します。

# 事前準備:

CIDと物性値のファイルを準備してください。物性値のファイルのサンプルは /sample_instance/input/ にあります。

# 予備実験:

使用法:

	python SINGLE_pre.py (.csv) (_values.txt) -l(lasso)/-ann(ANN)/-rf(RF)
サンプル:

	python SINGLE_pre.py ../Module_1/sample_instance/output/Bp_small_desc_norm.csv sample_instance/input/Bp_small_norm_values.txt -l


出力は次のようになります:

```
Bp_small	400	15	15	0	LASSO	0.549641087647406	0.41578203853379403	33.90013933181763 -l ./log/LLR_LASSO_Bp_small_desc_norm.csv 0.54 258 30 20
```

最後の項目（-annから始まるか、Lassoを使用する場合は-lから始まる）をコピーし、次のコマンドを使用して評価実験を実行してください。

# 評価実験:
使用法:
  
	python SINGLE_eval.py (.csv) (_values.txt) (-o (OUTPUT_PREFIX)) (-..., 予備実験の結果の最後の項目)
サンプル:

	Bp_small	400	15	15	0	LASSO	0.549641087647406	0.41578203853379403	33.90013933181763		-l ./log/LLR_LASSO_Bp_small_desc_norm.csv 0.54 258 30 20
