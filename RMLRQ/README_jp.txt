このモデルには４つのモデルがあり，モジュール１とモジュール４は先行研究と同じです．
詳細は"ku-dml/mol-infer/2LMM-LLR/"のモジュール１とモジュール４を確認してください．

このフォルダーはモジュール２とモジュール３の詳細を含んでいます．

モジュール2では2次記述子を生成しR-MLRを使用して予測関数を得ます．
このフォルダは2つのサブフォルダから構成されます．
	generate_quadratic_descritporsフォルダには2次記述子を生成するコードが含まれます．
	construct_prediction_functionsフォルダはMLRに基づく予測関数を構築するために使用される記述子を選択するBSP手順を使用するコードが含まれます．

モジュール3ではMILPを解いて化学グラフを生成します．
このフォルダには2つのサブフォルダがあります．
	MILPフォルダにはケミカルグラフを生成するコードが含まれています．
	詳細は"ku-dml/mol-infer/2LMM-LLR"のモジュール3を参照してください．
	GNSフォルダにはより多くのケミカルグラフを生成するためにグリッド近傍探索を使用するコードが含まれています．
	詳細はku-dml/mol-infer/Grid-neighbor-searchを参照してください．
	
各フォルダにはそれぞれreadme，コード，サンプルのインスタンスがあります．

参考文献: "Molecular Design Based on Integer Programming and Quadratic Descriptors in a Two-layered Model" by
Jianshen Zhu, Naveed Ahmed Azam, Shengjuan Cao, Ryota Ido, Kazuya Haraguchi, Liang Zhao, Hiroshi Nagamochi and Tatsuya Akutsu

