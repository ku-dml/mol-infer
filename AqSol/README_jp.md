<p align="center">
  <a href="/AqSol/README_en.md">English</a>
  ·
  <a href="/AqSol/README_jp.md">日本語</a>
</p>

このフォルダには、この論文で使用されたコードとデータセットが含まれています。

Muniba Batool, Naveed Ahmed Azam, Jianshen Zhu, Kazuya Haraguchi, Liang Zhao and Tatsuya Akutsu:
A Unified Model for Inferring Chemical Compounds with Given Aqueous Solubility.

このフォルダの構造は以下の通りです．
	
1. モジュール１
	a. 前処理
		 PubChemデータベースからSDFファイルを抽出．
		 化合物の除去（セクション4で説明）
	b. 記述子の生成 (セクション3.1.1で説明)
2. モジュール２
	MLR，LLR-ANN，LLR-LLR，FSP-MLR，MLR-LOO，LLR-ANN-LOO，FSP-MLR-LOO，FSP-LOO-MLRで予測関数の構築（セクション3.1.5で説明）
3. モジュール３
	推論のためのMILPの計算（セクション3.2で説明）

それぞれのフォルダの詳細な説明については、各フォルダを参照してください。