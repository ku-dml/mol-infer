<p align="center">
  <a href="/AqSol/README_en.md">English</a>
  ·
  <a href="/AqSol/README_jp.md">日本語</a>
</p>

This folder contains codes and datasets used in this paper.

Muniba Batool, Naveed Ahmed Azam, Jianshen Zhu, Kazuya Haraguchi, Liang Zhao and Tatsuya Akutsu:
A Unified Model for Inferring Chemical Compounds with Given Aqueous Solubility.

The folder structure is organized as follows:
	
1. Module 1
	a. Pre-processing
		 Extracting SDF's from Pubchem database.
		 Eliminating compounds (describe in section 4)
	b. Generating descriptors (describe in section 3.1.1)
2. Module 2
	Constructing prediction function by MLR, LLR-ANN, LLR-LLR, FSP-MLR, MLR-LOO, LLR-ANN-LOO, FSP-MLR-LOO, FSP-LOO-MLR (describe in section 3.1.5)
3. Module 3
	Solving an MILP for inference problem(describe in section 3.2)	

Please refer each folder for detailed descriptions.