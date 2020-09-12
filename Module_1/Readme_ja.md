[English](Readme.md)

# Module 1: SDFファイルから特徴ベクトルを計算する

最終更新日：2020年9月12日

## クィックスタート

モジュール１は、与えられた化合物の特徴ベクトル (FV) を計算する。
一個または複数個の化合物は、標準的なSDFファイルによって与えられる。
特徴ベクトルは、独自のFVフォーマットで出力される。
これらのフォーマットは、以下の節で説明する。
ここではまず使用例を示す。

コンパイル
```
g++ -std=c++11 -Wall -O3 -o fv4_in_ex fv4_in_ex.cpp
```
実行  
```
$ ./fv4_in_ex input.sdf output.csv
```
一番目と二番目の引数がそれぞれ入力のSDFファイルと出力のFVファイルを指定する。例えば、同梱の例に対しては次のようになる。
```
$ ./fv4_in_ex sample1.sdf sample1.csv
```

## SDF フォーマット

このプログラムは入力のフォーマットとして標準的な SDF (Structure Data File) を使用している。
https://www.chem-station.com/blog/2012/04/sdf.html などの解説が分かりやすい。
さらに、正確な定義書としては、公式資料 http://help.accelrysonline.com/ulm/onelab/1.0/content/ulm_pdfs/direct/reference/ctfileformats2016.pdf を参照するとよい。例として、sample1.sdf (https://pubchem.ncbi.nlm.nih.gov/compound/128703) を添付した。

## FV フォーマット

本プログラムの出力ファイルは、独自の FV (Feature Vector, 特徴ベクトル) フォーマットを採用している。
このテキストファイルは、カンマで区切った CSV ファイルであり、Excelなどの表計算ソフトで開くことができる。
具体的に、一行目には特徴ベクトルの構成要素が、二行目以降の各行には特徴ベクトルの数値データが記入される。
例として、sample1.sdf に対して計算した結果 sample1.csv を示す。それぞれの構成要素はそのあとに説明する。

```
CID,n,M,C_in,C_ex,O_in,O_ex,N_in,N_ex,H,C1O_in,C1O_ex,C2O_in,C2O_ex,C1N_in,C1N_ex,C1C_in,C1C_ex,#degree1_in,#degree1_ex,#degree2_in,#degree2_ex,#degree3_in,#degree3_ex,#degree4_in,#degree4_ex,#double_bond_in,#double_bond_ex,#triple_bond_in,#triple_bond_ex,Diameter,Bc_121_in,Bc_121_ex,Bc_122_in,Bc_122_ex,Bc_123_in,Bc_123_ex,Bc_131_in,Bc_131_ex,Bc_132_in,Bc_132_ex,Bc_141_in,Bc_141_ex,Bc_221_in,Bc_221_ex,Bc_222_in,Bc_222_ex,Bc_223_in,Bc_223_ex,Bc_231_in,Bc_231_ex,Bc_232_in,Bc_232_ex,Bc_241_in,Bc_241_ex,Bc_331_in,Bc_331_ex,Bc_332_in,Bc_332_ex,Bc_341_in,Bc_341_ex,Bc_441_in,Bc_441_ex,2-branch_height,2-branch_leaf_number
128703,24,130.833,14,3,0,6,1,0,0,0,5,0,1,2,0,12,3,0,7,10,2,5,0,0,0,0,1,0,0,0.75,0,2,0,0,0,0,0,4,0,1,0,0,9,1,0,0,0,0,1,1,0,0,0,0,4,0,0,0,0,0,0,0,1,2
```

構成要素の説明
+ CID:  PubChem (https://pubchem.ncbi.nlm.nih.gov) におけるCID。例えば sample1.sdf にある要素は、https://pubchem.ncbi.nlm.nih.gov/compound/128703 である。
+ n: 水素を除く原子の数
+ M: 平均分子質量
![M = \frac{1}{n}\sum_{a}\lfloor 10 \cdot mass(a)\rfloor](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+M+%3D+%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Ba%7D%5Clfloor+10+%5Ccdot+mass%28a%29%5Crfloor)
+ C_in, O_in, N_in: 内部原子の数
+ C_ex, O_ex, N_ex: 外部原子の数
+ H: 原子の数
+ C1O_in, C2O_in, C1N_in, C1C_in: 内部パスの数。例えば、C1O_in は、CとOの内部単結合の数、C2O_in はCとOの内部二重結合の数を示す。
+ C1O_ex, C2O_ex, C1N_ex, C1C_ex: 外部パスの数。例えば、C1O_in は、CとOの外部単結合の数、C2O_in はCとOの外部二重結合の数を示す。
+ #degree1_in,#degree2_in,#degree3_in,#degree4_in: それぞれの次数 (価数) を持つ内部原子の数
+ #degree1_ex,#degree2_ex,#degree3_ex,#degree4_ex: それぞれの次数 (価数) を持つ外部原子の数
+ #double_bond_in,#triple_bond_in: 内部二重結合と内部三重結合の数
+ #double_bond_ex,#triple_bond_ex: 外部二重結合と 外部三重結合の数
+ Diameter: 直径/n
+ Bc_xyz_in: 内部次数構成 (x, y, z)、ただし x ≤ y は z-重結合の両端点の次数を表す。
+ Bc_xyz_ex: 外部次数構成 (x, y, z)、ただし x ≤ y は z-重結合の両端点の次数を表す。
+ 2-branch_height: 2-branch-height bh_2。別途論文を参照されたい。
+ 2-branch_leaf_number: 2-branch-leaf-number bl_2。別途論文を参照されたい。

## プログラムノート
計算プログラム fv4_in_ex は、入力の SDF ファイルに対して、計算された特徴ベクトルを FV フォーマットで出力する。
ここではプログラム利用上の注意事項を示す。
1. 原子の質量は、プログラムの中にハードコード仕様になっている。執筆の時点では、次のようになっているが、
必要の場合、追加してご利用下さい。
```
対応している原子の質量 (function init_MassMap())
M["B"]  = 108;
M["C"]  = 120;
M["O"]  = 160;
M["N"]  = 140;
M["F"]  = 190;
M["Si"] = 280;
M["P"]  = 310;
M["S"]  = 320;
M["Cl"] = 355;
M["V"]  = 510;
M["Br"] = 800;
M["Cd"] = 1124;
M["I"]  = 1270;
M["Hg"] = 2006;
M["Pb"] = 2072;
M["Al"] = 269;
```
2. デバッグ出力をみるには、 fv4_in_ex.cppの```bool debug = true;``` (line 27) にしてコンパイルしてください。
