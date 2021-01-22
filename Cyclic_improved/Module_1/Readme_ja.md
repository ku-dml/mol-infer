[English](Readme.md)

# Module 1: SDFファイルから特徴ベクトルを計算する

最終更新日：2021年1月20日

モジュール１は、与えられた化合物の特徴ベクトルを計算する。


## クイックスタート

ここではまず使用例を示す。

化合物は標準的なSDFファイルによって与えられなければならない。

特徴ベクトルは、独自のFVフォーマットを持つファイル (csv形式) で出力される。
これらのフォーマットは、以下の節で説明する。

### SDFファイルからFVファイル (csv形式) を生成する

コンパイル
```
g++ -std=c++11 -Wall -O3 -o FV_ec fv_ec.cpp
```
fv_ec.cpp と common.cpp は同じフォルダに置かなければならない。

実行  
```
$ ./FV_ec INPUT.sdf OUTPUT.csv
```
一番目と二番目の引数には、
それぞれ入力のSDFファイルと出力のFVファイルを指定する。
例えば、同梱の例に対しては次のようになる。
```
$ ./FV_ec sample1.sdf sample1.csv

$ ./FV_ec sample2.sdf sample2.csv
```

引数を与えずに実行すれば (もしくは引数が適切に与えられなかった場合)、
指定すべき引数が出力されて終了する。

### SDFファイルから、すでに得られたFVファイル (csv形式) と同じ記述子を持つ特徴ベクトルを生成する

コンパイル
```
g++ -std=c++11 -Wall -O3 -o FV_proj fv_proj.cpp
```
fv_proj.cpp と common.cpp は同じフォルダに置かなければならない。


実行  
```
$ ./FV_proj DESCRIPTOR.csv INPUT.sdf OUTPUT.csv
```
一番目、二番目、三番目の引数には、それぞれ
- 入力のFVファイル (csv)
- 入力のSDFファイル
- 出力のFVファイル (csv)
を指定する。例えば、同梱の例に対しては次のようになる。
```
$ ./FV_ec sample2.csv sample1.sdf sample1_on_2.csv
```
sample1.sdfに含まれる化合物に関して、
sample2.csvに現れる記述子に限定された特徴ベクトルが
sample1_on_2.csvに出力される。

引数を与えずに実行すれば (もしくは引数が適切に与えられなかった場合)、
指定すべき引数が出力されて終了する。

### 2つのプログラムの役割
2つのプログラムは、sdfファイルに含まれる化合物を
Module 2 (ニューラルネットワークの学習) でどのように取り扱うかによって使い分ける必要がある。
+ FV_ec: INPUT.sdf ファイルに含まれる化合物から学習を行い、新しくニューラルネットワークを作りたい場合。
+ FV_proj: すでに DESCRIPTOR.csv からニューラルネットワークが作られていて、
そのニューラルネットワークを使って、INPUT.sdf に含まれる化合物の化学的性質の値を
推定したい場合。

### 重要な注意

この Cyclic_improved モデル、および前身の Cyclic モデルでは、入力される化合物にいくつかの仮定をしている。

+ <font color="red">すべての化合物は環を含んでいなければならない。</font>（すなわちグラフ理論で言うところの閉路である。）環を含まない化合物がSDFに存在する場合、特徴ベクトル生成器は無限ループに陥る可能性がある。
+ ***芳香辺 (aromatic edge) に対応していない。***必ずボンド幅が 1 もしくは 2 の辺に置き換えたSDFを用いること。
+ ***炭素数が3以下、電荷を帯びた原子を持つ化合物には対応していない。***また原子の質量および価数も、開発側が想定している値を取らなければならない。その値については fv_ec.cpp のコードを参照すること。

特徴ベクトル生成器 FV_ec には、SDFが上記仮定に沿ったものか否かを判定する機能はついていない。使用者は自分でSDFを整形する必要がある。


## SDF フォーマット		

このプログラムは入力のフォーマットとして標準的な SDF (Structure Data File) を使用している。
https://www.chem-station.com/blog/2012/04/sdf.html などの解説が分かりやすい。
さらに、正確な定義書としては、公式資料 http://help.accelrysonline.com/ulm/onelab/1.0/content/ulm_pdfs/direct/reference/ctfileformats2016.pdf を参照するとよい。例として、sample1.sdf (https://pubchem.ncbi.nlm.nih.gov/compound/6140) を添付した。

## FV フォーマット

本プログラムの出力ファイルは、独自の FV (Feature Vector, 特徴ベクトル) フォーマットを採用している。
このテキストファイルは、カンマで区切った CSV ファイルであり、Excelなどの表計算ソフトで開くことができる。
具体的に、一行目には特徴ベクトルの記述子が、二行目以降の各行には特徴ベクトルの数値データが記入される。
例として、sample1.sdf に対して FV_ec を実行した結果得られる sample1.csv を示す。それぞれの構成要素はそのあとに説明する。

```
CID,n,cs,ch,bl_2,ms,dg_co_1,dg_co_2,dg_co_3,dg_co_4,dg_nc_1,dg_nc_2,dg_nc_3,dg_nc_4,bd_co_2,bd_co_3,bd_in_2,bd_in_3,bd_ex_2,bd_ex_3,ns_co_C3,ns_co_C2,ns_nc_O1,ns_nc_N1,ns_nc_C2,ns_nc_C3,ec_co_C2_C3_2,ec_co_C2_C2_1,ec_co_C2_C3_1,ec_co_C2_C2_2,ec_in_C2_C3_1,ec_in_C3_C2_1,ec_ex_C3_N1_1,ec_ex_C3_C3_1,ec_ex_C3_O1_1,ec_ex_C3_O1_2,nsH
6140,12,6,4,1,128.333,0,5,1,0,3,1,2,0,3,0,0,0,1,0,1,5,2,1,1,2,1,2,1,2,1,1,1,1,1,1,11
```

構成要素の概要は以下の通り。詳細は論文を参照のこと。

+ CID:  PubChem (https://pubchem.ncbi.nlm.nih.gov) におけるCID。例えば sample1.sdf にある要素は、https://pubchem.ncbi.nlm.nih.gov/compound/6140 (フェニルアラニン; Phenylalanine) である。
+ n: 水素を除く原子の数
+ cs: コアに属する原子の数
+ ch: コアの高さ
+ bl: 2-leafの個数
+ ms: 平均分子質量  
![M = \frac{1}{n}\sum_{a}\lfloor 10 \cdot mass(a)\rfloor](https://render.githubusercontent.com/render/math?math=%5Cdisplaystyle+M+%3D+%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Ba%7D%5Clfloor+10+%5Ccdot+mass%28a%29%5Crfloor)
+ dg_co_1, dg_co_2, dg_co_3, dg_co_4: コアに属する原子で、次数が 1,2,3,4 のものの個数
+ dg_nc_1, dg_nc_2, dg_nc_3, dg_nc_4: コアに属さない原子で、次数が 1,2,3,4 のものの個数 
+ bd_co_2, bd_co_3: コア二重結合およびコア三重結合の数
+ bd_in_2, bd_in_3: 内部二重結合および内部三重結合の数
+ bd_ex_2, bd_ex_3: 外部二重結合および外部三重結合の数
+ ns_co_Xd: コアに属する元素記号 X の原子で、次数が d のものの個数。たとえば ns_co_C3 は、コアに属する炭素原子 C で次数が3のものの個数を示す。
+ ns_nc_Xd: コアに属さない元素記号 X の原子で、次数がdのものの個数。
+ ec_co_Xx_Yy_2, ec_co_Xx_Yy_3: コア二重結合およびコア三重結合であって (無向辺)、次数 x の原子 X と次数 y の原子 Y を結ぶものの本数。たとえば ec_co_C2_C3_2 は、次数2の炭素原子 C および次数3の炭素原子 C を結ぶコア二重結合の本数を示す。 
+ ec_in_Xx_Yy_2, ec_in_Xx_Yy_3: 内部二重結合および内部三重結合であって (親から子への有向辺)、次数 x の原子 X と次数 y の原子 Y を結ぶものの本数。
+ ec_ex_Xx_Yy_2, ec_ex_Xx_Yy_3: 外部二重結合および外部三重結合であって (親から子への有向辺)、次数 x の原子 X と次数 y の原子 Y を結ぶものの本数。 
+ nsH: 水素原子の個数。

※ ns,ecで始まる記述子は、入力SDF内の化合物に現れるもののみがFVファイルに出力される。

## プログラムノート

ここではプログラム利用上の注意事項を示す。

1. FV_ec の本体の C++ プログラムは fv_ec.cpp,
FV_proj の本体の C++ プログラムは fv_proj.cppとなっているが、
共通して用いられる関数は common.cpp に入っている。
両プログラムとも、コンパイル時にはソースファイルを同じディレクトリに置いておく必要がある。

2. 原子の質量は、プログラムの中にハードコード仕様になっている。common.cpp内の関数 init_MassMap() で以下のように定められているが、質量の変更や、ほかの原子が必要な場合には、編集して再度コンパイルすることで利用できる。
```
対応している原子の質量 (common.cpp にある関数 init_MassMap())
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

3. デバッグ出力をみるには、
fv_ec.cppの29行目もしくは
fv_proj.cpp の37行目を
```bool debug = true;```
としてコンパイルし直すこと。
