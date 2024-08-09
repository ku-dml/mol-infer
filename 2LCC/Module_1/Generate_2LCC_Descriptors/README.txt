This module corresponds to the second stage in the paper.

/***** 2LCCM_v030 *****/

- 2LCCM_v030 is different from 2LMM_v018 in that
  v030 has the CC ksi(C) of cycle C as a descriptor.

- 2LMM_v018 is different from 2LMM_v017 in that
  v018 has the rank r(G) of G as a descriptor. 

- 2LMM_v017 is different from 2LMM_v016 in the order of descriptors; 
  in this version, leaf_ac comes after fringe-configuration. 

- Haraguchi started 2LMM_v016
  by the following emails by the professor to him on Jun 6th, 2021. 
  (It is now called leaf_ac; written on July 1st, 2021)

-----
以前の特徴ベクトルとの比較実験を見ると，
BpやGapなどは明らかに2LM-Mのほうが成績が落ちているので，
葉の枝のacの情報は残すべきではないかと推察しています．
そこで，葉uとその隣接点vに対して，順序付き ac として
　 (alpha(u), alpha(v), beta(uv))
を定めます（これは (alpha(v), alpha(u), beta(uv))と区別する）．

そして，2LM-Mで葉（次数１）の特徴ベクトルに
葉の枝の 順序付き ac の頻度を記述子として加えて，
BpやGapに対し再度学習実験をしてもらえますか．  

この実験結果によって，今後，2LM-Mに 葉の枝の 順序付き ac を
追加するか決めたいと考えています．
-----
はい，水素無しの化学グラフでの話です．
cation/anionについてもとりあえず無視してください．
-----
