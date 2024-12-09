# 事前準備 
## fv_2LMM.cpp のコンパイル
```bash
cd src/Module_1
g++ -O2 -Wall -std=c++20 -o FV_2LMM_V019 fv_2LMM.cpp
```

## データの準備
eliminate.py で不要なデータを削除
```bash
python eliminate.py (data.sdf)
```

sample
```bash
python eliminate.py sample_instance/Bp_small.sdf
```

(optional) 
limit_atoms.py で原子種を制限
```bash
python limit_atoms.py (data.sdf) (原子の種類)
```

sample
```bash
python limit_atoms.py sample_instance/Bp_small.sdf C O N S Cl H
```

# 実行
## 特徴ベクトルの計算
サンプルスクリプト
```bash
./FV_2LMM_V019 sample_instance/Bp_small.sdf sample_instance/output/Bp_small 
```