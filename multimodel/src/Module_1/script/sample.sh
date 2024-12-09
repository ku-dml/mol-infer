# 物性
props=(Bp Fp Kow Lp Mp Sl)

# eliminate.py
for prop in ${props[@]}; do
  python3 eliminate.py sample_instance/input/${prop}_small.sdf
done

# FV_2LMM_V019
for prop in ${props[@]}; do
  ./FV_2LMM_V019 sample_instance/input/${prop}_small.sdf sample_instance/output/${prop}_small
done
