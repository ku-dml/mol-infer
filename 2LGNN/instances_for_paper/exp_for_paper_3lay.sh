#!/bin/sh


python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 3_homo_3lay_16_32_FT16 -val homo > qm9_3_homo.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 4_lumo_3lay_16_32_FT16 -val lumo > qm9_4_lumo.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 5_gap_3lay_16_32_FT16 -val gap > qm9_5_gap.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 1_mu_3lay_16_32_FT16 -val mu > qm9_1_mu.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 2_alpha_3lay_16_32_FT16 -val alpha > qm9_2_alpha.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 6_r2_3lay_16_32_FT16 -val r2 > qm9_6_r2.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 7_zpve_3lay_16_32_FT16 -val zpve > qm9_7_zpve.txt
# python3 GNN_main_qm9.py qm9_all.sdf qm9_complete_values.csv -o 8_u0_3lay_16_32_FT16 -val u0 > qm9_8_u0.txt
# python3 GNN_main_qm9.py qm9_all.sdf qm9_complete_values.csv -o 9_u298_3lay_16_32_FT16 -val u298 > qm9_9_u298.txt
# python3 GNN_main_qm9.py qm9_all.sdf qm9_complete_values.csv -o 10_h298_3lay_16_32_FT16 -val h298 > qm9_10_h298.txt
# python3 GNN_main_qm9.py qm9_all.sdf qm9_complete_values.csv -o 11_g298_3lay_16_32_FT16 -val g298 > qm9_11_g298.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 12_cv_3lay_16_32_FT16 -val cv > qm9_12_cv.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 13_u0_atom_3lay_16_32_FT16 -val u0_atom > qm9_13_u0_atom.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 14_u298_atom_3lay_16_32_FT16 -val u298_atom > qm9_14_u298_atom.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 15_h298_atom_3lay_16_32_FT16 -val h298_atom > qm9_15_h298_atom.txt
python3 GNN_main_qm9.py ./qm9_all.sdf ./qm9_complete_values.csv -o 16_g298_atom_3lay_16_32_FT16 -val g298_atom > qm9_16_g298_atom.txt