
import os

infer_name = "infer_2LMM_SEP.py"

prop_n = "At"
prop = f"./MILP_result/{prop_n}/{prop_n}"

tv_list = [("240", "250"), ("240", "250"), ("240", "250"), ("190", "200"),("300", "310"), ("360", "370"),("230", "240")]
abbr_list = ["a", "b1", "b2", "b3", "b4", "c", "d"]

for (instance_name_abbr, (tv_lb, tv_ub)) in zip(abbr_list, tv_list):
    instance_name = f"./MILP_result/{prop_n}/instance_{instance_name_abbr}_2LMM.txt"
    fringe_name = f"./MILP_result/{prop_n}/ins_{instance_name_abbr}_fringe_2LMM.txt"

    os.system(f"python {infer_name} {prop} {tv_lb} {tv_ub} {instance_name} {fringe_name} ./MILP_result/sdf/{prop_n}_{instance_name_abbr}_{tv_lb}_{tv_ub} -d1 ANN -d2 R-MLR > ./MILP_result/stdout/stdout_{prop_n}_{instance_name_abbr}_{tv_lb}_{tv_ub}.txt")


