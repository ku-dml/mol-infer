import sys
pwd = sys.path[0]
Module3_path = pwd[:pwd.rfind("MILP-range")+len("MILP-range")]
sys.path.append(Module3_path)

from Module_3.main_infer import infer

print(infer("config/config.yaml"))