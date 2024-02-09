# reduce size of D_i before BSP with THRESHOLD

############### import ###############
import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model

from src import LR, ALR, BSP, ANN, RF, DT

import time, sys, copy, itertools, math, warnings, argparse, json

warnings.simplefilter('ignore')
####################################################
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv_D1", help='input csv file D1')
    parser.add_argument("input_value_D1", help='input value file D1')
    parser.add_argument("input_csv_D2", help='input csv file D2')
    parser.add_argument("input_value_D2", help='input value file D2')
    # parser.add_argument("output_pre", help='input output directory for results of preliminary experiments')
    # parser.add_argument("output_eval", help='input output file for results of evalutation experiments')
    parser.add_argument('-l1', '--lasso1', nargs='*', type=str)
    parser.add_argument('-rbsp1', '--rbsp1', nargs='*', type=str)
    parser.add_argument('-ann1', '--ann1', nargs='*', type=str)
    parser.add_argument('-rf1', '--rf1', nargs='*', type=str)
    parser.add_argument('-alr1', '--alr1', nargs='*', type=str)
    parser.add_argument('-dt1', '--dt1', nargs='*', type=str)

    parser.add_argument('-l2', '--lasso2', nargs='*', type=str)
    parser.add_argument('-rbsp2', '--rbsp2', nargs='*', type=str)
    parser.add_argument('-ann2', '--ann2', nargs='*', type=str)
    parser.add_argument('-rf2', '--rf2', nargs='*', type=str)
    parser.add_argument('-alr2', '--alr2', nargs='*', type=str)
    parser.add_argument('-dt2', '--dt2', nargs='*', type=str)

    args = parser.parse_args()

    return(args) 

def read_dataset(data_csv, value_txt):
    
    ### read files ###
    # read the csv and the observed values
    fv = pd.read_csv(data_csv) 
    value = pd.read_csv(value_txt)      

    fv_list = list(fv.columns)
    fv_list = fv_list[1:]

    ### prepare data set ###
    # prepare CIDs
    CIDs = np.array(fv['CID'])
    # prepare target, train, test arrays
    target = np.array(value['a'])
    # construct dictionary: CID to feature vector
    fv_dict = {}
    for cid,row in zip(CIDs, fv.values[:,1:]):
        fv_dict[cid] = row
    # construct dictionary: CID to target value
    target_dict = {}
    for cid, val in zip(np.array(value['CID']), np.array(value['a'])):
        target_dict[cid] = val
    # check CIDs: target_values_filename should contain all CIDs that appear in descriptors_filename
    for cid in CIDs:
        if cid not in target_dict:
            sys.stderr.write('error: {} misses the target value of CID {}\n'.format(target_values_filename, cid))
            exit(1)
    # construct x and y so that the CIDs are ordered in ascending order
    CIDs.sort()
    x = np.array([fv_dict[cid] for cid in CIDs])
    y = np.array([target_dict[cid] for cid in CIDs])

    return (CIDs,x,y,fv_list)

def read_dataset2(data_csv, value_txt, CIDs):
    
    ### read files ###
    # read the csv and the observed values
    fv = pd.read_csv(data_csv, index_col=0) 
    value = pd.read_csv(value_txt)      

    fv_list = list(fv.columns)
    # fv_list = fv_list[1:]

    # ### prepare data set ###
    # # prepare CIDs
    # CIDs = np.array(fv['CID'])
    # # prepare target, train, test arrays
    # target = np.array(value['a'])
    # # construct dictionary: CID to feature vector
    # fv_dict = {}
    # for cid,row in zip(CIDs, fv.values[:,1:]):
    #     fv_dict[cid] = row
    # # construct dictionary: CID to target value
    # target_dict = {}
    # for cid, val in zip(np.array(value['CID']), np.array(value['a'])):
    #     target_dict[cid] = val
    # # check CIDs: target_values_filename should contain all CIDs that appear in descriptors_filename
    # for cid in CIDs:
    #     if cid not in target_dict:
    #         sys.stderr.write('error: {} misses the target value of CID {}\n'.format(target_values_filename, cid))
    #         exit(1)
    # # construct x and y so that the CIDs are ordered in ascending order
    # CIDs.sort()

    x = fv.loc[CIDs]

    # x = np.array([fv_dict[cid] for cid in CIDs])
    # y = np.array([target_dict[cid] for cid in CIDs])

    return x.values, fv_list

def output_dataset(CIDs, fv_list, x, D_best, selected_fv_filename):

    x_new = pd.DataFrame(index=CIDs)
    x_new.index.name = "CID"
    for i in D_best:
        x_new[fv_list[i]] = x[:, i]

    x_new.to_csv(selected_fv_filename)

def main(args):

    try:
        CIDs_D1, x_D1, y_D1, fv_list_D1 = read_dataset(args.input_csv_D1, args.input_value_D1)
        CIDs_D2, x_D2, y_D2, fv_list_D2 = read_dataset(args.input_csv_D2, args.input_value_D2)

        way_D1 = ""
        if args.lasso1 is not None:
            way_D1 = "Lasso"
        elif args.rbsp1 is not None:
            way_D1 = "BSP"
        elif args.ann1 is not None:
            way_D1 = "ANN"
        elif args.rf1 is not None:
            way_D1 = "RF"
        # elif args.alr1 is not None:
        #     way_D1 = "ALR"
        elif args.dt1 is not None:
            way_D1 = "DT"
        else:
            sys.stderr.write("Please specify one learning method for D1\n")
            exit(1)

        way_D2 = ""
        if args.lasso2 is not None:
            way_D2 = "Lasso"
        elif args.rbsp2 is not None:
            way_D2 = "BSP"
        elif args.ann2 is not None:
            way_D2 = "ANN"
        elif args.rf2 is not None:
            way_D2 = "RF"
        # elif args.alr2 is not None:
        #     way_D2 = "ALR"
        elif args.dt2 is not None:
            way_D2 = "DT"
        else:
            sys.stderr.write("Please specify one learning method for D2\n")
            exit(1)

    except:
        sys.stderr.write("usage: {} (input_csv_D1.csv)(input_value_D1.txt)(input_csv_D2.csv)(input_value_D2.txt)\n\n".format(sys.argv[0]))
        exit(1)

    x_D1 = x_D1.astype(np.float64)
    y_D1 = y_D1.astype(np.float64)
    x_D2 = x_D2.astype(np.float64)
    y_D2 = y_D2.astype(np.float64)

    # re-normalization, only for pre
    min_y = np.amin(y_D1)
    max_y = np.amax(y_D1)
    y_D1 = (y_D1 - min_y) / (max_y - min_y)

    min_y = np.amin(y_D2)
    max_y = np.amax(y_D2)
    y_D2 = (y_D2 - min_y) / (max_y - min_y)

    sp_loc = args.input_csv_D1.find('_desc_norm')
    sp_loc2 = args.input_csv_D1.rfind('/')
    if sp_loc2 == -1:
        prop_name_D1 = args.input_csv_D1[:sp_loc]
    else:
        prop_name_D1 = args.input_csv_D1[sp_loc2 + 1:sp_loc]

    sp_loc = args.input_csv_D2.find('_desc_norm')
    sp_loc2 = args.input_csv_D2.rfind('/')
    if sp_loc2 == -1:
        prop_name_D2 = args.input_csv_D2[:sp_loc]
    else:
        prop_name_D2 = args.input_csv_D2[sp_loc2 + 1:sp_loc]

    nonzero_linear_org_D1 = 0
    nonzero_quadratic_org_D1 = 0
    for fv in fv_list_D1:
        if "*" in fv:
            nonzero_quadratic_org_D1 += 1
        else:
            nonzero_linear_org_D1 += 1

    nonzero_linear_org_D2 = 0
    nonzero_quadratic_org_D2 = 0
    for fv in fv_list_D2:
        if "*" in fv:
            nonzero_quadratic_org_D2 += 1
        else:
            nonzero_linear_org_D2 += 1

    time_start = time.time()

    log_filename_D1 = "./log/SEP_pre_" + prop_name_D1 + "_" + way_D1 + "_log.txt"
    with open(log_filename_D1, 'w') as f:
        f.close()

    log_filename_D2 = "./log/SEP_pre_" + prop_name_D2 + "_" + way_D2 + "_log.txt"
    with open(log_filename_D2, 'w') as f:
        f.close()

    if way_D1 == "Lasso":
        time_start = time.time()
        best_lambda_D1, best_R2_train_D1, best_R2_test_D1 = LR.Lasso_pre(x_D1, y_D1, log_filename_D1)
        time_D1 = time.time() - time_start

        with open(log_filename_D1, 'a') as f:
            f.write(f"BEST LAMBDA: {best_lambda_D1}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D1}\n")
            f.write(f"R2 TEST: {best_R2_test_D1}\n")
            f.write(f"TIME: {time_D1}\n")
            f.close()  

    elif way_D1 == "BSP":
        time_start = time.time()
        D_best_D1, best_R2_train_D1, best_R2_test_D1, best_linear_before_D1, best_quadratic_before_D1 = BSP.LLR_ALR_pre(x_D1, y_D1, fv_list_D1, log_filename_D1)
        time_D1 = time.time() - time_start

        selected_fv_filename_D1 = "./log/MLR_based_BSP_" + prop_name_D1 + "_desc_norm.csv"
        output_dataset(CIDs_D1, fv_list_D1, x_D1, D_best_D1, selected_fv_filename_D1)

        with open(log_filename_D1, 'a') as f:
            f.write(f"D_BEST: {D_best_D1}\n")
            f.write(f"D_BEST_SIZE: {len(D_best_D1)}\n")
            f.write(f"D_BEST_LINEAR_BEFORE: {best_linear_before_D1}\n")
            f.write(f"D_BEST_QUAD_BEFORE: {best_quadratic_before_D1}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D1}\n")
            f.write(f"R2 TEST: {best_R2_test_D1}\n")
            f.write(f"TIME: {time_D1}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D1}\n")
            f.close() 

    elif way_D1 == "ANN":
        time_start = time.time()
        if "quadratic" in args.input_csv_D1:
            quad = True
        else:
            quad = False
        best, best_R2_train_D1, best_R2_test_D1 = ANN.LLR_ANN_pre(x_D1, y_D1, fv_list_D1, log_filename_D1, quad)
        time_D1 = time.time() - time_start

        _, arch_D1, rho_D1, maxitr_D1 = best[0]
        arr_D1 = best[1][1]
        selected_fv_filename_D1 = "./log/LLR_ANN_" + prop_name_D1 + "_desc_norm.csv"
        output_dataset(CIDs_D1, fv_list_D1, x_D1, arr_D1, selected_fv_filename_D1)

        with open(log_filename_D1, 'a') as f:
            f.write(f"BEST ARCH: {arch_D1}\n")
            f.write(f"BEST RHO: {rho_D1}\n")
            f.write(f"BEST MAXITR: {maxitr_D1}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D1}\n")
            f.write(f"R2 TEST: {best_R2_test_D1}\n")
            f.write(f"TIME: {time_D1}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D1}\n")
            f.close()  

    elif way_D1 == "RF":
        time_start = time.time()
        if "quadratic" in args.input_csv_D1:
            quad = True
        else:
            quad = False
        best_R2_train_D1, best_R2_test_D1, best_params = RF.hyper_param_tuning(x_D1, y_D1, log_filename_D1, quad)
        time_D1 = time.time() - time_start

        k_D1, params_D1 = best_params
        arr_D1 = RF.kbest(x_D1, y_D1, k_D1)
        selected_fv_filename_D1 = "./log/RF_" + prop_name_D1 + "_desc_norm.csv"
        output_dataset(CIDs_D1, fv_list_D1, x_D1, arr_D1, selected_fv_filename_D1)        

        with open(log_filename_D1, 'a') as f:
            f.write(f"BEST K: {k_D1}\n")
            f.write(f"BEST PARAMS: {json.dumps(params_D1)}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D1}\n")
            f.write(f"R2 TEST: {best_R2_test_D1}\n")
            f.write(f"TIME: {time_D1}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D1}\n")
            f.close()  

    # elif way_D1 == "ALR":
    #     time_start = time.time()
    #     best_lambda_D1, best_R2_train_D1, best_R2_test_D1 = ALR.ALR_pre(x_D1, y_D1, log_filename_D1)
    #     time_D1 = time.time() - time_start

    #     with open(log_filename_D1, 'a') as f:
    #         f.write(f"BEST LAMBDA: {best_lambda_D1}\n")
    #         f.write(f"R2 TRAIN: {best_R2_train_D1}\n")
    #         f.write(f"R2 TEST: {best_R2_test_D1}\n")
    #         f.write(f"TIME: {time_D1}\n")
    #         f.close() 

    elif way_D1 == "DT":
        time_start = time.time()
        if "quadratic" in args.input_csv_D1:
            quad = True
        else:
            quad = False
        best_R2_train_D1, best_R2_test_D1, best_params = DT.hyper_param_tuning(x_D1, y_D1, log_filename_D1, quad)
        time_D1 = time.time() - time_start

        k_D1, params_D1 = best_params
        arr_D1 = DT.kbest(x_D1, y_D1, k_D1)
        selected_fv_filename_D1 = "./log/DT_" + prop_name_D1 + "_desc_norm.csv"
        output_dataset(CIDs_D1, fv_list_D1, x_D1, arr_D1, selected_fv_filename_D1)        

        with open(log_filename_D1, 'a') as f:
            f.write(f"BEST K: {k_D1}\n")
            f.write(f"BEST PARAMS: {json.dumps(params_D1)}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D1}\n")
            f.write(f"R2 TEST: {best_R2_test_D1}\n")
            f.write(f"TIME: {time_D1}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D1}\n")
            f.close()  


    if way_D2 == "Lasso":
        time_start = time.time()
        best_lambda_D2, best_R2_train_D2, best_R2_test_D2 = LR.Lasso_pre(x_D2, y_D2, log_filename_D2)
        time_D2 = time.time() - time_start

        with open(log_filename_D2, 'a') as f:
            f.write(f"BEST LAMBDA: {best_lambda_D2}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D2}\n")
            f.write(f"R2 TEST: {best_R2_test_D2}\n")
            f.write(f"TIME: {time_D2}\n")
            f.close()  
            
    elif way_D2 == "BSP":
        time_start = time.time()
        D_best_D2, best_R2_train_D2, best_R2_test_D2, best_linear_before_D2, best_quadratic_before_D2 = BSP.LLR_ALR_pre(x_D2, y_D2, fv_list_D2, log_filename_D2)
        time_D2 = time.time() - time_start

        selected_fv_filename_D2 = "./log/MLR_based_BSP_" + prop_name_D2 + "_desc_norm.csv"
        output_dataset(CIDs_D2, fv_list_D2, x_D2, D_best_D2, selected_fv_filename_D2)

        with open(log_filename_D2, 'a') as f:
            f.write(f"D_BEST: {D_best_D2}\n")
            f.write(f"D_BEST_SIZE: {len(D_best_D2)}\n")
            f.write(f"D_BEST_LINEAR_BEFORE: {best_linear_before_D2}\n")
            f.write(f"D_BEST_QUAD_BEFORE: {best_quadratic_before_D2}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D2}\n")
            f.write(f"R2 TEST: {best_R2_test_D2}\n")
            f.write(f"TIME: {time_D2}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D2}\n")
            f.close() 

    elif way_D2 == "ANN":
        time_start = time.time()
        if "quadratic" in args.input_csv_D2:
            quad = True
        else:
            quad = False
        best, best_R2_train_D2, best_R2_test_D2 = ANN.LLR_ANN_pre(x_D2, y_D2, fv_list_D2, log_filename_D2, quad)
        time_D2 = time.time() - time_start

        _, arch_D2, rho_D2, maxitr_D2 = best[0]
        arr_D2 = best[1][1]
        selected_fv_filename_D2 = "./log/LLR_ANN_" + prop_name_D2 + "_desc_norm.csv"
        output_dataset(CIDs_D2, fv_list_D2, x_D2, arr_D2, selected_fv_filename_D2)

        with open(log_filename_D2, 'a') as f:
            f.write(f"BEST ARCH: {arch_D2}\n")
            f.write(f"BEST RHO: {rho_D2}\n")
            f.write(f"BEST MAXITR: {maxitr_D2}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D1}\n")
            f.write(f"R2 TEST: {best_R2_test_D2}\n")
            f.write(f"TIME: {time_D2}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D2}\n")
            f.close()  

    elif way_D2 == "RF":
        time_start = time.time()
        if "quadratic" in args.input_csv_D2:
            quad = True
        else:
            quad = False
        best_R2_train_D2, best_R2_test_D2, best_params = RF.hyper_param_tuning(x_D2, y_D2, log_filename_D2, quad)
        time_D2 = time.time() - time_start

        k_D2, params_D2 = best_params
        arr_D2 = RF.kbest(x_D2, y_D2, k_D2)
        selected_fv_filename_D2 = "./log/RF_" + prop_name_D2 + "_desc_norm.csv"
        output_dataset(CIDs_D2, fv_list_D2, x_D2, arr_D2, selected_fv_filename_D2)        

        with open(log_filename_D2, 'a') as f:
            f.write(f"BEST K: {k_D2}\n")
            f.write(f"BEST PARAMS: {json.dumps(params_D2)}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D2}\n")
            f.write(f"R2 TEST: {best_R2_test_D2}\n")
            f.write(f"TIME: {time_D2}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D2}\n")
            f.close()  

    # elif way_D2 == "ALR":
    #     time_start = time.time()
    #     best_lambda_D2, best_R2_train_D2, best_R2_test_D2 = ALR.ALR_pre(x_D2, y_D2, log_filename_D2)
    #     time_D2 = time.time() - time_start

    #     with open(log_filename_D2, 'a') as f:
    #         f.write(f"BEST LAMBDA: {best_lambda_D2}\n")
    #         f.write(f"R2 TRAIN: {best_R2_train_D2}\n")
    #         f.write(f"R2 TEST: {best_R2_test_D2}\n")
    #         f.write(f"TIME: {time_D2}\n")
    #         f.close() 

    elif way_D2 == "DT":
        time_start = time.time()
        if "quadratic" in args.input_csv_D2:
            quad = True
        else:
            quad = False
        best_R2_train_D2, best_R2_test_D2, best_params = DT.hyper_param_tuning(x_D2, y_D2, log_filename_D2, quad)
        time_D2 = time.time() - time_start

        k_D2, params_D2 = best_params
        arr_D2 = DT.kbest(x_D2, y_D2, k_D2)
        selected_fv_filename_D2 = "./log/DT_" + prop_name_D2 + "_desc_norm.csv"
        output_dataset(CIDs_D2, fv_list_D2, x_D2, arr_D2, selected_fv_filename_D2)        

        with open(log_filename_D2, 'a') as f:
            f.write(f"BEST K: {k_D2}\n")
            f.write(f"BEST PARAMS: {json.dumps(params_D2)}\n")
            f.write(f"R2 TRAIN: {best_R2_train_D2}\n")
            f.write(f"R2 TEST: {best_R2_test_D2}\n")
            f.write(f"TIME: {time_D2}\n")
            f.write(f"SELECT DESC FILENAME: {selected_fv_filename_D2}\n")
            f.close() 

    output_str = f"{prop_name_D1}\t{x_D1.shape[0]}\t{x_D1.shape[1]}\t{nonzero_linear_org_D1}\t{nonzero_quadratic_org_D1}\t"
    output_str += f"{prop_name_D2}\t{x_D2.shape[0]}\t{x_D2.shape[1]}\t{nonzero_linear_org_D2}\t{nonzero_quadratic_org_D2}\t"

    output_str += f"{way_D1}\t{best_R2_train_D1}\t{best_R2_test_D1}\t{time_D1}\t{way_D2}\t{best_R2_train_D2}\t{best_R2_test_D2}\t{time_D2}\t\t"
    if way_D1 == "Lasso":
        output_str += f"-l1 {best_lambda_D1}\t"
    elif way_D1 == "BSP":
        output_str += f"-rbsp1 {selected_fv_filename_D1}\t"
    elif way_D1 == "ANN":
        arch_str = ""
        for arc in arch_D1:
            arch_str += f"{arc} "
        output_str += f"-ann1 {selected_fv_filename_D1} {rho_D1} {maxitr_D1} {arch_str}\t"
    elif way_D1 == "RF":
        params_str = ""
        for (params_name, value) in params_D1.items():
            params_str += f"{params_name} {value} "
        output_str += f"-rf1 {selected_fv_filename_D1} {params_str}\t"
    # elif way_D1 == "ALR":
    #     output_str += f"-alr1 {best_lambda_D1}\t"
    elif way_D1 == "DT":
        params_str = ""
        for (params_name, value) in params_D1.items():
            params_str += f"{params_name} {value} "
        output_str += f"-dt1 {selected_fv_filename_D1} {params_str}\t"    

    if way_D2 == "Lasso":
        output_str += f"-l2 {best_lambda_D2}\t"
    elif way_D2 == "BSP":
        output_str += f"-rbsp2 {selected_fv_filename_D2}\t"
    elif way_D2 == "ANN":
        arch_str = ""
        for arc in arch_D2:
            arch_str += f"{arc} "
        output_str += f"-ann2 {selected_fv_filename_D2} {rho_D2} {maxitr_D2} {arch_str}\t"
    elif way_D2 == "RF":
        params_str = ""
        for (params_name, value) in params_D2.items():
            params_str += f"{params_name} {value} "
        output_str += f"-rf2 {selected_fv_filename_D2} {params_str}\t"
    # elif way_D2 == "ALR":
    #     output_str += f"-alr2 {best_lambda_D2}\t"
    elif way_D2 == "DT":
        params_str = ""
        for (params_name, value) in params_D2.items():
            params_str += f"{params_name} {value} "
        output_str += f"-dt2 {selected_fv_filename_D2} {params_str}\t"

    print(output_str)

if __name__ == "__main__":
    # prop = "Bp_large"
    # main((0, f"./data_var0/{prop}_var0_desc_norm.csv", f"./data_var0/{prop}_norm_values.txt"))
    # main((0, "./data_regression_var0/BHL_large_var0_theta0.001_D1_quadratic_h5000_desc_norm.csv", "./data_regression_var0/BHL_large_var0_theta0.001_D1_values.txt", "./data_regression_var0/BHL_large_var0_theta0.001_D2_quadratic_h5000_desc_norm.csv", "./data_regression_var0/BHL_large_var0_theta0.001_D2_values.txt"))
    # main(sys.argv)
    main(get_args())

