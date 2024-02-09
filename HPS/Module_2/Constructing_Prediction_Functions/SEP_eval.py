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

from src import LR, ALR, BSP, ANN, RF
from src import EVAL

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
    parser.add_argument('-l1', '--lasso1', nargs='+', type=str)
    parser.add_argument('-rbsp1', '--rbsp1', nargs='+', type=str)
    parser.add_argument('-ann1', '--ann1', nargs='+', type=str)
    parser.add_argument('-rf1', '--rf1', nargs='+', type=str)
    parser.add_argument('-alr1', '--alr1', nargs='+', type=str)
    parser.add_argument('-dt1', '--dt1', nargs='+', type=str)
    
    parser.add_argument('-l2', '--lasso2', nargs='+', type=str)
    parser.add_argument('-rbsp2', '--rbsp2', nargs='+', type=str)
    parser.add_argument('-ann2', '--ann2', nargs='+', type=str)
    parser.add_argument('-rf2', '--rf2', nargs='+', type=str)
    parser.add_argument('-alr2', '--alr2', nargs='+', type=str)
    parser.add_argument('-dt2', '--dt2', nargs='+', type=str)
    

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

    params_D1 = {}
    params_D2 = {}

    try:
        way_D1 = ""
        if args.lasso1 is not None:
            way_D1 = "Lasso"
            params_D1["lambda"] = float(args.lasso1[0])
        elif args.rbsp1 is not None:
            way_D1 = "BSP"
            params_D1["selected_fv_filename"] = args.rbsp1[0]
        elif args.ann1 is not None:
            way_D1 = "ANN"
            params_D1["selected_fv_filename"] = args.ann1[0]
            params_D1["rho"] = float(args.ann1[1])
            params_D1["maxitr"] = int(args.ann1[2])
            params_D1["arch"] = list(map(int, args.ann1[3:]))
        elif args.rf1 is not None:
            way_D1 = "RF"
            params_D1["selected_fv_filename"] = args.rf1[0]
            for i in range(1, len(args.rf1), 2):
                if args.rf1[i + 1] == "None":
                    params_D1[args.rf1[i]] = None
                elif args.rf1[i + 1] == "sqrt":
                    params_D1[args.rf1[i]] = "sqrt"
                elif "." in args.rf1[i + 1]:
                    params_D1[args.rf1[i]] = float(args.rf1[i + 1])
                else:
                    params_D1[args.rf1[i]] = int(args.rf1[i + 1])
        # elif args.alr1 is not None:
        #     way_D1 = "ALR"
        #     params_D1["lambda"] = float(args.alr1[0])
        elif args.dt1 is not None:
            way_D1 = "DT"
            params_D1["selected_fv_filename"] = args.dt1[0]
            for i in range(1, len(args.dt1), 2):
                if args.dt1[i + 1] == "None":
                    params_D1[args.dt1[i]] = None
                elif args.dt1[i + 1] == "sqrt":
                    params_D1[args.dt1[i]] = "sqrt"
                elif "." in args.dt1[i + 1]:
                    params_D1[args.dt1[i]] = float(args.dt1[i + 1])
                else:
                    params_D1[args.dt1[i]] = int(args.dt1[i + 1])
        else:
            sys.stderr.write("Please specify one learning method for D1\n")
            exit(1)

        way_D2 = ""
        if args.lasso2 is not None:
            way_D2 = "Lasso"
            params_D2["lambda"] = float(args.lasso2[0])
        elif args.rbsp2 is not None:
            way_D2 = "BSP"
            params_D2["selected_fv_filename"] = args.rbsp2[0]
        elif args.ann2 is not None:
            way_D2 = "ANN"
            params_D2["selected_fv_filename"] = args.ann2[0]
            params_D2["rho"] = float(args.ann2[1])
            params_D2["maxitr"] = int(args.ann2[2])
            params_D2["arch"] = list(map(int, args.ann2[3:]))
        elif args.rf2 is not None:
            way_D2 = "RF"
            params_D2["selected_fv_filename"] = args.rf2[0]
            for i in range(1, len(args.rf2), 2):
                if args.rf2[i + 1] == "None":
                    params_D2[args.rf2[i]] = None
                elif args.rf2[i + 1] == "sqrt":
                    params_D2[args.rf2[i]] = "sqrt"
                elif "." in args.rf2[i + 1]:
                    params_D2[args.rf2[i]] = float(args.rf2[i + 1])
                else:
                    params_D2[args.rf2[i]] = int(args.rf2[i + 1])
        # elif args.alr2 is not None:
        #     way_D2 = "ALR"
        #     params_D1["lambda"] = float(args.alr2[0])
        elif args.dt2 is not None:
            way_D2 = "DT"
            params_D2["selected_fv_filename"] = args.dt2[0]
            for i in range(1, len(args.dt2), 2):
                if args.dt2[i + 1] == "None":
                    params_D2[args.dt2[i]] = None
                elif args.dt2[i + 1] == "sqrt":
                    params_D2[args.dt2[i]] = "sqrt"
                elif "." in args.dt2[i + 1]:
                    params_D2[args.dt2[i]] = float(args.dt2[i + 1])
                else:
                    params_D2[args.dt2[i]] = int(args.dt2[i + 1])
        else:
            sys.stderr.write("Please specify one learning method for D2\n")
            exit(1)

        params_D1["way"] = way_D1
        params_D2["way"] = way_D2

        if "selected_fv_filename" in params_D1:
            CIDs_D1, x_D1, y_D1, fv_list_D1 = read_dataset(params_D1["selected_fv_filename"], args.input_value_D1)
        else:
            CIDs_D1, x_D1, y_D1, fv_list_D1 = read_dataset(args.input_csv_D1, args.input_value_D1)
        if "selected_fv_filename" in params_D2:
            CIDs_D2, x_D2, y_D2, fv_list_D2 = read_dataset(params_D2["selected_fv_filename"], args.input_value_D2)
        else:
            CIDs_D2, x_D2, y_D2, fv_list_D2 = read_dataset(args.input_csv_D2, args.input_value_D2)

    except:
        sys.stderr.write("usage: {} (input_csv_D1.csv)(input_value_D1.txt)(input_csv_D2.csv)(input_value_D2.txt)\n\n".format(sys.argv[0]))
        exit(1)

    x_D1 = x_D1.astype(np.float64)
    y_D1 = y_D1.astype(np.float64)
    x_D2 = x_D2.astype(np.float64)
    y_D2 = y_D2.astype(np.float64)

    # # re-normalization, only for pre
    # for eval, not here
    # min_y = np.amin(y_D1)
    # max_y = np.amax(y_D1)
    # y_D1 = (y_D1 - min_y) / (max_y - min_y)

    # min_y = np.amin(y_D2)
    # max_y = np.amax(y_D2)
    # y_D2 = (y_D2 - min_y) / (max_y - min_y)

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

    log_filename_D1 = "./log_eval/SEP_eval_" + prop_name_D1 + "_" + way_D1 + "_log.txt"
    with open(log_filename_D1, 'w') as f:
        f.close()

    log_filename_D2 = "./log_eval/SEP_eval_" + prop_name_D2 + "_" + way_D2 + "_log.txt"
    with open(log_filename_D2, 'w') as f:
        f.close()

    output_prefix_D1 = f"./linreg/{prop_name_D1}_"
    output_prefix_D2 = f"./linreg/{prop_name_D2}_"

    R2_train, R2_test, R2_train_D1, R2_test_D1, R2_train_D2, R2_test_D2, sum_time_D1, sum_time_D2 = \
        EVAL.eval_separate_plane(
            x_D1, y_D1, fv_list_D1, log_filename_D1, params_D1, output_prefix_D1,
            x_D2, y_D2, fv_list_D2, log_filename_D2, params_D2, output_prefix_D2
        )

    output_str = f"{prop_name_D1}\t{x_D1.shape[0]}\t{x_D1.shape[1]}\t{nonzero_linear_org_D1}\t{nonzero_quadratic_org_D1}\t"
    output_str += f"{prop_name_D2}\t{x_D2.shape[0]}\t{x_D2.shape[1]}\t{nonzero_linear_org_D2}\t{nonzero_quadratic_org_D2}\t"
    output_str += f"{way_D1}\t{np.median(R2_train_D1)}\t{np.amin(R2_train_D1)}\t{np.amax(R2_train_D1)}\t"
    output_str += f"{np.median(R2_test_D1)}\t{np.amin(R2_test_D1)}\t{np.amax(R2_test_D1)}\t{sum_time_D1}\t"
    output_str += f"{way_D2}\t{np.median(R2_train_D2)}\t{np.amin(R2_train_D2)}\t{np.amax(R2_train_D2)}\t"
    output_str += f"{np.median(R2_test_D2)}\t{np.amin(R2_test_D2)}\t{np.amax(R2_test_D2)}\t{sum_time_D2}\t"
    output_str += f"ALL\t{np.median(R2_train)}\t{np.amin(R2_train)}\t{np.amax(R2_train)}\t"
    output_str += f"{np.median(R2_test)}\t{np.amin(R2_test)}\t{np.amax(R2_test)}\t{sum_time_D1 + sum_time_D2}\t"

    print(output_str)

if __name__ == "__main__":
    # prop = "Bp_large"
    # main((0, f"./data_var0/{prop}_var0_desc_norm.csv", f"./data_var0/{prop}_norm_values.txt"))
    # main((0, "./data_regression_var0/BHL_large_var0_theta0.001_D1_quadratic_h5000_desc_norm.csv", "./data_regression_var0/BHL_large_var0_theta0.001_D1_values.txt", "./data_regression_var0/BHL_large_var0_theta0.001_D2_quadratic_h5000_desc_norm.csv", "./data_regression_var0/BHL_large_var0_theta0.001_D2_values.txt"))
    # main(sys.argv)
    main(get_args())

