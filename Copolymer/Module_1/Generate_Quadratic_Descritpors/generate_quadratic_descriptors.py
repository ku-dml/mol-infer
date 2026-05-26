############
#   -- deal with stuck issue, v4, 04/11

############### import ###############
import numpy as np
import pandas as pd

import pulp 

import time,sys,copy,itertools,math,warnings, os

from sklearn.model_selection import KFold
from sklearn.linear_model import Lasso
from sklearn.feature_selection import SelectKBest, f_regression

warnings.simplefilter('ignore')

####################################################
Times_pre = 1
Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 10000
ZERO_TOL = 0.000001

STEP = 0

MAX_NUM = 1000000

RANDOM_SEED = 11

# number of subset of descs and instances
k_p = 200
n_p = 200

# minimun dataset size to run
MIN_SIZE = 200

# list of values of lambda used in Lasso (preliminary)
lambda_list = [0.00001]
# lambda_list = [0, 0.000001, 0.00001, 0.0001, 0.001]
# lambda_list = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10, 25, 50, 100]
# lambda_list = [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 2, 5, 10, 25, 50, 100]

################################################
def learn_Lasso(x_train, y_train, x_test, y_test, a=1.0):
    lasso = Lasso(alpha=a, max_iter=10**5)
    lasso.fit(x_train, y_train)
    r2train = lasso.score(x_train,y_train)
    r2test = lasso.score(x_test,y_test)
    nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])
    return (lasso, nonzero, r2train, r2test)

def learn_Lasso_eval(x_train, y_train, fv_list, a=1.0):
    lasso = Lasso(alpha=a, max_iter=10**5)
    lasso.fit(x_train, y_train)
    # nonzero = len([w for w in lasso.coef_ if abs(w)>=ZERO_TOL])

    nonzero_linear = 0
    nonzero_quadratic = 0

    selected_fv = []

    for (i, w) in enumerate(lasso.coef_):
        if abs(w) >= ZERO_TOL:
            # print(fv_list[i])
            selected_fv.append(fv_list[i])
            if "*" in fv_list[i]:
                nonzero_quadratic += 1
            else:
                nonzero_linear += 1

    return (lasso, selected_fv, nonzero_linear, nonzero_quadratic)

################################################
def read_dataset(data_csv, value_txt):
    
    ### read files ###
    # read the csv and the observed values
    x = pd.read_csv(data_csv, index_col=0)
    value = pd.read_csv(value_txt)

    ### prepare data set ###
    # prepare CIDs
    CIDs = np.array(x.index)
    # prepare target, train, test arrays
    target = np.array(value['a'])
    # construct dictionary: CID to feature vector
    fv_dict = {}
    for cid,row in zip(CIDs, x.values):
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
    
    y = np.array([target_dict[cid] for cid in CIDs])

    return x, y

def Lasso_pre(x, y):

    # TsR = {lmd: [] for lmd in lambda_list}

    # for split_seed in range(1, Times_pre+1):
    #     # 10000 is added because this is for preliminary experiments
    #     kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

    #     for lmd in lambda_list:
    #         for train, test in kf.split(x):
    #             _, nonzero, r2train, r2test = learn_Lasso(x[train], y[train], x[test], y[test], a=lmd)
    #             TsR[lmd].append(r2test)

    # TsR_avg = {lmd: np.mean(np.array(TsR[lmd])) for lmd in TsR}

    # best_lambda = max(TsR_avg, key=TsR_avg.get)

    # # lambda_list_2 = list()
    # #
    # # if best_lambda == lambda_list[0]:
    # #     d1 = 0
    # #     d2 = lambda_list[0]
    # #     d3 = lambda_list[1]
    # #
    # #     d = (d3 - d2) / STEP
    # #     lambda_list_2 = [d2 + d * i for i in range(STEP + 1)]
    # #
    # # elif best_lambda == lambda_list[-1]:
    # #     d1 = lambda_list[-2]
    # #     d2 = lambda_list[-2]
    # #     d3 = lambda_list[-1]
    # #
    # #     d = (d3 - d2) / STEP
    # #     lambda_list_2 = [d2 + d * i for i in range(STEP + 1)]
    # # else:
    # #     idx = [i for (i, lmd) in enumerate(lambda_list) if lmd == best_lambda]
    # #     d1 = lambda_list[idx[0] - 1]
    # #     d2 = lambda_list[idx[0]]
    # #     d3 = lambda_list[idx[0] + 1]
    # #
    # #     d = (d2 - d1)/ (STEP / 2)
    # #     lambda_list_2 = [d1 + d * i for i in range(int(STEP / 2))]
    # #     d = (d3 - d2)/ (STEP / 2)
    # #     lambda_list_2_1 = [d2 + d * i for i in range(int(STEP / 2) + 1)]
    # #     lambda_list_2.extend(lambda_list_2_1)
    # #
    # # for lmd in lambda_list_2:
    # #     if lmd not in TsR:
    # #         TsR[lmd] = []
    # #
    # # for split_seed in range(1, Times_pre+1):
    # #     # 10000 is added because this is for preliminary experiments
    # #     kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)
    # #
    # #     for lmd in lambda_list_2:
    # #         if len(TsR[lmd]) == Times_pre * Fold:
    # #             continue
    # #         for train, test in kf.split(x):
    # #             _, nonzero, r2train, r2test = learn_Lasso(x[train], y[train], x[test], y[test], a=lmd)
    # #             TsR[lmd].append(r2test)

    # # TsR_avg = {lmd: np.mean(np.array(TsR[lmd])) for lmd in TsR}
    # #
    # # best_lambda = max(TsR_avg, key=TsR_avg.get)

    best_lambda = lambda_list[0]

    return best_lambda

def main(argv):
    try:
        x, y = read_dataset(argv[1], argv[2])
        h = int(argv[3])

    except:
        sys.stderr.write("usage: {} (input_data.csv)(input_values.txt)(h)\n\n".format(sys.argv[0]))
        exit(1)

    print(x.shape)

    min_y = np.amin(y)
    max_y = np.amax(y)
    y = (y - min_y) / (max_y - min_y)

    CIDs = np.array(x.index)
    N = len(CIDs)

    x_col = list(x.columns)

    desc_list = list()
    desc_list_tmp = list()
    desc_list_all = list()

    for i in range(len(x_col)):
        new_name = x_col[i]
        tmp = x[new_name]
        if abs(np.var(np.array(tmp))) >= ZERO_TOL:
            # print(i, new_name)
            desc_list_all.append(new_name)

    # xi*xj
    for i in range(len(x_col)):
        for j in range(i, len(x_col)):
            new_name = x_col[i] + "*" + x_col[j]

            tmp = np.multiply(x[x_col[i]], x[x_col[j]])
            if abs(np.var(np.array(tmp))) >= ZERO_TOL:
                # print(i, j, new_name)
                desc_list_all.append(new_name)

    # xi*(1-xj)
    for i in range(len(x_col)):
        for j in range(len(x_col)):
            new_name = x_col[i] + "*(1-" + x_col[j] + ")"
            
            tmp = np.multiply(x[x_col[i]], 1 - x[x_col[j]])
            if abs(np.var(np.array(tmp))) >= ZERO_TOL:
                # print(i, j, new_name)
                desc_list_all.append(new_name)

    print(f"S={len(desc_list_all)}")

    n_min = min(N, n_p)
    rng = np.random.default_rng(RANDOM_SEED)

    desc_list_org = desc_list_all[:]

    while len(desc_list_all) > h:
        # uncomment this if only use the last D' for kbest
        desc_list_org = desc_list_all[:]

        desc_list = list()
        
        selected_desc = rng.permutation(len(desc_list_all))

        selected_n = rng.permutation(N)
        selected_n = selected_n[:n_min]
        x_tmp = pd.DataFrame(index=x.index.take(selected_n))

        for i in selected_desc:

            fname = desc_list_all[i]
            if '*' not in fname:
                tmp = x[fname].take(selected_n)
            elif '(' not in fname:
                sp_loc = fname.find('*')
                x_name1 = fname[:sp_loc]
                x_name2 = fname[sp_loc + 1:]

                tmp = np.multiply(x[x_name1].take(selected_n), x[x_name2].take(selected_n))
            else:
                sp_loc = fname.find('*')
                x_name1 = fname[:sp_loc]
                x_name2 = fname[sp_loc + 4: -1]

                tmp = np.multiply(x[x_name1].take(selected_n), 1 - x[x_name2].take(selected_n))

            desc_list_tmp.append(fname)
            x_tmp[fname] = tmp

            if len(desc_list_tmp) == k_p:
                # time_start = time.time()
                best_lambda = Lasso_pre(x_tmp.values, y[selected_n])
                _, selected_fv, _, _ = learn_Lasso_eval(x_tmp.values, y[selected_n], desc_list_tmp, a=best_lambda)
                # print(time.time() - time_start)
                
                desc_list.extend(selected_fv)
                # print(selected_fv)

                selected_n = rng.permutation(N)
                selected_n = selected_n[:n_min]
                x_tmp = pd.DataFrame(index=x.index.take(selected_n))
                desc_list_tmp = list()

        if len(desc_list_tmp) != 0:
            best_lambda = Lasso_pre(x_tmp.values, y[selected_n])
            _, selected_fv, _, _ = learn_Lasso_eval(x_tmp.values, y[selected_n], desc_list_tmp, a=best_lambda)

            desc_list.extend(selected_fv)
            # print(selected_fv)

            desc_list_tmp = list()

        print(f"in while loop, |D'|={len(desc_list)}")

        if len(desc_list) == len(desc_list_all):
            desc_list_all = desc_list[:]
            break

        desc_list_all = desc_list[:]

    desc_list = desc_list_all[:]
    print(f"|D'|={len(desc_list)}")

    # filename = argv[1]
    # sp_loc = filename.find('_desc_norm')
    # sp_loc2 = filename_paper.find('.')
    # sp_str = str()
    # for a in des_indices:
    #     sp_str += '_' + str(a)
    # output_filename = filename[:sp_loc] + f"_indices_part1.txt"
    # with open(output_filename, 'w') as f:
    #     for _x in desc_list:
    #         f.write(f"{_x}\n")

    x_new = pd.DataFrame(index=x.index)
    for i in range(len(desc_list)):
        if '*' not in desc_list[i]:
            x_new[desc_list[i]] = x[desc_list[i]]
        elif '(' not in desc_list[i]:
            sp_loc = desc_list[i].find('*')
            x_name1 = desc_list[i][:sp_loc]
            x_name2 = desc_list[i][sp_loc + 1:]
            # tmp = list()
            # for cid in CIDs:
            #     tmp.append(x.loc[cid, x_name1] * x.loc[cid, x_name2])
            # tmp = np.array(tmp)
            tmp = np.multiply(x[x_name1], x[x_name2])
            x_new[desc_list[i]] = tmp
        else:
            sp_loc = desc_list[i].find('*')
            x_name1 = desc_list[i][:sp_loc]
            x_name2 = desc_list[i][sp_loc + 4: -1]
            # tmp = list()
            # for cid in CIDs:
            #     tmp.append(x.loc[cid, x_name1] * (1 - x.loc[cid, x_name2]))
            # tmp = np.array(tmp)
            tmp = np.multiply(x[x_name1], 1 - x[x_name2])
            x_new[desc_list[i]] = tmp

    if len(desc_list) >= h:
        print("|D'|>h")
        kbest = SelectKBest(f_regression, k=h)
        kbest.fit(x_new.values, y)
        arr = kbest.get_support(indices=True)
        # print(arr)
        selected_desc = [desc_list[i] for i in arr]
        x_new = x_new[selected_desc]
        # print(x_new[selected_desc])
    elif len(desc_list) < h:
        print("|D'|<h")
        desc_list_tmp = list()
        x_tmp = pd.DataFrame(index=x.index)

        if len(desc_list) != len(desc_list_org):
            # print(len(desc_list_org))

            for i in desc_list_org:
                fname = i
                if fname in desc_list:
                    continue
                if '*' not in fname:
                    tmp = x[fname]
                elif '(' not in fname:
                    sp_loc = fname.find('*')
                    x_name1 = fname[:sp_loc]
                    x_name2 = fname[sp_loc + 1:]

                    tmp = np.multiply(x[x_name1], x[x_name2])
                else:
                    sp_loc = fname.find('*')
                    x_name1 = fname[:sp_loc]
                    x_name2 = fname[sp_loc + 4: -1]

                    tmp = np.multiply(x[x_name1], 1 - x[x_name2])

                desc_list_tmp.append(fname)
                x_tmp[fname] = tmp

            # print("x_tmp:", x_tmp.values.shape)
            kbest = SelectKBest(f_regression, k=h-len(desc_list))
            kbest.fit(x_tmp.values, y)
            arr = kbest.get_support(indices=True)

            for _i, i in enumerate(arr):
                # print(_i + 1, desc_list_tmp[i])
                # if desc_list_tmp[i] in desc_list:
                #     print(_i, u, desc_list_tmp[i])
                x_new[desc_list_tmp[i]] = x_tmp.values[:, i]

    selected_linear = 0
    selected_quad = 0

    for fname in list(x_new.columns):
        if '*' in fname:
            selected_quad += 1
        else:
            selected_linear += 1

    print(f"#selected 1-D desc: {selected_linear}, #selected 2-D desc: {selected_quad}.")


    # # normalized_x = x

    # # print(normalized_x["dg_1"])

    filename = argv[1]
    sp_loc = filename.find('_desc_norm')
    # sp_loc2 = filename_paper.find('.')
    # sp_str = str()
    # for a in des_indices:
    #     sp_str += '_' + str(a)
    output_filename = filename[:sp_loc] + f"_quadratic_h{h}" + filename[sp_loc:]

    print(x_new.values.shape)

    x_new.to_csv(output_filename, sep=',')

if __name__ == "__main__":
    # R_DIR = "./classification"
    # fnames = os.listdir(path=R_DIR)
    # fnames.sort()
    # for fname in fnames:
    #     if '.csv' not in fname or '._' in fname:
    #         continue
    #     filename = os.path.join(R_DIR, fname)
    #     print(filename)
    #     main((0, filename))
    main(sys.argv)
    # prop = "At-Organic_org"
    # main((0, f"RefIdx_3elem_var0_desc_norm.csv", "RefIdx_norm_values.txt", 20000))

    