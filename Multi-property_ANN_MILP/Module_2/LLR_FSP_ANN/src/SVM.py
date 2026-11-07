### This part is devoted to only linear kernel ###

import numpy as np
import pandas as pd
from sklearn.svm import LinearSVR, LinearSVC
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.feature_selection import SelectKBest, f_regression

import time,sys,copy,itertools,math,warnings,datetime,json

from src import CLSF, FUN

C_list = list(np.arange(0.0001, 0.001, 0.0001)) + list(np.arange(0.001, 0.01, 0.001)) + \
         list(np.arange(0.01, 0.1, 0.01)) + list(np.arange(0.1, 1.0, 0.05)) + \
         list(np.arange(1.0, 10.0, 0.5)) + list(np.arange(10.0, 100.0, 10.0))

SEARCH_PARAMS = {
    'C': C_list,
}
STEP = 20

K_CANDIDATE = [1000, 5000]

Times_pre = 5
Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 10000
ZERO_TOL = 0.000001

warnings.simplefilter('ignore')

####################################################
def learn_SVM(
    x_train, y_train, x_test, y_test, params=LinearSVR().get_params(), metric='R2'
):
    if metric == 'R2':
        svm = LinearSVR()
        svm.set_params(**params)
        svm.fit(x_train, y_train)
    elif metric == 'BACC' or metric == 'AUCROC':
        svm = LinearSVC()
        svm.set_params(**params)
        svm.set_params(class_weight="balanced")
        svm.fit(x_train, y_train)

    if metric == 'R2':
        r2train = svm.score(x_train, y_train)
        r2test = svm.score(x_test, y_test)
    elif metric == 'BACC':
        r2train, r2test = CLSF.calc_bacc(svm, x_train, y_train, x_test, y_test)
    elif metric == 'AUCROC':
        r2train, r2test = CLSF.calc_aucroc(svm, x_train, y_train, x_test, y_test)

    return r2train, r2test, svm.predict(x_train), svm.predict(x_test)

# kbest returns selected_dataset x and mask
def kbest(x, y, k):
    selector = SelectKBest(score_func=f_regression, k=k)
    selector.fit(x, y)
    # x = selector.transform(x)
    return selector.get_support(indices=True)

def hyper_param_tuning(x, y, log_filename, quad=True, metric='R2', Sbar=FUN.Status_Bar()):
    best_train_score = None
    best_test_score = None
    best_params = None

    fv_size = x.shape[1]

    if quad:
        k_cand = [fv_size]
    else:
        k_cand = [math.floor(fv_size * i / 10) for i in range(5, 11)]

    keys = list(SEARCH_PARAMS)
    Sbar.set_total_times(len(k_cand) * len(list(itertools.product(*map(SEARCH_PARAMS.get, keys)))) * Times_pre * Fold)
    Sbar.output()
    
    for k in k_cand:
        arr = kbest(x, y, k)
        x_selected = x[:, arr]
        
        for values in itertools.product(*map(SEARCH_PARAMS.get, keys)):
            params = dict(zip(keys, values))

            R2_train = []
            R2_test = []

            with open(log_filename, 'a') as f:
                f.write(f"k: {k}\n")
                f.write(f"params: {json.dumps(params)}\n")
                f.close()

            time_start_org = time.time()

            for split_seed in range(1, Times_pre+1):
                kf = KFold(n_splits=Fold, shuffle=True, random_state=RANDOM_SEED_BASE+split_seed)

                for train, test in kf.split(x_selected):
                    time_start = time.time()
                    r2train, r2test, _, _ = learn_SVM(x_selected[train], y[train], x_selected[test], y[test], params=params, metric=metric)
                    time_end = time.time()
                    R2_train.append(r2train)
                    R2_test.append(r2test)

                    with open(log_filename, 'a') as f:
                        f.write(f"{split_seed}, {r2train}, {r2test}, {time_end - time_start}\n")
                        f.close()

                    Sbar.append(time_end - time_start)

            with open(log_filename, 'a') as f:
                f.write(f"median train score: {np.median(R2_train)}\n")
                f.write(f"median test score: {np.median(R2_test)}\n")
                f.write(f"time: {time.time() - time_start_org}\n\n")
                f.close()

            Sbar.output()

            if best_test_score is None:
                best_train_score = np.median(R2_train)
                best_test_score = np.median(R2_test)
                best_params = (k, params)
            elif np.median(R2_test) > best_test_score:
                best_train_score = np.median(R2_train)
                best_test_score = np.median(R2_test)
                best_params = (k, params)                

    return best_train_score, best_test_score, best_params 

            
