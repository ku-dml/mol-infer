import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.tree import DecisionTreeRegressor
from sklearn.feature_selection import SelectKBest, f_regression


import time,sys,copy,itertools,math,warnings,datetime,json

warnings.simplefilter('ignore')

####################################################
SEARCH_PARAMS = {
    'max_features': [1, 2, 3, 10, "sqrt", 1.0],
    'min_samples_split': [2, 3, 5, 10],
    'min_samples_leaf': [1, 2, 3, 5, 10],
    'max_depth': [2, 5, 10, 12, 15, 17, 20, None],
    'ccp_alpha': [0.000, 0.001, 0.002, 0.005]
}

K_CANDIDATE = [1000, 5000]

Times_pre = 5
Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 1000

####################################################
def learn_decision_tree(x_train, y_train, x_test, y_test, params=DecisionTreeRegressor().get_params()):
    dt = DecisionTreeRegressor()
    dt.set_params(**params)
    dt.fit(x_train, y_train)
    r2train = dt.score(x_train, y_train)
    r2test = dt.score(x_test, y_test)
    return r2train, r2test, dt.predict(x_train), dt.predict(x_test)

# kbest returns selected_dataset x and mask
def kbest(x, y, k):
    selector = SelectKBest(score_func=f_regression, k=k)
    selector.fit(x, y)
    # x = selector.transform(x)
    return selector.get_support(indices=True)

def hyper_param_tuning(x, y, log_filename, quad=True):
    best_train_score = None
    best_test_score = None
    best_params = None

    fv_size = x.shape[1]

    if quad:
        k_cand = [fv_size]
    else:
        k_cand = [math.floor(fv_size * i / 10) for i in range(5, 11)]

    for k in k_cand:
        arr = kbest(x, y, k)
        x_selected = x[:, arr]

        keys = list(SEARCH_PARAMS)
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
                    r2train, r2test, _, _ = learn_decision_tree(x_selected[train], y[train], x_selected[test], y[test], params=params)
                    time_end = time.time()
                    R2_train.append(r2train)
                    R2_test.append(r2test)

                    with open(log_filename, 'a') as f:
                        f.write(f"{split_seed}, {r2train}, {r2test}, {time_end - time_start}\n")
                        f.close()

            with open(log_filename, 'a') as f:
                f.write(f"median train score: {np.median(R2_train)}\n")
                f.write(f"median test score: {np.median(R2_test)}\n")
                f.write(f"time: {time.time() - time_start_org}\n\n")
                f.close()

            if best_test_score is None:
                best_train_score = np.median(R2_train)
                best_test_score = np.median(R2_test)
                best_params = (k, params)
            elif np.median(R2_test) > best_test_score:
                best_train_score = np.median(R2_train)
                best_test_score = np.median(R2_test)
                best_params = (k, params)                

    return best_train_score, best_test_score, best_params 

            
