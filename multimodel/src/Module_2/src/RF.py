import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_regression

from src import CLSF, FUN

import time,sys,copy,itertools,math,warnings,datetime,json

warnings.simplefilter('ignore')

####################################################
SEARCH_PARAMS = {
    'max_features': [1, 2, 3, 10, "sqrt", 1.0],
    'min_samples_split': [2, 3, 5, 10],
    'min_samples_leaf': [1, 2, 3, 5, 10],
    'max_depth': [2, 5, 10, 12, 15, 17, 20, None],
    'n_estimators': [5, 10, 25, 50, 75, 100]
}

K_CANDIDATE = [1000, 5000]

Times_pre = 5
Times_eva = 10
Fold = 5
RANDOM_SEED_BASE = 1000

####################################################
def learn_random_forest(
    x_train, y_train, x_test, y_test, params=RandomForestRegressor().get_params(), metric='R2'
):
    if metric == 'R2':
        randomforest = RandomForestRegressor()
        randomforest.set_params(**params)
        randomforest.fit(x_train, y_train)
    elif metric == 'BACC' or metric == 'AUCROC':
        randomforest = RandomForestClassifier()
        randomforest.set_params(**params)
        randomforest.set_params(class_weight="balanced")
        randomforest.fit(x_train, y_train)

    if metric == 'R2':
        r2train = randomforest.score(x_train, y_train)
        r2test = randomforest.score(x_test, y_test)
    elif metric == 'BACC':
        r2train, r2test = CLSF.calc_bacc(randomforest, x_train, y_train, x_test, y_test)
    elif metric == 'AUCROC':
        r2train, r2test = CLSF.calc_aucroc(randomforest, x_train, y_train, x_test, y_test)

    return randomforest, r2train, r2test, randomforest.predict(x_train), randomforest.predict(x_test)

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
                    _, r2train, r2test, _, _ = learn_random_forest(x_selected[train], y[train], x_selected[test], y[test], params=params, metric=metric)
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

def write_file(rf, output_filename):
    with open(output_filename, 'w') as f:
        f.write(f"{rf.get_params()['n_estimators']}\n")

        for tree in rf.estimators_:
            # set up information of the tree
            children_left = tree.tree_.children_left
            children_right = tree.tree_.children_right
            feature = tree.tree_.feature
            threshold = tree.tree_.threshold
            value = tree.tree_.value[:, 0, 0]
            # lists of nodes
            leaf_nodes = []
            inner_nodes = []
            stack = [0]
            while len(stack) > 0:
                node_id = stack.pop()
                left_node_id = children_left[node_id]
                right_node_id = children_right[node_id]
                if left_node_id != right_node_id:
                    inner_nodes.append(node_id)
                    stack.append(left_node_id)
                    stack.append(right_node_id)
                else:
                    leaf_nodes.append(node_id)

            # write inner nodes information
            f.write(f"{len(inner_nodes)}\n")
            while len(inner_nodes) > 0:
                node_id = inner_nodes.pop()
                left_node_id = children_left[node_id]
                right_node_id = children_right[node_id]
                f.write(
                    f'{node_id} {left_node_id} {right_node_id} {feature[node_id]} {threshold[node_id]}\n')

            # write leaf nodes information
            f.write(f"{len(leaf_nodes)}\n")
            while len(leaf_nodes) > 0:
                node_id = leaf_nodes.pop()
                f.write(f'{node_id} {value[node_id]}\n')




