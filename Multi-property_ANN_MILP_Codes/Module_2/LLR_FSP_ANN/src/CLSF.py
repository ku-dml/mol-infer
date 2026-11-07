import numpy as np
from sklearn.metrics import r2_score, roc_auc_score
import copy

ZERO_TOL = 0.000001

def get_theta(f, x, y, y_predict=None):
    if y_predict is None:
        R = f.predict(x)
    else:
        R = copy.deepcopy(y_predict)
    I = np.argsort(R)
    R.sort()
    ###
    C1 = len([cls for cls in y if cls==1])
    C0 = len(y)-C1
    Pos = [-1,]
    TP = C1
    TN = C0
    best_balacc = 0.5
    ###
    k = 0
    while k<len(y):
        K = [k,]
        while k<len(y)-1 and R[k]==R[k+1]:
            k += 1
            K.append(k)
        for l in K:
            if y[I[l]] == 1:
                TP += -1
            else:
                TN += 1
        balacc = 0.5 * (TP/C1 + TN/C0)
        if balacc > best_balacc:
            Pos = [K[-1],]
            best_balacc = balacc
        elif balacc == best_balacc:
            Pos.append(K[-1])
        k += 1
    ###
    # determine theta
    pos = Pos[len(Pos)//2]
    if pos == -1:
        theta = R[0] - ZERO_TOL
    elif pos == len(y)-1:
        theta = R[-1] + ZERO_TOL
    else:
        theta = (R[pos]+R[pos+1])/2

    '''
    ### the following verification schemes do not work for balanced accuracy without modification ###
    # verification I: whether acc == number of class 0
    zero_class = len([cls for cls in y if cls==0])
    if acc != zero_class:
        sys.stderr.write(f"error: zero_class {zero_class} != acc {acc}\n")
        exit(1)
        
    # verification II: whether theta works correctly
    verif_acc = 0
    for k,v in enumerate(R):
        i = I[k]
        if v<theta and y[i]==0:
            verif_acc += 1
        elif v>=theta and y[i]==1:
            verif_acc += 1
    if best_acc != verif_acc:
        sys.stderr.write(f"warning: best_acc {best_acc} != verif_acc {verif_acc}\n")
    '''
        
    return theta

def calc_bacc(reg, x_train, y_train, x_test, y_test, y_predict_train=None, y_predict_test=None):
    theta = get_theta(reg, x_train, y_train, y_predict_train)

    if y_predict_train is None:
        P = reg.predict(x_train)
    else:
        P = copy.deepcopy(y_predict_train)
    # aucroc_train = roc_auc_score(y_train, P)
    TParr = len([i for i in range(len(x_train)) if P[i]>=theta and y_train[i]==1])
    FNarr = len([i for i in range(len(x_train)) if P[i]< theta and y_train[i]==1])
    TNarr = len([i for i in range(len(x_train)) if P[i]< theta and y_train[i]==0])
    FParr = len([i for i in range(len(x_train)) if P[i]>=theta and y_train[i]==0])
    TPR_test = TParr / (TParr + FNarr)
    TNR_test = TNarr / (TNarr + FParr)
    # acc_test = (TParr + TNarr) / (TParr + FNarr + TNarr + FParr)
    bacctrain = (TPR_test + TNR_test) / 2

    if y_predict_test is None:
        P = reg.predict(x_test)
    else:
        P = copy.deepcopy(y_predict_test)
    # aucroc_test = roc_auc_score(y_test, P)
    TParr = len([i for i in range(len(x_test)) if P[i]>=theta and y_test[i]==1])
    FNarr = len([i for i in range(len(x_test)) if P[i]< theta and y_test[i]==1])
    TNarr = len([i for i in range(len(x_test)) if P[i]< theta and y_test[i]==0])
    FParr = len([i for i in range(len(x_test)) if P[i]>=theta and y_test[i]==0])
    TPR_test = TParr / (TParr + FNarr)
    TNR_test = TNarr / (TNarr + FParr)
    # acc_test = (TParr + TNarr) / (TParr + FNarr + TNarr + FParr)
    bacctest = (TPR_test + TNR_test) / 2

    return bacctrain, bacctest

def calc_aucroc(reg, x_train, y_train, x_test, y_test, y_predict_train=None, y_predict_test=None):
    if y_predict_train is None:
        P = reg.predict(x_train)
    else:
        P = copy.deepcopy(y_predict_train)
    aucroc_train = roc_auc_score(y_train, P)

    if y_predict_test is None:
        P = reg.predict(x_test)
    else:
        P = copy.deepcopy(y_predict_test)
    aucroc_test = roc_auc_score(y_test, P)    

    return aucroc_train, aucroc_test

