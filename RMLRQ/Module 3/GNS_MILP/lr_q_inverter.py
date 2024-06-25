# -*- coding: utf-8 -*-

"""
This file implements
 - a class to store the architecture of a linear regressor (LR)
 - functions to initialize vatiables and build a pulp MILP model for
   inverting a DT of a given architecture

Author:    Discrete Mathematics Lab,
           Department of Applied Mathematics and Physics,
           Graduate School of Informatics,
           Kyoto University
"""

# standard imports
import sys
import pandas as pd
import numpy as np
import subprocess

# import pulp for MILP modeling
# import pulp
import pulp_modified as pulp

# the splitter in the dt file
SPL = ' '

class LinearReg_q:
    def __init__(self):
        self.K = 0
        self.K_l = 0
        self.K_q = 0
        self.K_l_ind = list()
        self.K_q_ind = list()
        self.p = 6
        self.p_q = self.p + 3
        self.coef = list()
        self.bias = 0
        self.x = dict()
        self.y = dict()
        self.x_q = dict()
        self.z_q = dict()
        self.z_tmp = dict()
        self.z = dict()
        self.weight_var = dict()
        self.predict = dict()

        self.M = 10000

    def set_p(self, new_p):
        self.p = new_p
        self.p_q = self.p + 3

    def build_var(self, prop="def"):
        # self.x = {i: pulp.LpVariable(f"LRq_x({i})", 0, 1, cat=pulp.LpContinuous) for i in range(1, self.K_q + 1)}
        # self.y = {i: pulp.LpVariable(f"LRq_y({i})", 0, 1, cat=pulp.LpContinuous) for i in range(1, self.K_q + 1)}
        # self.x_q = {(i, q): pulp.LpVariable(f"LRq_xq({i},{q})", 0, 1, cat=pulp.LpBinary) 
        #                 for i in range(1, self.K_q + 1) for q in range(self.p + 1)}
        # self.z = {i: pulp.LpVariable(f"LRq_z({i})", 0, 1, cat=pulp.LpContinuous) for i in range(1, self.K_q + 1)}
        # self.z_q = {(i, q): pulp.LpVariable(f"LRq_zq({i},{q})", 0, 1, cat=pulp.LpContinuous) 
        #                 for i in range(1, self.K_q + 1) for q in range(self.p + 1)}

        # self.weight_var = {i: pulp.LpVariable(f"LR_x_{i}", cat=pulp.LpContinuous) for i in range(1, self.K_l + 1)}
        # self.predict = {0: pulp.LpVariable(f"LR_y", cat=pulp.LpContinuous)}

        self.x = {i: pulp.LpVariable(f"LRq_x({i})_{prop}", 0, cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.y = {i: pulp.LpVariable(f"LRq_y({i})_{prop}", 0, cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.x_q = {(i, q): pulp.LpVariable(f"LRq_xq({i},{q})_{prop}", 0, 1, cat=pulp.LpBinary) 
                        for i in self.K_q_ind for q in range(self.p_q + 1)}
        self.z_tmp = {i: pulp.LpVariable(f"LRq_z_tmp({i})_{prop}", 0, cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.z = {i: pulp.LpVariable(f"LRq_z({i})_{prop}", cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.z_q = {(i, q): pulp.LpVariable(f"LRq_zq({i},{q})_{prop}", 0, cat=pulp.LpContinuous) 
                        for i in self.K_q_ind for q in range(self.p_q + 1)}

        self.weight_var = {i: pulp.LpVariable(f"LR_x_{i}_{prop}", cat=pulp.LpContinuous) for i in range(1, self.K_l + 1)}
        self.weight_var[0] = 0
        # self.predict = {0: pulp.LpVariable(f"LR_y_{prop}", 0, 1, cat=pulp.LpContinuous)}
        self.predict = {0: pulp.LpVariable(f"LR_y_{prop}", cat=pulp.LpContinuous)}

        return self.weight_var, self.x, self.y, self.predict[0]

    def build_constraints(self, MILP, descriptors_q_x, descriptors_q_y, descriptors_q_x_minus, descriptors_q_y_minus, prop="def"):
        for i in range(1, self.K_q + 1):
            q_ind = self.K_q_ind[i - 1]
            MILP += (2**(self.p + 1) - 1) * self.x[q_ind] >= \
                    pulp.lpSum([2**q * self.x_q[(q_ind, q)] for q in range(self.p_q + 1)]) - 1, f"LRq-{i}-1_{prop}"
            MILP += (2**(self.p + 1) - 1) * self.x[q_ind] <= \
                    pulp.lpSum([2**q * self.x_q[(q_ind, q)] for q in range(self.p_q + 1)]), f"LRq-{i}-2_{prop}"
            for q in range(self.p_q + 1):
                MILP += self.z_q[(q_ind, q)] <= self.M * self.x_q[(q_ind, q)], f"LRq-{i}-3-{q}_{prop}"
                MILP += self.z_q[(q_ind, q)] >= self.y[q_ind] - self.M * (1 - self.x_q[(q_ind, q)]), f"LRq-{i}-4-{q}_{prop}"
                MILP += self.z_q[(q_ind, q)] <= self.y[q_ind] + self.M * (1 - self.x_q[(q_ind, q)]), f"LRq-{i}-5-{q}_{prop}"
            MILP += self.z_tmp[q_ind] * (2**(self.p + 1) - 1) == \
                    pulp.lpSum([2**q * self.z_q[(q_ind, q)] for q in range(self.p_q + 1)]), f"LRq-{i}-6_{prop}"

            if descriptors_q_x_minus[i] and descriptors_q_y_minus[i]:
                MILP += self.z[q_ind] == 1 - self.weight_var[descriptors_q_x[i]] - self.weight_var[descriptors_q_y[i]] + self.z_tmp[q_ind], f"LRq-{i}-7_{prop}"
            elif descriptors_q_x_minus[i]:
                MILP += self.z[q_ind] == self.weight_var[descriptors_q_y[i]] - self.z_tmp[q_ind], f"LRq-{i}-8_{prop}"
            elif descriptors_q_y_minus[i]:
                MILP += self.z[q_ind] == self.weight_var[descriptors_q_x[i]] - self.z_tmp[q_ind], f"LRq-{i}-9_{prop}"
            else:
                MILP += self.z[q_ind] == self.z_tmp[q_ind], f"LRq-{i}-10_{prop}"

        MILP += self.predict[0] == \
                    pulp.lpSum([self.coef[j - 1] * self.weight_var[i + 1] for (i, j) in enumerate(self.K_l_ind)]) + \
                    pulp.lpSum([self.coef[i - 1] * self.z[i] for i in self.K_q_ind]) + \
                    self.bias, f"LRq_lr_{prop}"

    def build_constraints_y(self, MILP, y_lb, y_ub, prop="def"):
        if f"LRq_ub_{prop}" in MILP.constraints:
            MILP.constraints[f"LRq_ub_{prop}"] = self.predict[0] <= y_ub
        else:
            MILP += self.predict[0] <= y_ub, f"LRq_ub_{prop}"

        if f"LRq_lb_{prop}" in MILP.constraints:
            MILP.constraints[f"LRq_lb_{prop}"] = self.predict[0] >= y_lb
        else:
            MILP += self.predict[0] >= y_lb, f"LRq_lb_{prop}"

        return MILP

    def _predict(self, des):
        y = 0
        if len(des) != self.K:
            print("""
            Error predicting the value.
            """)
            sys.exit()
        # else:
        y = sum([self.coef[i] * des[i] for i in range(0, self.K)])
        y += self.bias

        return y

    def _predict_alpha_star(self, des, descriptors_list, descriptors_q_x, descriptors_q_y, descriptors_q_x_minus, descriptors_q_y_minus):
        y = 0

        z_tmp = dict()
        z_l = list()
        for i in range(1, self.K_l + 1):
            fv_name = descriptors_list[i]
            if i <= len(self.K_l_ind):
                z_tmp[self.K_l_ind[i - 1]] = des[i].value()
            z_l.append(des[i].value())
            # print(f"x_hat[{i}] = {des[desc_dict[fv_name]]}")

        for i in range(1, self.K_q + 1):
            if descriptors_q_x[i] != 0:
                _x = des[descriptors_q_x[i]].value()
            else:
                _x = 0
            if descriptors_q_y[i] != 0:
                _y = des[descriptors_q_y[i]].value()
            else:
                _y = 0
            if descriptors_q_x_minus[i]:
                _x = 1 - _x
            if descriptors_q_y_minus[i]:
                _y = 1 - _y
            z_tmp[self.K_q_ind[i - 1]] = _x * _y
            # z.append(_x * _y)

        z = [z_tmp[j] for j in range(1, self.K + 1)]

        # if len(des) != self.K:
        #     print("""
        #     Error predicting the value.
        #     """)
        #     sys.exit()
        # else:
        y = sum([self.coef[i] * z[i] for i in range(0, self.K)])
        y += self.bias

        return y

########## readline except comments ##########
def readline_except_comment(fp):
    while 1:
        s = fp.readline()
        if s[0]!='#':
            s = s.replace('\n','')
            while len(s)>0 and s[-1] == ' ':
                s = s[:-1]
            if len(s)==0:
                sys.stderr.write("error: illegal line is found.\n")
                exit(1)
            break
    return s

def read_LRq(fp):
    LR = LinearReg_q()
    try:
        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(int, string.split(SPL)))
        ## read dimonsion K
        LR.K = int(arr[0])

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read weight w
        for i in range(len(arr)):
            if abs(arr[i]) > 100000:
                arr[i] = 0.0
        LR.coef = arr

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read bias b
        LR.bias = arr[0]

    except:
        sys.stderr.write("error: failed to read the linear regression file.\n")
        sys.exit(1)
        
    return LR

def create_random_LRq(LR_org, RANDOM_SEED=100):
    LR = LinearReg_q()

    LR.K = LR_org.K
    LR.K_l = LR_org.K_l
    LR.K_q = LR_org.K_q
    LR.K_l_ind = LR_org.K_l_ind[:]
    LR.K_q_ind = LR_org.K_q_ind[:]
    LR.p = LR_org.p
    LR.p_q = LR_org.p_q

    rng = np.random.default_rng(RANDOM_SEED)
    size_I_plus = rng.integers(1, LR.K + 1)
    size_I_minus = LR.K - size_I_plus
    # l_tmp = [i for i in range(LR.K)]
    l_tmp = rng.permutation(LR.K)
    I_plus = l_tmp[:size_I_plus]
    I_minus = l_tmp[size_I_plus:]

    for i in range(LR.K):
        if i in I_plus:
            LR.coef.append(1.0 / LR.K)
        else:
            LR.coef.append(-1.0 / LR.K)

    LR.bias = size_I_minus / LR.K

    return LR

def read_fv_descriptors(training_data_filename):
    """
    Read the textual names for the descriptors
    in the feature vector, as the first line from
    the supplied training data file.
    Return a list of these names as strings
    """
    try:
        data_frame = pd.read_csv(training_data_filename, sep=",")
    except:
        print("""
        Error reading the file {}
        with pandas.
        """.format(training_data_filename))
        sys.exit()
    return list(data_frame.columns)


def read_training_data(training_data_filename):
    """
    Given a set of feature vectors as training data,
    return a matrix that contains one feature vector
    per row
    """
    try:
        data_frame = pd.read_csv(training_data_filename, sep=",")
    except:
        print("""
        Error reading the file {}
        with pandas.
        """.format(training_data_filename))
        sys.exit()
    # print(data_frame.values) # testing

    try:
        table = data_frame.values  # numpy array
        table = table[:, 1:]
    except:
        print("""
        Exception in converting the dataframe
        to a numpy array, file
        {}
        """.format(training_data_filename))
        sys.exit()
    # Success
    return table        
 
def inspection(
    desc_filename, lr_q_filename, x_hat, std_eps,
    K_l, K_q, K_l_ind, K_q_ind,
    descriptors_list, descriptors_l_list, descriptors_q_x, descriptors_q_y, descriptors_q_x_minus, descriptors_q_y_minus
):
    fp = open(lr_q_filename)
    LR = read_LRq(fp)
    data = read_training_data(desc_filename)

    data_frame = pd.read_csv(desc_filename, sep=",")
    columns = data_frame.columns.tolist()

    # data_frame_test = pd.read_csv(desc_test_filename, sep=",")
    # columns_test = data_frame_test.columns.tolist()

    columns_dict = dict()

    desc_dict = dict()

    for fv_name in descriptors_l_list:
        i_dict = -1
        for j in range(len(columns)):
            if fv_name[0] != 'F':
                if columns[j] == fv_name:
                    i_dict = j
                    break
            else:
                ind_i = fv_name.find('_', 3)
                ind_j = columns[j].find('_', 3)
                if columns[j][ind_j:] == fv_name[ind_i:]:
                    i_dict = j
                    break
        desc_dict[fv_name] = i_dict - 1

    y = list()

    for des in data:
        z_tmp = dict()
        z_l = list()
        for i in range(1, K_l + 1):
            fv_name = descriptors_list[i]
            if i <= len(K_l_ind):
                z_tmp[K_l_ind[i - 1]] = des[desc_dict[fv_name]]
            z_l.append(des[desc_dict[fv_name]])
            # print(f"x_hat[{i}] = {des[desc_dict[fv_name]]}")

        for i in range(1, K_q + 1):
            if descriptors_q_x[i] != 0:
                fv_name1 = descriptors_list[descriptors_q_x[i]]
                _x = des[desc_dict[fv_name1]]
            else:
                _x = 0
            if descriptors_q_y[i] != 0:
                fv_name2 = descriptors_list[descriptors_q_y[i]]
                _y = des[desc_dict[fv_name2]]
            else:
                _y = 0
            if descriptors_q_x_minus[i]:
                _x = 1 - _x
            if descriptors_q_y_minus[i]:
                _y = 1 - _y
            z_tmp[K_q_ind[i - 1]] = _x * _y
            # z.append(_x * _y)

        z = [z_tmp[j] for j in range(1, LR.K + 1)]
        y.append(LR._predict(z))

        for i in range(1, K_l + 1):
            if abs(x_hat[i].value() - z_l[i - 1]) > 2*std_eps:
                print(f"i = {i}. [{descriptors_list[i]}] Difference in desc ({z_l[i - 1]}) and x_hat ({x_hat[i].value()}).")

    return y[0], list(data_frame.values[0, 1:])

