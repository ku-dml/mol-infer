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
import pulp

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
        self.p_q = self.p + 2
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

        self.M = 5

    def build_var(self):
        # self.x = {i: pulp.LpVariable(f"LRq_x({i})", 0, 1, cat=pulp.LpContinuous) for i in range(1, self.K_q + 1)}
        # self.y = {i: pulp.LpVariable(f"LRq_y({i})", 0, 1, cat=pulp.LpContinuous) for i in range(1, self.K_q + 1)}
        # self.x_q = {(i, q): pulp.LpVariable(f"LRq_xq({i},{q})", 0, 1, cat=pulp.LpBinary) 
        #                 for i in range(1, self.K_q + 1) for q in range(self.p + 1)}
        # self.z = {i: pulp.LpVariable(f"LRq_z({i})", 0, 1, cat=pulp.LpContinuous) for i in range(1, self.K_q + 1)}
        # self.z_q = {(i, q): pulp.LpVariable(f"LRq_zq({i},{q})", 0, 1, cat=pulp.LpContinuous) 
        #                 for i in range(1, self.K_q + 1) for q in range(self.p + 1)}

        # self.weight_var = {i: pulp.LpVariable(f"LR_x_{i}", cat=pulp.LpContinuous) for i in range(1, self.K_l + 1)}
        # self.predict = {0: pulp.LpVariable(f"LR_y", cat=pulp.LpContinuous)}

        self.x = {i: pulp.LpVariable(f"LRq_x({i})", 0, cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.y = {i: pulp.LpVariable(f"LRq_y({i})", 0, cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.x_q = {(i, q): pulp.LpVariable(f"LRq_xq({i},{q})", 0, 1, cat=pulp.LpBinary) 
                        for i in self.K_q_ind for q in range(self.p_q + 1)}
        self.z_tmp = {i: pulp.LpVariable(f"LRq_z_tmp({i})", 0, cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.z = {i: pulp.LpVariable(f"LRq_z({i})", cat=pulp.LpContinuous) for i in self.K_q_ind}
        self.z_q = {(i, q): pulp.LpVariable(f"LRq_zq({i},{q})", 0, cat=pulp.LpContinuous) 
                        for i in self.K_q_ind for q in range(self.p_q + 1)}

        self.weight_var = {i: pulp.LpVariable(f"LR_x_{i}", cat=pulp.LpContinuous) for i in range(1, self.K_l + 1)}
        self.weight_var[0] = 0
        self.predict = {0: pulp.LpVariable(f"LR_y", cat=pulp.LpContinuous)}

        return self.weight_var, self.x, self.y, self.predict[0]

    def build_constraints(self, MILP, y_lb, y_ub, descriptors_q_x, descriptors_q_y, descriptors_q_x_minus, descriptors_q_y_minus):
        for i in range(1, self.K_q + 1):
            q_ind = self.K_q_ind[i - 1]
            MILP += (2**(self.p + 1) - 1) * self.x[q_ind] >= \
                    pulp.lpSum([2**q * self.x_q[(q_ind, q)] for q in range(self.p_q + 1)]) - 1, f"LRq-{i}-1-1"
            MILP += (2**(self.p + 1) - 1) * self.x[q_ind] <= \
                    pulp.lpSum([2**q * self.x_q[(q_ind, q)] for q in range(self.p_q + 1)]), f"LRq-{i}-1-2"
            for q in range(self.p_q + 1):
                MILP += self.z_q[(q_ind, q)] <= self.M * self.x_q[(q_ind, q)], f"LRq-{i}-2-1-{q}"
                MILP += self.z_q[(q_ind, q)] >= self.y[q_ind] - self.M * (1 - self.x_q[(q_ind, q)]), f"LRq-{i}-2-2-{q}-1"
                MILP += self.z_q[(q_ind, q)] <= self.y[q_ind] + self.M * (1 - self.x_q[(q_ind, q)]), f"LRq-{i}-2-2-{q}-2"
            MILP += self.z_tmp[q_ind] * (2**(self.p + 1) - 1) == \
                    pulp.lpSum([2**q * self.z_q[(q_ind, q)] for q in range(self.p_q + 1)]), f"LRq-{i}-3"

            if descriptors_q_x_minus[i] and descriptors_q_y_minus[i]:
                MILP += self.z[q_ind] == 1 - self.weight_var[descriptors_q_x[i]] - self.weight_var[descriptors_q_y[i]] + self.z_tmp[q_ind], f"LRq-{i}-4"
            elif descriptors_q_x_minus[i]:
                MILP += self.z[q_ind] == self.weight_var[descriptors_q_y[i]] - self.z_tmp[q_ind], f"LRq-{i}-4"
            elif descriptors_q_y_minus[i]:
                MILP += self.z[q_ind] == self.weight_var[descriptors_q_x[i]] - self.z_tmp[q_ind], f"LRq-{i}-4"
            else:
                MILP += self.z[q_ind] == self.z_tmp[q_ind], f"LRq-{i}-4"

        MILP += self.predict[0] == \
                    pulp.lpSum([self.coef[i - 1] * self.weight_var[j + 1] for (i, j) in enumerate(self.K_l_ind)]) + \
                    pulp.lpSum([self.coef[i - 1] * self.z[i] for i in self.K_q_ind]) + \
                    self.bias, "LRq_lr"
        MILP += self.predict[0] <= y_ub, "LRq_ub"
        MILP += self.predict[0] >= y_lb, "LRq_lb"

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

    # for i in range(len(columns_test)):
    #     i_dict = -1
    #     for j in range(len(columns)):
    #         if columns_test[i][0] != 'F':
    #             if columns[j] == columns_test[i]:
    #                 i_dict = j
    #                 break
    #         else:
    #             ind_i = columns_test[i].find('_', 3)
    #             ind_j = columns[j].find('_', 3)
    #             if columns[j][ind_j:] == columns_test[i][ind_i:]:
    #                 i_dict = j
    #                 break
    #     columns_dict[i] = i_dict

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

    return y

# def inspection_no_x(desc_filename, desc_test_filename, lr_filename):
#     fp = open(lr_filename)
#     LR = read_LR(fp)
#     data = read_training_data(desc_filename)

#     data_frame = pd.read_csv(desc_filename, sep=",")
#     columns = data_frame.columns.tolist()

#     data_frame_test = pd.read_csv(desc_test_filename, sep=",")
#     columns_test = data_frame_test.columns.tolist()

#     columns_dict = dict()

#     for i in range(len(columns_test)):
#         i_dict = -1
#         for j in range(len(columns)):
#             if columns_test[i][0] != 'F':
#                 if columns[j] == columns_test[i]:
#                     i_dict = j
#                     break
#             else:
#                 ind_i = columns_test[i].find('_', 3)
#                 ind_j = columns[j].find('_', 3)
#                 if columns[j][ind_j:] == columns_test[i][ind_i:]:
#                     i_dict = j
#                     break
#         columns_dict[i] = i_dict

#     y = list()

#     for des in data:
#         y.append(LR._predict(des, columns_dict))

#         # for i in range(1, len(des) + 1):
#         #     if abs(x_hat[i].value() - des[columns_dict[i] - 1]) > 2*std_eps:
#         #         print(f"i = {i}. [{columns[i]}] Difference in desc ({des[columns_dict[i] - 1]}) and x_hat ({x_hat[i].value()}).")


#     return y

# def main(argv):
#     sdf_all_filename = argv[1]
#     sdf_filename = argv[2]
#     desc_filename = argv[3]
#     lr_filename = argv[4]
#     fv_gen_name = "./2LMM_v018/FV_2LMM_V018"
#     subprocess.run([fv_gen_name, sdf_all_filename, f"test", sdf_filename, "test_tmp"],
#                        stdout=subprocess.DEVNULL)
#     y = inspection_no_x("test_tmp_desc_norm.csv", desc_filename, lr_filename)

#     print(y)
#     print(y[0] * (117 + 165) - 117)

# if __name__=="__main__":
#     main((0, "OptR.sdf", "OptR_c_70_71.sdf", "OptR_desc_norm.csv", "OptR_linreg.txt"))

