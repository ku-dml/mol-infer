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
from ... import pulp_modified as pulp
# the splitter in the dt file
SPL = ' '

class LinearReg:
    def __init__(self):
        self.K = 0
        self.coef = list()
        self.bias = 0
        self.weight_var = dict()
        self.predict = dict()

    def build_weight_var(self, I_integer, I_nonneg, prop):
        self.weight_var = {i: pulp.LpVariable(f"LR_x_{i}_{prop}", cat=pulp.LpContinuous) for i in range(1, self.K + 1)}
        self.predict = {0: pulp.LpVariable(f"LR_y_{prop}", cat=pulp.LpContinuous)}

        # for i in range(1, self.K + 1):
        #     if i in I_integer:
        #         self.weight_var[i].cat = pulp.LpInteger
        #     if i in I_nonneg:
        #         self.weight_var[i].lowbound = 0

        return self.weight_var, self.predict[0]

    def build_constraints(self, MILP, y_lb, y_ub, prop):
        MILP += self.predict[0] == pulp.lpSum([self.coef[i - 1] * self.weight_var[i] for i in range(1, self.K + 1)]) + self.bias, "LR_lr_{}".format(prop)
        MILP += self.predict[0] <= y_ub, "LR_ub_{}".format(prop)
        MILP += self.predict[0] >= y_lb, "LR_lb_{}".format(prop)

        return MILP

    # def _predict(self, des, columns_dict):
    #     y = 0
    #     # if len(des) != self.K:
    #     #     print("""
    #     #     Error predicting the value.
    #     #     """)
    #     #     sys.exit()
    #     # else:
    #     y = sum([self.coef[i] * des[columns_dict[i + 1] - 1] for i in range(0, self.K)])
    #     y +=  self.bias

    #     return y
    
    def _predict(self, des):
        y = 0
        if len(des) != self.K:
            print("""
            Error predicting the value.
            """)
            sys.exit()

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

def read_LR(fp):
    LR = LinearReg()
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
        LR.bias = int(arr[0])

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
 
def inspection(desc_filename, LR, x_hat, std_eps,
    num_fv,
    descriptors_list):
    # fp = open(lr_filename)
    # LR = read_LR(fp)
    # data = read_training_data(desc_filename)

    # data_frame = pd.read_csv(desc_filename, sep=",")
    # columns = data_frame.columns.tolist()

    # data_frame_test = pd.read_csv(desc_test_filename, sep=",")
    # columns_test = data_frame_test.columns.tolist()

    # columns_dict = dict()

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

    data = read_training_data(desc_filename)

    data_frame = pd.read_csv(desc_filename, sep=",")
    columns = data_frame.columns.tolist()

    # data_frame_test = pd.read_csv(desc_test_filename, sep=",")
    # columns_test = data_frame_test.columns.tolist()

    columns_dict = dict()

    desc_dict = dict()

    for fv_name in descriptors_list:
        if fv_name == "CID":
            continue
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

    # y = list()

    # for des in data:
        
        
    #     print(des.shape)
    #     y.append(LR._predict(des, columns_dict))

    #     for i in range(1, len(des)):
    #         if abs(x_hat[i].value() - des[columns_dict[i] - 1]) > 2*std_eps:
    #             print(f"i = {i}. [{columns[i]}] Difference in desc ({des[columns_dict[i] - 1]}) and x_hat ({x_hat[i].value()}).")


    y = list()
    for des in data:
        z_tmp = dict()
        for i in range(1, num_fv):
            fv_name = descriptors_list[i]
            z_tmp[i] = des[desc_dict[fv_name]]
            # print(f"x_hat[{i}] = {des[desc_dict[fv_name]]}")

        z = [z_tmp[j] for j in range(1, num_fv)]
        y.append(LR._predict(z))
        
        for i in range(1, num_fv):
            if abs(x_hat[i].value() - z_tmp[i]) > 2*std_eps:
                print(f"i = {i}. [{descriptors_list[i]}] Difference in desc ({z_tmp[i]}) and x_hat ({x_hat[i].value()}).")


    return y

def inspection_no_x(desc_filename, desc_test_filename, lr_filename):
    fp = open(lr_filename)
    LR = read_LR(fp)
    data = read_training_data(desc_filename)

    data_frame = pd.read_csv(desc_filename, sep=",")
    columns = data_frame.columns.tolist()

    data_frame_test = pd.read_csv(desc_test_filename, sep=",")
    columns_test = data_frame_test.columns.tolist()

    columns_dict = dict()

    for i in range(len(columns_test)):
        i_dict = -1
        for j in range(len(columns)):
            if columns_test[i][0] != 'F':
                if columns[j] == columns_test[i]:
                    i_dict = j
                    break
            else:
                ind_i = columns_test[i].find('_', 3)
                ind_j = columns[j].find('_', 3)
                if columns[j][ind_j:] == columns_test[i][ind_i:]:
                    i_dict = j
                    break
        columns_dict[i] = i_dict

    y = list()

    for des in data:
        y.append(LR._predict(des, columns_dict))

        # for i in range(1, len(des) + 1):
        #     if abs(x_hat[i].value() - des[columns_dict[i] - 1]) > 2*std_eps:
        #         print(f"i = {i}. [{columns[i]}] Difference in desc ({des[columns_dict[i] - 1]}) and x_hat ({x_hat[i].value()}).")


    return y[0]

def main(argv):
    sdf_all_filename = argv[1]
    sdf_filename = argv[2]
    desc_filename = argv[3]
    lr_filename = argv[4]
    fv_gen_name = "./2LMM_v018/FV_2LMM_V018"
    subprocess.run([fv_gen_name, sdf_all_filename, f"test", sdf_filename, "test_tmp"],
                       stdout=subprocess.DEVNULL)
    y = inspection_no_x("test_tmp_desc_norm.csv", desc_filename, lr_filename)

    print(y)
    print(y[0] * (117 + 165) - 117)

if __name__=="__main__":
    main((0, "OptR.sdf", "OptR_c_70_71.sdf", "OptR_desc_norm.csv", "OptR_linreg.txt"))

