# -*- coding: utf-8 -*-

"""
This file implements
 - a class to store the architecture of an elastic linear regressor (LR)
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
import subprocess, math

# import pulp for MILP modeling
import pulp

# the splitter in the dt file
SPL = ' '

def scale_dcp(x, x_min, x_max):
    # print(x, x_min, x_max)
    if x_max == x_min:
        ans = 0
    else:
        ans = (x - x_min) / (x_max - x_min)
    return ans

class ElasticLinearReg:
    def __init__(self):
        self.K = 0
        self.coef = dict()
        self.coef_w = dict()
        self.bias = None
        self.y = dict()
        self.yy = None
        self.Delta = dict()   # value of function which is differenced, after being normalized
        self.Delta_Ms = dict()   # value of function which is differenced, after being normalized
        self.x_hat = None
        self.delta = None
        self.delta_atm = None
        self.delta_Ms = None
        self.c_min = dict()
        self.c_max = dict()
        self.mass_ind = None

        self.epsilon = dict()

        self.atm_LB = None
        self.atm_UB = None
        self.Ms_LB = None
        self.Ms_UB = None

        self.M = 1000000

    def y_original(self, y_new):
        ans = 0
        if self.y[1] == self.y[2]:
            ans = y_new / (2 * self.y[2] + self.y[0])
        else:
            if self.y[1] > self.y[2]:
                ans = (- (2 * self.y[2] + self.y[0]) + \
                        math.sqrt((2 * self.y[2] + self.y[0]) ** 2 + 4 * (self.y[1] - self.y[2]) * y_new)) / \
                        (2 * (self.y[1] - self.y[2]))
            else:
                ans = (- (2 * self.y[2] + self.y[0]) - \
                        math.sqrt((2 * self.y[2] + self.y[0]) ** 2 + 4 * (self.y[1] - self.y[2]) * y_new)) / \
                        (2 * (self.y[1] - self.y[2]))   
        return ans             

    def psi_j(self, j, x):
        return self.coef_w[(j, 0)] * x + self.coef_w[(j, 1)] * (x ** 2) + self.coef_w[(j, 2)] * (1 - (x - 1) ** 2)

    def build_Delta(self, mass_ind, min_dcp, max_dcp, fv_lb, fv_ub):
        self.mass_ind = mass_ind
        for i in range(1, self.K + 1):
            if i != mass_ind:
                # self.c_min[i] = min_dcp[i]
                # self.c_max[i] = max_dcp[i]
                self.c_min[i] = fv_lb[i]
                self.c_max[i] = fv_ub[i]

                self.epsilon[i] = 1000000

                self.Delta[(i, self.c_min[i])] = self.psi_j(i, scale_dcp(self.c_min[i], min_dcp[i], max_dcp[i]))
                self.epsilon[i] = min(self.epsilon[i], abs(self.Delta[(i, self.c_min[i])]))
                for j in range(self.c_min[i] + 1, self.c_max[i] + 1):
                    self.Delta[(i, j)] = self.psi_j(i, scale_dcp(j, min_dcp[i], max_dcp[i])) - \
                                         self.psi_j(i, scale_dcp(j - 1, min_dcp[i], max_dcp[i]))
                    self.epsilon[i] = min(self.epsilon[i], abs(self.Delta[(i, j)]))
                self.epsilon[i] /= 100000
                # print(i, self.epsilon[i])

    def calc_Ms_lb_ub(self, set_Lambda, n_LB, n_star, na_LB, na_UB, mass):
        min_mass_1_tmp = 1000000
        min_mass_2_tmp = 1000000
        max_mass_tmp = 0
        for atom in set_Lambda:
            if atom.valence == 1 and atom.mass < min_mass_1_tmp:
                min_mass_1_tmp = atom.mass
            elif atom.valence >= 2 and atom.mass < min_mass_2_tmp:
                min_mass_2_tmp = atom.mass
            if atom.mass > max_mass_tmp:
                max_mass_tmp = atom.mass
        self.Ms_LB = math.floor(min_mass_1_tmp * (n_LB * 0.75) + min_mass_2_tmp * (n_LB * 0.25) + mass["H1"] * na_LB["H1"])
        self.Ms_UB = n_star * max_mass_tmp + mass["H1"] * na_UB["H1"]
        # print(f"{self.Ms_LB}   {self.Ms_UB}")

    def build_Delta_Ms(self, mass_ind, min_dcp, max_dcp, n_LB, n_star, na_LB, na_UB):
        self.c_min[mass_ind] = self.Ms_LB
        self.c_max[mass_ind] = self.Ms_UB
        self.atm_LB = n_LB + na_LB["H1"]
        self.atm_UB = n_star + na_UB["H1"]

        # print(mass_ind, min_dcp[mass_ind], max_dcp[mass_ind])

        psi_j_tmp = dict()
        for s in range(self.Ms_LB, self.Ms_UB + 1):
            for i in range(self.atm_LB, self.atm_UB + 1):
                psi_j_tmp[(s, i)] = int(self.psi_j(mass_ind, scale_dcp(s / i, min_dcp[mass_ind], max_dcp[mass_ind])) * 1e10) / 1e10

        for s in range(self.Ms_LB, self.Ms_UB + 1):
            for i in range(self.atm_LB, self.atm_UB + 1):
                if s == self.Ms_LB:
                    self.Delta_Ms[(s, i)] = psi_j_tmp[(s, i)]
                else:
                    self.Delta_Ms[(s, i)] = psi_j_tmp[(s, i)] - psi_j_tmp[(s - 1, i)]

        # print(self.Delta_Ms)
        # print(min(self.Delta_Ms))
        min_abs = 1000000
        for x_tmp in self.Delta_Ms:
            if abs(self.Delta_Ms[x_tmp]) < min_abs:
                min_abs = abs(self.Delta_Ms[x_tmp])
        self.epsilon[mass_ind] = min_abs / 100000
        self.M = self.Ms_UB / self.atm_LB * 10

    def build_var(self, mass_ind):
        self.mass_ind = mass_ind
        self.x_hat = {i: pulp.LpVariable(f"x_hat({i})", cat=pulp.LpContinuous) 
                        for i in range(1, self.K + 1)}
        self.delta = {(j, s): pulp.LpVariable(f"ELR_delta({j},{s})", 0, 1, cat=pulp.LpBinary)
                        for j in range(1, self.K + 1) if j != mass_ind for s in range(self.c_min[j], self.c_max[j] + 1)}
        self.delta_atm = {i: pulp.LpVariable(f"ELR_delta_atom({i})", 0, 1, cat=pulp.LpBinary)
                            for i in range(self.atm_LB, self.atm_UB + 1)}
        self.delta_Ms = {s: pulp.LpVariable(f"ELR_delta_Ms({s})", 0, 1, cat=pulp.LpBinary)
                            for s in range(self.Ms_LB, self.Ms_UB + 1)}

        self.yy = pulp.LpVariable(f"y_ELR", cat=pulp.LpContinuous)

        return self.yy

    def build_constraints(self, MILP, descriptors, mass, y_LB, y_UB, n_G, na, na_ex, MASS):

        yy_LB = self.y[0] * y_LB + self.y[1] * (y_LB ** 2) + self.y[2] * (1 - ((y_LB - 1) ** 2))
        yy_UB = self.y[0] * y_UB + self.y[1] * (y_UB ** 2) + self.y[2] * (1 - ((y_UB - 1) ** 2))

        MILP += pulp.lpSum([self.coef[j] * self.x_hat[j]] for j in range(1, self.K + 1)) + self.bias == self.yy

        MILP += pulp.lpSum([self.coef[j] * self.x_hat[j]] for j in range(1, self.K + 1)) + self.bias >= yy_LB, f"milp-ELR-(98)-1"
        MILP += pulp.lpSum([self.coef[j] * self.x_hat[j]] for j in range(1, self.K + 1)) + self.bias <= yy_UB, f"milp-ELR-(98)-2"

        for j in range(1, self.K + 1):
            if j != self.mass_ind:
                MILP += pulp.lpSum(self.delta[(j, s)] for s in range(self.c_min[j], self.c_max[j] + 1)) + self.c_min[j] - 1 == descriptors[j], f"milp-ELR-(99)-1-{j}"
                for s in range(self.c_min[j], self.c_max[j]):
                    MILP += self.delta[(j, s)] >= self.delta[(j, s + 1)], f"milp-ELR-(99)-2-{j}-{s}"
                MILP += self.x_hat[j] >= pulp.lpSum(self.Delta[(j, s)] * self.delta[(j, s)]
                                for s in range(self.c_min[j], self.c_max[j] + 1)) - self.epsilon[j], f"milp-ELR-(99)-3-{j}-1"
                MILP += self.x_hat[j] <= pulp.lpSum(self.Delta[(j, s)] * self.delta[(j, s)]
                                for s in range(self.c_min[j], self.c_max[j] + 1)) + self.epsilon[j], f"milp-ELR-(99)-3-{j}-2"

        MILP += pulp.lpSum(self.delta_atm[i] for i in range(self.atm_LB, self.atm_UB + 1)) == 1, f"milp-ELR-(100)-1"
        MILP += pulp.lpSum(i * self.delta_atm[i] for i in range(self.atm_LB, self.atm_UB + 1)) == n_G + na_ex["H1"], f"milp-ELR-(100)-2"
        MILP += MASS == \
                pulp.lpSum(self.delta_Ms[s] for s in range(self.Ms_LB, self.Ms_UB + 1)) + self.Ms_LB - 1, f"milp-ELR-(100)-3"
        for s in range(self.Ms_LB, self.Ms_UB):
            MILP += self.delta_Ms[s] >= self.delta_Ms[s + 1], f"milp-ELR-(101)-1-{s}"
        for i in range(self.atm_LB, self.atm_UB + 1):
            MILP += self.x_hat[self.mass_ind] >= pulp.lpSum(self.Delta_Ms[(s, i)] * self.delta_Ms[s] 
                        for s in range(self.Ms_LB, self.Ms_UB + 1)) - self.M * (1 - self.delta_atm[i]) - self.epsilon[self.mass_ind], f"milp-ELR-(101)-2-{i}-1"
            MILP += self.x_hat[self.mass_ind] <= pulp.lpSum(self.Delta_Ms[(s, i)] * self.delta_Ms[s] 
                        for s in range(self.Ms_LB, self.Ms_UB + 1)) + self.M * (1 - self.delta_atm[i]) + self.epsilon[self.mass_ind], f"milp-ELR-(101)-2-{i}-2"

        return MILP

    def _predict(self, des, columns_dict):
        y = 0
        # if len(des) != self.K:
        #     print("""
        #     Error predicting the value.
        #     """)
        #     sys.exit()
        # else:

        des_new = dict()
        for i in range(1, self.K + 1):
            des_new[i] = self.psi_j(i, des[columns_dict[i] - 1]) 

        y = sum([self.coef[i] * des_new[i] for i in range(1, self.K + 1)])
        y += self.bias

        return y, des_new


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

def read_ELR(fp):
    ELR = ElasticLinearReg()
    try:
        ############
        # not totally determined
        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(int, string.split(SPL)))
        ## read dimonsion K
        ELR.K = int(arr[0])

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read weight w
        for i in range(1, ELR.K + 1):
            ELR.coef[i] = arr[i - 1]

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read weight w_0
        for i in range(1, ELR.K + 1):
            ELR.coef_w[(i, 0)] = arr[i - 1]

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read weight w_1
        for i in range(1, ELR.K + 1):
            ELR.coef_w[(i, 1)] = arr[i - 1]

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read weight w_2
        for i in range(1, ELR.K + 1):
            ELR.coef_w[(i, 2)] = arr[i - 1]

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read bias b
        ELR.bias = arr[0]

        string = readline_except_comment(fp)
        # print(string)
        arr = list(map(float, string.split(SPL)))
        ## read weight for y
        ELR.y[0] = arr[0]
        ELR.y[1] = arr[1]
        ELR.y[2] = arr[2]
        ##############

    except:
        sys.stderr.write("error: failed to read the linear regression file.\n")
        sys.exit(1)
        
    return ELR    


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
 
def inspection(desc_filename, desc_test_filename, elr_filename, x_hat, std_eps):
    fp = open(elr_filename)
    ELR = read_ELR(fp)
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
        y_ins, des_new = ELR._predict(des, columns_dict)
        y_ins = ELR.y_original(y_ins)
        y.append(y_ins)

        for i in range(1, len(des) + 1):
            if abs(x_hat[i].value() - des_new[i]) > 2*std_eps:
                print(f"i = {i}. [{columns[i]}] Difference in desc ({des_new[i]}) and x_hat ({x_hat[i].value()}).")

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

