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
from typing import List, Dict, Union
import numpy as np
import pandas as pd

# import pulp for MILP modeling
from Module_3.libs import pulp_modified as pulp

# the splitter in the dt file
SPL = ' '
M = 10.0

# LeafNode
class LeafNode:
    def __init__(self, idx, value):
        self.idx = idx
        self.value = value

# InnerNode
class InnerNode:
    def __init__(self, idx, left_ch, right_ch, feature_idx, threshold):
        self.idx = idx
        self.left_ch = left_ch
        self.right_ch = right_ch
        self.feature_idx = feature_idx
        self.threshold = threshold

class DecisionTree:
    def __init__(self, children_left, children_right, feature, threshold, value, tree_idx=0):
        # set up information of the tree
        self.tree_idx: int = tree_idx
        self.children_left: Dict[int, int] = children_left
        self.children_right: Dict[int, int] = children_right
        self.feature: Dict[int, int] = feature
        self.threshold: Dict[int, int] = threshold
        self.value: Dict[int, float] = value

    def get_node(self, idx) -> Union[LeafNode, InnerNode]:
        left = self.children_left.get(idx)
        right = self.children_right.get(idx)
        if left == right:
            # Node idx is leaf node
            return LeafNode(idx, self.value[idx])
        else:
            # Node idx is inner node
            return InnerNode(idx, left, right, self.feature[idx], self.threshold[idx])

    # left_leaf_children returns a list of its left leaf nodes

    def left_leaf_children(self, root_node_idx) -> List[LeafNode]:
        left_leaf_children_list = []
        stack = [self.children_left[root_node_idx]]
        while len(stack) > 0:
            idx = stack.pop()
            left = self.children_left.get(idx)
            right = self.children_right.get(idx)
            if left != right:
                stack.append(left)
                stack.append(right)
            else:
                left_leaf_children_list.append(LeafNode(idx, self.value[idx]))
        return left_leaf_children_list

    # right_leaf_children returns a list of its right leaf nodes
    def right_leaf_children(self, root_node_idx) -> List[LeafNode]:
        right_leaf_children_list = []
        stack = [self.children_right[root_node_idx]]
        while len(stack) > 0:
            idx = stack.pop()
            left = self.children_left.get(idx)
            right = self.children_right.get(idx)
            if left != right:
                stack.append(left)
                stack.append(right)
            else:
                right_leaf_children_list.append(LeafNode(idx, self.value[idx]))
        return right_leaf_children_list

    def predict(self, feature_vector):
        # set root node
        node = self.get_node(0)
        while type(node) != LeafNode:
            if feature_vector[node.feature_idx] <= node.threshold:
                node = self.get_node(node.left_ch)
            else:
                node = self.get_node(node.right_ch)
        return node.value

class RandomForst:
    def __init__(self, dt_list=[]):
        self.dt_list: List[DecisionTree] = dt_list

    def predict(self, feature_vector):
        prediction_values = list()
        for dt in self.dt_list:
            prediction_values.append(dt.predict(feature_vector))
        return sum(prediction_values)/len(prediction_values)

class DTRegInv:
    def __init__(self):
        self.prediction_val = pulp.LpVariable("temp")
        self.delta_list: Dict[int, pulp.LpVariable] = dict()

    def build_var(self, decision_tree_pulp_model: DecisionTree, property_name):
        self.prediction_val = pulp.LpVariable(
            f"DT_{decision_tree_pulp_model.tree_idx}_p_{property_name}", cat=pulp.LpContinuous)
        # the number of decision_tree_pulp_model.value is equal to that of leaves
        self.delta_list = {idx: pulp.LpVariable(
            f"DT_{decision_tree_pulp_model.tree_idx}_delta_{idx}_{property_name}", cat=pulp.LpBinary) for idx in decision_tree_pulp_model.value.keys()}

    # build_constraints builds constraints for a given decision tree
    def build_constraints(self, MILP, feature_vector, decision_tree_pulp_model: DecisionTree, property_name):
        # predictin value constraints (3)
        deltas = list()
        values = list()
        for idx in decision_tree_pulp_model.value.keys():
            deltas.append(self.delta_list[idx])
            values.append(decision_tree_pulp_model.value[idx])
        MILP += pulp.lpDot(deltas,
                           values) == self.prediction_val, f"DT_{decision_tree_pulp_model.tree_idx}_prediction_value_{property_name}"
        # deltas constraint (4)
        MILP += pulp.lpSum(
            self.delta_list.values()) == 1, f"DT_{decision_tree_pulp_model.tree_idx}_one_fv_visits_only_one_leaf_{property_name}"

        # inner_node constraints
        stack = [0]
        while len(stack) > 0:
            node = decision_tree_pulp_model.get_node(stack.pop())
            if type(node) == InnerNode:
                stack.append(node.left_ch)
                stack.append(node.right_ch)

                # inner_node constraints for left children leaves. (5)
                leaves = decision_tree_pulp_model.left_leaf_children(
                    node.idx)
                for leaf_node in leaves:
                    MILP += feature_vector[node.feature_idx] <= node.threshold + M*(
                        1-self.delta_list[leaf_node.idx]), f"DT_{decision_tree_pulp_model.tree_idx}_inner_{node.idx}_leaf_{leaf_node.idx}_{property_name}"
                # inner_node constraints for left children leaves. (6)
                leaves = decision_tree_pulp_model.right_leaf_children(node.idx)
                for leaf_node in leaves:
                    MILP += feature_vector[node.feature_idx] >= node.threshold + np.finfo(np.float64).eps - M*(
                        1-self.delta_list[leaf_node.idx]), f"DT_{decision_tree_pulp_model.tree_idx}_inner_{node.idx}_leaf_{leaf_node.idx}_{property_name}"
            else:
                continue
        return MILP

class RFRegInv:
    def __init__(self, n_estimators, property_name=""):
        self.ensembled_val: pulp.LpVariable = pulp.LpVariable(
            f"y_{property_name}", cat=pulp.LpContinuous)
        self.dt_reg_invs: List[DTRegInv] = [
            DTRegInv() for _ in range(n_estimators)]

    def build_var(self, random_forest: RandomForst, property_name = ""):
        for dt_reg_inv, tree in zip(self.dt_reg_invs, random_forest.dt_list):
            dt_reg_inv.build_var(tree, property_name)

    # build_constraints builds constraints for a given random forest
    def build_constraints(self, MILP, y_lb, y_ub, feature_vector, random_forest: RandomForst, property_name = ""):
        # constraints for lower and upper bounds. (1)
        MILP += self.ensembled_val >= y_lb
        MILP += self.ensembled_val <= y_ub
        # constrain for prediction values and the ensemble value. (2)
        prediction_value_list = list(
            map(lambda x: x.prediction_val, self.dt_reg_invs))
        MILP += self.ensembled_val == pulp.lpSum(prediction_value_list)/len(
            self.dt_reg_invs), "average-of-values-is-equal-to-y_{}".format(property_name)
        for dt_reg_inv, tree in zip(self.dt_reg_invs, random_forest.dt_list):
            # build constraints for each tree. (3), (4), (5), (6)
            MILP = dt_reg_inv.build_constraints(
                MILP, feature_vector, tree, property_name)
        return MILP

def read_rf(save_path):
    with open(save_path, 'r') as f:
        dt_list = []
        n_estimators = int(f.readline())
        for tree_idx in range(n_estimators):
            # set up tree information
            children_left = dict()
            children_right = dict()
            feature = dict()
            threshold = dict()
            value = dict()

            # read inner nodes information
            num_inner_nodes = int(f.readline())
            for _ in range(num_inner_nodes):
                idx, left_ch, right_ch, feature_idx, thresh = f.readline().split()
                idx = int(idx)
                left_ch = int(left_ch)
                right_ch = int(right_ch)
                feature_idx = int(feature_idx)
                thresh = float(thresh)
                children_left[idx] = left_ch
                children_right[idx] = right_ch
                feature[idx] = feature_idx
                threshold[idx] = thresh

            # read leave nodes information
            num_leaf_nodes = int(f.readline())
            for _ in range(num_leaf_nodes):
                idx, val = f.readline().split()
                idx = int(idx)
                val = float(val)
                value[idx] = val
            dt_list.append(DecisionTree(
                children_left, children_right, feature, threshold, value, tree_idx))
    return RandomForst(dt_list)

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
        # sys.exit()
        raise ValueError
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
        # sys.exit()
        raise ValueError
    # Success
    return table      

# def inspection(
#     desc_filename, rf_model, x_hat, std_eps,
#     K_l, K_q, K_l_ind, K_q_ind,
#     descriptors_list, descriptors_l_list, descriptors_q_x, descriptors_q_y, descriptors_q_x_minus, descriptors_q_y_minus,
#     num_fv
# ):
#     data = read_training_data(desc_filename)

#     data_frame = pd.read_csv(desc_filename, sep=",")
#     columns = data_frame.columns.tolist()

#     # data_frame_test = pd.read_csv(desc_test_filename, sep=",")
#     # columns_test = data_frame_test.columns.tolist()

#     columns_dict = dict()

#     desc_dict = dict()

#     for fv_name in descriptors_l_list:
#         i_dict = -1
#         for j in range(len(columns)):
#             if fv_name[0] != 'F':
#                 if columns[j] == fv_name:
#                     i_dict = j
#                     break
#             else:
#                 ind_i = fv_name.find('_', 3)
#                 ind_j = columns[j].find('_', 3)
#                 if columns[j][ind_j:] == fv_name[ind_i:]:
#                     i_dict = j
#                     break
#         desc_dict[fv_name] = i_dict - 1

#     # for i in range(len(columns_test)):
#     #     i_dict = -1
#     #     for j in range(len(columns)):
#     #         if columns_test[i][0] != 'F':
#     #             if columns[j] == columns_test[i]:
#     #                 i_dict = j
#     #                 break
#     #         else:
#     #             ind_i = columns_test[i].find('_', 3)
#     #             ind_j = columns[j].find('_', 3)
#     #             if columns[j][ind_j:] == columns_test[i][ind_i:]:
#     #                 i_dict = j
#     #                 break
#     #     columns_dict[i] = i_dict

#     y = list()

#     for des in data:
#         z_tmp = dict()
#         z_l = list()
#         for i in range(1, K_l + 1):
#             fv_name = descriptors_list[i]
#             if i <= len(K_l_ind):
#                 z_tmp[K_l_ind[i - 1]] = des[desc_dict[fv_name]]
#             z_l.append(des[desc_dict[fv_name]])
#             # print(f"x_hat[{i}] = {des[desc_dict[fv_name]]}")

#         for i in range(1, K_q + 1):
#             if descriptors_q_x[i] != 0:
#                 fv_name1 = descriptors_list[descriptors_q_x[i]]
#                 _x = des[desc_dict[fv_name1]]
#             else:
#                 _x = 0
#             if descriptors_q_y[i] != 0:
#                 fv_name2 = descriptors_list[descriptors_q_y[i]]
#                 _y = des[desc_dict[fv_name2]]
#             else:
#                 _y = 0
#             if descriptors_q_x_minus[i]:
#                 _x = 1 - _x
#             if descriptors_q_y_minus[i]:
#                 _y = 1 - _y
#             z_tmp[K_q_ind[i - 1]] = _x * _y
#             # z.append(_x * _y)

#         z = [z_tmp[j] for j in range(1, num_fv + 1)]
#         y.append(rf_model.predict(z))

#         for i in range(1, K_l + 1):
#             if abs(x_hat[i].value() - z_l[i - 1]) > 2*std_eps:
#                 print(f"i = {i}. [{descriptors_list[i]}] Difference in desc ({z_l[i - 1]}) and x_hat ({x_hat[i].value()}).")

#     return y

def inspection(
    desc_filename, rf_model, x_hat, std_eps,
    num_fv,
    descriptors_list
):
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

    y = list()
    for des in data:
        z_tmp = dict()
        for i in range(1, num_fv):
            fv_name = descriptors_list[i]
            z_tmp[i] = des[desc_dict[fv_name]]
            # print(f"x_hat[{i}] = {des[desc_dict[fv_name]]}")

        z = [z_tmp[j] for j in range(1, num_fv)]
        y.append(rf_model.predict(z))
        
        for i in range(1, num_fv):
            if abs(x_hat[i].value() - z_tmp[i]) > 2*std_eps:
                print(f"i = {i}. [{descriptors_list[i]}] Difference in desc ({z_tmp[i]}) and x_hat ({x_hat[i].value()}).")

    return y

def inspection_for_chi(
    desc_filename, rf_model, x_hat, std_eps,
    K_l, K_q, K_l_ind, K_q_ind,
    descriptors_list, descriptors_l_list, descriptors_q_x, descriptors_q_y, descriptors_q_x_minus, descriptors_q_y_minus,
    fv_other, num_fv
):
    data = read_training_data(desc_filename)

    data_frame = pd.read_csv(desc_filename, sep=",")
    columns = data_frame.columns.tolist()

    # data_frame_test = pd.read_csv(desc_test_filename, sep=",")
    # columns_test = data_frame_test.columns.tolist()

    columns_dict = dict()

    desc_dict = dict()

    for fv_name in descriptors_l_list:
        if fv_name in fv_other:
            desc_dict[fv_name] = fv_other[fv_name]
        else:
            i_dict = -1
            for j in range(len(columns)):
                if fv_name[0] != 'F':
                    if columns[j] + "_ch1" == fv_name:
                        i_dict = j
                        break
                else:
                    ind_i = fv_name.find('_', 3)
                    ind_j = columns[j].find('_', 3)
                    if columns[j][ind_j:] + "_ch1" == fv_name[ind_i:]:
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
            fv_val = des[desc_dict[fv_name]] if fv_name not in fv_other else desc_dict[fv_name]
            if i <= len(K_l_ind):
                z_tmp[K_l_ind[i - 1]] = fv_val
            z_l.append(fv_val)
            # print(f"x_hat[{i}] = {des[desc_dict[fv_name]]}")

        for i in range(1, K_q + 1):
            if descriptors_q_x[i] != 0:
                fv_name1 = descriptors_list[descriptors_q_x[i]]
                fv_val = des[desc_dict[fv_name1]] if fv_name1 not in fv_other else desc_dict[fv_name1]
                _x = fv_val
            else:
                _x = 0
            if descriptors_q_y[i] != 0:
                fv_name2 = descriptors_list[descriptors_q_y[i]]
                fv_val = des[desc_dict[fv_name2]] if fv_name2 not in fv_other else desc_dict[fv_name2]
                _y = fv_val
            else:
                _y = 0
            if descriptors_q_x_minus[i]:
                _x = 1 - _x
            if descriptors_q_y_minus[i]:
                _y = 1 - _y
            z_tmp[K_q_ind[i - 1]] = _x * _y
            # z.append(_x * _y)

        z = [z_tmp[j] for j in range(1, num_fv + 1)]
        y.append(rf_model.predict(z))

        for i in range(1, K_l + 1):
            if abs(x_hat[i].value() - z_l[i - 1]) > 2*std_eps:
                print(f"i = {i}. [{descriptors_list[i]}] Difference in desc ({z_l[i - 1]}) and x_hat ({x_hat[i].value()}).")

    return y
