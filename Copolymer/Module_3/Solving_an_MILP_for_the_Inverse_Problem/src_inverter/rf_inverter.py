from __future__ import annotations
from .inverter_base import *
# import pulp_modified as pulp
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, ClassVar, Union
import pandas as pd
import numpy as np
from src.utils import _read_clean_lines, _pop_int, _pop_ints, _pop_floats, _pop_float, _pop_strs, run_with_retry

@dataclass(frozen=False, kw_only=True)
class RFInverter(InverterBase):
	### An inverter class for random forest
	name: str
	_type: ClassVar[Literal["INV_RF"]] = "INV_RF"

	# LeafNode
	@dataclass(frozen=True)
	class LeafNode:
		idx: int
		value: float

	# InnerNode
	@dataclass(frozen=True)
	class InnerNode:
		idx: int
		left_ch: int
		right_ch: int
		feature_idx: int
		threshold: float

	class DecisionTree:
		def __init__(self, children_left, children_right, feature, threshold, value, tree_idx=0):
			# set up information of the tree
			self.tree_idx: int = tree_idx
			self.children_left: Dict[int, int] = children_left
			self.children_right: Dict[int, int] = children_right
			self.feature: Dict[int, int] = feature
			self.threshold: Dict[int, int] = threshold
			self.value: Dict[int, float] = value

		def get_node(self, idx) -> Union[RFInverter.LeafNode, RFInverter.InnerNode]:
			left = self.children_left.get(idx)
			right = self.children_right.get(idx)
			if left == right:
				# Node idx is leaf node
				return RFInverter.LeafNode(idx, self.value[idx])
			else:
				# Node idx is inner node
				return RFInverter.InnerNode(idx, left, right, self.feature[idx], self.threshold[idx])

		# left_leaf_children returns a list of its left leaf nodes
		def left_leaf_children(self, root_node_idx) -> List[RFInverter.LeafNode]:
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
					left_leaf_children_list.append(RFInverter.LeafNode(idx, self.value[idx]))
			return left_leaf_children_list

		# right_leaf_children returns a list of its right leaf nodes
		def right_leaf_children(self, root_node_idx) -> List[RFInverter.LeafNode]:
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
					right_leaf_children_list.append(RFInverter.LeafNode(idx, self.value[idx]))
			return right_leaf_children_list

		def predict(self, feature_vector):
			# set root node
			node = self.get_node(0)
			while type(node) != RFInverter.LeafNode:
				if feature_vector[node.feature_idx] <= node.threshold:
					node = self.get_node(node.left_ch)
				else:
					node = self.get_node(node.right_ch)
			return node.value

	class RandomForest:
		def __init__(self, dt_list=[]):
			self.dt_list: List[RFInverter.DecisionTree] = dt_list

		def predict(self, feature_vector):
			prediction_values = list()
			for dt in self.dt_list:
				prediction_values.append(dt.predict(feature_vector))
			return sum(prediction_values) / len(prediction_values)

	def __post_init__(self):
		self.output = dict()
		self.index_set = dict() # store the relation between the indices of fringe trees in fv files and in seed graph files, index_set[index_in_fv] = index_in_seed_graph
		self.x_hat = dict() # store the normalized feature vector value
		self.std_eps = 1e-5
		self.descriptors = dict()
		self.forbidden_node = dict()

		self.prediction_val = dict()
		self.delta_list = dict()

		self.M = 10

	def read_fringe_tree(self, fv_fringe_tree_filename):
		self.fv_fringe_tree = dict()
		with open(fv_fringe_tree_filename, 'r') as f:
			lines = f.readlines()
			if lines[0][0] == 'F':
				lines = lines[1:]
			for line in lines:
				if len(line.split(",")) < 4:
					continue
				ind = int(line.split(",")[0])
				str1 = line.split(",")[1]
				str2 = line.split(",")[2]
				str3 = line.split(",")[3].replace('\n', '')
				self.fv_fringe_tree[ind] = (str1, str2, str3)

	def prepare_fringe_tree_index_set(self, asgml_ind, strF, init=True):
		if asgml_ind not in self.index_set or init:
			self.index_set[asgml_ind] = dict()
		for (key, f) in strF.items():
			for (ind, fv_fringe_tree) in self.fv_fringe_tree.items():
				if fv_fringe_tree == f:
					self.index_set[asgml_ind][ind] = key
					break

	def set_up_inverter(self, invl):
		self.sdf_filename = invl.sdf_filename
		self.fv_gen_path = invl.fv_gen_path
		self._use_quadratic = invl._use_quadratic

		rf_filename = invl.rf_filename

		try:
			Q = _read_clean_lines(rf_filename)
			dt_list = []
			n_estimators = _pop_int(Q)
			for tree_idx in range(n_estimators):
				# set up tree information
				children_left = dict()
				children_right = dict()
				feature = dict()
				threshold = dict()
				value = dict()

				# read inner nodes information
				num_inner_nodes = _pop_int(Q)
				for _ in range(num_inner_nodes):
					idx, left_ch, right_ch, feature_idx, thresh = _pop_strs(Q)
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
				num_leaf_nodes = _pop_int(Q)
				for _ in range(num_leaf_nodes):
					idx, val = _pop_strs(Q)
					idx = int(idx)
					val = float(val)
					value[idx] = val
				dt_list.append(RFInverter.DecisionTree(children_left, children_right, feature, threshold, value, tree_idx))

			self.rf = RFInverter.RandomForest(dt_list)
		except:
			raise ValueError(f"Error: Failed to read the inverter {self.name}!")
	
		self.read_fringe_tree(invl.fringe_filename)
		self.prepare_fv_list(invl.desc_norm_selected_filename)
		self.prepare_desc_max_min(invl.desc_filename)
		self.prepare_value_max_min(invl.desc_norm_selected_filename, invl.values_filename)

		self.num_K = len(self.fv_list)
		
	def set_up_variables(self, asgml):
		self.prediction_val[asgml.ind] = dict()
		self.delta_list[asgml.ind] = dict()
		for dt in self.rf.dt_list:
			self.prediction_val[asgml.ind][dt.tree_idx] = pulp.LpVariable(f"RF_DT_{dt.tree_idx}_p_ASGM{asgml.ind}", cat=pulp.LpContinuous)
			self.delta_list[asgml.ind][dt.tree_idx] = {idx: pulp.LpVariable(f"RF_DT_{dt.tree_idx}_delta_{idx}_ASGM{asgml.ind}", 0, 1, cat=pulp.LpBinary) 
														for idx in dt.value.keys()}

		### The output variable uses self.output as the variable name
		self.output[asgml.ind] = pulp.LpVariable(f"RF_y_ASGM{asgml.ind}", cat=pulp.LpContinuous)

	def set_up_x_hat_variables(self, asgml):
		### the idx in scikit-learn for DT starts from 0, different from the default setting of 1
		self.x_hat[asgml.ind] = {i: pulp.LpVariable(f"x_hat({i})_ASGM{asgml.ind}") for i in range(1, self.num_K + 1)}
				
	def add_constraints_C1(self, MILP, asgml):
		prediction_value_list = [self.prediction_val[asgml.ind][dt.tree_idx] for dt in self.rf.dt_list]
		MILP += self.output[asgml.ind] == pulp.lpSum(prediction_value_list) / len(self.rf.dt_list), f"RF_y_ASGM{asgml.ind}"

		for dt in self.rf.dt_list:
			deltas = list()
			values = list()
			for idx in dt.value.keys():
				deltas.append(self.delta_list[asgml.ind][dt.tree_idx][idx])
				values.append(dt.value[idx])
			MILP += pulp.lpDot(deltas, values) == self.prediction_val[asgml.ind][dt.tree_idx], f"RF_DT_{dt.tree_idx}_prediction_value_ASGM{asgml.ind}"

			MILP += pulp.lpSum(self.delta_list[asgml.ind][dt.tree_idx].values()) == 1, f"RF_DT_{dt.tree_idx}_deltas_ASGM{asgml.ind}"

			stack = [0]
			while len(stack) > 0:
				node = dt.get_node(stack.pop())
				if type(node) == RFInverter.InnerNode:
					stack.append(node.left_ch)
					stack.append(node.right_ch)

					leaves = dt.left_leaf_children(node.idx)
					for leaf_node in leaves:
						### the +1 here deals with the difference of idx
						MILP += self.x_hat[asgml.ind][node.feature_idx + 1] <= node.threshold + self.M * (1 - self.delta_list[asgml.ind][dt.tree_idx][leaf_node.idx]), \
								f"RF_DT_{dt.tree_idx}_inner_{node.idx}_leaf_{leaf_node.idx}_ASGM{asgml.ind}"
					leaves = dt.right_leaf_children(node.idx)
					for leaf_node in leaves:
						### the +1 here deals with the difference of idx
						MILP += self.x_hat[asgml.ind][node.feature_idx + 1] >= node.threshold + self.std_eps - self.M * (1 - self.delta_list[asgml.ind][dt.tree_idx][leaf_node.idx]), \
								f"RF_DT_{dt.tree_idx}_inner_{node.idx}_leaf_{leaf_node.idx}_ASGM{asgml.ind}"
				else:
					continue

		return MILP

	def _inspection(self, fv_test):
		return self.rf.predict(fv_test)

	def inspection(self, asgml, fv_test, std_eps=None):
		if std_eps is None:
			std_eps = self.std_eps
		insp_value = self._inspection(fv_test)

		found_error = False
		for i in range(1, self.num_K + 1):
			if abs(self.x_hat[asgml.ind][i].value() - fv_test[i - 1]) > 2 * std_eps:
				print(f"Assignment: {asgml.ind}. [i = {i}, {self.fv_list[i - 1]}] Difference between fv_from_MILP ({self.x_hat[asgml.ind][i].value():.6f}) and fv_from_sdf ({fv_test[i - 1]:.6f}).")
				found_error = True
		if not found_error:
			print(f"Assignment: {asgml.ind}. No difference between feature vectors!")

		return insp_value

	