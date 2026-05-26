from .inverter_base import *
# import pulp_modified as pulp
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, ClassVar
import pandas as pd
import numpy as np
from src.utils import _read_clean_lines, _pop_int, _pop_ints, _pop_floats, _pop_float, run_with_retry

@dataclass(frozen=False, kw_only=True)
class ANNInverter(InverterBase):
	### An inverter class for ANN
	name: str
	_type: ClassVar[Literal["INV_ANN"]] = "INV_ANN"

	def __post_init__(self, name="ann"):
		self.output = dict()
		self.index_set = dict() # store the relation between the indices of fringe trees in fv files and in seed graph files, index_set[index_in_fv] = index_in_seed_graph
		self.x_hat = dict() # store the normalized feature vector value
		self.std_eps = 1e-5

		self.x = dict()
		self.y = dict()
		self.z = dict()

	@property
	def non_input_nodes(self):
		for layer in self.layers[1:]:
			for node in layer:
				yield node

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

	def _initialize(self, weight_tensor, bias_matrix):
		self.nodes = list()
		self.layers = list()
		self.weights = dict()
		self.biases = dict()
		self.in_neighbors = dict()

		for l, bias_layer in enumerate(bias_matrix):
			layer_nodes = list()
			for k, bias in enumerate(bias_layer):
				node = (l + 1, k + 1)
				self.nodes.append(node)
				self.biases[node] = bias
				layer_nodes.append(node)
			self.layers.append(layer_nodes)
		for l, layer in enumerate(self.layers):
			try:
				for node in self.layers[l+1]:
					self.in_neighbors[node] = layer
			except IndexError:
				pass
			except ValueError:
				raise ValueError("Received Value Error in distributing ANN nodes by layers")

		for l, weight_matrix in enumerate(weight_tensor):
			for i, weight_row in enumerate(weight_matrix):
				for j, weight in enumerate(weight_row):
					u = (l + 1, i + 1)
					v = (l + 2, j + 1)
					self.weights[(u, v)] = weight

	def _initialize_constants(self):
		# td = np.array(training_data)

		# if td.shape[1] != len(ann.input_layer):
		#     print("###############")
		#     print("""
		#     Vector length mismatch between training data
		#     and ANN input layer!
		#     """)
		#     sys.exit()

		# prepare dictionaries to store values
		a_low = dict()
		a_high = dict()
		a = dict()
		b_low = dict()
		b_high = dict()
		b_hat = dict()
		b = dict()
		c = dict()
		z13 = dict()
		# first, calculate b_low, b_high values for the
		# nodes of the input layer
		# for node, column in zip(ann.input_layer, td.T):
		#     b_low[node] = min(column)
		#     b_high[node] = 1000
		#     b[node] = (b_low[node], 0, b_high[node])

		for node in self.layers[0]:
			b_low[node] = -5
			b_high[node] = 5
			b[node] = (b_low[node], 0, b_high[node])

		# calculate the ranges for the remaining layers
		for v in self.non_input_nodes:
			a_low[v], a_high[v] = 0, 0
			for u in self.in_neighbors[v]:
				w_uv = self.weights[(u, v)]
				if w_uv > 0:
					a_low[v] += w_uv * b_low[u]
					a_high[v] += w_uv * b_high[u]
				else:
					a_high[v] += w_uv * b_low[u]
					a_low[v] += w_uv * b_high[u]
			if self.biases[v] >= 0:
				a_high[v] += self.biases[v]
			else:
				a_low[v] += self.biases[v]
			# a_high[v] += self.biases[v]
			# a_low[v] += self.biases[v]

			a[v] = (a_low[v], 0, a_high[v])

			b_low[v] = 0

			if v in self.layers[-1]:
				b_low[v] = -1000
			
			if a_high[v] > 0:
				b_high[v] = a_high[v]
			else:
				b_high[v] = 0
			
			b[v] = (b_low[v], 0, b_high[v])
			b_hat[v] = 2 * a_high[v] - a_low[v]
			c[v] = (0, 1)
			z13[v] = (1, 0)

		return a, b, b_hat, c, z13

	def set_up_inverter(self, invl):
		self.sdf_filename = invl.sdf_filename
		self.fv_gen_path = invl.fv_gen_path
		self._use_quadratic = invl._use_quadratic

		weight_filename = invl.weight_filename
		biases_filename = invl.biases_filename

		try:
			Q = _read_clean_lines(weight_filename)
			layer_sizes = _pop_ints(Q)
			self.layer_num = len(layer_sizes)
			self.num_K = layer_sizes[0]
			weight_tensor = list()
			for ell in layer_sizes[:-1]:
				weight_matrix = list()
				for _ in range(ell):
					weight_line = _pop_floats(Q)
					weight_matrix.append(weight_line)
				weight_tensor.append(weight_matrix)
			Q = _read_clean_lines(biases_filename)
			bias_matrix = list()
			# The input layer has no bias, and we define it to be 0
			bias_matrix.append([0] * layer_sizes[0])
			for ell in layer_sizes[1:]:
				bias_row = list()
				for _ in range(ell):
					bias_line = _pop_float(Q)
					bias_row.append(bias_line)
				bias_matrix.append(bias_row)
		except:
			raise ValueError(f"Error: Failed to read the inverter {self.name}!")

		self._initialize(weight_tensor, bias_matrix)
		self.ann_a, self.ann_b, self.ann_b_hat, self.ann_c, self.ann_z = self._initialize_constants()

		self.read_fringe_tree(invl.fringe_filename)
		self.prepare_fv_list(invl.desc_norm_selected_filename)
		self.prepare_desc_max_min(invl.desc_filename)
		self.prepare_value_max_min(invl.desc_norm_selected_filename, invl.values_filename)

	def set_up_inverter_for_GNN(self, weight_tensor, bias_matrix, layer_sizes):
		self.layer_num = len(layer_sizes)
		self._initialize(weight_tensor, bias_matrix)
		self.ann_a, self.ann_b, self.ann_b_hat, self.ann_c, self.ann_z = self._initialize_constants()
		
	def set_up_variables(self, asgml):
		a, b = self.ann_a, self.ann_b
		low = 0
		high = 2

		self.x[asgml.ind] = {node: pulp.LpVariable(f"ANN_x_{node}_ASGM{asgml.ind}", a[node][low], a[node][high], cat=pulp.LpContinuous) 
					for node in self.non_input_nodes}
		self.y[asgml.ind] = {node: pulp.LpVariable(f"ANN_y_{node}_ASGM{asgml.ind}", b[node][low], b[node][high], cat=pulp.LpContinuous) 
					for node in self.nodes}
		self.z[asgml.ind] = {node: pulp.LpVariable(f"ANN_z_{node}_ASGM{asgml.ind}", 0, 1, cat=pulp.LpBinary) 
					for node in self.non_input_nodes}

		### The output variable uses self.output as the variable name
		self.output[asgml.ind] = self.y[asgml.ind][(self.layer_num, 1)]

	def set_up_x_hat_variables(self, asgml):
		self.x_hat[asgml.ind] = {i: self.y[asgml.ind][(1, i)] for i in range(1, self.num_K + 1)}

	def add_constraints_C1(self, MILP, asgml):
		z = dict()
		for node in self.non_input_nodes:
			z[node] = (self.ann_z[node][0], self.z[asgml.ind][node], self.ann_z[node][1])

		for v in self.non_input_nodes:
			in_u = self.in_neighbors[v]
			# w_v = [self.weights[(u, v)] for u in in_u if u not in self.forbidden_node[asgml.ind]]
			# in_y = [self.y[asgml.ind][u] for u in in_u if u not in self.forbidden_node[asgml.ind]]
			w_v = [self.weights[(u, v)] for u in in_u]
			in_y = [self.y[asgml.ind][u] for u in in_u]

			MILP += self.x[asgml.ind][v] == pulp.lpDot(in_y, w_v) + self.biases[v], f"ANN_Output_node_{v}_ASGM{asgml.ind}"

			if v in self.layers[-1]:
				MILP += self.y[asgml.ind][v] == self.x[asgml.ind][v], f"ANN_ReLU_y_{v}_ASGM{asgml.ind}"
			else:
				MILP += self.x[asgml.ind][v] - self.ann_a[v][1] <= (self.ann_a[v][2] - self.ann_a[v][0]) * z[v][1], \
						f"ANN_ReLU_x_{v}_ASGM{asgml.ind}"
				if self.ann_a[v][2] > 0:
					MILP += self.x[asgml.ind][v] - self.ann_a[v][1] >= (self.ann_a[v][0] - self.ann_a[v][2]) * (1 - z[v][1]), \
							f"ANN_ReLU_x_{v}_lb_ASGM{asgml.ind}"
					for i in [0, 1]:
						MILP += self.y[asgml.ind][v] <= self.ann_c[v][i] * (self.x[asgml.ind][v] - self.ann_a[v][i]) + \
								self.ann_b[v][i] + self.ann_b_hat[v] * (1 + z[v][i + 1] - z[v][i]), \
								f"ANN_ReLU_y_{v}_{i}_ub_ASGM{asgml.ind}"
						MILP += self.y[asgml.ind][v] >= self.ann_c[v][i] * (self.x[asgml.ind][v] - self.ann_a[v][i]) + \
								self.ann_b[v][i] - self.ann_b_hat[v] * (1 + z[v][i + 1] - z[v][i]), \
								f"ANN_ReLU_y_{v}_{i}_lb_ASGM{asgml.ind}"
				else:
					MILP += z[v][1] == 0, f"ANN_ReLU_x_{v}_lb_ASGM{asgml.ind}"
					MILP += self.y[asgml.ind][v] == 0, f"ANN_ReLU_y_{v}_ASGM{asgml.ind}"

		return MILP

	def _activition(self, x):
		return max(0, x)  # ReLU

	def _inspection(self, fv_test):
		cur_layer = self.layers[0]
		vec_y = np.array(fv_test)
		for next_layer in self.layers[1:]:
			weight_matrix = np.array([[self.weights[(u, v)] for u in cur_layer] for v in next_layer])
			bias_vector = np.array([self.biases[v] for v in next_layer])
			vec_x = np.dot(weight_matrix, vec_y) + bias_vector

			if next_layer == self.layers[-1]:
				vec_y = np.array([x for x in vec_x])
			else:
				vec_y = np.array([self._activition(x) for x in vec_x])

			cur_layer = next_layer
		insp_value = vec_y[0]
		return insp_value

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

	