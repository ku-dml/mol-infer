from .inverter_base import *
# import pulp_modified as pulp
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, ClassVar
import pandas as pd
from src.utils import _read_clean_lines, _pop_int, _pop_floats, _pop_float, run_with_retry

@dataclass(frozen=False, kw_only=True)
class MLRInverter(InverterBase):
	### An inverter class for linear regression
	name: str
	_type: ClassVar[Literal["INV_MLR"]] = "INV_MLR"

	def __post_init__(self):
		self.output = dict()
		self.index_set = dict() # store the relation between the indices of fringe trees in fv files and in seed graph files, index_set[index_in_fv] = index_in_seed_graph
		self.x_hat = dict() # store the normalized feature vector value
		self.std_eps = 1e-5

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

		linreg_filename = invl.linreg_filename
		try:
			Q = _read_clean_lines(linreg_filename)
			self.num_K = _pop_int(Q) ## read dimension K
			self.weight = _pop_floats(Q) ## read weight w
			self.bias = _pop_float(Q) ## read bias b
		except:
			raise ValueError(f"Error: Failed to read the inverter {self.name}!")

		self.read_fringe_tree(invl.fringe_filename)
		self.prepare_fv_list(invl.desc_norm_selected_filename)
		self.prepare_desc_max_min(invl.desc_filename)
		self.prepare_value_max_min(invl.desc_norm_selected_filename, invl.values_filename)
		
	def set_up_variables(self, asgml):
		### The output variable uses self.output as the variable name
		self.output[asgml.ind] = pulp.LpVariable(f"LR_y_ASGM{asgml.ind}", cat=pulp.LpContinuous)

	def set_up_x_hat_variables(self, asgml):
		self.x_hat[asgml.ind] = {i: pulp.LpVariable(f"x_hat({i})_ASGM{asgml.ind}") for i in range(1, self.num_K + 1)}
				
	def add_constraints_C1(self, MILP, asgml):
		MILP += self.output[asgml.ind] == pulp.lpSum([self.weight[i - 1] * self.x_hat[asgml.ind][i] for i in range(1, self.num_K + 1)]) + self.bias, f"LR_lr_ASGM{asgml.ind}"
		return MILP

	def _inspection(self, fv_test):
		insp_value = sum([self.weight[i - 1] * fv_test[i - 1] for i in range(1, self.num_K + 1)]) + self.bias
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

	