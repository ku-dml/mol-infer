from dataclasses import dataclass, field
from .seed_graph_base import SeedGraphBase
from src.config import AssignmentLoader
from src_inverter.inverter_base import InverterBase
from abc import ABC, abstractmethod
from typing import Dict, List, Literal, Optional, ClassVar
from src.utils import run_with_retry
from src.config import AssignmentLoader
# import pulp_modified as pulp
import pandas as pd
import copy
import numpy as np

@dataclass(frozen=False, kw_only=True)
class AssignmentBase(ABC):
	sg: SeedGraphBase | List[SeedGraphBase]
	inv: InverterBase | List[InverterBase]
	asgml: AssignmentLoader
	_type: ClassVar[Literal["BASE"]] = "BASE"

	@property
	def type(self):
		return type(self)._type

	@abstractmethod
	def get_fv_variable(self):
		pass

	@abstractmethod
	def set_up_variables_norm(self):
		pass

	@abstractmethod
	def connect_fv_norm_x_hat(self):
		pass

	@abstractmethod
	def add_constraints_lb_ub(self):
		pass

	@abstractmethod
	def get_output(self):
		pass

@dataclass(frozen=False, kw_only=True)
class AssignmentNormal(AssignmentBase):
	sg: SeedGraphBase
	inv: InverterBase
	_type: ClassVar[Literal["NORMAL"]] = "NORMAL"
	_support_quadratic: ClassVar[bool] = False

	def get_fv_variable(self, fv_name):
		### TO DO, add other kind of descriptors

		sg = self.sg
		asgml = self.asgml
		inv = self.inv

		ans = None
		flag = 0
		if fv_name == "CID":
			flag = -1
			pass
		elif fv_name == "rank":
			ans = sg.rank_G
		elif fv_name == "n":
			ans = sg.n_G
		elif fv_name == "n_in":
			ans = sg.n_G_int
		elif fv_name == "link_edge":
			### polymer
			ans = sg.n_lnk
		elif fv_name == "ms":
			# ans = None# / (10 * n_G)
			ans = sg.mass_n
		elif fv_name[0] == "d":
			if fv_name[3] == "i":                
				j = int(fv_name[6])
				ans = sg.dg_int[j]
			else:
				j = int(fv_name[3])
				ans = sg.dg[j]
		elif fv_name == "bd_in_2":
			ans = sg.bd_int[2]
		elif fv_name == "bd_in_3":
			ans = sg.bd_int[3]
		elif fv_name[0] == "n" and fv_name[3] == "i":
			ele = fv_name[6:] # may need change, Zhu, 0109
			if ele in sg.Lambda_int:
				ans = sg.na_int[ele]
			else:
				ans = 0
		elif fv_name[0] == "n" and fv_name[3] == "e":
			ele = fv_name[6:]
			if ele in sg.Lambda_ex:
				ans = sg.na_ex[ele]
			else:
				ans = 0
		elif fv_name[0] == "F" and fv_name[1] == "C":
			index = int(fv_name.split("_")[1])
			if index in inv.index_set[asgml.ind].keys():
				ans = sg.fc[inv.index_set[asgml.ind][index]]
			else:
				ans = 0
				flag = -1
		elif fv_name[0] == "L":
			lines = fv_name.split("_")
			nu = (lines[1], lines[2], int(lines[3]))
			if nu in sg.Gamma_lf_ac:
				ans = sg.ac_lf[nu]
			else:
				ans = 0
				flag = -1
		elif fv_name[0] == "l":
			### polymer
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[pos[1] + 1: pos[2]]
			d1 = int(str_tmp[pos[2] + 1: pos[3]])
			ele2 = str_tmp[pos[3] + 1: pos[4]]
			d2 = int(str_tmp[pos[4] + 1: pos[5]])
			mul = int(str_tmp[-1])
			ec_tmp1 = ((ele1, d1), (ele2, d2), mul)
			ec_tmp2 = ((ele2, d2), (ele1, d1), mul)

			if ec_tmp1 in sg.Gamma_lnk_less or ec_tmp1 in sg.Gamma_lnk_equal:
				ans = sg.ec_lnk[ec_tmp1]
			elif ec_tmp2 in sg.Gamma_lnk_less or ec_tmp2 in sg.Gamma_lnk_equal:
				ans = sg.ec_lnk[ec_tmp2]
			else:
				ans = 0
		elif fv_name[0] == "C" and fv_name[1] == "S":
			### polymer
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[pos[2] + 1: pos[3]]
			d1 = int(str_tmp[-1])

			if (ele1, d1) in sg.Lambda_dg_int:
				# descriptors[i] = delta_cnt[(1, (ele1, d1))] + delta_cnt[(2, (ele1, d1))]
				ans = pulp.lpSum(sg.delta_cnt[(mu1, mu2, m)] 
					for (mu1, mu2, m) in (sg.Gamma_cnt_less + sg.Gamma_cnt_equal + sg.Gamma_cnt_great) if mu1 == (ele1, d1)) + \
					pulp.lpSum(sg.delta_cnt[(mu1, mu2, m)] 
					for (mu1, mu2, m) in (sg.Gamma_cnt_less + sg.Gamma_cnt_equal + sg.Gamma_cnt_great) if mu2 == (ele1, d1))
			else:
				ans = 0
		else:
			# ec
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[: pos[0]]
			d1 = int(str_tmp[pos[0] + 1: pos[1]])
			ele2 = str_tmp[pos[1] + 1: pos[2]]
			d2 = int(str_tmp[pos[2] + 1: pos[3]])
			mul = int(str_tmp[-1])
			ec_tmp1 = ((ele1, d1), (ele2, d2), mul)
			ec_tmp2 = ((ele2, d2), (ele1, d1), mul)

			if ec_tmp1 in sg.Gamma_int_less or ec_tmp1 in sg.Gamma_int_equal:
				ans = sg.ec_int[ec_tmp1]
			elif ec_tmp2 in sg.Gamma_int_less or ec_tmp2 in sg.Gamma_int_equal:
				ans = sg.ec_int[ec_tmp2]
			else:
				flag = -1
				ans = 0

		return ans, flag

	def set_up_variables_norm(self, asgml, num_K):
		self.fv_norm = {i: pulp.LpVariable(f"fv_norm({i})_ASGM{asgml.ind}") for i in range(1, num_K + 1)}

	def add_constraints_norm(self, 
		MILP, sg, fv_norm, 
		asgml, num_K, fv_list, descriptors, max_dcp, min_dcp, forbidden_node, std_eps
	):

		for i in range(1, num_K + 1):
			fv_name = fv_list[i - 1]
			if fv_name not in forbidden_node:
				if max_dcp[fv_name] == min_dcp[fv_name]:
					MILP += fv_norm[i] == 0, f"connect_norm_ASGM{asgml.ind}_{i}"
				else:
					MILP += fv_norm[i] >= (descriptors[fv_name] - std_eps - min_dcp[fv_name]) / (max_dcp[fv_name] - min_dcp[fv_name]), f"connect_norm_ASGM{asgml.ind}_{i}-1"
					MILP += fv_norm[i] <= (descriptors[fv_name] + std_eps - min_dcp[fv_name]) / (max_dcp[fv_name] - min_dcp[fv_name]), f"connect_norm_ASGM{asgml.ind}_{i}-2"
			else:
				MILP += fv_norm[i] == 0, f"connect_norm_ASGM{asgml.ind}_{i}"

		return MILP

	def connect_fv_norm_x_hat(self, MILP, asgml, fv_norm, num_K, x_hat):
		for i in range(1, num_K + 1):
			MILP += fv_norm[i] == x_hat[i], f"connect_ASGM{asgml.ind}_{i}"

		return MILP

	def add_constraints_lb_ub(self, MILP, output, asgml, y_min, y_max):
		### normalize the target value
		target_value_lb = (asgml.lb - y_min) / (y_max - y_min)
		target_value_ub = (asgml.ub - y_min) / (y_max - y_min)
		
		MILP += output >= target_value_lb, f"LB_ASGM{asgml.ind}"
		MILP += output <= target_value_ub, f"UB_ASGM{asgml.ind}"

		return MILP

	### connect C1 and C2
	def connect_C1_C2(self, MILP):
		self.inv.set_up_variables(self.asgml)
		self.output = self.inv.output[self.asgml.ind]
		self.inv.set_up_x_hat_variables(self.asgml)
		self.set_up_variables_norm(self.asgml, self.inv.num_K)
		self.inv.prepare_fringe_tree_index_set(self.asgml.ind, self.sg.strF)

		self.descriptors = dict()
		self.forbidden_node = list()
		for i in range(1, self.inv.num_K + 1):
			fv_name = self.inv.fv_list[i - 1]
			self.descriptors[fv_name], flag = self.get_fv_variable(fv_name)
			if flag == -1:
				self.forbidden_node.append(fv_name)

		MILP = self.add_constraints_norm(MILP, self.sg, self.fv_norm, self.asgml, self.inv.num_K, self.inv.fv_list, self.descriptors, self.inv.max_dcp, self.inv.min_dcp, self.forbidden_node, self.inv.std_eps)
		MILP = self.inv.add_constraints_C1(MILP, self.asgml)
		MILP = self.connect_fv_norm_x_hat(MILP, self.asgml, self.fv_norm, self.inv.num_K, self.inv.x_hat[self.asgml.ind])
		MILP = self.add_constraints_lb_ub(MILP, self.inv.output[self.asgml.ind], self.asgml, self.inv.y_min, self.inv.y_max)
		
		return MILP

	def _inspection_pre(self, milp_sdf_filename, test_prefix):
		cmd = [self.inv.fv_gen_path, self.inv.sdf_filename, f"test_tmp_ASGM{self.asgml.ind}", milp_sdf_filename, test_prefix]
		run_with_retry(cmd)
		test_desc_norm_filename = f"{test_prefix}_desc_norm.csv"
		df = pd.read_csv(test_desc_norm_filename, sep=",", index_col=0)
		idx = df.index[0]

		fv_test = []
		for fv_name in self.inv.fv_list:
			if not fv_name.startswith("F"):
				fv_test.append(float(df.at[idx, fv_name]))
			else:
				# fringe configuration
				ind_fv_name = fv_name.find('_', 3)
				match = None
				for fv_name_test in df.columns:
					ind2 = fv_name_test.find('_', 3)
					if ind2 != -1 and ind_fv_name != -1 and fv_name[ind_fv_name:] == fv_name_test[ind2:]:
						match = fv_name_test
						break
				fv_test.append(float(df.at[idx, match]) if match else 0.0)

		return fv_test

	def inspection(self, sdf_filename_prefix, test_prefix):
		sdf_filename = f"{sdf_filename_prefix}.sdf"
		fv_test = self._inspection_pre(sdf_filename, test_prefix)
		insp_value = self.inv.inspection(self.asgml, fv_test)
		insp_value_scaled = insp_value * (self.inv.y_max - self.inv.y_min) + self.inv.y_min

		return insp_value, insp_value_scaled

	def get_output(self):
		y_MILP = self.output.value()
		y_scale = y_MILP * (self.inv.y_max - self.inv.y_min) + self.inv.y_min
		
		return y_MILP, y_scale

@dataclass(frozen=False, kw_only=True)
class AssignmentQuad(AssignmentBase):
	sg: SeedGraphBase
	inv: InverterBase
	_type: ClassVar[Literal["QUAD"]] = "QUAD"
	_support_quadratic: ClassVar[bool] = True

	def get_fv_variable(self, fv_name):
		### TO DO, add other kind of descriptors

		sg = self.sg
		asgml = self.asgml
		inv = self.inv

		ans = None
		flag = 0
		if fv_name == "CID":
			flag = -1
			pass
		elif fv_name == "rank":
			ans = sg.rank_G
		elif fv_name == "n":
			ans = sg.n_G
		elif fv_name == "n_in":
			ans = sg.n_G_int
		elif fv_name == "link_edge":
			### polymer
			ans = sg.n_lnk
		elif fv_name == "ms":
			# ans = None# / (10 * n_G)
			ans = sg.mass_n
		elif fv_name[0] == "d":
			if fv_name[3] == "i":                
				j = int(fv_name[6])
				ans = sg.dg_int[j]
			else:
				j = int(fv_name[3])
				ans = sg.dg[j]
		elif fv_name == "bd_in_2":
			ans = sg.bd_int[2]
		elif fv_name == "bd_in_3":
			ans = sg.bd_int[3]
		elif fv_name[0] == "n" and fv_name[3] == "i":
			ele = fv_name[6:] # may need change, Zhu, 0109
			if ele in sg.Lambda_int:
				ans = sg.na_int[ele]
			else:
				ans = 0
		elif fv_name[0] == "n" and fv_name[3] == "e":
			ele = fv_name[6:]
			if ele in sg.Lambda_ex:
				ans = sg.na_ex[ele]
			else:
				ans = 0
		elif fv_name[0] == "F" and fv_name[1] == "C":
			index = int(fv_name.split("_")[1])
			if index in inv.index_set[asgml.ind].keys():
				ans = sg.fc[inv.index_set[asgml.ind][index]]
			else:
				ans = 0
				flag = -1
		elif fv_name[0] == "L":
			lines = fv_name.split("_")
			nu = (lines[1], lines[2], int(lines[3]))
			if nu in sg.Gamma_lf_ac:
				ans = sg.ac_lf[nu]
			else:
				ans = 0
				flag = -1
		elif fv_name[0] == "l":
			### polymer
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[pos[1] + 1: pos[2]]
			d1 = int(str_tmp[pos[2] + 1: pos[3]])
			ele2 = str_tmp[pos[3] + 1: pos[4]]
			d2 = int(str_tmp[pos[4] + 1: pos[5]])
			mul = int(str_tmp[-1])
			ec_tmp1 = ((ele1, d1), (ele2, d2), mul)
			ec_tmp2 = ((ele2, d2), (ele1, d1), mul)

			if ec_tmp1 in sg.Gamma_lnk_less or ec_tmp1 in sg.Gamma_lnk_equal:
				ans = sg.ec_lnk[ec_tmp1]
			elif ec_tmp2 in sg.Gamma_lnk_less or ec_tmp2 in sg.Gamma_lnk_equal:
				ans = sg.ec_lnk[ec_tmp2]
			else:
				ans = 0
		elif fv_name[0] == "C" and fv_name[1] == "S":
			### polymer
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[pos[2] + 1: pos[3]]
			d1 = int(str_tmp[-1])

			if (ele1, d1) in sg.Lambda_dg_int:
				# descriptors[i] = delta_cnt[(1, (ele1, d1))] + delta_cnt[(2, (ele1, d1))]
				ans = pulp.lpSum(sg.delta_cnt[(mu1, mu2, m)] 
					for (mu1, mu2, m) in (sg.Gamma_cnt_less + sg.Gamma_cnt_equal + sg.Gamma_cnt_great) if mu1 == (ele1, d1)) + \
					pulp.lpSum(sg.delta_cnt[(mu1, mu2, m)] 
					for (mu1, mu2, m) in (sg.Gamma_cnt_less + sg.Gamma_cnt_equal + sg.Gamma_cnt_great) if mu2 == (ele1, d1))
			else:
				ans = 0
		else:
			# ec
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[: pos[0]]
			d1 = int(str_tmp[pos[0] + 1: pos[1]])
			ele2 = str_tmp[pos[1] + 1: pos[2]]
			d2 = int(str_tmp[pos[2] + 1: pos[3]])
			mul = int(str_tmp[-1])
			ec_tmp1 = ((ele1, d1), (ele2, d2), mul)
			ec_tmp2 = ((ele2, d2), (ele1, d1), mul)

			if ec_tmp1 in sg.Gamma_int_less or ec_tmp1 in sg.Gamma_int_equal:
				ans = sg.ec_int[ec_tmp1]
			elif ec_tmp2 in sg.Gamma_int_less or ec_tmp2 in sg.Gamma_int_equal:
				ans = sg.ec_int[ec_tmp2]
			else:
				flag = -1
				ans = 0

		return ans, flag

	def set_up_variables_norm(self, asgml, num_K):
		self.fv_norm = {i: pulp.LpVariable(f"fv_norm({i})_ASGM{asgml.ind}") for i in range(1, num_K + 1)}

	def add_constraints_norm(self, 
		MILP, sg, fv_norm, 
		asgml, num_K, fv_list, descriptors, max_dcp, min_dcp, forbidden_node, std_eps
	):

		for i in range(1, num_K + 1):
			fv_name = fv_list[i - 1]
			if fv_name not in forbidden_node:
				if max_dcp[fv_name] == min_dcp[fv_name]:
					MILP += fv_norm[i] == 0, f"connect_norm_ASGM{asgml.ind}_{i}"
				else:
					MILP += fv_norm[i] >= (descriptors[fv_name] - std_eps - min_dcp[fv_name]) / (max_dcp[fv_name] - min_dcp[fv_name]), f"connect_norm_ASGM{asgml.ind}_{i}-1"
					MILP += fv_norm[i] <= (descriptors[fv_name] + std_eps - min_dcp[fv_name]) / (max_dcp[fv_name] - min_dcp[fv_name]), f"connect_norm_ASGM{asgml.ind}_{i}-2"
			else:
				MILP += fv_norm[i] == 0, f"connect_norm_ASGM{asgml.ind}_{i}"

		return MILP

	def connect_fv_norm_x_hat(self, MILP, asgml, fv_norm, inv_num_K, x_hat):
		self._p = 6
		self._p_q = self._p + 3
		self._M = 9
		self.x_q = {(i, q): pulp.LpVariable(f"Q_xq({i},{q})_ASGM{asgml.ind}", 0, 1, cat=pulp.LpBinary)
					for (i, info) in self.descriptors_info.items() if len(info) == 2 for q in range(self._p_q + 1)}
		self.z_q = {(i, q): pulp.LpVariable(f"Q_zq({i},{q})_ASGM{asgml.ind}", 0, cat=pulp.LpContinuous)
					for (i, info) in self.descriptors_info.items() if len(info) == 2 for q in range(self._p_q + 1)}
		self.z_tmp = {i: pulp.LpVariable(f"Q_z_tmp({i})_ASGM{asgml.ind}", 0, cat=pulp.LpContinuous)
					for (i, info) in self.descriptors_info.items() if len(info) == 2}

		for i in range(1, inv_num_K + 1):
			info = self.descriptors_info[i]
			if len(info) == 1:
				# linear descriptor
				fv_name = info[0]
				MILP += x_hat[i] == fv_norm[self.descriptors_ind[fv_name]], f"connect_ASGM{asgml.ind}_{i}"
			else:
				# quadratic descriptor
				fv_name1 = info[0][0]
				fv_name2 = info[1][0]
				q_ind_1 = self.descriptors_ind[fv_name1]
				q_ind_2 = self.descriptors_ind[fv_name2]
				
				MILP += (2**(self._p + 1) - 1) * fv_norm[q_ind_1] >= \
						pulp.lpSum([2**q * self.x_q[(i, q)] for q in range(self._p_q + 1)]) - 1, f"connect_Q-{i}-1-1_ASGM{asgml.ind}"
				MILP += (2**(self._p + 1) - 1) * fv_norm[q_ind_1] <= \
						pulp.lpSum([2**q * self.x_q[(i, q)] for q in range(self._p_q + 1)]), f"connect_Q-{i}-1-2_ASGM{asgml.ind}"
				for q in range(self._p_q + 1):
					MILP += self.z_q[(i, q)] <= self._M * self.x_q[(i, q)], f"connect_Q-{i}-2-1-{q}_ASGM{asgml.ind}"
					MILP += self.z_q[(i, q)] >= fv_norm[q_ind_2] - self._M * (1 - self.x_q[(i, q)]), f"connect_Q-{i}-2-2-{q}-1_ASGM{asgml.ind}"
					MILP += self.z_q[(i, q)] <= fv_norm[q_ind_2] + self._M * (1 - self.x_q[(i, q)]), f"connect_Q-{i}-2-2-{q}-2_ASGM{asgml.ind}"
				MILP += self.z_tmp[i] * (2**(self._p + 1) - 1) == \
					pulp.lpSum([2**q * self.z_q[(i, q)] for q in range(self._p_q + 1)]), f"connect_Q-{i}-3_ASGM{asgml.ind}"

				if info[0][1] == -1 and info[1][1] == -1:
					MILP += x_hat[i] == 1 - fv_norm[q_ind_1] - fv_norm[q_ind_2] + self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"
				elif info[0][1] == -1:
					MILP += x_hat[i] == fv_norm[q_ind_2] - self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"
				elif info[1][1] == -1:
					MILP += x_hat[i] == fv_norm[q_ind_1] - self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"
				else:
					MILP += x_hat[i] == self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"

		return MILP				

	def add_constraints_lb_ub(self, MILP, output, asgml, y_min, y_max):
		### normalize the target value
		target_value_lb = (asgml.lb - y_min) / (y_max - y_min)
		target_value_ub = (asgml.ub - y_min) / (y_max - y_min)
		
		MILP += output >= target_value_lb, f"LB_ASGM{asgml.ind}"
		MILP += output <= target_value_ub, f"UB_ASGM{asgml.ind}"

		return MILP

	def prepare_quadratic_desc(self):
		self.descriptors_info = dict()      # a dict stores all information of descriptors, including linear and quadratic
		self.descriptors_ind = dict()      # a dict stores the index of the linear descriptors
		self.descriptors = dict()          # a dict stores the completer linear descriptors
		self.descriptors_list = list()     # a list to store names of linear descritpors
		self.forbidden_node = list()

		for i in range(1, self.inv.num_K + 1):
			fv_name = self.inv.fv_list[i - 1]
			if "*" not in fv_name:
				# linear descriptor
				self.descriptors_list.append(fv_name)
				self.descriptors[fv_name], flag = self.get_fv_variable(fv_name)
				self.descriptors_info[i] = [fv_name]
				self.descriptors_ind[fv_name] = len(self.descriptors_list)
				if flag == -1:
					self.forbidden_node.append(fv_name)
			elif "*" in fv_name:
				# quadratic descriptor
				sp_loc1 = fv_name.find("*")
				fv_name1 = fv_name[:sp_loc1]
				fv_name2 = fv_name[sp_loc1 + 1:]
				if "1-" in fv_name1:
					minus_1 = -1
					fv_name1 = fv_name1[3:-1]
				else:
					minus_1 = 1
				if "1-" in fv_name2:
					minus_2 = -1
					fv_name2 = fv_name2[3:-1]
				else:
					minus_2 = 1
				
				if fv_name1 not in self.descriptors:
					self.descriptors_list.append(fv_name1)
					self.descriptors[fv_name1], flag = self.get_fv_variable(fv_name1)
					self.descriptors_ind[fv_name1] = len(self.descriptors_list)
					if flag == -1:
						self.forbidden_node.append(fv_name1)
				if fv_name2 not in self.descriptors:
					self.descriptors_list.append(fv_name2)
					self.descriptors[fv_name2], flag = self.get_fv_variable(fv_name2)
					self.descriptors_ind[fv_name2] = len(self.descriptors_list)
					if flag == -1:
						self.forbidden_node.append(fv_name2)
				self.descriptors_info[i] = [(fv_name1, minus_1), (fv_name2, minus_2)]

	### connect C1 and C2
	def connect_C1_C2(self, MILP):
		self.inv.set_up_variables(self.asgml)
		self.output = self.inv.output[self.asgml.ind]
		self.inv.set_up_x_hat_variables(self.asgml)
		self.inv.prepare_fringe_tree_index_set(self.asgml.ind, self.sg.strF)
		self.prepare_quadratic_desc()
		self.num_K = len(self.descriptors_list)
		self.set_up_variables_norm(self.asgml, self.num_K)
		
		MILP = self.add_constraints_norm(MILP, self.sg, self.fv_norm, self.asgml, self.num_K, self.descriptors_list, self.descriptors, self.inv.max_dcp, self.inv.min_dcp, self.forbidden_node, self.inv.std_eps)
		MILP = self.inv.add_constraints_C1(MILP, self.asgml)
		MILP = self.connect_fv_norm_x_hat(MILP, self.asgml, self.fv_norm, self.inv.num_K, self.inv.x_hat[self.asgml.ind])
		MILP = self.add_constraints_lb_ub(MILP, self.inv.output[self.asgml.ind], self.asgml, self.inv.y_min, self.inv.y_max)
		
		return MILP

	def _inspection_pre(self, milp_sdf_filename, test_prefix):
		cmd = [self.inv.fv_gen_path, self.inv.sdf_filename, f"test_tmp_ASGM{self.asgml.ind}", milp_sdf_filename, test_prefix]
		run_with_retry(cmd)
		test_desc_norm_filename = f"{test_prefix}_desc_norm.csv"
		df = pd.read_csv(test_desc_norm_filename, sep=",", index_col=0)
		idx = df.index[0]

		fv_test = []
		for fv_name in self.descriptors_list:
			if not fv_name.startswith("F"):
				fv_test.append(float(df.at[idx, fv_name]))
			else:
				# fringe configuration
				ind_fv_name = fv_name.find('_', 3)
				match = None
				for fv_name_test in df.columns:
					ind2 = fv_name_test.find('_', 3)
					if ind2 != -1 and ind_fv_name != -1 and fv_name[ind_fv_name:] == fv_name_test[ind2:]:
						match = fv_name_test
						break
				fv_test.append(float(df.at[idx, match]) if match else 0.0)

		return fv_test

	def inspection(self, sdf_filename_prefix, test_prefix):
		sdf_filename = f"{sdf_filename_prefix}.sdf"
		fv_test_linear = self._inspection_pre(sdf_filename, test_prefix)
		fv_test = []
		for i in range(1, self.inv.num_K + 1):
			info = self.descriptors_info[i]
			if len(info) == 1:
				# linear descriptor
				fv_name = info[0]
				fv_test.append(fv_test_linear[self.descriptors_ind[fv_name] - 1])
			else:
				# quadratic descriptor
				fv_name1 = info[0][0]
				fv_name2 = info[1][0]
				q_ind_1 = self.descriptors_ind[fv_name1]
				q_ind_2 = self.descriptors_ind[fv_name2]
				_x = fv_test_linear[q_ind_1 - 1]
				_y = fv_test_linear[q_ind_2 - 1]
				if info[0][1] == -1:
					_x = 1 - _x
				if info[1][1] == -1:
					_y = 1 - _y
				fv_test.append(_x * _y)

		insp_value = self.inv.inspection(self.asgml, fv_test, std_eps=1/(2**(self._p + 1) - 1))
		insp_value_scaled = insp_value * (self.inv.y_max - self.inv.y_min) + self.inv.y_min

		return insp_value, insp_value_scaled

	def get_output(self):
		y_MILP = self.output.value()
		y_scale = y_MILP * (self.inv.y_max - self.inv.y_min) + self.inv.y_min
		
		return y_MILP, y_scale

@dataclass(frozen=False, kw_only=True)
class AssignmentGNN(AssignmentBase):
	sg: SeedGraphBase
	inv: InverterBase
	_type: ClassVar[Literal["GNN"]] = "GNN"
	_support_quadratic: ClassVar[bool] = False

	def set_up_variables_norm(self):
		# no fv for GNN
		pass

	def add_constraints_norm(self):
		# no fv for GNN
		pass

	def connect_fv_norm_x_hat(self):
		# no fv for GNN
		pass

	def get_fv_variable(self):
		# no fv for GNN
		pass

	def add_constraints_lb_ub(self, MILP, output, asgml, y_mean, y_std):
		### standardize the target value
		target_value_lb = (asgml.lb - y_mean) / y_std
		target_value_ub = (asgml.ub - y_mean) / y_std
		
		MILP += output >= target_value_lb, f"LB_ASGM{asgml.ind}"
		MILP += output <= target_value_ub, f"UB_ASGM{asgml.ind}"

		return MILP

	### connect C1 and C2
	def connect_C1_C2(self, MILP):
		self.inv.build_ft_dataset(self.inv.GNN, self.sg.set_F)
		self.inv.set_up_variables(self.asgml, self.sg)
		self.output = self.inv.output[self.asgml.ind]
		
		MILP = self.inv.add_constraints_C1(MILP, self.asgml, self.sg)
		MILP = self.add_constraints_lb_ub(MILP, self.inv.output[self.asgml.ind], self.asgml, self.inv.y_mean, self.inv.y_std)
		
		return MILP

	def inspection(self, sdf_filename_prefix, test_prefix):
		sdf_filename = f"{sdf_filename_prefix}.sdf"
		insp_value = self.inv.inspection(sdf_filename)
		insp_value_scaled = insp_value * self.inv.y_std + self.inv.y_mean

		return insp_value, insp_value_scaled

	def get_output(self):
		y_MILP = self.output.value()
		y_scale = y_MILP * self.inv.y_std + self.inv.y_mean
		
		return y_MILP, y_scale

@dataclass(frozen=False, kw_only=True)
class AssignmentMixingVector(AssignmentBase):
	sg: List[SeedGraphBase]
	inv: InverterBase
	_type: ClassVar[Literal["MV"]] = "MV"
	_support_quadratic: bool = field(init=False)
	ratio: Dict[str, float] = field(init=False, default_factory=dict)

	def __post_init__(self):
		self._support_quadratic = self.inv._use_quadratic
		ratio_tmp = self.asgml.extra.get("ratio")
		if ratio_tmp is None:
			ratio_tmp = [1.0 / len(self.sg)] * len(self.sg)

		for (sg, r) in zip(self.sg, ratio_tmp):
			self.ratio[sg.name] = r
	
	def get_fv_variable(self, fv_name, sg):
		### TO DO, add other kind of descriptors

		asgml = self.asgml
		inv = self.inv

		ans = None
		flag = 0
		if fv_name == "CID":
			flag = -1
			pass
		elif fv_name == "rank":
			ans = sg.rank_G
		elif fv_name == "n":
			ans = sg.n_G
		elif fv_name == "n_in":
			ans = sg.n_G_int
		elif fv_name == "link_edge":
			### polymer
			ans = sg.n_lnk
		elif fv_name == "ms":
			# ans = None# / (10 * n_G)
			ans = sg.mass_n
		elif fv_name[0] == "d":
			if fv_name[3] == "i":                
				j = int(fv_name[6])
				ans = sg.dg_int[j]
			else:
				j = int(fv_name[3])
				ans = sg.dg[j]
		elif fv_name == "bd_in_2":
			ans = sg.bd_int[2]
		elif fv_name == "bd_in_3":
			ans = sg.bd_int[3]
		elif fv_name[0] == "n" and fv_name[3] == "i":
			ele = fv_name[6:] # may need change, Zhu, 0109
			if ele in sg.Lambda_int:
				ans = sg.na_int[ele]
			else:
				ans = 0
		elif fv_name[0] == "n" and fv_name[3] == "e":
			ele = fv_name[6:]
			if ele in sg.Lambda_ex:
				ans = sg.na_ex[ele]
			else:
				ans = 0
		elif fv_name[0] == "F" and fv_name[1] == "C":
			index = int(fv_name.split("_")[1])
			if index in inv.index_set[(asgml.ind, sg.name)].keys():
				ans = sg.fc[inv.index_set[(asgml.ind, sg.name)][index]]
			else:
				ans = 0
				flag = -1
		elif fv_name[0] == "L":
			lines = fv_name.split("_")
			nu = (lines[1], lines[2], int(lines[3]))
			if nu in sg.Gamma_lf_ac:
				ans = sg.ac_lf[nu]
			else:
				ans = 0
				flag = -1
		elif fv_name[0] == "l":
			### polymer
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[pos[1] + 1: pos[2]]
			d1 = int(str_tmp[pos[2] + 1: pos[3]])
			ele2 = str_tmp[pos[3] + 1: pos[4]]
			d2 = int(str_tmp[pos[4] + 1: pos[5]])
			mul = int(str_tmp[-1])
			ec_tmp1 = ((ele1, d1), (ele2, d2), mul)
			ec_tmp2 = ((ele2, d2), (ele1, d1), mul)

			if ec_tmp1 in sg.Gamma_lnk_less or ec_tmp1 in sg.Gamma_lnk_equal:
				ans = sg.ec_lnk[ec_tmp1]
			elif ec_tmp2 in sg.Gamma_lnk_less or ec_tmp2 in sg.Gamma_lnk_equal:
				ans = sg.ec_lnk[ec_tmp2]
			else:
				ans = 0
		elif fv_name[0] == "C" and fv_name[1] == "S":
			### polymer
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[pos[2] + 1: pos[3]]
			d1 = int(str_tmp[-1])

			if (ele1, d1) in sg.Lambda_dg_int:
				# descriptors[i] = delta_cnt[(1, (ele1, d1))] + delta_cnt[(2, (ele1, d1))]
				ans = pulp.lpSum(sg.delta_cnt[(mu1, mu2, m)] 
					for (mu1, mu2, m) in (sg.Gamma_cnt_less + sg.Gamma_cnt_equal + sg.Gamma_cnt_great) if mu1 == (ele1, d1)) + \
					pulp.lpSum(sg.delta_cnt[(mu1, mu2, m)] 
					for (mu1, mu2, m) in (sg.Gamma_cnt_less + sg.Gamma_cnt_equal + sg.Gamma_cnt_great) if mu2 == (ele1, d1))
			else:
				ans = 0
		else:
			# ec
			str_tmp = fv_name
			pos = list()
			posnum = -1
			for k in range(len(str_tmp)):
				if str_tmp[k] == '_':
					pos.append(k)
			ele1 = str_tmp[: pos[0]]
			d1 = int(str_tmp[pos[0] + 1: pos[1]])
			ele2 = str_tmp[pos[1] + 1: pos[2]]
			d2 = int(str_tmp[pos[2] + 1: pos[3]])
			mul = int(str_tmp[-1])
			ec_tmp1 = ((ele1, d1), (ele2, d2), mul)
			ec_tmp2 = ((ele2, d2), (ele1, d1), mul)

			if ec_tmp1 in sg.Gamma_int_less or ec_tmp1 in sg.Gamma_int_equal:
				ans = sg.ec_int[ec_tmp1]
			elif ec_tmp2 in sg.Gamma_int_less or ec_tmp2 in sg.Gamma_int_equal:
				ans = sg.ec_int[ec_tmp2]
			else:
				flag = -1
				ans = 0

		return ans, flag

	def set_up_variables_norm(self, asgml, num_K):
		self.fv_norm = {i: pulp.LpVariable(f"fv_norm({i})_ASGM{asgml.ind}") for i in range(1, num_K + 1)}

	def add_constraints_norm(self, 
		MILP, sg, fv_norm, 
		asgml, num_K, fv_list, descriptors, max_dcp, min_dcp, forbidden_node, std_eps
	):
		# This is where MV model differs form normal model
		# The fv is weighted *BEFORE* normalization
		for i in range(1, num_K + 1):
			fv_name = fv_list[i - 1]
			if fv_name not in forbidden_node:
				if max_dcp[fv_name] == min_dcp[fv_name]:
					MILP += fv_norm[i] == 0, f"connect_norm_ASGM{asgml.ind}_{i}"
				else:
					MILP += fv_norm[i] >= (pulp.lpSum([self.ratio[sg.name] * descriptors[(sg.name, fv_name)] for sg in self.sg]) 
						- std_eps - min_dcp[fv_name]) / (max_dcp[fv_name] - min_dcp[fv_name]), f"connect_norm_ASGM{asgml.ind}_{i}-1"
					MILP += fv_norm[i] <= (pulp.lpSum([self.ratio[sg.name] * descriptors[(sg.name, fv_name)] for sg in self.sg]) 
						+ std_eps - min_dcp[fv_name]) / (max_dcp[fv_name] - min_dcp[fv_name]), f"connect_norm_ASGM{asgml.ind}_{i}-2"
			else:
				MILP += fv_norm[i] == 0, f"connect_norm_ASGM{asgml.ind}_{i}"

		return MILP

	def connect_fv_norm_x_hat(self, MILP, asgml, fv_norm, inv_num_K, x_hat):
		self._p = 6
		self._p_q = self._p + 3
		self._M = 9
		self.x_q = {(i, q): pulp.LpVariable(f"Q_xq({i},{q})_ASGM{asgml.ind}", 0, 1, cat=pulp.LpBinary)
					for (i, info) in self.descriptors_info.items() if len(info) == 2 for q in range(self._p_q + 1)}
		self.z_q = {(i, q): pulp.LpVariable(f"Q_zq({i},{q})_ASGM{asgml.ind}", 0, cat=pulp.LpContinuous)
					for (i, info) in self.descriptors_info.items() if len(info) == 2 for q in range(self._p_q + 1)}
		self.z_tmp = {i: pulp.LpVariable(f"Q_z_tmp({i})_ASGM{asgml.ind}", 0, cat=pulp.LpContinuous)
					for (i, info) in self.descriptors_info.items() if len(info) == 2}

		for i in range(1, inv_num_K + 1):
			info = self.descriptors_info[i]
			if len(info) == 1:
				# linear descriptor
				fv_name = info[0]
				MILP += x_hat[i] == fv_norm[self.descriptors_ind[fv_name]], f"connect_ASGM{asgml.ind}_{i}"
			else:
				# quadratic descriptor
				fv_name1 = info[0][0]
				fv_name2 = info[1][0]
				q_ind_1 = self.descriptors_ind[fv_name1]
				q_ind_2 = self.descriptors_ind[fv_name2]
				
				MILP += (2**(self._p + 1) - 1) * fv_norm[q_ind_1] >= \
						pulp.lpSum([2**q * self.x_q[(i, q)] for q in range(self._p_q + 1)]) - 1, f"connect_Q-{i}-1-1_ASGM{asgml.ind}"
				MILP += (2**(self._p + 1) - 1) * fv_norm[q_ind_1] <= \
						pulp.lpSum([2**q * self.x_q[(i, q)] for q in range(self._p_q + 1)]), f"connect_Q-{i}-1-2_ASGM{asgml.ind}"
				for q in range(self._p_q + 1):
					MILP += self.z_q[(i, q)] <= self._M * self.x_q[(i, q)], f"connect_Q-{i}-2-1-{q}_ASGM{asgml.ind}"
					MILP += self.z_q[(i, q)] >= fv_norm[q_ind_2] - self._M * (1 - self.x_q[(i, q)]), f"connect_Q-{i}-2-2-{q}-1_ASGM{asgml.ind}"
					MILP += self.z_q[(i, q)] <= fv_norm[q_ind_2] + self._M * (1 - self.x_q[(i, q)]), f"connect_Q-{i}-2-2-{q}-2_ASGM{asgml.ind}"
				MILP += self.z_tmp[i] * (2**(self._p + 1) - 1) == \
					pulp.lpSum([2**q * self.z_q[(i, q)] for q in range(self._p_q + 1)]), f"connect_Q-{i}-3_ASGM{asgml.ind}"

				if info[0][1] == -1 and info[1][1] == -1:
					MILP += x_hat[i] == 1 - fv_norm[q_ind_1] - fv_norm[q_ind_2] + self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"
				elif info[0][1] == -1:
					MILP += x_hat[i] == fv_norm[q_ind_2] - self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"
				elif info[1][1] == -1:
					MILP += x_hat[i] == fv_norm[q_ind_1] - self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"
				else:
					MILP += x_hat[i] == self.z_tmp[i], f"connect_Q-{i}-4_ASGM{asgml.ind}"

		return MILP				

	def add_constraints_lb_ub(self, MILP, output, asgml, y_min, y_max):
		### normalize the target value
		target_value_lb = (asgml.lb - y_min) / (y_max - y_min)
		target_value_ub = (asgml.ub - y_min) / (y_max - y_min)
		
		MILP += output >= target_value_lb, f"LB_ASGM{asgml.ind}"
		MILP += output <= target_value_ub, f"UB_ASGM{asgml.ind}"

		return MILP

	def prepare_quadratic_desc(self):
		self.descriptors_info = dict()      # a dict stores all information of descriptors, including linear and quadratic
		self.descriptors_ind = dict()      # a dict stores the index of the linear descriptors
		self.descriptors = dict()          # a dict stores the completer linear descriptors
		self.descriptors_list = list()     # a list to store names of linear descritpors
		self.forbidden_node = list()

		for i in range(1, self.inv.num_K + 1):
			fv_name = self.inv.fv_list[i - 1]
			if "*" not in fv_name:
				# linear descriptor
				self.descriptors_list.append(fv_name)
				to_forbid = True
				for sg in self.sg:
					self.descriptors[(sg.name, fv_name)], flag = self.get_fv_variable(fv_name, sg)
					if flag != -1:
						to_forbid = False
				self.descriptors_info[i] = [fv_name]
				self.descriptors_ind[fv_name] = len(self.descriptors_list)
				if to_forbid:
					self.forbidden_node.append(fv_name)
			elif "*" in fv_name:
				# quadratic descriptor
				sp_loc1 = fv_name.find("*")
				fv_name1 = fv_name[:sp_loc1]
				fv_name2 = fv_name[sp_loc1 + 1:]
				if "1-" in fv_name1:
					minus_1 = -1
					fv_name1 = fv_name1[3:-1]
				else:
					minus_1 = 1
				if "1-" in fv_name2:
					minus_2 = -1
					fv_name2 = fv_name2[3:-1]
				else:
					minus_2 = 1
				
				if fv_name1 not in self.descriptors:
					self.descriptors_list.append(fv_name1)
					to_forbid = True
					for sg in self.sg:
						self.descriptors[(sg.name, fv_name1)], flag = self.get_fv_variable(fv_name1, sg)
						if flag != -1:
							to_forbid = False
					self.descriptors_ind[fv_name1] = len(self.descriptors_list)
					if to_forbid:
						self.forbidden_node.append(fv_name1)
				if fv_name2 not in self.descriptors:
					self.descriptors_list.append(fv_name2)
					to_forbid = True
					for sg in self.sg:
						self.descriptors[(sg.name, fv_name2)], flag = self.get_fv_variable(fv_name2, sg)
						if flag != -1:
							to_forbid = False
					self.descriptors_ind[fv_name2] = len(self.descriptors_list)
					if to_forbid:
						self.forbidden_node.append(fv_name2)
				self.descriptors_info[i] = [(fv_name1, minus_1), (fv_name2, minus_2)]

	### connect C1 and C2
	def connect_C1_C2(self, MILP):
		self.inv.set_up_variables(self.asgml)
		self.output = self.inv.output[self.asgml.ind]
		self.inv.set_up_x_hat_variables(self.asgml)
		for sg in self.sg:
			self.inv.prepare_fringe_tree_index_set((self.asgml.ind, sg.name), sg.strF)
		self.prepare_quadratic_desc()
		self.num_K = len(self.descriptors_list)
		self.set_up_variables_norm(self.asgml, self.num_K)
		
		MILP = self.add_constraints_norm(MILP, self.sg, self.fv_norm, self.asgml, self.num_K, self.descriptors_list, self.descriptors, self.inv.max_dcp, self.inv.min_dcp, self.forbidden_node, self.inv.std_eps)
		MILP = self.inv.add_constraints_C1(MILP, self.asgml)
		MILP = self.connect_fv_norm_x_hat(MILP, self.asgml, self.fv_norm, self.inv.num_K, self.inv.x_hat[self.asgml.ind])
		MILP = self.add_constraints_lb_ub(MILP, self.inv.output[self.asgml.ind], self.asgml, self.inv.y_min, self.inv.y_max)
		
		return MILP

	def _inspection_pre(self, milp_sdf_filename_list, test_prefix):
		fv_test = np.zeros(len(self.descriptors_list), dtype=float)
		for (sg_name, milp_sdf_filename) in milp_sdf_filename_list.items():
			cmd = [self.inv.fv_gen_path, self.inv.sdf_filename, f"test_tmp_ASGM{self.asgml.ind}", milp_sdf_filename, test_prefix]
			run_with_retry(cmd)

			test_desc_filename = f"{test_prefix}_desc.csv"

			df = pd.read_csv(test_desc_filename, sep=",", index_col=0)
			idx = df.index[0]

			fv_test_tmp = []
			for fv_name in self.descriptors_list:
				if not fv_name.startswith("F"):
					fv_test_tmp.append(float(df.at[idx, fv_name]))
				else:
					# fringe configuration
					ind_fv_name = fv_name.find('_', 3)
					match = None
					for fv_name_test in df.columns:
						ind2 = fv_name_test.find('_', 3)
						if ind2 != -1 and ind_fv_name != -1 and fv_name[ind_fv_name:] == fv_name_test[ind2:]:
							match = fv_name_test
							break
					fv_test_tmp.append(float(df.at[idx, match]) if match else 0.0)
			fv_test += np.asarray(fv_test_tmp, dtype=float) * self.ratio[sg_name]

		# normalization
		fv_test_norm = []
		for (fv_name, fv_desc) in zip(self.descriptors_list, list(fv_test)):
			match = None
			if not fv_name.startswith("F"):
				match = fv_name
			else:
				# fringe configuration
				ind_fv_name = fv_name.find('_', 3)
				match = None
				for fv_name_test in df.columns:
					ind2 = fv_name_test.find('_', 3)
					if ind2 != -1 and ind_fv_name != -1 and fv_name[ind_fv_name:] == fv_name_test[ind2:]:
						match = fv_name_test
						break
			if match:
				fv_max = self.inv.max_dcp[match]
				fv_min = self.inv.min_dcp[match]
				if abs(fv_max - fv_min) < 1e-6:
					fv_test_norm.append(0.0)
				else:
					fv_test_norm.append((fv_desc - fv_min) / (fv_max - fv_min))
			else:
				fv_test_norm.append(0.0)

		return list(fv_test_norm)

	def inspection(self, sdf_filename_prefix, test_prefix):
		sdf_filename = {sg.name: f"{sdf_filename_prefix}_{sg.name}.sdf" for sg in self.sg}
		fv_test_linear = self._inspection_pre(sdf_filename, test_prefix)
		fv_test = []
		for i in range(1, self.inv.num_K + 1):
			info = self.descriptors_info[i]
			if len(info) == 1:
				# linear descriptor
				fv_name = info[0]
				fv_test.append(fv_test_linear[self.descriptors_ind[fv_name] - 1])
			else:
				# quadratic descriptor
				fv_name1 = info[0][0]
				fv_name2 = info[1][0]
				q_ind_1 = self.descriptors_ind[fv_name1]
				q_ind_2 = self.descriptors_ind[fv_name2]
				_x = fv_test_linear[q_ind_1 - 1]
				_y = fv_test_linear[q_ind_2 - 1]
				if info[0][1] == -1:
					_x = 1 - _x
				if info[1][1] == -1:
					_y = 1 - _y
				fv_test.append(_x * _y)

		insp_value = self.inv.inspection(self.asgml, fv_test, std_eps=1/(2**(self._p + 1) - 1) + self.inv.std_eps)
		insp_value_scaled = insp_value * (self.inv.y_max - self.inv.y_min) + self.inv.y_min

		return insp_value, insp_value_scaled

	def get_output(self):
		y_MILP = self.output.value()
		y_scale = y_MILP * (self.inv.y_max - self.inv.y_min) + self.inv.y_min
		
		return y_MILP, y_scale

def create_assignment(asgml, sg, inv):
	_t = asgml.type.upper()

	if _t == "NORMAL":
		return AssignmentNormal(asgml=asgml, sg=sg, inv=inv)
	elif _t == "QUAD":
		return AssignmentQuad(asgml=asgml, sg=sg, inv=inv)
	elif _t == "GNN":
		return AssignmentGNN(asgml=asgml, sg=sg, inv=inv)
	elif _t == "MV":
		return AssignmentMixingVector(asgml=asgml, sg=sg, inv=inv)
	else:
		raise ValueError(f"Unknown assignment type: {asgml.type} !")

def check_assignment(asgm):
	if asgm.inv._use_quadratic and not asgm._support_quadratic:
		raise ValueError(
			f"Assignment {asgm.asgml.ind} does not support quadratic descriptors "
			f"while a learning model ({asgm.inv.name}) using quadratic descriptors is assigned!"
		)
	if asgm.type == "GNN" and not asgm.inv.type == "INV_GNN":
		raise ValueError(
			f"Assignment {asgm.asgml.ind} can only be used for GNN models!"
		)


