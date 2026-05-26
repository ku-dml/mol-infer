import sys, re
from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, ClassVar
# import pulp_modified as pulp
import pandas as pd
from .seed_graph_base import *
from src.utils import _read_clean_lines, _pop_int, _pop_ints, _pop_strs

_prevent_double_bound = True

@dataclass(frozen=False, kw_only=True)
class SeedGraph_2LMM_Polymer(Rdkit2DPostprocessMixin, SeedGraphBase):
	### a class to save the information about seed graphs
	### 2LMM model for polymer
	name: str
	_type: ClassVar[Literal["SEED_GRAPH_2LMM_POLYMER"]] = "SEED_GRAPH_2LMM_POLYMER"

	def __post_init__(self):
		self.MAX_BOND = 3
		self.MAX_VAL = 4

	def read_seed_graph(self, filename):
		Q = _read_clean_lines(filename)

		### read V_C ###
		self.num_V_C = _pop_int(Q)
		self.V_C = tuple(range(1, self.num_V_C + 1))
		
		### read E_C ###
		self.num_E_C = _pop_int(Q)
		self.E_C = {}
		for e in range(self.num_E_C):
			arr = _pop_ints(Q)
			self.E_C[arr[0]] = (arr[0], arr[1], arr[2])   # Add arr[0] to distinguish two edges with same starting and ending vertices, by Zhu
			
		### read ell_LB and ell_UB ###
		self.ell_LB = {}
		self.ell_UB = {}
		for e in range(self.num_E_C):
			arr = _pop_ints(Q)
			self.ell_LB[arr[0]] = arr[1]
			self.ell_UB[arr[0]] = arr[2]

		### compute E_ge_two, E_ge_one, E_zero_one, E_equal_one ###
		self.E_ge_two = []
		self.E_ge_one = []
		self.E_zero_one = []
		self.E_equal_one = []
		self.I_ge_two = []
		self.I_ge_one = []
		self.I_zero_one = []
		self.I_equal_one = []
		for e in self.E_C:
			if self.ell_LB[e] >= 2:
				self.E_ge_two.append(self.E_C[e])
				self.I_ge_two.append(e)
			elif self.ell_LB[e] == 1 and self.ell_UB[e] >= 2:
				self.E_ge_one.append(self.E_C[e])
				self.I_ge_one.append(e)
			elif self.ell_LB[e] == 0 and self.ell_UB[e] == 1:
				self.E_zero_one.append(self.E_C[e])
				self.I_zero_one.append(e)
			elif self.ell_LB[e] == 1 and self.ell_UB[e] == 1:
				self.E_equal_one.append(self.E_C[e])
				self.I_equal_one.append(e)
			else:
				raise ValueError(f"Strange enge bounds: eid={e}, lb={self.ell_LB[e]}, ub={self.ell_UB[e]}.")
				
		### read n_LB_int and n_UB_int ###
		self.n_LB_int = _pop_int(Q)
		self.n_UB_int = _pop_int(Q)

		### read n_LB_lnk and n_UB_lnk ###  #cao
		self.n_LB_lnk = _pop_int(Q)
		self.n_UB_lnk = _pop_int(Q)

		# read n_LB and n_star
		self.n_LB = _pop_int(Q)
		self.n_star = _pop_int(Q)

		# read rho
		self.rho = _pop_int(Q)

		### read ch_LB and ch_UB ###
		self.ch_LB = {}
		self.ch_UB = {}
		for v in range(self.num_V_C):
			arr = _pop_ints(Q)
			self.ch_LB[arr[0]] = arr[1]
			self.ch_UB[arr[0]] = arr[2]
		for e in range(len(self.E_ge_two + self.E_ge_one)):
			arr = _pop_ints(Q)
			self.ch_LB[self.E_C[arr[0]]] = arr[1]
			self.ch_UB[self.E_C[arr[0]]] = arr[2]

		### read bl_LB and bl_UB ###
		self.bl_LB = {}
		self.bl_UB = {}
		for v in range(self.num_V_C):
			arr = _pop_ints(Q)
			self.bl_LB[arr[0]] = arr[1]
			self.bl_UB[arr[0]] = arr[2]
		for e in range(len(self.E_ge_two + self.E_ge_one)):
			arr = _pop_ints(Q)
			self.bl_LB[self.E_C[arr[0]]] = arr[1]
			self.bl_UB[self.E_C[arr[0]]] = arr[2]

		# read Lambda
		self.Lambda = _pop_strs(Q)

		# read Lambda_dg_int
		num = _pop_int(Q)
		self.Lambda_dg_int = list()
		for i in range(num):
			arr = _pop_strs(Q)
			self.Lambda_dg_int.append((arr[0], int(arr[1])))

		# read Gamma_int_ac
		num = _pop_int(Q)
		self.Gamma_int_ac = list()
		nu_int = list()
		for i in range(num):
			arr = _pop_strs(Q)
			tmp_1 = (arr[0], arr[1], int(arr[2]))
			tmp_2 = (arr[1], arr[0], int(arr[2]))
			nu_int.append(tmp_1)
			if tmp_1 not in self.Gamma_int_ac:
				self.Gamma_int_ac.append(tmp_1)
			if tmp_2 not in self.Gamma_int_ac:
				self.Gamma_int_ac.append(tmp_2)

		# read Gamma_int
		num = _pop_int(Q)
		self.Gamma_int = list()
		gam_int = list()
		for i in range(num):
			arr = _pop_strs(Q)
			tmp_1 = ((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4]))
			tmp_2 = ((arr[2], int(arr[3])), (arr[0], int(arr[1])), int(arr[4]))
			gam_int.append(tmp_1)
			if tmp_1 not in self.Gamma_int:
				self.Gamma_int.append(tmp_1)
			if tmp_2 not in self.Gamma_int:
				self.Gamma_int.append(tmp_2)

		# read Gamma_ac_lnk #cao
		num = _pop_int(Q)
		self.Gamma_ac_lnk = list()
		nu_lnk = list()
		for i in range(num):
			arr = _pop_strs(Q)
			tmp_1 = (arr[0], arr[1], int(arr[2]))
			tmp_2 = (arr[1], arr[0], int(arr[2]))
			nu_lnk.append(tmp_1)
			if tmp_1 not in self.Gamma_ac_lnk:
				self.Gamma_ac_lnk.append(tmp_1)
			if tmp_2 not in self.Gamma_ac_lnk:
				self.Gamma_ac_lnk.append(tmp_2)

		# read Gamma_lnk #cao
		num = _pop_int(Q)
		self.Gamma_lnk = list()
		gam_lnk = list()
		for i in range(num):
			arr = _pop_strs(Q)
			tmp_1 = ((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4]))
			tmp_2 = ((arr[2], int(arr[3])), (arr[0], int(arr[1])), int(arr[4]))
			gam_lnk.append(tmp_1)
			if tmp_1 not in self.Gamma_lnk:
				self.Gamma_lnk.append(tmp_1)
			if tmp_2 not in self.Gamma_lnk:
				self.Gamma_lnk.append(tmp_2)

		# read Lambda_star
		self.Lambda_star = {i: set() for i in range(1, self.num_V_C + 1)}
		for i in range(1, self.num_V_C + 1):
			arr = _pop_strs(Q)
			ind = int(arr[0])
			arr.pop(0)
			for a in arr:
				self.Lambda_star[ind].add(a)
		
		# read na_LB and na_UB
		num = _pop_int(Q)
		self.na_LB = {}
		self.na_UB = {}
		for i in range(num):
			arr = _pop_strs(Q)
			self.na_LB[arr[0]] = int(arr[1])
			self.na_UB[arr[0]] = int(arr[2])

		self.Lambda_int = list()
		# read na_LB_int and na_UB_int
		num = _pop_int(Q)
		self.na_LB_int = {}
		self.na_UB_int = {}
		for i in range(num):
			arr = _pop_strs(Q)
			self.na_LB_int[arr[0]] = int(arr[1])
			self.na_UB_int[arr[0]] = int(arr[2])
			self.Lambda_int.append(arr[0])

		# read ns_LB_int and ns_UB_int
		num = _pop_int(Q)
		self.ns_LB_int = {}
		self.ns_UB_int = {}
		for i in range(num):
			arr = _pop_strs(Q)
			self.ns_LB_int[(arr[0], int(arr[1]))] = int(arr[2])
			self.ns_UB_int[(arr[0], int(arr[1]))] = int(arr[3])

		# read ns_LB_cnt and ns_UB_cnt
		num = _pop_int(Q)
		self.ns_LB_cnt = {}
		self.ns_UB_cnt = {}
		for i in range(num):
			arr = _pop_strs(Q)
			self.ns_LB_cnt[(arr[0], int(arr[1]))] = int(arr[2])
			self.ns_UB_cnt[(arr[0], int(arr[1]))] = int(arr[3])

		# read ac_LB_int and ac_UB_int
		num = _pop_int(Q)
		self.ac_LB_int = {}
		self.ac_UB_int = {}
		for i in range(num):
			arr = _pop_strs(Q)
			a1, a2, m = nu_int[int(arr[0]) - 1]
			self.ac_LB_int[(a1, a2, m)] = int(arr[1])
			self.ac_LB_int[(a2, a1, m)] = int(arr[1])
			self.ac_UB_int[(a1, a2, m)] = int(arr[2])
			self.ac_UB_int[(a2, a1, m)] = int(arr[2])

		# read ec_LB_int and ec_UB_int
		num = _pop_int(Q)
		self.ec_LB_int = {}
		self.ec_UB_int = {}
		for i in range(num):
			arr = _pop_strs(Q)
			a1, a2, m = gam_int[int(arr[0]) - 1]
			self.ec_LB_int[(a1, a2, m)] = int(arr[1])
			self.ec_LB_int[(a2, a1, m)] = int(arr[1])
			self.ec_UB_int[(a1, a2, m)] = int(arr[2])
			self.ec_UB_int[(a2, a1, m)] = int(arr[2])

		# read ac_LB_lnk and ac_UB_lnk #cao
		num = _pop_int(Q)
		self.ac_LB_lnk = {}
		self.ac_UB_lnk = {}
		for i in range(num):
			arr = _pop_strs(Q)
			a1, a2, m = nu_lnk[int(arr[0]) - 1]
			self.ac_LB_lnk[(a1, a2, m)] = int(arr[1])
			self.ac_LB_lnk[(a2, a1, m)] = int(arr[1])
			self.ac_UB_lnk[(a1, a2, m)] = int(arr[2])
			self.ac_UB_lnk[(a2, a1, m)] = int(arr[2])

		# read ec_LB_lnk and ec_UB_lnk #cao
		num = _pop_int(Q)
		self.ec_LB_lnk = {}
		self.ec_UB_lnk = {}
		for i in range(num):
			arr = _pop_strs(Q)
			a1, a2, m = gam_lnk[int(arr[0]) - 1]
			self.ec_LB_lnk[(a1, a2, m)] = int(arr[1])
			self.ec_LB_lnk[(a2, a1, m)] = int(arr[1])
			self.ec_UB_lnk[(a1, a2, m)] = int(arr[2])
			self.ec_UB_lnk[(a2, a1, m)] = int(arr[2])

		# read E_C_lnk #cao
		self.num_E_C_lnk = _pop_int(Q)
		self.E_C_lnk = _pop_ints(Q)

		# read bd2_LB and bd2_UB
		self.bd2_LB = {}
		self.bd2_UB = {}
		for e in range(len(self.E_C)):
			arr = _pop_ints(Q)
			self.bd2_LB[self.E_C[arr[0]]] = arr[1]
			self.bd2_UB[self.E_C[arr[0]]] = arr[2]

		# read bd3_LB and bd3_UB
		self.bd3_LB = {}
		self.bd3_UB = {}
		for e in range(len(self.E_C)):
			arr = _pop_ints(Q)
			self.bd3_LB[self.E_C[arr[0]]] = arr[1]
			self.bd3_UB[self.E_C[arr[0]]] = arr[2]

		# read ac_LB_lf and ac_UB_lf
		num = _pop_int(Q)
		self.ac_LB_lf = dict()
		self.ac_UB_lf = dict()
		for e in range(num):
			arr = _pop_strs(Q)
			self.ac_LB_lf[(arr[0], arr[1], int(arr[2]))] = int(arr[3])
			self.ac_UB_lf[(arr[0], arr[1], int(arr[2]))] = int(arr[4])
		arr = _pop_ints(Q)
		self.ac_LB_lf_common = arr[0]
		self.ac_UB_lf_common = arr[1]   
			
		####################################
		# # Undefined constants for instances but used in MILP
		self.r_GC = self.num_E_C - (self.num_V_C - 1)
		self.dg_LB = [0, 0, 0, 0, 0]
		self.dg_UB = [self.n_star, self.n_star, self.n_star, self.n_star, self.n_star]

	def read_fringe_trees(self, fringe_filename):
		self.set_F = list()
		self.strF = dict()

		self.fc_LB = dict()
		self.fc_UB = dict()
		
		with open(fringe_filename,'r') as f:
			lines = f.readlines()
			for line in lines:
				if len(line.split(",")) < 4:
					continue
				ind = int(line.split(",")[0])
				str1 = line.split(",")[1]
				str2 = line.split(",")[2]
				str3 = line.split(",")[3].replace('\n', '')
				if len(line.split(",")) > 4:
					LB_tmp = line.split(",")[4].replace('\n', '')
					LB_tmp = LB_tmp.replace(' ', '')
					self.fc_LB[ind] = int(LB_tmp)
					UB_tmp = line.split(",")[5].replace('\n', '')
					UB_tmp = UB_tmp.replace(' ', '')
					self.fc_UB[ind] = int(UB_tmp)
				else:
					self.fc_LB[ind] = 0
					# self.fc_UB[ind] = 10
					self.fc_UB[ind] = self.num_V_C
				
				psi = FringeTree()
				psi.read_from_seq(str1, str2, str3)
				psi.index = ind
				flag = True
				for (a, d) in psi.vertex:
					if a not in self.Lambda:
						flag = False
						break
				if flag:
					self.strF[ind] = (str1, str2, str3)
					self.set_F.append(psi)

		self.Lambda_ex = list()
		for psi in self.set_F:
			for (a, d) in psi.vertex[1:]:
				if a not in self.Lambda_ex and a in self.Lambda:
					self.Lambda_ex.append(a)

		# return set_F, Lambda_ex, strF, fc_LB, fc_UB

	# def prepare_lambda_pre_do(self, ann_training_data_filename, Lambda, Lambda_dg_int, Lambda_int, Gamma_int, Gamma_int_ac):
	# 	''' A function to prepare some variables about the descriptors 
	# 		from the feature vector '''

	# 	descriptors = dict()      # get variable from fv file
	# 	fv = pd.read_csv(ann_training_data_filename, sep=",")
	# 	num_fv = len(list(fv.columns))

	# 	fv_Lambda = list()

	# 	for i, fv_name in enumerate(list(fv.columns)):
	# 		try:
	# 			if fv_name[0] == "n" and fv_name[3] == "i":
	# 				ele = fv_name[6:] # may need change, Zhu, 0109
	# 				if ele not in fv_Lambda:
	# 					fv_Lambda.append(ele)
	# 			elif fv_name[0] == "n" and fv_name[3] == "e":
	# 				ele = fv_name[6:]
	# 				if ele not in fv_Lambda:
	# 					fv_Lambda.append(ele)
	# 		except:
	# 			pass

	# 	new_Lambda = list()
	# 	new_Lambda_dg_int = list()
	# 	new_Lambda_int = list()
	# 	new_Gamma_int = list()
	# 	new_Gamma_int_ac = list()

	# 	new_Lambda.append("H1")

	# 	for ele in Lambda:
	# 		if ele in fv_Lambda:
	# 			new_Lambda.append(ele)
	# 	for (ele, d) in Lambda_dg_int:
	# 		if ele in fv_Lambda:
	# 			new_Lambda_dg_int.append((ele, d))
	# 	for ele in Lambda_int:
	# 		if ele in fv_Lambda:
	# 			new_Lambda_int.append(ele)
	# 	for ((a, d1), (b, d2), m) in Gamma_int:
	# 		if a in fv_Lambda and b in fv_Lambda:
	# 			new_Gamma_int.append(((a, d1), (b, d2), m))
	# 	for (a, b, m) in Gamma_int_ac:
	# 		if a in fv_Lambda and b in fv_Lambda:
	# 			new_Gamma_int_ac.append((a, b, m))

	# 	return new_Lambda, new_Lambda_dg_int, new_Lambda_int, new_Gamma_int, new_Gamma_int_ac

	# def prepare_lambda_pre(self, desc_filename):
	# 	self.Lambda, self.Lambda_dg_int, self.Lambda_int, self.Gamma_int, self.Gamma_int_ac = \
	# 		self.prepare_lambda_pre_do(desc_filename, self.Lambda, self.Lambda_dg_int, self.Lambda_int, self.Gamma_int, self.Gamma_int_ac)

	def prepare_Lambda(self, set_Lambda, Lambda):
		new_set_Lambda = list()
		for atom in set_Lambda:
			if atom.symbol in Lambda:
				new_set_Lambda.append(atom)
			else:
				for atom2 in Lambda:
					if re.sub(r'[0-9]', '', atom2) == atom.symbol:
						name = re.sub(r'[A-Z]', '', atom2)
						name = re.sub(r'[a-z]', '', name)

						new_set_Lambda.append(CG_element(atom2, int(name), atom.mass, atom.ele_neg, atom.ion_energy, atom.ele_affinity))

		return new_set_Lambda

	def prepare_Lambda_dg(self, 
		set_Lambda, Lambda, Lambda_int, Lambda_dg_int, Lambda_ex, Gamma_int, Gamma_int_ac
	):
		#   A function to prepare the set Lambda_dg and Gamma
		set_Lambda_dg = list()
		for atom in set_Lambda:
			if atom.symbol != "H" and atom.symbol != "H1":
				for i in range(1, atom.valence + 4):
					if (atom.symbol, i) in Lambda_dg_int:
						set_Lambda_dg.append((atom, i))

		Lambda_dg = [(ele.symbol, i) for (ele, i) in set_Lambda_dg]  # The set Lambda_dg
		epsilon = "e"
		epsilon_dg = ("e", 0)
		Code_Lambda_dg = {ele_dg : i + 1 for i, ele_dg in enumerate(Lambda_dg)}
		Code_Lambda_dg[epsilon_dg] = 0
		Code_Lambda = {ele : i + 1 for i, ele in enumerate(Lambda)}
		Code_Lambda[epsilon] = 0

		Code_Lambda_int = {ele : i + 1 for i, ele in enumerate(Lambda_int)}
		Code_Lambda_int[epsilon] = 0

		Code_Lambda_ex = {ele : i + 1 for i, ele in enumerate(Lambda_ex)}
		Code_Lambda_ex[epsilon] = 0

		MAX_CODE = len(Lambda)
		MAX_CODE_dg = len(Lambda_dg)

		MAX_CODE_int = len(Lambda_int)
		MAX_CODE_ex = len(Lambda_ex)

		Code_Gamma_ec_int = {gamma : i + 1 for i, gamma in enumerate(Gamma_int)}
		Code_Gamma_ac_int = {gamma : i + 1 for i, gamma in enumerate(Gamma_int_ac)}

		val = {ele.symbol: ele.valence for ele in set_Lambda}
		val_temp = {(ele.symbol, i): ele.valence for (ele, i) in set_Lambda_dg}
		val.update(val_temp)
		mass = {ele.symbol: ele.mass for ele in set_Lambda}
		ele_neg = {ele.symbol: ele.ele_neg for ele in set_Lambda}
		ion_energy = {ele.symbol: ele.ion_energy for ele in set_Lambda}
		ele_affinity = {ele.symbol: ele.ele_affinity for ele in set_Lambda}

		return set_Lambda_dg, Lambda_dg, epsilon, epsilon_dg, Code_Lambda_dg, \
			Code_Lambda, Code_Lambda_int, Code_Lambda_ex, \
			MAX_CODE, MAX_CODE_dg, MAX_CODE_int, MAX_CODE_ex, \
			Code_Gamma_ec_int, Code_Gamma_ac_int, val, mass, ele_neg, ion_energy, ele_affinity

	def prepare_Gamma_ac(self, 
		Gamma_int, Gamma_int_ac, Code_Lambda, Code_Lambda_dg, 
		Gamma_lnk, Gamma_ac_lnk, ns_LB_cnt, ns_UB_cnt
	):
		# A function to prepare the sets Gamma_int and Gamma_int_ac
		Gamma_int_less = list()
		Gamma_int_equal = list()
		Gamma_int_great = list()

		Gamma_lnk_less = list()
		Gamma_lnk_equal = list()
		Gamma_lnk_great = list()

		for ((a1, d1), (a2, d2), m) in Gamma_int:
			if Code_Lambda_dg[(a1, d1)] < Code_Lambda_dg[(a2, d2)]:
				Gamma_int_less.append(((a1, d1), (a2, d2), m))
			elif Code_Lambda_dg[(a1, d1)] == Code_Lambda_dg[(a2, d2)]:
				Gamma_int_equal.append(((a1, d1), (a2, d2), m))
			elif Code_Lambda_dg[(a1, d1)] > Code_Lambda_dg[(a2, d2)]:
				Gamma_int_great.append(((a1, d1), (a2, d2), m))

		for ((a1, d1), (a2, d2), m) in Gamma_lnk:
			if Code_Lambda_dg[(a1, d1)] < Code_Lambda_dg[(a2, d2)]:
				Gamma_lnk_less.append(((a1, d1), (a2, d2), m))
			elif Code_Lambda_dg[(a1, d1)] == Code_Lambda_dg[(a2, d2)]:
				Gamma_lnk_equal.append(((a1, d1), (a2, d2), m))
			elif Code_Lambda_dg[(a1, d1)] > Code_Lambda_dg[(a2, d2)]:
				Gamma_lnk_great.append(((a1, d1), (a2, d2), m))

		Gamma_int_ac_less = list()
		Gamma_int_ac_equal = list()
		Gamma_int_ac_great = list()

		Gamma_ac_lnk_less = list()  #IDO
		Gamma_ac_lnk_equal = list()
		Gamma_ac_lnk_great = list()
	 
		for ((a1, d1), (a2, d2), m) in Gamma_int:
			if Code_Lambda[a1] < Code_Lambda[a2]:
				if (a1, a2, m) not in Gamma_int_ac_less:
					Gamma_int_ac_less.append((a1, a2, m))
				if (a2, a1, m) not in Gamma_int_ac_great:
					Gamma_int_ac_great.append((a2, a1, m))
			elif Code_Lambda[a1] ==  Code_Lambda[a2]:
				if (a1, a2, m) not in Gamma_int_ac_equal:
					Gamma_int_ac_equal.append((a1, a2, m))
			elif Code_Lambda[a1] > Code_Lambda[a2]:
				if (a1, a2, m) not in Gamma_int_ac_great:
					Gamma_int_ac_great.append((a1, a2, m))  
				if (a2, a1, m) not in Gamma_int_ac_less:
					Gamma_int_ac_less.append((a2, a1, m))

		for ((a1, d1), (a2, d2), m) in Gamma_lnk:
			if Code_Lambda[a1] < Code_Lambda[a2]:
				if (a1, a2, m) not in Gamma_ac_lnk_less:
					Gamma_ac_lnk_less.append((a1, a2, m))
				if (a2, a1, m) not in Gamma_ac_lnk_great:
					Gamma_ac_lnk_great.append((a2, a1, m))
			elif Code_Lambda[a1] ==  Code_Lambda[a2]:
				if (a1, a2, m) not in Gamma_ac_lnk_equal:
					Gamma_ac_lnk_equal.append((a1, a2, m))
			elif Code_Lambda[a1] > Code_Lambda[a2]:
				if (a1, a2, m) not in Gamma_ac_lnk_great:
					Gamma_ac_lnk_great.append((a1, a2, m))  
				if (a2, a1, m) not in Gamma_ac_lnk_less:
					Gamma_ac_lnk_less.append((a2, a1, m))

		Gamma_tilde_ac_C = Gamma_int_ac_equal + Gamma_int_ac_great + Gamma_int_ac_less
		Gamma_tilde_ac_T = Gamma_int_ac_equal + Gamma_int_ac_great + Gamma_int_ac_less
		Gamma_tilde_ac_CT = Gamma_int_ac_equal + Gamma_int_ac_great + Gamma_int_ac_less
		Gamma_tilde_ac_TC = Gamma_int_ac_equal + Gamma_int_ac_great + Gamma_int_ac_less
		Gamma_tilde_ac_F = Gamma_int_ac_equal + Gamma_int_ac_great + Gamma_int_ac_less
		Gamma_tilde_ac_CF = Gamma_int_ac_equal + Gamma_int_ac_great + Gamma_int_ac_less
		Gamma_tilde_ac_TF = Gamma_int_ac_equal + Gamma_int_ac_great + Gamma_int_ac_less

		Gamma_tilde_ec_C = Gamma_int_equal + Gamma_int_less + Gamma_int_great
		Gamma_tilde_ec_T = Gamma_int_equal + Gamma_int_less + Gamma_int_great
		Gamma_tilde_ec_CT = Gamma_int_equal + Gamma_int_less + Gamma_int_great
		Gamma_tilde_ec_TC = Gamma_int_equal + Gamma_int_less + Gamma_int_great
		Gamma_tilde_ec_F = Gamma_int_equal + Gamma_int_less + Gamma_int_great
		Gamma_tilde_ec_CF = Gamma_int_equal + Gamma_int_less + Gamma_int_great
		Gamma_tilde_ec_TF = Gamma_int_equal + Gamma_int_less + Gamma_int_great

		Gamma_cnt_less = list()
		Gamma_cnt_equal = list()
		Gamma_cnt_great = list()

		for (mu1, mu2, m) in Gamma_lnk_less:
			if m == 1 and ns_LB_cnt[mu1] <= 1 and ns_UB_cnt[mu1] >= 1 and ns_LB_cnt[mu2] <= 1 and ns_UB_cnt[mu2] >= 1:
				Gamma_cnt_less.append((mu1, mu2, m))
		for (mu1, mu2, m) in Gamma_lnk_great:
			if m == 1 and ns_LB_cnt[mu1] <= 1 and ns_UB_cnt[mu1] >= 1 and ns_LB_cnt[mu2] <= 1 and ns_UB_cnt[mu2] >= 1:
				Gamma_cnt_great.append((mu1, mu2, m))
		for (mu1, mu2, m) in Gamma_lnk_equal:
			if m == 1 and ns_UB_cnt[mu1] == 2:
				Gamma_cnt_equal.append((mu1, mu2, m))

		return Gamma_int_less, Gamma_int_equal, Gamma_int_great, \
			Gamma_int_ac_less, Gamma_int_ac_equal, Gamma_int_ac_great, \
			Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, \
			Gamma_tilde_ac_F, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, \
			Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, \
			Gamma_tilde_ec_F, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF, \
			Gamma_lnk_less, Gamma_lnk_equal, Gamma_lnk_great, \
			Gamma_ac_lnk_less, Gamma_ac_lnk_equal, Gamma_ac_lnk_great, \
			Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great

	def prepare_Gamma_lf_ac(self, 
		set_Lambda, ac_LB_lf, ac_UB_lf, ac_LB_lf_common, ac_UB_lf_common
	):
		# A function to prepare the set Gamma_lf_ac
		Gamma_lf_ac = list()
		num = len(set_Lambda)
		for a1 in range(num):
			for a2 in range(num):
				for m in range(1, min(set_Lambda[a1].valence, set_Lambda[a2].valence) + 1):
					if m > 3:
						continue
					ac_lf_tmp = (set_Lambda[a1].symbol, set_Lambda[a2].symbol, m)
					Gamma_lf_ac.append(ac_lf_tmp)
					if ac_lf_tmp not in ac_LB_lf:
						ac_LB_lf[ac_lf_tmp] = ac_LB_lf_common
						ac_UB_lf[ac_lf_tmp] = ac_UB_lf_common

		return Gamma_lf_ac, ac_LB_lf, ac_UB_lf

	def prepare_Lambda_tilde(self, 
		Lambda_dg_int
	):

		Lambda_tilde_dg_C = Lambda_dg_int[:]
		Lambda_tilde_dg_T = Lambda_dg_int[:]
		Lambda_tilde_dg_F = Lambda_dg_int[:]

		return Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F

	def prepare_dataset_for_scheme_graph(self,
		n_star, rho, V_C, E_C, E_C_lnk, E_ge_two, E_ge_one, n_LB_int, n_UB_int, bl_LB, bl_UB, ch_LB, ch_UB, ell_UB,
		I_ge_two, I_ge_one, I_equal_one, I_zero_one, bd2_LB, bd3_LB, val, Lambda_star
	):
		# A function to prepare constants of the seed graph
		t_C = len(V_C)
		t_C_tilde = sum(1 for i in V_C if bl_UB[i] >= 1)
		# print(t_C_tilde)
		t_T = n_UB_int - len(V_C)
		t_F = n_star - n_LB_int

		k_C = len(E_ge_two + E_ge_one)
		k_C_tilde = len(E_ge_two)
		c_F = t_C_tilde + t_T
		m_C = len(E_C)
		head_C = {i : E_C[i][2] for i in range(1, m_C + 1)}
		tail_C = {i : E_C[i][1] for i in range(1, m_C + 1)}
		tail_F = {i : i if i <= t_C_tilde else i - t_C for i in range(1, c_F + 1)}
		E_C_plus = [[] for _ in range(t_C + 1)]
		E_C_minus = [[] for _ in range(t_C + 1)]
		I_ge_one_plus = [[] for _ in range(t_C + 1)]
		I_ge_one_minus = [[] for _ in range(t_C + 1)]
		I_ge_two_plus = [[] for _ in range(t_C + 1)]
		I_ge_two_minus = [[] for _ in range(t_C + 1)]
		I_zero_one_plus = [[] for _ in range(t_C + 1)]
		I_zero_one_minus = [[] for _ in range(t_C + 1)]
		I_equal_one_plus = [[] for _ in range(t_C + 1)]
		I_equal_one_minus = [[] for _ in range(t_C + 1)]
		
		for i in range(1, m_C + 1):
			_, tail, head = E_C[i]
			E_C_plus[tail].append(i)
			E_C_minus[head].append(i)
			if i in I_ge_one:
				I_ge_one_plus[tail].append(i)
				I_ge_one_minus[head].append(i)
			elif i in I_ge_two:
				I_ge_two_plus[tail].append(i)
				I_ge_two_minus[head].append(i)
			elif i in I_zero_one:
				I_zero_one_plus[tail].append(i)
				I_zero_one_minus[head].append(i)
			elif i in I_equal_one:
				I_equal_one_plus[tail].append(i)
				I_equal_one_minus[head].append(i)

		beta_star = {i: 0 for i in V_C}
		for i in V_C:
			beta_star[i] = sum([1 for e in I_ge_two + I_ge_one if E_C[e][1] == i or E_C[e][2] == i]) + \
						   sum([1 + bd2_LB[E_C[e]] + 2 * bd3_LB[E_C[e]] 
								for e in I_equal_one if E_C[e][1] == i or E_C[e][2] == i]) + \
						   bl_LB[i]
		delta_i = {i: 0 for i in V_C}
		for i in V_C:
			# if i in V_C_star:
			#     delta_i[i] = val[alpha_star[i]] - beta_star[i]
			# else:
			#     delta_i[i] = 4 - beta_star[i]
			delta_i[i] = max([val[a] for a in Lambda_star[i]]) - beta_star[i]

		I_lnk = list()
		
		for e in E_C_lnk:
			I_lnk.append(e)
		
		n_lnk_equal_one = len(set(I_lnk) & set(I_equal_one))  #IDO
			
		return t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, c_F, m_C, \
			head_C, tail_C, tail_F, E_C_plus, E_C_minus, \
			I_lnk, I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, \
			I_zero_one_plus, I_zero_one_minus, I_equal_one_plus, I_equal_one_minus, \
			delta_i, n_lnk_equal_one

	def prepare_fringe_tree(self,
		set_F, V_C, t_T, t_F, val, Lambda_int, Lambda_ex, Code_Lambda_int, Gamma_lf_ac
	):
		# A function to prepare constants for fringe trees
		set_F_E = set_F
		set_F_v = {v : set_F for v in V_C}

		F_C = set_F_v
		F_T = {i : set_F_E for i in range(1, t_T + 1)}
		F_F = {i : set_F_E for i in range(1, t_F + 1)}
		n_C = max(len(psi.vertex) - 1 for v in V_C for psi in set_F_v[v])
		n_T = max(len(psi.vertex) - 1 for psi in set_F_E)
		n_F = max(len(psi.vertex) - 1 for psi in set_F_E)

		Code_F = {psi : psi.index for psi in set_F}

		# number of non-root vertices
		n_psi_H = {Code_F[psi] : sum([1 for v in psi.vertex if v[0] != "H1"]) - 1 for psi in set_F}

		# degree of root
		deg_r_H = {Code_F[psi] : sum([1 for v in psi.adj[0] if psi.vertex[v][0] != "H1"]) for psi in set_F}
		deg_r_hyd = {Code_F[psi] : sum([1 for v in psi.adj[0] if psi.vertex[v][0] == "H1"]) for psi in set_F}

		# element of root
		atom_r = {Code_F[psi] : psi.root[0] for psi in set_F}
		alpha_r = {i : Code_Lambda_int[atom_r[i]] for i in atom_r}

		# sum of mul of edges incident to root
		beta_r = {Code_F[psi] : sum(psi.beta[0][v] for v in psi.adj[0]) for psi in set_F}
		beta_r_2 = {Code_F[psi]: sum([1 if psi.beta[0][v] == 2 else 0 for v in psi.adj[0]]) for psi in set_F}

		ht_H = {Code_F[psi] : psi.height for psi in set_F}

		v_ion = {Code_F[psi]: psi.chg[0] for psi in set_F}

		na_alpha_ex = {ele : {Code_F[psi] : 0} for psi in set_F for ele in Lambda_ex}
		for psi in set_F:
			na_ex_tmp = {ele : 0 for ele in Lambda_ex}
			for u, (ele, dep) in enumerate(psi.vertex[1:]):
				if ele in Lambda_ex:
					na_ex_tmp[ele] += 1

			for ele, d in na_alpha_ex.items():
				d[Code_F[psi]] = na_ex_tmp[ele]

		deg_fr = {(Code_F[psi], i): 0 for i in range(1, self.MAX_VAL + 1) for psi in set_F}
		for psi in set_F:
			# deg_tmp = {i: 0 for i in range(len(psi.vertex))}
			for i in range(1, len(psi.vertex)):
				if psi.vertex[i][0] == "H1":
					continue
				tmp = len([1 for j in psi.adj[i] if psi.vertex[j][0] != "H1"])
				# deg_tmp[i] = tmp
				deg_fr[(Code_F[psi], tmp)] += 1

		F_Cp = dict()
		for i in F_C.keys():
			F_Cp_tmp = dict()
			F_Cp_tmp[0] = set()
			F_Cp_tmp[1] = set()
			F_Cp_tmp[2] = set()
			
			for psi in F_C[i]:
				p = ht_H[Code_F[psi]]
				if p not in F_Cp_tmp.keys():
					F_Cp_tmp[p] = {psi}
				else:
					F_Cp_tmp[p].add(psi)
			F_Cp[i] = F_Cp_tmp

		F_Tp = dict()
		for i in F_T.keys():
			F_Tp_tmp = dict()
			F_Tp_tmp[0] = set()
			F_Tp_tmp[1] = set()
			F_Tp_tmp[2] = set()

			for psi in F_T[i]:
				p = ht_H[Code_F[psi]]
				if p not in F_Tp_tmp.keys():
					F_Tp_tmp[p] = {psi}
				else:
					F_Tp_tmp[p].add(psi)
			F_Tp[i] = F_Tp_tmp

		F_Fp = dict()
		for i in F_F.keys():
			F_Fp_tmp = dict()
			F_Fp_tmp[0] = set()
			F_Fp_tmp[1] = set()
			F_Fp_tmp[2] = set()

			for psi in F_F[i]:
				p = ht_H[Code_F[psi]]
				if p not in F_Fp_tmp.keys():
					F_Fp_tmp[p] = {psi}
				else:
					F_Fp_tmp[p].add(psi)
			F_Fp[i] = F_Fp_tmp

		ac_psi_lf = {nu : {Code_F[psi] : 0} for psi in set_F for nu in Gamma_lf_ac}
		for psi in set_F:
			ac_psi_lf_tmp = {nu : 0 for nu in Gamma_lf_ac}
			for i in range(1, len(psi.vertex)):
				if psi.vertex[i][0] == "H1":
					continue
				tmp = len([1 for j in psi.adj[i] if psi.vertex[j][0] != "H1"])
				if tmp == 1:
					symbol1 = psi.vertex[i][0]
					symbol2 = ""
					m = 0
					for j in psi.adj[i]:
						if psi.vertex[j][0] != "H1":
							symbol2 = psi.vertex[j][0]
							m = psi.beta[i][j]
							break
					ac_psi_lf_tmp[(symbol1, symbol2, m)] += 1

			for nu, d in ac_psi_lf.items():
				d[Code_F[psi]] = ac_psi_lf_tmp[nu]

		return Code_F, n_psi_H, deg_r_H, deg_r_hyd, beta_r, atom_r, ht_H, F_C, F_T, F_F, \
				n_C, n_T, n_F, F_Cp, F_Tp, F_Fp, set_F_E, set_F_v, \
				na_alpha_ex, alpha_r, deg_fr, v_ion, ac_psi_lf, beta_r_2

	def prepare_variables_selecting_core(self,
		# Constants
		t_C, k_C_tilde, k_C, t_T, m_C, n_LB_int, n_UB_int, n_LB_lnk, n_UB_lnk, ell_LB, ell_UB, bl_LB, bl_UB
	):
		# Define e_C
		e_C = {i : pulp.LpVariable(f"e_C({i})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
				for i in range(1, m_C + 1)}

		# Define v_T
		v_T = {i : pulp.LpVariable(f"v_T({i})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
				for i in range(1, t_T + 1)}

		# Define e_T    
		e_T = {i : pulp.LpVariable(f"e_T({i})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
				for i in range(1, (t_T + 1) + 1)}

		# Define chi_T
		chi_T = {i : pulp.LpVariable(f"chi_T({i})_SG{self.name}", 0, k_C, cat = pulp.LpInteger)
				for i in range(1, t_T + 1)}

		# Define clr_T
		clr_T = {k : pulp.LpVariable(f"clr_T({k})_SG{self.name}", ell_LB[k] - 1, ell_UB[k] - 1, cat = pulp.LpInteger)
				for k in range(1, k_C + 1)}
		clr_T[0] = pulp.LpVariable(f"clr_T(0)_SG{self.name}", 0, t_T, cat = pulp.LpInteger)

		# Define delta_chi_T
		delta_chi_T = {k : pulp.LpVariable(f"delta_chi_T({k})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
				for k in range(0, k_C + 1)} 
		chi_T_tmp = {(i, k) : pulp.LpVariable(f"chi_T({i},{k})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
				for i in range(1, t_T + 1) for k in range(0, k_C + 1)}
		chi_T.update(chi_T_tmp)

		# Define deg_tilde_C_plus
		deg_tilde_C_plus = {i : pulp.LpVariable(f"deg_tilde_C_plus({i})_SG{self.name}", 0, self.MAX_VAL, cat = pulp.LpInteger)
				for i in range(1, t_C + 1)}

		# Define deg_tilde_C_minus
		deg_tilde_C_minus = {i : pulp.LpVariable(f"deg_tilde_C_minus({i})_SG{self.name}", 0, self.MAX_VAL, cat = pulp.LpInteger)
				for i in range(1, t_C + 1)}

		# Define n_G_lnk  IDO
		n_lnk = pulp.LpVariable(f"n_lnk_SG{self.name}", n_LB_lnk, n_UB_lnk, cat = pulp.LpInteger)

		# Define rank_G
		rank_G = pulp.LpVariable(f"rank_G_SG{self.name}", 0, cat=pulp.LpInteger)

		return e_C, v_T, e_T, chi_T, clr_T, delta_chi_T, deg_tilde_C_plus, deg_tilde_C_minus, n_lnk, rank_G

	def add_constraints_selecting_core(self,
		# Model
		MILP,
		# Constants
		t_C, k_C_tilde, k_C, t_T, m_C, n_LB_lnk, n_UB_lnk, ell_LB, ell_UB, bl_LB, bl_UB,
		I_equal_one, I_equal_one_minus, I_equal_one_plus, I_ge_two, I_ge_one, I_ge_one_minus,
		I_ge_one_plus, I_zero_one, I_zero_one_minus, I_zero_one_plus, n_lnk_equal_one, I_lnk, r_GC,
		# Binary Variables
		e_C, v_T, e_T, delta_chi_T,
		# Integer Variables
		chi_T, clr_T, deg_tilde_C_plus, deg_tilde_C_minus, n_lnk, rank_G
	):
		# -------- Constraint (0) --------
		MILP += rank_G == r_GC - pulp.lpSum(1 - e_C[i] for i in I_zero_one), f"milp-2LMM-(0)-rank_SG{self.name}"

		# -------- Constraint (1) --------
		# Constraint of e_C, for i in I_(=1)
		for i in I_equal_one:
			MILP += e_C[i] == 1, f"milp-2LMH-(1)-{i}_SG{self.name}"

		# -------- Constraint (2) --------
		# Constraint of e_C, clr_T, for i in I_(\ge 2)
		for i in I_ge_two:
			MILP += e_C[i] == 0, f"milp-2LMH-(2)-1-{i}_SG{self.name}"
			MILP += clr_T[i] >= 1, f"milp-2LMH-(2)-2-{i}_SG{self.name}"

		# -------- Constraint (3) --------
		# Constraint of e_C, clr_T, for i in I_(\ge 1)
		for i in I_ge_one:
			MILP += e_C[i] + clr_T[i] >= 1, f"milp-2LMH-(3)-1-{i}_SG{self.name}"
			MILP += clr_T[i] <= t_T * (1 - e_C[i]), f"milp-2LMH-(3)-2-{i}_SG{self.name}"

		# -------- Constraint (4) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(e_C[c] for c in (I_ge_one_minus[i] + I_zero_one_minus[i] + I_equal_one_minus[i])) == \
					deg_tilde_C_minus[i], f"milp-2LMH-(4)-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(e_C[c] for c in (I_ge_one_plus[i] + I_zero_one_plus[i] + I_equal_one_plus[i])) == \
					deg_tilde_C_plus[i], f"milp-2LMH-(4)-2-{i}_SG{self.name}"

		# -------- Constraint (5) --------
		# Constraint of delta_chi_T, for i in [1, t_T]
		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(chi_T[(i, k)]
					for k in range(0, k_C + 1)) == 1, f"milp-2LMH-(5)-2-{i}_SG{self.name}"        
			MILP += pulp.lpSum(k * chi_T[(i, k)]
					for k in range(0, k_C + 1)) == chi_T[i], f"milp-2LMH-(5)-3-{i}_SG{self.name}"
			MILP += chi_T[(i, 0)] == 1 - v_T[i], f"milp-2LMH-(5)-1-{i}_SG{self.name}"

		# -------- Constraint (6) --------
		# Constraint of delta_chi_T, clr_T, for k in [0, k_C]
		for k in range(0, k_C + 1):
			MILP += pulp.lpSum(chi_T[(i, k)]
					for i in range(1, t_T + 1)) == clr_T[k], f"milp-2LMH-(6)-1-{k}_SG{self.name}"
			MILP += pulp.lpSum(chi_T[(i, k)]
					for i in range(1, t_T + 1)) <= t_T * delta_chi_T[k], f"milp-2LMH-(6)-2-{k}-1_SG{self.name}"
			MILP += pulp.lpSum(chi_T[(i, k)]
					for i in range(1, t_T + 1)) >= delta_chi_T[k], f"milp-2LMH-(6)-2-{k}-2_SG{self.name}"

		# -------- Constraint (7) --------
		# Constraint of e_T, v_T, chi_T, for i in [2, t_T]
		for i in range(2, t_T + 1):
			MILP += v_T[i - 1] >= v_T[i], f"milp-2LMH-(7)-1-{i}_SG{self.name}"
			MILP += chi_T[i - 1] - chi_T[i] <= k_C * (v_T[i - 1] - e_T[i]), f"milp-2LMH-(7)-2-{i}-1_SG{self.name}"
			MILP += chi_T[i - 1] - chi_T[i] >= v_T[i - 1] - e_T[i], f"milp-2LMH-(7)-2-{i}-2_SG{self.name}"

		# -------- Constraint (8) --------  IDO
		# Constraint of crl_T, n_G_int
		MILP += pulp.lpSum(clr_T[k] + 1 for k in range(1, k_C + 1) if k in I_lnk) + n_lnk_equal_one == n_lnk, f"mlp-2LMM-(8)_SG{self.name}"

		return MILP

	def prepare_variables_internal_vertices_and_edges(self,
		t_F, t_T, t_C, t_C_tilde, c_F, n_LB_int, n_UB_int, V_C, E_ge_one, E_ge_two, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C
	):
		# Define n_G_int
		n_G_int = pulp.LpVariable(f"n_G_int_SG{self.name}", n_LB_int, n_UB_int, cat=pulp.LpInteger)

		# Define v_F
		v_F = {i: pulp.LpVariable(f"v_F({i})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
					for i in range(1, t_F + 1)}
		# Define e_F
		e_F = {i : pulp.LpVariable(f"e_F({i})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
					for i in range(1, (t_F + 1) + 1)}
		# Define chi_F
		chi_F = {i : pulp.LpVariable(f"chi_F({i})_SG{self.name}", 0, c_F, cat = pulp.LpInteger)
					for i in range(1, t_F + 1)}
		# Define clr_F
		clr_F = {c : pulp.LpVariable(f"clr_F({c})_SG{self.name}", 0, t_F, cat = pulp.LpInteger)
					for c in range(0, c_F + 1)}
		# Define delta_chi_F
		delta_chi_F = {c : pulp.LpVariable(f"delta_chi_F({c})_SG{self.name}", bl_LB[c], 1, cat = pulp.LpInteger)
					for c in range(1, t_C_tilde + 1)}
		delta_chi_F[0] = pulp.LpVariable(f"delta_chi_F(0)_SG{self.name}", 0, 1, cat=pulp.LpBinary)
		delta_chi_F_tmp = {c: pulp.LpVariable(f"delta_chi_F({c})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for c in range(t_C_tilde + 1, c_F + 1)}
		delta_chi_F.update(delta_chi_F_tmp)
		chi_F_tmp = {(i, c) : pulp.LpVariable(f"chi_F({i},{c})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
					for i in range(1, t_F + 1) for c in range(0, c_F + 1)}
		chi_F.update(chi_F_tmp)
		# Define bl
		bl = {(k, i) : pulp.LpVariable(f"bl({k},{i})_SG{self.name}", 0, 1, cat = pulp.LpBinary)
					for k in (I_ge_two + I_ge_one) for i in range(1, t_T + 1)}

		return n_G_int, v_F, e_F, chi_F, clr_F, delta_chi_F, bl

	def add_constraints_internal_vertices_and_edges(self,
		# Model
		MILP,
		# Constants
		t_T, t_F, t_C, t_C_tilde, c_F, tail_F, I_ge_one, I_ge_two, bl_LB, bl_UB, E_C,
		# Binary Variables
		delta_chi_F, chi_T, e_F, v_F, v_T,
		# Integer Variables
		chi_F, clr_F, bl, n_G_int
	):
		# -------- Constraint (9) --------
		# Constraint of delta_chi_F, chi_F, for i in [1, t_F]
		for i in range(1, t_F + 1):
			MILP += chi_F[(i, 0)] == 1 - v_F[i], f"milp-2LMH-(9)-1-{i}_SG{self.name}"        
			MILP += pulp.lpSum(chi_F[(i, c)]
					for c in range(0, c_F + 1)) == 1, f"milp-2LMH-(9)-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(c * chi_F[(i, c)]
					for c in range(0, c_F + 1)) == chi_F[i], f"milp-2LMH-(9)-3-{i}_SG{self.name}"

		# -------- Constraint (10) --------
		# Constraint of delta_chi_F, clr_F, for c in [0, c_F]
		for c in range(0, c_F + 1):
			MILP += pulp.lpSum(chi_F[(i, c)]
					for i in range(1, t_F + 1)) == clr_F[c], f"milp-2LMH-(10)-1-{c}_SG{self.name}"
			MILP += pulp.lpSum(chi_F[(i, c)]
					for i in range(1, t_F + 1)) <= t_F * delta_chi_F[c], f"milp-2LMH-(10)-2-{c}-1_SG{self.name}"
			MILP += pulp.lpSum(chi_F[(i, c)]
					for i in range(1, t_F + 1)) >= delta_chi_F[c], f"milp-2LMH-(11)-2-{c}-2_SG{self.name}"

		# -------- Constraint (11) --------
		# Constraint of e_F
		MILP += e_F[1] == 0, f"milp-2LMH-(11)-1_SG{self.name}"
		MILP += e_F[t_F + 1] == 0, f"milp-2LMH-(11)-2_SG{self.name}"

		# -------- Constraint (12) --------
		# Constraint of e_F, chi_F, v_F, for i in [2, t_F]
		for i in range(2, t_F + 1):
			MILP += v_F[i - 1] >= v_F[i], f"milp-2LMH-(12)-1-{i}_SG{self.name}"
			MILP += chi_F[i - 1] - chi_F[i] <= c_F * (v_F[i - 1] - e_F[i]), f"milp-2LMH-(12)-2-{i}-1_SG{self.name}"
			MILP += chi_F[i - 1] - chi_F[i] >= v_F[i - 1] - e_F[i], f"milp-2LMH-(12)-2-{i}-2_SG{self.name}"

		# -------- Constraint (13) --------
		for k in (I_ge_two + I_ge_one):
			for i in range(1, t_T + 1):
				MILP += bl[(k, i)] >= delta_chi_F[t_C_tilde + i] + chi_T[(i, k)] - 1, f"milp-2LMH-(13)-{k}-{i}_SG{self.name}"

		# -------- Constraint (14) --------
		MILP += pulp.lpSum(bl[(k, i)] for k in (I_ge_one + I_ge_two) for i in range(1, t_T + 1)) <= \
				pulp.lpSum(delta_chi_F[t_C_tilde + i] for i in range(1, t_T + 1)), f"milp-2LMH-(14)_SG{self.name}"

		# -------- Constraint (15) --------
		for k in (I_ge_two + I_ge_one):
			MILP += pulp.lpSum(bl[(k, i)] for i in range(1, t_T + 1)) >= bl_LB[E_C[k]], f"milp-2LMH-(15)-{k}-1_SG{self.name}"
			MILP += pulp.lpSum(bl[(k, i)] for i in range(1, t_T + 1)) <= bl_UB[E_C[k]], f"milp-2LMH-(15)-{k}-2_SG{self.name}"

		# -------- Constraint (16) --------
		# Constraint of v_T (for specify sigma_co)
		MILP += t_C + pulp.lpSum(v_T[i] for i in range(1, t_T + 1)) + pulp.lpSum(v_F[i] for i in range(1, t_F + 1)) == n_G_int, f"milp-2LMH-(16)_SG{self.name}"

		return MILP

	def prepare_variables_fringe_trees(self,
		n_LB, n_star, rho, ch_LB, ch_UB, t_T, t_C, t_F, n_T, n_C, n_F, delta_i, I_ge_two, I_ge_one, V_C,
		E_ge_one, E_ge_two, v_T, v_F, F_C, F_T, F_F, Code_F, Gamma_lf_ac, ac_LB_lf, ac_UB_lf
	):
		# Define n_G
		n_G = pulp.LpVariable(f"n_G_SG{self.name}", n_LB, n_star, cat=pulp.LpInteger)

		# Define h_T, h_C, h_F
		h_T = {i: pulp.LpVariable(f"h_T({i})_SG{self.name}", 0, rho, cat=pulp.LpInteger)
				for i in range(1, t_T + 1)}
		h_C = {i: pulp.LpVariable(f"h_C({i})_SG{self.name}", 0, rho, cat=pulp.LpInteger)
				for i in range(1, t_C + 1)}
		h_F = {i: pulp.LpVariable(f"h_F({i})_SG{self.name}", 0, rho, cat=pulp.LpInteger)
				for i in range(1, t_F + 1)}

		# Define delta_fr_C, delta_fr_T, delta_fr_F
		delta_fr_C = {(i, Code_F[psi]): pulp.LpVariable(f"delta_fr_C({i},{Code_F[psi]})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_C + 1) for psi in F_C[i]}
		delta_fr_T = {(i, Code_F[psi]): pulp.LpVariable(f"delta_fr_T({i},{Code_F[psi]})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_T + 1) for psi in F_T[i]}
		delta_fr_F = {(i, Code_F[psi]): pulp.LpVariable(f"delta_fr_F({i},{Code_F[psi]})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_F + 1) for psi in F_F[i]}

		# Define deg_ex_C, deg_ex_T, deg_ex_F
		deg_ex_C = {i: pulp.LpVariable(f"deg_ex_C({i})_SG{self.name}", 0, self.MAX_VAL - 1, cat=pulp.LpInteger)
					for i in range(1, t_C + 1)}
		deg_ex_T = {i: pulp.LpVariable(f"deg_ex_T({i})_SG{self.name}", 0, self.MAX_VAL - 1, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		deg_ex_F = {i: pulp.LpVariable(f"deg_ex_F({i})_SG{self.name}", 0, self.MAX_VAL - 1, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}    

		# Define hyddeg_C, hyddeg_T, hyddeg_F
		hyddeg_C = {i: pulp.LpVariable(f"hyddeg_C({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_C + 1)}
		hyddeg_T = {i: pulp.LpVariable(f"hyddeg_T({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		hyddeg_F = {i: pulp.LpVariable(f"hyddeg_F({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}   

		# Define eledeg_C, eledeg_T, eledeg_F
		eledeg_C = {i: pulp.LpVariable(f"eledeg_C({i})_SG{self.name}", -3, 3, cat=pulp.LpInteger)
					for i in range(1, t_C + 1)}
		eledeg_T = {i: pulp.LpVariable(f"eledeg_T({i})_SG{self.name}", -3, 3, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		eledeg_F = {i: pulp.LpVariable(f"eledeg_F({i})_SG{self.name}", -3, 3, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}    

		# Define sigma
		sigma = {(k, i): pulp.LpVariable(f"sigma({k},{i})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for k in (I_ge_two + I_ge_one) for i in range(1, t_T + 1)} 

		# Define ac_lf
		ac_lf = {nu: pulp.LpVariable(f"ac_lf({nu})_SG{self.name}", ac_LB_lf[nu], ac_UB_lf[nu], cat=pulp.LpInteger)
					for nu in Gamma_lf_ac}

		return n_G, h_T, h_C, h_F, delta_fr_C, delta_fr_T, delta_fr_F, \
			deg_ex_C, deg_ex_T, deg_ex_F, hyddeg_C, hyddeg_T, hyddeg_F, eledeg_C, eledeg_T, eledeg_F, sigma, ac_lf

	def add_constraints_fringe_trees(self,
		# Model
		MILP,
		# Constants
		t_T, t_F, t_C, t_C_tilde, c_F, I_ge_one, I_ge_two, rho, n_star, n_T, n_C, n_F, ch_LB, ch_UB,
		E_C, F_C, F_T, F_F, F_Cp, F_Tp, F_Fp, Code_F, val, n_psi_H, deg_r_H, deg_r_hyd, ht_H, v_ion,
		ac_psi_lf, Gamma_lf_ac,
		# Binary Variables
		delta_chi_F, delta_chi_T, e_F, v_F, v_T, h_T, h_C, h_F, sigma, delta_fr_C, delta_fr_T, delta_fr_F,
		# Integer Variables
		chi_T, chi_F, clr_F, n_G, deg_ex_C, deg_ex_T, deg_ex_F, hyddeg_C, hyddeg_T, hyddeg_F,
		eledeg_C, eledeg_T, eledeg_F, ac_lf
	):
		# -------- Constraint (18) & (19) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(delta_fr_C[(i, Code_F[psi])] 
						for psi in F_C[i]) == 1, f"milp-2LMH-(18)-C-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(deg_r_H[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
						for psi in F_C[i]) == deg_ex_C[i], f"milp-2LMH-(19)-C-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(deg_r_hyd[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
						for psi in F_C[i]) == hyddeg_C[i], f"milp-2LMH-(19)-C-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(v_ion[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
						for psi in F_C[i]) == eledeg_C[i], f"milp-2LMH-(19)-C-3-{i}_SG{self.name}"

		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(delta_fr_T[(i, Code_F[psi])] 
						for psi in F_T[i]) == v_T[i], f"milp-2LMH-(18)-T-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(deg_r_H[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
						for psi in F_T[i]) == deg_ex_T[i], f"milp-2LMH-(19)-T-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(deg_r_hyd[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
						for psi in F_T[i]) == hyddeg_T[i], f"milp-2LMH-(19)-T-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(v_ion[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
						for psi in F_T[i]) == eledeg_T[i], f"milp-2LMH-(19)-T-3-{i}_SG{self.name}"

		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(delta_fr_F[(i, Code_F[psi])] 
						for psi in F_F[i]) == v_F[i], f"milp-2LMH-(18)-F-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(deg_r_H[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
						for psi in F_F[i]) == deg_ex_F[i], f"milp-2LMH-(19)-F-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(deg_r_hyd[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
						for psi in F_F[i]) == hyddeg_F[i], f"milp-2LMH-(19)-F-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(v_ion[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
						for psi in F_F[i]) == eledeg_F[i], f"milp-2LMH-(19)-F-3-{i}_SG{self.name}"

		# -------- Constraint (20) --------
		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(delta_fr_F[(i, Code_F[psi])]
						for psi in F_Fp[i][rho]) >= v_F[i] - e_F[i + 1], f"milp-2LMH-(20)-{i}_SG{self.name}"

		# -------- Constraint (21) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(ht_H[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] 
						for psi in F_C[i]) == h_C[i], f"milp-2LMH-(21)-C-{i}_SG{self.name}"

		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(ht_H[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] 
						for psi in F_T[i]) == h_T[i], f"milp-2LMH-(21)-T-{i}_SG{self.name}"

		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(ht_H[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] 
						for psi in F_F[i]) == h_F[i], f"milp-2LMH-(21)-F-{i}_SG{self.name}"

		# -------- Constraint (22) --------
		MILP += pulp.lpSum(n_psi_H[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
					for i in range(1, t_C + 1) for psi in F_C[i]) + \
				pulp.lpSum(n_psi_H[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
					for i in range(1, t_T + 1) for psi in F_T[i]) + \
				pulp.lpSum(n_psi_H[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
					for i in range(1, t_F + 1) for psi in F_F[i]) + \
				pulp.lpSum(v_T[i] for i in range(1, t_T + 1)) + \
				pulp.lpSum(v_F[i] for i in range(1, t_F + 1)) + t_C == n_G, f"milp-2LMH-(22)_SG{self.name}"

		# -------- Constraint (ac_lf) -------- 
		for nu in Gamma_lf_ac:
			MILP += pulp.lpSum(ac_psi_lf[nu][Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
						for i in range(1, t_C + 1) for psi in F_C[i]) + \
					pulp.lpSum(ac_psi_lf[nu][Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
						for i in range(1, t_T + 1) for psi in F_T[i]) + \
					pulp.lpSum(ac_psi_lf[nu][Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
						for i in range(1, t_F + 1) for psi in F_F[i]) == ac_lf[nu], f"milp-2LMM-ac-lf-{nu}_SG{self.name}"

		# -------- Constraint (23) --------
		for i in range(1, t_C_tilde + 1):
			MILP += h_C[i] >= ch_LB[i] - n_star * delta_chi_F[i], f"milp-2LMH-(23)-1-{i}_SG{self.name}"
			MILP += clr_F[i] + rho >= ch_LB[i], f"milp-2LMH-(23)-2-{i}_SG{self.name}"
			MILP += h_C[i] <= ch_UB[i], f"milp-2LMH-(23)-3-{i}_SG{self.name}"
			MILP += clr_F[i] + rho <= ch_UB[i] + n_star * (1 - delta_chi_F[i]), f"milp-2LMH-(23)-4-{i}_SG{self.name}"

		# -------- Constraint (24) --------
		for i in range(t_C_tilde + 1, t_C + 1):
			MILP += h_C[i] >= ch_LB[i], f"milp-2LMH-(24)-1-{i}_SG{self.name}"
			MILP += h_C[i] <= ch_UB[i], f"milp-2LMH-(24)-2-{i}_SG{self.name}" 

		# -------- Constraint (25) --------
		for i in range(1, t_T + 1):
			for k in (I_ge_one + I_ge_two):
				MILP += h_T[i] <= ch_UB[E_C[k]] + \
					n_star * (delta_chi_F[t_C_tilde + i] + 1 - chi_T[(i, k)]), f"milp-2LMH-(25)-1-{i}-{k}_SG{self.name}"
				MILP += clr_F[t_C_tilde + i] + rho <= ch_UB[E_C[k]] + \
					n_star * (2 - delta_chi_F[t_C_tilde + i] - chi_T[(i, k)]), f"milp-2LMH-(25)-2-{i}-{k}_SG{self.name}"

		# -------- Constraint (26) --------
		for k in (I_ge_one + I_ge_two):
			MILP += pulp.lpSum(sigma[(k, i)] for i in range(1, t_T + 1)) == delta_chi_T[k], f"milp-2LMH-(26)-{k}_SG{self.name}"

		# -------- Constraint (27) --------
		for i in range(1, t_T + 1):
			for k in (I_ge_one + I_ge_two):
				MILP += chi_T[(i, k)] >= sigma[(k, i)], f"milp-2LMH-(27)-1-{i}-{k}_SG{self.name}"
				MILP += h_T[i] >= ch_LB[E_C[k]] - \
					n_star * (delta_chi_F[t_C_tilde + i] + 1 - sigma[(k, i)]), f"milp-2LMH-(27)-2-{i}-{k}_SG{self.name}"
				MILP += clr_F[t_C_tilde + i] + rho >= ch_LB[E_C[k]] - \
					n_star * (2 - delta_chi_F[t_C_tilde + i] - sigma[(k, i)]), f"milp-2LMH-(27)-3-{i}-{k}_SG{self.name}"

		return MILP

	def prepare_variables_degree(self,
		t_C, t_T, t_F, n_C, n_T, n_F, delta_i, n_star, dg_LB, dg_UB, rho
	):
		# Define deg_C, deg_T, deg_F
		deg_C = {i: pulp.LpVariable(f"deg_C({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
				for i in range(1, t_C + 1)}
		deg_T = {i: pulp.LpVariable(f"deg_T({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
				for i in range(1, t_T + 1)}
		deg_F = {i: pulp.LpVariable(f"deg_F({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
				for i in range(1, t_F + 1)}

		# Define deg_CT, deg_TC
		deg_CT = {i: pulp.LpVariable(f"deg_CT({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
				for i in range(1, t_C + 1)}
		deg_TC = {i: pulp.LpVariable(f"deg_TC({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
				for i in range(1, t_C + 1)}

		# Define delta_dg_C, delta_dg_T, delta_dg_F
		delta_dg_C = {(i, d): pulp.LpVariable(f"delta_dg_C({i},{d})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(1, t_C + 1) for d in range(1, self.MAX_VAL + 1)}
		delta_dg_T = {(i, d): pulp.LpVariable(f"delta_dg_T({i},{d})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(1, t_T + 1) for d in range(0, self.MAX_VAL + 1)}
		delta_dg_F = {(i, d): pulp.LpVariable(f"delta_dg_F({i},{d})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(1, t_F + 1) for d in range(0, self.MAX_VAL + 1)}

		# Define dg
		dg = {d: pulp.LpVariable(f"dg({d})_SG{self.name}", dg_LB[d], dg_UB[d], cat=pulp.LpInteger) for d in range(1, self.MAX_VAL + 1)} 

		# Define deg_int_C, deg_int_T, deg_int_F
		deg_int_C = {i: pulp.LpVariable(f"deg_int_C({i})_SG{self.name}", 1, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_C + 1)}
		deg_int_T = {i: pulp.LpVariable(f"deg_int_T({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		deg_int_F = {i: pulp.LpVariable(f"deg_int_F({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}

		#Define delta_int_dg_C, delta_int_dg_T, delta_int_dg_F
		delta_int_dg_C = {(i, d): pulp.LpVariable(f"delta_int_dg_C({i},{d})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_C + 1) for d in range(1, self.MAX_VAL + 1)}
		delta_int_dg_T = {(i, d): pulp.LpVariable(f"delta_int_dg_T({i},{d})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_T + 1) for d in range(0, self.MAX_VAL + 1)}
		delta_int_dg_F = {(i, d): pulp.LpVariable(f"delta_int_dg_F({i},{d})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_F + 1) for d in range(0, self.MAX_VAL + 1)}

		# Define dg_int
		dg_int = {d: pulp.LpVariable(f"dg_int({d})_SG{self.name}", dg_LB[d], dg_UB[d], cat=pulp.LpInteger) for d in range(1, self.MAX_VAL + 1)}

		return deg_C, deg_T, deg_F, deg_CT, deg_TC, \
			delta_dg_C, delta_dg_T, delta_dg_F, dg, \
			deg_int_C, deg_int_T, deg_int_F, \
			delta_int_dg_C, delta_int_dg_T, delta_int_dg_F, \
			dg_int

	def add_constraints_degree(self,
		# Model
		MILP,
		# Constants
		t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, delta_i, I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus,
		rho, F_C, F_T, F_F, F_Cp, F_Tp, F_Fp, Code_F, deg_fr,
		# Binary Variables
		delta_dg_C, delta_dg_T, delta_dg_F, delta_int_dg_C, delta_int_dg_T, delta_int_dg_F,
		e_T, e_F, v_F, v_T, delta_chi_T, delta_chi_F, delta_fr_C, delta_fr_T, delta_fr_F,
		# Integer Variables
		deg_C, deg_T, deg_F, deg_CT, deg_TC, dg, deg_int_C, deg_int_T, deg_int_F, dg_int, 
		deg_tilde_C_minus, deg_tilde_C_plus, deg_ex_C, deg_ex_T, deg_ex_F, hyddeg_C, hyddeg_T, hyddeg_F
	):
		# -------- Constraint (28) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(delta_chi_T[k] 
						for k in (I_ge_two_plus[i] + I_ge_one_plus[i])) == deg_CT[i], f"milp-2LMH-(28)-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_chi_T[k] 
						for k in (I_ge_two_minus[i] + I_ge_one_minus[i])) == deg_TC[i], f"milp-2LMH-(28)-2-{i}_SG{self.name}"

		# -------- Constraint (29) --------
		for i in range(1, t_C_tilde + 1):
			MILP += deg_tilde_C_minus[i] + deg_tilde_C_plus[i] + \
					deg_CT[i] + deg_TC[i] + \
					delta_chi_F[i] == deg_int_C[i], f"milp-2LMH-(29)-{i}_SG{self.name}"

		# -------- Constraint (30) --------
		for i in range(t_C_tilde + 1, t_C + 1):
			MILP += deg_tilde_C_minus[i] + deg_tilde_C_plus[i] + \
					deg_CT[i] + deg_TC[i] == deg_int_C[i], f"milp-2LMH-(30)-{i}_SG{self.name}"

		# -------- Constraint (31) --------
		for i in range(1, t_C + 1):
			MILP += deg_int_C[i] + deg_ex_C[i] == deg_C[i], f"milp-2LMH-(31)-{i}_SG{self.name}"

		# -------- Constraint (32) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(delta_fr_C[(i, Code_F[psi])]
						for psi in F_Cp[i][rho]) >= 2 - deg_int_C[i], f"milp-(32)-{i}_SG{self.name}"

		# -------- Constraint (33) --------
		MILP += e_T[1] == 0, f"milp-2LMH-(33)-T1_SG{self.name}"
		MILP += e_T[t_T + 1] == 0, f"milp-2LMH-(33)-T2_SG{self.name}"
		for i in range(1, t_T + 1):
			MILP += 2 * v_T[i] + delta_chi_F[t_C_tilde + i] == \
					deg_int_T[i], f"milp-2LMH-(33)-1-{i}_SG{self.name}"
			MILP += deg_int_T[i] + deg_ex_T[i] == deg_T[i], f"milp-2LMH-(33)-2-{i}_SG{self.name}"

		# -------- Constraint (34) --------
		for i in range(1, t_F + 1):
			MILP += v_F[i] + e_F[i + 1] == deg_int_F[i], f"milp-2LMH-(34)-1-{i}_SG{self.name}"
			MILP += deg_int_F[i] + deg_ex_F[i] == deg_F[i], f"milp-2LMH-(34)-2-{i}_SG{self.name}"

		# -------- Constraint (35) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(delta_dg_C[(i, d)] 
					for d in range(1, self.MAX_VAL + 1)) == 1, f"milp-2LMH-(35)-C-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_dg_C[(i, d)] 
					for d in range(1, self.MAX_VAL + 1)) == deg_C[i], f"milp-2LMH-(35)-C-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_int_dg_C[(i, d)]
					for d in range(1, self.MAX_VAL + 1)) == 1, f"milp-2LMH-(35)-C-3-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_int_dg_C[(i, d)] 
					for d in range(1, self.MAX_VAL + 1)) == deg_int_C[i], f"milp-2LMH-(35)-C-4-{i}_SG{self.name}"

		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(delta_dg_T[(i, d)] 
					for d in range(0, self.MAX_VAL + 1)) == 1, f"milp-2LMH-(35)-T-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_dg_T[(i, d)]
					for d in range(1, self.MAX_VAL + 1)) == deg_T[i], f"milp-2LMH-(35)-T-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_int_dg_T[(i, d)]
					for d in range(0, self.MAX_VAL + 1)) == 1, f"milp-2LMH-(35)-T-3-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_int_dg_T[(i, d)] 
					for d in range(1, self.MAX_VAL + 1)) == deg_int_T[i], f"milp-2LMH-(35)-T-4-{i}_SG{self.name}"

		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(delta_dg_F[(i, d)] 
					for d in range(0, self.MAX_VAL + 1)) == 1, f"milp-2LMH-(35)-F-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_dg_F[(i, d)]
					for d in range(1, self.MAX_VAL + 1)) == deg_F[i], f"milp-2LMH-(35)-F-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_int_dg_F[(i, d)]
					for d in range(0, self.MAX_VAL + 1)) == 1, f"milp-2LMH-(35)-F-3-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_int_dg_F[(i, d)] 
					for d in range(1, self.MAX_VAL + 1)) == deg_int_F[i], f"milp-2LMH-(35)-F-4-{i}_SG{self.name}"

		# -------- Constraint (36) --------
		for d in range(1, self.MAX_VAL + 1):
			MILP += pulp.lpSum(delta_dg_C[(i, d)] for i in range(1, t_C + 1)) + \
					pulp.lpSum(delta_dg_T[(i, d)] for i in range(1, t_T + 1)) + \
					pulp.lpSum(delta_dg_F[(i, d)] for i in range(1, t_F + 1)) == dg[d], f"milp-2LMH-(36)-1-{d}_SG{self.name}"

		for d in range(1, self.MAX_VAL + 1):
			MILP += pulp.lpSum(delta_int_dg_C[(i, d)] for i in range(1, t_C + 1)) + \
					pulp.lpSum(delta_int_dg_T[(i, d)] for i in range(1, t_T + 1)) + \
					pulp.lpSum(delta_int_dg_F[(i, d)] for i in range(1, t_F + 1)) == dg_int[d], f"milp-2LMH-(36)-2-{d}_SG{self.name}"

		return MILP

	def prepare_variables_multiplicity(self,
		t_C, t_T, t_F, n_C, n_T, n_F, c_F, delta_i, I_equal_one, I_zero_one, I_ge_one, I_ge_two,
		E_C, bd2_LB, bd2_UB, bd3_LB, bd3_UB, n_UB_int, n_star, rho
	):
		# Define beta_T, beta_F
		beta_T = {i: pulp.LpVariable(f"beta_T({i})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
				for i in range(1, t_T + 2)} # use [1, t_T + 1] just for convenience
		beta_F = {i: pulp.LpVariable(f"beta_F({i})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
				for i in range(1, t_F + 2)} # use [1, t_F + 1] just for convenience

		# Define beta_C
		beta_C = {i: pulp.LpVariable(f"beta_C({i})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
				for i in (I_equal_one + I_zero_one + I_ge_one)}

		# Define beta_plus, beta_minus
		beta_plus = {k: pulp.LpVariable(f"beta_plus({k})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
					for k in (I_ge_two + I_ge_one)}
		beta_minus = {k: pulp.LpVariable(f"beta_minus({k})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
					for k in (I_ge_two + I_ge_one)}

		# Define beta_in
		beta_in = {c: pulp.LpVariable(f"beta_in({c})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
					for c in range(1, c_F + 1)}

		# Define beta_ex_C, beta_ex_T, beta_ex_F
		beta_ex_C = {i: pulp.LpVariable(f"beta_ex_C({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_C + 1)}
		beta_ex_T = {i: pulp.LpVariable(f"beta_ex_T({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		beta_ex_F = {i: pulp.LpVariable(f"beta_ex_F({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}    

		# Define delta_beta_T, delta_beta_F
		delta_beta_T = {(i, m): pulp.LpVariable(f"delta_beta_T({i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(2, t_T + 1) for m in range(self.MAX_BOND + 1)}
		delta_beta_F = {(i, m): pulp.LpVariable(f"delta_beta_F({i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(2, t_F + 1) for m in range(self.MAX_BOND + 1)}

		# Define delta_beta_C
		delta_beta_C = {(i, m): pulp.LpVariable(f"delta_beta_C({i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in (I_equal_one + I_zero_one + I_ge_one) for m in range(self.MAX_BOND + 1)}

		# Define delta_beta_plus, delta_beta_in
		delta_beta_plus = {(k, m): pulp.LpVariable(f"delta_beta_plus({k},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
							for k in (I_ge_two + I_ge_one) for m in range(self.MAX_BOND + 1)}
		delta_beta_minus = {(k, m): pulp.LpVariable(f"delta_beta_minus({k},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
							for k in (I_ge_two + I_ge_one) for m in range(self.MAX_BOND + 1)}
		
		# Define delta_beta_in
		delta_beta_in = {(c, m): pulp.LpVariable(f"delta_beta_in({c},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
							for c in range(1, c_F + 1) for m in range(self.MAX_BOND + 1)}
		
		# Define bd_int
		bd_int = {m: pulp.LpVariable(f"bd_int({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}

		# Define bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF
		bd_T = {m: pulp.LpVariable(f"bd_T({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}
		bd_C = {m: pulp.LpVariable(f"bd_C({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}
		bd_CT = {m: pulp.LpVariable(f"bd_CT({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}
		bd_TC = {m: pulp.LpVariable(f"bd_TC({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}
		bd_F = {m: pulp.LpVariable(f"bd_F({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}
		bd_CF = {m: pulp.LpVariable(f"bd_CF({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}
		bd_TF = {m: pulp.LpVariable(f"bd_TF({m})_SG{self.name}", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, self.MAX_BOND + 1)}
		
		return beta_T, beta_F, beta_C, beta_plus, beta_minus, beta_in, beta_ex_C, beta_ex_T, beta_ex_F, \
			delta_beta_T, delta_beta_F, delta_beta_C, \
			delta_beta_plus, delta_beta_minus, delta_beta_in, bd_int, \
			bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF

	def add_constraints_multiplicity(self,
		# Model
		MILP,
		# Constants
		t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, c_F, delta_i, head_C, tail_C, I_equal_one, I_zero_one,
		I_ge_one, I_ge_two, E_C, bd2_LB, bd2_UB, bd3_LB, bd3_UB, rho, F_C, F_T, F_F, Code_F, beta_r,
		# Binary Variables
		e_C, e_T, e_F, v_F, v_T, delta_chi_T, delta_chi_F, delta_beta_C, delta_beta_T, delta_beta_F,
		delta_beta_plus, delta_beta_minus, delta_beta_in, delta_fr_C, delta_fr_T, delta_fr_F,
		# Integer Variables
		beta_C, beta_T, beta_F, beta_plus, beta_minus, beta_in, bd_int, bd_C, bd_T, bd_CT, bd_TC,
		bd_F, bd_CF, bd_TF, beta_ex_C, beta_ex_T, beta_ex_F
	):
		# -------- Constraint (37) --------
		for i in (I_equal_one + I_zero_one + I_ge_one):
			MILP += beta_C[i] >= e_C[i], f"milp-2LMH-(37)-{i}-1_SG{self.name}"
			MILP += beta_C[i] <= self.MAX_BOND * e_C[i], f"milp-2LMH-(37)-{i}-2_SG{self.name}"

		# -------- Constraint (38) --------
		for i in range(2, t_T + 1):
			MILP += beta_T[i] >= e_T[i], f"milp-2LMH-(38)-T-{i}-1_SG{self.name}"
			MILP += beta_T[i] <= self.MAX_BOND * e_T[i], f"milp-2LMH-(38)-T-{i}-2_SG{self.name}"
		for i in range(2, t_F + 1):
			MILP += beta_F[i] >= e_F[i], f"milp-2LMH-(38)-F-{i}-1_SG{self.name}"
			MILP += beta_F[i] <= self.MAX_BOND * e_F[i], f"milp-2LMH-(38)-F-{i}-2_SG{self.name}"

		# -------- Constraint (39) --------
		for k in (I_ge_two + I_ge_one):
			MILP += beta_plus[k] >= delta_chi_T[k], f"milp-2LMH-(39)-1-{k}-1_SG{self.name}"
			MILP += beta_plus[k] <= self.MAX_BOND * delta_chi_T[k], f"milp-2LMH-(39)-1-{k}-2_SG{self.name}"
			MILP += beta_minus[k] >= delta_chi_T[k], f"milp-2LMH-(39)-2-{k}-1_SG{self.name}"
			MILP += beta_minus[k] <= self.MAX_BOND * delta_chi_T[k], f"milp-2LMH-(39)-2-{k}-2_SG{self.name}"

		# -------- Constraint (40) --------
		for c in range(1, c_F + 1):
			MILP += beta_in[c] >= delta_chi_F[c], f"milp-2LMH-(40)-{c}-1_SG{self.name}"
			MILP += beta_in[c] <= self.MAX_BOND * delta_chi_F[c], f"milp-2LMH-(40)-{c}-2_SG{self.name}"

		# -------- Constraint (41) --------
		for i in range(2, t_T + 1):
			MILP += pulp.lpSum(delta_beta_T[(i, m)] 
				for m in range(self.MAX_BOND + 1)) == 1, f"milp-2LMH-(41)-T-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_beta_T[(i, m)] 
				for m in range(self.MAX_BOND + 1)) == beta_T[i], f"milp-2LMH-(41)-T-2-{i}_SG{self.name}"
		for i in range(2, t_F + 1):
			MILP += pulp.lpSum(delta_beta_F[(i, m)] 
				for m in range(self.MAX_BOND + 1)) == 1, f"milp-2LMH-(41)-F-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_beta_F[(i, m)] 
				for m in range(self.MAX_BOND + 1)) == beta_F[i], f"milp-2LMH-(41)-F-2-{i}_SG{self.name}"

		# -------- Constraint (42) --------
		for i in (I_equal_one + I_zero_one + I_ge_one):
			MILP += pulp.lpSum(delta_beta_C[(i, m)] 
				for m in range(self.MAX_BOND + 1)) == 1, f"milp-2LMH-(42)-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_beta_C[(i, m)] 
				for m in range(self.MAX_BOND + 1)) == beta_C[i], f"milp-2LMH-(42)-2-{i}_SG{self.name}"

		# -------- Constraint (43) --------
		for k in (I_ge_two + I_ge_one):
			MILP += pulp.lpSum(delta_beta_plus[(k, m)]
				for m in range(self.MAX_BOND + 1)) == 1, f"milp-2LMH-(43)-1-{k}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_beta_plus[(k, m)]
				for m in range(self.MAX_BOND + 1)) == beta_plus[k], f"milp-2LMH-(43)-2-{k}_SG{self.name}"
		for k in (I_ge_two + I_ge_one):
			MILP += pulp.lpSum(delta_beta_minus[(k, m)]
				for m in range(self.MAX_BOND + 1)) == 1, f"milp-2LMH-(43)-3-{k}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_beta_minus[(k, m)]
				for m in range(self.MAX_BOND + 1)) == beta_minus[k], f"milp-2LMH-(43)-4-{k}_SG{self.name}"
		for c in range(1, c_F + 1):
			MILP += pulp.lpSum(delta_beta_in[(c, m)]
				for m in range(self.MAX_BOND + 1)) == 1, f"milp-2LMH-(43)-5-{c}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_beta_in[(c, m)]
				for m in range(self.MAX_BOND + 1)) == beta_in[c], f"milp-2LMH-(43)-6-{c}_SG{self.name}"

		# -------- Constraint (44) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(beta_r[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] 
						for psi in F_C[i]) == beta_ex_C[i], f"milp-2LMH-(44)-C-{i}_SG{self.name}"
		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(beta_r[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] 
						for psi in F_T[i]) == beta_ex_T[i], f"milp-2LMH-(44)-T-{i}_SG{self.name}"
		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(beta_r[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] 
						for psi in F_F[i]) == beta_ex_F[i], f"milp-2LMH-(44)-F-{i}_SG{self.name}"

		# -------- Constraint (45) --------
		for m in range(1, self.MAX_BOND + 1):
			MILP += pulp.lpSum(delta_beta_C[(i, m)] for i in (I_equal_one + I_zero_one + I_ge_one)) == bd_C[m], f"milp-2LMH-(45)-1-{m}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)) == bd_T[m], f"milp-2LMH-(45)-2-{m}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_plus[(k, m)] for k in (I_ge_two + I_ge_one)) == bd_CT[m], f"milp-2LMH-(45)-3-{m}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_minus[(k, m)] for k in (I_ge_two + I_ge_one)) == bd_TC[m], f"milp-2LMH-(45)-4-{m}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)) == bd_F[m], f"milp-2LMH-(45)-5-{m}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)) == bd_CF[m], f"milp-2LMH-(45)-6-{m}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)) == bd_TF[m], f"milp-2LMH-(45)-7-{m}_SG{self.name}"
			MILP += bd_C[m] + bd_T[m] + bd_F[m] + \
					bd_CT[m] + bd_TC[m] + bd_TF[m] + bd_CF[m] == bd_int[m], f"milp-2LMH-(45)-8-{m}_SG{self.name}"

		return MILP

	def prepare_variables_chemical_elements(self,
		t_C, t_T, t_F, n_C, n_T, n_F, delta_i, n_star, Lambda, epsilon, na_LB, na_UB, rho, na_LB_int, na_UB_int,
		Lambda_int, Lambda_ex, MAX_CODE, MAX_CODE_int, MAX_CODE_ex, Code_Lambda_int, Code_Lambda_ex,
		F_C, F_T, F_F, Code_F, alpha_r, na_alpha_ex
	):
		# Define beta_CT, beta_TC
		beta_CT = {i: pulp.LpVariable(f"beta_CT({i})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		beta_TC = {i: pulp.LpVariable(f"beta_TC({i})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}

		# Define beta_CF, beta_TF
		beta_CF = {i: pulp.LpVariable(f"beta_CF({i})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}
		beta_TF = {i: pulp.LpVariable(f"beta_TF({i})_SG{self.name}", 0, self.MAX_BOND, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}

		# Define alpha_C, alpha_T, alpha_F
		alpha_C = {i: pulp.LpVariable(f"alpha_C({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
					for i in range(1, t_C + 1)}
		alpha_T = {i: pulp.LpVariable(f"alpha_T({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		alpha_F = {i: pulp.LpVariable(f"alpha_F({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
					for i in range(1, t_F + 1)}

		# # Define delta_alpha_C, delta_alpha_T, delta_alpha_F
		delta_alpha_C = {(i, mu): pulp.LpVariable(f"delta_alpha_C({i},{mu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_C + 1) for mu in (Lambda_int + [epsilon])}
		delta_alpha_T = {(i, mu): pulp.LpVariable(f"delta_alpha_T({i},{mu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_T + 1) for mu in (Lambda_int + [epsilon])}
		delta_alpha_F = {(i, mu): pulp.LpVariable(f"delta_alpha_F({i},{mu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_F + 1) for mu in (Lambda_int + [epsilon])}

		# Define MASS
		MASS = pulp.LpVariable(f"Mass_SG{self.name}", cat=pulp.LpInteger)

		# Define na
		na = {atom: pulp.LpVariable(f"na({atom})_SG{self.name}", na_LB[atom], na_UB[atom], cat=pulp.LpInteger)
				for atom in Lambda}

		# Define na_int, na_C, na_T, na_F
		na_int = {atom: pulp.LpVariable(f"na_int({atom})_SG{self.name}", na_LB_int[atom], na_UB_int[atom], cat=pulp.LpInteger) 
				for atom in Lambda_int}
		na_C = {atom: pulp.LpVariable(f"na_C({atom})_SG{self.name}", 0, na_UB_int[atom], cat=pulp.LpInteger) 
				for atom in Lambda_int}
		na_T = {atom: pulp.LpVariable(f"na_T({atom})_SG{self.name}", 0, na_UB_int[atom], cat=pulp.LpInteger)
				for atom in Lambda_int}
		na_F = {atom: pulp.LpVariable(f"na_F({atom})_SG{self.name}", 0, na_UB_int[atom], cat=pulp.LpInteger)
				for atom in Lambda_int}

		# Define na_ex_C, na_ex_T, na_ex_F, na_ex
		na_ex_C = {atom: pulp.LpVariable(f"na_ex_C({atom})_SG{self.name}", 0, na_UB[atom], cat=pulp.LpInteger)
				for atom in Lambda_ex}
		na_ex_T = {atom: pulp.LpVariable(f"na_ex_T({atom})_SG{self.name}", 0, na_UB[atom], cat=pulp.LpInteger)
				for atom in Lambda_ex}
		na_ex_F = {atom: pulp.LpVariable(f"na_ex_F({atom})_SG{self.name}", 0, na_UB[atom], cat=pulp.LpInteger)
				for atom in Lambda_ex}
		na_ex = {atom: pulp.LpVariable(f"na_ex({atom})_SG{self.name}", 0, na_UB[atom], cat=pulp.LpInteger)
				for atom in Lambda_ex}

		return beta_CT, beta_TC, beta_CF, beta_TF, alpha_C, alpha_T, alpha_F, MASS, \
			delta_alpha_C, delta_alpha_T, delta_alpha_F, na, na_int, na_C, na_T, na_F, \
			na_ex_C, na_ex_T, na_ex_F, na_ex

	def add_constraints_chemical_elements(self,
		# Model
		MILP,
		# Constants
		t_T, t_F, t_C, t_C_tilde, n_T, n_C, n_F, c_F, delta_i, E_C, I_equal_one, I_zero_one, I_ge_one, I_ge_one_plus, I_ge_one_minus,
		I_ge_two, I_ge_two_plus, I_ge_two_minus, val, mass, Lambda, Code_Lambda, Lambda_star, epsilon, rho,
		na_LB_int, na_UB_int, Lambda_int, Lambda_ex, MAX_CODE, MAX_CODE_int, MAX_CODE_ex, Code_Lambda_int, Code_Lambda_ex,
		F_C, F_T, F_F, Code_F, alpha_r, na_alpha_ex,
		# Binary Variables
		v_T, v_F, e_T, e_F, delta_chi_T, delta_chi_F, chi_T, chi_F, delta_alpha_C, delta_alpha_T, delta_alpha_F,
		delta_fr_C, delta_fr_T, delta_fr_F,
		# Integer Variables
		deg_C, deg_T, deg_F, beta_C, beta_T, beta_F, beta_CT, beta_TC, beta_CF, beta_TF, beta_plus, beta_minus, beta_in,
		alpha_C, alpha_T, alpha_F, MASS, na, na_int, na_C, na_T, na_F, na_ex_C, na_ex_T, na_ex_F, beta_ex_C, beta_ex_T, beta_ex_F,
		na_ex, bd_int, eledeg_C, eledeg_T, eledeg_F
	):
		# -------- Constraint (46) --------
		for k in (I_ge_one + I_ge_two):
			for i in range(1, t_T + 1):
				MILP += beta_CT[i] >= beta_plus[k] - self.MAX_BOND * (e_T[i] - chi_T[(i, k)] + 1), f"milp-2LMH-(46)-1-{k}-{i}-1_SG{self.name}"
				MILP += beta_CT[i] <= beta_plus[k] + self.MAX_BOND * (e_T[i] - chi_T[(i, k)] + 1), f"milp-2LMH-(46)-1-{k}-{i}-2_SG{self.name}"
			for i in range(1, t_T + 1):
				MILP += beta_TC[i] >= beta_minus[k] - self.MAX_BOND * (e_T[i + 1] - chi_T[(i, k)] + 1), f"milp-2LMH-(46)-2-{k}-{i}-1_SG{self.name}"
				MILP += beta_TC[i] <= beta_minus[k] + self.MAX_BOND * (e_T[i + 1] - chi_T[(i, k)] + 1), f"milp-2LMH-(46)-2-{k}-{i}-2_SG{self.name}"

		MILP += pulp.lpSum(beta_CT[i] for i in range(1, t_T + 1)) == pulp.lpSum(beta_plus[k] for k in (I_ge_one + I_ge_two)), f"milp-2LMH-(46)-3_SG{self.name}"
		MILP += pulp.lpSum(beta_TC[i] for i in range(1, t_T + 1)) == pulp.lpSum(beta_minus[k] for k in (I_ge_one + I_ge_two)), f"milp-2LMH-(46)-4_SG{self.name}"

		# -------- Constraint (47) --------
		for c in range(1, t_C_tilde + 1):
			for i in range(1, t_F + 1):
				MILP += beta_CF[i] >= beta_in[c] - self.MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-2LMH-(47)-1-{c}-{i}-1_SG{self.name}"
				MILP += beta_CF[i] <= beta_in[c] + self.MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-2LMH-(47)-1-{c}-{i}-2_SG{self.name} "
		for c in range(t_C_tilde + 1, c_F + 1):
			for i in range(1, t_F + 1):
				MILP += beta_TF[i] >= beta_in[c] - self.MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-2LMH-(47)-2-{c}-{i}-1_SG{self.name}"
				MILP += beta_TF[i] <= beta_in[c] + self.MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-2LMH-(47)-2-{c}-{i}-2_SG{self.name}"

		MILP += pulp.lpSum(beta_CF[i] for i in range(1, t_F + 1)) == pulp.lpSum(beta_in[c] for c in range(1, t_C_tilde + 1)), f"milp-2LMH-(47)-3_SG{self.name}" 
		MILP += pulp.lpSum(beta_TF[i] for i in range(1, t_F + 1)) == pulp.lpSum(beta_in[c] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-2LMH-(47)-4_SG{self.name}" 

		# -------- Constraint (48) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(delta_alpha_C[(i, atom)] for atom in Lambda_int) == 1, \
					f"milp-2LMH-(48)-C-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Lambda_int[atom] * delta_alpha_C[(i, atom)]
						for atom in Lambda_int) == alpha_C[i], f"milp-2LMH-(48)-C-2-{i}_SG{self.name}"

		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(delta_alpha_T[(i, atom)] for atom in Lambda_int) == v_T[i], \
					f"milp-(48)-T-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Lambda_int[atom] * delta_alpha_T[(i, atom)]
						for atom in Lambda_int) == alpha_T[i], f"milp-2LMH-(48)-T-2-{i}_SG{self.name}"

		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(delta_alpha_F[(i, atom)] for atom in Lambda_int) == v_F[i], \
					f"milp-(48)-F-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Lambda_int[atom] * delta_alpha_F[(i, atom)]
						for atom in Lambda_int) == alpha_F[i], f"milp-2LMH-(48)-F-2-{i}_SG{self.name}"

		# -------- Constraint (49) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(alpha_r[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] 
						for psi in F_C[i]) == alpha_C[i], f"milp-2LMH-(49)-C-{i}_SG{self.name}"
		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(alpha_r[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] 
						for psi in F_T[i]) == alpha_T[i], f"milp-2LMH-(49)-T-{i}_SG{self.name}"
		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(alpha_r[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] 
						for psi in F_F[i]) == alpha_F[i], f"milp-2LMH-(49)-F-{i}_SG{self.name}"

		# -------- Constraint (50) --------
		for i in range(1, t_C_tilde + 1):
			MILP += pulp.lpSum(beta_C[j] for j in range(1, len(E_C) + 1) 
						if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
					pulp.lpSum(beta_plus[k] for k in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
					pulp.lpSum(beta_minus[k] for k in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
					beta_in[i] + beta_ex_C[i] - eledeg_C[i] == \
					pulp.lpSum(val[atom] * delta_alpha_C[(i, atom)] for atom in Lambda_int), f"milp-2LMH-(50)-{i}_SG{self.name}"

		# -------- Constraint (51) --------
		for i in range(t_C_tilde + 1, t_C + 1):
			MILP += pulp.lpSum(beta_C[j] for j in range(1, len(E_C) + 1) 
						if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
					pulp.lpSum(beta_plus[c] for c in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
					pulp.lpSum(beta_minus[c] for c in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
					beta_ex_C[i] - eledeg_C[i] == \
					pulp.lpSum(val[atom] * delta_alpha_C[(i, atom)] for atom in Lambda_int), f"milp-2LMH-(51)-{i}_SG{self.name}" 

		# -------- Constraint (52) --------
		MILP += beta_T[1] == 0, f"milp-2LMH-(52)-T1_SG{self.name}"
		MILP += beta_T[t_T + 1] == 0, f"milp-2LMH-(52)-T2_SG{self.name}"
		for i in range(1, t_T + 1):
			MILP += beta_T[i] + beta_T[i + 1] + beta_ex_T[i] + \
					beta_CT[i] + beta_TC[i] + beta_in[t_C_tilde + i] - eledeg_T[i] == \
					pulp.lpSum(val[atom] * delta_alpha_T[(i, atom)] for atom in Lambda_int), f"milp-2LMH-(52)-{i}_SG{self.name}"

		# -------- Constraint (53) --------
		MILP += beta_F[1] == 0, f"milp-2LMH-(53)-F1_SG{self.name}"
		MILP += beta_F[t_F + 1] == 0, f"milp-2LMH-(53)-F2_SG{self.name}"
		for i in range(1, t_F + 1):
			MILP += beta_F[i] + beta_F[i + 1] + beta_ex_F[i] + \
					beta_CF[i] + beta_TF[i] - eledeg_F[i] == \
					pulp.lpSum(val[atom] * delta_alpha_F[(i, atom)] for atom in Lambda_int), f"milp-2LMH-(53)-{i}_SG{self.name}"

		# -------- Constraint (54) --------
		for atom in Lambda_int:
			MILP += pulp.lpSum(delta_alpha_C[(i, atom)] for i in range(1, t_C + 1)) == na_C[atom], f"milp-2LMH-(54)-C-{atom}_SG{self.name}"
			MILP += pulp.lpSum(delta_alpha_T[(i, atom)] for i in range(1, t_T + 1)) == na_T[atom], f"milp-2LMH-(54)-T-{atom}_SG{self.name}"
			MILP += pulp.lpSum(delta_alpha_F[(i, atom)] for i in range(1, t_F + 1)) == na_F[atom], f"milp-2LMH-(54)-F-{atom}_SG{self.name}"

		# -------- Constraint (55) --------
		for atom in Lambda_ex:
			MILP += pulp.lpSum(na_alpha_ex[atom][Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
						for i in range(1, t_C + 1) for psi in F_C[i]) == na_ex_C[atom], f"milp-2LMH-(55)-C-{atom}_SG{self.name}"
			MILP += pulp.lpSum(na_alpha_ex[atom][Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
						for i in range(1, t_T + 1) for psi in F_T[i]) == na_ex_T[atom], f"milp-2LMH-(55)-T-{atom}_SG{self.name}"
			MILP += pulp.lpSum(na_alpha_ex[atom][Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
						for i in range(1, t_F + 1) for psi in F_F[i]) == na_ex_F[atom], f"milp-2LMH-(55)-F-{atom}_SG{self.name}"

		# -------- Constraint (56) --------
		for atom in Lambda_int:
			MILP += na_C[atom] + na_T[atom] + na_F[atom] == na_int[atom], f"milp-2LMH-(56)-1-{atom}_SG{self.name}"
		for atom in Lambda_ex:
			MILP += na_ex_C[atom] + na_ex_T[atom] + na_ex_F[atom] == na_ex[atom], f"milp-2LMH-(56)-2-{atom}_SG{self.name}"
		for atom in Lambda:
			if atom in Lambda_int and atom in Lambda_ex:
				MILP += na_int[atom] + na_ex[atom] == na[atom], f"milp-2LMH-(56)-3-{atom}_SG{self.name}"
			elif atom in Lambda_int:
				MILP += na_int[atom] == na[atom], f"milp-2LMH-(56)-4-{atom}_SG{self.name}"
			elif atom in Lambda_ex:
				MILP += na_ex[atom] == na[atom], f"milp-2LMH-(56)-5-{atom}_SG{self.name}"

		# -------- Constraint (57) --------
		MILP += pulp.lpSum(mass[atom] * na[atom] for atom in Lambda) == MASS, f"milp-2LMH-(57)_SG{self.name}"

		# -------- Constraint (58) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(delta_alpha_C[(i, atom)] for atom in Lambda_star[i]) == 1, f"milp-2LMH-(58)-{i}_SG{self.name}"

		return MILP

	def prepare_variables_number_of_bounds(self,
		k_C, t_T, bd_T
	): 
		# Define bd_T
		bd_T_temp = {(k, i, m): pulp.LpVariable(f"bd_T({k},{i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for k in range(1, k_C + 1) for i in range(2, t_T + 1) for m in {2, 3}}
		bd_T.update(bd_T_temp)

		return bd_T

	def add_constraints_number_of_bounds(self,
		# Model
		MILP,
		# Constants
		t_T, k_C, E_C, I_equal_one, I_zero_one, bd2_LB, bd2_UB, bd3_LB, bd3_UB,
		# Binary Variables
		chi_T, delta_beta_C, delta_beta_T, delta_beta_plus, delta_beta_minus, bd_T
		# Integer Variables
	):
		# -------- Constraint (59) --------
		for i in (I_equal_one + I_zero_one):
			MILP += delta_beta_C[(i, 2)] >= bd2_LB[E_C[i]], f"milp-2LMH-(59)-{i}-2-1_SG{self.name}"
			MILP += delta_beta_C[(i, 2)] <= bd2_UB[E_C[i]], f"milp-2LMH-(59)-{i}-2-2_SG{self.name}"
			MILP += delta_beta_C[(i, 3)] >= bd3_LB[E_C[i]], f"milp-2LMH-(59)-{i}-3-1_SG{self.name}"
			MILP += delta_beta_C[(i, 3)] <= bd3_UB[E_C[i]], f"milp-2LMH-(59)-{i}-3-2_SG{self.name}"

		# -------- Constraint (60) --------
		for k in range(1, k_C + 1):
			for i in range(2, t_T + 1):
				for m in {2, 3}:
					MILP += bd_T[(k, i, m)] >= \
						delta_beta_T[(i, m)] + chi_T[(i, k)] - 1, f"milp-2LMH-(60)-{k}-{i}-{m}_SG{self.name}"

		# -------- Constraint (61) --------
		for m in {2, 3}:
			MILP += pulp.lpSum(delta_beta_T[(j, m)] for j in range(2, t_T + 1)) >= \
					pulp.lpSum(bd_T[(k, i, m)] for k in range(1, k_C + 1) 
									for i in range(2, t_T + 1)), f"milp-2LMH-(61)-{m}_SG{self.name}"

		# -------- Constraint (62) --------
		for k in range(1, k_C + 1):
			MILP += pulp.lpSum(bd_T[(k, i, 2)] for i in range(2, t_T + 1)) + \
					delta_beta_plus[(k, 2)] + delta_beta_minus[(k, 2)] >= bd2_LB[E_C[k]], f"milp-2LMH-(62)-{k}-2-1_SG{self.name}"
			MILP += pulp.lpSum(bd_T[(k, i, 2)] for i in range(2, t_T + 1)) + \
					delta_beta_plus[(k, 2)] + delta_beta_minus[(k, 2)] <= bd2_UB[E_C[k]], f"milp-2LMH-(62)-{k}-2-2_SG{self.name}"
			MILP += pulp.lpSum(bd_T[(k, i, 3)] for i in range(2, t_T + 1)) + \
					delta_beta_plus[(k, 3)] + delta_beta_minus[(k, 3)] >= bd3_LB[E_C[k]], f"milp-2LMH-(62)-{k}-3-1_SG{self.name}"
			MILP += pulp.lpSum(bd_T[(k, i, 3)] for i in range(2, t_T + 1)) + \
					delta_beta_plus[(k, 3)] + delta_beta_minus[(k, 3)] <= bd3_UB[E_C[k]], f"milp-2LMH-(62)-{k}-3-2_SG{self.name}"

		return MILP

	def prepare_variables_adjacency_configuration(self,
		t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F, n_T, n_C, n_F, delta_i, rho,
		Gamma_int_ac, Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC,
		Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, ac_LB_int, ac_UB_int, Lambda_int, Lambda_ex, MAX_CODE_int, MAX_CODE_ex,
		Gamma_ac_lnk, ac_LB_lnk, ac_UB_lnk
	):
		# Define ac_int
		ac_int = {nu: pulp.LpVariable(f"ac_int{nu}_SG{self.name}", ac_LB_int[nu], ac_UB_int[nu], cat=pulp.LpInteger)
				for nu in Gamma_int_ac}

		# Define ac_C, ac_T, ac_F
		ac_C = {nu: pulp.LpVariable(f"ac_C({nu})_SG{self.name}", 0, m_C, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_C}
		ac_T = {nu: pulp.LpVariable(f"ac_T({nu})_SG{self.name}", 0, t_T, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_T}
		ac_F = {nu: pulp.LpVariable(f"ac_F({nu})_SG{self.name}", 0, t_F, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_F}

		# Define ac_CT, ac_TC
		ac_CT = {nu: pulp.LpVariable(f"ac_CT({nu})_SG{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
				for nu in Gamma_tilde_ac_CT}
		ac_TC = {nu: pulp.LpVariable(f"ac_TC({nu})_SG{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
				for nu in Gamma_tilde_ac_TC} 

		# Define ac_CF, ac_TF
		ac_CF = {nu: pulp.LpVariable(f"ac_CF({nu})_SG{self.name}", 0, t_C_tilde, cat=pulp.LpInteger)
				for nu in Gamma_tilde_ac_CF}
		ac_TF = {nu: pulp.LpVariable(f"ac_TF({nu})_SG{self.name}", 0, t_T, cat=pulp.LpInteger)
				for nu in Gamma_tilde_ac_TF}

		# Define delta_ac_C
		delta_ac_C = {(i, nu): pulp.LpVariable(f"delta_ac_C({i},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(k_C_tilde + 1, m_C + 1) for nu in Gamma_tilde_ac_C}

		# Define delta_ac_T
		delta_ac_T = {(i, nu): pulp.LpVariable(f"delta_ac_T({i},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(2, t_T + 1) for nu in Gamma_tilde_ac_T}
		
		# Define delta_ac_F
		delta_ac_F = {(i, nu): pulp.LpVariable(f"delta_ac_F({i},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(2, t_F + 1) for nu in Gamma_tilde_ac_F}

		# Define delta_ac_CT, delta_ac_TC
		delta_ac_CT = {(k, nu): pulp.LpVariable(f"delta_ac_CT({k},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for k in range(1, k_C + 1) for nu in Gamma_tilde_ac_CT}
		delta_ac_TC = {(k, nu): pulp.LpVariable(f"delta_ac_TC({k},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for k in range(1, k_C + 1) for nu in Gamma_tilde_ac_TC} 

		# Define delta_ac_CF, delta_ac_TF
		delta_ac_CF = {(c, nu): pulp.LpVariable(f"delta_ac_CF({c},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for c in range(1, t_C_tilde + 1) for nu in Gamma_tilde_ac_CF}
		delta_ac_TF = {(i, nu): pulp.LpVariable(f"delta_ac_TF({i},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(1, t_T + 1) for nu in Gamma_tilde_ac_TF}

		# Define alpha_CT, alpha_TC
		alpha_CT = {k: pulp.LpVariable(f"alpha_CT({k})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
					for k in range(1, k_C + 1)}
		alpha_TC = {k: pulp.LpVariable(f"alpha_TC({k})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
					for k in range(1, k_C + 1)}
		
		# Define alpha_CF
		alpha_CF = {c: pulp.LpVariable(f"alpha_CF({c})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
					for c in range(1, t_C_tilde + 1)}

		# Define alpha_TF
		alpha_TF = {i: pulp.LpVariable(f"alpha_TF({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}
		
		# Define Delta_ac_C_plus, Delta_ac_C_minus
		Delta_ac_C_plus = {i: pulp.LpVariable(f"Delta_ac_C_plus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(k_C_tilde + 1, m_C + 1)}
		Delta_ac_C_minus = {i: pulp.LpVariable(f"Delta_ac_C_minus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(k_C_tilde + 1, m_C + 1)}
		
		# Define Delta_ac_T_plus, Delta_ac_T_minus
		Delta_ac_T_plus = {i: pulp.LpVariable(f"Delta_ac_T_plus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(2, t_T + 1)}
		Delta_ac_T_minus = {i: pulp.LpVariable(f"Delta_ac_T_minus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(2, t_T + 1)}

		# Define Delta_ac_F_plus, Delta_ac_F_minus
		Delta_ac_F_plus = {i: pulp.LpVariable(f"Delta_ac_F_plus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(2, t_F + 1)}
		Delta_ac_F_minus = {i: pulp.LpVariable(f"Delta_ac_F_minus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(2, t_F + 1)}

		# Define Delta_ac_CT_plus, Delta_ac_CT_minus
		Delta_ac_CT_plus = {k: pulp.LpVariable(f"Delta_ac_CT_plus({k})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}
		Delta_ac_CT_minus = {k: pulp.LpVariable(f"Delta_ac_CT_minus({k})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}

		# Define Delta_ac_TC_plus, Delta_ac_TC_minus
		Delta_ac_TC_plus = {k: pulp.LpVariable(f"Delta_ac_TC_plus({k})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}
		Delta_ac_TC_minus = {k: pulp.LpVariable(f"Delta_ac_TC_minus({k})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}
		
		# Define Delta_ac_CF_plus, Delta_ac_CF_minus
		Delta_ac_CF_plus = {c: pulp.LpVariable(f"Delta_ac_CF_plus({c})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for c in range(1, t_C_tilde + 1)}
		Delta_ac_CF_minus = {c: pulp.LpVariable(f"Delta_ac_CF_minus({c})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for c in range(1, t_C_tilde + 1)}

		# Define Delta_ac_TF_plus, Delta_ac_TF_minus
		Delta_ac_TF_plus = {i: pulp.LpVariable(f"Delta_ac_TF_plus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(1, t_T + 1)}
		Delta_ac_TF_minus = {i: pulp.LpVariable(f"Delta_ac_TF_minus({i})_SG{self.name}", 0, MAX_CODE_int, cat=pulp.LpInteger)
						for i in range(1, t_T + 1)}

		# Define ac_lnk   IDO
		ac_lnk = {nu: pulp.LpVariable(f"ac_lnk({nu})_SF{self.name}", ac_LB_lnk[nu], ac_UB_lnk[nu], cat=pulp.LpInteger)
					for nu in Gamma_ac_lnk}
		
		# Define ac_C_lnk, ac_T_lnk   IDO
		ac_C_lnk = {nu: pulp.LpVariable(f"ac_C_lnk({nu})_SF{self.name}", 0, m_C, cat=pulp.LpInteger)
					for nu in Gamma_ac_lnk}
		ac_T_lnk = {nu: pulp.LpVariable(f"ac_T_lnk({nu})_SF{self.name}", 0, m_C, cat=pulp.LpInteger)
					for nu in Gamma_ac_lnk}

		# Define ac_CT_lnk, ac_TC_lnk   IDO
		ac_CT_lnk = {nu: pulp.LpVariable(f"ac_CT_lnk({nu})_SF{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
					for nu in Gamma_ac_lnk}
		ac_TC_lnk = {nu: pulp.LpVariable(f"ac_TC_lnk({nu})_SF{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
					for nu in Gamma_ac_lnk}
		
		 # Define delta_ac_T_lnk   IDO
		delta_ac_T_lnk = {(i, nu): pulp.LpVariable(f"delta_ac_T_lnk({i},{nu})_SF{self.name}", 0, 1, cat=pulp.LpBinary)
					for i in range(2, t_T+1) for nu in Gamma_ac_lnk}

		return ac_int, ac_C, ac_T, ac_F, ac_CT, ac_TC, ac_CF, ac_TF, \
			delta_ac_C, delta_ac_T, delta_ac_F, \
			alpha_CT, alpha_TC, alpha_CF, alpha_TF, \
			Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, \
			Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus, \
			Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus, Delta_ac_TF_minus, \
			delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF, \
			ac_lnk, ac_C_lnk, ac_T_lnk, ac_CT_lnk, ac_TC_lnk, delta_ac_T_lnk

	def add_constraints_adjacency_configuration(self,
		# Model
		MILP,
		# Constants
		t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, k_C, k_C_tilde, m_C, c_F, tail_C, head_C,
		Gamma_int_ac, Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F, Gamma_tilde_ac_CT,
		Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, Gamma_int_ac_less, Gamma_int_ac_equal,
		ac_LB_int, ac_UB_int, Lambda_int, Lambda_ex, Code_Lambda_int, Code_Lambda_ex, MAX_CODE_int, MAX_CODE_ex,
		rho, delta_i, Gamma_ac_lnk, Gamma_ac_lnk_less, Gamma_ac_lnk_equal,
		# Binary Variables
		delta_alpha_C, delta_alpha_T, delta_alpha_F, delta_beta_C, delta_beta_T, delta_beta_F, 
		delta_beta_plus, delta_beta_minus, delta_beta_in, delta_chi_T, delta_chi_F,
		chi_T, chi_F, delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF,
		e_C, e_T, e_F, v_T, v_F,
		Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus,
		Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus, Delta_ac_CF_plus, Delta_ac_CF_minus,
		Delta_ac_TF_plus, Delta_ac_TF_minus, delta_ac_T_lnk,
		# Integer Variables
		clr_T, alpha_C, alpha_T, alpha_F, alpha_CT, alpha_TC, alpha_CF, alpha_TF, beta_C, beta_T, beta_F, 
		beta_plus, beta_minus, beta_in, ac_C, ac_T, ac_F, ac_CT, ac_TC, ac_CF, ac_TF, ac_int,
		ac_lnk, ac_C_lnk, ac_T_lnk, ac_CT_lnk, ac_TC_lnk, I_lnk, n_lnk
	):
		# -------- Constraint (63) --------
		for nu in Gamma_int_ac:
			if nu not in Gamma_tilde_ac_C:
				MILP += ac_C[nu] == 0, f"milp-2LMH-(63)-1-{nu}_SG{self.name}"
			if nu not in Gamma_tilde_ac_T:
				MILP += ac_T[nu] == 0, f"milp-2LMH-(63)-2-{nu}_SG{self.name}"
			if nu not in Gamma_tilde_ac_F:
				MILP += ac_F[nu] == 0, f"milp-2LMH-(63)-3-{nu}_SG{self.name}"
			if nu not in Gamma_tilde_ac_CT:
				MILP += ac_CT[nu] == 0, f"milp-2LMH-(63)-4-{nu}_SG{self.name}"
			if nu not in Gamma_tilde_ac_TC:
				MILP += ac_TC[nu] == 0, f"milp-2LMH-(63)-5-{nu}_SG{self.name}"
			if nu not in Gamma_tilde_ac_CF:
				MILP += ac_CF[nu] == 0, f"milp-2LMH-(63)-6-{nu}_SG{self.name}"
			if nu not in Gamma_tilde_ac_TF:
				MILP += ac_TF[nu] == 0, f"milp-2LMH-(63)-7-{nu}_SG{self.name}"

		# -------- Constraint (64) --------
		for m in range(1, self.MAX_BOND + 1):
			MILP += pulp.lpSum(ac_C[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
					pulp.lpSum(delta_beta_C[(i, m)] for i in range(k_C_tilde + 1, m_C + 1)), f"milp-2LMH-(64)-1-{m}_SG{self.name}"
			MILP += pulp.lpSum(ac_T[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
					pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)), f"milp-2LMH-(64)-2-{m}_SG{self.name}"
			MILP += pulp.lpSum(ac_F[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
					pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)), f"milp-2LMH-(64)-3-{m}_SG{self.name}"
			MILP += pulp.lpSum(ac_CT[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
					pulp.lpSum(delta_beta_plus[(k, m)] for k in range(1, k_C + 1)), f"milp-2LMH-(64)-4-{m}_SG{self.name}"
			MILP += pulp.lpSum(ac_TC[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
					pulp.lpSum(delta_beta_minus[(k, m)] for k in range(1, k_C + 1)), f"milp-2LMH-(64)-5-{m}_SG{self.name}"
			MILP += pulp.lpSum(ac_CF[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
					pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)), f"milp-2LMH-(64)-6-{m}_SG{self.name}"
			MILP += pulp.lpSum(ac_TF[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
					pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-2LMH-(64)-7-{m}_SG{self.name}"

		# -------- Constraint (65) --------
		for i in range(k_C_tilde + 1, m_C + 1):
			MILP += pulp.lpSum(m * delta_ac_C[(i, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_C) == beta_C[i], f"milp-2LMH-(65)-1-{i}_SG{self.name}"
			MILP += Delta_ac_C_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_C[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_C) == alpha_C[tail_C[i]], f"milp-2LMH-(65)-2-{i}_SG{self.name}"
			MILP += Delta_ac_C_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_C[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_C) == alpha_C[head_C[i]], f"milp-2LMH-(65)-3-{i}_SG{self.name}"
			MILP += Delta_ac_C_plus[i] + Delta_ac_C_minus[i] <= 2 * MAX_CODE_int * (1 - e_C[i]), f"milp-2LMH-(65)-4-{i}_SG{self.name}"
		for nu in Gamma_tilde_ac_C:
			MILP += pulp.lpSum(delta_ac_C[(i, nu)] 
					for i in range(k_C_tilde + 1, m_C + 1)) == ac_C[nu], f"milp-2LMH-(65)-5-{nu}_SG{self.name}"

		# -------- Constraint (65) -------- IDO
		for nu in Gamma_ac_lnk:
			MILP += pulp.lpSum(delta_ac_C[(i, nu)]
				  for i in (set(I_lnk) & set(range(k_C + 1, m_C + 1)))) == ac_C_lnk[nu], f"milp-2LMM-(66)-{nu}_SG{self.name}"

		# -------- Constraint (66) --------
		for i in range(2, t_T + 1):
			MILP += pulp.lpSum(m * delta_ac_T[(i, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_T) == beta_T[i], f"milp-2LMH-(66)-1-{i}_SG{self.name}"
			MILP += Delta_ac_T_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_T[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_T) == alpha_T[i - 1], f"milp-2LMH-(66)-2-{i}_SG{self.name}"
			MILP += Delta_ac_T_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_T[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_T) == alpha_T[i], f"milp-2LMH-(66)-3-{i}_SG{self.name}"
			MILP += Delta_ac_T_plus[i] + Delta_ac_T_minus[i] <= 2 * MAX_CODE_int * (1 - e_T[i]), f"milp-2LMH-(66)-4-{i}_SG{self.name}"
		for nu in Gamma_tilde_ac_T:
			MILP += pulp.lpSum(delta_ac_T[(i, nu)] 
					for i in range(2, t_T + 1)) == ac_T[nu], f"milp-2LMH-(66)-5-{nu}_SG{self.name}"

		# -------- Constraint (67) --------
		for i in range(2, t_F + 1):
			MILP += pulp.lpSum(m * delta_ac_F[(i, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_F) == beta_F[i], f"milp-2LMH-(67)-1-{i}_SG{self.name}"
			MILP += Delta_ac_F_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_F[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[i - 1], f"milp-2LMH-(67)-2-{i}_SG{self.name}"
			MILP += Delta_ac_F_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_F[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[i], f"milp-2LMH-(67)-3-{i}_SG{self.name}"
			MILP += Delta_ac_F_plus[i] + Delta_ac_F_minus[i] <= 2 * MAX_CODE_ex * (1 - e_F[i]), f"milp-2LMH-(67)-4-{i}_SG{self.name}"
		for nu in Gamma_tilde_ac_F:
			MILP += pulp.lpSum(delta_ac_F[(i, nu)] 
					for i in range(2, t_F + 1)) == ac_F[nu], f"milp-2LMH-(67)-5-{nu}_SG{self.name}"

		# -------- Constraint (68) --------  IDO
		for i in range(2, t_T + 1):
			for nu in Gamma_ac_lnk:
				MILP += delta_ac_T[(i, nu)] + pulp.lpSum(chi_T[(i, k)] 
						for k in (set(I_lnk) & set(range(1, k_C + 1)))) >= 2 * delta_ac_T_lnk[(i, nu)], f"milp-2LMH-(68)-1-{i}-{nu}-1_SG{self.name}"
				MILP += delta_ac_T_lnk[(i, nu)] >= delta_ac_T[(i, nu)] + pulp.lpSum(chi_T[(i, k)] 
							for k in (set(I_lnk) & set(range(1, k_C + 1)))) - 1, f"milp-2LMH-(68)-1-{i}-{nu}-2_SG{self.name}" 
		for nu in Gamma_ac_lnk:
			MILP += pulp.lpSum(delta_ac_T_lnk[(i, nu)] for i in range(2, t_T + 1)) == ac_T_lnk[nu], f"milp-2LMH-(68)-2-{nu}_SG{self.name}"
		
		# -------- Constraint (69) --------
		for i in range(2, t_F + 1):
			MILP += pulp.lpSum(m * delta_ac_F[(i, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_F) == beta_F[i], f"milp-2LMH-(69)-1-{i}_SG{self.name}"
			MILP += Delta_ac_F_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_F[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[i - 1], f"milp-2LMH-(69)-2-{i}_SG{self.name}"
			MILP += Delta_ac_F_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_F[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[i], f"milp-2LMH-(69)-3-{i}_SG{self.name}"
			MILP += Delta_ac_F_plus[i] + Delta_ac_F_minus[i] <= 2 * MAX_CODE_ex * (1 - e_F[i]), f"milp-2LMH-(69)-4-{i}_SG{self.name}"
		for nu in Gamma_tilde_ac_F:
			MILP += pulp.lpSum(delta_ac_F[(i, nu)] 
					for i in range(2, t_F + 1)) == ac_F[nu], f"milp-2LMH-(69)-5-{nu}_SG{self.name}"

		# -------- Constraint (70) --------
		for k in range(1, k_C + 1):
			for i in range(1, t_T + 1):
				MILP += alpha_T[i] + MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i]) >= \
						alpha_CT[k], f"milp-2LMH-(70)-1-{k}-{i}_SG{self.name}"
				MILP += alpha_CT[k] >= alpha_T[i] - \
						MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i]), f"milp-2LMH-(70)-2-{k}-{i}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_ac_CT[(k, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_CT) == beta_plus[k], f"milp-2LMH-(70)-3-{k}_SG{self.name}"
			MILP += Delta_ac_CT_plus[k] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_CT[(k, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_CT) == alpha_C[tail_C[k]], f"milp-2LMH-(70)-4-{k}_SG{self.name}"
			MILP += Delta_ac_CT_minus[k] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_CT[(k, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_CT) == alpha_CT[k], f"milp-2LMH-(70)-5-{k}_SG{self.name}"
			MILP += Delta_ac_CT_plus[k] + Delta_ac_CT_minus[k] <= 2 * MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-2LMH-(70)-6-{k}_SG{self.name}"
		for nu in Gamma_tilde_ac_CT:
			MILP += pulp.lpSum(delta_ac_CT[(k, nu)] for k in range(1, k_C + 1)) == ac_CT[nu], f"milp-2LMH-(70)-7-{nu}_SG{self.name}"

		# -------- Constraint (71) --------   IDO
		for nu in Gamma_ac_lnk:
			MILP += pulp.lpSum(delta_ac_CT[(i, nu)] for i in (set(I_lnk) & set(range(1, k_C + 1)))) == ac_CT_lnk[nu], f"milp-2LMM-(71)-{nu}_SG{self.name}"

		# -------- Constraint (72) --------
		for k in range(1, k_C + 1):
			for i in range(1, t_T + 1):
				MILP += alpha_T[i] + MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i + 1]) >= \
						alpha_TC[k], f"milp-2LMH-(72)-1-{k}-{i}_SG{self.name}"
				MILP += alpha_TC[k] >= alpha_T[i] - \
						MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i + 1]), f"milp-2LMH-(72)-2-{k}-{i}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_ac_TC[(k, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_TC) == beta_minus[k], f"milp-2LMH-(72)-3-{k}_SG{self.name}"
			MILP += Delta_ac_TC_plus[k] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_TC[(k, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_TC) == alpha_TC[k], f"milp-2LMH-(72)-4-{k}_SG{self.name}"
			MILP += Delta_ac_TC_minus[k] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_TC[(k, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_TC) == alpha_C[head_C[k]], f"milp-2LMH-(72)-5-{k}_SG{self.name}"
			MILP += Delta_ac_TC_plus[k] + Delta_ac_TC_minus[k] <= 2 * MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-2LMH-(72)-6-{k}_SG{self.name}"
		for nu in Gamma_tilde_ac_TC:
			MILP += pulp.lpSum(delta_ac_TC[(k, nu)] for k in range(1, k_C + 1)) == ac_TC[nu], f"milp-2LMH-(72)-7-{nu}_SG{self.name}"

		# -------- Constraint (73) --------  IDO
		for nu in Gamma_ac_lnk:
			MILP += pulp.lpSum(delta_ac_TC[(i, nu)] for i in (set(I_lnk) & set(range(1, k_C + 1)))) == ac_TC_lnk[nu], f"milp-2LMM-(73)-{nu}_SG{self.name}"

		# -------- Constraint (74) --------
		for c in range(1, t_C_tilde + 1):
			for i in range(1, t_F + 1):
				MILP += alpha_F[i] + MAX_CODE_int * (1 - chi_F[(i, c)] + e_F[i]) >= \
						alpha_CF[c], f"milp-2LMH-(74)-1-{c}-{i}_SG{self.name}"
				MILP += alpha_CF[c] >= alpha_F[i] - \
						MAX_CODE_int * (1 - chi_F[(i, c)] + e_F[i]), f"milp-2LMH-(74)-2-{c}-{i}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_ac_CF[(c, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_CF) == beta_in[c], f"milp-2LMH-(74)-3-{c}_SG{self.name}"
			MILP += Delta_ac_CF_plus[c] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_CF[(c, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_CF) == alpha_C[c], f"milp-2LMH-(74)-4-{c}_SG{self.name}"
			MILP += Delta_ac_CF_minus[c] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_CF[(c, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_CF) == alpha_CF[c], f"milp-2LMH-(74)-5-{c}_SG{self.name}"
			MILP += Delta_ac_CF_plus[c] + Delta_ac_CF_minus[c] <= \
					2 * max(MAX_CODE_int, MAX_CODE_int) * (1 - delta_chi_F[c]), f"milp-2LMH-(74)-6-{c}_SG{self.name}"
		for nu in Gamma_tilde_ac_CF:
			MILP += pulp.lpSum(delta_ac_CF[(c, nu)] for c in range(1, t_C_tilde + 1)) == ac_CF[nu], f"milp-2LMH-(74)-7-{nu}_SG{self.name}"

		# -------- Constraint (75) --------
		for i in range(1, t_T + 1):
			for j in range(1, t_F + 1):
				MILP += alpha_F[j] + MAX_CODE_int * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]) >= \
						alpha_TF[i], f"milp-2LMH-(75)-1-{i}-{j}_SG{self.name}"
				MILP += alpha_TF[i] >= alpha_F[j] - \
						MAX_CODE_int * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]), f"milp-2LMH-(75)-2-{i}-{j}_SG{self.name}"
			MILP += pulp.lpSum(m * delta_ac_TF[(i, (a, b, m))] 
				for (a, b, m) in Gamma_tilde_ac_TF) == beta_in[i + t_C_tilde], f"milp-2LMH-(75)-3-{i}_SG{self.name}"
			MILP += Delta_ac_TF_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_TF[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_TF) == alpha_T[i], f"milp-2LMH-(75)-4-{i}_SG{self.name}"
			MILP += Delta_ac_TF_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_TF[(i, (a, b, m))]
						for (a, b, m) in Gamma_tilde_ac_TF) == alpha_TF[i], f"milp-2LMH-(75)-5-{i}_SG{self.name}"
			MILP += Delta_ac_TF_plus[i] + Delta_ac_TF_minus[i] <= \
					2 * max(MAX_CODE_int, MAX_CODE_int) * (1 - delta_chi_F[i + t_C_tilde]), f"milp-2LMH-(75)-6-{i}_SG{self.name}"
		for nu in Gamma_tilde_ac_TF:
			MILP += pulp.lpSum(delta_ac_TF[(i, nu)] for i in range(1, t_T + 1)) == ac_TF[nu], f"milp-2LMH-(75)-7-{nu}_SG{self.name}"

		# -------- Constraint (76) --------
		for (a, b, m) in Gamma_int_ac_less:
			MILP += ac_C[(a, b, m)] + ac_C[(b, a, m)] + \
					ac_T[(a, b, m)] + ac_T[(b, a, m)] + \
					ac_F[(a, b, m)] + ac_F[(b, a, m)] + \
					ac_CT[(a, b, m)] + ac_CT[(b, a, m)] + \
					ac_TC[(a, b, m)] + ac_TC[(b, a, m)] + \
					ac_CF[(a, b, m)] + ac_CF[(b, a, m)] + \
					ac_TF[(a, b, m)] + ac_TF[(b, a, m)] == ac_int[(a, b, m)], f"milp-2LMH-(76)-1-{(a, b, m)}_SG{self.name}"
		for nu in Gamma_int_ac_equal:
			MILP += ac_C[nu] + ac_T[nu] + ac_F[nu] + \
					ac_CT[nu] + ac_TC[nu] + ac_CF[nu] + ac_TF[nu] == ac_int[nu], f"milp-2LMH-(76)-2-{nu}_SG{self.name}"

		# -------- Constraint (77) --------  IDO
		for (a, b, m) in Gamma_ac_lnk_less:
			MILP += ac_C_lnk[(a, b, m)] + ac_C_lnk[(b, a, m)] + \
					ac_T_lnk[(a, b, m)] + ac_T_lnk[(b, a, m)] + \
					ac_CT_lnk[(a, b, m)] + ac_CT_lnk[(b, a, m)] + \
					ac_TC_lnk[(a, b, m)] + ac_TC_lnk[(b, a, m)] == ac_lnk[(a, b, m)], f"milp-2LMM-(77)-1-{(a, b, m)}_SG{self.name}"
		for nu in Gamma_ac_lnk_equal:
			MILP += ac_C_lnk[nu] + ac_T_lnk[nu] + \
					ac_CT_lnk[nu] + ac_TC_lnk[nu] == ac_lnk[nu], f"milp-2LMM-(77)-2-{nu}_SG{self.name}"

		MILP += pulp.lpSum(ac_lnk[gamma] for gamma in Gamma_ac_lnk_less) + \
				pulp.lpSum(ac_lnk[gamma] for gamma in Gamma_ac_lnk_equal) == n_lnk, f"milp-2LMM-polymer-lnk_ac_SG{self.name}"

		return MILP

	def prepare_variables_chemical_symbols(self,
		t_C, t_T, t_F, n_C, n_T, n_F, delta_i, n_star, n_LB_int, n_UB_int, ns_LB_int, ns_UB_int, Lambda_dg_int, epsilon_dg
	):
		# Define ns_int
		ns_int = {mu: pulp.LpVariable(f"ns_int({mu})_SG{self.name}", ns_LB_int[mu], ns_UB_int[mu], cat=pulp.LpInteger) for mu in Lambda_dg_int}

		# Define delta_ns_C, delta_ns_T, delta_ns_F
		delta_ns_C = {(i, mu): pulp.LpVariable(f"delta_ns_C({i},{mu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_C + 1) for mu in (Lambda_dg_int + [epsilon_dg])}
		delta_ns_T = {(i, mu): pulp.LpVariable(f"delta_ns_T({i},{mu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_T + 1) for mu in (Lambda_dg_int + [epsilon_dg])}
		delta_ns_F = {(i, mu): pulp.LpVariable(f"delta_ns_F({i},{mu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
						for i in range(1, t_F + 1) for mu in (Lambda_dg_int + [epsilon_dg])}
		
		return ns_int, delta_ns_C, delta_ns_T, delta_ns_F

	def add_constraints_chemical_symbols(self,
		# Model
		MILP,
		# Constants
		t_C, t_T, t_F, n_C, n_T, n_F, Lambda_dg_int, Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F,
		Code_Lambda, Code_Lambda_int, Code_Lambda_ex, epsilon_dg, delta_i, rho,
		# Binary Variables
		delta_ns_C, delta_ns_T, delta_ns_F,
		# Integer Variables
		alpha_C, alpha_T, alpha_F, deg_C, deg_T, deg_F, ns_int
	):
		# -------- Constraint (78) --------
		for i in range(1, t_C + 1):
			MILP += pulp.lpSum(delta_ns_C[(i, mu)] 
					for mu in (Lambda_tilde_dg_C + [epsilon_dg])) == 1, f"milp-2LMH-(78)-C-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Lambda_int[a] * delta_ns_C[(i, (a, d))] 
					for (a, d) in Lambda_tilde_dg_C) == alpha_C[i], f"milp-2LMH-(78)-C-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_ns_C[(i, (a, d))]
					for (a, d) in Lambda_tilde_dg_C) == deg_C[i], f"milp-2LMH-(78)-C-3-{i}_SG{self.name}"
		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(delta_ns_T[(i, mu)] 
					for mu in (Lambda_tilde_dg_T + [epsilon_dg])) == 1, f"milp-2LMH-(78)-T-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Lambda_int[a] * delta_ns_T[(i, (a, d))] 
					for (a, d) in Lambda_tilde_dg_T) == alpha_T[i], f"milp-2LMH-(78)-T-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_ns_T[(i, (a, d))]
					for (a, d) in Lambda_tilde_dg_T) == deg_T[i], f"milp-2LMH-(78)-T-3-{i}_SG{self.name}"
		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(delta_ns_F[(i, mu)] 
					for mu in (Lambda_tilde_dg_F + [epsilon_dg])) == 1, f"milp-2LMH-(78)-F-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Lambda_int[a] * delta_ns_F[(i, (a, d))] 
					for (a, d) in Lambda_tilde_dg_F) == alpha_F[i], f"milp-2LMH-(78)-F-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(d * delta_ns_F[(i, (a, d))]
					for (a, d) in Lambda_tilde_dg_F) == deg_F[i], f"milp-2LMH-(78)-F-3-{i}_SG{self.name}"
		
		# -------- Constraint (79) --------
		for mu in Lambda_dg_int:
			MILP += pulp.lpSum(delta_ns_C[(i, mu)] for i in range(1, t_C + 1)) + \
					pulp.lpSum(delta_ns_T[(i, mu)] for i in range(1, t_T + 1)) + \
					pulp.lpSum(delta_ns_F[(i, mu)] for i in range(1, t_F + 1)) == ns_int[mu], f"milp-2LMH-(79)-{mu}_SG{self.name}"

		return MILP

	def prepare_variables_edge_configuration(self,
		t_C, t_C_tilde, t_T, t_F, m_C, k_C, k_C_tilde, c_F, n_T, n_C, n_F, delta_i, rho,
		Gamma_int, Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_F, Gamma_tilde_ec_CT, Gamma_tilde_ec_TC,
		Gamma_tilde_ec_CF, Gamma_tilde_ec_TF, ec_LB_int, ec_UB_int,
		ec_LB_lnk, ec_UB_lnk, Gamma_lnk, Lambda_dg_int, Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great
	):
		# Define ec_int
		ec_int = {gamma: pulp.LpVariable(f"ec_int{gamma}_SG{self.name}", ec_LB_int[gamma], ec_UB_int[gamma], cat=pulp.LpInteger)
				for gamma in Gamma_int}

		# Define ec_C, ec_T, ec_F
		ec_C = {gamma: pulp.LpVariable(f"ec_C({gamma})_SG{self.name}", 0, m_C, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_C}
		ec_T = {gamma: pulp.LpVariable(f"ec_T({gamma})_SG{self.name}", 0, t_T, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_T}
		ec_F = {gamma: pulp.LpVariable(f"ec_F({gamma})_SG{self.name}", 0, t_F, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_F}

		# Define ec_CT, ec_TC
		ec_CT = {gamma: pulp.LpVariable(f"ec_CT({gamma})_SG{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
				for gamma in Gamma_tilde_ec_CT}
		ec_TC = {gamma: pulp.LpVariable(f"ec_TC({gamma})_SG{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
				for gamma in Gamma_tilde_ec_TC} 

		# Define ec_CF, ec_TF
		ec_CF = {gamma: pulp.LpVariable(f"ec_CF({gamma})_SG{self.name}", 0, t_C_tilde, cat=pulp.LpInteger)
				for gamma in Gamma_tilde_ec_CF}
		ec_TF = {gamma: pulp.LpVariable(f"ec_TF({gamma})_SG{self.name}", 0, t_T, cat=pulp.LpInteger)
				for gamma in Gamma_tilde_ec_TF}

		# Define delta_ec_C
		delta_ec_C = {(i, gamma): pulp.LpVariable(f"delta_ec_C({i},{gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(k_C_tilde + 1, m_C + 1) for gamma in Gamma_tilde_ec_C}

		# Define delta_ec_T
		delta_ec_T = {(i, gamma): pulp.LpVariable(f"delta_ec_T({i},{gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(2, t_T + 1) for gamma in Gamma_tilde_ec_T}
		
		# Define delta_ec_F
		delta_ec_F = {(i, gamma): pulp.LpVariable(f"delta_ec_F({i},{gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(2, t_F + 1) for gamma in Gamma_tilde_ec_F}

		# Define delta_ec_CT, delta_ec_TC
		delta_ec_CT = {(k, gamma): pulp.LpVariable(f"delta_ec_CT({k},{gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for k in range(1, k_C + 1) for gamma in Gamma_tilde_ec_CT}
		delta_ec_TC = {(k, gamma): pulp.LpVariable(f"delta_ec_TC({k},{gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for k in range(1, k_C + 1) for gamma in Gamma_tilde_ec_TC} 

		# Define delta_ec_CF, delta_ec_TF
		delta_ec_CF = {(c, gamma): pulp.LpVariable(f"delta_ec_CF({c},{gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for c in range(1, t_C_tilde + 1) for gamma in Gamma_tilde_ec_CF}
		delta_ec_TF = {(i, gamma): pulp.LpVariable(f"delta_ec_TF({i},{gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
				for i in range(1, t_T + 1) for gamma in Gamma_tilde_ec_TF}
		
		# Define deg_T_CT, deg_T_TC
		deg_T_CT = {k: pulp.LpVariable(f"deg_T_CT({k})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for k in range(1, k_C + 1)}
		deg_T_TC = {k: pulp.LpVariable(f"deg_T_TC({k})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for k in range(1, k_C + 1)}

		# Define deg_F_CF
		deg_F_CF = {c: pulp.LpVariable(f"deg_F_CF({c})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for c in range(1, t_C_tilde + 1)}

		# Define deg_F_TF
		deg_F_TF = {i: pulp.LpVariable(f"deg_F_TF({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
					for i in range(1, t_T + 1)}

		# Define Delta_ec_C_plus, Delta_ec_C_minus
		Delta_ec_C_plus = {i: pulp.LpVariable(f"Delta_ec_C_plus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(k_C_tilde + 1, m_C + 1)}
		Delta_ec_C_minus = {i: pulp.LpVariable(f"Delta_ec_C_minus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(k_C_tilde + 1, m_C + 1)}
		
		# Define Delta_ec_T_plus, Delta_ec_T_minus
		Delta_ec_T_plus = {i: pulp.LpVariable(f"Delta_ec_T_plus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(2, t_T + 1)}
		Delta_ec_T_minus = {i: pulp.LpVariable(f"Delta_ec_T_minus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(2, t_T + 1)}

		# Define Delta_ec_F_plus, Delta_ec_F_minus
		Delta_ec_F_plus = {i: pulp.LpVariable(f"Delta_ec_F_plus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(2, t_F + 1)}
		Delta_ec_F_minus = {i: pulp.LpVariable(f"Delta_ec_F_minus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(2, t_F + 1)}
		
		# Define Delta_ec_CT_plus, Delta_ec_CT_minus
		Delta_ec_CT_plus = {k: pulp.LpVariable(f"Delta_ec_CT_plus({k})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}
		Delta_ec_CT_minus = {k: pulp.LpVariable(f"Delta_ec_CT_minus({k})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}

		# Define Delta_ec_TC_plus, Delta_ec_TC_minus
		Delta_ec_TC_plus = {k: pulp.LpVariable(f"Delta_ec_TC_plus({k})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}
		Delta_ec_TC_minus = {k: pulp.LpVariable(f"Delta_ec_TC_minus({k})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for k in range(1, k_C + 1)}
		
		# Define Delta_ec_CF_plus, Delta_ec_CF_minus
		Delta_ec_CF_plus = {c: pulp.LpVariable(f"Delta_ec_CF_plus({c})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for c in range(1, t_C_tilde + 1)}
		Delta_ec_CF_minus = {c: pulp.LpVariable(f"Delta_ec_CF_minus({c})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for c in range(1, t_C_tilde + 1)}

		# Define Delta_ec_TF_plus, Delta_ec_TF_minus
		Delta_ec_TF_plus = {i: pulp.LpVariable(f"Delta_ec_TF_plus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(1, t_T + 1)}
		Delta_ec_TF_minus = {i: pulp.LpVariable(f"Delta_ec_TF_minus({i})_SG{self.name}", 0, self.MAX_VAL, cat=pulp.LpInteger)
						for i in range(1, t_T + 1)}

		# Define ec_lnk  IDO
		ec_lnk = {gamma: pulp.LpVariable(f"ec_lnk({gamma})_SG{self.name}", ec_LB_lnk[gamma], ec_UB_lnk[gamma], cat=pulp.LpInteger) for gamma in Gamma_lnk}

		# Define ec_C_lnk, ec_T_lnk IDO
		ec_C_lnk = {gamma: pulp.LpVariable(f"ec_C_lnk({gamma})_SG{self.name}", 0, m_C, cat=pulp.LpInteger) for gamma in Gamma_lnk}
		ec_T_lnk = {gamma: pulp.LpVariable(f"ec_T_lnk({gamma})_SG{self.name}", 0, m_C, cat=pulp.LpInteger) for gamma in Gamma_lnk}


		# Define ec_CT_lnk, ec_TC_lnk   IDO
		ec_CT_lnk = {nu: pulp.LpVariable(f"ec_CT_lnk({nu})_SG{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
					for nu in Gamma_lnk}
		ec_TC_lnk = {nu: pulp.LpVariable(f"ec_TC_lnk({nu})_SG{self.name}", 0, min(k_C, t_T), cat=pulp.LpInteger)
					for nu in Gamma_lnk}
		
		 # Define delta_ec_T_lnk   IDO
		delta_ec_T_lnk = {(i, nu): pulp.LpVariable(f"delta_ec_T_lnk({i},{nu})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for i in range(2, t_T+1) for nu in Gamma_lnk}

		# Define delta_cnt
		delta_cnt = {gamma: pulp.LpVariable(f"delta_cnt({gamma})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for gamma in (Gamma_cnt_less + Gamma_cnt_equal + Gamma_cnt_great)}
		
		return ec_int, ec_C, ec_T, ec_F, ec_CT, ec_TC, ec_CF, ec_TF, \
			delta_ec_C, delta_ec_T, delta_ec_F, \
			delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF, \
			deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF, \
			Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, Delta_ec_F_plus, Delta_ec_F_minus, \
			Delta_ec_CT_plus, Delta_ec_CT_minus, Delta_ec_TC_plus, Delta_ec_TC_minus, \
			Delta_ec_CF_plus, Delta_ec_CF_minus, Delta_ec_TF_plus, Delta_ec_TF_minus, \
			ec_lnk, ec_C_lnk, ec_T_lnk, ec_CT_lnk, ec_TC_lnk, delta_ec_T_lnk, delta_cnt

	def add_constraints_edge_configuration(self,
		# Model
		MILP,
		# Constants
		t_C, t_C_tilde, t_T, t_F, n_C, n_T, n_F, k_C, k_C_tilde, m_C, c_F, tail_C, head_C,
		Gamma_int, Gamma_int_less, Gamma_int_equal, Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_F,
		Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF, 
		Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_F, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF,
		rho, delta_i, Code_Gamma_ac_int, Code_Gamma_ec_int, 
		Gamma_lnk, Gamma_lnk_less, Gamma_lnk_equal, Lambda_dg_int, ns_LB_cnt, ns_UB_cnt, I_lnk,
		Gamma_cnt_less, Gamma_cnt_equal, Gamma_cnt_great,
		# Binary Variables
		delta_dg_C, delta_dg_T, delta_dg_F, delta_chi_T, delta_chi_F, chi_T, chi_F, delta_beta_C, delta_beta_T, delta_beta_F,
		delta_beta_plus, delta_beta_minus, delta_beta_in, delta_ac_C, delta_ac_T, delta_ac_F,
		delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF, delta_ec_C, delta_ec_T, delta_ec_F,
		delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF,
		e_C, e_T, e_F, v_T, v_F,
		delta_ec_T_lnk, delta_cnt,
		# Integer Variables
		clr_T, deg_C, deg_T, deg_F, deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF, ec_C, ec_T, ec_F, ec_CT, ec_TC, ec_CF, ec_TF, ec_int,
		Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, Delta_ec_F_plus, Delta_ec_F_minus,
		Delta_ec_CT_plus, Delta_ec_CT_minus, Delta_ec_TC_plus, Delta_ec_TC_minus,
		Delta_ec_CF_plus, Delta_ec_CF_minus, Delta_ec_TF_plus, Delta_ec_TF_minus,
		ec_lnk, ec_T_lnk, ec_C_lnk, ec_CT_lnk, ec_TC_lnk, n_lnk
	):
		# -------- Constraint (80) --------
		for gamma in Gamma_int:
			if gamma not in Gamma_tilde_ec_C:
				MILP += ec_C[gamma] == 0, f"milp-2LMH-(80)-1-{gamma}_SG{self.name}"
			if gamma not in Gamma_tilde_ec_T:
				MILP += ec_T[gamma] == 0, f"milp-2LMH-(80)-2-{gamma}_SG{self.name}"
			if gamma not in Gamma_tilde_ec_F:
				MILP += ec_F[gamma] == 0, f"milp-2LMH-(80)-3-{gamma}_SG{self.name}"
			if gamma not in Gamma_tilde_ec_CT:
				MILP += ec_CT[gamma] == 0, f"milp-2LMH-(80)-4-{gamma}_SG{self.name}"
			if gamma not in Gamma_tilde_ec_TC:
				MILP += ec_TC[gamma] == 0, f"milp-2LMH-(80)-5-{gamma}_SG{self.name}"
			if gamma not in Gamma_tilde_ec_CF:
				MILP += ec_CF[gamma] == 0, f"milp-2LMH-(80)-6-{gamma}_SG{self.name}"
			if gamma not in Gamma_tilde_ec_TF:
				MILP += ec_TF[gamma] == 0, f"milp-2LMH-(80)-7-{gamma}_SG{self.name}"

		# -------- Constraint (81) --------
		for m in range(1, self.MAX_BOND + 1):
			MILP += pulp.lpSum(ec_C[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
					pulp.lpSum(delta_beta_C[(i, m)] for i in range(k_C_tilde + 1, m_C + 1)), f"milp-2LMH-(81)-1-{m}_SG{self.name}"
			MILP += pulp.lpSum(ec_T[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
					pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)), f"milp-2LMH-(81)-2-{m}_SG{self.name}"
			MILP += pulp.lpSum(ec_F[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
					pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)), f"milp-2LMH-(81)-3-{m}_SG{self.name}"
			MILP += pulp.lpSum(ec_CT[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
					pulp.lpSum(delta_beta_plus[(k, m)] for k in range(1, k_C + 1)), f"milp-2LMH-(81)-4-{m}_SG{self.name}"
			MILP += pulp.lpSum(ec_TC[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
					pulp.lpSum(delta_beta_minus[(k, m)] for k in range(1, k_C + 1)), f"milp-2LMH-(81)-5-{m}_SG{self.name}"
			MILP += pulp.lpSum(ec_CF[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
					pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)), f"milp-2LMH-(81)-6-{m}_SG{self.name}"
			MILP += pulp.lpSum(ec_TF[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
					pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-2LMH-(81)-8-{m}_SG{self.name}"

		# -------- Constraint (82) --------
		for i in range(k_C_tilde + 1, m_C + 1):
			MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_C[(i, ((a, d1), (b, d2), m))] 
						for ((a, d1), (b, d2), m) in Gamma_tilde_ec_C) == \
					pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_C[(i, nu)] 
						for nu in Gamma_tilde_ac_C), f"milp-2LMH-(82)-1-{i}_SG{self.name}"
			MILP += Delta_ec_C_plus[i] + pulp.lpSum(d * delta_ec_C[(i, ((a, d), b, m))]
						for ((a, d), b, m) in Gamma_tilde_ec_C) == deg_C[tail_C[i]], f"milp-2LMH-(82)-2-{i}_SG{self.name}"
			MILP += Delta_ec_C_minus[i] + pulp.lpSum(d * delta_ec_C[(i, (a, (b, d), m))]
						for (a, (b, d), m) in Gamma_tilde_ec_C) == deg_C[head_C[i]], f"milp-2LMH-(82)-3-{i}_SG{self.name}"
			MILP += Delta_ec_C_plus[i] + Delta_ec_C_minus[i] <= 8 * (1 - e_C[i]), f"milp-2LMH-(82)-4-{i}_SG{self.name}"
		for gamma in Gamma_tilde_ec_C:
			MILP += pulp.lpSum(delta_ec_C[(i, gamma)] 
					for i in range(k_C_tilde + 1, m_C + 1)) == ec_C[gamma], f"milp-2LMH-(82)-5-{gamma}_SG{self.name}"
		
		# -------- Constraint (83) --------
		for nu in Gamma_lnk:
			 MILP += pulp.lpSum(delta_ec_C[(i, nu)]
				  for i in (set(I_lnk) & set(range(k_C + 1, m_C + 1)))) == ec_C_lnk[nu], f"milp-2LMM-(83)-{nu}_SG{self.name}"

		# -------- Constraint (84) --------
		for i in range(2, t_T + 1):
			MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_T[(i, ((a, d1), (b, d2), m))] 
						for ((a, d1), (b, d2), m) in Gamma_tilde_ec_T) == \
					pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_T[(i, nu)] 
						for nu in Gamma_tilde_ac_T), f"milp-2LMH-(84)-1-{i}_SG{self.name}"
			MILP += Delta_ec_T_plus[i] + pulp.lpSum(d * delta_ec_T[(i, ((a, d), b, m))]
						for ((a, d), b, m) in Gamma_tilde_ec_T) == deg_T[i - 1], f"milp-2LMH-(84)-2-{i}_SG{self.name}"
			MILP += Delta_ec_T_minus[i] + pulp.lpSum(d * delta_ec_T[(i, (a, (b, d), m))]
						for (a, (b, d), m) in Gamma_tilde_ec_T) == deg_T[i], f"milp-2LMH-(84)-3-{i}_SG{self.name}"
			MILP += Delta_ec_T_plus[i] + Delta_ec_T_minus[i] <= 8 * (1 - e_T[i]), f"milp-2LMH-(84)-4-{i}_SG{self.name}"
		for gamma in Gamma_tilde_ec_T:
			MILP += pulp.lpSum(delta_ec_T[(i, gamma)] 
					for i in range(2, t_T + 1)) == ec_T[gamma], f"milp-2LMH-(84)-5-{gamma}_SG{self.name}"

		# -------- Constraint (85) --------  IDO
		for i in range(2, t_T + 1):
			for nu in Gamma_lnk:
				MILP += delta_ec_T[(i, nu)] + pulp.lpSum(chi_T[(i, k)] 
						for k in (set(I_lnk) & set(range(1, k_C + 1)))) >= 2 * delta_ec_T_lnk[(i, nu)], f"milp-2LMH-(85)-1-{i}-{nu}-1_SG{self.name}"
				MILP += delta_ec_T_lnk[(i, nu)] >= delta_ec_T[(i, nu)] + pulp.lpSum(chi_T[(i, k)] 
							for k in (set(I_lnk) & set(range(1, k_C + 1)))) - 1, f"milp-2LMH-(85)-1-{i}-{nu}-2_SG{self.name}"
		for nu in Gamma_lnk:
			MILP += pulp.lpSum(delta_ec_T_lnk[(i, nu)] for i in range(2, t_T + 1)) == ec_T_lnk[nu], f"milp-2LMH-(85)-2-{nu}_SG{self.name}"  

		# -------- Constraint (86) --------
		for i in range(2, t_F + 1):
			MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_F[(i, ((a, d1), (b, d2), m))]
							   for ((a, d1), (b, d2), m) in Gamma_tilde_ec_F) == \
					pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_F[(i, nu)]
							   for nu in Gamma_tilde_ac_F), f"milp-2LMH-(86)-1-{i}_SG{self.name}"
			MILP += Delta_ec_F_plus[i] + pulp.lpSum(d * delta_ec_F[(i, ((a, d), b, m))]
						for ((a, d), b, m) in Gamma_tilde_ec_F) == deg_F[i - 1], f"milp-2LMH-(86)-2-{i}_SG{self.name}"
			MILP += Delta_ec_F_minus[i] + pulp.lpSum(d * delta_ec_F[(i, (a, (b, d), m))]
						for (a, (b, d), m) in Gamma_tilde_ec_F) == deg_F[i], f"milp-2LMH-(86)-3-{i}_SG{self.name}"
			MILP += Delta_ec_F_plus[i] + Delta_ec_F_minus[i] <= 8 * (1 - e_F[i]), f"milp-2LMH-(86)-4-{i}_SG{self.name}"
		for gamma in Gamma_tilde_ec_F:
			MILP += pulp.lpSum(delta_ec_F[(i, gamma)]
							   for i in range(2, t_F + 1)) == ec_F[gamma], f"milp-2LMH-(86)-5-{gamma}_SG{self.name}"

		# -------- Constraint (87) --------
		for k in range(1, k_C + 1):
			for i in range(1, t_T + 1):
				MILP += deg_T[i] + self.MAX_VAL * (1 - chi_T[(i, k)] + e_T[i]) >= \
						deg_T_CT[k], f"milp-2LMH-(87)-1-{k}-{i}_SG{self.name}"
				MILP += deg_T_CT[k] >= deg_T[i] - \
						self.MAX_VAL * (1 - chi_T[(i, k)] + e_T[i]), f"milp-2LMH-(87)-2-{k}-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_CT[(k, ((a, d1), (b, d2), m))] 
						for ((a, d1), (b, d2), m) in Gamma_tilde_ec_CT) == \
					pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_CT[(k, nu)] 
						for nu in Gamma_tilde_ac_CT), f"milp-2LMH-(87)-3-{k}_SG{self.name}"
			MILP += Delta_ec_CT_plus[k] + pulp.lpSum(d * delta_ec_CT[(k, ((a, d), b, m))]
						for ((a, d), b, m) in Gamma_tilde_ec_CT) == deg_C[tail_C[k]], f"milp-2LMH-(87)-4-{k}_SG{self.name}"
			MILP += Delta_ec_CT_minus[k] + pulp.lpSum(d * delta_ec_CT[(k, (a, (b, d), m))]
						for (a, (b, d), m) in Gamma_tilde_ec_CT) == deg_T_CT[k], f"milp-2LMH-(87)-5-{k}_SG{self.name}"
			MILP += Delta_ec_CT_plus[k] + Delta_ec_CT_minus[k] <= 8 * (1 - delta_chi_T[k]), f"milp-2LMH-(87)-6-{k}_SG{self.name}"
		for gamma in Gamma_tilde_ec_CT:
			MILP += pulp.lpSum(delta_ec_CT[(k, gamma)] for k in range(1, k_C + 1)) == ec_CT[gamma], f"milp-2LMH-(87)-7-{gamma}_SG{self.name}"

		# -------- Constraint (88) --------   IDO
		for nu in Gamma_lnk:
			MILP += pulp.lpSum(delta_ec_CT[(i, nu)] for i in (set(I_lnk) & set(range(1, k_C + 1)))) == ec_CT_lnk[nu], f"milp-2LMM-(88)-{nu}_SG{self.name}"

		# -------- Constraint (89) --------
		for k in range(1, k_C + 1):
			for i in range(1, t_T + 1):
				MILP += deg_T[i] + self.MAX_VAL * (1 - chi_T[(i, k)] + e_T[i + 1]) >= \
						deg_T_TC[k], f"milp-2LMH-(89)-1-{k}-{i}_SG{self.name}"
				MILP += deg_T_TC[k] >= deg_T[i] - \
						self.MAX_VAL * (1 - chi_T[(i, k)] + e_T[i + 1]), f"milp-2LMH-(89)-2-{k}-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_TC[(k, ((a, d1), (b, d2), m))]
						for ((a, d1), (b, d2), m) in Gamma_tilde_ec_TC) == \
					pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_TC[(k, nu)]
						for nu in Gamma_tilde_ac_TC), f"milp-2LMH-(89)-3-{k}_SG{self.name}"
			MILP += Delta_ec_TC_plus[k] + pulp.lpSum(d * delta_ec_TC[(k, ((a, d), b, m))]
						for ((a, d), b, m) in Gamma_tilde_ec_TC) == deg_T_TC[k], f"milp-2LMH-(89)-4-{k}_SG{self.name}"
			MILP += Delta_ec_TC_minus[k] + pulp.lpSum(d * delta_ec_TC[(k, (a, (b, d), m))]
						for (a, (b, d), m) in Gamma_tilde_ec_TC) == deg_C[head_C[k]], f"milp-2LMH-(89)-5-{k}_SG{self.name}"
			MILP += Delta_ec_TC_plus[k] + Delta_ec_TC_minus[k] <= 8 * (1 - delta_chi_T[k]), f"milp-2LMH-(89)-6-{k}_SG{self.name}"
		for gamma in Gamma_tilde_ec_TC:
			MILP += pulp.lpSum(delta_ec_TC[(k, gamma)] for k in range(1, k_C + 1)) == ec_TC[gamma], f"milp-2LMH-(89)-7-{gamma}_SG{self.name}"

		# -------- Constraint (90) --------  IDO
		for nu in Gamma_lnk:
			MILP += pulp.lpSum(delta_ec_TC[(i, nu)] for i in (set(I_lnk) & set(range(1, k_C + 1)))) == ec_TC_lnk[nu], f"milp-2LMM-(90)-{nu}_SG{self.name}"

		# -------- Constraint (91) --------
		for c in range(1, t_C_tilde + 1):
			for i in range(1, t_F + 1):
				MILP += deg_F[i] + self.MAX_VAL * (1 - chi_F[(i, c)] + e_F[i]) >= \
						deg_F_CF[c], f"milp-2LMH-(91)-1-{c}-{i}_SG{self.name}"
				MILP += deg_F_CF[c] >= deg_F[i] - \
						self.MAX_VAL * (1 - chi_F[(i, c)] + e_F[i]), f"milp-2LMH-(91)-2-{c}-{i}_SG{self.name}"
			MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_CF[(c, ((a, d1), (b, d2), m))] 
						for ((a, d1), (b, d2), m) in Gamma_tilde_ec_CF) == \
					pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_CF[(c, nu)] 
						for nu in Gamma_tilde_ac_CF), f"milp-2LMH-(91)-3-{c}_SG{self.name}"
			MILP += Delta_ec_CF_plus[c] + pulp.lpSum(d * delta_ec_CF[(c, ((a, d), b, m))]
						for ((a, d), b, m) in Gamma_tilde_ec_CF) == deg_C[c], f"milp-2LMH-(91)-4-{c}_SG{self.name}"
			MILP += Delta_ec_CF_minus[c] + pulp.lpSum(d * delta_ec_CF[(c, (a, (b, d), m))]
						for (a, (b, d), m) in Gamma_tilde_ec_CF) == deg_F_CF[c], f"milp-2LMH-(91)-5-{c}_SG{self.name}"
			MILP += Delta_ec_CF_plus[c] + Delta_ec_CF_minus[c] <= 8 * (1 - delta_chi_F[c]), f"milp-2LMH-(91)-6-{c}_SG{self.name}"
		for gamma in Gamma_tilde_ec_CF:
			MILP += pulp.lpSum(delta_ec_CF[(c, gamma)] for c in range(1, t_C_tilde + 1)) == ec_CF[gamma], f"milp-2LMH-(91)-7-{gamma}_SG{self.name}"

		# -------- Constraint (92) --------
		for i in range(1, t_T + 1):
			for j in range(1, t_F + 1):
				MILP += deg_F[j] + self.MAX_VAL * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]) >= \
						deg_F_TF[i], f"milp-2LMH-(92)-1-{i}-{j}_SG{self.name}"
				MILP += deg_F_TF[i] >= deg_F[j] - \
						self.MAX_VAL * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]), f"milp-2LMH-(92)-2-{i}-{j}_SG{self.name}"
			MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_TF[(i, ((a, d1), (b, d2), m))] 
						for ((a, d1), (b, d2), m) in Gamma_tilde_ec_TF) == \
					pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_TF[(i, nu)] 
						for nu in Gamma_tilde_ac_TF), f"milp-2LMH-(92)-3-{i}_SG{self.name}"
			MILP += Delta_ec_TF_plus[i] + pulp.lpSum(d * delta_ec_TF[(i, ((a, d), b, m))]
						for ((a, d), b, m) in Gamma_tilde_ec_TF) == deg_T[i], f"milp-2LMH-(92)-4-{i}_SG{self.name}"
			MILP += Delta_ec_TF_minus[i] + pulp.lpSum(d * delta_ec_TF[(i, (a, (b, d), m))]
						for (a, (b, d), m) in Gamma_tilde_ec_TF) == deg_F_TF[i], f"milp-2LMH-(92)-5-{i}_SG{self.name}"
			MILP += Delta_ec_TF_plus[i] + Delta_ec_TF_minus[i] <= 8 * (1 - delta_chi_F[i + t_C_tilde]), f"milp-2LMH-(92)-6-{i}_SG{self.name}"
		for gamma in Gamma_tilde_ec_TF:
			MILP += pulp.lpSum(delta_ec_TF[(i, gamma)] for i in range(1, t_T + 1)) == ec_TF[gamma], f"milp-2LMH-(92)-7-{gamma}_SG{self.name}"
		
		# -------- Constraint (93) --------
		for (mu1, mu2, m) in Gamma_int_less:
			MILP += ec_C[(mu1, mu2, m)] + ec_C[(mu2, mu1, m)] + \
					ec_T[(mu1, mu2, m)] + ec_T[(mu2, mu1, m)] + \
					ec_F[(mu1, mu2, m)] + ec_F[(mu2, mu1, m)] + \
					ec_CT[(mu1, mu2, m)] + ec_CT[(mu2, mu1, m)] + \
					ec_TC[(mu1, mu2, m)] + ec_TC[(mu2, mu1, m)] + \
					ec_CF[(mu1, mu2, m)] + ec_CF[(mu2, mu1, m)] + \
					ec_TF[(mu1, mu2, m)] + ec_TF[(mu2, mu1, m)] == ec_int[(mu1, mu2, m)], f"milp-2LMH-(93)-1-{(mu1, mu2, m)}_SG{self.name}"
		for gamma in Gamma_int_equal:
			MILP += ec_C[gamma] + ec_T[gamma] + ec_F[gamma] + \
					ec_CT[gamma] + ec_TC[gamma] + \
					ec_CF[gamma] + ec_TF[gamma] == ec_int[gamma], f"milp-2LMH-(93)-2-{gamma}_SG{self.name}"

		# -------- Constraint (94) --------  IDO
		for (a, b, m) in Gamma_lnk_less:
			MILP += ec_C_lnk[(a, b, m)] + ec_C_lnk[(b, a, m)] + \
					ec_T_lnk[(a, b, m)] + ec_T_lnk[(b, a, m)] + \
					ec_CT_lnk[(a, b, m)] + ec_CT_lnk[(b, a, m)] + \
					ec_TC_lnk[(a, b, m)] + ec_TC_lnk[(b, a, m)]  == ec_lnk[(a, b, m)], f"milp-2LMM-(94)-1-{(a, b, m)}_SG{self.name}"
		for nu in Gamma_lnk_equal:
			MILP += ec_C_lnk[nu] + ec_T_lnk[nu] + \
					ec_CT_lnk[nu] + ec_TC_lnk[nu]  == ec_lnk[nu], f"milp-2LMM-(94)-2-{nu}_SG{self.name}"
		MILP += pulp.lpSum(ec_lnk[gamma] for gamma in Gamma_lnk_less) + \
				pulp.lpSum(ec_lnk[gamma] for gamma in Gamma_lnk_equal) == n_lnk, f"milp-2LMM-polymer-lnk_ec_SG{self.name}"

		# -------- Constraint (101) --------
		MILP += pulp.lpSum(delta_cnt[gamma] for gamma in (Gamma_cnt_less + Gamma_cnt_equal + Gamma_cnt_great)) == 1, f"milp-2LMM-polymer-(101)-1_SG{self.name}"

		# -------- Constraint (102) --------
		for (mu1, mu2, m) in (Gamma_cnt_less + Gamma_cnt_equal):
			MILP += delta_cnt[(mu1, mu2, m)] <= ec_lnk[(mu1, mu2, m)], f"milp-2LMM-polymer-(101)-2-{(mu1, mu2, m)}_SG{self.name}"
		for (mu1, mu2, m) in Gamma_cnt_great:
			MILP += delta_cnt[(mu1, mu2, m)] <= ec_lnk[(mu2, mu1, m)], f"milp-2LMM-polymer-(101)-3-{(mu1, mu2, m)}_SG{self.name}"

		return MILP

	def prepare_variables_fringe_configuration(self,
		t_C, t_T, t_F, set_F, Code_F, fc_LB, fc_UB
	):
		#Define fc
		fc = {Code_F[psi]: pulp.LpVariable(f"fc({Code_F[psi]})_SG{self.name}", fc_LB[Code_F[psi]], fc_UB[Code_F[psi]], cat=pulp.LpInteger)
			for psi in set_F}

		return fc

	def add_constraints_fringe_configuration(self,
		# Model
		MILP,
		# Constants
		t_C, t_T, t_F, set_F, Code_F,
		# Binary Variables
		delta_fr_C, delta_fr_T, delta_fr_F,
		# Integer Variables
		fc
	):
		# -------- Constraint (85) --------
		for psi in set_F:
			MILP += pulp.lpSum(delta_fr_C[(i, Code_F[psi])] for i in range(1, t_C + 1)) + \
					pulp.lpSum(delta_fr_T[(i, Code_F[psi])] for i in range(1, t_T + 1)) + \
					pulp.lpSum(delta_fr_F[(i, Code_F[psi])] for i in range(1, t_F + 1)) == fc[Code_F[psi]], f"milp-2LMH-(85)-{Code_F[psi]}_SG{self.name}"

		return MILP

	def add_constraints_mass_n(self,
		MILP,
		n_LB, n_star, n_G, na_UB, na_ex, MASS
	):
		delta_n = {n: pulp.LpVariable(f"delta_n({n})_SG{self.name}".format(n), 0, 1, cat=pulp.LpBinary) for n in range(n_LB, n_star + na_UB["H1"] + 1)}

		MILP += pulp.lpSum(delta_n[n] for n in range(n_LB, n_star + na_UB["H1"] + 1)) == 1, f"milp-n-1_SG{self.name}"
		MILP += pulp.lpSum(n * delta_n[n] for n in range(n_LB, n_star + na_UB["H1"] + 1)) == n_G + na_ex["H1"], f"milp-n-2_SG{self.name}"

		mass_n = pulp.LpVariable(f"mass_n_SG{self.name}", 0, 100000)
		for n in range(n_LB, n_star + na_UB["H1"] + 1):
			MILP += mass_n <= MASS * (1 / n) + 100000 * (1 - delta_n[n]), f"milp-mass_n-{n}-1_SG{self.name}"
			MILP += mass_n >= MASS * (1 / n) - 100000 * (1 - delta_n[n]), f"milp-mass_n-{n}-2_SG{self.name}"

		return MILP, delta_n, mass_n

	def add_constraints_prevent_continuous_double_bond_complete(self,
		MILP,
		t_C, t_C_tilde, t_T, t_F, beta_r_2, Code_F, F_C, F_T, F_F, E_C,
		I_equal_one, I_zero_one, I_ge_one, I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus,
		delta_alpha_C, delta_alpha_T, delta_alpha_F, delta_fr_C, delta_fr_T, delta_fr_F,
		delta_beta_C, delta_beta_T, delta_beta_F, delta_beta_plus, delta_beta_minus, delta_beta_in,
		beta_CT, beta_TC, beta_CF, beta_TF
	):
		delta_beta_CT = {(i, m): pulp.LpVariable(f"delta_beta_CT({i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for i in range(1, t_T + 1) for m in range(0, self.MAX_BOND + 1)}
		delta_beta_TC = {(i, m): pulp.LpVariable(f"delta_beta_TC({i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for i in range(1, t_T + 1) for m in range(0, self.MAX_BOND + 1)}
		delta_beta_CF = {(i, m): pulp.LpVariable(f"delta_beta_CF({i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for i in range(1, t_F + 1) for m in range(0, self.MAX_BOND + 1)}
		delta_beta_TF = {(i, m): pulp.LpVariable(f"delta_beta_TF({i},{m})_SG{self.name}", 0, 1, cat=pulp.LpBinary)
					for i in range(1, t_F + 1) for m in range(0, self.MAX_BOND + 1)}
		
		for i in range(1, t_T + 1):
			MILP += pulp.lpSum(delta_beta_CT[(i, m)] for m in range(0, self.MAX_BOND + 1)) == 1, f"milp-2LMM-prevent-continuous-double-bound_CT-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_CT[(i, m)] * m for m in range(0, self.MAX_BOND + 1)) == beta_CT[i], f"milp-2LMM-prevent-continuous-double-bound_CT-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_TC[(i, m)] for m in range(0, self.MAX_BOND + 1)) == 1, f"milp-2LMM-prevent-continuous-double-bound_TC-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_TC[(i, m)] * m for m in range(0, self.MAX_BOND + 1)) == beta_TC[i], f"milp-2LMM-prevent-continuous-double-bound_TC-2-{i}_SG{self.name}"
		for i in range(1, t_F + 1):
			MILP += pulp.lpSum(delta_beta_CF[(i, m)] for m in range(0, self.MAX_BOND + 1)) == 1, f"milp-2LMM-prevent-continuous-double-bound_CF-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_CF[(i, m)] * m for m in range(0, self.MAX_BOND + 1)) == beta_CF[i], f"milp-2LMM-prevent-continuous-double-bound_CF-2-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_TF[(i, m)] for m in range(0, self.MAX_BOND + 1)) == 1, f"milp-2LMM-prevent-continuous-double-bound_TF-1-{i}_SG{self.name}"
			MILP += pulp.lpSum(delta_beta_TF[(i, m)] * m for m in range(0, self.MAX_BOND + 1)) == beta_TF[i], f"milp-2LMM-prevent-continuous-double-bound_TF-2-{i}_SG{self.name}"

		for i in range(1, t_C_tilde + 1):
			MILP += pulp.lpSum(delta_beta_C[(j, 2)] for j in range(1, len(E_C) + 1) 
						if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
					pulp.lpSum(delta_beta_plus[(k, 2)] for k in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
					pulp.lpSum(delta_beta_minus[(k, 2)] for k in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
					delta_beta_in[(i, 2)] + pulp.lpSum(beta_r_2[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] for psi in F_C[i]) <= \
					4 - 3 * delta_alpha_C[(i, "C4")], f"milp-2LMM-prevent-continuous-double-bound_C-{i}_SG{self.name}"

		for i in range(t_C_tilde + 1, t_C + 1):
			MILP += pulp.lpSum(delta_beta_C[(j, 2)] for j in range(1, len(E_C) + 1) 
						if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
					pulp.lpSum(delta_beta_plus[(c, 2)] for c in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
					pulp.lpSum(delta_beta_minus[(c, 2)] for c in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
					pulp.lpSum(beta_r_2[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] for psi in F_C[i]) <= \
					4 - 3 * delta_alpha_C[(i, "C4")], f"milp-2LMM-prevent-continuous-double-bound_C-{i}_SG{self.name}"

		for i in range(1, t_T + 1):
			MILP += (delta_beta_T[(i, 2)] if i != 1 else 0) + (delta_beta_T[(i + 1, 2)] if i != t_T else 0) + pulp.lpSum(beta_r_2[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] for psi in F_T[i]) + \
					delta_beta_CT[(i, 2)] + delta_beta_TC[(i, 2)] + delta_beta_in[(t_C_tilde + i, 2)] <= \
					4 - 3 * delta_alpha_T[(i, "C4")], f"milp-2LMM-prevent-continuous-double-bound_T-{i}_SG{self.name}"

		for i in range(1, t_F + 1):
			MILP += (delta_beta_F[(i, 2)] if i != 1 else 0) + (delta_beta_F[(i + 1, 2)] if i != t_F else 0)+ pulp.lpSum(beta_r_2[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] for psi in F_F[i]) + \
					delta_beta_CF[(i, 2)] + delta_beta_TF[(i, 2)] <= \
					4 - 3 * delta_alpha_F[(i, "C4")], f"milp-2LMM-prevent-continuous-double-bound_F-{i}_SG{self.name}"

		return MILP

	def add_constraints_prevent_continuous_double_bond(self,
		MILP,
		t_C, t_T, t_F,
		delta_alpha_C, delta_alpha_T, delta_alpha_F,
		deg_C, deg_T, deg_F, hyddeg_C, hyddeg_T, hyddeg_F
	):
		for i in range(1, t_C + 1):
			MILP += deg_C[i] + hyddeg_C[i] >= 3 * delta_alpha_C[(i, "C4")], f"milp-2LMM-prevent-continuous-double-bound_C-{i}_SG{self.name}"
		for i in range(1, t_T + 1):
			MILP += deg_T[i] + hyddeg_T[i] >= 3 * delta_alpha_T[(i, "C4")], f"milp-2LMM-prevent-continuous-double-bound_T-{i}_SG{self.name}"
		for i in range(1, t_F + 1):
			MILP += deg_F[i] + hyddeg_F[i] >= 3 * delta_alpha_F[(i, "C4")], f"milp-2LMM-prevent-continuous-double-bound_F-{i}_SG{self.name}"
		
		return MILP

	def set_up_seed_graph(self):
		set_Lambda = prepare_CG_element_info()

		self.set_Lambda = self.prepare_Lambda(set_Lambda, self.Lambda)
		self.set_Lambda_dg, self.Lambda_dg, self.epsilon, self.epsilon_dg, self.Code_Lambda_dg, \
			self.Code_Lambda, self.Code_Lambda_int, self.Code_Lambda_ex, \
			self.MAX_CODE, self.MAX_CODE_dg, self.MAX_CODE_int, self.MAX_CODE_ex, \
			self.Code_Gamma_ec_int, self.Code_Gamma_ac_int, self.val, self.mass, \
			self.ele_neg, self.ion_energy, self.ele_affinity = \
				self.prepare_Lambda_dg( 
					self.set_Lambda, self.Lambda, self.Lambda_int, self.Lambda_dg_int, self.Lambda_ex, self.Gamma_int, self.Gamma_int_ac
				)
		self.Gamma_int_less, self.Gamma_int_equal, self.Gamma_int_great, \
			self.Gamma_int_ac_less, self.Gamma_int_ac_equal, self.Gamma_int_ac_great, \
			self.Gamma_tilde_ac_C, self.Gamma_tilde_ac_T, self.Gamma_tilde_ac_CT, self.Gamma_tilde_ac_TC, \
			self.Gamma_tilde_ac_F, self.Gamma_tilde_ac_CF, self.Gamma_tilde_ac_TF, \
			self.Gamma_tilde_ec_C, self.Gamma_tilde_ec_T, self.Gamma_tilde_ec_CT, self.Gamma_tilde_ec_TC, \
			self.Gamma_tilde_ec_F, self.Gamma_tilde_ec_CF, self.Gamma_tilde_ec_TF , \
			self.Gamma_lnk_less, self.Gamma_lnk_equal, self.Gamma_lnk_great, \
			self.Gamma_ac_lnk_less, self.Gamma_ac_lnk_equal, self.Gamma_ac_lnk_great, \
			self.Gamma_cnt_less, self.Gamma_cnt_equal, self.Gamma_cnt_great = \
				self.prepare_Gamma_ac(
					self.Gamma_int, self.Gamma_int_ac, self.Code_Lambda, self.Code_Lambda_dg,
					self.Gamma_lnk, self.Gamma_ac_lnk, self.ns_LB_cnt, self.ns_UB_cnt
				)
		self.Gamma_lf_ac, self.ac_LB_lf, self.ac_UB_lf = \
			self.prepare_Gamma_lf_ac(
				self.set_Lambda, self.ac_LB_lf, self.ac_UB_lf, self.ac_LB_lf_common, self.ac_UB_lf_common
			)
		self.Lambda_tilde_dg_C, self.Lambda_tilde_dg_T, self.Lambda_tilde_dg_F = \
			self.prepare_Lambda_tilde(self.Lambda_dg_int)
		self.t_C, self.t_C_tilde, self.t_T, self.t_F, self.k_C, self.k_C_tilde, self.c_F, self.m_C, \
			self.head_C, self.tail_C, self.tail_F, self.E_C_plus, self.E_C_minus, \
			self.I_lnk, self.I_ge_one_plus, self.I_ge_one_minus, self.I_ge_two_plus, self.I_ge_two_minus, \
			self.I_zero_one_plus, self.I_zero_one_minus, self.I_equal_one_plus, self.I_equal_one_minus, \
			self.delta_i, self.n_lnk_equal_one = \
				self.prepare_dataset_for_scheme_graph(
					self.n_star, self.rho, self.V_C, self.E_C, self.E_C_lnk, self.E_ge_two, self.E_ge_one,
					self.n_LB_int, self.n_UB_int, self.bl_LB, self.bl_UB, self.ch_LB, self.ch_UB, 
					self.ell_UB, self.I_ge_two, self.I_ge_one, self.I_equal_one, self.I_zero_one, 
					self.bd2_LB, self.bd3_LB, self.val, self.Lambda_star
				)
		self.Code_F, self.n_psi_H, self.deg_r_H, self.deg_r_hyd, self.beta_r, self.atom_r, self.ht_H, \
			self.F_C, self.F_T, self.F_F, self.n_C, self.n_T, self.n_F, self.F_Cp, self.F_Tp, self.F_Fp, \
			self.set_F_E, self.set_F_v, self.na_alpha_ex, self.alpha_r, self.deg_fr, self.v_ion, self.ac_psi_lf, self.beta_r_2 = \
				self.prepare_fringe_tree(
					self.set_F, self.V_C, self.t_T, self.t_F, self.val, self.Lambda_int, 
					self.Lambda_ex, self.Code_Lambda_int, self.Gamma_lf_ac
				)

	def set_up_variables(self):
		self.e_C, self.v_T, self.e_T, self.chi_T, self.clr_T, self.delta_chi_T, self.deg_tilde_C_plus, \
			self.deg_tilde_C_minus, self.n_lnk, self.rank_G = \
				self.prepare_variables_selecting_core(
					self.t_C, self.k_C_tilde, self.k_C, self.t_T, self.m_C, self.n_LB_int, self.n_UB_int, 
					self.n_LB_lnk, self.n_UB_lnk, self.ell_LB, self.ell_UB, self.bl_LB, self.bl_UB
				)		
		self.n_G_int, self.v_F, self.e_F, self.chi_F, self.clr_F, self.delta_chi_F, self.bl = \
			self.prepare_variables_internal_vertices_and_edges(
				self.t_F, self.t_T, self.t_C, self.t_C_tilde, self.c_F, self.n_LB_int, self.n_UB_int, self.V_C, 
				self.E_ge_one, self.E_ge_two, self.I_ge_one, self.I_ge_two, self.bl_LB, self.bl_UB, self.E_C
			)
		self.n_G, self.h_T, self.h_C, self.h_F, self.delta_fr_C, self.delta_fr_T, self.delta_fr_F, \
			self.deg_ex_C, self.deg_ex_T, self.deg_ex_F, self.hyddeg_C, self.hyddeg_T, self.hyddeg_F, \
			self.eledeg_C, self.eledeg_T, self.eledeg_F, self.sigma, self.ac_lf = \
				self.prepare_variables_fringe_trees(
					self.n_LB, self.n_star, self.rho, self.ch_LB, self.ch_UB, self.t_T, self.t_C, self.t_F,
					self.n_T, self.n_C, self.n_F, self.delta_i, self.I_ge_two, self.I_ge_one, self.V_C, 
					self.E_ge_one, self.E_ge_two, self.v_T, self.v_F, self.F_C, self.F_T, self.F_F, 
					self.Code_F, self.Gamma_lf_ac, self.ac_LB_lf, self.ac_UB_lf
				)
		self.deg_C, self.deg_T, self.deg_F, self.deg_CT, self.deg_TC, \
			self.delta_dg_C, self.delta_dg_T, self.delta_dg_F, self.dg, \
			self.deg_int_C, self.deg_int_T, self.deg_int_F, \
			self.delta_int_dg_C, self.delta_int_dg_T, self.delta_int_dg_F, self.dg_int = \
				self.prepare_variables_degree(
					self.t_C, self.t_T, self.t_F, self.n_C, self.n_T, self.n_F,
					self.delta_i, self.n_star, self.dg_LB, self.dg_UB, self.rho
				)
		self.beta_T, self.beta_F, self.beta_C, self.beta_plus, self.beta_minus, self.beta_in, \
			self.beta_ex_C, self.beta_ex_T, self.beta_ex_F, self.delta_beta_T, self.delta_beta_F, self.delta_beta_C, \
			self.delta_beta_plus, self.delta_beta_minus, self.delta_beta_in, self.bd_int, \
			self.bd_C, self.bd_T, self.bd_CT, self.bd_TC, self.bd_F, self.bd_CF, self.bd_TF = \
				self.prepare_variables_multiplicity(
					self.t_C, self.t_T, self.t_F, self.n_C, self.n_T, self.n_F, self.c_F, self.delta_i,
					self.I_equal_one, self.I_zero_one, self.I_ge_one, self.I_ge_two, self.E_C,
					self.bd2_LB, self.bd2_UB, self.bd3_LB, self.bd3_UB, self.n_UB_int, self.n_star, self.rho
				)
		self.beta_CT, self.beta_TC, self.beta_CF, self.beta_TF, self.alpha_C, self.alpha_T, self.alpha_F, self.MASS, \
			self.delta_alpha_C, self.delta_alpha_T, self.delta_alpha_F, self.na, self.na_int, self.na_C, self.na_T, self.na_F, \
			self.na_ex_C, self.na_ex_T, self.na_ex_F, self.na_ex = \
				self.prepare_variables_chemical_elements(
					self.t_C, self.t_T, self.t_F, self.n_C, self.n_T, self.n_F, self.delta_i,
					self.n_star, self.Lambda, self.epsilon, self.na_LB, self.na_UB, self.rho, 
					self.na_LB_int, self.na_UB_int, self.Lambda_int, self.Lambda_ex, self.MAX_CODE,
					self.MAX_CODE_int, self.MAX_CODE_ex, self.Code_Lambda_int, self.Code_Lambda_ex, 
					self.F_C, self.F_T, self.F_F, self.Code_F, self.alpha_r, self.na_alpha_ex
				)
		self.bd_T = self.prepare_variables_number_of_bounds(self.k_C, self.t_T, self.bd_T)
		self.ac_int, self.ac_C, self.ac_T, self.ac_F, self.ac_CT, self.ac_TC, self.ac_CF, self.ac_TF, \
			self.delta_ac_C, self.delta_ac_T, self.delta_ac_F, \
			self.alpha_CT, self.alpha_TC, self.alpha_CF, self.alpha_TF, \
			self.Delta_ac_C_plus, self.Delta_ac_C_minus, self.Delta_ac_T_plus, self.Delta_ac_T_minus, \
			self.Delta_ac_F_plus, self.Delta_ac_F_minus, \
			self.Delta_ac_CT_plus, self.Delta_ac_CT_minus, self.Delta_ac_TC_plus, self.Delta_ac_TC_minus, \
			self.Delta_ac_CF_plus, self.Delta_ac_CF_minus, self.Delta_ac_TF_plus, self.Delta_ac_TF_minus, \
			self.delta_ac_CT, self.delta_ac_TC, self.delta_ac_CF, self.delta_ac_TF, \
			self.ac_lnk, self.ac_C_lnk, self.ac_T_lnk, self.ac_CT_lnk, self.ac_TC_lnk, self.delta_ac_T_lnk = \
				self.prepare_variables_adjacency_configuration(
					self.t_C, self.t_C_tilde, self.t_T, self.t_F, self.m_C, self.k_C, self.k_C_tilde, self.c_F, 
					self.n_T, self.n_C, self.n_F, self.delta_i, self.rho, self.Gamma_int_ac, self.Gamma_tilde_ac_C, 
					self.Gamma_tilde_ac_T, self.Gamma_tilde_ac_F, self.Gamma_tilde_ac_CT, self.Gamma_tilde_ac_TC, 
					self.Gamma_tilde_ac_CF, self.Gamma_tilde_ac_TF, 
					self.ac_LB_int, self.ac_UB_int, self.Lambda_int, self.Lambda_ex, self.MAX_CODE_int, self.MAX_CODE_ex,
					self.Gamma_ac_lnk, self.ac_LB_lnk, self.ac_UB_lnk
				)
		self.ns_int, self.delta_ns_C, self.delta_ns_T, self.delta_ns_F = \
			self.prepare_variables_chemical_symbols(
				self.t_C, self.t_T, self.t_F, self.n_C, self.n_T, self.n_F, self.delta_i,
				self.n_star, self.n_LB_int, self.n_UB_int, self.ns_LB_int, self.ns_UB_int, self.Lambda_dg_int, self.epsilon_dg
			)
		self.ec_int, self.ec_C, self.ec_T, self.ec_F, self.ec_CT, self.ec_TC, self.ec_CF, self.ec_TF, \
			self.delta_ec_C, self.delta_ec_T, self.delta_ec_F, \
			self.delta_ec_CT, self.delta_ec_TC, self.delta_ec_CF, self.delta_ec_TF, \
			self.deg_T_CT, self.deg_T_TC, self.deg_F_CF, self.deg_F_TF, \
			self.Delta_ec_C_plus, self.Delta_ec_C_minus, self.Delta_ec_T_plus, self.Delta_ec_T_minus, \
			self.Delta_ec_F_plus, self.Delta_ec_F_minus, \
			self.Delta_ec_CT_plus, self.Delta_ec_CT_minus, self.Delta_ec_TC_plus, self.Delta_ec_TC_minus, \
			self.Delta_ec_CF_plus, self.Delta_ec_CF_minus, self.Delta_ec_TF_plus, self.Delta_ec_TF_minus, \
			self.ec_lnk, self.ec_C_lnk, self.ec_T_lnk, self.ec_CT_lnk, self.ec_TC_lnk, self.delta_ec_T_lnk, self.delta_cnt = \
				self.prepare_variables_edge_configuration(
					self.t_C, self.t_C_tilde, self.t_T, self.t_F, self.m_C, self.k_C, self.k_C_tilde, self.c_F,
					self.n_T, self.n_C, self.n_F, self.delta_i, self.rho, self.Gamma_int,
					self.Gamma_tilde_ec_C, self.Gamma_tilde_ec_T, self.Gamma_tilde_ec_F,
					self.Gamma_tilde_ec_CT, self.Gamma_tilde_ec_TC, self.Gamma_tilde_ec_CF, self.Gamma_tilde_ec_TF,
					self.ec_LB_int, self.ec_UB_int, self.ec_LB_lnk, self.ec_UB_lnk, self.Gamma_lnk,
					self.Lambda_dg_int, self.Gamma_cnt_less, self.Gamma_cnt_equal, self.Gamma_cnt_great
				)
		self.fc = self.prepare_variables_fringe_configuration(
			self.t_C, self.t_T, self.t_F, self.set_F, self.Code_F, self.fc_LB, self.fc_UB
		)

	def add_constraints_C2(self, MILP):
		MILP = self.add_constraints_selecting_core(
			MILP,
			self.t_C, self.k_C_tilde, self.k_C, self.t_T, self.m_C, self.n_LB_lnk, self.n_UB_lnk, self.ell_LB, self.ell_UB, 
			self.bl_LB, self.bl_UB, self.I_equal_one, self.I_equal_one_minus, self.I_equal_one_plus,
			self.I_ge_two, self.I_ge_one, self.I_ge_one_minus, self.I_ge_one_plus, self.I_zero_one, 
			self.I_zero_one_minus, self.I_zero_one_plus, self.n_lnk_equal_one, self.I_lnk, self.r_GC,
			self.e_C, self.v_T, self.e_T, self.delta_chi_T, self.chi_T, self.clr_T, 
			self.deg_tilde_C_plus, self.deg_tilde_C_minus, self.n_lnk, self.rank_G
		)
		MILP = self.add_constraints_internal_vertices_and_edges(
			MILP,
			self.t_T, self.t_F, self.t_C, self.t_C_tilde, self.c_F, self.tail_F, 
			self.I_ge_one, self.I_ge_two, self.bl_LB, self.bl_UB, self.E_C, self.delta_chi_F, self.chi_T, 
			self.e_F, self.v_F, self.v_T, self.chi_F, self.clr_F, self.bl, self.n_G_int
		)
		MILP = self.add_constraints_fringe_trees(
			MILP,
			self.t_T, self.t_F, self.t_C, self.t_C_tilde, self.c_F, self.I_ge_one, self.I_ge_two, 
			self.rho, self.n_star, self.n_T, self.n_C, self.n_F, self.ch_LB, self.ch_UB, self.E_C, 
			self.F_C, self.F_T, self.F_F, self.F_Cp, self.F_Tp, self.F_Fp,
			self.Code_F, self.val, self.n_psi_H, self.deg_r_H, self.deg_r_hyd, self.ht_H, 
			self.v_ion, self.ac_psi_lf, self.Gamma_lf_ac,
			self.delta_chi_F, self.delta_chi_T, self.e_F, self.v_F, self.v_T, self.h_T, self.h_C, self.h_F,
			self.sigma, self.delta_fr_C, self.delta_fr_T, self.delta_fr_F,
			self.chi_T, self.chi_F, self.clr_F, self.n_G, self.deg_ex_C, self.deg_ex_T, self.deg_ex_F, 
			self.hyddeg_C, self.hyddeg_T, self.hyddeg_F, self.eledeg_C, self.eledeg_T, self.eledeg_F, self.ac_lf
		)
		MILP = self.add_constraints_degree(
			MILP,
			self.t_T, self.t_F, self.t_C, self.t_C_tilde, self.n_T, self.n_C, self.n_F, self.delta_i,
			self.I_ge_one_plus, self.I_ge_one_minus, self.I_ge_two_plus, self.I_ge_two_minus, self.rho,
			self.F_C, self.F_T, self.F_F, self.F_Cp, self.F_Tp, self.F_Fp, self.Code_F, self.deg_fr,
			self.delta_dg_C, self.delta_dg_T, self.delta_dg_F, self.delta_int_dg_C, self.delta_int_dg_T, self.delta_int_dg_F,
			self.e_T, self.e_F, self.v_F, self.v_T, self.delta_chi_T, self.delta_chi_F, self.delta_fr_C, self.delta_fr_T, self.delta_fr_F,
			self.deg_C, self.deg_T, self.deg_F, self.deg_CT, self.deg_TC, self.dg, self.deg_int_C, self.deg_int_T, self.deg_int_F,
			self.dg_int, self.deg_tilde_C_minus, self.deg_tilde_C_plus, self.deg_ex_C, self.deg_ex_T, self.deg_ex_F, 
			self.hyddeg_C, self.hyddeg_T, self.hyddeg_F
		)
		MILP = self.add_constraints_multiplicity(
			MILP,
			self.t_T, self.t_F, self.t_C, self.t_C_tilde, self.n_T, self.n_C, self.n_F, self.c_F, self.delta_i,
			self.head_C, self.tail_C, self.I_equal_one, self.I_zero_one, self.I_ge_one, self.I_ge_two,
			self.E_C, self.bd2_LB, self.bd2_UB, self.bd3_LB, self.bd3_UB, self.rho, self.F_C, self.F_T, self.F_F,
			self.Code_F, self.beta_r, self.e_C, self.e_T, self.e_F, self.v_F, self.v_T,
			self.delta_chi_T, self.delta_chi_F, self.delta_beta_C, self.delta_beta_T, self.delta_beta_F,
			self.delta_beta_plus, self.delta_beta_minus, self.delta_beta_in,
			self.delta_fr_C, self.delta_fr_T, self.delta_fr_F,
			self.beta_C, self.beta_T, self.beta_F, self.beta_plus, self.beta_minus, self.beta_in,
			self.bd_int, self.bd_C, self.bd_T, self.bd_CT, self.bd_TC, self.bd_F, self.bd_CF, self.bd_TF,
			self.beta_ex_C, self.beta_ex_T, self.beta_ex_F
		)
		MILP = self.add_constraints_chemical_elements(
			MILP,
			self.t_T, self.t_F, self.t_C, self.t_C_tilde, self.n_T, self.n_C, self.n_F, self.c_F, self.delta_i, self.E_C,
			self.I_equal_one, self.I_zero_one, self.I_ge_one, self.I_ge_one_plus, self.I_ge_one_minus,
			self.I_ge_two, self.I_ge_two_plus, self.I_ge_two_minus, self.val, self.mass,
			self.Lambda, self.Code_Lambda, self.Lambda_star, self.epsilon, self.rho, self.na_LB_int, self.na_UB_int,
			self.Lambda_int, self.Lambda_ex, self.MAX_CODE, self.MAX_CODE_int, self.MAX_CODE_ex,
			self.Code_Lambda_int, self.Code_Lambda_ex, self.F_C, self.F_T, self.F_F, self.Code_F, self.alpha_r, self.na_alpha_ex,
			self.v_T, self.v_F, self.e_T, self.e_F, self.delta_chi_T, self.delta_chi_F, self.chi_T, self.chi_F, 
			self.delta_alpha_C, self.delta_alpha_T, self.delta_alpha_F, self.delta_fr_C, self.delta_fr_T, self.delta_fr_F, 
			self.deg_C, self.deg_T, self.deg_F, self.beta_C, self.beta_T, self.beta_F, self.beta_CT, self.beta_TC, self.beta_CF, self.beta_TF,
			self.beta_plus, self.beta_minus, self.beta_in, self.alpha_C, self.alpha_T, self.alpha_F, self.MASS, self.na, self.na_int,
			self.na_C, self.na_T, self.na_F, self.na_ex_C, self.na_ex_T, self.na_ex_F, self.beta_ex_C, self.beta_ex_T, self.beta_ex_F, self.na_ex, self.bd_int,
			self.eledeg_C, self.eledeg_T, self.eledeg_F
		)
		MILP = self.add_constraints_number_of_bounds(
			MILP,
			self.t_T, self.k_C, self.E_C, self.I_equal_one, self.I_zero_one, self.bd2_LB, self.bd2_UB, 
			self.bd3_LB, self.bd3_UB, self.chi_T, self.delta_beta_C, self.delta_beta_T, self.delta_beta_plus, 
			self.delta_beta_minus, self.bd_T
		)
		MILP = self.add_constraints_adjacency_configuration(
			MILP,
			self.t_C, self.t_C_tilde, self.t_T, self.t_F, self.n_C, self.n_T, self.n_F, self.k_C, self.k_C_tilde, 
			self.m_C, self.c_F, self.tail_C, self.head_C,
			self.Gamma_int_ac, self.Gamma_tilde_ac_C, self.Gamma_tilde_ac_T, self.Gamma_tilde_ac_F, self.Gamma_tilde_ac_CT,
			self.Gamma_tilde_ac_TC, self.Gamma_tilde_ac_CF, self.Gamma_tilde_ac_TF, self.Gamma_int_ac_less, self.Gamma_int_ac_equal,
			self.ac_LB_int, self.ac_UB_int, self.Lambda_int, self.Lambda_ex,
			self.Code_Lambda_int, self.Code_Lambda_ex, self.MAX_CODE_int, self.MAX_CODE_ex, self.rho, self.delta_i,
			self.Gamma_ac_lnk, self.Gamma_ac_lnk_less, self.Gamma_ac_lnk_equal,
			self.delta_alpha_C, self.delta_alpha_T, self.delta_alpha_F, self.delta_beta_C, self.delta_beta_T, self.delta_beta_F,
			self.delta_beta_plus, self.delta_beta_minus, self.delta_beta_in, self.delta_chi_T, self.delta_chi_F, self.chi_T, self.chi_F,
			self.delta_ac_C, self.delta_ac_T, self.delta_ac_F, self.delta_ac_CT, self.delta_ac_TC, self.delta_ac_CF, self.delta_ac_TF,
			self.e_C, self.e_T, self.e_F, self.v_T, self.v_F, self.Delta_ac_C_plus, self.Delta_ac_C_minus, self.Delta_ac_T_plus,
			self.Delta_ac_T_minus, self.Delta_ac_F_plus, self.Delta_ac_F_minus, self.Delta_ac_CT_plus, self.Delta_ac_CT_minus,
			self.Delta_ac_TC_plus, self.Delta_ac_TC_minus, self.Delta_ac_CF_plus, self.Delta_ac_CF_minus, self.Delta_ac_TF_plus,
			self.Delta_ac_TF_minus, self.delta_ac_T_lnk, self.clr_T, self.alpha_C, self.alpha_T, self.alpha_F, self.alpha_CT, self.alpha_TC, self.alpha_CF, self.alpha_TF, 
			self.beta_C, self.beta_T, self.beta_F, self.beta_plus, self.beta_minus, self.beta_in, 
			self.ac_C, self.ac_T, self.ac_F, self.ac_CT, self.ac_TC, self.ac_CF, self.ac_TF, self.ac_int,
			self.ac_lnk, self.ac_C_lnk, self.ac_T_lnk, self.ac_CT_lnk, self.ac_TC_lnk, self.I_lnk, self.n_lnk
		)
		MILP = self.add_constraints_chemical_symbols(
			MILP,
			self.t_C, self.t_T, self.t_F, self.n_C, self.n_T, self.n_F,
			self.Lambda_dg_int, self.Lambda_tilde_dg_C, self.Lambda_tilde_dg_T, self.Lambda_tilde_dg_F,
			self.Code_Lambda, self.Code_Lambda_int, self.Code_Lambda_ex, self.epsilon_dg, self.delta_i, self.rho,
			self.delta_ns_C, self.delta_ns_T, self.delta_ns_F,
			self.alpha_C, self.alpha_T, self.alpha_F, self.deg_C, self.deg_T, self.deg_F, self.ns_int
		)
		MILP = self.add_constraints_edge_configuration(
			MILP,
			self.t_C, self.t_C_tilde, self.t_T, self.t_F, self.n_C, self.n_T, self.n_F, self.k_C, self.k_C_tilde, self.m_C, self.c_F,
			self.tail_C, self.head_C, self.Gamma_int, self.Gamma_int_less, self.Gamma_int_equal,
			self.Gamma_tilde_ec_C, self.Gamma_tilde_ec_T, self.Gamma_tilde_ec_F, self.Gamma_tilde_ec_CT,
			self.Gamma_tilde_ec_TC, self.Gamma_tilde_ec_CF, self.Gamma_tilde_ec_TF,
			self.Gamma_tilde_ac_C, self.Gamma_tilde_ac_T, self.Gamma_tilde_ac_F, self.Gamma_tilde_ac_CT,
			self.Gamma_tilde_ac_TC, self.Gamma_tilde_ac_CF, self.Gamma_tilde_ac_TF,
			self.rho, self.delta_i, self.Code_Gamma_ac_int, self.Code_Gamma_ec_int,
			self.Gamma_lnk, self.Gamma_lnk_less, self.Gamma_lnk_equal, self.Lambda_dg_int,
			self.ns_LB_cnt, self.ns_UB_cnt, self.I_lnk, self.Gamma_cnt_less, self.Gamma_cnt_equal, self.Gamma_cnt_great,
			self.delta_dg_C, self.delta_dg_T, self.delta_dg_F, self.delta_chi_T, self.delta_chi_F, self.chi_T, self.chi_F,
			self.delta_beta_C, self.delta_beta_T, self.delta_beta_F, self.delta_beta_plus, self.delta_beta_minus, self.delta_beta_in,
			self.delta_ac_C, self.delta_ac_T, self.delta_ac_F, self.delta_ac_CT, self.delta_ac_TC, self.delta_ac_CF, self.delta_ac_TF,
			self.delta_ec_C, self.delta_ec_T, self.delta_ec_F, self.delta_ec_CT, self.delta_ec_TC, self.delta_ec_CF, self.delta_ec_TF,
			self.e_C, self.e_T, self.e_F, self.v_T, self.v_F, self.delta_ec_T_lnk, self.delta_cnt,
			self.clr_T, self.deg_C, self.deg_T, self.deg_F,
			self.deg_T_CT, self.deg_T_TC, self.deg_F_CF, self.deg_F_TF,
			self.ec_C, self.ec_T, self.ec_F, self.ec_CT, self.ec_TC, self.ec_CF, self.ec_TF,
			self.ec_int, self.Delta_ec_C_plus, self.Delta_ec_C_minus, self.Delta_ec_T_plus, self.Delta_ec_T_minus,
			self.Delta_ec_F_plus, self.Delta_ec_F_minus, self.Delta_ec_CT_plus, self.Delta_ec_CT_minus,
			self.Delta_ec_TC_plus, self.Delta_ec_TC_minus, self.Delta_ec_CF_plus, self.Delta_ec_CF_minus,
			self.Delta_ec_TF_plus, self.Delta_ec_TF_minus, self.ec_lnk, self.ec_T_lnk, self.ec_C_lnk,
			self.ec_CT_lnk, self.ec_TC_lnk, self.n_lnk
		)
		MILP = self.add_constraints_fringe_configuration(
			MILP,
			self.t_C, self.t_T, self.t_F, self.set_F, self.Code_F, 
			self.delta_fr_C, self.delta_fr_T, self.delta_fr_F, self.fc
		)
		MILP, self.delta_n, self.mass_n = self.add_constraints_mass_n(
			MILP, 
			self.n_LB, self.n_star, self.n_G, self.na_UB, self.na_ex, self.MASS
		)

		if _prevent_double_bound:
			# MILP = self.add_constraints_prevent_continuous_double_bond(
			# 	MILP,
			# 	self.t_C, self.t_T, self.t_F,
			# 	self.delta_alpha_C, self.delta_alpha_T, self.delta_alpha_F,
			# 	self.deg_C, self.deg_T, self.deg_F, self.hyddeg_C, self.hyddeg_T, self.hyddeg_F
			# )
			MILP = self.add_constraints_prevent_continuous_double_bond_complete(
				MILP,
				self.t_C, self.t_C_tilde, self.t_T, self.t_F, self.beta_r_2, self.Code_F, self.F_C, self.F_T, self.F_F, self.E_C,
				self.I_equal_one, self.I_zero_one, self.I_ge_one, self.I_ge_one_plus, self.I_ge_one_minus, self.I_ge_two_plus, self.I_ge_two_minus,
				self.delta_alpha_C, self.delta_alpha_T, self.delta_alpha_F, self.delta_fr_C, self.delta_fr_T, self.delta_fr_F,
				self.delta_beta_C, self.delta_beta_T, self.delta_beta_F, self.delta_beta_plus, self.delta_beta_minus, self.delta_beta_in,
				self.beta_CT, self.beta_TC, self.beta_CF, self.beta_TF
			)

		return MILP

	def print_sdf_file(self,
		t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, Lambda_int, Lambda_ex,
		head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, I_lnk, set_F, Code_F, E_C,
		# Variables
		n_G, v_T, v_F, alpha_C, alpha_T, alpha_F, beta_C, beta_T, beta_F,
		beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, delta_cnt,

		outputfilename
		):

		''' A function to output sdf file.  '''
		n_tmp = 1000
		n = 0

		graph_col = {i : " " for i in range(1, n_tmp + 1)}
		graph_adj = {(i, j) : 0 for i in range(1, n_tmp + 1) for j in range(1, n_tmp + 1)}

		graph_index = 0
		m = 0

		index_T = {(i, j) : 0 for i in range(1, t_T + 1) for j in range(n_T + 1)}
		index_C = {(i, j) : 0 for i in range(1, t_C + 1) for j in range(n_C + 1)}
		index_F = {(i, j) : 0 for i in range(1, t_F + 1) for j in range(n_F + 1)}

		chg = list()
		chg_ele = dict()

		for i in range(1, t_T + 1):
			if round(v_T[i].value()) != 0:
				graph_index += 1
				index_T[(i, 0)] = graph_index
				graph_col[graph_index] = Lambda_int[round(alpha_T[i].value()) - 1]

				psi = FringeTree()
				for psi_tmp in set_F:
					if round(delta_fr_T[i, Code_F[psi_tmp]].value()) != 0:
						psi = psi_tmp

				for j in range(1, len(psi.vertex)):
					a_tmp, h_tmp = psi.vertex[j]
					graph_index += 1
					index_T[(i, j)] = graph_index
					graph_col[graph_index] = a_tmp

				for j in range(0, len(psi.vertex)):
					if psi.chg[j] != 0:
						chg.append(index_T[i, j])
						chg.append(psi.chg[j])
						chg_ele[index_T[i, j]] = 4 - psi.chg[j]
					else:
						chg_ele[index_T[i, j]] = 0

				for j1 in range(len(psi.vertex)):
					for j2 in range(j1 + 1, len(psi.vertex)):
						if psi.beta[j1][j2] != 0:
							m += 1
							ind1 = index_T[(i, j1)]
							ind2 = index_T[(i, j2)]
							graph_adj[(ind1, ind2)] = psi.beta[j1][j2]
							graph_adj[(ind2, ind1)] = psi.beta[j1][j2]

		for i in range(1, t_C + 1):
			graph_index += 1
			index_C[(i, 0)] = graph_index
			graph_col[graph_index] = Lambda_int[round(alpha_C[i].value()) - 1]

			psi = FringeTree()
			for psi_tmp in set_F:
				if round(delta_fr_C[i, Code_F[psi_tmp]].value()) != 0:
					psi = psi_tmp

			for j in range(1, len(psi.vertex)):
				a_tmp, h_tmp = psi.vertex[j]
				graph_index += 1
				index_C[(i, j)] = graph_index
				graph_col[graph_index] = a_tmp

			for j in range(0, len(psi.vertex)):
				if psi.chg[j] != 0:
					chg.append(index_C[i, j])
					chg.append(psi.chg[j])
					chg_ele[index_C[i, j]] = 4 - psi.chg[j]
				else:
					chg_ele[index_C[i, j]] = 0

			for j1 in range(len(psi.vertex)):
				for j2 in range(j1 + 1, len(psi.vertex)):
					if psi.beta[j1][j2] != 0:
						m += 1
						ind1 = index_C[(i, j1)]
						ind2 = index_C[(i, j2)]
						graph_adj[(ind1, ind2)] = psi.beta[j1][j2]
						graph_adj[(ind2, ind1)] = psi.beta[j1][j2]

		for i in range(1, t_F + 1):
			if round(v_F[i].value()) != 0:
				graph_index += 1
				index_F[(i, 0)] = graph_index
				graph_col[graph_index] = Lambda_int[round(alpha_F[i].value()) - 1]

				psi = FringeTree()
				for psi_tmp in set_F:
					if round(delta_fr_F[i, Code_F[psi_tmp]].value()) != 0:
						psi = psi_tmp

				for j in range(1, len(psi.vertex)):
					a_tmp, h_tmp = psi.vertex[j]
					graph_index += 1
					index_F[(i, j)] = graph_index
					graph_col[graph_index] = a_tmp

				for j in range(0, len(psi.vertex)):
					if psi.chg[j] != 0:
						chg.append(index_F[i, j])
						chg.append(psi.chg[j])
						chg_ele[index_F[i, j]] = 4 - psi.chg[j]
					else:
						chg_ele[index_F[i, j]] = 0

				for j1 in range(len(psi.vertex)):
					for j2 in range(j1 + 1, len(psi.vertex)):
						if psi.beta[j1][j2] != 0:
							m += 1
							ind1 = index_F[(i, j1)]
							ind2 = index_F[(i, j2)]
							graph_adj[(ind1, ind2)] = psi.beta[j1][j2]
							graph_adj[(ind2, ind1)] = psi.beta[j1][j2]

		for i in range(2, t_T + 1):
			mul = round(beta_T[i].value())
			if mul != 0:
				m += 1
				ind1 = index_T[(i - 1, 0)]
				ind2 = index_T[(i, 0)]
				graph_adj[(ind1, ind2)] = mul
				graph_adj[(ind2, ind1)] = mul

		for i in range(2, t_F + 1):
			mul = round(beta_F[i].value())
			if mul != 0:
				m += 1
				ind1 = index_F[(i - 1, 0)]
				ind2 = index_F[(i, 0)]
				graph_adj[(ind1, ind2)] = mul
				graph_adj[(ind2, ind1)] = mul

		for i in (I_equal_one + I_zero_one + I_ge_one):
			mul = round(beta_C[i].value())
			if mul != 0:
				m += 1
				ind1 = index_C[(head_C[i], 0)]
				ind2 = index_C[(tail_C[i], 0)]
				graph_adj[(ind1, ind2)] = mul
				graph_adj[(ind2, ind1)] = mul

		for i in range(1, t_T + 1):
			mul = round(beta_CT[i].value())
			ei1 = round(e_T[i].value())
			ei2 = round(e_T[i + 1].value())
			if mul != 0 and ei1 == 0:
				for c in (I_ge_two + I_ge_one):
					if round(chi_T[(i, c)].value()) == 1:
						m += 1
						ind1 = index_C[(tail_C[c], 0)]
						ind2 = index_T[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul 
						break
		for i in range(1, t_T + 1):
			mul = round(beta_TC[i].value())
			ei1 = round(e_T[i].value())
			ei2 = round(e_T[i + 1].value())
			if mul != 0 and ei2 == 0:
				for c in (I_ge_two + I_ge_one):
					if round(chi_T[(i, c)].value()) == 1:
						m += 1
						ind1 = index_C[(head_C[c], 0)]
						ind2 = index_T[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul
						break
		for i in range(1, t_F + 1):
			mul = round(beta_CF[i].value())
			ei1 = round(e_F[i].value())
			if mul != 0 and ei1 == 0:
				for c in range(1, t_C_tilde + 1):
					if round(chi_F[(i, c)].value()) == 1:
						m += 1
						ind1 = index_C[(c, 0)]
						ind2 = index_F[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul 
						break
		for i in range(1, t_F + 1):
			mul = round(beta_TF[i].value())
			ei1 = round(e_F[i].value())
			if mul != 0 and ei1 == 0:
				for c in range(t_C_tilde + 1, c_F + 1):
					if round(chi_F[(i, c)].value()) == 1:
						m += 1
						ind1 = index_T[(c - t_C_tilde, 0)]
						ind2 = index_F[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul 
						break

		n = graph_index

		base_vertices = list()
		base_edges = list()
		base_edges_color = list()
		visited = set()
		core_set = set()
		link_edges = list()
		for i in range(1, t_C + 1):
			base_vertices.append(index_C[(i, 0)])
			core_set.add(index_C[(i, 0)])
		for i in range(1, t_T + 1):
			if index_T[(i, 0)] != 0:
				core_set.add(index_T[(i, 0)])    

		deg_tmp = {i: 0 for i in range(1, n + 1)}
		for i in range(1, n + 1):
			for j in range(1, n + 1):
				if i != j and graph_col[j] != "H1" and graph_adj[(i, j)] != 0:
					deg_tmp[i] += 1

		for i in range(t_C):
			u = base_vertices[i]
			visited.add(u)
			for v in range(1, n + 1):
				if graph_adj[(u, v)] != 0 and v in core_set:
					tmp = list()
					tmp.append(u)
					tmp.append(v)
					uu = u
					vv = v
					while vv not in visited and vv not in base_vertices:
						visited.add(vv)
						for ww in range(1, n + 1):
							if graph_adj[(vv, ww)] != 0 and ww != uu and ww in core_set:
								tmp.append(ww)
								uu = vv
								vv = ww
								break
					if len(tmp) > 2 or v not in visited:
						base_edges.append(tmp)
						if len(tmp) == 2:
							for c in (I_equal_one + I_zero_one + I_ge_one):
								e = E_C[c]
								if (index_C[(e[1], 0)] == tmp[0] and index_C[(e[2], 0)] == tmp[1]) or \
										(index_C[(e[2], 0)] == tmp[0] and index_C[(e[1], 0)] == tmp[1]):
									base_edges_color.append(c)
									break
						elif len(tmp) > 2:
							ind_T = 0
							for j in range(1, t_T + 1):
								if index_T[(j, 0)] == tmp[1]:
									ind_T = j
									break
							c = round(chi_T[ind_T].value())
							base_edges_color.append(c)
						if len(link_edges) == 0 and c in I_lnk:
							for i in range(len(tmp) - 1):
								uu = tmp[i]
								vv = tmp[i + 1]
								if graph_adj[(uu, vv)] == 1:
									ec1 = ((graph_col[uu], deg_tmp[uu]), (graph_col[vv], deg_tmp[vv]), graph_adj[(uu, vv)])
									ec2 = ((graph_col[vv], deg_tmp[vv]), (graph_col[uu], deg_tmp[uu]), graph_adj[(uu, vv)])
									if delta_cnt[ec1].value() > 1e-6 or delta_cnt[ec2].value() > 1e-6:
										link_edges.append((uu, vv))
										break

		n = graph_index
		with open(outputfilename, "w") as f:
			f.write("1\n")
			f.write("mol-infer\n")
			f.write("2LMM model for polymer\n")
			f.write("{:3}{:3}  0  0  0  0  0  0  0  0999 V2000 \n".format(n, m))
			
			for i in range(1, n + 1):
				atom_tmp = re.sub(r'[0-9]', '', graph_col[i])
				f.write("    0.0000    0.0000    0.0000 {:2}  0  {:1}  0  0  0  0  0  0  0  0  0  0\n".format(atom_tmp, chg_ele[i]))

			for i in range(1, n + 1):
				for j in range(i + 1, n + 1):
					if (i, j) not in link_edges and (j, i) not in link_edges and graph_adj[(i, j)] > 0:
						f.write("{:3}{:3}{:3}  0  0  0  0\n".format(i, j, graph_adj[(i, j)]))
			u, v = link_edges[0]
			f.write("{:3}{:3}{:3}  0  0  0  0\n".format(u, v, graph_adj[(u, v)]))

			if len(chg) != 0:
				f.write("M  CHG{:3}".format(int(len(chg) / 2)))
				for tmp in chg:
					f.write("{:4}".format(tmp))
				f.write("\n")

			f.write("M  END\n")
			f.write("$$$$\n")
			f.close()

		print(f"{outputfilename} generated.")

		return index_C, index_T, graph_adj, graph_index

	def print_sdf_file_with_e(self,
		t_C, t_T, t_F, t_C_tilde, n_C, n_T, n_F, c_F, Lambda, Lambda_int, Lambda_ex,
		head_C, tail_C, I_equal_one, I_zero_one, I_ge_one, I_ge_two, I_lnk, set_F, Code_F, E_C,
		# Variables
		n_G, v_T, v_F, alpha_C, alpha_T, alpha_F, beta_C, beta_T, beta_F,
		beta_CT, beta_TC, beta_CF, beta_TF, chi_T, chi_F, e_T, e_F, delta_fr_C, delta_fr_T, delta_fr_F, delta_cnt,

		outputfilename
		):

		''' A function to output sdf file.  '''
		n_tmp = 1000
		n = 0

		graph_col = {i : " " for i in range(1, n_tmp + 1)}
		graph_adj = {(i, j) : 0 for i in range(1, n_tmp + 1) for j in range(1, n_tmp + 1)}

		graph_index = 0
		m = 0

		index_T = {(i, j) : 0 for i in range(1, t_T + 1) for j in range(n_T + 1)}
		index_C = {(i, j) : 0 for i in range(1, t_C + 1) for j in range(n_C + 1)}
		index_F = {(i, j) : 0 for i in range(1, t_F + 1) for j in range(n_F + 1)}

		chg = list()
		chg_ele = dict()

		for i in range(1, t_T + 1):
			if round(v_T[i].value()) != 0:
				graph_index += 1
				index_T[(i, 0)] = graph_index
				graph_col[graph_index] = Lambda_int[round(alpha_T[i].value()) - 1]

				psi = FringeTree()
				for psi_tmp in set_F:
					if round(delta_fr_T[i, Code_F[psi_tmp]].value()) != 0:
						psi = psi_tmp

				for j in range(1, len(psi.vertex)):
					a_tmp, h_tmp = psi.vertex[j]
					graph_index += 1
					index_T[(i, j)] = graph_index
					graph_col[graph_index] = a_tmp

				for j in range(0, len(psi.vertex)):
					if psi.chg[j] != 0:
						chg.append(index_T[i, j])
						chg.append(psi.chg[j])
						chg_ele[index_T[i, j]] = 4 - psi.chg[j]
					else:
						chg_ele[index_T[i, j]] = 0

				for j1 in range(len(psi.vertex)):
					for j2 in range(j1 + 1, len(psi.vertex)):
						if psi.beta[j1][j2] != 0:
							m += 1
							ind1 = index_T[(i, j1)]
							ind2 = index_T[(i, j2)]
							graph_adj[(ind1, ind2)] = psi.beta[j1][j2]
							graph_adj[(ind2, ind1)] = psi.beta[j1][j2]

		for i in range(1, t_C + 1):
			graph_index += 1
			index_C[(i, 0)] = graph_index
			graph_col[graph_index] = Lambda_int[round(alpha_C[i].value()) - 1]

			psi = FringeTree()
			for psi_tmp in set_F:
				if round(delta_fr_C[i, Code_F[psi_tmp]].value()) != 0:
					psi = psi_tmp

			for j in range(1, len(psi.vertex)):
				a_tmp, h_tmp = psi.vertex[j]
				graph_index += 1
				index_C[(i, j)] = graph_index
				graph_col[graph_index] = a_tmp

			for j in range(0, len(psi.vertex)):
				if psi.chg[j] != 0:
					chg.append(index_C[i, j])
					chg.append(psi.chg[j])
					chg_ele[index_C[i, j]] = 4 - psi.chg[j]
				else:
					chg_ele[index_C[i, j]] = 0

			for j1 in range(len(psi.vertex)):
				for j2 in range(j1 + 1, len(psi.vertex)):
					if psi.beta[j1][j2] != 0:
						m += 1
						ind1 = index_C[(i, j1)]
						ind2 = index_C[(i, j2)]
						graph_adj[(ind1, ind2)] = psi.beta[j1][j2]
						graph_adj[(ind2, ind1)] = psi.beta[j1][j2]

		for i in range(1, t_F + 1):
			if round(v_F[i].value()) != 0:
				graph_index += 1
				index_F[(i, 0)] = graph_index
				graph_col[graph_index] = Lambda_int[round(alpha_F[i].value()) - 1]

				psi = FringeTree()
				for psi_tmp in set_F:
					if round(delta_fr_F[i, Code_F[psi_tmp]].value()) != 0:
						psi = psi_tmp

				for j in range(1, len(psi.vertex)):
					a_tmp, h_tmp = psi.vertex[j]
					graph_index += 1
					index_F[(i, j)] = graph_index
					graph_col[graph_index] = a_tmp

				for j in range(0, len(psi.vertex)):
					if psi.chg[j] != 0:
						chg.append(index_F[i, j])
						chg.append(psi.chg[j])
						chg_ele[index_F[i, j]] = 4 - psi.chg[j]
					else:
						chg_ele[index_F[i, j]] = 0

				for j1 in range(len(psi.vertex)):
					for j2 in range(j1 + 1, len(psi.vertex)):
						if psi.beta[j1][j2] != 0:
							m += 1
							ind1 = index_F[(i, j1)]
							ind2 = index_F[(i, j2)]
							graph_adj[(ind1, ind2)] = psi.beta[j1][j2]
							graph_adj[(ind2, ind1)] = psi.beta[j1][j2]

		graph_index += 1
		graph_col[graph_index] = "R"
		chg_ele[graph_index] = 0
		graph_index += 1
		graph_col[graph_index] = "R"
		chg_ele[graph_index] = 0

		for i in range(2, t_T + 1):
			mul = round(beta_T[i].value())
			if mul != 0:
				m += 1
				ind1 = index_T[(i - 1, 0)]
				ind2 = index_T[(i, 0)]
				graph_adj[(ind1, ind2)] = mul
				graph_adj[(ind2, ind1)] = mul

		for i in range(2, t_F + 1):
			mul = round(beta_F[i].value())
			if mul != 0:
				m += 1
				ind1 = index_F[(i - 1, 0)]
				ind2 = index_F[(i, 0)]
				graph_adj[(ind1, ind2)] = mul
				graph_adj[(ind2, ind1)] = mul

		for i in (I_equal_one + I_zero_one + I_ge_one):
			mul = round(beta_C[i].value())
			if mul != 0:
				m += 1
				ind1 = index_C[(head_C[i], 0)]
				ind2 = index_C[(tail_C[i], 0)]
				graph_adj[(ind1, ind2)] = mul
				graph_adj[(ind2, ind1)] = mul

		for i in range(1, t_T + 1):
			mul = round(beta_CT[i].value())
			ei1 = round(e_T[i].value())
			ei2 = round(e_T[i + 1].value())
			if mul != 0 and ei1 == 0:
				for c in (I_ge_two + I_ge_one):
					if round(chi_T[(i, c)].value()) == 1:
						m += 1
						ind1 = index_C[(tail_C[c], 0)]
						ind2 = index_T[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul 
						break
		for i in range(1, t_T + 1):
			mul = round(beta_TC[i].value())
			ei1 = round(e_T[i].value())
			ei2 = round(e_T[i + 1].value())
			if mul != 0 and ei2 == 0:
				for c in (I_ge_two + I_ge_one):
					if round(chi_T[(i, c)].value()) == 1:
						m += 1
						ind1 = index_C[(head_C[c], 0)]
						ind2 = index_T[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul
						break
		for i in range(1, t_F + 1):
			mul = round(beta_CF[i].value())
			ei1 = round(e_F[i].value())
			if mul != 0 and ei1 == 0:
				for c in range(1, t_C_tilde + 1):
					if round(chi_F[(i, c)].value()) == 1:
						m += 1
						ind1 = index_C[(c, 0)]
						ind2 = index_F[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul 
						break
		for i in range(1, t_F + 1):
			mul = round(beta_TF[i].value())
			ei1 = round(e_F[i].value())
			if mul != 0 and ei1 == 0:
				for c in range(t_C_tilde + 1, c_F + 1):
					if round(chi_F[(i, c)].value()) == 1:
						m += 1
						ind1 = index_T[(c - t_C_tilde, 0)]
						ind2 = index_F[(i, 0)]
						graph_adj[(ind1, ind2)] = mul
						graph_adj[(ind2, ind1)] = mul 
						break

		n = graph_index

		base_vertices = list()
		base_edges = list()
		base_edges_color = list()
		visited = set()
		core_set = set()
		link_edges = list()
		for i in range(1, t_C + 1):
			base_vertices.append(index_C[(i, 0)])
			core_set.add(index_C[(i, 0)])
		for i in range(1, t_T + 1):
			if index_T[(i, 0)] != 0:
				core_set.add(index_T[(i, 0)])    

		deg_tmp = {i: 0 for i in range(1, n + 1)}
		for i in range(1, n + 1):
			for j in range(1, n + 1):
				if i != j and graph_col[j] != "H1" and graph_adj[(i, j)] != 0:
					deg_tmp[i] += 1

		for i in range(t_C):
			u = base_vertices[i]
			visited.add(u)
			for v in range(1, n + 1):
				if graph_adj[(u, v)] != 0 and v in core_set:
					tmp = list()
					tmp.append(u)
					tmp.append(v)
					uu = u
					vv = v
					while vv not in visited and vv not in base_vertices:
						visited.add(vv)
						for ww in range(1, n + 1):
							if graph_adj[(vv, ww)] != 0 and ww != uu and ww in core_set:
								tmp.append(ww)
								uu = vv
								vv = ww
								break
					if len(tmp) > 2 or v not in visited:
						base_edges.append(tmp)
						if len(tmp) == 2:
							for c in (I_equal_one + I_zero_one + I_ge_one):
								e = E_C[c]
								if (index_C[(e[1], 0)] == tmp[0] and index_C[(e[2], 0)] == tmp[1]) or \
										(index_C[(e[2], 0)] == tmp[0] and index_C[(e[1], 0)] == tmp[1]):
									base_edges_color.append(c)
									break
						elif len(tmp) > 2:
							ind_T = 0
							for j in range(1, t_T + 1):
								if index_T[(j, 0)] == tmp[1]:
									ind_T = j
									break
							c = round(chi_T[ind_T].value())
							base_edges_color.append(c)
						if len(link_edges) == 0 and c in I_lnk:
							for i in range(len(tmp) - 1):
								uu = tmp[i]
								vv = tmp[i + 1]
								if graph_adj[(uu, vv)] == 1:
									ec1 = ((graph_col[uu], deg_tmp[uu]), (graph_col[vv], deg_tmp[vv]), graph_adj[(uu, vv)])
									ec2 = ((graph_col[vv], deg_tmp[vv]), (graph_col[uu], deg_tmp[uu]), graph_adj[(uu, vv)])
									if delta_cnt[ec1].value() > 1e-6 or delta_cnt[ec2].value() > 1e-6:
										link_edges.append((uu, vv))
										graph_adj[(uu, graph_index - 1)] = 1
										graph_adj[(graph_index - 1, uu)] = 1
										graph_adj[(vv, graph_index)] = 1
										graph_adj[(graph_index, vv)] = 1
										m += 1
										break

		n = graph_index
		with open(outputfilename, "w") as f:
			f.write("1\n")
			f.write("mol-infer\n")
			f.write("2LMM model for polymer\n")
			f.write("{:3}{:3}  0  0  0  0  0  0  0  0999 V2000 \n".format(n, m))
			
			for i in range(1, n + 1):
				atom_tmp = re.sub(r'[0-9]', '', graph_col[i])
				f.write("    0.0000    0.0000    0.0000 {:2}  0  {:1}  0  0  0  0  0  0  0  0  0  0\n".format(atom_tmp, chg_ele[i]))

			for i in range(1, n + 1):
				for j in range(i + 1, n + 1):
					if (i, j) not in link_edges and (j, i) not in link_edges and graph_adj[(i, j)] > 0:
						f.write("{:3}{:3}{:3}  0  0  0  0\n".format(i, j, graph_adj[(i, j)]))

			if len(chg) != 0:
				f.write("M  CHG{:3}".format(int(len(chg) / 2)))
				for tmp in chg:
					f.write("{:4}".format(tmp))
				f.write("\n")

			f.write("M  END\n")
			f.write("$$$$\n")
			f.close()

		print(f"{outputfilename} generated.")

		return index_C, index_T, graph_adj, graph_index

	def _print_sdf(self, sdf_filename_prefix):
		sdf_filename_no_e = f"{sdf_filename_prefix}.sdf"
		self.index_C, self.index_T, self.graph_adj, self.graph_ind = self.print_sdf_file(
			self.t_C, self.t_T, self.t_F, self.t_C_tilde, self.n_C, self.n_T, self.n_F, self.c_F, 
			self.Lambda, self.Lambda_int, self.Lambda_ex,
			self.head_C, self.tail_C, self.I_equal_one, self.I_zero_one, self.I_ge_one, self.I_ge_two, 
			self.I_lnk, self.set_F, self.Code_F, self.E_C,
			self.n_G, self.v_T, self.v_F, self.alpha_C, self.alpha_T, self.alpha_F,
			self.beta_C, self.beta_T, self.beta_F, self.beta_CT, self.beta_TC, self.beta_CF, self.beta_TF, self.chi_T, self.chi_F,
			self.e_T, self.e_F, self.delta_fr_C, self.delta_fr_T, self.delta_fr_F, self.delta_cnt,
			sdf_filename_no_e
		)

		sdf_filename_with_e = f"{sdf_filename_prefix}_e.sdf"
		self.index_C, self.index_T, self.graph_adj, self.graph_ind = self.print_sdf_file_with_e(
			self.t_C, self.t_T, self.t_F, self.t_C_tilde, self.n_C, self.n_T, self.n_F, self.c_F, 
			self.Lambda, self.Lambda_int, self.Lambda_ex,
			self.head_C, self.tail_C, self.I_equal_one, self.I_zero_one, self.I_ge_one, self.I_ge_two, 
			self.I_lnk, self.set_F, self.Code_F, self.E_C,
			self.n_G, self.v_T, self.v_F, self.alpha_C, self.alpha_T, self.alpha_F,
			self.beta_C, self.beta_T, self.beta_F, self.beta_CT, self.beta_TC, self.beta_CF, self.beta_TF, self.chi_T, self.chi_F,
			self.e_T, self.e_F, self.delta_fr_C, self.delta_fr_T, self.delta_fr_F, self.delta_cnt,
			sdf_filename_with_e
		)

		return [sdf_filename_no_e, sdf_filename_with_e]

	def _print_gstar_file(self, partition_filename,
		graph_ind, chi_T,
		t_C, t_T, index_C, index_T, graph_adj, E_C, ch_LB, ch_UB, I_ge_one, I_ge_two, I_zero_one, I_equal_one,
		set_F_v, set_F_E
	):
		n = graph_ind

		base_vertices = list()
		base_edges = list()
		base_edges_color = list()
		visited = set()
		core_set = set()
		for i in range(1, t_C + 1):
			base_vertices.append(index_C[(i, 0)])
			core_set.add(index_C[(i, 0)])
		for i in range(1, t_T + 1):
			if index_T[(i, 0)] != 0:
				core_set.add(index_T[(i, 0)])    
		for i in range(t_C):
			u = base_vertices[i]
			visited.add(u)
			for v in range(1, n + 1):
				if graph_adj[(u, v)] != 0 and v in core_set:
					tmp = list()
					tmp.append(u)
					tmp.append(v)
					uu = u
					vv = v
					while vv not in visited and vv not in base_vertices:
						visited.add(vv)
						for ww in range(1, n + 1):
							if graph_adj[(vv, ww)] != 0 and ww != uu and ww in core_set:
								tmp.append(ww)
								uu = vv
								vv = ww
								break
					if len(tmp) > 2 or v not in visited:
						base_edges.append(tmp)
						if len(tmp) == 2:
							for c in (I_equal_one + I_zero_one + I_ge_one):
								e = E_C[c]
								if (index_C[(e[1], 0)] == tmp[0] and index_C[(e[2], 0)] == tmp[1]) or \
										(index_C[(e[2], 0)] == tmp[0] and index_C[(e[1], 0)] == tmp[1]):
									base_edges_color.append(c)
									# print(c)
									break
						elif len(tmp) > 2:
							ind_T = 0
							for j in range(1, t_T + 1):
								if index_T[(j, 0)] == tmp[1]:
									ind_T = j
									break
							c = round(chi_T[ind_T].value())
							base_edges_color.append(c)

		with open(partition_filename, "w") as f:
			f.write(f"{t_C}\n")
			for i in range(t_C):
				f.write(f"{base_vertices[i]}\n")
				f.write(f"{ch_LB[i + 1]} {ch_UB[i + 1]}\n")
				for j in set_F_v[i + 1]:
					f.write(f"{j.index} ")
				f.write("\n")

			f.write(f"{len(base_edges)}\n")
			for i in range(len(base_edges)):
				f.write(f"{base_edges[i][0]}")
				for j in range(1, len(base_edges[i])):
					f.write(f" {base_edges[i][j]}")
				f.write("\n")
				c = base_edges_color[i]
				if c in (I_ge_two + I_ge_one):
					f.write(f"{ch_LB[E_C[c]]} {ch_UB[E_C[c]]}\n")
				else:
					f.write("0 0\n")
				for j in set_F_E:
					f.write(f"{j.index} ")
				f.write("\n")

	def print_partition_file(self, partition_filename):
		if not hasattr(self, 'graph_ind'):
			raise ValueError("Please run print_sdf() before print_partition_file().")

		self._print_gstar_file(partition_filename,
			self.graph_ind, self.chi_T,
			self.t_C, self.t_T, self.index_C, self.index_T, self.graph_adj, self.E_C, self.ch_LB, self.ch_UB,
			self.I_ge_one, self.I_ge_two, self.I_zero_one, self.I_equal_one,
			self.set_F_v, self.set_F_E
		)







