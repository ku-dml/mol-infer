from abc import ABC, abstractmethod
from typing import Dict, List, Literal, Optional, ClassVar
from dataclasses import dataclass, field
import pandas as pd
import numpy as np

@dataclass(frozen=False, kw_only=True)
class InverterBase(ABC):
	_type: ClassVar[Literal["BASE"]] = "BASE"

	@property
	def type(self):
		return type(self)._type

	@abstractmethod
	def set_up_inverter(self):
		pass

	@abstractmethod
	def set_up_variables(self):
		pass

	@abstractmethod
	def set_up_x_hat_variables(self):
		pass

	@abstractmethod
	def add_constraints_C1(self):
		pass

	@abstractmethod
	def inspection(self):
		pass

	def prepare_fv_list(self, desc_norm_selected_filename):
		fv = pd.read_csv(desc_norm_selected_filename, sep=",")
		self.fv_list = list(fv.columns)[1:] ### remove "CID" column

	def prepare_desc_max_min(self, desc_filename):
		fv = pd.read_csv(desc_filename, sep=",")
		# prepare for normalization and standardization
		self.max_dcp = dict()
		self.min_dcp = dict()

		for fv_name in fv.columns:
			if fv_name == "CID":
				continue
			s = fv[fv_name]
			self.max_dcp[fv_name] = float(s.max())
			self.min_dcp[fv_name] = float(s.min())

	def prepare_value_max_min(self, desc_norm_selected_filename, values_filename):
		fv = pd.read_csv(desc_norm_selected_filename)
		value = pd.read_csv(values_filename)

		CIDs = np.array(fv['CID'])
		target = np.array(value['a'])
		fv_dict = {}
		for cid,row in zip(CIDs, fv.values[:,1:]):
			fv_dict[cid] = row
		target_dict = {}
		for cid, val in zip(np.array(value['CID']), np.array(value['a'])):
			target_dict[cid] = val
		for cid in CIDs:
			if cid not in target_dict:
				raise ValueError(f"Error: {values_filename} misses the value of CID {cid}")

		y = np.array([target_dict[cid] for cid in CIDs])

		self.y_min = np.amin(y)
		self.y_max = np.amax(y)


