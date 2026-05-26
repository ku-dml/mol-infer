from dataclasses import dataclass, field
import itertools
from pathlib import Path
from typing import Dict, List, Literal, Optional, ClassVar
from abc import ABC, abstractmethod
import yaml, time
import pandas as pd
import numpy as np
import pulp_modified as pulp

def _opt_path(d, key):
	v = d.get(key)
	return Path(v).expanduser().resolve() if v else None

def _normalize_config_input(v, _type=str):
	if isinstance(v, _type):
		return v
	if isinstance(v, (list, tuple)):
		return list(v)
	raise TypeError(f"Invalid input: {v!r}")

@dataclass(frozen=True)
class AssignmentLoader:
	_counter = itertools.count(1)
	ind: int = field(init=False)
	type: str
	seed_graph: str | List[str]
	learning_model: str
	lb: float
	ub: float
	extra: dict[str, object] = field(default_factory=dict)

	def __post_init__(self):
		object.__setattr__(self, "ind", next(self._counter))

	@classmethod
	def from_dict(cls, d):
		d = dict(d)
		return cls(
			seed_graph=_normalize_config_input(d["seed_graph"], _type=str), 
			type=d.get("type", "NORMAL"),
			learning_model=d["learning_model"],
			lb=float(d["lb"]), 
			ub=float(d["ub"]),
			extra=d
		)

@dataclass(frozen=True)
class SeedGraphLoader:
	name: str
	type: str
	filename: Path
	fringe_filename: Path
	assignment: List[AssignmentLoader]

	@classmethod
	def from_dict(cls, d):
		return cls(
			name=d["name"], 
			type=d["type"], 
			filename=Path(d["filename"]), 
			fringe_filename=Path(d["fringe_filename"]),
			assignment=list()
		)

	def validate_files(self):
		_files_to_check = [self.filename, self.fringe_filename]
		for p in _files_to_check:
			if not p.exists():
				raise FileNotFoundError(f"[SeedGraph:{self.name}] Not found: {p}")

	def __post_init__(self):
		self.validate_files()

@dataclass(frozen=False, kw_only=True)
class LearningModelLoader(ABC):
	name: str
	_type: ClassVar[Literal["BASE"]] = "BASE"
	fv_gen_path: Optional[Path] = None
	sdf_filename: Optional[Path] = None
	desc_filename: Optional[Path] = None
	desc_norm_selected_filename: Optional[Path] = None
	fringe_filename: Optional[Path] = None
	values_filename: Optional[Path] = None

	@abstractmethod
	def validate_files(self):
		pass

	@property
	def type(self):
		return type(self)._type

	def __post_init__(self):
		self.validate_files()

@dataclass(frozen=False, kw_only=True)
class MLRModelLoader(LearningModelLoader):
	_type: ClassVar[Literal["MLR"]] = "MLR"
	linreg_filename: Path
	_use_quadratic: bool = False

	@classmethod
	def from_dict(cls, d):
		flags = d.get("flags", []) or []

		return cls(name=d["name"], 
			desc_filename=Path(d["desc_filename"]).expanduser().resolve(), 
			desc_norm_selected_filename=Path(d["desc_norm_selected_filename"]).expanduser().resolve(),
			fringe_filename=Path(d["fringe_filename"]).expanduser().resolve(),
			values_filename=Path(d["values_filename"]).expanduser().resolve(),
			linreg_filename=Path(d["linreg_filename"]).expanduser().resolve(),
			fv_gen_path=_opt_path(d, "fv_gen_path"),
			sdf_filename=_opt_path(d, "sdf_filename"),
			_use_quadratic=("use_quadratic" in flags)
		)

	def validate_files(self):
		_files_to_check = [self.desc_filename, self.desc_norm_selected_filename, self.fringe_filename, self.values_filename, self.linreg_filename]
		if self.sdf_filename is not None:
			_files_to_check.append(self.sdf_filename)

		missing = [p for p in _files_to_check if p is None or not p.exists()]
		if missing:
			names = ", ".join(str(p) for p in missing if p is not None)
			raise FileNotFoundError(f"[LearningModel:{self.name}] Not found: {names}")

@dataclass(frozen=False, kw_only=True)
class ANNModelLoader(LearningModelLoader):
	_type: ClassVar[Literal["ANN"]] = "ANN"
	weight_filename: Path
	biases_filename: Path
	_use_quadratic: bool = False

	@classmethod
	def from_dict(cls, d):
		flags = d.get("flags", []) or []

		return cls(name=d["name"], 
			desc_filename=Path(d["desc_filename"]).expanduser().resolve(),
			desc_norm_selected_filename=Path(d["desc_norm_selected_filename"]).expanduser().resolve(),
			fringe_filename=Path(d["fringe_filename"]).expanduser().resolve(),
			values_filename=Path(d["values_filename"]).expanduser().resolve(),
			weight_filename=Path(d["weight_filename"]).expanduser().resolve(), 
			biases_filename=Path(d["biases_filename"]).expanduser().resolve(),
			fv_gen_path=_opt_path(d, "fv_gen_path"),
			sdf_filename=_opt_path(d, "sdf_filename"),
			_use_quadratic=("use_quadratic" in flags)
		)

	def validate_files(self):
		_files_to_check = [self.desc_filename, self.desc_norm_selected_filename, self.fringe_filename, self.values_filename, self.weight_filename, self.biases_filename]
		if self.sdf_filename is not None:
			_files_to_check.append(self.sdf_filename)
		
		missing = [p for p in _files_to_check if p is None or not p.exists()]
		if missing:
			names = ", ".join(str(p) for p in missing if p is not None)
			raise FileNotFoundError(f"[LearningModel:{self.name}] Not found: {names}")


@dataclass(frozen=False, kw_only=True)
class RFModelLoader(LearningModelLoader):
	_type: ClassVar[Literal["RF"]] = "RF"
	rf_filename: Path
	_use_quadratic: bool = False

	@classmethod
	def from_dict(cls, d):
		flags = d.get("flags", []) or []

		return cls(name=d["name"], 
			desc_filename=Path(d["desc_filename"]).expanduser().resolve(),
			desc_norm_selected_filename=Path(d["desc_norm_selected_filename"]).expanduser().resolve(),
			fringe_filename=Path(d["fringe_filename"]).expanduser().resolve(),
			values_filename=Path(d["values_filename"]).expanduser().resolve(),
			rf_filename=Path(d["rf_filename"]).expanduser().resolve(),
			fv_gen_path=_opt_path(d, "fv_gen_path"),
			sdf_filename=_opt_path(d, "sdf_filename"),
			_use_quadratic=("use_quadratic" in flags)
		)

	def validate_files(self):
		_files_to_check = [self.desc_filename, self.desc_norm_selected_filename, self.fringe_filename, self.values_filename, self.rf_filename]
		if self.sdf_filename is not None:
			_files_to_check.append(self.sdf_filename)

		missing = [p for p in _files_to_check if p is None or not p.exists()]
		if missing:
			names = ", ".join(str(p) for p in missing if p is not None)
			raise FileNotFoundError(f"[LearningModel:{self.name}] Not found: {names}")

@dataclass(frozen=False, kw_only=True)
class GNNModelLoader(LearningModelLoader):
	_type: ClassVar[Literal["GNN"]] = "GNN"
	model_filename: Path
	config_filename: Path
	_use_quadratic: bool = False

	@classmethod
	def from_dict(cls, d):
		flags = d.get("flags", []) or []

		return cls(name=d["name"], 
			desc_filename=None,
			desc_norm_selected_filename=None,
			fringe_filename=None,
			values_filename=None,
			model_filename=Path(d["model_filename"]).expanduser().resolve(),
			config_filename=Path(d["config_filename"]).expanduser().resolve(),
			_use_quadratic=("use_quadratic" in flags)
		)

	def validate_files(self):
		_files_to_check = [self.model_filename, self.config_filename]
		if self.sdf_filename is not None:
			_files_to_check.append(self.sdf_filename)

		missing = [p for p in _files_to_check if p is None or not p.exists()]
		if missing:
			names = ", ".join(str(p) for p in missing if p is not None)
			raise FileNotFoundError(f"[LearningModel:{self.name}] Not found: {names}")

@dataclass(frozen=True)
class SolverLoader:
	type: str

@dataclass(frozen=True)
class CPLEXSolverLoader(SolverLoader):
	type: Literal["CPLEX"]
	CPLEX_PATH: Path
	CPLEX_MSG: bool
	CPLEX_TIMELIMIT: float

	@classmethod
	def from_dict(cls, d):
		return cls(type="CPLEX", 
			CPLEX_PATH=Path(d["CPLEX_PATH"]), 
			CPLEX_MSG=d.get("CPLEX_MSG", False), 
			CPLEX_TIMELIMIT=d.get("CPLEX_TIMELIMIT", 0.0)
		)

# @dataclass(frozen=True)
# class GurobiSolverLoader(SolverLoader):
# # 	type: Literal["Gurobi"]
	
# 	@classmethod
# 	def from_dict(cls, d):


class ConfigLoader:
	def __init__(self, path: Path):
		self.path = Path(path)
		with self.path.open("r") as f:
			self.raw = yaml.safe_load(f) or {}

		self._debug = self.raw.get("DEBUG", False)

		if "solver" not in self.raw:
			raise ValueError(f"Solver information not found in {self.path}!")

		s = self.raw.get("solver", [])
		if s["type"] == "CPLEX":
			self.solver = CPLEXSolverLoader.from_dict(s)
		else:
			raise ValueError(f"[Solver:{s['type']}] Not Support!")

		if "seed_graphs" not in self.raw or len(self.raw.get("seed_graphs")) == 0:
			raise ValueError(f"Seed graph information not found in {self.path}!")

		self.seed_graphs_list = [SeedGraphLoader.from_dict(g) for g in self.raw.get("seed_graphs", [])]

		if "learning_models" not in self.raw or len(self.raw.get("learning_models")) == 0:
			raise ValueError(f"Learning model information not found in {self.path}!")		

		self.learning_models_list = list()
		for g in self.raw.get("learning_models", []):
			if g["type"] == "mlr":
				self.learning_models_list.append(MLRModelLoader.from_dict(g))
			elif g["type"] == "ann":
				self.learning_models_list.append(ANNModelLoader.from_dict(g))
			elif g["type"] == "rf":
				self.learning_models_list.append(RFModelLoader.from_dict(g))
			elif g["type"] == "gnn":
				self.learning_models_list.append(GNNModelLoader.from_dict(g))
			else:
				# not supporting type
				raise ValueError(f"[LearningModel:{g['name']}, {g['type']}] Not Support!")

		self.graphs = {g.name: g for g in self.seed_graphs_list}
		self.learning_models = {g.name: g for g in self.learning_models_list}

		if "assignments" not in self.raw or len(self.raw.get("assignments")) == 0:
			raise ValueError(f"Assignment information not found in {self.path}!")			

		self.assignments_list = [AssignmentLoader.from_dict(a) for a in self.raw.get("assignments", [])]

		# self._validate_files()
		# self._validate()
				
	def exp_start(self):
		self.start_time = time.time()

	def exp_end(self):
		self.end_time = time.time()

	def init_end(self):
		self.init_end_time = time.time()

	def prepare_solver(self):
		if self.solver.type == "CPLEX":
			if self.solver.CPLEX_TIMELIMIT > 0:
				CPLEX = pulp.CPLEX(path=self.solver.CPLEX_PATH, msg=self.solver.CPLEX_MSG, timeLimit=self.solver.CPLEX_TIMELIMIT)
			else:
				CPLEX = pulp.CPLEX(path=self.solver.CPLEX_PATH, msg=self.solver.CPLEX_MSG)
		return CPLEX

	def _validate(self):
		for a in self.assignments_list:
			if a.seed_graph not in self.graphs:
				raise ValueError(f"Undefined seed graph: {a.seed_graph}")
			if a.learning_model not in self.learning_models:
				raise ValueError(f"Undefined learning model: {a.learning_model}")
			self.graphs[a.seed_graph].assignment.append(a)

	def __post_init__(self):
		self._validate()


