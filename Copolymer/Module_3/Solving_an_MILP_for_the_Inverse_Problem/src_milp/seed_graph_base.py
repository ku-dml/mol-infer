from collections import namedtuple
from dataclasses import dataclass
from abc import ABC, abstractmethod
from typing import Dict, List, Literal, Optional, ClassVar
from rdkit import Chem
from rdkit.Chem import AllChem
import os, tempfile, shutil, subprocess
from pathlib import Path

@dataclass(frozen=True)
class CG_element:
	symbol: str
	valence: int
	mass: int
	ele_neg: float
	ion_energy: float
	ele_affinity: float

def prepare_CG_element_info():
	''' A function to prepare information of chemical elements. '''

	set_Lambda = list()
	set_Lambda.append(CG_element("B", 3, 108, 2.0, 0.801, 0.24))
	set_Lambda.append(CG_element("C", 4, 120, 2.6, 1.086, 1.27))
	set_Lambda.append(CG_element("O", 2, 160, 3.4, 1.314, 1.465))
	set_Lambda.append(CG_element("N", 3, 140, 3.0, 1.402, 0.000))
	set_Lambda.append(CG_element("F", 1, 190, 4.0, 1.681, 3.34))
	set_Lambda.append(CG_element("Si", 4, 280, 1.9, 0.787, 1.24))
	set_Lambda.append(CG_element("P", 4, 309, 2.2, 1.012, 0.77))
	set_Lambda.append(CG_element("S", 2, 320, 2.6, 1.000, 2.08))
	set_Lambda.append(CG_element("Cl", 1, 354, 3.2, 1.251, 3.61))
	set_Lambda.append(CG_element("V", 3, 509, 1.6, 0.651, -0.051))
	set_Lambda.append(CG_element("Br", 1, 800, 3.0, 1.140, 3.36))
	set_Lambda.append(CG_element("Cd", 2, 1124, 1.7, 0.868, 0))
	set_Lambda.append(CG_element("I", 1, 1270, 2.7, 1.008, 0.35))
	set_Lambda.append(CG_element("Hg", 2, 2006, 1.9, 1.007, 0))
	set_Lambda.append(CG_element("Pb", 2, 2072, 1.8, 0.716, 0))
	set_Lambda.append(CG_element("Al", 3, 270, 1.6, 0.578, 0.46))
	set_Lambda.append(CG_element("H", 1, 10, 2.2, 1.312, 0.754))
	set_Lambda.append(CG_element("e*", 2, 0, 0, 0, 0))

	return set_Lambda

class FringeTree():
	def __init__(self):
		self.root = ("e*", 0)
		self.index = 0
		self.vertex = []
		self.adj = []
		self.alpha = []
		self.beta = []
		self.height = 0
		self.chg = []

	def read_from_seq(self, str1, str2, str3):
		seq1 = str1.split()
		seq2 = [int(mul) for mul in str2.split()]
		seq3 = [int(chg) for chg in str3.split()]
		# self.index = ind
		self.vertex = [(seq1[j], int(seq1[j + 1])) for j in range(0, len(seq1), 2)]
		self.root = self.vertex[0]
		self.height = max(self.vertex[v][1] for v in range(len(self.vertex)) if self.vertex[v][0] != "H1")
		self.adj = [set() for _ in range(len(self.vertex))]
		self.beta = [[0 for _ in range(len(self.vertex))] for _ in range(len(self.vertex))]
		self.chg = [chg for chg in seq3]
		for j in range(len(seq2)):
			cld = j + 1
			prt = max(v for v in range(j + 1) if self.vertex[v][1] == self.vertex[cld][1] - 1)   
			self.adj[prt].add(cld)
			self.adj[cld].add(prt)
			self.beta[prt][cld] = seq2[j]
			self.beta[cld][prt] = seq2[j]
			# print(str(prt) + " " + str(cld) + " " + str(j) + " " + str(seq2[j]))


@dataclass(frozen=False, kw_only=True)
class SeedGraphBase(ABC):
	_type: ClassVar[Literal["BASE"]] = "BASE"

	@property
	def type(self):
		return type(self)._type
	
	@abstractmethod
	def read_seed_graph(self):
		pass

	@abstractmethod
	def read_fringe_trees(self):
		pass

	@abstractmethod
	def set_up_seed_graph(self):
		pass

	@abstractmethod
	def set_up_variables(self):
		pass

	@abstractmethod
	def add_constraints_C2(self):
		pass

	@abstractmethod
	def _print_sdf(self):
		pass

	def print_sdf(self, *args, **kwargs):
		paths = list(self._print_sdf(*args, **kwargs))
		self._post_print_sdf(paths)

	def _post_print_sdf(self, sdf_filename_list):
		pass

class Rdkit2DPostprocessMixin:
	### A mixin to add a postprocess when outputing sdf files
	def _post_print_sdf(self, sdf_filename_list):
		for p in sdf_filename_list:
			# self._compute_2d_inplace(p)
			self._sdf_to_svg_obabel(p)

	@staticmethod
	def _compute_2d_inplace(sdf_path):
		suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
		fd, tmp_path = tempfile.mkstemp(suffix=".sdf", prefix="rdk2d_")
		os.close(fd)

		wrote = 0
		try:
			writer = Chem.SDWriter(tmp_path)
			writer.SetKekulize(True)
			for mol in suppl:
				if mol is None:
					continue
				flags = (Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
				Chem.SanitizeMol(mol, sanitizeOps=flags)
				
				AllChem.Compute2DCoords(mol)
				writer.write(mol)
				wrote += 1
			writer.close()
			os.replace(tmp_path, sdf_path)
		except Exception:
			try: os.remove(tmp_path)
			except OSError: pass
			raise
		return wrote

	@staticmethod
	def _sdf_to_svg_obabel(sdf_path, obabel_cmd=None, overwrite=True):
		sdf_path = Path(sdf_path)
		svg_path = sdf_path.with_suffix(".svg")

		if obabel_cmd is None:
			obabel_cmd = shutil.which("obabel") or shutil.which("obabel.exe")
		if obabel_cmd is None:
			print(f"Cannot find Obabel. Skip.")
			return None

		if svg_path.exists():
			if overwrite:
				try:
					svg_path.unlink()
				except OSError:
					pass
			else:
				return str(svg_path)
		
		cmd = [
			obabel_cmd,
			"-i", "sdf", str(sdf_path),
			"-o", "svg",
			"-O", str(svg_path),
			"-d",
			"-xd"
		]

		try:
			subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			print(f"{svg_path} generated.")
		except subprocess.CalledProcessError as e:
			raise RuntimeError(f"Fail to generate svg from sdf file: {e}") from e

		return str(svg_path)

