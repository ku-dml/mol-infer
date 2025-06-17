###
# define the class ChemicalGraph
###

import sys, copy, math
import numpy as np
import pandas as pd
from collections import deque, defaultdict
from operator import itemgetter, attrgetter

from rdkit import Chem

from src import constants

class Atom:
    def __init__(self):
        self.vertex_name = ""
        self.atom_name = ""
        self.atom_type = ""
        self.index = None
        self.status = False
        self.charge = 0
        self.rad = 0
        self.label = None

        self.adj = list()
        self.adj_atom_type = dict()
        self.beta = dict()
        self.fringe_tree_TS = None


    def __str__(self):
        return self.atom_name

    def copy_from(self, u):
        self.vertex_name = u.vertex_name
        self.atom_name = u.atom_name
        self.atom_type = u.atom_type
        self.index = u.index
        self.status = u.status
        self.charge = u.charge
        self.rad = u.rad
        self.label = u.label

    def copy(self):
        #### need to use this when copying!!!
        repulica = copy.deepcopy(self)

        return repulica        

    @property
    def degree(self):
        return len(self.adj)

    @property
    def valence(self):
        return sum([mul for mul in self.beta.values()])

    @property
    def hyd(self):
        return sum([1 for a in self.adj_atom_type.values() if a == 'H'])

    def _label(self):
        self.label = self.valence * constants.CONST._atom_num_upper + constants.CONST.atomic_number[self.atom_type]
        return self.label

    def change_atom_name(self):
        val = self.valence
        self.atom_name += f"{val}"
        self._label()

    def remove_index_in_adj(self, index_to_remove):
        if index_to_remove in self.adj:
            self.adj.remove(index_to_remove)
            self.beta.pop(index_to_remove, None)

    def add_neighbour(self, v, mul, v_atom_type):
        self.adj.append(v)
        self.beta[v] = mul
        self.adj_atom_type[v] = v_atom_type

    def check_atom_type(self, atom_type_name):
        return 1 if self.atom_type == atom_type_name else 0

    def check_adj_bond(self, bond_type):
        return 1 if bond_type in self.beta.values() else 0

    def get_feature(self):
        feature = list()
        #####

        # feature.append(constants.CONST.atomic_number[self.atom_type])
        # feature.append(constants.CONST.atom_mass[self.atom_type] / 10)
        # feature.append(constants.CONST.atom_electronegativity[self.atom_type])
        # feature.append(constants.CONST.atom_ion_energy[self.atom_type])
        # feature.append(constants.CONST.atom_electron_affinity[self.atom_type])

        feature.append(self.check_atom_type('C'))
        feature.append(self.check_atom_type('O'))
        feature.append(self.check_atom_type('N'))

        # feature.append(constants.CONST.atom_mass[self.atom_type] / (10 * constants.CONST.atomic_number[self.atom_type]))
        
        feature.append(self.degree)
        feature.append(self.valence)
        feature.append(self.hyd)
        feature.append(self.charge)

        return np.array(feature)

    def get_feature_for_FT(self, root=None):
        ### only for the part of GNN of embedding FT, possibly different from the above one

        feature = list()
        #####
        if root is not None:
            feature.append(1 if self.index == root else 0)
        feature.append(self.check_atom_type('C'))
        feature.append(self.check_atom_type('O'))
        feature.append(self.check_atom_type('N'))
        feature.append(self.check_atom_type('H'))
        feature.append(self.check_adj_bond(2))
        feature.append(self.check_adj_bond(3))
        

        feature.append(constants.CONST.atomic_number[self.atom_type])
        # feature.append(constants.CONST.atom_mass[self.atom_type] / 10)
        # feature.append(constants.CONST.atom_electronegativity[self.atom_type])
        feature.append(constants.CONST.atom_ion_energy[self.atom_type])
        feature.append(constants.CONST.atom_electron_affinity[self.atom_type])
        
        feature.append(self.degree)
        feature.append(self.valence)
        feature.append(self.hyd)
        feature.append(self.charge)

        # m = len(feature)
        # for i in range(m):
        #     for j in range(i, m):
        #         feature.append(feature[i] * feature[j])

        return np.array(feature)

    def get_feature_FT(self, TS_list):
        feature = list()

        if self.fringe_tree_TS is not None:
            # # feature.append(self.fringe_tree_TS.atom_num)
            # # feature.append(self.fringe_tree_TS.height)
            # # feature.append(self.fringe_tree_TS.num_of_H)
            # # feature.append(self.fringe_tree_TS.average_mass / 100)
            # # feature.append(self.fringe_tree_TS.sum_positive)
            # # feature.append(self.fringe_tree_TS.sum_negative)
            # feature.append(self.fringe_tree_TS.num_double_bond)
            # feature.append(self.fringe_tree_TS.num_triple_bond)
            # feature.append(self.fringe_tree_TS.ZagrebIndex2(-0.5))
            # feature.append(self.fringe_tree_TS.WienerIndex())
            # feature.append(self.fringe_tree_TS.ABCIndex())
            # feature.append(self.fringe_tree_TS.ZagrebIndex1(-0.5)) 

            # feature = self.fringe_tree_TS.get_FT_feature()  
            feature = TS_list[self.fringe_tree_TS].get_FT_feature()

        return np.array(feature)     

class ChemicalGraph:
    def __init__(self, cg=None):
        if cg is None:
            self.CID = ""
            self.SDF = list()
            self.vertex_list = dict()
        else:
            self.CID = cg.CID
            self.SDF = cg.SDF[:]
            self.vertex_list = copy.deepcopy(cg.vertex_list)
        self.clear_tmp_variables()

    def __getitem__(self, u):
        if u in self.vertex_list.keys():
            return self.vertex_list[u]
        else:
            return None

    def __len__(self):
        return self.atom_num

    def clear_tmp_variables(self):
        ### temporary variables in order to shorten running time
        self._num_double_bond = None
        self._num_triple_bond = None
        self._core_indices = None
        self._atom_num = None
        self._edge_num = None
        self._average_mass = None
        self._interior_indices = None

    def copy(self):
        #### need to use this when copying!!!
        repulica = copy.deepcopy(self)
        repulica.clear_tmp_variables()

        return repulica

    def construct_from_SDF_file(self, f, store_SDF=False):
        flag = True
        ans = 1
        line_all = []
        while True:
            line = f.readline()
            if not line:
                ans = 0
                break
            line_all.append(line)

            if flag:
                self.CID = line.rstrip('\n')
                line = f.readline()
                line_all.append(line)
                line = f.readline()
                line_all.append(line)
                line = f.readline()
                line_all.append(line)
                if not line:
                    ans = 0
                    return 0
                n, m = int(line[0:3]), int(line[3:6])
                for vertex in range(n):
                    atom = Atom()
                    line = f.readline()
                    line_all.append(line)
                    line_split = line.split()

                    atom.index = vertex
                    atom.atom_type = line_split[3]
                    atom.atom_name = line_split[3]
                    atom.status = True

                    self.vertex_list[vertex] = atom

                for i in range(m):
                    line = f.readline()
                    line_all.append(line)
                    v1, v2, mul = int(line[0:3]), int(line[3:6]), int(line[6:9])
                    self.vertex_list[v1 - 1].add_neighbour(v2 - 1, mul, self.vertex_list[v2 - 1].atom_type)
                    self.vertex_list[v2 - 1].add_neighbour(v1 - 1, mul, self.vertex_list[v1 - 1].atom_type) 

                flag = False
            else:
                if line[0:4] == "$$$$":
                    break
                # if line[0:6] == "M  RAG" and no_radical:
                #     ans = -2
                if line[0:6] == "M  CHG":
                    # if no_ion:
                    #     ans = -1
                    # else:
                    line_split = line.split()
                    n_chg = int(line_split[2])
                    for i in range(n_chg):
                        v = int(line_split[2 * i + 3])
                        chg = int(line_split[2 * i + 4])
                        self.vertex_list[v - 1].charge = chg

        if store_SDF:
            self.SDF = line_all[:]

        return ans

    def construct_from_mol(self, mol, no_ion=True, no_radical=True, store_SDF=False):
        self.CID = mol.GetProp('_Name')
        n, m = mol.GetNumAtoms(), mol.GetNumBonds()

        mol_idx = dict()

        for vertex, a in enumerate(mol.GetAtoms()):
            atom = Atom()
            atom.index = vertex

            mol_idx[a.GetIdx()] = vertex
            atom.atom_type = a.GetSymbol()
            atom.atom_name = a.GetSymbol()
            atom.status = True

            atom.charge = a.GetFormalCharge()
            atom.rad = a.GetNumRadicalElectrons()

            self.vertex_list[vertex] = atom

        for b in mol.GetBonds():
            v1, v2, mul = b.GetBeginAtom().GetIdx(), b.GetEndAtom().GetIdx(), b.GetBondTypeAsDouble()
            atom1, atom2 = self.vertex_list[mol_idx[v1]], self.vertex_list[mol_idx[v2]]
            atom1.add_neighbour(mol_idx[v2], mul, atom2.atom_type)
            atom2.add_neighbour(mol_idx[v1], mul, atom1.atom_type)


        if store_SDF:
            self.SDF = copy.deepcopy(mol)

    # def output_SDF(self, filename, add=True):
    #     r = 'a' if add else 'w'
    #     with open(filename, r) as f:
    #         for line in self.SDF:
    #             f.write(line)

    def output_SDF(self, filename, add=True):
        r = 'a' if add else 'w'
        f = open(filename, r)
        writer = Chem.SDWriter(f)
        writer.write(self.SDF)
    
    def change_atom_name(self):
        for vertex in self.vertex_list.values():
            vertex.change_atom_name()

    @property
    def atom_num(self):
        if self._atom_num is None:
            self._atom_num = len(self.vertex_list)

        return self._atom_num

    @property
    def edge_num(self):
        if self._edge_num is None:
            num = 0
            for vertex in self.vertex_list.values():
                num += vertex.degree

            self._edge_num = int(num / 2)

        return self._edge_num

    @property
    def rank(self):
        return self.edge_num - self.atom_num + 1

    def degree(self, u):
        return self.vertex_list[u].degree

    def remove_vertex(self, vertex_index):
        ## a sub method used by other methods
        for vertex in self.vertex_list.values():
            vertex.remove_index_in_adj(vertex_index)

        self.vertex_list.pop(vertex_index, None)

    def remove_atom_type(self, atom_type):
        vertex_index_to_remove = [vertex.index for vertex in self.vertex_list.values() 
                    if vertex.atom_type == atom_type]
        
        for index in vertex_index_to_remove:
            self.remove_vertex(index)

        self.clear_tmp_variables()

    def remove_leaf(self):
        vertex_index_to_remove = [vertex.index for vertex in self.vertex_list.values() 
                    if vertex.degree == 1]
        
        for index in vertex_index_to_remove:
            self.remove_vertex(index)

        self.clear_tmp_variables()

        return len(vertex_index_to_remove)

    @property
    def core_indices(self):
        if self._core_indices is None:
            cg_copy = self.copy()

            while len(cg_copy) != 2 and cg_copy.remove_leaf() != 0:
                pass

            self._core_indices = list(cg_copy.vertex_list.keys())

        return self._core_indices

    def interior_indices(self, k=2):
        if self._interior_indices is None:
            cg_copy = self.copy()

            cg_copy.remove_atom_type('H')

            for _ in range(k):
                if len(cg_copy) == 2:
                    break
                cg_copy.remove_leaf()

            self._interior_indices = list(cg_copy.vertex_list.keys())

        return self._interior_indices

    @property
    def num_double_bond(self):
        if self._num_double_bond is None:
            num = 0
            for vertex in self.vertex_list.values():
                num += len([1 for mul in vertex.beta.values() if mul == 2])

            self._num_double_bond = int(num / 2)

        return self._num_double_bond

    @property
    def num_triple_bond(self):
        if self._num_triple_bond is None:
            num = 0
            for vertex in self.vertex_list.values():
                num += len([1 for mul in vertex.beta.values() if mul == 3])

            self._num_triple_bond = int(num / 2)

        return self._num_triple_bond

    @property
    def average_mass(self):
        if self._average_mass is None:
            num = sum([constants.CONST.atom_mass[vertex.atom_type] for vertex in self.vertex_list.values()])
            num /= self.atom_num()

            self._average_mass = num

        return self._average_mass

    def add_vertex(self, v, u, mul, index_for_u=None):
        #### v is the index to attach, u is a Atom object
        if index_for_u is None:
            index_for_u = max(self.vertex_list.keys()) + 1
            u.index = index_for_u

        self.vertex_list[index_for_u] = u
        self.vertex_list[v].add_neighbour(index_for_u, mul, u.atom_type)
        self.vertex_list[index_for_u].add_neighbour(v, mul, self.vertex_list[v].atom_type)

        return index_for_u

class TreeSignature():
    def __init__(self):
        self.delta = list()
        self.mul = list()
        self.charge = list()

    def __eq__(self, other):
        return (self.delta == other.delta and self.mul == other.mul and self.charge == other.charge)

    def __str__(self):
        delta_str = [constants.CONST.get_chemical_symbol(d) if d > constants.CONST._atom_num_upper 
                        else d for d in self.delta]
        return f"({delta_str}, {self.mul}, {self.charge})"

    def __len__(self):
        return self.atom_num()

    @property
    def atom_num(self):
        return len(self.charge)

    @property
    def height(self):
        return max([d for d in self.delta if d <= constants.CONST._atom_num_upper])

    @property
    def num_of_H(self):
        return self.delta.count(constants.CONST._atom_num_upper + 1)

    @property
    def average_mass(self):
        mass = 0
        for d in self.delta:
            if d <= constants.CONST._atom_num_upper:
                continue
            val = math.floor(d / constants.CONST._atom_num_upper)
            atomic_number = d - val * constants.CONST._atom_num_upper
            atom = [a for a in constants.CONST.atomic_number.keys() if constants.CONST.atomic_number[a] == atomic_number]
            mass += constants.CONST.atom_mass[atom[0]]

        mass /= self.atom_num

        return mass

    @property
    def sum_positive(self):
        return sum([c for c in self.charge if c > 0])

    @property
    def sum_negative(self):
        return sum([c for c in self.charge if c < 0])

    @property
    def num_double_bond(self):
        return self.mul.count(2)

    @property
    def num_triple_bond(self):
        return self.mul.count(3)

    def ZagrebIndex1(self, lmd=-0.5):
        zi1 = 0

        deg, _ = self._get_degree_vector()
        for d in deg:
            if d != 0:
                zi1 += math.pow(d, lmd)

        return zi1

    def ZagrebIndex2(self, lmd=-0.5):
        zi2 = 0

        deg, edge_set = self._get_degree_vector()
        for (u, v) in edge_set:
            zi2 += math.pow(deg[u] * deg[v], lmd)

        return zi2

    def ABCIndex(self):
        ABCI = 0

        deg, edge_set = self._get_degree_vector()

        for (u, v) in edge_set:
            ABCI += math.pow((deg[u] + deg[v] - 2) / (deg[u] * deg[v]), 0.5)

        return ABCI

    def WienerIndex(self, normalized=True):
        _, edge_set = self._get_degree_vector()
        N = self.atom_num
        D = np.ones((N, N)) * (N + 1)
        for (u, v) in edge_set:
            D[u, v] = 1
        for i in range(N):
            D[i, i] = 0

        for k in range(N):
            for i in range(N):
                for j in range(N):
                    if D[i, k] + D[k, j] < D[i, j]:
                        D[i, j] = D[i, k] + D[k, j]

        WI = 0.0
        for i in range(N):
            for j in range(i + 1, N):
                WI += D[i, j]

        if normalized:
            WI /= math.pow(N, 2)

        return WI


    def _get_degree_vector(self):
        # This implementation used the fact that the height is at most 3
        deg = [0 for i in range(self.atom_num)]

        ind_parent = [0, 0, 0]
        edge_set = list()
        for i in range(len(self.mul)):
            if self.delta[2 * i + 3] == 1:
                deg[ind_parent[0]] += 1
                deg[i + 1] += 1
                ind_parent[1] = i + 1
                edge_set.append((0, i + 1))
            elif self.delta[2 * i + 3] == 2:
                deg[ind_parent[1]] += 1
                deg[i + 1] = 1
                ind_parent[2] = i + 1
                edge_set.append((ind_parent[1], i + 1))
            elif self.delta[2 * i + 3] == 3:
                deg[ind_parent[2]] += 1
                deg[i + 1] += 1
                edge_set.append((ind_parent[2], i + 1))

        return deg, edge_set

    def _get_depth_vector(self):
        depth = [self.delta[i] for i in range(1, len(self.delta), 2)]
        return depth

    def extend(self, other):
        self.delta.extend(other.delta)
        self.mul.extend(other.mul)
        self.charge.extend(other.charge)

    def get_FT_feature(self):
        # define default fv for a fringe tree:¥
        feature = []
        # feature.append(self.fringe_tree_TS.atom_num)
        # feature.append(self.fringe_tree_TS.height)
        # feature.append(self.fringe_tree_TS.num_of_H)
        # feature.append(self.fringe_tree_TS.average_mass / 100)
        # feature.append(self.fringe_tree_TS.sum_positive)
        # feature.append(self.fringe_tree_TS.sum_negative)
        # feature.append(self.num_double_bond)
        # feature.append(self.num_triple_bond)
        # feature.append(self.ZagrebIndex2(-0.5))
        # feature.append(self.WienerIndex())
        # feature.append(self.ABCIndex())
        # feature.append(self.ZagrebIndex1(-0.5))

        return feature

class RootedTree(ChemicalGraph):
    def __init__(self, cg=None, root=None):
        super().__init__(cg)
        self.root = root
        self.parent = dict()

    def get_Fringe_tree(self, cg, root, k=2):
        # core_indices = cg.core_indices
        interior_indices = cg.interior_indices(k)
        self.root = root
        self.vertex_list[root] = Atom()
        self.vertex_list[root].copy_from(cg.vertex_list[root])

        visited = set()

        u = root
        self.parent[u] = None

        queue = deque([u])
        while len(queue) > 0:
            u = queue.popleft()
            visited.add(u)
            for v in cg.vertex_list[u].adj:
                if v in interior_indices or v in visited:
                    continue
                new_v = Atom()
                new_v.copy_from(cg.vertex_list[v])

                self.add_vertex(u, new_v, cg.vertex_list[u].beta[v], index_for_u=v)
                self.parent[v] = u

                queue.append(v)

    def _get_TS(self, u, depth):
        child_num = self.vertex_list[u].degree - 1
        if u == self.root:
            child_num += 1

        if child_num <= 0:
            TS = TreeSignature()
            TS.delta = [self.vertex_list[u].label, depth]
            TS.charge = [self.vertex_list[u].charge]

            return TS

        tmp_TS_list = list()
        for v in self.vertex_list[u].adj:
            if v == self.parent[u]:
                continue

            tmp_TS = self._get_TS(v, depth + 1)
            tmp_TS.mul.insert(0, self.vertex_list[u].beta[v])
            tmp_TS_list.append(tmp_TS)


        tmp_TS_list.sort(key=attrgetter('delta', 'mul', 'charge'))

        TS = TreeSignature()
        TS.delta = [self.vertex_list[u].label, depth]
        TS.charge = [self.vertex_list[u].charge]

        for tmp_TS in tmp_TS_list:
            TS.extend(tmp_TS)

        return TS

    def get_TS(self):
        return self._get_TS(self.root, 0)

class ChemicalGraph_list():
    def __init__(self, extra_fv_used=False):
        self.cg_dict = dict()  # CID -> index
        self.property_dict = dict() # CID -> value

        if extra_fv_used:
            self.extra_fv_dict = dict() # CID -> fv
        else:
            self.extra_fv_dict = None 

        self.fringe_tree_list = list() # fringe tree list of the whole dataset
        self.TS_list = list()          # TS of fringe tree

    def __iter__(self):
        return iter(self.cg_dict.values())

    def __next__(self):
        return next(self.cg_dict.values())

    def __str__(self):
        return str(self.cg_dict.keys())

    def __getitem__(self, CID):
        return self.at(CID)      

    def __len__(self):
        return len(self.cg_dict)  

    def append(self, cg):
        self.cg_dict[cg.CID] = cg
        self.property_dict[cg.CID] = None

    def set_value(self, CID, value):
        if CID in self.property_dict.keys():
            self.property_dict[CID] = value

    def get_value(self, CID):
        return self.property_dict[CID]

    # def normalization(self, y_min=0, y_max=1):
    #     if y_min > y_max:
    #         print("y_min must be smaller than y_max!!")
    #     else:
    #         max_y = max(self.property_dict.values())
    #         min_y = min(self.property_dict.values())

    #         for CID in self.property_dict.keys():
    #             value = self.property_dict[CID]
    #             self.property_dict[CID] = (value - min_y) / (max_y - min_y) * (y_max - y_min) + y_min

    def standardization(self):
        values = np.array([val for val in self.property_dict.values()])
        y_mean = values.mean(axis=0)
        y_std = values.std(axis=0)

        # if abs(y_std) > 1e-6:
        for CID in self.property_dict.keys():
            value = self.property_dict[CID]
            self.property_dict[CID] = (value - y_mean) / y_std

    def y_one_hot(self, n_class=2):
        self.class_size = defaultdict(int)
        for CID in self.property_dict.keys():
            value = int(self.property_dict[CID])
            self.class_size[value] += 1
            self.property_dict[CID] = [1 if i == value else 0 for i in range(n_class)]

    def at(self, CID):
        if CID in self.cg_dict:
            return self.cg_dict[CID]
        else:
            return None

    def read_value_file(self, filename):
        with open(filename, 'r') as f:
            line = f.readline()

            line = f.readline()
            while line:
                line_split = line.split(',')
                self.set_value(line_split[0], [float(val) for val in line_split[1:]])

                line = f.readline()

    def read_value_file_csv(self, filename, val):
        fp = pd.read_csv(filename, index_col=0)
        fp = fp[val]

        tmp_filename = "./tmp.txt"
        fp.to_csv(tmp_filename, sep=',')

        with open(tmp_filename, 'r') as f:
            line = f.readline()

            line = f.readline()
            while line:
                line_split = line.split(',')
                self.set_value(line_split[0], [float(val) for val in line_split[1:]])

                line = f.readline()

    def read_extra_fv_file(self, filename):
        x = pd.read_csv(filename, index_col=0)
        CIDs = list(x.index)
        CIDs_str = list(map(str, CIDs))
        length = None

        for CID, CID_str in zip(CIDs, CIDs_str):
            if CID_str in self.cg_dict.keys():
                self.extra_fv_dict[CID_str] = np.array(x.loc[CID])
                if length is None:
                    length = len(np.array(x.loc[CID]))

        return length

    def get_fringe_tree(self, k=2, store=True):
        for cg in self:
            interior_indices = cg.interior_indices(k)

            for u in interior_indices:
                rt = RootedTree()
                rt.get_Fringe_tree(cg, u)

                rt_TS = rt.get_TS()

                try:
                    rt_index = self.TS_list.index(rt_TS)
                except:
                    self.fringe_tree_list.append(rt)
                    self.TS_list.append(rt_TS)
                    rt_index = len(self.TS_list) - 1

                if store:
                    # cg.vertex_list[u].fringe_tree_TS = copy.deepcopy(rt_TS)
                    cg.vertex_list[u].fringe_tree_TS = rt_index

        return self.fringe_tree_list

def construct_from_SDF_file(filename, extra_fv_used=False, store_SDF=False, use_rdkit=False, check_qm9=False):
    if not use_rdkit:
        if check_qm9:
            qm9_uncharacterized = pd.read_csv("./qm9_uncharacterized.csv", index_col=0)
            qm9_uncharacterized = [f"gdb_{idx}" for idx in list(qm9_uncharacterized.index)]

        with open(filename, 'r') as f:

            cg_list = ChemicalGraph_list(extra_fv_used)
            while True:
                g = ChemicalGraph()
                ans = g.construct_from_SDF_file(f, store_SDF=store_SDF)
                if ans == 0:
                    break
                if ans == 1:
                    g.change_atom_name()

                    if not check_qm9 or (check_qm9 and g.CID not in qm9_uncharacterized):
                        cg_list.append(g)

    else:
        ### use rdkit ###
        suppl = Chem.SDMolSupplier(filename)
        cg_list = ChemicalGraph_list(extra_fv_used)

        for mol in suppl:
            print(mol.GetProp('_Name'))
            try:
                mol = Chem.AddHs(mol)
            except:
                filename = "ERROR.txt"
                writer = Chem.SDWriter(filename)
                writer.write(mol)
                break
            g = ChemicalGraph()
            g.construct_from_mol(mol, store_SDF=store_SDF)
            g.change_atom_name()
            cg_list.append(g)

    return cg_list

def create_H_suppresed_chemical_graph(cg_list):
    cg_list_without_H = ChemicalGraph_list()

    for cg in cg_list:
        cg_tmp = cg.copy()
        cg_tmp.remove_atom_type('H')

        cg_list_without_H.append(cg_tmp)
        cg_list_without_H.set_value(cg_tmp.CID, cg_list.property_dict[cg_tmp.CID])

    return cg_list_without_H

# def read_value_file(cg_list, filename):
#     with open(filename, 'r') as f:
#         line = f.readline()

#         line = f.readline()
#         while line:
#             line_split = line.split(',')
#             cg_list.set_value(line_split[0], float(line_split[1]))

#             line = f.readline()


