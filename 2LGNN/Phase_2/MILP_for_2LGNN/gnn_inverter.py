import torch
from torch_geometric.loader import DataLoader

from src import chemicalgraph, network, dataset

import pandas as pd

from read_instance_2layer_2LMM_GNN import *

import ann_inverter

# import pulp
import pulp_modified as pulp

def read_GNN_config(config_file):
    params = dict()
    with open(config_file, 'r') as f:
        line = f.readline()
        while line:
            line_all = [l.rstrip(',\n') for l in line.split(' ')]
            param_name = line_all[0]
            if 'ANN' in param_name:
                arr = tuple([int(n) for n in line_all[1:]])
            else:
                arr = int(line_all[1])
            params[param_name] = arr
            line = f.readline()

    if 'extra_x_length' not in params:
        params['extra_x_length'] = None

    return params

def cRT_to_RT(psi):
    RT = chemicalgraph.RootedTree()
    RT.CID = psi.index
    for (i, v) in enumerate(psi.vertex):
        atom = chemicalgraph.Atom()
        atom.index = i
        atom.atom_type = v[0].rstrip('1234567890')
        atom.atom_name = v[0]
        atom.status = True
        atom.chg = psi.chg[i]

        RT.vertex_list[i] = atom

    n = len(psi.vertex)
    for i in range(n):
        for j in range(i + 1, n):
            if psi.beta[i][j] != 0:
                RT.parent[j] = i
                v1, v2, mul = i, j, psi.beta[i][j]
                RT.vertex_list[v1].add_neighbour(v2, mul, RT.vertex_list[v2].atom_type)
                RT.vertex_list[v2].add_neighbour(v1, mul, RT.vertex_list[v1].atom_type) 

    RT.root = 0

    return RT

class GNN_inverter:
    def __init__(self, model_path=None):
        self.device = torch.device('cpu')

        self.model_path = model_path
        self.GNN = None
        self.ft_embedding = None
        

    def _build_GNN(self, params):
        self.GNN = network.GNN_2L_no_edge_type(
                in_channels_FT=params['in_channels_FT'],
                hidden_channels_FT=params['hidden_channels_FT'],
                num_layers_FT=params['num_layers_FT'],
                p_channels_FT=params['p_channels_FT'],
                ANN_channels_FT=params['ANN_channels_FT'],
                in_channels_IN=params['in_channels_IN'],
                hidden_channels_IN=params['hidden_channels_IN'],
                num_layers_IN=params['num_layers_IN'],
                p_channels_IN=params['p_channels_IN'],
                ANN_channels_IN=params['ANN_channels_IN'],
                extra_x_length=params['extra_x_length'],
                device=self.device
            )

        self.learned_weights = torch.load(self.model_path, map_location=self.device)
        self.GNN.load_state_dict(self.learned_weights)

        self.p_cv = params['p_channels_IN']
        self.num_layers = params['num_layers_IN']
        self.num_node_features = params['in_channels_IN'] + params['ANN_channels_FT'][-1] 
        self.hidden_layer_size = params['hidden_channels_IN']
        self.ann_layers = [self.p_cv] + params['ANN_channels_IN']

    def build_ft_fataset(self, GNN, set_F):
        fringe_tree_list = list()
        ft_index_map = dict()
        for i, psi in enumerate(set_F):
            rt = cRT_to_RT(psi)
            fringe_tree_list.append(rt)
            ft_index_map[psi.index] = i
        ft_dataset = dataset.FringeTree_Datasets(root='./data')
        ft_dataset.construct_from_FringeTree_list(fringe_tree_list)

        ft_dataloader = DataLoader(ft_dataset, batch_size=len(ft_dataset))
        for ft_dataset_batch in ft_dataloader:
            break

        ft_dataset_batch.to(self.device)
        GNN.eval()
        self.ft_embedding = GNN.GNN_FT(ft_dataset_batch).detach().numpy()
        self.ft_embedding = {psi.index: self.ft_embedding[ft_index_map[psi.index], :] for psi in set_F}

    def _get_weights(self, learned_weights):
        self.weight_w = dict()
        self.weight_bias = dict()
        self.weight_w_cv = dict()

        p_cv = self.p_cv
        l_max = self.num_layers
        z_max = self.num_node_features
        h_size = self.hidden_layer_size

        for l in range(0, l_max):
            for z_p in range(1, (z_max if l == 0 else h_size) + 1):
                for z in range(1, h_size + 1):
                    self.weight_w[(0, l, z_p, z)] = learned_weights[f'GNN_2L_int.convs.{l}.lin_l_self.weight'][z-1, z_p-1].item()
                    self.weight_w[(1, l, z_p, z)] = learned_weights[f'GNN_2L_int.convs.{l}.lin_l.weight'][z-1, z_p-1].item()
            for z in range(1, h_size + 1):
                self.weight_bias[(0, z, l)] = learned_weights[f'GNN_2L_int.convs.{l}.bias'][z-1].item()
        for p in range(1, p_cv + 1):
            for z in range(1, h_size + 1):
                self.weight_w_cv[(z, p)] = learned_weights['GNN_2L_int.convs_to_p.weight'][p-1, z-1].item()

        self.ann_weight_tensor = list()
        self.ann_bias_matrix = list()

        self.ann_bias_matrix.append([0] * self.ann_layers[0])

        for i in range(len(self.ann_layers) - 1):
            tmp_array = learned_weights[f'GNN_2L_int.Linears.{i}.weight'].numpy().T
            self.ann_weight_tensor.append(np.array(tmp_array))
            tmp_array = learned_weights[f'GNN_2L_int.Linears.{i}.bias'].numpy()
            self.ann_bias_matrix.append(np.array(tmp_array))

    def prepare_variables(self,
        t_C,
        t_C_tilde,
        t_T,
        t_F,
        k_C,
        k_C_tilde,
        m_C
    ):
        p_cv = self.p_cv
        l_max = self.num_layers
        z_max = self.num_node_features
        h_size = self.hidden_layer_size

        # Define theta_cv
        theta_cv_tmp = {p: pulp.LpVariable(f"theta_cv_tmp({p})") for p in range(1, p_cv + 1)}
        delta_cv_tmp = {p: pulp.LpVariable(f"delta_cv_tmp({p})", 0, 1, cat=pulp.LpBinary) for p in range(1, p_cv + 1)}
        theta_cv = {p: pulp.LpVariable(f"theta_cv({p})") for p in range(1, p_cv + 1)}

        # Define theta_C
        theta_C_tmp = {(i, z, ll): pulp.LpVariable(f"theta_C_tmp({i},{z},{ll})") for i in range(1, t_C + 1)
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        delta_C_tmp = {(i, z, ll): pulp.LpVariable(f"delta_C_tmp({i},{z},{ll})", 0, 1, cat=pulp.LpBinary) for i in range(1, t_C + 1)
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        theta_C = {(i, z, ll): pulp.LpVariable(f"theta_C({i},{z},{ll})") for i in range(1, t_C + 1)
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        
        # Define theta_T
        theta_T_tmp = {(i, z, ll): pulp.LpVariable(f"theta_T_tmp({i},{z},{ll})") for i in range(0, t_T + 2)  # for convenience
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        delta_T_tmp = {(i, z, ll): pulp.LpVariable(f"delta_T_tmp({i},{z},{ll})", 0, 1, cat=pulp.LpBinary) for i in range(0, t_T + 2)  # for convenience
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        theta_T = {(i, z, ll): pulp.LpVariable(f"theta_T({i},{z},{ll})") for i in range(0, t_T + 2)  # for convenience
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        
        # Define theta_F
        theta_F_tmp = {(i, z, ll): pulp.LpVariable(f"theta_F_tmp({i},{z},{ll})") for i in range(0, t_F + 2)  # for convenience
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        delta_F_tmp = {(i, z, ll): pulp.LpVariable(f"delta_F_tmp({i},{z},{ll})", 0, 1, cat=pulp.LpBinary) for i in range(0, t_F + 2)  # for convenience
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}
        theta_F = {(i, z, ll): pulp.LpVariable(f"theta_F({i},{z},{ll})") for i in range(0, t_F + 2)  # for convenience
                    for ll in range(0, l_max + 1) for z in range(1, (z_max if ll == 0 else h_size) + 1)}

        # Define theta_C_minus, theta_C_plus
        theta_C_minus = {(k, z, ll): pulp.LpVariable(f"theta_C_minus({k},{z},{ll})") for k in range(k_C_tilde + 1, m_C + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}
        theta_C_plus = {(k, z, ll): pulp.LpVariable(f"theta_C_plus({k},{z},{ll})") for k in range(k_C_tilde + 1, m_C + 1) # instead of [1, k_C], it should be [k_C_tilde + 1, m_C]
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}

        # Define theta_T_minus, theta_T_plus
        theta_T_minus = {(i, z, ll): pulp.LpVariable(f"theta_T_minus({i},{z},{ll})") for i in range(1, t_T + 2)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}
        theta_T_plus = {(i, z, ll): pulp.LpVariable(f"theta_T_plus({i},{z},{ll})") for i in range(1, t_T + 2)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}

        # Define theta_F_minus, theta_F_plus
        theta_F_minus = {(i, z, ll): pulp.LpVariable(f"theta_F_minus({i},{z},{ll})") for i in range(1, t_F + 2)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}
        theta_F_plus = {(i, z, ll): pulp.LpVariable(f"theta_F_plus({i},{z},{ll})") for i in range(1, t_F + 2)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}
        
        # Define theta_CT_T, theta_TC_T
        theta_CT_T = {(k, z, ll): pulp.LpVariable(f"theta_CT_T({k},{z},{ll})") for k in range(1, k_C + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}
        theta_TC_T = {(k, z, ll): pulp.LpVariable(f"theta_TC_T({k},{z},{ll})") for k in range(1, k_C + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}

        # Define theta_CT_C, theta_TC_C, theta_TF_F
        theta_CT_C = {(i, z, ll): pulp.LpVariable(f"theta_CT_C({i},{z},{ll})") for i in range(1, t_T + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)} # why R+
        theta_TC_C = {(i, z, ll): pulp.LpVariable(f"theta_TC_C({i},{z},{ll})") for i in range(1, t_T + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}
        theta_TF_F = {(i, z, ll): pulp.LpVariable(f"theta_TF_F({i},{z},{ll})") for i in range(1, t_T + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}

        # Define theta_CF_F
        theta_CF_F = {(c, z, ll): pulp.LpVariable(f"theta_CF_F({c},{z},{ll})") for c in range(1, t_C_tilde + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}

        # Define theta_CF_C, theta_TF_T
        theta_CF_C = {(i, z, ll): pulp.LpVariable(f"theta_CF_C({i},{z},{ll})") for i in range(1, t_F + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}
        theta_TF_T = {(i, z, ll): pulp.LpVariable(f"theta_TF_T({i},{z},{ll})") for i in range(1, t_F + 1)
                    for ll in range(0, l_max) for z in range(1, h_size + 1)}

        return theta_cv, theta_cv_tmp, delta_cv_tmp, theta_C, theta_T, theta_F, theta_C_tmp, theta_T_tmp, theta_F_tmp, delta_C_tmp, delta_T_tmp, delta_F_tmp, \
            theta_C_minus, theta_C_plus, theta_T_minus, theta_T_plus, theta_F_minus, theta_F_plus, \
            theta_CT_T, theta_TC_T, theta_CT_C, theta_TC_C, theta_TF_F, theta_CF_F, theta_CF_C, theta_TF_T

    def add_constraints_graph_convolution(
        self,
        # Model
        MILP,
        # Constants
        t_C,
        t_C_tilde,
        t_T,
        t_F,
        k_C,
        k_C_tilde,
        m_C,
        val,
        Code_F,
        Lambda_int,
        I_a_minus,
        I_a_plus,
        I_b_minus,
        I_b_plus,
        F_C,
        F_T,
        F_F,
        head_C,
        tail_C,
        # Binary Variables
        delta_alpha_C,
        delta_alpha_T,
        delta_alpha_F,
        delta_fr_C,
        delta_fr_T,
        delta_fr_F,
        delta_beta_C,
        delta_beta_T,
        delta_beta_F,
        delta_beta_plus,
        delta_beta_minus,
        delta_beta_in,
        e_C,
        e_T,
        e_F,
        delta_chi_T,
        delta_chi_F,
        chi_T,
        chi_F,
        v_T,
        v_F,
        # Integer Variables
        alpha_C,
        alpha_T,
        alpha_F,
        deg_C,
        deg_T,
        deg_F,
        eledeg_C,
        eledeg_T,
        eledeg_F,
        hyddeg_C,
        hyddeg_T,
        hyddeg_F,
        # General Variables
        theta_cv, 
        theta_cv_tmp,
        delta_cv_tmp,
        theta_C, 
        theta_T, 
        theta_F,
        theta_C_tmp, 
        theta_T_tmp, 
        theta_F_tmp,
        delta_C_tmp,
        delta_T_tmp,
        delta_F_tmp,
        theta_C_minus, 
        theta_C_plus, 
        theta_T_minus, 
        theta_T_plus, 
        theta_F_minus, 
        theta_F_plus,
        theta_CT_T, 
        theta_TC_T, 
        theta_CT_C, 
        theta_TC_C, 
        theta_TF_F, 
        theta_CF_F, 
        theta_CF_C, 
        theta_TF_T
    ):
        p_cv = self.p_cv
        l_max = self.num_layers
        z_max = self.num_node_features  
        h_size = self.hidden_layer_size
        M = {l: 10000 for l in range(0, l_max + 1)}   
        kappa = 0.1

        ft_embedding = self.ft_embedding

        w, bias, w_cv = self.weight_w, self.weight_bias, self.weight_w_cv

        # -------- Constraint (69) --------
        for i in range(1, t_C + 1):
            MILP += theta_C[(i, 1, 0)] == delta_alpha_C[(i, "C4")], f"MILP-2LM-M-GC-(65)-C-1-{i}"
            MILP += theta_C[(i, 2, 0)] == delta_alpha_C[(i, "O2")], f"MILP-2LM-M-GC-(65)-C-2-{i}"
            MILP += theta_C[(i, 3, 0)] == delta_alpha_C[(i, "N3")], f"MILP-2LM-M-GC-(65)-C-3-{i}"
            MILP += theta_C[(i, 4, 0)] == deg_C[i] + hyddeg_C[i], f"MILP-2LM-M-GC-(65)-C-4-{i}"
            MILP += theta_C[(i, 5, 0)] == pulp.lpSum(val[a] * delta_alpha_C[(i, a)] for a in Lambda_int) - eledeg_C[i], f"MILP-2LM-M-GC-(65)-C-5-{i}"
            MILP += theta_C[(i, 6, 0)] == hyddeg_C[i], f"MILP-2LM-M-GC-(65)-C-6-{i}"
            MILP += theta_C[(i, 7, 0)] == eledeg_C[i], f"MILP-2LM-M-GC-(65)-C-7-{i}"

            for j in range(8, z_max + 1):
                MILP += theta_C[(i, j, 0)] == pulp.lpSum(ft_embedding[Code_F[psi]][j - 8] * delta_fr_C[(i, Code_F[psi])] for psi in F_C[i]), f"MILP-2LM-M-GC-(65)-C-{j}-{i}"

            # for j in range(1, z_max + 1):
            #     MILP += theta_C[(i, j, 0)] == 1, f"test-C-{i}-{j}"

        for i in range(1, t_T + 1):
            MILP += theta_T[(i, 1, 0)] == delta_alpha_T[(i, "C4")], f"MILP-2LM-M-GC-(65)-T-1-{i}"
            MILP += theta_T[(i, 2, 0)] == delta_alpha_T[(i, "O2")], f"MILP-2LM-M-GC-(65)-T-2-{i}"
            MILP += theta_T[(i, 3, 0)] == delta_alpha_T[(i, "N3")], f"MILP-2LM-M-GC-(65)-T-3-{i}"
            MILP += theta_T[(i, 4, 0)] == deg_T[i] + hyddeg_T[i], f"MILP-2LM-M-GC-(65)-T-4-{i}"
            MILP += theta_T[(i, 5, 0)] == pulp.lpSum(val[a] * delta_alpha_T[(i, a)] for a in Lambda_int) - eledeg_T[i], f"MILP-2LM-M-GC-(65)-T-5-{i}"
            MILP += theta_T[(i, 6, 0)] == hyddeg_T[i], f"MILP-2LM-M-GC-(65)-T-6-{i}"
            MILP += theta_T[(i, 7, 0)] == eledeg_T[i], f"MILP-2LM-M-GC-(65)-T-7-{i}"

            for j in range(8, z_max + 1):
                MILP += theta_T[(i, j, 0)] == pulp.lpSum(ft_embedding[Code_F[psi]][j - 8] * delta_fr_T[(i, Code_F[psi])] for psi in F_T[i]), f"MILP-2LM-M-GC-(65)-T-{j}-{i}"

            # for j in range(1, z_max + 1):
            #     MILP += theta_T[(i, j, 0)] == v_T[i], f"test-T-{i}-{j}"

        for i in range(1, t_F + 1):
            MILP += theta_F[(i, 1, 0)] == delta_alpha_F[(i, "C4")], f"MILP-2LM-M-GC-(65)-F-1-{i}"
            MILP += theta_F[(i, 2, 0)] == delta_alpha_F[(i, "O2")], f"MILP-2LM-M-GC-(65)-F-2-{i}"
            MILP += theta_F[(i, 3, 0)] == delta_alpha_F[(i, "N3")], f"MILP-2LM-M-GC-(65)-F-3-{i}"
            MILP += theta_F[(i, 4, 0)] == deg_F[i] + hyddeg_F[i], f"MILP-2LM-M-GC-(65)-F-4-{i}"
            MILP += theta_F[(i, 5, 0)] == pulp.lpSum(val[a] * delta_alpha_F[(i, a)] for a in Lambda_int) - eledeg_F[i], f"MILP-2LM-M-GC-(65)-F-5-{i}"
            MILP += theta_F[(i, 6, 0)] == hyddeg_F[i], f"MILP-2LM-M-GC-(65)-F-6-{i}"
            MILP += theta_F[(i, 7, 0)] == eledeg_F[i], f"MILP-2LM-M-GC-(65)-F-7-{i}"

            for j in range(8, z_max + 1):
                MILP += theta_F[(i, j, 0)] == pulp.lpSum(ft_embedding[Code_F[psi]][j - 8] * delta_fr_F[(i, Code_F[psi])] for psi in F_F[i]), f"MILP-2LM-M-GC-(65)-F-{j}-{i}"
    
            # for j in range(1, z_max + 1):
            #     MILP += theta_F[(i, j, 0)] == v_F[i], f"test-F-{i}-{j}"

        # -------- Constraint (70) --------
        for z in range(1, h_size + 1):  # missing
            for l in range(0, (l_max - 1) + 1):  # missing
                for i in range(1, t_C_tilde + 1):
                    MILP += pulp.lpSum(w[(0, l, z_p, z)] * theta_C[(i, z_p, l)]
                                for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                            pulp.lpSum(theta_C_minus[(k, z, l)] for k in I_a_minus[i]) + \
                            pulp.lpSum(theta_C_plus[(k, z, l)] for k in I_a_plus[i]) + \
                            pulp.lpSum(theta_CT_T[(k, z, l)] for k in I_b_plus[i]) + \
                            pulp.lpSum(theta_TC_T[(k, z, l)] for k in I_b_minus[i]) + \
                            theta_CF_F[(i, z, l)] + bias[(0, z, l)] == theta_C_tmp[(i, z, l + 1)], f"MILP-2LM-M-GC-(66)-{z}-{l}-{i}-1"

                    MILP += theta_C[(i, z, l + 1)] - kappa * theta_C_tmp[(i, z, l + 1)] >= -M[l + 1] * delta_C_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-C-1-{z}-{l}-{i}-1"
                    MILP += theta_C[(i, z, l + 1)] - kappa * theta_C_tmp[(i, z, l + 1)] <= M[l + 1] * delta_C_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-C-1-{z}-{l}-{i}-2"
                    MILP += theta_C[(i, z, l + 1)] - theta_C_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_C_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-C-2-{z}-{l}-{i}-1"
                    MILP += theta_C[(i, z, l + 1)] - theta_C_tmp[(i, z, l + 1)] <= M[l + 1] * (1 - delta_C_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-C-2-{z}-{l}-{i}-2"
                    MILP += theta_C_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_C_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-C-3-{z}-{l}-{i}-1"
                    MILP += theta_C_tmp[(i, z, l + 1)] <= M[l + 1] * delta_C_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-C-3-{z}-{l}-{i}-2"

                    MILP += theta_C[(i, z, l + 1)] <= M[l + 1], f"MILP-2LM-M-GC-(66)-{z}-{l}-{i}-2"
                    MILP += theta_C[(i, z, l + 1)] >= -M[l + 1], f"MILP-2LM-M-GC-(66)-{z}-{l}-{i}-3" # by Zhu

        # -------- Constraint (71) --------
        for z in range(1, h_size + 1):
            for l in range(0, (l_max - 1) + 1):
                for i in range(t_C_tilde + 1, t_C + 1):
                    MILP += pulp.lpSum(w[(0, l, z_p, z)] * theta_C[(i, z_p, l)]
                                for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                            pulp.lpSum(theta_C_minus[(k, z, l)] for k in I_a_minus[i]) + \
                            pulp.lpSum(theta_C_plus[(k, z, l)] for k in I_a_plus[i]) + \
                            pulp.lpSum(theta_CT_T[(k, z, l)] for k in I_b_plus[i]) + \
                            pulp.lpSum(theta_TC_T[(k, z, l)] for k in I_b_minus[i]) + \
                            bias[(0, z, l)] == theta_C_tmp[(i, z, l + 1)], f"MILP-2LM-M-GC-(67)-{z}-{l}-{i}-1"

                    MILP += theta_C[(i, z, l + 1)] - kappa * theta_C_tmp[(i, z, l + 1)] >= -M[l + 1] * delta_C_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-C-1-{z}-{l}-{i}-1"
                    MILP += theta_C[(i, z, l + 1)] - kappa * theta_C_tmp[(i, z, l + 1)] <= M[l + 1] * delta_C_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-C-1-{z}-{l}-{i}-2"
                    MILP += theta_C[(i, z, l + 1)] - theta_C_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_C_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-C-2-{z}-{l}-{i}-1"
                    MILP += theta_C[(i, z, l + 1)] - theta_C_tmp[(i, z, l + 1)] <= M[l + 1] * (1 - delta_C_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-C-2-{z}-{l}-{i}-2"
                    MILP += theta_C_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_C_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-C-3-{z}-{l}-{i}-1"
                    MILP += theta_C_tmp[(i, z, l + 1)] <= M[l + 1] * delta_C_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-C-3-{z}-{l}-{i}-2"

                    MILP += theta_C[(i, z, l + 1)] <= M[l + 1], f"MILP-2LM-M-GC-(67)-{z}-{l}-{i}-2" # by Zhu 
                    MILP += theta_C[(i, z, l + 1)] >= -M[l + 1], f"MILP-2LM-M-GC-(67)-{z}-{l}-{i}-3" # by Zhu 

        # tmp_T = {(i, z, l): pulp.LpVariable(f"tmp_T({i},{z},{l})") for i in range(1, t_T + 1) for z in range(1, h_size + 1) for l in range(0, l_max)}
        # tmp_F = {(i, z, l): pulp.LpVariable(f"tmp_F({i},{z},{l})") for i in range(1, t_F + 1) for z in range(1, h_size + 1) for l in range(0, l_max)}

        # -------- Constraint (72) --------
        for z in range(1, h_size + 1):
            for l in range(0, (l_max - 1) + 1):
                for i in range(1, t_T + 1):
                    MILP += theta_T_tmp[(i, z, l + 1)] - (pulp.lpSum(w[(0, l, z_p, z)] * theta_T[(i, z_p, l)]
                                for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                            theta_T_minus[(i, z, l)] + \
                            theta_T_plus[(i, z, l)] + \
                            theta_CT_C[(i, z, l)] + \
                            theta_TC_C[(i, z, l)] + \
                            theta_TF_F[(i, z, l)] + \
                            bias[(0, z, l)]) <= M[l + 1] * (1 - v_T[i]), f"MILP-2LM-M-GC-(68)-{z}-{l}-{i}-1-1"
                    MILP += theta_T_tmp[(i, z, l + 1)] - (pulp.lpSum(w[(0, l, z_p, z)] * theta_T[(i, z_p, l)]
                                for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                            theta_T_minus[(i, z, l)] + \
                            theta_T_plus[(i, z, l)] + \
                            theta_CT_C[(i, z, l)] + \
                            theta_TC_C[(i, z, l)] + \
                            theta_TF_F[(i, z, l)] + \
                            bias[(0, z, l)]) >= -M[l + 1] * (1 - v_T[i]), f"MILP-2LM-M-GC-(68)-{z}-{l}-{i}-1-2"

                    MILP += theta_T[(i, z, l + 1)] - kappa * theta_T_tmp[(i, z, l + 1)] >= -M[l + 1] * delta_T_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-T-1-{z}-{l}-{i}-1"
                    MILP += theta_T[(i, z, l + 1)] - kappa * theta_T_tmp[(i, z, l + 1)] <= M[l + 1] * delta_T_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-T-1-{z}-{l}-{i}-2"
                    MILP += theta_T[(i, z, l + 1)] - theta_T_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_T_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-T-2-{z}-{l}-{i}-1"
                    MILP += theta_T[(i, z, l + 1)] - theta_T_tmp[(i, z, l + 1)] <= M[l + 1] * (1 - delta_T_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-T-2-{z}-{l}-{i}-2"
                    MILP += theta_T_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_T_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-T-3-{z}-{l}-{i}-1"
                    MILP += theta_T_tmp[(i, z, l + 1)] <= M[l + 1] * delta_T_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-T-3-{z}-{l}-{i}-2"

                    MILP += theta_T[(i, z, l + 1)] <= M[l + 1] * v_T[i] , f"MILP-2LM-M-GC-(68)-{z}-{l}-{i}-2" # by Zhu
                    MILP += theta_T[(i, z, l + 1)] >= -M[l + 1] * v_T[i], f"MILP-2LM-M-GC-(68)-{z}-{l}-{i}-3" # by Zhu
                if (0, z, l) in theta_T:
                    MILP += theta_T[(0, z, l)] == 0, f"MILP-2LM-M-GC-(68)-{z}-{l}-3"
                if (t_T + 1, z, l) in theta_T:
                    MILP += theta_T[(t_T + 1, z, l)] == 0, f"MILP-2LM-M-GC-(68)-{z}-{l}-4"
        # for i in range(1, t_T + 1):
        #     for z in range(1, h_size + 1):
        #         for l in range(0, l_max):
        #             MILP += tmp_T[(i, z, l)] <= M[l + 1] * (1 - v_T[i]), f"MILP-2LM-M-GC-(68)-{z}-{l}-{i}-5" # by Zhu
        #             MILP += tmp_T[(i, z, l)] >= -M[l + 1] * (1 - v_T[i]), f"MILP-2LM-M-GC-(68)-{z}-{l}-{i}-6" # by Zhu

        # -------- Constraint (73) --------
        for z in range(1, h_size + 1):
            for l in range(0, (l_max - 1) + 1):
                for i in range(1, t_F + 1):
                    MILP += theta_F_tmp[(i, z, l + 1)] - (pulp.lpSum(w[(0, l, z_p, z)] * theta_F[(i, z_p, l)]
                                for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                            theta_F_minus[(i, z, l)] + \
                            theta_F_plus[(i, z, l)] + \
                            theta_CF_C[(i, z, l)] + \
                            theta_TF_T[(i, z, l)] + \
                            bias[(0, z, l)]) <= M[l + 1] * (1 - v_F[i]), f"MILP-2LM-M-GC-(69)-{z}-{l}-{i}-1-1"
                    MILP += theta_F_tmp[(i, z, l + 1)] - (pulp.lpSum(w[(0, l, z_p, z)] * theta_F[(i, z_p, l)]
                                for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                            theta_F_minus[(i, z, l)] + \
                            theta_F_plus[(i, z, l)] + \
                            theta_CF_C[(i, z, l)] + \
                            theta_TF_T[(i, z, l)] + \
                            bias[(0, z, l)]) >= -M[l + 1] * (1 - v_F[i]), f"MILP-2LM-M-GC-(69)-{z}-{l}-{i}-1-2"

                    MILP += theta_F[(i, z, l + 1)] - kappa * theta_F_tmp[(i, z, l + 1)] >= -M[l + 1] * delta_F_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-F-1-{z}-{l}-{i}-1"
                    MILP += theta_F[(i, z, l + 1)] - kappa * theta_F_tmp[(i, z, l + 1)] <= M[l + 1] * delta_F_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-F-1-{z}-{l}-{i}-2"
                    MILP += theta_F[(i, z, l + 1)] - theta_F_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_F_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-F-2-{z}-{l}-{i}-1"
                    MILP += theta_F[(i, z, l + 1)] - theta_F_tmp[(i, z, l + 1)] <= M[l + 1] * (1 - delta_F_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-F-2-{z}-{l}-{i}-2"
                    MILP += theta_F_tmp[(i, z, l + 1)] >= -M[l + 1] * (1 - delta_F_tmp[(i, z, l + 1)]), f"MILP-GC-LeakyReLU-F-3-{z}-{l}-{i}-1"
                    MILP += theta_F_tmp[(i, z, l + 1)] <= M[l + 1] * delta_F_tmp[(i, z, l + 1)], f"MILP-GC-LeakyReLU-C-F-{z}-{l}-{i}-2"

                    MILP += theta_F[(i, z, l + 1)] <= M[l + 1] * v_F[i], f"MILP-2LM-M-GC-(69)-{z}-{l}-{i}-2"
                    MILP += theta_F[(i, z, l + 1)] >= -M[l + 1] * v_F[i], f"MILP-2LM-M-GC-(69)-{z}-{l}-{i}-3" # by Zhu
                if (0, z, l) in theta_F:
                    MILP += theta_F[(0, z, l)] == 0, f"MILP-2LM-M-GC-(69)-{z}-{l}-3"
                if (t_F + 1, z, l) in theta_F:
                    MILP += theta_F[(t_F + 1, z, l)] == 0, f"MILP-2LM-M-GC-(69)-{z}-{l}-4"
        # for i in range(1, t_F + 1):
        #     for z in range(1, h_size + 1):
        #         for l in range(0, l_max):
        #             MILP += tmp_F[(i, z, l)] <= M[l + 1] * (1 - v_F[i]), f"MILP-2LM-M-GC-(69)-{z}-{l}-{i}-5" # by Zhu
        #             MILP += tmp_F[(i, z, l)] >= -M[l + 1] * (1 - v_F[i]), f"MILP-2LM-M-GC-(69)-{z}-{l}-{i}-6" # by Zhu

        # -------- Constraint (74) --------
        for z in range(1, h_size + 1):  
            for l in range(0, l_max):
                for k in range(k_C_tilde + 1, m_C + 1): # range modified
                    MILP += theta_C_minus[(k, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(tail_C[k], z_p, l)]
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) <= \
                                M[l] * (1 - e_C[k]), f"MILP-2LM-M-GC-(70)-1-{z}-{l}-{k}-1"
                    MILP += theta_C_minus[(k, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(tail_C[k], z_p, l)]
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) >= \
                                -M[l] * (1 - e_C[k]), f"MILP-2LM-M-GC-(70)-1-{z}-{l}-{k}-2"
                    MILP += theta_C_plus[(k, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(head_C[k], z_p, l)]
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) <= \
                                M[l] * (1 - e_C[k]), f"MILP-2LM-M-GC-(70)-2-{z}-{l}-{k}-1"
                    MILP += theta_C_plus[(k, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(head_C[k], z_p, l)]
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) >= \
                                -M[l] * (1 - e_C[k]), f"MILP-2LM-M-GC-(70)-2-{z}-{l}-{k}-2"
        for k in range(k_C_tilde + 1, m_C + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_C_minus[(k, z, l)] <= M[l] * e_C[k], f"MILP-2LM-M-GC-(70)-3-{k}-{z}-{l}-1"
                    MILP += theta_C_plus[(k, z, l)] <= M[l] * e_C[k], f"MILP-2LM-M-GC-(70)-3-{k}-{z}-{l}-2"
                    MILP += theta_C_minus[(k, z, l)] >= -M[l] * e_C[k], f"MILP-2LM-M-GC-(70)-4-{k}-{z}-{l}-1" # by Zhu
                    MILP += theta_C_plus[(k, z, l)] >= -M[l] * e_C[k], f"MILP-2LM-M-GC-(70)-4-{k}-{z}-{l}-2" # by Zhu

        # -------- Constraint (75) --------
        for z in range(1, h_size + 1):
            for l in range(0, l_max):
                for i in range(2, t_T + 1):
                    MILP += theta_T_minus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i - 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) <= \
                                M[l] * (1 - e_T[i]), f"MILP-2LM-M-GC-(71)-1-{z}-{l}-{i}-1"
                    MILP += theta_T_minus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i - 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) >= \
                                -M[l] * (1 - e_T[i]), f"MILP-2LM-M-GC-(71)-1-{z}-{l}-{i}-2"
                for i in range(1, t_T):
                    MILP += theta_T_plus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i + 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) <= \
                                M[l] * (1 - e_T[i + 1]), f"MILP-2LM-M-GC-(71)-2-{z}-{l}-{i}-1"
                    MILP += theta_T_plus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i + 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) >= \
                                -M[l] * (1 - e_T[i + 1]), f"MILP-2LM-M-GC-(71)-2-{z}-{l}-{i}-2"
        for i in range(1, t_T + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_T_minus[(i, z, l)] <= M[l] * e_T[i], f"MILP-2LM-M-GC-(71)-3-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_T_plus[(i, z, l)] <= M[l] * e_T[i + 1], f"MILP-2LM-M-GC-(71)-3-{i}-{z}-{l}-2" # by Zhu
                    MILP += theta_T_minus[(i, z, l)] >= -M[l] * e_T[i], f"MILP-2LM-M-GC-(71)-4-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_T_plus[(i, z, l)] >= -M[l] * e_T[i + 1], f"MILP-2LM-M-GC-(71)-4-{i}-{z}-{l}-2" # by Zhu              

        # -------- Constraint (76) --------
        for z in range(1, h_size + 1):
            for l in range(0, l_max):
                for i in range(2, t_F + 1):
                    MILP += theta_F_minus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(i - 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) <= \
                                M[l] * (1 - e_F[i]), f"MILP-2LM-M-GC-(72)-1-{z}-{l}-{i}-1"
                    MILP += theta_F_minus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(i - 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) >= \
                                -M[l] * (1 - e_F[i]), f"MILP-2LM-M-GC-(72)-1-{z}-{l}-{i}-2"
                for i in range(1, t_F):
                    MILP += theta_F_plus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(i + 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) <= \
                                M[l] * (1 - e_F[i + 1]), f"MILP-2LM-M-GC-(72)-2-{z}-{l}-{i}-1"
                    MILP += theta_F_plus[(i, z, l)] - pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(i + 1, z_p, l)] 
                                    for z_p in range(1, (z_max if l == 0 else h_size) + 1)) >= \
                                -M[l] * (1 - e_F[i + 1]), f"MILP-2LM-M-GC-(72)-2-{z}-{l}-{i}-2"
        for i in range(1, t_F + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_F_minus[(i, z, l)] <= M[l] * e_F[i], f"MILP-2LM-M-GC-(72)-3-{i}-{z}-{l}-1"
                    MILP += theta_F_plus[(i, z, l)] <= M[l] * e_F[i + 1], f"MILP-2LM-M-GC-(72)-3-{i}-{z}-{l}-2" # by Zhu
                    MILP += theta_F_minus[(i, z, l)] >= -M[l] * e_F[i], f"MILP-2LM-M-GC-(72)-4-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_F_plus[(i, z, l)] >= -M[l] * e_F[i + 1], f"MILP-2LM-M-GC-(72)-4-{i}-{z}-{l}-2" # by Zhu

        # -------- Constraint (77) --------
        for k in range(1, k_C + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    for i in range(1, t_T + 1):
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i]) <= \
                                    theta_CT_T[(k, z, l)], f"MILP-2LM-M-GC-(73)-1-{k}-{z}-{l}-{i}-1"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i]) >= \
                                    theta_CT_T[(k, z, l)], f"MILP-2LM-M-GC-(73)-1-{k}-{z}-{l}-{i}-2"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i + 1]) <= \
                                    theta_TC_T[(k, z, l)], f"MILP-2LM-M-GC-(73)-2-{k}-{z}-{l}-{i}-1"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(i, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i + 1]) >= \
                                    theta_TC_T[(k, z, l)], f"MILP-2LM-M-GC-(73)-2-{k}-{z}-{l}-{i}-2"
        for k in range(1, k_C + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_CT_T[(k, z, l)] <= M[l] * delta_chi_T[k], f"MILP-2LM-M-GC-(73)-3-{k}-{z}-{l}-1"
                    MILP += theta_TC_T[(k, z, l)] <= M[l] * delta_chi_T[k], f"MILP-2LM-M-GC-(73)-3-{k}-{z}-{l}-2"
                    MILP += theta_CT_T[(k, z, l)] >= -M[l] * delta_chi_T[k], f"MILP-2LM-M-GC-(73)-4-{k}-{z}-{l}-1" # by Zhu
                    MILP += theta_TC_T[(k, z, l)] >= -M[l] * delta_chi_T[k], f"MILP-2LM-M-GC-(73)-4-{k}-{z}-{l}-2" # by Zhu

        # -------- Constraint (78) --------
        for i in range(1, t_T + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    for k in range(1, k_C + 1):
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(tail_C[k], z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i]) <= \
                                    theta_CT_C[(i, z, l)], f"MILP-2LM-M-GC-(74)-1-{i}-{z}-{l}-{k}-1"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(tail_C[k], z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i]) >= \
                                    theta_CT_C[(i, z, l)], f"MILP-2LM-M-GC-(74)-1-{i}-{z}-{l}-{k}-2"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(head_C[k], z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i + 1]) <= \
                                    theta_TC_C[(i, z, l)], f"MILP-2LM-M-GC-(74)-2-{i}-{z}-{l}-{k}-1"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(head_C[k], z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_T[(i, k)] + e_T[i + 1]) >= \
                                    theta_TC_C[(i, z, l)], f"MILP-2LM-M-GC-(74)-2-{i}-{z}-{l}-{k}-2"
        for i in range(1, t_T + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_CT_C[(i, z, l)] <= M[l] * (1 - e_T[i]), f"MILP-2LM-M-GC-(74)-3-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_TC_C[(i, z, l)] <= M[l] * (1 - e_T[i + 1]), f"MILP-2LM-M-GC-(74)-3-{i}-{z}-{l}-2" # by Zhu
                    MILP += theta_CT_C[(i, z, l)] >= -M[l] * (1 - e_T[i]), f"MILP-2LM-M-GC-(74)-5-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_TC_C[(i, z, l)] >= -M[l] * (1 - e_T[i + 1]), f"MILP-2LM-M-GC-(74)-5-{i}-{z}-{l}-2" # by Zhu

                    MILP += theta_CT_C[(i, z, l)] <= M[l] * v_T[i], f"MILP-2LM-M-GC-(74)-4-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_TC_C[(i, z, l)] <= M[l] * v_T[i], f"MILP-2LM-M-GC-(74)-4-{i}-{z}-{l}-2" # by Zhu
                    MILP += theta_CT_C[(i, z, l)] >= -M[l] * v_T[i], f"MILP-2LM-M-GC-(74)-6-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_TC_C[(i, z, l)] >= -M[l] * v_T[i], f"MILP-2LM-M-GC-(74)-6-{i}-{z}-{l}-2" # by Zhu
                    
        # -------- Constraint (79) --------
        for c in range(1, t_C_tilde + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    for i in range(1, t_F + 1):
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(i, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_F[(i, c)] + e_F[i]) <= \
                                    theta_CF_F[(c, z, l)], f"MILP-2LM-M-GC-(75)-1-{c}-{z}-{l}-{i}-1"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(i, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_F[(i, c)] + e_F[i]) >= \
                                    theta_CF_F[(c, z, l)], f"MILP-2LM-M-GC-(75)-1-{c}-{z}-{l}-{i}-2"
        for c in range(1, t_C_tilde + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_CF_F[(c, z, l)] <= M[l] * delta_chi_F[c], f"MILP-2LM-M-GC-(75)-2-{c}-{z}-{l}-1"
                    MILP += theta_CF_F[(c, z, l)] >= -M[l] * delta_chi_F[c], f"MILP-2LM-M-GC-(75)-2-{c}-{z}-{l}-2" # by Zhu

        # -------- Constraint (80) --------
        for i in range(1, t_F + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    for c in range(1, t_C_tilde + 1):
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(c, z_p, l)] # should be c
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_F[(i, c)] + e_F[i]) <= \
                                    theta_CF_C[(i, z, l)], f"MILP-2LM-M-GC-(76)-1-{i}-{z}-{l}-{c}-1"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_C[(c, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_F[(i, c)] + e_F[i]) >= \
                                    theta_CF_C[(i, z, l)], f"MILP-2LM-M-GC-(76)-1-{i}-{z}-{l}-{c}-2"
        for i in range(1, t_F + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_CF_C[(i, z, l)] <= M[l] * pulp.lpSum(chi_F[(i, c)] 
                                for c in range(1, t_C_tilde + 1)), f"MILP-2LM-M-GC-(76)-2-{i}-{z}-{l}-1"
                    MILP += theta_CF_C[(i, z, l)] >= -M[l] * pulp.lpSum(chi_F[(i, c)] 
                                for c in range(1, t_C_tilde + 1)), f"MILP-2LM-M-GC-(76)-2-{i}-{z}-{l}-2" # by Zhu
                    MILP += theta_CF_C[(i, z, l)] <= M[l] * (1 - e_F[i]), f"MILP-2LM-M-GC-(76)-3-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_CF_C[(i, z, l)] >= -M[l] * (1 - e_F[i]), f"MILP-2LM-M-GC-(76)-3-{i}-{z}-{l}-2" # by Zhu

        # -------- Constraint (81) --------
        for i in range(1, t_T + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    for j in range(1, t_F + 1):
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(j, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_F[(j, t_C_tilde + i)] + e_F[j]) <= \
                                    theta_TF_F[(i, z, l)], f"MILP-2LM-M-GC-(77)-1-{i}-{z}-{l}-{j}-1"
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_F[(j, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_F[(j, t_C_tilde + i)] + e_F[j]) >= \
                                    theta_TF_F[(i, z, l)], f"MILP-2LM-M-GC-(77)-1-{i}-{z}-{l}-{j}-2"
        for i in range(1, t_T + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_TF_F[(i, z, l)] <= M[l] * delta_chi_F[t_C_tilde + i], f"MILP-2LM-M-GC-(77)-2-{i}-{z}-{l}-1"
                    MILP += theta_TF_F[(i, z, l)] >= -M[l] * delta_chi_F[t_C_tilde + i], f"MILP-2LM-M-GC-(77)-2-{i}-{z}-{l}-2" # by Zhu

        # -------- Constraint (82) --------
        for i in range(1, t_F + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    for j in range(1, t_T + 1):
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(j, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) - \
                                M[l] * (1 - chi_F[(i, t_C_tilde + j)] + e_F[i]) <= \
                                    theta_TF_T[(i, z, l)], f"MILP-2LM-M-GC-(78)-1-{i}-{z}-{l}-{j}-1" # by Zhu
                        MILP += pulp.lpSum(w[(1, l, z_p, z)] * theta_T[(j, z_p, l)]
                                        for z_p in range(1, (z_max if l == 0 else h_size) + 1)) + \
                                M[l] * (1 - chi_F[(i, t_C_tilde + j)] + e_F[i]) >= \
                                    theta_TF_T[(i, z, l)], f"MILP-2LM-M-GC-(78)-1-{i}-{z}-{l}-{j}-2" # by Zhu
        for i in range(1, t_F + 1):
            for z in range(1, h_size + 1):
                for l in range(0, l_max):
                    MILP += theta_TF_T[(i, z, l)] <= M[l] * pulp.lpSum(chi_F[(i, t_C_tilde + j)]
                                for j in range(1, t_T + 1)), f"MILP-2LM-M-GC-(78)-2-{i}-{z}-{l}-1"
                    MILP += theta_TF_T[(i, z, l)] >= -M[l] * pulp.lpSum(chi_F[(i, t_C_tilde + j)]
                                for j in range(1, t_T + 1)), f"MILP-2LM-M-GC-(78)-2-{i}-{z}-{l}-2" # by Zhu
                    MILP += theta_TF_T[(i, z, l)] <= M[l] * (1 - e_F[i]), f"MILP-2LM-M-GC-(78)-3-{i}-{z}-{l}-1" # by Zhu
                    MILP += theta_TF_T[(i, z, l)] >= -M[l] * (1 - e_F[i]), f"MILP-2LM-M-GC-(78)-3-{i}-{z}-{l}-2" # by Zhu

        # -------- Constraint (83) --------
        for p in range(1, p_cv + 1):
            MILP += pulp.lpSum(w_cv[(z, p)] * theta_C[(i, z, l_max)] 
                        for i in range(1, t_C + 1) for z in range(1, h_size + 1)) + \
                    pulp.lpSum(w_cv[(z, p)] * theta_T[(i, z, l_max)] 
                        for i in range(1, t_T + 1) for z in range(1, h_size + 1)) + \
                    pulp.lpSum(w_cv[(z, p)] * theta_F[(i, z, l_max)] 
                        for i in range(1, t_F + 1) for z in range(1, h_size + 1)) == \
                    theta_cv_tmp[p], f"MILP-2LM-M-GC-(79)-{p}"

            MILP += theta_cv[p] - kappa * theta_cv_tmp[p] >= -M[l_max] * delta_cv_tmp[p], f"MILP-GC-LeakyReLU-p-1-{p}-1"
            MILP += theta_cv[p] - kappa * theta_cv_tmp[p] <= M[l_max] * delta_cv_tmp[p], f"MILP-GC-LeakyReLU-p-1-{p}-2"
            MILP += theta_cv[p] - theta_cv_tmp[p] >= -M[l_max] * (1 - delta_cv_tmp[p]), f"MILP-GC-LeakyReLU-p-2-{p}-1"
            MILP += theta_cv[p] - theta_cv_tmp[p] <= M[l_max] * (1 - delta_cv_tmp[p]), f"MILP-GC-LeakyReLU-p-2-{p}-2"
            MILP += theta_cv_tmp[p] >= -M[l_max] * (1 - delta_cv_tmp[p]), f"MILP-GC-LeakyReLU-p-3-{p}-1"
            MILP += theta_cv_tmp[p] <= M[l_max] * delta_cv_tmp[p], f"MILP-GC-LeakyReLU-p-3-{p}-2"
                    
        return MILP

    def build_MILP(self, model_path, params, MILP, GNN, learned_weights, set_F,
            t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, m_C, val, Code_F, Lambda_int, I_a_minus, I_a_plus, I_b_minus, I_b_plus, F_C, F_T, F_F, head_C, tail_C,
            delta_alpha_C, delta_alpha_T, delta_alpha_F, delta_fr_C, delta_fr_T, delta_fr_F, delta_beta_C, delta_beta_T, delta_beta_F, 
            delta_beta_plus, delta_beta_minus, delta_beta_in, e_C, e_T, e_F, delta_chi_T, delta_chi_F, chi_T, chi_F, v_T, v_F,
            alpha_C, alpha_T, alpha_F, deg_C, deg_T, deg_F, eledeg_C, eledeg_T, eledeg_F, hyddeg_C, hyddeg_T, hyddeg_F,
            target_value_lb, target_value_ub
        ):
        self.model_path = model_path

        self.p_cv = params['p_channels_IN']
        self.num_layers = params['num_layers_IN']
        self.num_node_features = params['in_channels_IN'] + params['ANN_channels_FT'][-1] 
        self.hidden_layer_size = params['hidden_channels_IN']
        self.ann_layers = [self.p_cv] + list(params['ANN_channels_IN'])

        # self._build_GNN(params)
        self.build_ft_fataset(GNN, set_F)
        self._get_weights(learned_weights)

        theta_cv, theta_cv_tmp, delta_cv_tmp, theta_C, theta_T, theta_F, theta_C_tmp, theta_T_tmp, theta_F_tmp, delta_C_tmp, delta_T_tmp, delta_F_tmp, \
            theta_C_minus, theta_C_plus, theta_T_minus, theta_T_plus, theta_F_minus, theta_F_plus, \
            theta_CT_T, theta_TC_T, theta_CT_C, theta_TC_C, theta_TF_F, theta_CF_F, theta_CF_C, theta_TF_T = \
                self.prepare_variables(t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, m_C)

        MILP = self.add_constraints_graph_convolution(
                MILP,
                t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, m_C, val, Code_F, Lambda_int, I_a_minus, I_a_plus, I_b_minus, I_b_plus, F_C, F_T, F_F, head_C, tail_C,
                delta_alpha_C, delta_alpha_T, delta_alpha_F, delta_fr_C, delta_fr_T, delta_fr_F, delta_beta_C, delta_beta_T, delta_beta_F, 
                delta_beta_plus, delta_beta_minus, delta_beta_in, e_C, e_T, e_F, delta_chi_T, delta_chi_F, chi_T, chi_F, v_T, v_F,
                alpha_C, alpha_T, alpha_F, deg_C, deg_T, deg_F, eledeg_C, eledeg_T, eledeg_F, hyddeg_C, hyddeg_T, hyddeg_F,
        
                theta_cv, theta_cv_tmp, delta_cv_tmp, theta_C, theta_T, theta_F, theta_C_tmp, theta_T_tmp, theta_F_tmp,
                delta_C_tmp, delta_T_tmp, delta_F_tmp, theta_C_minus, theta_C_plus, theta_T_minus, theta_T_plus, theta_F_minus, theta_F_plus,
                theta_CT_T, theta_TC_T, theta_CT_C, theta_TC_C, theta_TF_F, theta_CF_F, theta_CF_C, theta_TF_T
            )

        ann = ann_inverter.ANN(self.ann_weight_tensor, self.ann_bias_matrix)
        milp_ann_constants = ann_inverter.initialize_constants(ann)
        ann_a, ann_b, ann_b_hat, ann_c, ann_z = milp_ann_constants
        milp_ann_vars = ann_inverter.initialize_lp_variables(ann, ann_a, ann_b)
        MILP = ann_inverter.build_MILP_ReLU(
            MILP,
            ann,
            milp_ann_vars,
            milp_ann_constants,
            target_value_lb,
            target_value_ub,
            # eps,
            list()
        )

        _, y, _ = milp_ann_vars       

        for i in range(1, self.p_cv + 1):
            MILP += y[(1, i)] == theta_cv[i], f"ann-connecting-{i}"

        return MILP, y[ann.output_layer[0]]

    def inspection(self, GNN, test_sdf, device):
        cg_list = chemicalgraph.construct_from_SDF_file(test_sdf, use_rdkit=False)
        cg_list.get_fringe_tree()
        cg_list.set_value("1", 0.0)
        cg_dataset = dataset.ChemicalGraph_Datasets(root='./data/')
        cg_dataset.construct_from_ChemicalGraph_list(cg_list, FT_use=True)
        ft_dataset = dataset.FringeTree_Datasets(root='./data')
        ft_dataset.construct_from_FringeTree_list(cg_list.fringe_tree_list)

        loader = DataLoader(cg_dataset, batch_size=len(cg_dataset))
        ft_dataloader = DataLoader(ft_dataset, batch_size=len(ft_dataset))
        for ft_dataset_batch in ft_dataloader:
            break

        GNN.eval()
        for data in loader:
            data = data.to(device)
            ft_dataset_batch.to(device)
            output = GNN(data, ft_dataset_batch)
            y_pred = output.cpu().detach().numpy()



        return y_pred[0][0]





