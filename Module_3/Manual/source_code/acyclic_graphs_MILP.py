"""
Copyright (C) 2020
by Discrete Mathematics Lab, 
   Department of Applied Mathematics and Physics, 
   Graduate School of Informatics,
   Kyoto University

Licensed under the MIT license, see License.txt
"""

"""
acyclic_graphs_MILP.py

This file implements functions that given 
files that contain the weight and bias values of
a trained ANN, construct an MILP model for
solving the inverse problem to
get a vector of resources for constructing 
an acyclic graph
"""

import pulp
import time
from collections import namedtuple
import sys
import ann_inverter
import pandas as pd
import math

CG_element = namedtuple("CG_element", "symbol, valence, mass")

def prepare_CG_element_info():
    ''' A function to prepare information of chemical elements. '''

    set_Lambda = list()
    set_Lambda.append(CG_element("B", 3, 108))
    set_Lambda.append(CG_element("C", 4, 120))
    set_Lambda.append(CG_element("O", 2, 160))
    set_Lambda.append(CG_element("N", 3, 140))
    set_Lambda.append(CG_element("F", 1, 190))
    set_Lambda.append(CG_element("Si", 4, 280))
    set_Lambda.append(CG_element("P", 4, 310))
    set_Lambda.append(CG_element("S", 2, 320))
    set_Lambda.append(CG_element("Cl", 1, 355))
    set_Lambda.append(CG_element("V", 3, 510))
    set_Lambda.append(CG_element("Br", 1, 800))
    set_Lambda.append(CG_element("Cd", 2, 1124))
    set_Lambda.append(CG_element("I", 1, 1270))
    set_Lambda.append(CG_element("Hg", 2, 2006))
    set_Lambda.append(CG_element("Pb", 2, 2076))
    set_Lambda.append(CG_element("Al", 3, 269))

    return set_Lambda

def prepare_CG_from_fv(set_Lambda, ann_training_data_filename):
    fv = pd.read_csv(ann_training_data_filename, sep=",")

    set_Lambda_new = list()
    for atom in set_Lambda:
        symbol_in = atom.symbol + "_in"
        symbol_ex = atom.symbol + "_ex"
        if symbol_in in fv.columns or symbol_ex in fv.columns:
            set_Lambda_new.append(atom)

    return set_Lambda_new

def prepare_fv(ann_training_data_filename,
    n_star, dia_star, k_star, d_max, bl_star, bh_star, 
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    delta_alpha, delta_beta_tilde, delta_beta_hat, ce_in, ce_ex, ce, 
    Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex,
    delta_tau, delta_tau_hat, ac_in, ac_ex, ac, n_H, 
    bc, bc_in, bc_ex, delta_dc, delta_dc_hat 
    ):
    ''' A function to prepare some variables about the descriptors 
        from the feature vector '''

    descriptors = dict()      # get variable from fv file
    stringoutput = dict()     # string used for output
    fv = pd.read_csv(ann_training_data_filename, sep=",")
    num_fv = len(list(fv.columns))

    # some variables used for AD constraint
    ceset = set()
    acset = set()
    desc_index_v = set()
    desc_index_e = set()

    for i, fv_name in enumerate(list(fv.columns)):
        if fv_name == "CID":
            pass
        elif fv_name == "n":
            descriptors[i] = n_star
            stringoutput[i] = "               n : "
        elif fv_name == "M":
            descriptors[i] = Mass# / (10 * n_star)
            stringoutput[i] = "             M/n : "
        elif fv_name == "#double_bond_in":
            descriptors[i] = b_in[2]
            stringoutput[i] = "    #double_bond_in : "
            desc_index_e.add(i)
        elif fv_name == "#double_bond_ex":
            descriptors[i] = b_ex[2]
            stringoutput[i] = "    #double_bond_ex : "
            desc_index_e.add(i)
        elif fv_name == "#triple_bond_in":
            descriptors[i] = b_in[3]
            stringoutput[i] = "    #triple_bond_in : "
            desc_index_e.add(i)
        elif fv_name == "#triple_bond_ex":
            descriptors[i] = b_ex[3]
            stringoutput[i] = "    #triple_bond_ex : "
            desc_index_e.add(i)
        elif fv_name == "H":
            descriptors[i] = n_H
            stringoutput[i] = "               H : "
        elif fv_name[0] == "#" and fv_name[1] == "d" and fv_name[-1] == "n":
            j = int(fv_name[7])
            descriptors[i] = dg_in[j]
            desc_index_v.add(i)
            stringoutput[i] = "        #degree{}_in : ".format(j)
        elif fv_name[0] == "#" and fv_name[1] == "d" and fv_name[-1] == "x":
            j = int(fv_name[7])
            descriptors[i] = dg_ex[j]
            desc_index_v.add(i)
            stringoutput[i] = "        #degree{}_ex : ".format(j)
        elif fv_name[0] == "B" and fv_name[1] == "c" and fv_name[-1] == "n":
            j1 = int(fv_name[3])
            j2 = int(fv_name[4])
            j3 = int(fv_name[5])
            descriptors[i] = bc_in[(j1, j2, j3)]
            desc_index_e.add(i)
            stringoutput[i] = "          Bc_{}{}{}_in : ".format(j1, 
                j2, j3)
        elif fv_name[0] == "B" and fv_name[1] == "c" and fv_name[-1] == "x":
            j1 = int(fv_name[3])
            j2 = int(fv_name[4])
            j3 = int(fv_name[5])
            descriptors[i] = bc_ex[(j1, j2, j3)]
            desc_index_e.add(i)
            stringoutput[i] = "          Bc_{}{}{}_ex : ".format(j1, 
                j2, j3)
        elif fv_name == "Diameter":
            descriptors[i] = dia_star / n_star
            stringoutput[i] = "      Diameter/n : "
        elif fv_name[-1] == "t":
            descriptors[i] = bh_star
            stringoutput[i] = "   branch-height : "
        elif fv_name[-1] == "r":
            descriptors[i] = bl_star
            stringoutput[i] = "branch-leaf-number : "
        else:
            str_tmp = fv_name
            posnum = -1
            for k in range(len(str_tmp)):
                if str_tmp[k] in "0123456789":
                    posnum = k
                    break
            if posnum == -1:
                ce_tmp = str_tmp[:-3]
                if str_tmp[-1] == "n":
                    descriptors[i] = ce_in[ce_tmp]
                elif str_tmp[-1] == "x":
                    descriptors[i] = ce_ex[ce_tmp]
                stringoutput[i] = \
                    "               {} : ".format(str_tmp)
                ceset.add(ce_tmp)
                desc_index_v.add(i)
            else:
                str_atom1 = str_tmp[0:k]
                str_atom2 = str_tmp[k + 1:-3]
                mul = int(str_tmp[k])
                gg = (str_atom1, str_atom2, mul)
                if gg not in Gamma_less:
                    gg = (str_atom2, str_atom1, mul)
                if str_tmp[-1] == "n":
                    descriptors[i] = ac_in[gg]
                elif str_tmp[-1] == "x":
                    descriptors[i] = ac_ex[gg]
                stringoutput[i] = "             {} : ".format(str_tmp)
                acset.add(gg)
                desc_index_e.add(i)

    return descriptors, stringoutput, num_fv, ceset, acset, \
        desc_index_v, desc_index_e


def prepare_dataset_for_scheme_graph(
    n_star, d_max, dia_star, k_star, bl_star, bh_star): 
    ''' A function for preparing some variables of the scheme graph '''
    
    a = d_max
    b = d_max - 1
    c = bh_star

    root_num = 1

    if bh_star == 2:
        s_star = 10
        # print("s* = " + str(s_star) + '\n')
    
        E_B_left = {1, 4}
        E_B_right = {3, 9}
        V_B_left = {2, 5}
        V_B_right = {4, 10}
        V_B_bh_star = {5, 6, 7, 8, 9, 10}
        V_B_s = [set() for _ in range(s_star + 1)]
        V_B_s[5] = {2, 5}
        V_B_s[6] = {2, 6}
        V_B_s[7] = {3, 7}
        V_B_s[8] = {3, 8}
        V_B_s[9] = {4, 9}
        V_B_s[10] = {4, 10}
        E_B_s = [set() for _ in range(s_star + 1)]
        E_B_s[5] = {1, 4}
        E_B_s[6] = {1, 5}
        E_B_s[7] = {2, 6}
        E_B_s[8] = {2, 7}
        E_B_s[9] = {3, 8}
        E_B_s[10] = {3, 9}
        
        E_B = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        E1_plus = [set() for _ in range(s_star + 1)]
        E1_plus[1] = {1, 2, 3}
        E1_plus[2] = {4, 5}
        E1_plus[3] = {6, 7}
        E1_plus[4] = {8, 9}
        E1_plus[5] = set()
        E1_plus[6] = set()
        E1_plus[7] = set()
        E1_plus[8] = set()
        E1_plus[9] = set()
        E1_plus[10] = set()
        E1_minus = [set() for _ in range(s_star + 1)]
        #E1_minus[1] = set()
        E1_minus[2] = {1}
        E1_minus[3] = {2}
        E1_minus[4] = {3}
        E1_minus[5] = {4}
        E1_minus[6] = {5}
        E1_minus[7] = {6}
        E1_minus[8] = {7}
        E1_minus[9] = {8}
        E1_minus[10] = {9}
        tail = [0, 1, 1, 1, 2, 2, 3, 3, 4, 4]
        head = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        Cld_B = [set() for _ in range(s_star + 1)]
        Cld_B[1] = {2, 3, 4}
        Cld_B[2] = {5, 6}
        Cld_B[3] = {7, 8}
        Cld_B[4] = {9, 10}
        Cld_B[5] = set()
        Cld_B[6] = set()
        Cld_B[7] = set()
        Cld_B[8] = set()
        Cld_B[9] = set()
        Cld_B[10] = set()

    elif bh_star == 1:
        s_star = 4
        # print("s* = " + str(s_star) + '\n')
    
        E_B_left = {1}
        E_B_right = {3}
        V_B_left = {2}
        V_B_right = {4}
        V_B_bh_star = {2, 3, 4}
        V_B_s = [set() for _ in range(s_star + 1)]
        V_B_s[2] = {2}
        V_B_s[3] = {3}
        V_B_s[4] = {4}
        E_B_s = [set() for _ in range(s_star + 1)]
        E_B_s[2] = {1}
        E_B_s[3] = {2}
        E_B_s[4] = {3}
        
        E_B = [1, 2, 3]
        E1_plus = [set() for _ in range(s_star + 1)]
        E1_plus[1] = {1, 2, 3}
        E1_minus = [set() for _ in range(s_star + 1)]
        E1_minus[2] = {1}
        E1_minus[3] = {2}
        E1_minus[4] = {3}
        tail = [0, 1, 1, 1]
        head = [0, 2, 3, 4]

        Cld_B = [set() for _ in range(s_star + 1)]
        Cld_B[1] = {2, 3, 4}

    c_star = s_star - 1
    t_star = n_star - (bh_star - 1) - (k_star + 1) * bl_star

    ###########4.3 scheme graph##########
    n_S_tree = int(1 + (d_max - 1) * 
        ((d_max - 1) ** k_star - 1)/(d_max - 2)) 
    n_T_tree = int(1 + (d_max - 2) * 
        ((d_max - 1) ** k_star - 1)/(d_max - 2)) 

    if d_max == 3:
        Cld_S = [set() for _ in range(n_S_tree + 1)]
        Cld_S[1] = {2, 3}
        Cld_S[2] = {4, 5}
        Cld_S[3] = {6, 7}

        Cld_T = [set() for _ in range(n_T_tree + 1)]
        Cld_T[1] = {2}
        Cld_T[2] = {3, 4}

        prt_S = [0, 0, 1, 1, 2, 2, 3, 3]
        prt_T = [0, 0, 1, 2, 2]

        P_T_prc = [set() for _ in range(n_T_tree + 1)]
        P_T_prc[1] = {2}
        P_T_prc[2] = {3}
        P_T_prc[3] = {4}
        P_T_prc[4] = set()

        P_S_prc = [set() for _ in range(n_S_tree + 1)]
        P_S_prc[1] = {2}
        P_S_prc[2] = {3, 4}
        P_S_prc[3] = {6}
        P_S_prc[4] = {5}
        P_S_prc[5] = set()
        P_S_prc[6] = {7}
        P_S_prc[7] = set()

        Dsn_S = [set() for _ in range(k_star + 1)]
        Dsn_S[2] = {4, 5, 6, 7}
        Dsn_S[1] = {2, 3}
        Dsn_S[0] = {1}
    else:   # dmax = 4
        Cld_S = [set() for _ in range(n_S_tree + 1)]
        Cld_S[1] = {2, 3, 4}
        Cld_S[2] = {5, 6, 7}
        Cld_S[3] = {8, 9, 10}
        Cld_S[4] = {11, 12, 13}

        Cld_T = [set() for _ in range(n_T_tree + 1)]
        Cld_T[1] = {2, 3}
        Cld_T[2] = {4, 5, 6}
        Cld_T[3] = {7, 8, 9}

        prt_S = [0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
        prt_T = [0, 0, 1, 1, 2, 2, 2, 3, 3, 3]

        P_T_prc = [set() for _ in range(n_T_tree + 1)]
        P_T_prc[1] = {2}
        P_T_prc[2] = {3, 4}
        P_T_prc[3] = {7}
        P_T_prc[4] = {5}
        P_T_prc[5] = {6}
        P_T_prc[7] = {8}
        P_T_prc[8] = {9}

        P_S_prc = [set() for _ in range(n_S_tree + 1)]
        P_S_prc[1] = {2}
        P_S_prc[2] = {3, 5}
        P_S_prc[3] = {4, 8}
        P_S_prc[4] = {11}
        P_S_prc[5] = {6}
        P_S_prc[6] = {7}
        P_S_prc[8] = {9}
        P_S_prc[9] = {10}
        P_S_prc[11] = {12}
        P_S_prc[12] = {13}

        Dsn_S = [set() for _ in range(k_star + 1)]
        Dsn_S[2] = {5, 6, 7, 8, 9, 10, 11, 12, 13}
        Dsn_S[1] = {2, 3, 4}
        Dsn_S[0] = {1}

    
    Bcset = [(1, 2, 1), (1, 2, 2), (1, 2, 3), (1, 3, 1), (1, 3, 2), 
             (1, 4 ,1), (2, 2, 1), (2, 2, 2), (2, 2, 3), (2, 3, 1), 
             (2, 3, 2), (2, 4, 1), (3, 3 ,1), (3, 3, 2), (3, 4, 1), 
             (4, 4, 1)]

    return dia_star, d_max, k_star, bh_star, t_star, s_star,\
    c_star, root_num, bl_star, E_B, E1_plus, E1_minus, tail, head,\
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right, \
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S, \
    P_S_prc, P_T_prc, Bcset

def prepare_variables_4_2(
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):
    ''' A function to prepare variables used in section 4.2 '''
    ############### # 4.2 Variables from here ############### 
    a = {i: pulp.LpVariable("a({})".format(i), 0, 1, 
        cat=pulp.LpBinary) for i in (E_B)}

    e_st = {(s, t): pulp.LpVariable("e_st({},{})".format(s, t), 0, 1, 
        cat=pulp.LpBinary) for s in range(1, s_star + 1)
        for t in range(1, t_star + 1)}

    e_ts = {(t, s): pulp.LpVariable("e_ts({},{})".format(t, s), 0, 1, 
        cat=pulp.LpBinary) for s in range(1, s_star + 1) 
        for t in range(1, t_star + 1)}

    chi = {t: pulp.LpVariable("chi({})".format(t), 0, c_star, 
        cat=pulp.LpInteger) for t in range(1, t_star + 1)}

    delta_clr = {(t, c): pulp.LpVariable("delta_clr({},{})".format(t, 
        c), 0, 1, cat=pulp.LpBinary) for t in range(1, t_star + 1) 
        for c in range(0, c_star + 1)}

    clr = {c: pulp.LpVariable("clr({})".format(c), 0, t_star, 
        cat=pulp.LpInteger) for c in range(0, c_star + 1)}

    degbt_plus = {s: pulp.LpVariable("deg_bt_plus({})".format(s), 
        0, 4, cat=pulp.LpInteger) for s in range(1, s_star + 1)}

    degbt_minus = {s: pulp.LpVariable("deg_bt_minus({})".format(s), 
        0, 4, cat=pulp.LpInteger) for s in range(1, s_star + 1)}

    return a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus

def prepare_variables_4_3(
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):
    ''' A function to prepare variables used in section 4.3 '''
    # ############### 4.3 Variables from here ############### 

    tilde = {s: pulp.LpVariable("tilde({})".format(s), 0, 1, 
        cat=pulp.LpBinary) for s in range(1, s_star + 1)}
    
    u = {(s, i): pulp.LpVariable("u({},{})".format(s, i), 0, 1, 
        cat=pulp.LpBinary) for s in range(1, s_star + 1) 
        for i in range(1, n_S_tree + 1)}
    
    v = {(t, i): pulp.LpVariable("v({},{})".format(t, i), 0, 1, 
        cat=pulp.LpBinary) for t in range(1, t_star + 1) 
        for i in range(1, n_T_tree + 1)}
    
    e = {t: pulp.LpVariable("e({})".format(t), 0, 1, 
        cat=pulp.LpBinary) for t in range(1, t_star + 2)}

    return tilde, u, v, e

def prepare_variables_4_4(
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):
    ''' A function to prepare variables used in section 4.4 '''
    # ############### 4.4 Variables from here ###############
    alpha_tilde = {(p, i): \
        pulp.LpVariable("alpha_tilde({},{})".format(p, i), 0, 
            MAX_CODE, cat=pulp.LpInteger)
            for p in range(1, s_star + t_star + 1) 
            for i in range(1, n_S_tree + 1)}

    return  alpha_tilde

def prepare_variables_4_5(    
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):

    ''' A function to prepare variables used in section 4.5 '''
    # # # ############### 4.5 Variables from here ###############

    beta_tilde = {i: pulp.LpVariable("beta_tilde({})".format(i), 0, 
        MAX_BOND, cat=pulp.LpInteger) for i in (E_B)}
    beta_tilde_temp = {(p, i): \
        pulp.LpVariable("beta_tilde({},{})".format(p, i), 0, 
            MAX_BOND, cat=pulp.LpInteger)
            for p in range(1, s_star + t_star + 1) 
            for i in range(2, n_S_tree + 1)}
    beta_tilde.update(beta_tilde_temp)
    beta_tilde_temp = {(t, 1): \
        pulp.LpVariable("beta_tilde({},1)".format(t), 0, MAX_BOND, 
            cat=pulp.LpInteger) for t in range(1, t_star + 2)}
    beta_tilde.update(beta_tilde_temp)
    beta_hat = {(s, t): \
        pulp.LpVariable("beta_hat({},{})".format(s, t), 0, MAX_BOND, 
            cat=pulp.LpInteger)
            for s in range(1, s_star + 1) 
            for t in range(1, t_star + 1)}

    return beta_tilde, beta_hat

def prepare_variables_4_6(    
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):

    ''' A function to prepare variables used in section 4.6 '''
    # ############### 4.6 Variables from here ###############
    delta_alpha = {(p, i, atom):\
        pulp.LpVariable("delta_alpha({},{},{})".format(p, i, atom), 
            0, 1, cat=pulp.LpBinary)
            for p in range(1, s_star + t_star + 1) 
            for i in range(1, n_S_tree + 1) for atom in Extra_Lambda}
    
    delta_beta_tilde = {(i, k):\
        pulp.LpVariable("delta_beta_tilde({},{})".format(i, k), 0, 1,
        cat=pulp.LpBinary) for i in (E_B) 
        for k in range(0, MAX_BOND + 1)}
    
    delta_beta_tilde_temp = {(p, i, k):\
        pulp.LpVariable("delta_beta_tilde({},{},{})".format(p, i, k),
        0, 1, cat=pulp.LpBinary)
        for p in range(1, s_star + t_star + 1) 
        for i in range(2, n_S_tree + 1) 
        for k in range(0, MAX_BOND + 1)}
    
    delta_beta_tilde.update(delta_beta_tilde_temp)
    
    delta_beta_tilde_temp = {(t, 1, k):\
        pulp.LpVariable("delta_beta_tilde({},1,{}".format(t, k), 
            0, 1, cat=pulp.LpBinary)
            for t in range(1, t_star + 2) 
            for k in range(0, MAX_BOND + 1)}
    
    delta_beta_tilde.update(delta_beta_tilde_temp)
    
    delta_beta_hat = {(s, t, k):\
        pulp.LpVariable("delta_beta_hat({},{},{})".format(s, t, k), 
            0, 1, cat=pulp.LpBinary)
            for s in range(1, s_star + 1)
            for t in range(1, t_star + 1) 
            for k in range(0, MAX_BOND + 1)}
    
    return delta_alpha, delta_beta_tilde, delta_beta_hat

def prepare_variables_4_7(    
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):

    ''' A function to prepare variables used in section 4.7 '''
    # ############## 4.7 Variables from here ###############
    ce = {atom: pulp.LpVariable("ce({})".format(atom), 0, n_star, 
        cat=pulp.LpInteger) for atom in Lambda}

    # in/ex 07/23
    ce_in = {atom: pulp.LpVariable("ce_in({})".format(atom), 
        0, n_star, cat=pulp.LpInteger) for atom in Lambda}
    
    ce_ex = {atom: pulp.LpVariable("ce_ex({})".format(atom),
        0, n_star, cat=pulp.LpInteger) for atom in Lambda}
    
    Mass = pulp.LpVariable("mass")
    # in/ex 07/23
    
    b_in = {k: pulp.LpVariable("b_in({})".format(k), 0, 2 * n_star, 
        cat=pulp.LpInteger) for k in range(1, MAX_BOND + 1)}
    
    b_ex = {k: pulp.LpVariable("b_ex({})".format(k), 0, 2 * n_star, 
        cat=pulp.LpInteger) for k in range(1, MAX_BOND + 1)}
    
    b = {k: pulp.LpVariable("b({})".format(k), 0, 2 * n_star, 
        cat=pulp.LpInteger) for k in range(1, MAX_BOND + 1)}
    n_H = pulp.LpVariable('n_H', 0, 4 * n_star, cat=pulp.LpInteger)

    return ce_in, ce_ex, ce, Mass, b_in, b_ex, b, n_H

def prepare_variables_4_8(   
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):

    ''' A function to prepare variables used in section 4.8 '''
    # ############### 4.8 Variables from here ###############
    deg = {(p, i): pulp.LpVariable("deg({},{})".format(p, i), 0, 4, 
        cat=pulp.LpInteger) for p in range(1, s_star + t_star + 1) 
        for i in range(1, n_S_tree + 1)}

    delta_deg = {(p, i, d):\
        pulp.LpVariable("delta_deg({},{},{})".format(p, i, d), 
            0, 1, cat=pulp.LpBinary)
            for p in range(1, s_star + t_star + 1) 
            for i in range(1, n_S_tree + 1) 
            for d in range(MAX_VAL + 1)}

    # modified 20/07/23
    dg_in = {d: pulp.LpVariable("dg_in({})".format(d), 0, n_star, 
        cat=pulp.LpInteger) for d in range(1, 5)}

    # modified 20/07/23
    dg_ex = {d: pulp.LpVariable("dg_ex({})".format(d), 0, n_star, 
        cat=pulp.LpInteger) for d in range(1, 5)}

    return deg, delta_deg, dg_in, dg_ex

def prepare_variables_4_9(    
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):

    ''' A function to prepare variables used in section 4.9 '''
    # ############### 4.9 Variables from here ###############
    delta_tau = {(i, gamma):\
        pulp.LpVariable("delta_tau({},{})".format(i, gamma), 0, 1, 
            cat=pulp.LpBinary) for i in (E_B) for gamma in Gamma_ALL}
    delta_tau_temp = {(t, 1, gamma):\
        pulp.LpVariable("delta_tau({},{},{})".format(t, 1, gamma), 
            0, 1, cat=pulp.LpBinary)
            for t in range(2, t_star + 1) for gamma in Gamma_ALL}
    delta_tau.update(delta_tau_temp)
    delta_tau_temp = {(p, i, gamma):\
        pulp.LpVariable("delta_tau({},{},{})".format(p, i, gamma), 
            0, 1, cat=pulp.LpBinary)
            for p in range(1, s_star + t_star + 1) 
            for i in range(2, n_S_tree + 1) for gamma in Gamma_ALL}
    delta_tau.update(delta_tau_temp)
    delta_tau_hat = {(s, t, gamma):\
        pulp.LpVariable("delta_tau_hat({},{},{})".format(s, t, gamma),
        0, 1, cat=pulp.LpBinary)
        for s in range(1, s_star + 1) for t in range(1, t_star + 1) 
        for gamma in Gamma_ALL}
    ac_in = {gamma: pulp.LpVariable("ac_in({})".format(gamma), 
        0, n_star, cat=pulp.LpInteger)
        for gamma in (Gamma_less | Gamma_equal)}
    # for var in ac_in:
    #     print(var)
    ac_ex = {gamma: pulp.LpVariable("ac_ex({})".format(gamma), 
        0, n_star, cat=pulp.LpInteger)
        for gamma in (Gamma_less | Gamma_equal)}
    ac = {gamma: pulp.LpVariable("ac({})".format(gamma), 
        0, n_star, cat=pulp.LpInteger)
        for gamma in (Gamma_less | Gamma_equal)}

    return delta_tau, delta_tau_hat, ac_in, ac_ex, ac

def prepare_variables_4_10(    
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL
    ):

    ''' A function to prepare variables used in section 4.10 '''
    # ############### 4.10 Variables from here ###############
    bc = {(d1, d2, k):\
        pulp.LpVariable("bc({},{},{})".format(d1, d2, k), 
            0, n_star - 1, cat=pulp.LpInteger)
            for (d1, d2, k) in Bcset}

    bc_in = {(d1, d2, k):\
        pulp.LpVariable("bc_in({},{},{})".format(d1, d2, k), 
            0, n_star - 1, cat=pulp.LpInteger)
            for (d1, d2, k) in Bcset}

    bc_ex = {(d1, d2, k):\
        pulp.LpVariable("bc_ex({},{},{})".format(d1, d2, k), 
            0, n_star - 1, cat=pulp.LpInteger)
            for (d1, d2, k) in Bcset}
    
    delta_dc = {(i, d1, d2, k):\
        pulp.LpVariable("delta_dc({},{},{},{})".format(i, d1, d2, k),
            0, 1, cat=pulp.LpBinary)
            for i in (E_B) for d1 in range(MAX_VAL + 1) 
            for d2 in range(MAX_VAL + 1) for k in range(5)}
    
    delta_dc_temp = {(t, 1, d1, d2, k):\
        pulp.LpVariable("delta_dc({},{},{},{},{})".format(t, 1, 
            d1, d2, k), 0, 1, cat=pulp.LpBinary)
            for t in range(2, t_star + 1) 
            for d1 in range(MAX_VAL + 1) for d2 in range(MAX_VAL + 1)
            for k in range(5)}
    
    delta_dc.update(delta_dc_temp)
    
    delta_dc_temp = {(p, i, d1, d2, k):\
        pulp.LpVariable("delta_dc({},{},{},{},{})".format(p, i, d1, 
            d2, k), 0, 1, cat=pulp.LpBinary)
            for p in range(1, s_star + t_star + 1) 
            for i in range(2, n_S_tree + 1) 
            for d1 in range(MAX_VAL + 1)
            for d2 in range(MAX_VAL + 1) for k in range(5)}
    
    delta_dc.update(delta_dc_temp)
    
    delta_dc_hat = {(s, t, d1, d2, k):\
            pulp.LpVariable("delta_dc_hat({},{},{},{},{})".format(s, 
                t, d1, d2, k), 0, 1, cat=pulp.LpBinary)
                for s in range(1, s_star + 1) 
                for t in range(1, t_star + 1) 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1) for k in range(5)}

    return bc, bc_in, bc_ex, delta_dc, delta_dc_hat

def prepare_variables_B_8(
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset, 
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a_3, a_1_3, Ncset, Ncset_e, ce_nc, dg_nc, ml
    ):
    ''' A function to prepare variables used in section B.8 '''
    ############### # B.8 Variables from here ###############
    nc = {nu: pulp.LpVariable("nc({})".format(nu), 0, n_star, 
        cat=pulp.LpInteger) for nu in Ncset_e}

    delta_nc = {(p, i, nu):\
        pulp.LpVariable("delta_nc({},{},{})".format(p, i, nu),
            0, 1, cat=pulp.LpBinary)
            for p in range(1, s_star + t_star + 1) 
            for i in range(1, n_S_tree + 1) 
            for nu in Ncset_e}

    return nc, delta_nc

def add_constraints_4_2(
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus):

    ''' A function to add constraints in section 4.2 '''
     ############### 4.2 constraints from here ###############
    # (1) checked 200305
    for t in range(1, t_star + 1):
        MILP += \
            pulp.lpSum([delta_clr[(t, c)] 
                for c in range(0, c_star + 1)]) == 1, \
            "milp-(8)-(1)-{}".format(t)

    # (2) checked 200305
    for t in range(1, t_star + 1):
        MILP += \
            pulp.lpSum([c * delta_clr[(t, c)] 
                for c in range(0, c_star + 1)]) == chi[t], \
            "milp-(8)-(2)-{}".format(t)
    
    # (3) checked 200305
    for c in range(0, c_star + 1):
        MILP += \
            pulp.lpSum([delta_clr[(t, c)] 
                for t in range(1, t_star + 1)]) == clr[c], \
            "milp-(9)-{}".format(c)
    
    # (4) 
    for i in range(1, c_star + 1):
        MILP += t_star * (1 - a[i]) >= clr[i], "milp-(10)-{}".format(i)
    
    # (5) checked 200305
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += e_st[(s, t)] + e_ts[(t, s)] <= 1,\
            "milp-(11)-{}-{}".format(s, t)
    
    # (6) checked 200305
    for c in range(1, c_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([e_ts[(t, s)] 
                for s in range(1, s_star + 1) 
                if s != head[c]]) <= 1 - delta_clr[(t, c)], \
            "milp-(12)-(1)-{}-{}".format(t, c)
    
    # (7) checked 200305
    for c in range(1, c_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([e_st[(s, t)] 
                for s in range(1, s_star + 1) 
                if s != tail[c]]) <= 1 - delta_clr[(t, c)],\
                "milp-(12)-(2)-{}-{}".format(t, c)
    
    # (8) checked 200305
    for s in range(1, s_star + 1):
        MILP += \
            pulp.lpSum([a[i] for i in E1_minus[s]]) + \
            pulp.lpSum([e_ts[(t, s)]for t in range(1, t_star + 1)]) \
            == degbt_minus[s], "milp-(13)-(1)-{}".format(s)
    
    #(9) checked 200305
    for s in range(1, s_star + 1):
        MILP += \
            pulp.lpSum([a[i] for i in E1_plus[s]]) + \
            pulp.lpSum([e_st[(s, t)] for t in range(1, t_star + 1)]) \
            == degbt_plus[s], "milp-(13)-(2)-{}".format(s)
    
    # (10) checked 200305
    for s in range(1, s_star + 1):
        MILP += degbt_minus[s] + degbt_plus[s] <= d_max,\
        "milp-(13)-(3)-{}".format(s)

    return MILP
    
def add_constraints_4_3(    
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e
    ):

    ''' A function to add constraints in section 4.3 '''
     ############### 4.3 constraints from here ###############
    # (11)  P_S_prc
    for s in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            for j in P_S_prc[i]:
                MILP += u[(s, i)] >= u[(s, j)],\
                "milp-(14)-{}-{}-{}".format(s, i, j)
    
    # (12)  P_T_prc
    for t in range(1, t_star + 1):
        for i in range(1, n_T_tree + 1):
            for j in P_T_prc[i]:
                MILP += v[(t, i)] >= v[(t, j)], \
                "milp-(15)-{}-{}-{}".format(t, i, j)
    
    # (13)
    MILP += pulp.lpSum([u[(s, i)] for s in range(1, s_star + 1) 
                for i in range(1, n_S_tree + 1)]) + \
            pulp.lpSum([v[(t, i)] for t in range(1, t_star + 1) 
                for i in range(1, n_T_tree + 1)]) == n_star, \
            "milp-(16)"
    
    # (14)
    for s in range(1, s_star + 1):
        MILP += pulp.lpSum([u[(s, i)] 
            for i in range(1, n_S_tree + 1)]) <=\
                2 + 2 * pulp.lpSum([u[(s, j)] 
                    for j in Cld_S[1]]), "milp-(17)-{}".format(s)
    
    # (15)
    for t in range(1, t_star + 1):
        MILP += pulp.lpSum([v[(t, i)] 
            for i in range(1, n_T_tree + 1)]) <=\
                2 + 2 * pulp.lpSum([v[(t, j)] 
                    for j in Cld_T[1]]), "milp-(18)-{}".format(t)
    
    
    # (16)
    MILP += e[1] == 0, "milp-(19)-(1)"
    MILP += e[t_star + 1] == 0, "milp-(19)-(2)"
    
    # (17)
    for t in range(1, t_star + 1):
        MILP += \
            e[t + 1] + pulp.lpSum([e_ts[(t, s)] 
                for s in range(1, s_star + 1)]) == v[(t, 1)], \
            "milp-(19)-(3)-{}".format(t)
    
    # (18)
    for t in range(1, t_star + 1):
        MILP += \
            e[t] + pulp.lpSum([e_st[(s, t)] 
                for s in range(1, s_star + 1)]) == v[(t, 1)], \
            "milp-(19)-4-{}".format(t)

    # new 20/07/23
    for t in range(1, t_star + 1):
        MILP += pulp.lpSum([delta_clr[(t, c)]
                for c in range(1, c_star + 1)]) == v[(t, 1)], \
                "milp-(20)-{}".format(t)
    
    #(21)
    for t in range(1, t_star):
        MILP += chi[t] - chi[t + 1] >= v[(t, 1)] - e[t + 1] , \
        "milp-(21)-{}-lb".format(t)

    # (21)
    for t in range(1, t_star):
        MILP += c_star * (1 - e[t + 1]) >= chi[t] - chi[t + 1], \
        "milp-(21)-{}-ub".format(t)

    # (22) 20/08/06
    for i in range(1, c_star + 1):
        MILP += a[i] + pulp.lpSum([e_ts[(t, i + 1)]
                for t in range(1, t_star + 1)]) == u[(i + 1, 1)], \
                "milp-(22)-{}".format(i)
    

    # (23) only for root s = 1 diameter is even number
    MILP += tilde[1] == 1, "milp-(23)-(1)"
    MILP += u[(1, 1)] == 1, "milp-(23)-(2)"

    # (24)
    for s in range(1, s_star + 1):
        MILP += tilde[s] <= u[(s, 1)], "milp-(24)-{}".format(s)

    # (25)only for root s = 1 diameter is even number
    for s in range(2, s_star + 1):
        expr = pulp.lpSum(u[(s_dush, 1)] for s_dush in Cld_B[s])
        MILP += (d_max - 1) * tilde[s] >= expr,  \
            "milp-(25)-(1)-{}-ub".format(s)
        MILP += expr >= 2 * tilde[s], "milp-(25)-(1)-{}-lb".format(s)

    #(26)
    for s in range(2, s_star + 1):
        MILP += pulp.lpSum(u[(s, i)] 
            for i in Dsn_S[k_star]) >= u[(s, 1)] - tilde[s],\
        "milp-(25)-(2)-{}".format(s)

    # (27) modified 20/07/22
    MILP += pulp.lpSum(u[(s, 1)] - tilde[s] 
        for s in range(2, s_star + 1)) == bl_star, \
        "milp-(26)-(1)"

    #(28)
    MILP += pulp.lpSum(u[(s, 1)] for s in V_B_bh_star) >= 1, \
        "milp-(26)-(2)"

    #(29)
    MILP += pulp.lpSum(u[(s, 1)] for s in V_B_left) + \
        pulp.lpSum(clr[i] for i in E_B_left) \
        == math.ceil(dia_star/2) - k_star, "milp-(27)-(1)"
    
    #(30)
    MILP += pulp.lpSum(u[(s, 1)] for s in V_B_right) + \
        pulp.lpSum(clr[i] for i in E_B_right) \
        == math.floor(dia_star/2) - k_star, "milp-(27)-(2)"

    # modified 20/07/28
    if bh_star == 2:
        MILP += clr[4] <= clr[5], "milp-(28)-(1)"
        MILP += clr[9] <= clr[8], "milp-(28)-(2)"
        MILP += clr[6] <= clr[7], "milp-(28)-(3)"
    
        MILP += u[(3, 1)] + u[(6, 1)] + clr[2] + clr[6] <= \
            math.floor(dia_star/2) - k_star, "milp-(28)-(4)"
    elif bh_star == 1:
        MILP += u[(3, 1)] + clr[2] <= \
            math.floor(dia_star/2) - k_star, "milp-(28)"
        # for s in V_B_bh_star - V_B_left - V_B_right:
        #     MILP += pulp.lpSum(u[(s1, 1)]
        #         for s1 in V_B_s[s]) + \
        #         pulp.lpSum(clr[i] for i in E_B_s[s]) <= \
        #         math.floor(dia_star/2) - k_star, "milp-(28)-{}".format(s)

    return MILP

def add_constraints_4_5(   
    MILP, 
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat
    ):

    # A function to add constraints in section 4.5
    # (31)
    for i in (E_B):
        MILP += beta_tilde[i] >= a[i], "milp-(29)-{}-lb".format(i)
        MILP += beta_tilde[i] <= 3 * a[i], "milp-(29)-{}-ub".format(i)
    
    # (32)
    for s in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += beta_tilde[(s, i)] >= u[(s, i)],\
                "milp-(30)-{}-{}-lb".format(s, i)
            MILP += beta_tilde[(s, i)] <= 3 * u[(s, i)],\
                "milp-(30)-{}-{}-ub".format(s, i)
    
    # (33)
    for t in range(1, t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += beta_tilde[(s_star + t, i)] >= v[(t, i)],\
                "milp-(31)-{}-{}-lb".format(t, i)
            MILP += beta_tilde[(s_star + t, i)] <= 3 * v[(t, i)], \
                "milp-(31)-{}-{}-ub".format(t, i)    
    
    # (34)
    for t in range(1, t_star + 2):
        MILP += beta_tilde[(t, 1)] >= e[t], \
            "milp-(32)-{}-lb".format(t)
        MILP += beta_tilde[(t, 1)] <= 3 * e[t], \
            "milp-(32)-{}-ub".format(t)    
    
    # (35)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += beta_hat[(s, t)] >= e_st[(s, t)] + e_ts[(t, s)], \
                "milp-(33)-{}-{}-lb".format(s, t)
            MILP += beta_hat[(s, t)] <= 3 * (e_st[(s, t)] +
                e_ts[(t, s)]), "milp-(33)-{}-{}-ub".format(s, t)

    return MILP

#B.5
def add_constraints_4_6(
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    delta_alpha, delta_beta_tilde, delta_beta_hat
    ):

    ''' A function to add constraints in section 4.6 '''
    # (36)
    for p in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            MILP += pulp.lpSum([delta_alpha[(p, i, atom)] 
                for atom in Extra_Lambda]) == 1, \
            "milp-(34)-{}-{}".format(p, i)
    for p in range(s_star + 1,  s_star + t_star + 1):
        for i in range(1, n_T_tree + 1):
            MILP += pulp.lpSum([delta_alpha[(p, i, atom)] 
                for atom in Extra_Lambda]) == 1, \
            "milp-(34)-{}-{}".format(p, i)
    
    # (37)
    for p in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            MILP += pulp.lpSum([Code_Lambda[atom] * 
                delta_alpha[(p, i, atom)] for atom in Extra_Lambda]) \
                == alpha_tilde[(p, i)], "milp-(35)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(1, n_T_tree + 1):
            MILP += pulp.lpSum([Code_Lambda[atom] * 
                delta_alpha[(p, i, atom)] for atom in Extra_Lambda]) \
                == alpha_tilde[(p, i)], "milp-(35)-{}-{}".format(p, i)

    # (38)
    for i in (E_B):
        MILP += pulp.lpSum([delta_beta_tilde[(i, k)] 
            for k in range(0, MAX_BOND + 1)]) == 1, \
        "milp-(36)-(1)-{}".format(i)
    
    # (39)
    for i in (E_B):
        MILP += pulp.lpSum([k * delta_beta_tilde[(i, k)] 
            for k in range(1, MAX_BOND + 1)]) == \
                beta_tilde[i], "milp-(36)-(2)-{}".format(i)
    
    # (40)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([delta_beta_tilde[(p, i, k)] 
                for k in range(0, MAX_BOND + 1)]) == 1, \
            "milp-(37)-(1)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([delta_beta_tilde[(p, i, k)] 
                for k in range(0, MAX_BOND + 1)]) == 1, \
            "milp-(37)-(1)-{}-{}".format(p, i)
    
    # (41)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([k * delta_beta_tilde[(p, i, k)]
                for k in range(1, MAX_BOND + 1)]) == \
                    beta_tilde[(p, i)], "milp-(37)-(2)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([k * delta_beta_tilde[(p, i, k)]
                for k in range(1, MAX_BOND + 1)]) == \
                    beta_tilde[(p, i)], "milp-(37)-(2)-{}-{}".format(p, i)
    
    # (42)
    for t in range(1, t_star + 2):
        MILP += pulp.lpSum([delta_beta_tilde[(t, 1, k)] 
            for k in range(0, MAX_BOND + 1)]) == 1, \
        "milp-(38)-(1)-{}".format(t)
    
    # (43)
    for t in range(1, t_star + 2):
        MILP += pulp.lpSum([k * delta_beta_tilde[(t, 1, k)] 
            for k in range(1, MAX_BOND + 1)]) == \
                beta_tilde[(t, 1)], "milp-(38)-(2)-{}".format(t)
    
    # (44)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([delta_beta_hat[(s, t, k)] 
                for k in range(0, MAX_BOND + 1)]) == 1, \
            "milp-(39)-(1)-{}-{}".format(s, t)
    
    # (45)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([k * delta_beta_hat[(s, t, k)] 
                for k in range(0, MAX_BOND + 1)]) == \
                    beta_hat[(s, t)], "milp-(39)-(2)-{}-{}".format(s, t)

    # (46) #infeasible
    for s in range(1, s_star + 1):
        MILP += pulp.lpSum([beta_tilde[i] 
                    for i in (E1_plus[s] | E1_minus[s])]) + \
                pulp.lpSum([beta_hat[(s, t)] 
                    for t in range(1, t_star + 1)]) + \
                pulp.lpSum([beta_tilde[(s, j)] 
                    for j in Cld_S[1]]) <= \
                pulp.lpSum([val[atom] * delta_alpha[(s, 1, atom)] 
                    for atom in Lambda]), "milp-(40)-{}".format(s)
        MILP += MAX_VAL * (pulp.lpSum([beta_tilde[i] 
                    for i in (E1_plus[s] | E1_minus[s])]) +
                pulp.lpSum([beta_hat[(s, t)] 
                    for t in range(1, t_star + 1)]) +
                pulp.lpSum([beta_tilde[(s, j)] 
                    for j in Cld_S[1]])) >= \
                pulp.lpSum([val[atom] * delta_alpha[(s, 1, atom)] 
                    for atom in Lambda]), "milp-(40)-{}-max".format(s)
    
    # (47) #infeasible
    for t in range(1, t_star + 1):
        MILP += pulp.lpSum([beta_hat[(s, t)] 
                    for s in range(1, s_star + 1)]) + \
                beta_tilde[(t, 1)] + beta_tilde[(t + 1, 1)] + \
                pulp.lpSum([beta_tilde[(s_star + t, j)] 
                    for j in Cld_T[1]]) <= \
                pulp.lpSum([val[atom] * 
                    delta_alpha[(s_star + t, 1, atom)] 
                    for atom in Lambda]), "milp-(41)-{}".format(t)
        MILP += MAX_VAL * (pulp.lpSum([beta_hat[(s, t)] 
                    for s in range(1, s_star + 1)]) +
                beta_tilde[(t, 1)] + beta_tilde[(t + 1, 1)] +
                pulp.lpSum([beta_tilde[(s_star + t, j)] 
                    for j in Cld_T[1]])) >= \
                pulp.lpSum([val[atom] * 
                    delta_alpha[(s_star + t, 1, atom)] 
                    for atom in Lambda]), "milp-(41)-{}-max".format(t)
    
    # (48) #infeasible
    for s in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += beta_tilde[(s, i)] +\
                        pulp.lpSum([beta_tilde[(s, j)] 
                            for j in Cld_S[i]]) <= \
                    pulp.lpSum([val[atom] * delta_alpha[(s, i, atom)] 
                        for atom in Lambda]), \
                    "milp-(42)-{}-{}".format(s, i)
            MILP += MAX_VAL * (beta_tilde[(s, i)] + 
                        pulp.lpSum([beta_tilde[(s, j)] 
                            for j in Cld_S[i]])) >= \
                    pulp.lpSum([val[atom] * delta_alpha[(s, i, atom)] 
                        for atom in Lambda]), \
                    "milp-(42)-{}-{}-max".format(s, i)
    
    # (49) #infeasible
    for t in range(1, t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += beta_tilde[(s_star + t, i)] + \
                    pulp.lpSum([beta_tilde[(s_star + t, j)] 
                        for j in Cld_T[i]]) <= \
                    pulp.lpSum([val[atom] * 
                        delta_alpha[(s_star + t, i, atom)] 
                        for atom in Lambda]),\
                        "milp-(43)-{}-{}".format(t, i)
            MILP += MAX_VAL * (beta_tilde[(s_star + t, i)] + 
                        pulp.lpSum([beta_tilde[(s_star + t, j)] 
                            for j in Cld_T[i]])) >= \
                    pulp.lpSum([val[atom] * 
                        delta_alpha[(s_star + t, i, atom)] 
                        for atom in Lambda]), \
                    "milp-(43)-{}-{}-max".format(t, i)
    
    return MILP

#B.6 
def add_constraints_4_7(
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    delta_alpha, delta_beta_tilde, delta_beta_hat,
    ce_in, ce_ex, ce, Mass, b_in, b_ex, b, n_H
    ):

    ''' A function to add constraints in section 4.7 '''
    # (50) modified 20/07/23
    for atom in Lambda:
        MILP += pulp.lpSum([delta_alpha[(p, 1, atom)] 
            for p in range(1, s_star + t_star + 1)]) == \
                ce_in[atom], "milp-(44)-(1)-{}".format(atom)
    
    # (51) modified 20/07/23
    for atom in Lambda:
        MILP += pulp.lpSum([delta_alpha[(p, i, atom)] 
            for p in range(1, s_star + 1) 
            for i in range(2, n_S_tree + 1)]) + \
            pulp.lpSum([delta_alpha[(p, i, atom)] 
            for p in range(s_star + 1, s_star + t_star + 1) 
            for i in range(2, n_T_tree + 1)]) == \
                ce_ex[atom], "milp-(44)-(2)-{}".format(atom)

    # for atom in Lambda:
    #     MILP += ce[atom] == ce_in[atom] + ce_ex[atom], \
    #         "milp-(60)-2-{}".format(atom)
    
    # (52) modified 07/23
    MILP += pulp.lpSum([mass[atom] * (ce_in[atom] + ce_ex[atom]) 
        for atom in Lambda]) == Mass, "milp-(45)"
    
    # (53) modified 07/23
    for k in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum([delta_beta_tilde[(i, k)] 
                    for i in (E_B)]) + \
                pulp.lpSum([delta_beta_hat[(s, t, k)] 
                    for s in range(1, s_star + 1) 
                    for t in range(1, t_star + 1)]) + \
                pulp.lpSum([delta_beta_tilde[(t, 1, k)] 
                    for t in range(2, t_star + 1)]) == b_in[k], \
                "milp-(46)-{}".format(k)
    
    # (54) modified 07/23
    for k in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum([delta_beta_tilde[(p, i, k)] 
                    for p in range(1, s_star + 1) 
                    for i in range(2, n_S_tree + 1)]) + \
                pulp.lpSum([delta_beta_tilde[(p, i, k)] 
                    for p in range(s_star + 1, s_star + t_star + 1) 
                    for i in range(2, n_T_tree + 1)])  == \
                b_ex[k], "milp-(47)-{}".format(k)
    
    # # (55) 
    # for k in range(1, MAX_BOND + 1):
    #     MILP += b_in[k] + b_ex[k] == b[k], "milp-(63)-2-{}".format(k)

    # 20/07/23
    MILP += pulp.lpSum([val[atom] * (ce_in[atom] + ce_ex[atom])
            for atom in Lambda]) - \
            2 * (n_star - 1) - 2 * b_in[2] - 2 * b_ex[2] - \
            2 * b_in[3] - 2 * b_ex[3] == n_H, \
            "milp-(48)".format(k)

    return MILP

# B.7
def add_constraints_4_8(
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    delta_alpha, delta_beta_tilde, delta_beta_hat,
    ce_in, ce_ex, ce, Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex
    ):

    ''' A function to add constraints in section 4.8 '''
    # # (56)
    for s in range(1, s_star + 1):
        MILP += pulp.lpSum([a[i] 
                    for i in (E1_minus[s] | E1_plus[s])]) + \
                pulp.lpSum([(e_st[(s, t)] + e_ts[(t, s)]) 
                    for t in range(1, t_star + 1)]) + \
                pulp.lpSum([u[(s, j)] 
                    for j in Cld_S[1]]) == deg[(s, 1)], \
                "milp-(49)-{}".format(s)
    
    # (57)
    for s in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += u[(s, i)] + pulp.lpSum([u[(s, j)] 
                for j in Cld_S[i]]) == deg[(s, i)], \
            "milp-(50)-{}-{}".format(s, i)
    
    # (58)
    for t in range(1, t_star + 1):
        MILP += 2 * v[(t, 1)] + pulp.lpSum([v[(t, j)] 
            for j in Cld_T[1]]) == deg[(s_star + t, 1)], \
        "milp-(51)-{}".format(t)
    
    # (59)
    for t in range(1, t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += v[(t, i)] + pulp.lpSum([v[(t, j)] 
                for j in Cld_T[i]]) == deg[(s_star + t, i)],\
            "milp-(52)-{}-{}".format(t, i)
    
    # (60)
    for p in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            MILP += pulp.lpSum([delta_deg[(p, i, d)] 
                for d in range(0, MAX_VAL + 1)]) == 1, \
            "milp-(53)-(1)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(1, n_T_tree + 1):
            MILP += pulp.lpSum([delta_deg[(p, i, d)] 
                for d in range(0, MAX_VAL + 1)]) == 1, \
            "milp-(53)-(1)-{}-{}".format(p, i)
    
    # (61)
    for p in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            MILP += pulp.lpSum([d * delta_deg[(p, i, d)] 
                for d in range(1, MAX_VAL + 1)]) == \
                    deg[(p, i)], "milp-(53)-(2)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(1, n_T_tree + 1):
            MILP += pulp.lpSum([d * delta_deg[(p, i, d)] 
                for d in range(1, MAX_VAL + 1)]) == \
                    deg[(p, i)], "milp-(53)-(2)-{}-{}".format(p, i)
    
    # modified 20/07/23
    for d in range(1, MAX_VAL + 1):
        MILP += pulp.lpSum([delta_deg[(p, 1, d)] 
            for p in range(1, s_star + t_star + 1)])== \
                dg_in[d], "milp-(54)-(1)-{}".format(d)

    # modified 20/07/23
    for d in range(1, MAX_VAL + 1):
        MILP += pulp.lpSum([delta_deg[(p, i, d)] 
            for p in range(1, s_star + 1) 
            for i in range(2, n_S_tree + 1)]) + \
            pulp.lpSum([delta_deg[(p, i, d)] 
            for p in range(s_star + 1, s_star + t_star + 1) 
            for i in range(2, n_T_tree + 1)]) == \
                dg_ex[d], "milp-(54)-(2)-{}".format(d)
    
    # 20/07/23
    if d_max == 3:
        MILP += dg_in[4] + dg_ex[4] == 0, "milp-(55)"
    else:
        MILP += dg_in[4] + dg_ex[4] >= 1, "milp-(55)"

    return MILP


# B.8
def add_constraints_4_9(
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    delta_alpha, delta_beta_tilde, delta_beta_hat,
    ce_in, ce_ex, ce, Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex,
    delta_tau, delta_tau_hat, ac_in, ac_ex, ac, n_H 
    ):

    ''' A function to add constraints in section 4.9 '''
    #  ############### 4.9 constraints from here ###############
    # (64)
    for i in (E_B):
        MILP += pulp.lpSum([delta_tau[(i, gamma)] 
            for gamma in Gamma_ALL]) == 1, "milp-(56)-(1)-{}".format(i)
    
    # (65)
    for i in (E_B):
        MILP += pulp.lpSum([Code_Lambda[atom1] *
            delta_tau[(i, (atom1, atom2, k))] 
            for (atom1, atom2, k) in Gamma_ALL]) == \
                alpha_tilde[(tail[i], 1)], "milp-(56)-(2)-{}".format(i)
    
    #(66)
    for i in (E_B):
        MILP += pulp.lpSum([Code_Lambda[atom2] * 
            delta_tau[(i, (atom1, atom2, k))] 
            for (atom1, atom2, k) in Gamma_ALL]) == \
                alpha_tilde[(head[i], 1)], "milp-(56)-(3)-{}".format(i)
    
    # (67)
    for i in (E_B):
        MILP += pulp.lpSum([k * delta_tau[(i, (atom1, atom2, k))] 
            for (atom1, atom2, k) in Gamma_ALL]) == \
                beta_tilde[i], "milp-(56)-(4)-{}".format(i)
    
    # (68)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([delta_tau[(t, 1, gamma)] 
            for gamma in Gamma_ALL]) == 1, "milp-(57)-(1)-{}".format(t)
    
    # (69)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([Code_Lambda[atom1] * 
            delta_tau[(t, 1, (atom1, atom2, k))] 
            for (atom1, atom2, k) in Gamma_ALL]) \
                == alpha_tilde[(s_star + t - 1, 1)], \
                "milp-(57)-(2)-{}".format(t)
    
    # (70)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([Code_Lambda[atom2] * 
            delta_tau[(t, 1, (atom1, atom2, k))] 
            for (atom1, atom2, k) in Gamma_ALL]) \
                == alpha_tilde[(s_star + t, 1)], \
                "milp-(57)-(3)-{}".format(t)
    
    # (71)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([k * delta_tau[(t, 1, (atom1, atom2, k))] 
            for (atom1, atom2, k) in Gamma_ALL]) \
                == beta_tilde[(t, 1)], "milp-(57)-(4)-{}".format(t)
    
    # (72)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([delta_tau[(p, i, gamma)] 
                for gamma in Gamma_ALL]) == 1, \
            "milp-(58)-(1)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([delta_tau[(p, i, gamma)] 
                for gamma in Gamma_ALL]) == 1, \
            "milp-(58)-(1)-{}-{}".format(p, i)
    
    # (73)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([Code_Lambda[atom1] * 
                delta_tau[(p, i, (atom1, atom2, k))] 
                for (atom1, atom2, k) in Gamma_ALL]) == \
            alpha_tilde[(p, prt_S[i])], "milp-(58)-(2)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([Code_Lambda[atom1] * 
                delta_tau[(p, i, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            alpha_tilde[(p, prt_T[i])], "milp-(58)-(2)-{}-{}".format(p, i)
    
    # (74)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([Code_Lambda[atom2] * 
                delta_tau[(p, i, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            alpha_tilde[(p, i)], "milp-(58)-(3)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([Code_Lambda[atom2] * 
                delta_tau[(p, i, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            alpha_tilde[(p, i)], "milp-(58)-(3)-{}-{}".format(p, i)
    
    # (75) #infeasible
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([k * 
                delta_tau[(p, i, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            beta_tilde[(p, i)], "milp-(58)-(4)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([k * 
                delta_tau[(p, i, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            beta_tilde[(p, i)], "milp-(58)-(4)-{}-{}".format(p, i)
    
    # (76)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([delta_tau_hat[s, t, gamma] 
                for gamma in Gamma_ALL]) == 1, \
            "milp-(59)-(1)-{}-{}".format(s, t)
    
    # (77)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([Code_Lambda[atom1] * 
                delta_tau_hat[(s, t, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            alpha_tilde[(s, 1)], "milp-(59)-(2)-{}-{}".format(s, t)
    
    # (78) #infeasible
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([Code_Lambda[atom2] * 
                delta_tau_hat[(s, t, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            alpha_tilde[(s_star + t, 1)], \
            "milp-(59)-(3)-{}-{}".format(s, t)
    
    # (79) #infeasible
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([k * 
                delta_tau_hat[(s, t, (atom1, atom2, k))]
                for (atom1, atom2, k) in Gamma_ALL]) == \
            beta_hat[(s, t)], "milp-(59)-(4)-{}-{}".format(s, t)
    
    # (80) modified 20/07/23
    for (atom1, atom2, k) in Gamma_less:
        MILP += pulp.lpSum([delta_tau[(i, (atom1, atom2, k))] + 
                    delta_tau[(i, (atom2, atom1, k))] 
                    for i in E_B]) + \
                pulp.lpSum([delta_tau_hat[(s, t, (atom1, atom2, k))] + 
                    delta_tau_hat[(s, t, (atom2, atom1, k))]
                    for s in range(1, s_star + 1) 
                    for t in range(1, t_star + 1)]) + \
                pulp.lpSum([delta_tau[(t, 1, (atom1, atom2, k))] + 
                    delta_tau[(t, 1, (atom2, atom1, k))]
                    for t in range(2, t_star + 1)]) == \
                ac_in[(atom1, atom2, k)], \
                "milp-(60)-{}".format((atom1, atom2, k))
    
    # (81) modified 20/07/23
    for gamma in Gamma_equal:
        MILP += pulp.lpSum([delta_tau[(i, gamma)] for i in E_B]) + \
                pulp.lpSum([delta_tau_hat[(s, t, gamma)] 
                    for s in range(1, s_star + 1) 
                    for t in range(1, t_star + 1)]) + \
                pulp.lpSum([delta_tau[(t, 1, gamma)] 
                    for t in range(2, t_star + 1)]) == \
                ac_in[gamma], "milp-(61)-{}".format(gamma)
    
    # (82) modified 20/07/23
    for (atom1, atom2, k) in Gamma_less:
        MILP += pulp.lpSum([delta_tau[(p, i, (atom1, atom2, k))] + \
                    delta_tau[(p, i, (atom2, atom1, k))]
                    for p in range(1, s_star + 1) 
                    for i in range(2, n_S_tree + 1)]) + \
                pulp.lpSum([delta_tau[(p, i, (atom1, atom2, k))] + \
                    delta_tau[(p, i, (atom2, atom1, k))]
                    for p in range(s_star + 1, s_star + t_star + 1) 
                    for i in range(2, n_T_tree + 1)]) == \
                ac_ex[(atom1, atom2, k)], \
                "milp-(62)-{}".format((atom1, atom2, k))
    
    # (83)
    for gamma in Gamma_equal:
        MILP += pulp.lpSum([delta_tau[(p, i, gamma)] 
            for p in range(1, s_star + 1) 
            for i in range(2, n_S_tree + 1)]) + \
                pulp.lpSum([delta_tau[(p, i, gamma)] 
            for p in range(s_star + 1, s_star + t_star + 1) 
            for i in range(2, n_T_tree + 1)]) == \
                ac_ex[gamma], "milp-(63)-{}".format(gamma)
    
    return MILP


#B.9
def add_constraints_4_10(
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    delta_alpha, delta_beta_tilde, delta_beta_hat,
    ce_in, ce_ex, ce, Mass, b_in, b_ex, b, deg, delta_deg, dg_in, dg_ex,
    delta_tau, delta_tau_hat, ac_in, ac_ex, ac, n_H,
    bc, bc_in, bc_ex, delta_dc, delta_dc_hat 
    ):

    ''' A function to add constraints in section 4.10 '''
    # (85)
    for i in (E_B):
        MILP += pulp.lpSum([delta_dc[(i, d1, d2, k)] 
            for d1 in range(MAX_VAL + 1) 
            for d2 in range(MAX_VAL + 1)
            for k in range(4)]) == 1, "milp-(64)-(1)-{}".format(i)
    
    # (86)
    for i in (E_B):
        MILP += pulp.lpSum([k * delta_dc[(i, d1, d2, k)] 
            for d1 in range(MAX_VAL + 1) 
            for d2 in range(MAX_VAL + 1)
            for k in range(4)]) == beta_tilde[i], \
        "milp-(64)-(2)-{}".format(i)
    
    # (87)
    for i in (E_B):
        MILP += pulp.lpSum([d1 * delta_dc[(i, d1, d2, k)] 
            for d1 in range(1, MAX_VAL + 1) 
            for d2 in range(MAX_VAL + 1)\
            for k in range(4)]) == deg[(tail[i], 1)], \
        "milp-(64)-(3)-{}".format(i)
    
    # (88)
    for i in (E_B):
        MILP += pulp.lpSum([d2 * delta_dc[(i, d1, d2, k)] 
            for d1 in range(MAX_VAL + 1) 
            for d2 in range(1, MAX_VAL + 1)
            for k in range(4)]) == deg[head[i], 1], \
        "milp-(64)-(4)-{}".format(i)
    
    # (89)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([delta_dc[(t, 1, d1, d2, k)] 
            for d1 in range(MAX_VAL + 1) 
            for d2 in range(MAX_VAL + 1)
            for k in range(4)]) == 1, "milp-(65)-(1)-{}".format(t)
    
    # (90)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([k * delta_dc[(t, 1, d1, d2, k)] 
            for d1 in range(MAX_VAL + 1)
            for d2 in range(MAX_VAL + 1) 
            for k in range(4)]) == beta_tilde[(t, 1)], \
        "milp-(65)-(2)-{}".format(t)
    
    # (91)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([d1 * delta_dc[(t, 1, d1, d2, k)] 
            for d1 in range(1, MAX_VAL + 1)
            for d2 in range(MAX_VAL + 1) 
            for k in range(4)]) == deg[s_star + t - 1, 1], \
        "milp-(65)-(3)-{}".format(t)
    
    # (92)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([d2 * delta_dc[(t, 1, d1, d2, k)] 
            for d1 in range(MAX_VAL + 1)
            for d2 in range(1, MAX_VAL + 1) for k in range(4)]) == \
        deg[s_star + t, 1], "milp-(65)-(4)-{}".format(t)
    
    
    # (93)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([delta_dc[(p, i, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1)
                for k in range(4)]) == 1, \
            "milp-(66)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([delta_dc[(p, i, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1)
                for k in range(4)]) == 1, \
            "milp-(66)-{}-{}".format(p, i)
    
    # (94)
    for s in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([k * delta_dc[(s, i, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1) 
                for d2 in range(MAX_VAL + 1)
                for k in range(4)]) == beta_tilde[(s, i)],\
            "milp-(67)-{}-{}".format(s, i)
    
    # (95)
    for t in range(1, t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([k * delta_dc[(s_star + t, 
                i, d1, d2, k)] for d1 in range(MAX_VAL + 1)
                                for d2 in range(MAX_VAL + 1) 
                                for k in range(4)]) == \
            beta_tilde[(s_star + t, i)], \
            "milp-(68)-{}-{}".format(t, i)
    
    # (96)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([d1 * delta_dc[(p, i, d1, d2, k)] 
                for d1 in range(1, MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1) for k in range(4)]) == \
                deg[(p, prt_S[i])], "milp-(69)-(1)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([d1 * delta_dc[(p, i, d1, d2, k)] 
                for d1 in range(1, MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1) for k in range(4)]) == \
                deg[(p, prt_T[i])], "milp-(69)-(1)-{}-{}".format(p, i)
    
    # (97)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([d2 * delta_dc[(p, i, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(1, MAX_VAL + 1) 
                for k in range(4)]) == \
                    deg[(p, i)], "milp-(69)-(2)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([d2 * delta_dc[(p, i, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(1, MAX_VAL + 1) 
                for k in range(4)]) == \
                    deg[(p, i)], "milp-(69)-(2)-{}-{}".format(p, i)
    
    # (98)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([delta_dc_hat[(s, t, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1) 
                for k in range(4)]) == 1, \
            "milp-(70)-(1)-{}-{}".format(s, t)
    
    # (99)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([k * delta_dc_hat[(s, t, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1) 
                for k in range(4)]) == \
                    beta_hat[(s, t)], "milp-(70)-(2)-{}-{}".format(s, t)

    # (100)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([d1 * delta_dc_hat[(s, t, d1, d2, k)] 
                for d1 in range(1, MAX_VAL + 1)
                for d2 in range(MAX_VAL + 1) for k in range(4)]) == \
            deg[(s, 1)], "milp-(70)-(3)-{}-{}".format(s, t)

    # (101)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([d2 * delta_dc_hat[(s, t, d1, d2, k)] 
                for d1 in range(MAX_VAL + 1)
                for d2 in range(1, MAX_VAL + 1) 
                for k in range(4)]) == \
                deg[(s_star + t, 1)], "milp-(70)-(4)-{}-{}".format(s, t)

    # modified 20/07/28
    for (d1, d2, k) in Bcset:
        if d1 != d2:
            MILP += pulp.lpSum([delta_dc[(i, d1, d2, k)] 
                        for i in E_B]) + \
                    pulp.lpSum([delta_dc[(i, d2, d1, k)] 
                        for i in E_B]) + \
                    pulp.lpSum([delta_dc[(t, 1, d1, d2, k)] 
                        for t in range(2, t_star + 1)]) + \
                    pulp.lpSum([delta_dc[(t, 1, d2, d1, k)] 
                        for t in range(2, t_star + 1)]) + \
                    pulp.lpSum([delta_dc_hat[(s, t, d1, d2, k)] 
                        for s in range(1, s_star + 1)
                        for t in range(1, t_star + 1)]) + \
                    pulp.lpSum([delta_dc_hat[(s, t, d2, d1, k)] 
                        for s in range(1, s_star + 1)
                        for t in range(1, t_star + 1)]) == \
                    bc_in[(d1, d2, k)], \
                    "milp-(71)-(1)-{}-{}-{}".format(d1, d2, k)

            MILP += pulp.lpSum([delta_dc[(p, i, d1, d2, k)] 
                        for p in range(1, s_star + 1)
                        for i in range(2, n_S_tree + 1)]) + \
                    pulp.lpSum([delta_dc[(p, i, d1, d2, k)] 
                        for p in range(s_star + 1, s_star + t_star + 1)
                        for i in range(2, n_T_tree + 1)]) + \
                    pulp.lpSum([delta_dc[(p, i, d2, d1, k)] 
                        for p in range(1, s_star + 1)
                        for i in range(2, n_S_tree + 1)]) + \
                    pulp.lpSum([delta_dc[(p, i, d2, d1, k)] 
                        for p in range(s_star + 1, s_star + t_star + 1)
                        for i in range(2, n_T_tree + 1)]) == \
                    bc_ex[(d1, d2, k)], \
                    "milp-(71)-(2)-{}-{}-{}".format(d1, d2, k)

        elif d1 == d2:
            MILP += pulp.lpSum([delta_dc[(i, d1, d2, k)] 
                        for i in E_B]) + \
                    pulp.lpSum([delta_dc[(t, 1, d1, d2, k)] 
                        for t in range(2, t_star + 1)]) + \
                    pulp.lpSum([delta_dc_hat[(s, t, d1, d2, k)] 
                        for s in range(1, s_star + 1)
                        for t in range(1, t_star + 1)]) == \
                    bc_in[(d1, d2, k)], \
                    "milp-(72)-(1)-{}-{}-{}".format(d1, d2, k)

            MILP += pulp.lpSum([delta_dc[(p, i, d1, d2, k)] 
                        for p in range(1, s_star + 1)
                        for i in range(2, n_S_tree + 1)]) + \
                    pulp.lpSum([delta_dc[(p, i, d1, d2, k)] 
                        for p in range(s_star + 1, s_star + t_star + 1)
                        for i in range(2, n_T_tree + 1)]) == \
                    bc_ex[(d1, d2, k)], \
                    "milp-(72)-(2)-{}-{}-{}".format(d1, d2, k)

    return MILP

def add_constraints_B_8(
    MILP,
    n_star, dia_star, k_star, d_max, bl_star, bh_star,
    s_star, t_star, c_star, root_num, E_B, E1_plus, E1_minus, tail, head,
    n_S_tree, n_T_tree, E_B_left, E_B_right, V_B_left, V_B_right,
    V_B_bh_star, V_B_s, E_B_s, Cld_S, Cld_T, Cld_B, prt_S, prt_T, Dsn_S,  
    P_S_prc, P_T_prc, Bcset,
    MAX_BOND, MAX_VAL, MAX_CODE,
    Lambda, Extra_Lambda, Code_Lambda, val, mass, 
    Gamma_equal, Gamma_less, Gamma_great, Gamma_zero, Gamma_ALL,
    a_3, a_1_3, Ncset, Ncset_e, ce_nc, dg_nc, ml,
    a, e_st, e_ts, chi, delta_clr, clr, degbt_plus, degbt_minus, 
    tilde, u, v, e, alpha_tilde, beta_tilde, beta_hat,
    delta_alpha, delta_beta_tilde, delta_beta_hat,
    ce_bt, ce_ft, ce, Mass, b_bt, b_ft, b, deg, delta_deg, dg,
    delta_tau, delta_tau_hat, ac_bt, ac_ft, ac, n_H,
    dc, delta_dc, delta_dc_hat,
    nc, delta_nc
    ):
    ''' A function to add constraints in section B.8 '''
    # (74)
    for p in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            MILP += pulp.lpSum([delta_nc[(p, i, nu)] 
                        for nu in Ncset_e]) == 1, \
                    "milp-(112)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(1, n_T_tree + 1):
            MILP += pulp.lpSum([delta_nc[(p, i, nu)] 
                        for nu in Ncset_e]) == 1, \
                    "milp-(112)-{}-{}".format(p, i)

    # (75)
    for p in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            MILP += pulp.lpSum([ce_nc[nu] * delta_nc[(p, i, nu)] 
                        for nu in Ncset_e]) == \
                    alpha_tilde[(p, i)], \
                    "milp-(113)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(1, n_T_tree + 1):
            MILP += pulp.lpSum([ce_nc[nu] * delta_nc[(p, i, nu)] 
                        for nu in Ncset_e]) == \
                    alpha_tilde[(p, i)], \
                    "milp-(113)-{}-{}".format(p, i)

    # (76)
    for p in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            MILP += pulp.lpSum([dg_nc[nu] * delta_nc[(p, i, nu)] 
                        for nu in Ncset_e]) == \
                    deg[(p, i)], "milp-(114)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(1, n_T_tree + 1):
            MILP += pulp.lpSum([dg_nc[nu] * delta_nc[(p, i, nu)] 
                        for nu in Ncset_e]) == \
                    deg[(p, i)], "milp-(114)-{}-{}".format(p, i)

    # (77)
    for s in range(1, s_star + 1):
        MILP += pulp.lpSum([beta_tilde[i] 
                    for i in (E1_plus[s] | E1_minus[s])]) + \
                pulp.lpSum([beta_hat[(s, t)] 
                    for t in range(1, t_star + 1)]) + \
                pulp.lpSum([beta_tilde[(s, j)] 
                    for j in Cld_S[1]]) == \
                pulp.lpSum([ml[nu] * delta_nc[(s, 1, nu)] 
                    for nu in Ncset_e]), "milp-(115)-{}".format(s)

    # (78)
    for t in range(1, t_star + 1):
        MILP += pulp.lpSum([beta_hat[(s, t)] 
                    for s in range(1, s_star + 1)]) + \
                beta_tilde[(t, 1)] + beta_tilde[(t + 1, 1)] + \
                pulp.lpSum([beta_tilde[(s_star + t, j)] 
                    for j in Cld_T[1]]) == \
                pulp.lpSum([ml[nu] * 
                    delta_nc[(s_star + t, 1, nu)] 
                    for nu in Ncset_e]), "milp-(116)-{}".format(t)

    # (79)
    for s in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += beta_tilde[(s, i)] +\
                    pulp.lpSum([beta_tilde[(s, j)] 
                        for j in Cld_S[i]]) == \
                    pulp.lpSum([ml[nu] * delta_nc[(s, i, nu)] 
                        for nu in Ncset_e]), \
                    "milp-(117)-{}-{}".format(s, i)
    
    # (80)
    for t in range(1, t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += beta_tilde[(s_star + t, i)] + \
                    pulp.lpSum([beta_tilde[(s_star + t, j)] 
                        for j in Cld_T[i]]) == \
                    pulp.lpSum([ml[nu] * 
                        delta_nc[(s_star + t, i, nu)] 
                        for nu in Ncset_e]),\
                        "milp-(118)-{}-{}".format(t, i)

    # (81)
    for i in (E_B):
        MILP += pulp.lpSum([delta_nc[(tail[i], 1, nu)]
                    for nu in (a_3 | a_1_3)]) >= \
                beta_tilde[i] - 2, "milp-(119)-{}".format(i)

    # (82)
    for i in (E_B):
        MILP += pulp.lpSum([delta_nc[(head[i], 1, nu)]
                    for nu in (a_3 | a_1_3)]) >= \
                beta_tilde[i] - 2, "milp-(120)-{}".format(i)

    # (83)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([delta_nc[(s_star + t - 1, 1, nu)]
                    for nu in (a_3 | a_1_3)]) >= \
                beta_tilde[(t, 1)] - 2, "milp-(121)-{}".format(t)

    # (84)
    for t in range(2, t_star + 1):
        MILP += pulp.lpSum([delta_nc[(s_star + t, 1, nu)]
                    for nu in (a_3 | a_1_3)]) >= \
                beta_tilde[(t, 1)] - 2, "milp-(122)-{}".format(t)

    # (85)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([delta_nc[(p, prt_S[i], nu)]
                        for nu in (a_3 | a_1_3)]) >= \
                    beta_tilde[(p, i)] - 2, \
                    "milp-(123)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([delta_nc[(p, prt_T[i], nu)]
                        for nu in (a_3 | a_1_3)]) >= \
                    beta_tilde[(p, i)] - 2, \
                    "milp-(123)-{}-{}".format(p, i)

    # (86)
    for p in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            MILP += pulp.lpSum([delta_nc[(p, i, nu)]
                        for nu in (a_3 | a_1_3)]) >= \
                    beta_tilde[(p, i)] - 2, \
                    "milp-(124)-{}-{}".format(p, i)
    for p in range(s_star + 1, s_star + t_star + 1):
        for i in range(2, n_T_tree + 1):
            MILP += pulp.lpSum([delta_nc[(p, i, nu)]
                        for nu in (a_3 | a_1_3)]) >= \
                    beta_tilde[(p, i)] - 2, \
                    "milp-(124)-{}-{}".format(p, i)

    # (87)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([delta_nc[(s, 1, nu)]
                        for nu in (a_3 | a_1_3)]) >= \
                    beta_hat[(s, t)] - 2, \
                    "milp-(125)-{}-{}".format(s, t)

    # (88)
    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            MILP += pulp.lpSum([delta_nc[(s_star + t, 1, nu)]
                        for nu in (a_3 | a_1_3)]) >= \
                    beta_hat[(s, t)] - 2, \
                    "milp-(126)-{}-{}".format(s, t)

    # (89)
    for nu in Ncset:
        MILP += pulp.lpSum([delta_nc[(p, i, nu)]
                    for p in range(1, s_star + 1)
                    for i in range(1, n_S_tree + 1)]) + \
                pulp.lpSum([delta_nc[(p, i, nu)]
                    for p in range(s_star + 1, s_star + t_star + 1)
                    for i in range(1, n_T_tree + 1)]) == \
                nc[nu], "milp-(127)-{}".format(nu)

    return MILP


def add_constraints_ANN(MILP,
                        descriptors,
                        num_fv,
                        y, n_star):
    ''' A function to add constraints used in ANN '''
    for i in range(1, num_fv):
        if i == 2:  # Mass
            MILP += y[(1, i)] == descriptors[i] * (1 / n_star),\
            "ann_input_{}".format(i)
        else:
            MILP += y[(1, i)] == descriptors[i], \
            "ann_input_{}".format(i)

    return MILP

def add_constraints_AD(MILP,
                       n_star,
                       num_fv, ceset, acset,
                       ce_in, ce_ex, ac_in, ac_ex, y,
                       Lambda,Gamma_equal, Gamma_less,
                       desc_index_v, desc_index_e,
                       ann_training_data_filename,
                       ad_min_filename,
                       ad_max_filename):
    ''' A function to add constraints for AD '''
    fv_raw = pd.read_csv(ann_training_data_filename)
    fv = fv_raw.values

    for atom in Lambda:
        if atom not in ceset:
            MILP += ce_in[atom] == 0, "ad-ce-{}-in".format(atom)
            MILP += ce_ex[atom] == 0, "ad-ce-{}-ex".format(atom)
    for gamma in (Gamma_equal | Gamma_less):
        if gamma not in acset:
            MILP += ac_in[gamma] == 0, "ad-ac-{}-in".format(gamma)
            MILP += ac_ex[gamma] == 0, "ad-ac-{}-ex".format(gamma)

    domain_lower = [float('inf') for _ in range(fv.shape[1] + 1)]
    domain_upper = [-float('inf') for _ in range(fv.shape[1] + 1)]

    for i in range(fv.shape[0]):
        for j in range(fv.shape[1]):
            if j in desc_index_v:
                temp = fv[i][j] / fv[i][1]
                if temp < domain_lower[j]:
                    domain_lower[j] = temp
                elif temp > domain_upper[j]:
                    domain_upper[j] = temp
            elif j in desc_index_e:
                temp = fv[i][j] / (fv[i][1] - 1)
                if temp < domain_lower[j]:
                    domain_lower[j] = temp
                elif temp > domain_upper[j]:
                    domain_upper[j] = temp
            else:
                temp = fv[i][j]
                if temp < domain_lower[j]:
                    domain_lower[j] = temp
                elif temp > domain_upper[j]:
                    domain_upper[j] = temp


    fvmin = pd.read_csv(ad_min_filename, index_col="n")
    fvmax = pd.read_csv(ad_max_filename, index_col="n")

    i = 0
    for des in fv_raw.columns:
        if des in fvmin.columns:
            if i in desc_index_v:
                if fvmin.at[n_star, des] / n_star < domain_lower[i]:
                    domain_lower[i] = fvmin.at[n_star, des] / n_star
                if fvmax.at[n_star, des] / n_star > domain_lower[i]:
                    domain_upper[i] = fvmax.at[n_star, des] / n_star
            elif i in desc_index_e:
                if fvmin.at[n_star, des] / (n_star - 1) < domain_lower[i]:
                    domain_lower[i] = fvmin.at[n_star, des] / (n_star - 1)
                if fvmax.at[n_star, des] / (n_star - 1) > domain_lower[i]:
                    domain_upper[i] = fvmax.at[n_star, des] / (n_star - 1)
        i += 1

    for i in range(1, num_fv):
        if i in desc_index_v:
            MILP += y[(1, i)] * (1 / n_star) >= domain_lower[i], \
                "ad-{}-lb".format(i)
            MILP += y[(1, i)] * (1 / n_star) <= domain_upper[i], \
                "ad-{}-ub".format(i)
        elif i in desc_index_e:
            MILP += y[(1, i)] * (1 / (n_star - 1)) >= \
                domain_lower[i], "ad-{}-lb".format(i)
            MILP += y[(1, i)] * (1 / (n_star - 1)) <= \
                domain_upper[i], "ad-{}-ub".format(i)

    return MILP

def print_gstar_file(
    dia_star, d_max,
    Bcset,
    set_Lambda,
    Gamma_equal, Gamma_less,
    ce_in, ce_ex, ac_in, ac_ex, 
    bc_in, bc_ex, dg_in, dg_ex,
    outputfilename):

    ''' A function to output input file for g* program '''

    with open(outputfilename, "w") as f:
        f.write(str(round(len(set_Lambda))) + "\n")
        for atom in set_Lambda:
            f.write("{} {} {} {} {}\n".format(atom.symbol, 
                atom.mass, atom.valence, 
                str(round(ce_in[atom.symbol].value())),
                str(round(ce_ex[atom.symbol].value()))))

        f.write(str(len(Gamma_equal) + len(Gamma_less)) + "\n")
        for (atom1, atom2, k) in (Gamma_less | Gamma_equal):
            f.write("{} {} {} {} {}\n".format(atom1, 
                atom2, k, 
                str(round(ac_in[(atom1, atom2, k)].value())),
                str(round(ac_ex[(atom1, atom2, k)].value()))))

        for i in range(1, 5):
            f.write("{} {}\n".format(
                str(round(dg_in[i].value())),
                str(round(dg_ex[i].value()))))

        f.write(str(dia_star) + "\n")
        f.write(str(d_max) + "\n")

        for (d1, d2, k) in Bcset:
            f.write("{} {} {} {} {}\n".format(d1, d2, k, 
                str(round(bc_in[(d1, d2, k)].value())),
                str(round(bc_ex[(d1, d2, k)].value()))))

        f.write("\n")
        f.close()

    return

def print_sdf_file(
    n_star, E_B, tail, head,
    s_star, t_star, n_S_tree, n_T_tree,
    prt_S, prt_T,
    Lambda,
    u, v, alpha_tilde, beta_tilde, beta_hat,
    outputfilename
    ):

    ''' A function to output sdf file.  '''

    graph_col = {i : " " for i in range(1, n_star + 1)}
    graph_adj = {(i, j) : 0 for i in range(1, n_star + 1) for j in range(1, n_star + 1)}

    graph_index = 0

    index_u = {(s, i) : 0 for s in range(1, s_star + 1) for i in range(1, n_S_tree + 1)}
    index_v = {(t, i) : 0 for t in range(1, t_star + 1) for i in range(1, n_T_tree + 1)}

    for s in range(1, s_star + 1):
        for i in range(1, n_S_tree + 1):
            if round(u[(s, i)].value()) != 0:
                graph_index += 1
                index_u[(s, i)] = graph_index
                graph_col[graph_index] = Lambda[round(alpha_tilde[(s, i)].value()) - 1]

    for t in range(1, t_star + 1):
        for i in range(1, n_T_tree + 1):
            if round(v[(t, i)].value()) != 0:
                graph_index += 1
                index_v[(t, i)] = graph_index
                graph_col[graph_index] = Lambda[round(alpha_tilde[(s_star + t, i)].value()) - 1]

    for i in E_B:
        if round(beta_tilde[i].value()) != 0:
            ind1 = index_u[(tail[i], 1)]
            ind2 = index_u[(head[i], 1)]
            mul = round(beta_tilde[i].value())
            graph_adj[(ind1, ind2)] = mul
            graph_adj[(ind2, ind1)] = mul

    for s in range(1, s_star + 1):
        for i in range(2, n_S_tree + 1):
            if round(beta_tilde[(s, i)].value()) != 0:
                ind1 = index_u[(s, i)]
                ind2 = index_u[(s, prt_S[i])]
                mul = round(beta_tilde[(s, i)].value())
                graph_adj[(ind1, ind2)] = mul
                graph_adj[(ind2, ind1)] = mul

    for t in range(1, t_star + 1):
        for i in range(2, n_T_tree + 1):
            if round(beta_tilde[(s_star + t, i)].value()) != 0:
                ind1 = index_v[(t, i)]
                ind2 = index_v[(t, prt_T[i])]
                mul = round(beta_tilde[(s_star + t, i)].value())
                graph_adj[(ind1, ind2)] = mul
                graph_adj[(ind2, ind1)] = mul

    for t in range(1, t_star + 2):
        if round(beta_tilde[(t, 1)].value()) != 0:
            ind1 = index_v[(t, 1)]
            ind2 = index_v[(t - 1, 1)]
            mul = round(beta_tilde[(t, 1)].value())
            graph_adj[(ind1, ind2)] = mul
            graph_adj[(ind2, ind1)] = mul

    for s in range(1, s_star + 1):
        for t in range(1, t_star + 1):
            if round(beta_hat[(s, t)].value()) != 0:
                ind1 = index_u[(s, 1)]
                ind2 = index_v[(t, 1)]
                mul = round(beta_hat[(s, t)].value())
                graph_adj[(ind1, ind2)] = mul
                graph_adj[(ind2, ind1)] = mul

    # for i in range(1, n_star+1):
    #     for j in range(1, n_star+1):
    #         if graph_adj[(i, j)] != 0:
    #             print("{}_{}:{}".format(i, j, graph_adj[(i, j)]))

    with open(outputfilename, "w") as f:
        f.write("1\n")
        f.write("MILP\n")
        f.write("acyclic graph\n")
        f.write("{:3}{:3}  0  0  0  0  0  0  0  0999 V2000 \n".format(n_star, n_star - 1))
        
        for i in range(1, n_star + 1):
            f.write("    0.0000    0.0000    0.0000 {:2}  0  0  0  0  0  0  0  0  0  0  0  0\n".format(graph_col[i]))

        for i in range(1, n_star + 1):
            for j in range(i + 1, n_star + 1):
                if graph_adj[(i, j)] > 0:
                    f.write("{:3}{:3}{:3}  0  0  0  0\n".format(i, j, graph_adj[(i, j)]))

        f.write("M  END\n")
        f.write("$$$$\n")
        f.close()

    return
