import pulp
import time
from collections import namedtuple
import sys
import pandas as pd

CG_element = namedtuple("CG_element", "symbol, valence, mass")
MAX_BOND = 3
MAX_VAL = 4

class chemicalRootedTree():
    def __init__(self):
        self.root = ("e", 0)
        self.index = 0
        self.vertex = []
        self.adj = []
        self.alpha = []
        self.beta = []
        self.height = 0

##############################
#     Table of Variables     #
##############################
"""
# -------- Global ---------
    Gamma_co : \Gamma^{\rm co}
    Gamma_in : \Gamma^{\rm in}
    Gamma_ex : \Gamma^{\rm ex}
    Gamma_co_less : \Gamma^{\rm co}_{<}
    Gamma_co_equal : \Gamma^{\rm co}_{=}
    Gamma_co_great : \Gamma^{\rm co}_{>}
    Gamma_co_ac : \Gamma^{\rm co}_{\rm ac}
    Gamma_in_ac : \Gamma^{\rm in}_{\rm ac}
    Gamma_ex_ac : \Gamma^{\rm ex}_{\rm ac}
    Gamma_co_ac_less : \Gamma^{\rm co}_{\rm ac, <}
    Gamma_co_ac_equal : \Gamma^{\rm co}_{\rm ac, =}
    Gamma_co_ac_great : \Gamma^{\rm co}_{\rm ac, >}
    Gamma_tilde_ac_C : \tilde{\Gamma_{\rm ac}^{\rm C}}
    Gamma_tilde_ac_T : \tilde{\Gamma_{\rm ac}^{\rm T}}
    Gamma_tilde_ac_CT : \tilde{\Gamma_{\rm ac}^{\rm CT}}
    Gamma_tilde_ac_TC : \tilde{\Gamma_{\rm ac}^{\rm TC}}
    Gamma_tilde_ac_F : \tilde{\Gamma_{\rm ac}^{\rm F}}
    Gamma_tilde_ac_CF : \tilde{\Gamma_{\rm ac}^{\rm CF}}
    Gamma_tilde_ac_TF : \tilde{\Gamma_{\rm ac}^{\rm TF}}
    Gamma_tilde_ac_X : \tilde{\Gamma_{\rm ac}^{\rm X}}
    Lambda : \Lambda
    Lambda_dg_co : \Lambda_{\rm dg}^{\rm co}
    Lambda_dg_nc : \Lambda_{\rm dg}^{\rm nc}
    Lambda_tilde_dg_C : \tilde{\Lambda_{\rm dg}^{\rm C}}
    Lambda_tilde_dg_T : \tilde{\Lambda_{\rm dg}^{\rm C}}
    Lambda_tilde_dg_F : \tilde{\Lambda_{\rm dg}^{\rm C}}
    Gamma_tilde_ec_C : \tilde{\Gamma_{\rm ec}^{\rm C}}
    Gamma_tilde_ec_T : \tilde{\Gamma_{\rm ec}^{\rm T}}
    Gamma_tilde_ec_CT : \tilde{\Gamma_{\rm ec}^{\rm CT}}
    Gamma_tilde_ec_TC : \tilde{\Gamma_{\rm ec}^{\rm TC}}
    Gamma_tilde_ec_F : \tilde{\Gamma_{\rm ec}^{\rm F}}
    Gamma_tilde_ec_CF : \tilde{\Gamma_{\rm ec}^{\rm CF}}
    Gamma_tilde_ec_TF : \tilde{\Gamma_{\rm ec}^{\rm TF}}
    Gamma_tilde_ec_X : \tilde{\Gamma_{\rm ec}^{\rm X}}
# -------- A.1 Selecting Core-vertices and Core-edges   --------
    e_C : e^{\rm C}
    v_T : v^{^rm T}
    e_T : e^{\rm T}
    chi_T : \chi^{\rm T}
    clr_T : {\rm clr}^{\rm T}
    delta_chi_T : \delta_\chi^{\rm T}
    deg_tilde_C_plus : \tilde{{\rm deg}_{\rm C}}^{+}
    deg_tilde_C_minus : \tilde{{\rm deg}_{\rm C}}^{-}
# -------- A.2 Constraints for Including Internal Vertices and Edges --------
    bl_G : {\rm bl}^{\rm G} 
    v_F : v^{\rm F}
    e_F : e^{\rm F}
    chi_F : \chi^{\rm F}
    clr_F : {\rm clr}^{\rm F}
    delta_chi_F :\delta_\chi^{\rm F}
    sigma : \sigma
    sigma_C : \sigma^{\rm C}
    sigma_T : \sigma^{\rm T}
    sigma_F : \sigma^{\rm F}
    bl : {\rm bl}
# -------- A.3 Constraints for Including Fringe-trees --------
    n_G : n_G
    ch_G : {\rm cg}_G
    v_T : v^{\rm T}
    v_C : v^{\rm C}
    v_F : v^{\rm F}
tilde_C    h_T : h^{\rm T}
    h_C : h^{\rm C}
    h_F : h^{\rm F}
    sigma : \sigma
# -------- A.4 Descriptor for the Number of Specified Degree --------
    deg_C : {\rm deg}^{\rm C}
    deg_T : {\rm deg}^{\rm T}
    deg_F : {\rm deg}^{\rm F}
    deg_CT : {\rm deg}^{\rm CT}
    deg_TC : {\rm deg}^{\rm TC}
    delta_dg_C : \delta_{\rm dg}^{\rm C}
    delta_dg_T : \delta_{\rm dg}^{\rm T}
    delta_dg_F : \delta_{\rm dg}^{\rm F}
    dg : {\rm dg}
    dg_co : {\rm dg}^{\rm co}
    dg_C : {\rm dg}_{\rm C}
    dg_T : {\rm dg}_{\rm T}
    dg_nc : {\rm dg}^{\rm nc}
    dg_in : {\rm dg}^{\rm in}
    dg_ex : {\rm dg}^{\rm ex}
    dg_Cp : {\rm dg}_{\rm C(p)}
    dg_Tp : {\rm dg}_{\rm T(p)}
    dg_Fp : {\rm dg}_{\rm F(p)}
# -------- A.5 Assigning Multiplicity --------
    beta_T : \beta^{\rm T}
    beta_F : \beta^{\rm F}
    beta_C : \beta^{\rm C}
    beta_plus : \beta^+
    beta_minus : \beta^-
    beta_in : \beta^{\rm in}
    delta_beta_T : \delta_\beta^{\rm T}
    delta_beta_F : \delta_\beta^{\rm F}
    delta_beta_C : \delta_\beta^{\rm C}
    delta_beta_plus : \delta_\beta^+
    delta_beta_minus : \delta_\beta^-
    delta_beta_in : \delta_\beta^{\rm in}
    bd : {\rm bd}
    bd_co : {\rm bd}^{\rm co}
    bd_in : {\rm bd}^{\rm in}
    bd_ex : {\rm bd}^{\rm ex}
    bd_C : {\rm bd}_{\rm C}
    bd_T : {\rm bd}_{\rm T}
    bd_CT : {\rm bd}_{\rm CT}
    bd_TC : {\rm bd}_{\rm TC}
    bd_F : {\rm bd}_{\rm F}
    bd_CF : {\rm bd}_{\rm CF}
    bd_TF : {\rm bd}_{\rm TF}
    bd_Cp : {\rm bd}_{\rm C(p)}
    bd_Tp : {\rm bd}_{\rm T(p)}
    bd_Fp : {\rm bd}_{\rm F(@)}
# -------- A.6 Assigning Chemical Elements and Valence Condition --------
    beta_CT : \beta^{\rm CT}
    beta_TC : \beta^{\rm TC}
    beta_CF : \beta^{\rm CF}
    beta_TF : \beta^{\rm TF}
    alpha_T : \alpha^{\rm T}
    alpha_F : \alpha^{\rm F}
    alpha_C : \alpha^{\rm C}
    delta_alpha_T : \delta_\alpha^{\rm T}
    delta_alpha_F : \delta_\alpha^{\rm F}
    delta_alpha_C : \delta_\alpha^{\rm C}
    MASS : {\rm Mass}
    n_H : n_{\rm H}
    na : {\rm na}
    na_co : {\rm na}^{\rm co}
    na_C : {\rm na}_{\rm C}
    na_T : {\rm na}_{\rm T}
    na_nc : {\rm na}^{\rm nc}
    na_in : {\rm na}^{\rm in}
    na_Cp : {\rm na}_{\rm C(p)}
    na_Tp : {\rm na}_{\rm T(p)}
    na_Fp : {\rm na}_{\rm F(p)}
# -------- A.7 Constraints for Bounds on the Number of Bonds --------
    bd_T : {\rm bd}_{\rm T}
# -------- A.8 Descriptors for the Number of Adjacency-configuration --------
    delta_alpha_CT : \delta_\alpha^{\rm CT}
    delta_alpha_TC : \delta_\alpha^{\rm TC}
    delta_alpha_TF : \delta_\alpha^{\rm TF}
    delta_alpha_CF : \delta_\alpha^{\rm CF}
    ac_C : {\rm ac}_{\rm C}
    ac_T : {\rm ac}_{\rm T}
    ac_F : {\rm ac}_{\rm F}
    ac_Cp : {\rm ac}_{\rm C(p)}
    ac_Tp : {\rm ac}_{\rm T(p)}
    ac_Fp : {\rm ac}_{\rm F(p)}
    ac_CT : {\rm ac}_{\rm CT}
    ac_TC : {\rm ac}_{\rm TC}
    ac_CF : {\rm ac}_{\rm CF}
    ac_TF : {\rm ac}_{\rm TF}
    delta_ac_C: \delta_{\rm ac}^{\rm C}
    delta_ac_T: \delta_{\rm ac}^{\rm T}
    delta_ac_F: \delta_{\rm ac}^{\rm F}
    delta_ac_Cp: \delta_{\rm ac}^{\rm C(p)}
    delta_ac_Tp: \delta_{\rm ac}^{\rm T(p)}
    delta_ac_Fp: \delta_{\rm ac}^{\rm F(p)}
    delta_ac_CT: \delta_{\rm ac}^{\rm CT}
    delta_ac_TC: \delta_{\rm ac}^{\rm TC}
    delta_ac_CF: \delta_{\rm ac}^{\rm CF}
    delta_ac_TF: \delta_{\rm ac}^{\rm TF}
    ac_co : {\rm ac}^{\rm co}
    ac_in : {\rm ac}^{\rm in}
    ac_ex : {\rm ac}^{\rm ex}
    ac_nc : {\rm ac}^{\rm nc}
# -------- A.9 Descriptor for the Number of Chemical Symbols --------
    ns_co : {\rm ns}^{\rm co}
    ns_nc : {\rm ns}^{\rm nc}
    delta_ns_C : \delta_{\rm ns}^{\rm C}
    delta_ns_T : \delta_{\rm ns}^{\rm T}
    delta_ns_F : \delta_{\rm ns}^{\rm F}
# -------- A.10 Descriptor for the Number of Edge-configurations --------
    delta_dg_CT : \delta_{\rm dg}^{\rm CT}
    delta_dg_TC : \delta_{\rm dg}^{\rm TC}
    delta_dg_CF : \delta_{\rm dg}^{\rm CF}
    delta_dg_TF : \delta_{\rm dg}^{\rm TF}
    ec_C : {\rm ec}_{\rm C}
    ec_T : {\rm ec}_{\rm T}
    ec_F : {\rm ec}_{\rm F}
    ec_Cp : {\rm ec}_{\rm C(p)}
    ec_Tp : {\rm ec}_{\rm T(p)}
    ec_Fp : {\rm ec}_{\rm F(p)}
    ec_CT : {\rm ec}_{\rm CT}
    ec_TC : {\rm ec}_{\rm TC}
    ec_CF : {\rm ec}_{\rm CF}
    ec_TF : {\rm ec}_{\rm TF}
    delta_ec_C: \delta_{\rm ec}^{\rm C}
    delta_ec_T: \delta_{\rm ec}^{\rm T}
    delta_ec_F: \delta_{\rm ec}^{\rm F}
    delta_ec_Cp: \delta_{\rm ec}^{\rm C(p)}
    delta_ec_Tp: \delta_{\rm ec}^{\rm T(p)}
    delta_ec_Fp: \delta_{\rm ec}^{\rm F(p)}
    delta_ec_CT: \delta_{\rm ec}^{\rm CT}
    delta_ec_TC: \delta_{\rm ec}^{\rm TC}
    delta_ec_CF: \delta_{\rm ec}^{\rm CF}
    delta_ec_TF: \delta_{\rm ec}^{\rm TF}
    ec_co : {\rm ec}^{\rm co}
    ec_in : {\rm ec}^{\rm in}
    ec_ex : {\rm ec}^{\rm ex}
    ec_nc : {\rm ec}^{\rm nc}
     
"""

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

def prepare_Lambda(
    set_Lambda,
    Lambda
):
    new_set_Lambda = list()
    for atom in set_Lambda:
        if atom.symbol in Lambda:
            new_set_Lambda.append(atom)

    return new_set_Lambda

def prepare_Lambda_dg(
    set_Lambda,
    Lambda,
    Lambda_int,
    Lambda_ex,
    Gamma_int,
    Gamma_int_ac
):
    #   A function to prepare the set Lambda_dg and Gamma
    set_Lambda_dg = list()
    for atom in set_Lambda:
        for i in range(1, atom.valence + 1):
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

    return set_Lambda_dg, Lambda_dg, epsilon, epsilon_dg, Code_Lambda_dg, \
        Code_Lambda, Code_Lambda_int, Code_Lambda_ex, \
        MAX_CODE, MAX_CODE_dg, MAX_CODE_int, MAX_CODE_ex, \
        Code_Gamma_ec_int, Code_Gamma_ac_int, val, mass

# A function to prepare sets used in A.8 and A.10
def prepare_Gamma_ac(
    Gamma_int,
    Gamma_int_ac,
    Code_Lambda,
    Code_Lambda_dg
):
    Gamma_int_less = list()
    Gamma_int_equal = list()
    Gamma_int_great = list()

    for ((a1, d1), (a2, d2), m) in Gamma_int:
        if Code_Lambda_dg[(a1, d1)] < Code_Lambda_dg[(a2, d2)]:
            Gamma_int_less.append(((a1, d1), (a2, d2), m))
        elif Code_Lambda_dg[(a1, d1)] == Code_Lambda_dg[(a2, d2)]:
            Gamma_int_equal.append(((a1, d1), (a2, d2), m))
        elif Code_Lambda_dg[(a1, d1)] > Code_Lambda_dg[(a2, d2)]:
            Gamma_int_great.append(((a1, d1), (a2, d2), m))

    Gamma_int_ac_less = list()
    Gamma_int_ac_equal = list()
    Gamma_int_ac_great = list()
 
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

    return Gamma_int_less, Gamma_int_equal, Gamma_int_great, \
        Gamma_int_ac_less, Gamma_int_ac_equal, Gamma_int_ac_great, \
        Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, \
        Gamma_tilde_ac_F, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, \
        Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, \
        Gamma_tilde_ec_F, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF

# A function to prepare sets used in A.9
def prepare_Lambda_tilde(
    Lambda_dg_int
):

    Lambda_tilde_dg_C = Lambda_dg_int[:]
    Lambda_tilde_dg_T = Lambda_dg_int[:]
    Lambda_tilde_dg_F = Lambda_dg_int[:]

    return Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F

def prepare_F_Lambda(
    Lambda
):
    F_Lambda = list()
    for ele in Lambda:
        f_tmp = chemicalRootedTree()
        f_tmp.root = (ele[0], 0)
        F_Lambda.append(f_tmp)

    return F_Lambda

# Prepare Constants of Seed Graph
def prepare_dataset_for_scheme_graph(
    n_star,
    rho,
    V_C,
    E_C,
    E_ge_two,
    E_ge_one,
    n_LB_int,
    n_UB_int,
    bl_LB,
    bl_UB,
    ch_LB,
    ch_UB,
    ell_UB,
    I_ge_two,
    I_ge_one,
    I_equal_one,
    I_zero_one,
    bd2_LB,
    bd3_LB,
    val,
    Lambda_star
):
    t_C = len(V_C)
    t_C_tilde = sum(1 for i in V_C if bl_UB[i] >= 1)
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
        
    return t_C, t_C_tilde, t_T, t_F, k_C, k_C_tilde, c_F, m_C, \
        head_C, tail_C, tail_F, E_C_plus, E_C_minus, \
        I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, \
        I_zero_one_plus, I_zero_one_minus, I_equal_one_plus, I_equal_one_minus, \
        delta_i

def prepare_fringe_tree(
    set_F,
    F_Lambda,
    V_C,
    t_T,
    t_F,
    val,
    Lambda_int,
    Lambda_ex,
    Code_Lambda
    ):
    
    set_F_E = set_F
    set_F_v = {v : set_F for v in V_C}

    F_C = set_F_v
    F_T = {i : set_F_E for i in range(1, t_T + 1)}
    F_F = {i : set_F_E for i in range(1, t_F + 1)}
    n_C = max(len(psi.vertex) - 1 for v in V_C for psi in set_F_v[v])
    n_T = max(len(psi.vertex) - 1 for psi in set_F_E)
    n_F = max(len(psi.vertex) - 1 for psi in set_F_E)

    Code_F = {psi : i + 1 for i, psi in enumerate(set_F + F_Lambda)}

    # number of non-root vertices
    n_psi = {Code_F[psi] : len(psi.vertex) - 1 for psi in set_F}

    # degree of root
    deg_r = {Code_F[psi] : len(psi.adj[0]) for psi in set_F}

    # element of root
    atom_r = {Code_F[psi] : psi.root[0] for psi in (set_F + F_Lambda)} 
    alpha_r = {i : Code_Lambda[atom_r[i]] for i in atom_r}

    # sum of mul of edges incident to root
    beta_r = {Code_F[psi] : sum(psi.beta[0][v] for v in psi.adj[0]) for psi in set_F}

    ht = {Code_F[psi] : psi.height for psi in set_F}

    for psi in F_Lambda:
        n_psi[Code_F[psi]] = 0
        deg_r[Code_F[psi]] = 0
        beta_r[Code_F[psi]] = 0
        ht[Code_F[psi]] = 0

    n_H = dict()
    na_alpha_ex = {ele : {i + 1 : 0} for i in range(len(set_F)) for ele in Lambda_ex}
    for i, psi in enumerate(set_F):
        n_H_tmp = {d : 0 for d in range(MAX_VAL)}
        na_ex_tmp = {ele : 0 for ele in Lambda_ex}
        for u, (ele, dep) in enumerate(psi.vertex[1:]):
            beta_tmp = 0
            na_ex_tmp[ele] += 1
            for v in psi.adj[u + 1]:
                beta_tmp += psi.beta[u + 1][v]
            d_tmp = val[ele] - beta_tmp
            n_H_tmp[d_tmp] += 1

        for ele, d in na_alpha_ex.items():
            d[i + 1] = na_ex_tmp[ele]

        n_H[i + 1] = n_H_tmp

    F_Cp = dict()
    for i in F_C.keys():
        F_Cp_tmp = dict()
        for psi in F_C[i]:
            p = ht[Code_F[psi]]
            if p not in F_Cp_tmp.keys():
                F_Cp_tmp[p] = {psi}
            else:
                F_Cp_tmp[p].add(psi)
        F_Cp[i] = F_Cp_tmp

    F_Tp = dict()
    for i in F_T.keys():
        F_Tp_tmp = dict()
        for psi in F_T[i]:
            p = ht[Code_F[psi]]
            if p not in F_Tp_tmp.keys():
                F_Tp_tmp[p] = {psi}
            else:
                F_Tp_tmp[p].add(psi)
        F_Tp[i] = F_Tp_tmp

    F_Fp = dict()
    for i in F_F.keys():
        F_Fp_tmp = dict()
        for psi in F_F[i]:
            p = ht[Code_F[psi]]
            if p not in F_Fp_tmp.keys():
                F_Fp_tmp[p] = {psi}
            else:
                F_Fp_tmp[p].add(psi)
        F_Fp[i] = F_Fp_tmp

    return Code_F, n_psi, deg_r, beta_r, atom_r, ht, F_C, F_T, F_F, \
            n_C, n_T, n_F, F_Cp, F_Tp, F_Fp, set_F_E, set_F_v, n_H, \
            na_alpha_ex, alpha_r

def read_fringe_tree(fv_fringe_tree_filename, strF):

    fv_fringe_tree = dict()
    index_set = dict()
    with open(fv_fringe_tree_filename,'r') as f:
        lines = f.readlines()[1:]
        for line in lines:
            ind = int(line.split(",")[0])
            str1 = line.split(",")[1]
            str2 = line.split(",")[2].replace('\n', '')
            fv_fringe_tree[ind] = (str1, str2)
            for key, f in strF.items():
                if (str1, str2) == f:
                    index_set[ind] = key

    # print(fv_fringe_tree)
    # print(index_set)

    return fv_fringe_tree, index_set

def prepare_fv( 
    ann_training_data_filename,
    Lambda_int,
    Lambda_ex,
    Gamma_int,
    fv_fringe_tree,
    index_set,
    n_G,
    n_G_int,
    MASS,
    dg,
    dg_int,
    hydg,
    bd_int,
    na_int,
    na_ex,
    ec_int,
    fc
):
    ''' A function to prepare some variables about the descriptors 
        from the feature vector '''

    descriptors = dict()      # get variable from fv file
    # stringoutput = dict()     # string used for output
    fv = pd.read_csv(ann_training_data_filename, sep=",")
    num_fv = len(list(fv.columns))

    # prepare for normalization and standardization
    max_dcp = dict()
    min_dcp = dict()
    avg_dcp = dict()
    sd_dcp = dict()
    for i in range(1,num_fv):
        max_dcp[i] = fv.iloc[:, i].max()
        min_dcp[i] = fv.iloc[:, i].min()
        avg_dcp[i] = fv.iloc[:, i].mean()
        sd_dcp[i] = fv.iloc[:, i].std()

    mass_ind = -1

    forbidden_node = list()

    for i, fv_name in enumerate(list(fv.columns)):
        if fv_name == "CID":
            pass
        elif fv_name == "n":
            descriptors[i] = n_G
            # stringoutput[i] = "               n : "
        elif fv_name == "n_in":
            descriptors[i] = n_G_int
        elif fv_name == "ms":
            # descriptors[i] = MASS# / (10 * n_G)
            mass_ind = i
            # stringoutput[i] = "             M/n : "
        elif fv_name[0] == "d":
            if fv_name[3] == "i":                
                j = int(fv_name[6])
                descriptors[i] = dg_int[j]
            elif fv_name[3] == "h":
                j = int(fv_name[7])
                descriptors[i] = hydg[j]
            else:
                j = int(fv_name[3])
                descriptors[i] = dg[j]
        elif fv_name == "bd_in_2":
            descriptors[i] = bd_int[2]
            # stringoutput[i] = "    #double_bond_in : "
        elif fv_name == "bd_in_3":
            descriptors[i] = bd_int[3]
            # stringoutput[i] = "    #double_bond_in : "
        elif fv_name[0] == "n" and fv_name[3] == "i":
            ele = fv_name[6:] # may need change, Zhu, 0109
            if ele in Lambda_int:
                descriptors[i] = na_int[ele]
            else:
                descriptors[i] = 0
            # stringoutput[i] = "        #degree{}_ex : ".format(j)
        elif fv_name[0] == "n" and fv_name[3] == "e":
            ele = fv_name[6:]
            if ele in Lambda_ex:
                descriptors[i] = na_ex[ele]
            else:
                descriptors[i] = 0
            # stringoutput[i] = "        #degree{}_ex : ".format(j)
        elif fv_name[0] == "F" and fv_name[1] == "C":
            index = int(fv_name.split("_")[1])
            if index in index_set.keys():
                descriptors[i] = fc[index_set[index]]
            else:
                descriptors[i] = 0
                forbidden_node.append((1, i))
        else:
            str_tmp = fv_name
            pos = list()
            posnum = -1
            for k in range(len(str_tmp)):
                if str_tmp[k] == '_':
                    pos.append(k)
            ele1 = str_tmp[: pos[0] - 1]
            d1 = int(str_tmp[pos[0] - 1])
            ele2 = str_tmp[pos[0] + 1: pos[1] - 1]
            d2 = int(str_tmp[pos[1] - 1])
            mul = int(str_tmp[-1])
            ec_tmp = ((ele1, d1), (ele2, d2), mul)
            if ec_tmp in Gamma_int:
                descriptors[i] = ec_int[ec_tmp]
            else:
                descriptors[i] = 0

    return descriptors, num_fv, mass_ind, max_dcp, min_dcp, avg_dcp, sd_dcp, forbidden_node

def add_constraints_mass_n(MILP,
    n_LB, n_star, n_G, MASS):
    
    delta_n = {n: pulp.LpVariable(f"delta_n({n})".format(n), 0, 1, cat=pulp.LpBinary) for n in range(n_LB, n_star + 1)}

    MILP += pulp.lpSum(delta_n[n] for n in range(n_LB, n_star + 1)) == 1, f"milp-ann-n-1"
    MILP += pulp.lpSum(n * delta_n[n] for n in range(n_LB, n_star + 1)) == n_G, f"milp-ann-n-2"

    mass_n = pulp.LpVariable("mass_n", 0, 100000)
    for n in range(n_LB, n_star + 1):
        MILP += mass_n <= MASS * (1 / n) + 100000 * (1 - delta_n[n]), f"milp-mass_n-{n}-1"
        MILP += mass_n >= MASS * (1 / n) - 100000 * (1 - delta_n[n]), f"milp-mass_n-{n}-2"

    return MILP, mass_n

def add_constraints_ANN(MILP,
                        x_hat,
                        num_fv,
                        y, mass_ind,
                        mass_n, 
                        forbidden_node,
                        prop = "def"):
    ''' A function to add constraints used in ANN '''

    for i in range(1, num_fv):
        # if i == mass_ind:  # Mass
        #     MILP += y[(1, i)] == x_hat[i],\
        #     "ann_input_{}_{}".format(i, prop)
        # elif (1, i) not in forbidden_node:
        #     MILP += y[(1, i)] == x_hat[i], \
        #     "ann_input_{}_{}".format(i, prop)
        if (1, i) not in forbidden_node:
            MILP += y[(1, i)] == x_hat[i], \
            "ann_input_{}_{}".format(i, prop)

    return MILP

# -------- A.12 Constraints for Normalization or Standardization of Feature Vectors --------
def prepare_variables_nor_std_fv(num_fv):

    x_hat = {i: pulp.LpVariable(f"x_hat({i})") for i in range(1, num_fv)}
    x_tilde = {i: pulp.LpVariable(f"x_tilde({i})") for i in range(1, num_fv)}

    return x_hat, x_tilde

def add_constraints_nor_std_fv(
    MILP,
    num_fv,
    mass_ind,
    descriptors,
    mass_n,
    max_dcp,
    min_dcp,
    avg_dcp,
    sd_dcp,
    x_hat,
    x_tilde,
    eps,
    forbidden_node
):
    for i in range(1, num_fv):
        if i == mass_ind:
            if max_dcp[i] == min_dcp[i]:
                MILP += x_hat[i] >= (1 - eps) * (mass_n - min_dcp[i]), f"milp-(88)-{i}-1-mass"
                MILP += x_hat[i] <= (1 + eps) * (mass_n - min_dcp[i]), f"milp-(88)-{i}-2-mass"
            else:
                MILP += x_hat[i] >= (1 - eps) * (mass_n - min_dcp[i]) / (max_dcp[i] - min_dcp[i]), f"milp-(88)-{i}-1-mass"
                MILP += x_hat[i] <= (1 + eps) * (mass_n - min_dcp[i]) / (max_dcp[i] - min_dcp[i]), f"milp-(88)-{i}-2-mass"
            # MILP += x_tilde[i] >= (mass_n - eps - avg_dcp[i]) / sd_dcp[i], f"milp-(89)-{i}-1-mass"
            # MILP += x_tilde[i] <= (mass_n + eps - avg_dcp[i]) / sd_dcp[i], f"milp-(89)-{i}-2-mass"
        elif (1, i) not in forbidden_node:
            if max_dcp[i] == min_dcp[i]:
                MILP += x_hat[i] >= (1 - eps) * (descriptors[i] - min_dcp[i]), f"milp-(88)-{i}-1(nor)"
                MILP += x_hat[i] <= (1 + eps) * (descriptors[i] - min_dcp[i]), f"milp-(88)-{i}-2(nor)"
            else:
                MILP += x_hat[i] >= (1 - eps) * (descriptors[i] - min_dcp[i]) / (max_dcp[i] - min_dcp[i]), f"milp-(88)-{i}-1(nor)"
                MILP += x_hat[i] <= (1 + eps) * (descriptors[i] - min_dcp[i]) / (max_dcp[i] - min_dcp[i]), f"milp-(88)-{i}-2(nor)"
            # MILP += x_tilde[i] >= (descriptors[i] - eps - avg_dcp[i]) / sd_dcp[i], f"milp-(89)-{i}-1(std)"
            # MILP += x_tilde[i] <= (descriptors[i] + eps - avg_dcp[i]) / sd_dcp[i], f"milp-(89)-{i}-2(std)"

    return MILP

# -------- A.1 Selecting Core-vertices and Core-edges   -------- 
def prepare_variables_selecting_core(
    # Constants
    t_C,
    k_C_tilde,
    k_C,
    t_T,
    m_C,
    ell_LB,
    ell_UB,
    bl_LB,
    bl_UB
):
    # Define e_C
    e_C = {i : pulp.LpVariable(f"e_C({i})", 0, 1, cat = pulp.LpBinary)
            for i in range(1, m_C + 1)}

    # Define v_T
    v_T = {(i, 0) : pulp.LpVariable(f"v_T({i},0)", 0, 1, cat = pulp.LpBinary)
            for i in range(1, t_T + 1)}

    # Define e_T    
    e_T = {i : pulp.LpVariable(f"e_T({i})", 0, 1, cat = pulp.LpBinary)
            for i in range(1, (t_T + 1) + 1)}

    # Define chi_T
    chi_T = {i : pulp.LpVariable(f"chi_T({i})", 0, k_C, cat = pulp.LpInteger)
            for i in range(1, t_T + 1)}

    # Define clr_T
    clr_T = {k : pulp.LpVariable(f"clr_T({k})", ell_LB[k] - 1, ell_UB[k] - 1, cat = pulp.LpInteger)
            for k in range(1, k_C + 1)}
    clr_T[0] = pulp.LpVariable("clr_T(0)", 0, t_T, cat = pulp.LpInteger)

    # Define delta_chi_T
    delta_chi_T = {k : pulp.LpVariable(f"delta_chi_T({k})", 0, 1, cat = pulp.LpBinary)
            for k in range(0, k_C + 1)} 
    chi_T_tmp = {(i, k) : pulp.LpVariable(f"chi_T({i},{k})", 0, 1, cat = pulp.LpBinary)
            for i in range(1, t_T + 1) for k in range(0, k_C + 1)}
    chi_T.update(chi_T_tmp)

    # Define deg_tilde_C_plus
    deg_tilde_C_plus = {i : pulp.LpVariable(f"deg_tilde_C_plus({i})", 0, MAX_VAL, cat = pulp.LpInteger)
            for i in range(1, t_C + 1)}

    # Define deg_tilde_C_minus
    deg_tilde_C_minus = {i : pulp.LpVariable(f"deg_tilde_C_minus({i})", 0, MAX_VAL, cat = pulp.LpInteger)
            for i in range(1, t_C + 1)}

    return e_C, v_T, e_T, chi_T, clr_T, delta_chi_T, deg_tilde_C_plus, deg_tilde_C_minus

def add_constraints_selecting_core(
    # Model
    MILP,
    # Constants
    t_C,
    k_C_tilde,
    k_C,
    t_T,
    m_C,
    n_LB_int,
    n_UB_int,
    ell_LB,
    ell_UB,
    bl_LB,
    bl_UB,
    I_equal_one,
    I_equal_one_minus,
    I_equal_one_plus,
    I_ge_two,
    I_ge_one,
    I_ge_one_minus,
    I_ge_one_plus,
    I_zero_one_minus,
    I_zero_one_plus,
    # Binary Variables
    e_C,
    v_T,
    e_T,
    delta_chi_T,
    # Integer Variables
    chi_T,
    clr_T,
    deg_tilde_C_plus,
    deg_tilde_C_minus
)->pulp.LpProblem:
    # -------- Constraint (3) --------
    # Constraint of e_C, for i in I_(=1)
    for i in I_equal_one:
        MILP += e_C[i] == 1, f"milp-(3)-{i}"

    # -------- Constraint (4) --------
    # Constraint of e_C, clr_T, for i in I_(\ge 2)
    for i in I_ge_two:
        MILP += e_C[i] == 0, f"milp-(4)-1-{i}"
        MILP += clr_T[i] >= 1, f"milp-(4)-2-{i}"

    # -------- Constraint (5) --------
    # Constraint of e_C, clr_T, for i in I_(\ge 1)
    for i in I_ge_one:
        MILP += e_C[i] + clr_T[i] >= 1, f"milp-(5)-1-{i}"
        MILP += clr_T[i] <= t_T * (1 - e_C[i]), f"milp-(5)-2-{i}"

    # -------- Constraint (6) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(e_C[c] for c in (I_ge_one_minus[i] + I_zero_one_minus[i] + I_equal_one_minus[i])) == \
                deg_tilde_C_minus[i], f"milp-(6)-1-{i}"
        MILP += pulp.lpSum(e_C[c] for c in (I_ge_one_plus[i] + I_zero_one_plus[i] + I_equal_one_plus[i])) == \
                deg_tilde_C_plus[i], f"milp-(6)-2-{i}"

    # -------- Constraint (7) --------
    # Constraint of delta_chi_T, for i in [1, t_T]
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(chi_T[(i, k)]
                for k in range(0, k_C + 1)) == 1, f"milp-(7)-2-{i}"        
        MILP += pulp.lpSum(k * chi_T[(i, k)]
                for k in range(0, k_C + 1)) == chi_T[i], f"milp-(7)-3-{i}"
        MILP += chi_T[(i, 0)] == 1 - v_T[(i, 0)], f"milp-(7)-1-{i}"

    # -------- Constraint (8) --------
    # Constraint of delta_chi_T, clr_T, for k in [0, k_C]
    for k in range(0, k_C + 1):
        MILP += pulp.lpSum(chi_T[(i, k)]
                for i in range(1, t_T + 1)) == clr_T[k], f"milp-(8)-1-{k}"
        MILP += pulp.lpSum(chi_T[(i, k)]
                for i in range(1, t_T + 1)) <= t_T * delta_chi_T[k], f"milp-(8)-2-{k}-1"
        MILP += pulp.lpSum(chi_T[(i, k)]
                for i in range(1, t_T + 1)) >= delta_chi_T[k], f"milp-(8)-2-{k}-2"

    # -------- Constraint (9) --------
    # Constraint of e_T, v_T, chi_T, for i in [2, t_T]
    for i in range(2, t_T + 1):
        MILP += v_T[(i - 1, 0)] >= v_T[(i, 0)], f"milp-(9)-1-{i}"
        MILP += chi_T[i - 1] - chi_T[i] <= k_C * (v_T[(i - 1, 0)] - e_T[i]), f"milp-(9)-2-{i}-1"
        MILP += chi_T[i - 1] - chi_T[i] >= v_T[(i - 1, 0)] - e_T[i], f"milp-(9)-2-{i}-2"


    return MILP

# -------- A.2 Constraints for Adding Internal Vertices and Edges --------
def prepare_variables_internal_vertices_and_edges(
    t_F,
    t_T,
    t_C,
    t_C_tilde,
    c_F,
    V_C,
    E_ge_one,
    E_ge_two,
    I_ge_one,
    I_ge_two,
    bl_LB,
    bl_UB,
    E_C,
    n_LB_int,
    n_UB_int
):
    # Define bl_G
    # Here max function should be sum, by Zhu.
    if len(E_ge_one + E_ge_two) == 0:
        bl_G_UB = max(bl_UB[v] for v in V_C)
    else:
        bl_G_UB = max(max(bl_UB[v] for v in V_C), max(bl_UB[e] for e in (E_ge_one + E_ge_two)))
    bl_G = pulp.LpVariable("bl_G", 0, bl_G_UB, cat = pulp.LpInteger)
    # Define v_F
    v_F = {(i, 0) : pulp.LpVariable(f"v_F({i},0)", 0, 1, cat = pulp.LpBinary)
                for i in range(1, t_F + 1)}
    # Define e_F
    e_F = {i : pulp.LpVariable(f"e_F({i})", 0, 1, cat = pulp.LpBinary)
                for i in range(1, (t_F + 1) + 1)}
    # Define chi_F
    chi_F = {i : pulp.LpVariable(f"chi_F({i})", 0, c_F, cat = pulp.LpInteger)
                for i in range(1, t_F + 1)}
    # Define clr_F
    clr_F = {c : pulp.LpVariable(f"clr_F({c})", 0, t_F, cat = pulp.LpInteger)
                for c in range(0, c_F + 1)}
    # Define delta_chi_F
    delta_chi_F = {c : pulp.LpVariable(f"delta_chi_F({c})", bl_LB[c], 1, cat = pulp.LpInteger)
                for c in range(1, t_C_tilde + 1)}
    delta_chi_F[0] = pulp.LpVariable(f"delta_chi_F(0)", 0, 1, cat=pulp.LpBinary)
    delta_chi_F_tmp = {c: pulp.LpVariable(f"delta_chi_F({c})", 0, 1, cat=pulp.LpBinary)
                for c in range(t_C_tilde + 1, c_F + 1)}
    delta_chi_F.update(delta_chi_F_tmp)
    chi_F_tmp = {(i, c) : pulp.LpVariable(f"chi_F({i},{c})", 0, 1, cat = pulp.LpBinary)
                for i in range(1, t_F + 1) for c in range(0, c_F + 1)}
    chi_F.update(chi_F_tmp)
    # Define bl
    bl = {(k, i) : pulp.LpVariable(f"bl({k},{i})", 0, 1, cat = pulp.LpBinary)
                for k in (I_ge_two + I_ge_one) for i in range(1, t_T + 1)}

    # Define n_G_int
    n_G_int = pulp.LpVariable("n_G_int", n_LB_int, n_UB_int, cat = pulp.LpInteger)

    return bl_G, v_F, e_F, chi_F, clr_F, delta_chi_F, bl, n_G_int

def add_constraints_internal_vertices_and_edges(
    # Model
    MILP,
    # Constants
    t_T,
    t_F,
    t_C,
    t_C_tilde,
    c_F,
    tail_F,
    I_ge_one,
    I_ge_two,
    bl_LB,
    bl_UB,
    E_C,
    # Binary Variables
    delta_chi_F,
    chi_T,
    e_F,
    v_F,
    v_T,
    # Integer Variables
    chi_F,
    clr_F,
    bl,
    bl_G,
    n_G_int
)->pulp.LpProblem:
    # -------- Constraint (10) --------
    # Constraint of delta_chi_F, chi_F, for i in [1, t_F]
    for i in range(1, t_F + 1):
        MILP += chi_F[(i, 0)] == 1 - v_F[(i, 0)], f"milp-(10)-1-{i}"        
        MILP += pulp.lpSum(chi_F[(i, c)]
                for c in range(0, c_F + 1)) == 1, f"milp-(10)-2-{i}"
        MILP += pulp.lpSum(c * chi_F[(i, c)]
                for c in range(0, c_F + 1)) == chi_F[i], f"milp-(10)-3-{i}"

    # -------- Constraint (11) --------
    # Constraint of delta_chi_F, clr_F, for c in [0, c_F]
    for c in range(0, c_F + 1):
        MILP += pulp.lpSum(chi_F[(i, c)]
                for i in range(1, t_F + 1)) == clr_F[c], f"milp-(11)-1-{c}"
        MILP += pulp.lpSum(chi_F[(i, c)]
                for i in range(1, t_F + 1)) <= t_F * delta_chi_F[c], f"milp-(11)-2-{c}-1"
        MILP += pulp.lpSum(chi_F[(i, c)]
                for i in range(1, t_F + 1)) >= delta_chi_F[c], f"milp-(11)-2-{c}-2"

    # -------- Constraint (12) --------
    # Constraint of e_F
    MILP += e_F[1] == 0, f"milp-(12)-1"
    MILP += e_F[t_F + 1] == 0, f"milp-(12)-2"

    # -------- Constraint (13) --------
    # Constraint of e_F, chi_F, v_F, for i in [2, t_F]
    for i in range(2, t_F + 1):
        MILP += v_F[(i - 1, 0)] >= v_F[(i, 0)], f"milp-(13)-1-{i}"
        MILP += chi_F[i - 1] - chi_F[i] <= c_F * (v_F[(i - 1, 0)] - e_F[i]), f"milp-(13)-2-{i}-1"
        MILP += chi_F[i - 1] - chi_F[i] >= v_F[(i - 1, 0)] - e_F[i], f"milp-(13)-2-{i}-2"

    # -------- Constraint (14) --------
    MILP += t_C + pulp.lpSum(v_T[(i, 0)] for i in range(1, t_T + 1)) + \
            pulp.lpSum(v_F[(i, 0)] for i in range(1, t_F + 1)) == n_G_int, f"milp-(14)"

    # -------- Constraint (15) --------
    MILP += pulp.lpSum(delta_chi_F[c] for c in range(1, c_F + 1)) == bl_G, f"milp-(15)"

    # -------- Constraint (16) --------
    for k in (I_ge_two + I_ge_one):
        for i in range(1, t_T + 1):
            MILP += bl[(k, i)] >= delta_chi_F[t_C_tilde + i] + chi_T[(i, k)] - 1, f"milp-(16)-{k}-{i}"

    # -------- Constraint (17) --------
    MILP += pulp.lpSum(bl[(k, i)] for k in (I_ge_one + I_ge_two) for i in range(1, t_T + 1)) <= \
            pulp.lpSum(delta_chi_F[t_C_tilde + i] for i in range(1, t_T + 1)), f"milp-(17)"

    # -------- Constraint (18) --------
    for k in (I_ge_two + I_ge_one):
        MILP += pulp.lpSum(bl[(k, i)] for i in range(1, t_T + 1)) >= bl_LB[E_C[k]], f"milp-(18)-{k}-1"
        MILP += pulp.lpSum(bl[(k, i)] for i in range(1, t_T + 1)) <= bl_UB[E_C[k]], f"milp-(18)-{k}-2"
    
    return MILP

# -------- A.3 Constraints for Including Fringe-trees --------
def prepare_variables_fringe_trees(
    n_LB, 
    n_star,
    rho,
    ch_LB,
    ch_UB,
    t_T,
    t_C, 
    t_F,
    n_T,
    n_C,
    n_F,
    delta_i,
    I_ge_two,
    I_ge_one,
    V_C,
    E_ge_one,
    E_ge_two,
    v_T,
    v_F,
    F_C,
    F_T,
    F_F,
    Code_F,
    ht,
    F_Lambda
):
    # Define n_G
    n_G = pulp.LpVariable(f"n_G", n_LB, n_star, cat=pulp.LpInteger)

    # Define h_T, h_C, h_F
    h_T = {i: pulp.LpVariable(f"h_T({i})", 0, rho, cat=pulp.LpInteger)
            for i in range(1, t_T + 1)}
    h_C = {i: pulp.LpVariable(f"h_C({i})", 0, rho, cat=pulp.LpInteger)
            for i in range(1, t_C + 1)}
    h_F = {i: pulp.LpVariable(f"h_F({i})", 0, rho, cat=pulp.LpInteger)
            for i in range(1, t_F + 1)}

    # Define delta_fr_C, delta_fr_T, delta_fr_F
    delta_fr_C = {(i, Code_F[psi]): pulp.LpVariable(f"delta_fr_C({i},{Code_F[psi]})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_C + 1) for psi in (F_C[i] + F_Lambda)}    
    delta_fr_T = {(i, Code_F[psi]): pulp.LpVariable(f"delta_fr_T({i},{Code_F[psi]})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_T + 1) for psi in (F_T[i] + F_Lambda)}
    delta_fr_F = {(i, Code_F[psi]): pulp.LpVariable(f"delta_fr_F({i},{Code_F[psi]})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_F + 1) for psi in (F_F[i] + F_Lambda)}  

    # Define deg_ex_C, deg_ex_T, deg_ex_F
    deg_ex_C = {i: pulp.LpVariable(f"deg_ex_C({i})", 0, MAX_VAL - 1, cat=pulp.LpInteger)
                for i in range(1, t_C + 1)}
    deg_ex_T = {i: pulp.LpVariable(f"deg_ex_T({i})", 0, MAX_VAL - 1, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    deg_ex_F = {i: pulp.LpVariable(f"deg_ex_F({i})", 0, MAX_VAL - 1, cat=pulp.LpInteger)
                for i in range(1, t_F + 1)}    

    # Define sigma
    sigma = {(k, i): pulp.LpVariable(f"sigma({k},{i})", 0, 1, cat=pulp.LpBinary)
                    for k in (I_ge_two + I_ge_one) for i in range(1, t_T + 1)} 

    return n_G, h_T, h_C, h_F, delta_fr_C, delta_fr_T, delta_fr_F, \
        deg_ex_C, deg_ex_T, deg_ex_F, sigma

def add_constraints_fringe_trees(
    # Model
    MILP,
    # Constants
    t_T,
    t_F,
    t_C,
    t_C_tilde,
    c_F,
    I_ge_one,
    I_ge_two,
    rho,
    n_star,
    n_T,
    n_C,
    n_F,
    ch_LB,
    ch_UB,
    E_C,
    F_C,
    F_T,
    F_F,
    F_Cp,
    F_Tp,
    F_Fp,
    F_Lambda,
    Code_F,
    val,
    n_psi,
    deg_r,
    ht,
    # Binary Variables
    delta_chi_F,
    delta_chi_T,
    e_F,
    v_F,
    v_T,
    h_T,
    h_C,
    h_F,
    sigma,
    delta_fr_C,
    delta_fr_T,
    delta_fr_F,
    # Integer Variables
    chi_T,
    chi_F,
    clr_F,
    n_G,
    deg_ex_C,
    deg_ex_T,
    deg_ex_F
)->pulp.LpProblem:
    # -------- Constraint (19) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_fr_C[(i, Code_F[psi])] 
                    for psi in (F_C[i] + F_Lambda)) == 1, f"milp-(19)-C-1-{i}"
        MILP += pulp.lpSum(deg_r[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
                    for psi in (F_C[i] + F_Lambda)) == deg_ex_C[i], f"milp-(19)-C-2-{i}"

    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(delta_fr_T[(i, Code_F[psi])] 
                    for psi in (F_T[i] + F_Lambda)) == v_T[(i, 0)], f"milp-(19)-T-1-{i}"
        MILP += pulp.lpSum(deg_r[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
                    for psi in (F_T[i] + F_Lambda)) == deg_ex_T[i], f"milp-(19)-T-2-{i}"

    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_fr_F[(i, Code_F[psi])] 
                    for psi in (F_F[i] + F_Lambda)) == v_F[(i, 0)], f"milp-(19)-F-1-{i}"
        MILP += pulp.lpSum(deg_r[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
                    for psi in (F_F[i] + F_Lambda)) == deg_ex_F[i], f"milp-(19)-F-2-{i}"

    # -------- Constraint (20) --------
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_fr_F[(i, Code_F[psi])]
                    for psi in F_Fp[i][rho]) >= v_F[(i, 0)] - e_F[i + 1], f"milp-(20)-{i}"

    # -------- Constraint (21) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(ht[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] 
                    for psi in F_C[i]) == h_C[i], f"milp-(21)-C-{i}"

    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(ht[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] 
                    for psi in F_T[i]) == h_T[i], f"milp-(21)-T-{i}"

    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(ht[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] 
                    for psi in F_F[i]) == h_F[i], f"milp-(21)-F-{i}"

    # -------- Constraint (22) --------
    MILP += pulp.lpSum(n_psi[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])]
                for i in range(1, t_C + 1) for psi in F_C[i]) + \
            pulp.lpSum(n_psi[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])]
                for i in range(1, t_T + 1) for psi in F_T[i]) + \
            pulp.lpSum(n_psi[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])]
                for i in range(1, t_F + 1) for psi in F_F[i]) + \
            pulp.lpSum(v_T[(i, 0)] for i in range(1, t_T + 1)) + \
            pulp.lpSum(v_F[(i, 0)] for i in range(1, t_F + 1)) + t_C == n_G, f"milp-(22)"

    # -------- Constraint (23) --------
    for i in range(1, t_C_tilde + 1):
        MILP += h_C[i] >= ch_LB[i] - n_star * delta_chi_F[i], f"milp-(23)-1-{i}"
        MILP += clr_F[i] + rho >= ch_LB[i], f"milp-(23)-2-{i}"
        MILP += h_C[i] <= ch_UB[i], f"milp-(23)-3-{i}"
        MILP += clr_F[i] + rho <= ch_UB[i] + n_star * (1 - delta_chi_F[i]), f"milp-(23)-4-{i}"

    # -------- Constraint (24) --------
    for i in range(t_C_tilde + 1, t_C + 1):
        MILP += h_C[i] >= ch_LB[i], f"milp-(24)-1-{i}"
        MILP += h_C[i] <= ch_UB[i], f"milp-(24)-2-{i}" 

    # -------- Constraint (25) --------
    for i in range(1, t_T + 1):
        for k in (I_ge_one + I_ge_two):
            MILP += h_T[i] <= ch_UB[E_C[k]] + \
                n_star * (delta_chi_F[t_C_tilde + i] + 1 - chi_T[(i, k)]), f"milp-(25)-1-{i}-{k}"
            MILP += clr_F[t_C_tilde + i] + rho <= ch_UB[E_C[k]] + \
                n_star * (2 - delta_chi_F[t_C_tilde + i] - chi_T[(i, k)]), f"milp-(25)-2-{i}-{k}"

    # -------- Constraint (26) --------
    for k in (I_ge_one + I_ge_two):
        MILP += pulp.lpSum(sigma[(k, i)] for i in range(1, t_T + 1)) == delta_chi_T[k], f"milp-(26)-{k}"

    # -------- Constraint (27) --------
    for i in range(1, t_T + 1):
        for k in (I_ge_one + I_ge_two):
            MILP += chi_T[(i, k)] >= sigma[(k, i)], f"milp-(27)-1-{i}-{k}"
            MILP += h_T[i] >= ch_LB[E_C[k]] - \
                n_star * (delta_chi_F[t_C_tilde + i] + 1 - sigma[(k, i)]), f"milp-(27)-2-{i}-{k}"
            MILP += clr_F[t_C_tilde + i] + rho >= ch_LB[E_C[k]] - \
                n_star * (2 - delta_chi_F[t_C_tilde + i] - sigma[(k, i)]), f"milp-(27)-3-{i}-{k}"

    return MILP

# -------- A.4 Descriptor for the Number of Specified Degree --------
def prepare_variables_degree(
    t_C,
    t_T,
    t_F,
    n_C,
    n_T,
    n_F,
    delta_i,
    n_star,
    dg_LB,
    dg_UB,
    rho
):
    # Define deg_C, deg_T, deg_F
    deg_C = {(i, j): pulp.LpVariable(f"deg_C({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_C + 1) for j in range(n_C + 1)}
    deg_T = {(i, j): pulp.LpVariable(f"deg_T({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_T + 1) for j in range(n_T + 1)}
    deg_F = {(i, j): pulp.LpVariable(f"deg_F({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_F + 1) for j in range(n_F + 1)}

    # Define deg_CT, deg_TC
    deg_CT = {i: pulp.LpVariable(f"deg_CT({i})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_C + 1)}
    deg_TC = {i: pulp.LpVariable(f"deg_TC({i})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_C + 1)}

    # Define delta_dg_C, delta_dg_T, delta_dg_F
    delta_dg_C = {(i, 0, d): pulp.LpVariable(f"delta_dg_C({i},0,{d})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_C + 1) for d in range(1, MAX_VAL + 1)}
    delta_dg_T = {(i, 0, d): pulp.LpVariable(f"delta_dg_T({i},0,{d})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_T + 1) for d in range(0, MAX_VAL + 1)}
    delta_dg_F = {(i, 0, d): pulp.LpVariable(f"delta_dg_F({i},0,{d})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_F + 1) for d in range(0, MAX_VAL + 1)}

    # Define dg
    dg = {d: pulp.LpVariable(f"dg({d})", dg_LB[d], dg_UB[d], cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)} 

    # Define deg_int_C, deg_int_T, deg_int_F
    deg_int_C = {i: pulp.LpVariable(f"deg_int_C({i})", 1, MAX_VAL, cat=pulp.LpInteger)
                for i in range(1, t_C + 1)}
    deg_int_T = {i: pulp.LpVariable(f"deg_int_T({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    deg_int_F = {i: pulp.LpVariable(f"deg_int_F({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                for i in range(1, t_F + 1)}

    #Define delta_int_dg_C, delta_int_dg_T, delta_int_dg_F
    delta_int_dg_C = {(i, d): pulp.LpVariable(f"delta_int_dg_C({i},{d})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_C + 1) for d in range(1, MAX_VAL + 1)}
    delta_int_dg_T = {(i, d): pulp.LpVariable(f"delta_int_dg_T({i},{d})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_T + 1) for d in range(0, MAX_VAL + 1)}
    delta_int_dg_F = {(i, d): pulp.LpVariable(f"delta_int_dg_F({i},{d})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_F + 1) for d in range(0, MAX_VAL + 1)}

    # Define dg_int
    dg_int = {d: pulp.LpVariable(f"dg_int({d})", dg_LB[d], dg_UB[d], cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}

    return deg_C, deg_T, deg_F, deg_CT, deg_TC, \
        delta_dg_C, delta_dg_T, delta_dg_F, dg, \
        deg_int_C, deg_int_T, deg_int_F, \
        delta_int_dg_C, delta_int_dg_T, delta_int_dg_F, \
        dg_int

def add_constraints_degree(
    # Model
    MILP,
    # Constants
    t_T,
    t_F,
    t_C,
    t_C_tilde,
    n_T,
    n_C,
    n_F,
    delta_i,
    I_ge_one_plus,
    I_ge_one_minus,
    I_ge_two_plus,
    I_ge_two_minus,
    rho,
    F_C,
    F_T,
    F_F,
    F_Cp,
    F_Tp,
    F_Fp,
    Code_F,
    deg_r,
    # Binary Variables
    delta_dg_C,
    delta_dg_T,
    delta_dg_F,
    delta_int_dg_C,
    delta_int_dg_T,
    delta_int_dg_F,
    e_T,
    e_F,
    v_F,
    v_T,
    delta_chi_T,
    delta_chi_F,
    delta_fr_C,
    # Integer Variables
    deg_C,
    deg_T,
    deg_F,
    deg_CT,
    deg_TC,
    dg,
    deg_int_C,
    deg_int_T,
    deg_int_F,
    dg_int,
    deg_tilde_C_minus,
    deg_tilde_C_plus,
    deg_ex_C,
    deg_ex_T,
    deg_ex_F
):
    # -------- Constraint (28) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_chi_T[k] 
                    for k in (I_ge_two_plus[i] + I_ge_one_plus[i])) == deg_CT[i], f"milp-(28)-1-{i}"
        MILP += pulp.lpSum(delta_chi_T[k] 
                    for k in (I_ge_two_minus[i] + I_ge_one_minus[i])) == deg_TC[i], f"milp-(28)-2-{i}"

    # -------- Constraint (29) --------
    for i in range(1, t_C_tilde + 1):
        MILP += deg_tilde_C_minus[i] + deg_tilde_C_plus[i] + \
                deg_CT[i] + deg_TC[i] + \
                delta_chi_F[i] == deg_int_C[i], f"milp-(29)-{i}"

    # -------- Constraint (30) --------
    for i in range(t_C_tilde + 1, t_C + 1):
        MILP += deg_tilde_C_minus[i] + deg_tilde_C_plus[i] + \
                deg_CT[i] + deg_TC[i] == deg_int_C[i], f"milp-(30)-{i}"

    # -------- Constraint (31) --------
    for i in range(1, t_C + 1):
        MILP += deg_int_C[i] + deg_ex_C[i] == deg_C[(i, 0)], f"milp-(31)-{i}"

    # -------- Constraint (32) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_fr_C[(i, Code_F[psi])]
                    for psi in F_Cp[i][rho]) >= 2 - deg_int_C[i], f"milp-(32)-{i}"

    # -------- Constraint (33) --------
    MILP += e_T[1] == 0, f"milp-(33)-T1"
    MILP += e_T[t_T + 1] == 0, f"milp-(33)-T2"
    for i in range(1, t_T + 1):
        MILP += 2 * v_T[(i, 0)] + delta_chi_F[t_C_tilde + i] == \
                deg_int_T[i], f"milp-(33)-1-{i}"
        MILP += deg_int_T[i] + deg_ex_T[i] == deg_T[(i, 0)], f"milp-(33)-2-{i}"

    # -------- Constraint (34) --------
    for i in range(1, t_F + 1):
        MILP += v_F[(i, 0)] + e_F[i + 1] == deg_int_F[i], f"milp-(34)-1-{i}"
        MILP += deg_int_F[i] + deg_ex_F[i] == deg_F[(i, 0)], f"milp-(34)-2-{i}"

    # -------- Constraint (35) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_dg_C[(i, 0, d)] 
                for d in range(1, MAX_VAL + 1)) == 1, f"milp-(35)-C-1-{i}"
        MILP += pulp.lpSum(d * delta_dg_C[(i, 0, d)] 
                for d in range(1, MAX_VAL + 1)) == deg_C[(i, 0)], f"milp-(35)-C-2-{i}"
        MILP += pulp.lpSum(delta_int_dg_C[(i, d)]
                for d in range(1, MAX_VAL + 1)) == 1, f"milp-(35)-C-3-{i}"
        MILP += pulp.lpSum(d * delta_int_dg_C[(i, d)] 
                for d in range(1, MAX_VAL + 1)) == deg_int_C[i], f"milp-(35)-C-4-{i}"

    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(delta_dg_T[(i, 0, d)] 
                for d in range(0, MAX_VAL + 1)) == 1, f"milp-(35)-T-1-{i}"
        MILP += pulp.lpSum(d * delta_dg_T[(i, 0, d)] 
                for d in range(1, MAX_VAL + 1)) == deg_T[(i, 0)], f"milp-(35)-T-2-{i}"
        MILP += pulp.lpSum(delta_int_dg_T[(i, d)]
                for d in range(0, MAX_VAL + 1)) == 1, f"milp-(35)-T-3-{i}"
        MILP += pulp.lpSum(d * delta_int_dg_T[(i, d)] 
                for d in range(1, MAX_VAL + 1)) == deg_int_T[i], f"milp-(35)-T-4-{i}"

    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_dg_F[(i, 0, d)] 
                for d in range(0, MAX_VAL + 1)) == 1, f"milp-(35)-F-1-{i}"
        MILP += pulp.lpSum(d * delta_dg_F[(i, 0, d)] 
                for d in range(1, MAX_VAL + 1)) == deg_F[(i, 0)], f"milp-(35)-F-2-{i}"
        MILP += pulp.lpSum(delta_int_dg_F[(i, d)]
                for d in range(0, MAX_VAL + 1)) == 1, f"milp-(35)-F-3-{i}"
        MILP += pulp.lpSum(d * delta_int_dg_F[(i, d)] 
                for d in range(1, MAX_VAL + 1)) == deg_int_F[i], f"milp-(35)-F-4-{i}"

    # -------- Constraint (36) --------
    for d in range(1, MAX_VAL + 1):
        MILP += pulp.lpSum(delta_dg_C[(i, 0, d)] for i in range(1, t_C + 1)) + \
                pulp.lpSum(delta_dg_T[(i, 0, d)] for i in range(1, t_T + 1)) + \
                pulp.lpSum(delta_dg_F[(i, 0, d)] for i in range(1, t_F + 1)) == dg[d], f"milp-(36)-1-{d}"

    for d in range(1, MAX_VAL + 1):
        MILP += pulp.lpSum(delta_int_dg_C[(i, d)] for i in range(1, t_C + 1)) + \
                pulp.lpSum(delta_int_dg_T[(i, d)] for i in range(1, t_T + 1)) + \
                pulp.lpSum(delta_int_dg_F[(i, d)] for i in range(1, t_F + 1)) == dg_int[d], f"milp-(36)-2-{d}"

    return MILP

# -------- A.5 Assigning Multiplicity --------
def prepare_variables_multiplicity(
    t_C,
    t_T,
    t_F,
    n_C,
    n_T,
    n_F,
    c_F,
    delta_i,
    I_equal_one,
    I_zero_one,
    I_ge_one,
    I_ge_two,
    E_C,
    bd2_LB,
    bd2_UB,
    bd3_LB,
    bd3_UB,
    n_UB_int, 
    n_star,
    rho
):
    # Define beta_T, beta_F
    beta_T = {i: pulp.LpVariable(f"beta_T({i})", 0, MAX_BOND, cat=pulp.LpInteger)
            for i in range(1, t_T + 2)} # use [1, t_T + 1] just for convenience
    beta_F = {i: pulp.LpVariable(f"beta_F({i})", 0, MAX_BOND, cat=pulp.LpInteger)
            for i in range(1, t_F + 2)} # use [1, t_F + 1] just for convenience

    # Define beta_C
    beta_C = {i: pulp.LpVariable(f"beta_C({i})", 0, MAX_BOND, cat=pulp.LpInteger)
            for i in (I_equal_one + I_zero_one + I_ge_one)}

    # Define beta_plus, beta_minus
    beta_plus = {k: pulp.LpVariable(f"beta_plus({k})", 0, MAX_BOND, cat=pulp.LpInteger)
                for k in (I_ge_two + I_ge_one)}
    beta_minus = {k: pulp.LpVariable(f"beta_minus({k})", 0, MAX_BOND, cat=pulp.LpInteger)
                for k in (I_ge_two + I_ge_one)}

    # Define beta_in
    beta_in = {c: pulp.LpVariable(f"beta_in({c})", 0, MAX_BOND, cat=pulp.LpInteger)
                for c in range(1, c_F + 1)}

    # Define beta_ex_C, beta_ex_T, beta_ex_F
    beta_ex_C = {i: pulp.LpVariable(f"beta_ex_C({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                for i in range(1, t_C + 1)}
    beta_ex_T = {i: pulp.LpVariable(f"beta_ex_T({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    beta_ex_F = {i: pulp.LpVariable(f"beta_ex_F({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                for i in range(1, t_F + 1)}    

    # Define delta_beta_T, delta_beta_F
    delta_beta_T = {(i, m): pulp.LpVariable(f"delta_beta_T({i},{m})", 0, 1, cat=pulp.LpBinary)
                    for i in range(2, t_T + 1) for m in range(MAX_BOND + 1)}
    delta_beta_F = {(i, m): pulp.LpVariable(f"delta_beta_F({i},{m})", 0, 1, cat=pulp.LpBinary)
                    for i in range(2, t_F + 1) for m in range(MAX_BOND + 1)}

    # Define delta_beta_C
    delta_beta_C = {(i, m): pulp.LpVariable(f"delta_beta_C({i},{m})", 0, 1, cat=pulp.LpBinary)
                    for i in (I_equal_one + I_zero_one + I_ge_one) for m in range(MAX_BOND + 1)}

    # # Define delta_beta_C, delta_beta_T, delta_beta_F
    # delta_beta_C_temp = {(i, j, m): pulp.LpVariable(f"delta_beta_C({i},{j},{m})", 0, 1, cat=pulp.LpBinary)
    #                     for i in range(1, t_C + 1) for j in range(n_C + 1) for m in range(MAX_BOND + 1)}
    # delta_beta_C.update(delta_beta_C_temp)
    # delta_beta_T_temp = {(i, j, m): pulp.LpVariable(f"delta_beta_T({i},{j},{m})", 0, 1, cat=pulp.LpBinary)
    #                     for i in range(1, t_T + 1) for j in range(n_T + 1) for m in range(MAX_BOND + 1)}
    # delta_beta_T.update(delta_beta_T_temp)
    # delta_beta_F_temp = {(i, j, m): pulp.LpVariable(f"delta_beta_F({i},{j},{m})", 0, 1, cat=pulp.LpBinary)
    #                     for i in range(1, t_F + 1) for j in range(n_F + 1) for m in range(MAX_BOND + 1)}
    # delta_beta_F.update(delta_beta_F_temp)

    # Define delta_beta_plus, delta_beta_in
    # The range should be [0, 1] instead of [0, 3], by Zhu
    delta_beta_plus = {(k, m): pulp.LpVariable(f"delta_beta_plus({k},{m})", 0, 1, cat=pulp.LpBinary)
                        for k in (I_ge_two + I_ge_one) for m in range(MAX_BOND + 1)}
    delta_beta_minus = {(k, m): pulp.LpVariable(f"delta_beta_minus({k},{m})", 0, 1, cat=pulp.LpBinary)
                        for k in (I_ge_two + I_ge_one) for m in range(MAX_BOND + 1)}
    
    # Define delta_beta_in
    # The range should be [0, 1] instead of [0, 3], by Zhu
    delta_beta_in = {(c, m): pulp.LpVariable(f"delta_beta_in({c},{m})", 0, 1, cat=pulp.LpBinary)
                        for c in range(1, c_F + 1) for m in range(MAX_BOND + 1)}
    
    # Define bd_int
    bd_int = {m: pulp.LpVariable(f"bd_int({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}

    # Define bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF
    bd_T = {m: pulp.LpVariable(f"bd_T({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_C = {m: pulp.LpVariable(f"bd_C({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_CT = {m: pulp.LpVariable(f"bd_CT({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_TC = {m: pulp.LpVariable(f"bd_TC({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_F = {m: pulp.LpVariable(f"bd_F({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_CF = {m: pulp.LpVariable(f"bd_CF({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_TF = {m: pulp.LpVariable(f"bd_TF({m})", 0, 2 * n_UB_int, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    
    return beta_T, beta_F, beta_C, beta_plus, beta_minus, beta_in, beta_ex_C, beta_ex_T, beta_ex_F, \
        delta_beta_T, delta_beta_F, delta_beta_C, \
        delta_beta_plus, delta_beta_minus, delta_beta_in, bd_int, \
        bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF

def add_constraints_multiplicity(
    # Model
    MILP,
    # Constants
    t_T,
    t_F,
    t_C,
    t_C_tilde,
    n_T,
    n_C,
    n_F,
    c_F,
    delta_i,
    head_C,
    tail_C,
    I_equal_one,
    I_zero_one,
    I_ge_one,
    I_ge_two,
    E_C,
    bd2_LB,
    bd2_UB,
    bd3_LB,
    bd3_UB,
    rho,
    F_C,
    F_T,
    F_F,
    Code_F,
    beta_r,
    # Binary Variables
    e_C,
    e_T,
    e_F,
    v_F,
    v_T,
    delta_chi_T,
    delta_chi_F,
    delta_beta_C,
    delta_beta_T,
    delta_beta_F,
    delta_beta_plus,
    delta_beta_minus,
    delta_beta_in,
    delta_fr_C,
    delta_fr_T,
    delta_fr_F,
    # Integer Variables
    beta_C,
    beta_T,
    beta_F,
    beta_plus, 
    beta_minus, 
    beta_in,
    bd_int,
    bd_C,
    bd_T,
    bd_CT,
    bd_TC,
    bd_F,
    bd_CF,
    bd_TF,
    beta_ex_C, 
    beta_ex_T, 
    beta_ex_F
):
    # -------- Constraint (37) --------
    for i in (I_equal_one + I_zero_one + I_ge_one):
        MILP += beta_C[i] >= e_C[i], f"milp-(37)-{i}-1"
        MILP += beta_C[i] <= MAX_BOND * e_C[i], f"milp-(37)-{i}-2"

    # -------- Constraint (38) --------
    for i in range(2, t_T + 1):
        MILP += beta_T[i] >= e_T[i], f"milp-(38)-T-{i}-1"
        MILP += beta_T[i] <= MAX_BOND * e_T[i], f"milp-(38)-T-{i}-2"
    for i in range(2, t_F + 1):
        MILP += beta_F[i] >= e_F[i], f"milp-(38)-F-{i}-1"
        MILP += beta_F[i] <= MAX_BOND * e_F[i], f"milp-(38)-F-{i}-2"

    # -------- Constraint (39) --------
    for k in (I_ge_two + I_ge_one):
        MILP += beta_plus[k] >= delta_chi_T[k], f"milp-(39)-1-{k}-1"
        MILP += beta_plus[k] <= MAX_BOND * delta_chi_T[k], f"milp-(39)-1-{k}-2"
        MILP += beta_minus[k] >= delta_chi_T[k], f"milp-(39)-2-{k}-1"
        MILP += beta_minus[k] <= MAX_BOND * delta_chi_T[k], f"milp-(39)-2-{k}-2"

    # -------- Constraint (40) --------
    for c in range(1, c_F + 1):
        MILP += beta_in[c] >= delta_chi_F[c], f"milp-(40)-{c}-1"
        MILP += beta_in[c] <= MAX_BOND * delta_chi_F[c], f"milp-(40)-{c}-2"

    # -------- Constraint (41) --------
    for i in range(2, t_T + 1):
        MILP += pulp.lpSum(delta_beta_T[(i, m)] 
            for m in range(MAX_BOND + 1)) == 1, f"milp-(41)-T-1-{i}"
        MILP += pulp.lpSum(m * delta_beta_T[(i, m)] 
            for m in range(MAX_BOND + 1)) == beta_T[i], f"milp-(41)-T-2-{i}"
    for i in range(2, t_F + 1):
        MILP += pulp.lpSum(delta_beta_F[(i, m)] 
            for m in range(MAX_BOND + 1)) == 1, f"milp-(41)-F-1-{i}"
        MILP += pulp.lpSum(m * delta_beta_F[(i, m)] 
            for m in range(MAX_BOND + 1)) == beta_F[i], f"milp-(41)-F-2-{i}"

    # -------- Constraint (42) --------
    for i in (I_equal_one + I_zero_one + I_ge_one):
        MILP += pulp.lpSum(delta_beta_C[(i, m)] 
            for m in range(MAX_BOND + 1)) == 1, f"milp-(42)-1-{i}"
        MILP += pulp.lpSum(m * delta_beta_C[(i, m)] 
            for m in range(MAX_BOND + 1)) == beta_C[i], f"milp-(42)-2-{i}"

    # -------- Constraint (43) --------
    for k in (I_ge_two + I_ge_one):
        MILP += pulp.lpSum(delta_beta_plus[(k, m)]
            for m in range(MAX_BOND + 1)) == 1, f"milp-(43)-1-{k}"
        MILP += pulp.lpSum(m * delta_beta_plus[(k, m)]
            for m in range(MAX_BOND + 1)) == beta_plus[k], f"milp-(43)-2-{k}"
    for k in (I_ge_two + I_ge_one):
        MILP += pulp.lpSum(delta_beta_minus[(k, m)]
            for m in range(MAX_BOND + 1)) == 1, f"milp-(43)-3-{k}"
        MILP += pulp.lpSum(m * delta_beta_minus[(k, m)]
            for m in range(MAX_BOND + 1)) == beta_minus[k], f"milp-(43)-4-{k}"
    for c in range(1, c_F + 1):
        MILP += pulp.lpSum(delta_beta_in[(c, m)]
            for m in range(MAX_BOND + 1)) == 1, f"milp-(43)-5-{c}"
        MILP += pulp.lpSum(m * delta_beta_in[(c, m)]
            for m in range(MAX_BOND + 1)) == beta_in[c], f"milp-(43)-6-{c}"

    # -------- Constraint (44) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(beta_r[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] 
                    for psi in F_C[i]) == beta_ex_C[i], f"milp-(44)-C-{i}"
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(beta_r[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] 
                    for psi in F_T[i]) == beta_ex_T[i], f"milp-(44)-T-{i}"
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(beta_r[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] 
                    for psi in F_F[i]) == beta_ex_F[i], f"milp-(44)-F-{i}"

    # -------- Constraint (45) --------
    for m in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum(delta_beta_C[(i, m)] for i in (I_equal_one + I_zero_one + I_ge_one)) == bd_C[m], f"milp-(45)-1-{m}"
        MILP += pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)) == bd_T[m], f"milp-(45)-2-{m}"
        MILP += pulp.lpSum(delta_beta_plus[(k, m)] for k in (I_ge_two + I_ge_one)) == bd_CT[m], f"milp-(45)-3-{m}"
        MILP += pulp.lpSum(delta_beta_minus[(k, m)] for k in (I_ge_two + I_ge_one)) == bd_TC[m], f"milp-(45)-4-{m}"
        MILP += pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)) == bd_F[m], f"milp-(45)-5-{m}"
        MILP += pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)) == bd_CF[m], f"milp-(45)-6-{m}"
        MILP += pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)) == bd_TF[m], f"milp-(45)-7-{m}"
        MILP += bd_C[m] + bd_T[m] + bd_F[m] + \
                bd_CT[m] + bd_TC[m] + bd_TF[m] + bd_CF[m] == bd_int[m], f"milp-(45)-8-{m}"

    return MILP

# -------- A.6 Assigning Chemical Elements and Valence Condition --------
def prepare_variables_chemical_elements(
    t_C,
    t_T,
    t_F,
    n_C,
    n_T,
    n_F,
    delta_i,
    n_star,
    Lambda,
    epsilon,
    na_LB,
    na_UB,
    rho,
    na_LB_int,
    na_UB_int,
    Lambda_int,
    Lambda_ex,
    MAX_CODE,
    MAX_CODE_int,
    MAX_CODE_ex,
    Code_Lambda_int,
    Code_Lambda_ex,
    F_C,
    F_T,
    F_F,
    Code_F,
    alpha_r,
    na_alpha_ex,
    n_H
):
    # Define beta_CT, beta_TC
    beta_CT = {i: pulp.LpVariable(f"beta_CT({i})", 0, MAX_BOND, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    beta_TC = {i: pulp.LpVariable(f"beta_TC({i})", 0, MAX_BOND, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}

    # Define beta_CF, beta_TF
    beta_CF = {i: pulp.LpVariable(f"beta_CF({i})", 0, MAX_BOND, cat=pulp.LpInteger)
                for i in range(1, t_F + 1)}
    beta_TF = {i: pulp.LpVariable(f"beta_TF({i})", 0, MAX_BOND, cat=pulp.LpInteger)
                for i in range(1, t_F + 1)}

    # Define alpha_C, alpha_T, alpha_F
    alpha_C = {(i, 0): pulp.LpVariable(f"alpha_C({i},{0})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                for i in range(1, t_C + 1)}
    alpha_T = {(i, 0): pulp.LpVariable(f"alpha_T({i},{0})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    alpha_F = {(i, 0): pulp.LpVariable(f"alpha_F({i},{0})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                for i in range(1, t_F + 1)}

    # # Define delta_alpha_C, delta_alpha_T, delta_alpha_F
    delta_alpha_C = {(i, 0, mu): pulp.LpVariable(f"delta_alpha_C({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_C + 1) for mu in (Lambda_int + [epsilon])}
    delta_alpha_T = {(i, 0, mu): pulp.LpVariable(f"delta_alpha_T({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_T + 1) for mu in (Lambda_int + [epsilon])}
    delta_alpha_F = {(i, 0, mu): pulp.LpVariable(f"delta_alpha_F({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_F + 1) for mu in (Lambda_int + [epsilon])}

    # Define MASS
    MASS = pulp.LpVariable(f"Mass", cat=pulp.LpInteger)

    # Define na
    na = {atom: pulp.LpVariable(f"na({atom})", na_LB[atom], na_UB[atom], cat=pulp.LpInteger)
            for atom in Lambda}

    # Define na_int, na_C, na_T, na_F
    na_int = {atom: pulp.LpVariable(f"na_int({atom})", na_LB_int[atom], na_UB_int[atom], cat=pulp.LpInteger) 
            for atom in Lambda}
    na_C = {atom: pulp.LpVariable(f"na_C({atom})", 0, na_UB_int[atom], cat=pulp.LpInteger) 
            for atom in Lambda}
    na_T = {atom: pulp.LpVariable(f"na_T({atom})", 0, na_UB_int[atom], cat=pulp.LpInteger)
            for atom in Lambda}
    na_F = {atom: pulp.LpVariable(f"na_F({atom})", 0, na_UB_int[atom], cat=pulp.LpInteger)
            for atom in Lambda}

    # Define na_ex_C, na_ex_T, na_ex_F
    na_ex_C = {atom: pulp.LpVariable(f"na_ex_C({atom})", 0, na_UB[atom], cat=pulp.LpInteger)
            for atom in Lambda_ex}
    na_ex_T = {atom: pulp.LpVariable(f"na_ex_T({atom})", 0, na_UB[atom], cat=pulp.LpInteger)
            for atom in Lambda_ex}
    na_ex_F = {atom: pulp.LpVariable(f"na_ex_F({atom})", 0, na_UB[atom], cat=pulp.LpInteger)
            for atom in Lambda_ex}

    na_ex = {atom: pulp.LpVariable(f"na_ex({atom})", 0, na_UB[atom], cat=pulp.LpInteger)
            for atom in Lambda_ex}

    # Define delta_hyd_C, T, F
    delta_hyd_C = {(i, d): pulp.LpVariable(f"delta_hyd_C({i},{d})", 0, 1, cat=pulp.LpBinary)
                for i in range(1, t_C + 1) for d in range(MAX_VAL)}
    delta_hyd_T = {(i, d): pulp.LpVariable(f"delta_hyd_T({i},{d})", 0, 1, cat=pulp.LpBinary)
                for i in range(1, t_T + 1) for d in range(MAX_VAL)}
    delta_hyd_F = {(i, d): pulp.LpVariable(f"delta_hyd_F({i},{d})", 0, 1, cat=pulp.LpBinary)
                for i in range(1, t_F + 1) for d in range(MAX_VAL)}

    # Define hydg
    hydg = {d: pulp.LpVariable(f"hydg({d})", 0, n_star, cat=pulp.LpInteger) for d in range(MAX_VAL)}

    return beta_CT, beta_TC, beta_CF, beta_TF, alpha_C, alpha_T, alpha_F, MASS, \
        delta_alpha_C, delta_alpha_T, delta_alpha_F, na, na_int, na_C, na_T, na_F, \
        na_ex_C, na_ex_T, na_ex_F, delta_hyd_C, delta_hyd_T, delta_hyd_F, hydg, \
        na_ex

def add_constraints_chemical_elements(
    # Model
    MILP,
    # Constants
    t_T,
    t_F,
    t_C,
    t_C_tilde,
    n_T,
    n_C,
    n_F,
    c_F,
    delta_i,
    E_C,
    I_equal_one,
    I_zero_one,
    I_ge_one,
    I_ge_one_plus,
    I_ge_one_minus,
    I_ge_two,
    I_ge_two_plus,
    I_ge_two_minus,
    val,
    mass,
    Lambda,
    Code_Lambda,
    Lambda_star,
    epsilon,
    rho,
    na_LB_int,
    na_UB_int,
    Lambda_int,
    Lambda_ex,
    MAX_CODE,
    MAX_CODE_int,
    MAX_CODE_ex,
    Code_Lambda_int,
    Code_Lambda_ex,
    F_C,
    F_T,
    F_F,
    F_Lambda,
    Code_F,
    alpha_r,
    na_alpha_ex,
    n_H,
    # Binary Variables
    v_T,
    v_F,
    e_T,
    e_F,
    delta_chi_T,
    delta_chi_F,
    chi_T,
    chi_F,
    delta_alpha_C,
    delta_alpha_T,
    delta_alpha_F,
    delta_fr_C,
    delta_fr_T,
    delta_fr_F,
    delta_hyd_C,
    delta_hyd_T,
    delta_hyd_F,
    # Integer Variables
    deg_C,
    deg_T,
    deg_F,
    beta_C,
    beta_T,
    beta_F,
    beta_CT,
    beta_TC,
    beta_CF,
    beta_TF,
    beta_plus,
    beta_minus,
    beta_in,
    alpha_C,
    alpha_T,
    alpha_F,
    MASS,
    na,
    na_int,
    na_C,
    na_T,
    na_F,
    na_ex_C,
    na_ex_T,
    na_ex_F,
    beta_ex_C,
    beta_ex_T,
    beta_ex_F,
    na_ex,
    bd_int,
    hydg
):
    # -------- Constraint (46) --------
    for k in (I_ge_one + I_ge_two):
        for i in range(1, t_T + 1):
            MILP += beta_CT[i] >= beta_plus[k] - MAX_BOND * (e_T[i] - chi_T[(i, k)] + 1), f"milp-(46)-1-{k}-{i}-1"
            MILP += beta_CT[i] <= beta_plus[k] + MAX_BOND * (e_T[i] - chi_T[(i, k)] + 1), f"milp-(46)-1-{k}-{i}-2"
        for i in range(1, t_T + 1):
            MILP += beta_TC[i] >= beta_minus[k] - MAX_BOND * (e_T[i + 1] - chi_T[(i, k)] + 1), f"milp-(46)-2-{k}-{i}-1"
            MILP += beta_TC[i] <= beta_minus[k] + MAX_BOND * (e_T[i + 1] - chi_T[(i, k)] + 1), f"milp-(46)-2-{k}-{i}-2"

    # -------- Constraint (47) --------
    for c in range(1, t_C_tilde + 1):
        for i in range(1, t_F + 1):
            MILP += beta_CF[i] >= beta_in[c] - MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-(47)-1-{c}-{i}-1"
            MILP += beta_CF[i] <= beta_in[c] + MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-(47)-1-{c}-{i}-2"
    for c in range(t_C_tilde + 1, c_F + 1):
        for i in range(1, t_F + 1):
            MILP += beta_TF[i] >= beta_in[c] - MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-(47)-2-{c}-{i}-1"
            MILP += beta_TF[i] <= beta_in[c] + MAX_BOND * (e_F[i] - chi_F[(i, c)] + 1), f"milp-(47)-2-{c}-{i}-2"

    # -------- Constraint (48) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_alpha_C[(i, 0, atom)] for atom in Lambda_int) == 1, \
                f"milp-(48)-C-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_int[atom] * delta_alpha_C[i, 0, atom]
                    for atom in Lambda_int) == alpha_C[(i, 0)], f"milp-(48)-C-2-{i}"

    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(delta_alpha_T[(i, 0, atom)] for atom in Lambda_int) == v_T[(i, 0)], \
                f"milp-(48)-T-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_int[atom] * delta_alpha_T[i, 0, atom]
                    for atom in Lambda_int) == alpha_T[(i, 0)], f"milp-(48)-T-2-{i}"

    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_alpha_F[(i, 0, atom)] for atom in Lambda_int) == v_F[(i, 0)], \
                f"milp-(48)-F-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_int[atom] * delta_alpha_F[i, 0, atom]
                    for atom in Lambda_int) == alpha_F[(i, 0)], f"milp-(48)-F-2-{i}"

    # -------- Constraint (49) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(alpha_r[Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] 
                    for psi in (F_C[i] + F_Lambda)) == alpha_C[(i, 0)], f"milp-(49)-C-{i}"
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(alpha_r[Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] 
                    for psi in (F_T[i] + F_Lambda)) == alpha_T[(i, 0)], f"milp-(49)-T-{i}"
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(alpha_r[Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] 
                    for psi in (F_F[i] + F_Lambda)) == alpha_F[(i, 0)], f"milp-(49)-F-{i}"

    # -------- Constraint (50) --------
    for i in range(1, t_C_tilde + 1):
        MILP += pulp.lpSum(beta_C[j] for j in range(1, len(E_C) + 1) 
                    if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
                pulp.lpSum(beta_plus[k] for k in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
                pulp.lpSum(beta_minus[k] for k in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
                beta_in[i] + beta_ex_C[i] + \
                pulp.lpSum(d * delta_hyd_C[(i, d)] for d in range(MAX_VAL)) == \
                pulp.lpSum(val[atom] * delta_alpha_C[(i, 0, atom)] for atom in Lambda_int), f"milp-(50)-{i}"

    # -------- Constraint (51) --------
    for i in range(t_C_tilde + 1, t_C + 1):
        MILP += pulp.lpSum(beta_C[j] for j in range(1, len(E_C) + 1) 
                    if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
                pulp.lpSum(beta_plus[c] for c in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
                pulp.lpSum(beta_minus[c] for c in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
                beta_ex_C[i] + pulp.lpSum(d * delta_hyd_C[(i, d)] for d in range(MAX_VAL)) == \
                pulp.lpSum(val[atom] * delta_alpha_C[(i, 0, atom)] for atom in Lambda_int), f"milp-(51)-{i}" 

    # -------- Constraint (52) --------
    MILP += beta_T[1] == 0, f"milp-(52)-T1"
    MILP += beta_T[t_T + 1] == 0, f"milp-(52)-T2"
    for i in range(1, t_T + 1):
        MILP += beta_T[i] + beta_T[i + 1] + beta_ex_T[i] + \
                beta_CT[i] + beta_TC[i] + beta_in[t_C_tilde + i] + \
                pulp.lpSum(d * delta_hyd_T[(i, d)] for d in range(MAX_VAL)) == \
                pulp.lpSum(val[atom] * delta_alpha_T[(i, 0, atom)] for atom in Lambda_int), f"milp-vector-(52)-{i}"

    # -------- Constraint (53) --------
    MILP += beta_F[1] == 0, f"milp-vector-(53)-F1"
    MILP += beta_F[t_F + 1] == 0, f"milp-vector-(53)-F2"
    for i in range(1, t_F + 1):
        MILP += beta_F[i] + beta_F[i + 1] + beta_ex_F[i] + \
                beta_CF[i] + beta_TF[i] + \
                pulp.lpSum(d * delta_hyd_F[(i, d)] for d in range(MAX_VAL)) == \
                pulp.lpSum(val[atom] * delta_alpha_F[(i, 0, atom)] for atom in Lambda_int), f"milp-vector-(53)-{i}"

    # -------- Constraint (54) --------
    for atom in Lambda_int:
        MILP += pulp.lpSum(delta_alpha_C[(i, 0, atom)] for i in range(1, t_C + 1)) == na_C[atom], f"milp-(54)-C-{atom}"
        MILP += pulp.lpSum(delta_alpha_T[(i, 0, atom)] for i in range(1, t_T + 1)) == na_T[atom], f"milp-(54)-T-{atom}"
        MILP += pulp.lpSum(delta_alpha_F[(i, 0, atom)] for i in range(1, t_F + 1)) == na_F[atom], f"milp-(54)-F-{atom}"

    # -------- Constraint (55) --------
    for atom in Lambda_ex: 
        MILP += pulp.lpSum(na_alpha_ex[atom][Code_F[psi]] * delta_fr_C[(i, Code_F[psi])] 
                    for i in range(1, t_C + 1) for psi in F_C[i]) == na_ex_C[atom], f"milp-(55)-C-{atom}"
        MILP += pulp.lpSum(na_alpha_ex[atom][Code_F[psi]] * delta_fr_T[(i, Code_F[psi])] 
                    for i in range(1, t_T + 1) for psi in F_T[i]) == na_ex_T[atom], f"milp-(55)-T-{atom}"
        MILP += pulp.lpSum(na_alpha_ex[atom][Code_F[psi]] * delta_fr_F[(i, Code_F[psi])] 
                    for i in range(1, t_F + 1) for psi in F_F[i]) == na_ex_F[atom], f"milp-(55)-F-{atom}"

    # -------- Constraint (56) --------
    for atom in Lambda_int:
        MILP += na_C[atom] + na_T[atom] + na_F[atom] == na_int[atom], f"milp-(56)-1-{atom}"
    for atom in Lambda_ex:
        MILP += na_ex_C[atom] + na_ex_T[atom] + na_ex_F[atom] == na_ex[atom], f"milp-(56)-2-{atom}"
    for atom in Lambda:
        if atom in Lambda_int and atom in Lambda_ex:
            MILP += na_int[atom] + na_ex[atom] == na[atom], f"milp-(56)-3-{atom}"
        elif atom in Lambda_int:
            MILP += na_int[atom] == na[atom], f"milp-(56)-4-{atom}"
        elif atom in Lambda_ex:
            MILP += na_ex[atom] == na[atom], f"milp-(56)-5-{atom}"

    # -------- Constraint (57) --------
    MILP += pulp.lpSum(mass[atom] * na[atom] for atom in Lambda) == MASS, f"milp-(57)"

    # -------- Constraint (58) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_hyd_C[(i, d)] for d in range(MAX_VAL)) == 1, f"milp-(58)-C-{i}"
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(delta_hyd_T[(i, d)] for d in range(MAX_VAL)) == v_T[(i, 0)], f"milp-(58)-T-{i}"
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_hyd_F[(i, d)] for d in range(MAX_VAL)) == v_F[(i, 0)], f"milp-(58)-F-{i}"

    # -------- Constraint (59) --------
    for d in range(MAX_VAL):
        MILP += pulp.lpSum(delta_hyd_C[(i, d)] for i in range(1, t_C + 1)) + \
                pulp.lpSum(delta_hyd_T[(i, d)] for i in range(1, t_T + 1)) + \
                pulp.lpSum(delta_hyd_F[(i, d)] for i in range(1, t_F + 1)) + \
                pulp.lpSum(n_H[Code_F[psi]][d] * delta_fr_C[(i, Code_F[psi])]
                    for i in range(1, t_C + 1) for psi in F_C[i]) + \
                pulp.lpSum(n_H[Code_F[psi]][d] * delta_fr_T[(i, Code_F[psi])]
                    for i in range(1, t_T + 1) for psi in F_T[i]) + \
                pulp.lpSum(n_H[Code_F[psi]][d] * delta_fr_F[(i, Code_F[psi])]
                    for i in range(1, t_F + 1) for psi in F_F[i]) == hydg[d], f"milp-(59)-{d}"

    # -------- Constraint (60) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_alpha_C[(i, 0, atom)] for atom in Lambda_star[i]) == 1, f"milp-(60)-{i}"

    return MILP

# -------- A.7 Constraints for Bounds on the Number of Bonds --------
def prepare_variables_number_of_bounds(
    k_C,
    t_T,
    bd_T
): 
    # Define bd_T
    bd_T_temp = {(k, i, m): pulp.LpVariable(f"bd_T({k},{i},{m})", 0, 1, cat=pulp.LpBinary)
                for k in range(1, k_C + 1) for i in range(2, t_T + 1) for m in {2, 3}}
    bd_T.update(bd_T_temp)

    return bd_T
    
def add_constraints_number_of_bounds(
    # Model
    MILP,
    # Constants
    t_T,
    k_C,
    E_C,
    I_equal_one,
    I_zero_one,
    bd2_LB,
    bd2_UB,
    bd3_LB,
    bd3_UB,
    # Binary Variables
    chi_T,
    delta_beta_C,
    delta_beta_T,
    delta_beta_plus,
    delta_beta_minus,
    bd_T
    # Integer Variables
):
    # -------- Constraint (61) --------
    for i in (I_equal_one + I_zero_one):
        MILP += delta_beta_C[(i, 2)] >= bd2_LB[E_C[i]], f"milp-(61)-{i}-2-1"
        MILP += delta_beta_C[(i, 2)] <= bd2_UB[E_C[i]], f"milp-(61)-{i}-2-2"
        MILP += delta_beta_C[(i, 3)] >= bd3_LB[E_C[i]], f"milp-(61)-{i}-3-1"
        MILP += delta_beta_C[(i, 3)] <= bd3_UB[E_C[i]], f"milp-(61)-{i}-3-2"

    # -------- Constraint (62) --------
    for k in range(1, k_C + 1):
        for i in range(2, t_T + 1):
            for m in {2, 3}:
                MILP += bd_T[(k, i, m)] >= \
                    delta_beta_T[(i, m)] + chi_T[(i, k)] - 1, f"milp-(62)-{k}-{i}-{m}"

    # -------- Constraint (63) --------
    for m in {2, 3}:
        MILP += pulp.lpSum(delta_beta_T[(j, m)] for j in range(2, t_T + 1)) >= \
                pulp.lpSum(bd_T[(k, i, m)] for k in range(1, k_C + 1) 
                                for i in range(2, t_T + 1)), f"milp-(63)-{m}"

    # -------- Constraint (64) --------
    for k in range(1, k_C + 1):
        MILP += pulp.lpSum(bd_T[(k, i, 2)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 2)] + delta_beta_minus[(k, 2)] >= bd2_LB[E_C[k]], f"milp-(64)-{k}-2-1"
        MILP += pulp.lpSum(bd_T[(k, i, 2)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 2)] + delta_beta_minus[(k, 2)] <= bd2_UB[E_C[k]], f"milp-(64)-{k}-2-2"
        MILP += pulp.lpSum(bd_T[(k, i, 3)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 3)] + delta_beta_minus[(k, 3)] >= bd3_LB[E_C[k]], f"milp-(64)-{k}-3-1"
        MILP += pulp.lpSum(bd_T[(k, i, 3)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 3)] + delta_beta_minus[(k, 3)] <= bd3_UB[E_C[k]], f"milp-(64)-{k}-3-2"

    return MILP

# -------- A.8 Descriptors for the Number of Adjacency-configuration --------
def prepare_variables_adjacency_configuration(
    t_C,
    t_C_tilde,
    t_T,
    t_F,
    m_C,
    k_C,
    k_C_tilde,
    c_F,
    n_T,
    n_C,
    n_F,
    delta_i,
    rho,
    Gamma_int_ac,
    Gamma_tilde_ac_C,
    Gamma_tilde_ac_T,
    Gamma_tilde_ac_F,
    Gamma_tilde_ac_CT,
    Gamma_tilde_ac_TC,
    Gamma_tilde_ac_CF,
    Gamma_tilde_ac_TF,
    ac_LB_int,
    ac_UB_int,
    Lambda_int,
    Lambda_ex,
    MAX_CODE_int,
    MAX_CODE_ex
):
    # Define ac_int
    ac_int = {nu: pulp.LpVariable(f"ac_int{nu}", ac_LB_int[nu], ac_UB_int[nu], cat=pulp.LpInteger)
            for nu in Gamma_int_ac}

    # Define ac_C, ac_T, ac_F
    ac_C = {nu: pulp.LpVariable(f"ac_C({nu})", 0, m_C, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_C}
    ac_T = {nu: pulp.LpVariable(f"ac_T({nu})", 0, t_T, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_T}
    ac_F = {nu: pulp.LpVariable(f"ac_F({nu})", 0, t_F, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_F}

    # Define ac_CT, ac_TC
    ac_CT = {nu: pulp.LpVariable(f"ac_CT({nu})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for nu in Gamma_tilde_ac_CT}
    ac_TC = {nu: pulp.LpVariable(f"ac_TC({nu})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for nu in Gamma_tilde_ac_TC} # Gamma_tilde_ac_TC instead of Gamma_tilde_ac_CT

    # Define ac_CF, ac_TF
    ac_CF = {nu: pulp.LpVariable(f"ac_CF({nu})", 0, t_C_tilde, cat=pulp.LpInteger)
            for nu in Gamma_tilde_ac_CF}
    ac_TF = {nu: pulp.LpVariable(f"ac_TF({nu})", 0, t_T, cat=pulp.LpInteger)
            for nu in Gamma_tilde_ac_TF}

    # Define delta_ac_C
    delta_ac_C = {(i, nu): pulp.LpVariable(f"delta_ac_C({i},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(k_C_tilde + 1, m_C + 1) for nu in Gamma_tilde_ac_C}

    # Define delta_ac_T
    delta_ac_T = {(i, nu): pulp.LpVariable(f"delta_ac_T({i},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_T + 1) for nu in Gamma_tilde_ac_T}
    
    # Define delta_ac_F
    delta_ac_F = {(i, nu): pulp.LpVariable(f"delta_ac_F({i},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_F + 1) for nu in Gamma_tilde_ac_F}

    # Define delta_ac_CT, delta_ac_TC
    delta_ac_CT = {(k, nu): pulp.LpVariable(f"delta_ac_CT({k},{nu})", 0, 1, cat=pulp.LpBinary)
            for k in range(1, k_C + 1) for nu in Gamma_tilde_ac_CT}
    delta_ac_TC = {(k, nu): pulp.LpVariable(f"delta_ac_TC({k},{nu})", 0, 1, cat=pulp.LpBinary)
            for k in range(1, k_C + 1) for nu in Gamma_tilde_ac_TC} # Gamma_tilde_ac_TC instead of Gamma_tilde_ac_CT

    # Define delta_ac_CF, delta_ac_TF
    delta_ac_CF = {(c, nu): pulp.LpVariable(f"delta_ac_CF({c},{nu})", 0, 1, cat=pulp.LpBinary)
            for c in range(1, t_C_tilde + 1) for nu in Gamma_tilde_ac_CF}
    delta_ac_TF = {(i, nu): pulp.LpVariable(f"delta_ac_TF({i},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_T + 1) for nu in Gamma_tilde_ac_TF}

    # Define alpha_CT, alpha_TC
    alpha_CT = {k: pulp.LpVariable(f"alpha_CT({k})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                for k in range(1, k_C + 1)}
    alpha_TC = {k: pulp.LpVariable(f"alpha_TC({k})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                for k in range(1, k_C + 1)}
    
    # Define alpha_CF
    alpha_CF = {c: pulp.LpVariable(f"alpha_CF({c})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                for c in range(1, t_C_tilde + 1)}

    # Define alpha_TF
    alpha_TF = {i: pulp.LpVariable(f"alpha_TF({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    
    # Define Delta_ac_C_plus, Delta_ac_C_minus
    Delta_ac_C_plus = {i: pulp.LpVariable(f"Delta_ac_C_plus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(k_C_tilde + 1, m_C + 1)}
    Delta_ac_C_minus = {i: pulp.LpVariable(f"Delta_ac_C_minus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(k_C_tilde + 1, m_C + 1)}
    
    # Define Delta_ac_T_plus, Delta_ac_T_minus
    Delta_ac_T_plus = {i: pulp.LpVariable(f"Delta_ac_T_plus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(2, t_T + 1)}
    Delta_ac_T_minus = {i: pulp.LpVariable(f"Delta_ac_T_minus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(2, t_T + 1)}

    # Define Delta_ac_F_plus, Delta_ac_F_minus
    Delta_ac_F_plus = {i: pulp.LpVariable(f"Delta_ac_F_plus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(2, t_F + 1)}
    Delta_ac_F_minus = {i: pulp.LpVariable(f"Delta_ac_F_minus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(2, t_F + 1)}

    # Define Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus
    # Delta_ac_C_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_C_plus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
    #                 for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    # Delta_ac_C_plus.update(Delta_ac_C_plus_temp)
    # Delta_ac_C_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_C_minus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
    #                 for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    # Delta_ac_C_minus.update(Delta_ac_C_minus_temp)
    # Delta_ac_T_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_T_plus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
    #                 for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    # Delta_ac_T_plus.update(Delta_ac_T_plus_temp)
    # Delta_ac_T_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_T_minus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
    #                 for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    # Delta_ac_T_minus.update(Delta_ac_T_minus_temp)
    # Delta_ac_F_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_F_plus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
    #                 for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    # Delta_ac_F_plus.update(Delta_ac_F_plus_temp)
    # Delta_ac_F_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_F_minus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
    #                 for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    # Delta_ac_F_minus.update(Delta_ac_F_minus_temp)
    
    # Define Delta_ac_CT_plus, Delta_ac_CT_minus
    Delta_ac_CT_plus = {k: pulp.LpVariable(f"Delta_ac_CT_plus({k})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    Delta_ac_CT_minus = {k: pulp.LpVariable(f"Delta_ac_CT_minus({k})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}

    # Define Delta_ac_TC_plus, Delta_ac_TC_minus
    Delta_ac_TC_plus = {k: pulp.LpVariable(f"Delta_ac_TC_plus({k})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    Delta_ac_TC_minus = {k: pulp.LpVariable(f"Delta_ac_TC_minus({k})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    
    # Define Delta_ac_CF_plus, Delta_ac_CF_minus
    Delta_ac_CF_plus = {c: pulp.LpVariable(f"Delta_ac_CF_plus({c})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for c in range(1, t_C_tilde + 1)}
    Delta_ac_CF_minus = {c: pulp.LpVariable(f"Delta_ac_CF_minus({c})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for c in range(1, t_C_tilde + 1)}

    # Define Delta_ac_TF_plus, Delta_ac_TF_minus
    Delta_ac_TF_plus = {i: pulp.LpVariable(f"Delta_ac_TF_plus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1)}
    Delta_ac_TF_minus = {i: pulp.LpVariable(f"Delta_ac_TF_minus({i})", 0, MAX_CODE_int, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1)}

    # # Define ac_co, ac_in, ac_ex
    # ac_co = {nu: pulp.LpVariable(f"ac_co{nu}", ac_LB_co[nu], ac_UB_co[nu], cat=pulp.LpInteger)
    #         for nu in Gamma_co_ac}
    # ac_in = {nu: pulp.LpVariable(f"ac_in{nu}", ac_LB_in[nu], ac_UB_in[nu], cat=pulp.LpInteger)
    #         for nu in Gamma_in_ac}
    # ac_ex = {nu: pulp.LpVariable(f"ac_ex{nu}", ac_LB_ex[nu], ac_UB_ex[nu], cat=pulp.LpInteger)
    #         for nu in Gamma_ex_ac}

    # # Define ac_nc
    # ac_nc = {nu: pulp.LpVariable(f"ac_nc{nu}", 0, m_UB, cat=pulp.LpBinary) for nu in (set(Gamma_in_ac) | set(Gamma_ex_ac))}
    

    return ac_int, ac_C, ac_T, ac_F, ac_CT, ac_TC, ac_CF, ac_TF, \
        delta_ac_C, delta_ac_T, delta_ac_F, \
        alpha_CT, alpha_TC, alpha_CF, alpha_TF, \
        Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, \
        Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus, \
        Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus, Delta_ac_TF_minus, \
        delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF
        # delta_alpha_CT, delta_alpha_TC, delta_alpha_CF, delta_alpha_TF, \


def add_constraints_adjacency_configuration(
    # Model
    MILP,
    # Constants
    t_C,
    t_C_tilde,
    t_T,
    t_F,
    n_C,
    n_T,
    n_F,
    k_C,
    k_C_tilde,
    m_C,
    c_F,
    tail_C,
    head_C,
    Gamma_int_ac,
    Gamma_tilde_ac_C,
    Gamma_tilde_ac_T,
    Gamma_tilde_ac_F,
    Gamma_tilde_ac_CT,
    Gamma_tilde_ac_TC,
    Gamma_tilde_ac_CF,
    Gamma_tilde_ac_TF,
    Gamma_int_ac_less,
    Gamma_int_ac_equal,
    ac_LB_int,
    ac_UB_int,
    Lambda_int,
    Lambda_ex,
    Code_Lambda_int,
    Code_Lambda_ex,
    MAX_CODE_int,
    MAX_CODE_ex,
    rho,
    delta_i,
    # Binary Variables
    delta_alpha_C,
    delta_alpha_T,
    delta_alpha_F,
    delta_beta_C,
    delta_beta_T,
    delta_beta_F,
    delta_beta_plus,
    delta_beta_minus,
    delta_beta_in,
    delta_chi_T,
    delta_chi_F,
    chi_T,
    chi_F,
    delta_ac_C,
    delta_ac_T,
    delta_ac_F,
    delta_ac_CT,
    delta_ac_TC,
    delta_ac_CF,
    delta_ac_TF,
    e_C,
    e_T,
    e_F,
    v_T,
    v_F,
    Delta_ac_C_plus,
    Delta_ac_C_minus,
    Delta_ac_T_plus,
    Delta_ac_T_minus,
    Delta_ac_F_plus,
    Delta_ac_F_minus,
    Delta_ac_CT_plus,
    Delta_ac_CT_minus,
    Delta_ac_TC_plus,
    Delta_ac_TC_minus,
    Delta_ac_CF_plus,
    Delta_ac_CF_minus,
    Delta_ac_TF_plus,
    Delta_ac_TF_minus,
    # Integer Variables
    alpha_C,
    alpha_T,
    alpha_F,
    alpha_CT,
    alpha_TC,
    alpha_CF,
    alpha_TF,
    beta_C,
    beta_T,
    beta_F,
    beta_plus,
    beta_minus,
    beta_in,
    ac_C,
    ac_T,
    ac_F,
    ac_CT,
    ac_TC,
    ac_CF,
    ac_TF,
    ac_int
):
    # -------- Constraint (65) --------
    for nu in Gamma_int_ac:
        if nu not in Gamma_tilde_ac_C:
            MILP += ac_C[nu] == 0, f"milp-(65)-1-{nu}"
        if nu not in Gamma_tilde_ac_T:
            MILP += ac_T[nu] == 0, f"milp-(65)-2-{nu}"
        if nu not in Gamma_tilde_ac_F:
            MILP += ac_F[nu] == 0, f"milp-(65)-3-{nu}"
        if nu not in Gamma_tilde_ac_CT:
            MILP += ac_CT[nu] == 0, f"milp-(65)-4-{nu}"
        if nu not in Gamma_tilde_ac_TC:
            MILP += ac_TC[nu] == 0, f"milp-(65)-5-{nu}"
        if nu not in Gamma_tilde_ac_CF:
            MILP += ac_CF[nu] == 0, f"milp-(65)-6-{nu}"
        if nu not in Gamma_tilde_ac_TF:
            MILP += ac_TF[nu] == 0, f"milp-(65)-7-{nu}"

    # -------- Constraint (66) --------
    for m in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum(ac_C[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
                pulp.lpSum(delta_beta_C[(i, m)] for i in range(k_C_tilde + 1, m_C + 1)), f"milp-(66)-1-{m}"
        MILP += pulp.lpSum(ac_T[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
                pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)), f"milp-(66)-2-{m}"
        MILP += pulp.lpSum(ac_F[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
                pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)), f"milp-(66)-3-{m}"
        MILP += pulp.lpSum(ac_CT[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
                pulp.lpSum(delta_beta_plus[(k, m)] for k in range(1, k_C + 1)), f"milp-(66)-4-{m}"
        MILP += pulp.lpSum(ac_TC[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
                pulp.lpSum(delta_beta_minus[(k, m)] for k in range(1, k_C + 1)), f"milp-(66)-5-{m}"
        MILP += pulp.lpSum(ac_CF[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)), f"milp-(66)-6-{m}"
        MILP += pulp.lpSum(ac_TF[(a, b, m1)] for (a, b, m1) in Gamma_int_ac if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-(66)-7-{m}"

    # -------- Constraint (67) --------
    for i in range(k_C_tilde + 1, m_C + 1):
        MILP += pulp.lpSum(m * delta_ac_C[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_C) == beta_C[i], f"milp-(67)-1-{i}"
        MILP += Delta_ac_C_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_C[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_C) == alpha_C[(tail_C[i], 0)], f"milp-(67)-2-{i}"
        MILP += Delta_ac_C_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_C[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_C) == alpha_C[(head_C[i], 0)], f"milp-(67)-3-{i}"
        MILP += Delta_ac_C_plus[i] + Delta_ac_C_minus[i] <= 2 * MAX_CODE_int * (1 - e_C[i]), f"milp-(67)-4-{i}"
    for nu in Gamma_tilde_ac_C:
        MILP += pulp.lpSum(delta_ac_C[(i, nu)] 
                for i in range(k_C_tilde + 1, m_C + 1)) == ac_C[nu], f"milp-(67)-5-{nu}"

    # -------- Constraint (68) --------
    for i in range(2, t_T + 1):
        MILP += pulp.lpSum(m * delta_ac_T[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_T) == beta_T[i], f"milp-(68)-1-{i}"
        MILP += Delta_ac_T_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_T[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_T) == alpha_T[(i - 1, 0)], f"milp-(68)-2-{i}"
        MILP += Delta_ac_T_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_T[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_T) == alpha_T[(i, 0)], f"milp-(68)-3-{i}"
        MILP += Delta_ac_T_plus[i] + Delta_ac_T_minus[i] <= 2 * MAX_CODE_int * (1 - e_T[i]), f"milp-(68)-4-{i}"
    for nu in Gamma_tilde_ac_T:
        MILP += pulp.lpSum(delta_ac_T[(i, nu)] 
                for i in range(2, t_T + 1)) == ac_T[nu], f"milp-(68)-5-{nu}"

    # -------- Constraint (69) --------
    for i in range(2, t_F + 1):
        MILP += pulp.lpSum(m * delta_ac_F[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_F) == beta_F[i], f"milp-(69)-1-{i}"
        MILP += Delta_ac_F_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_F[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[(i - 1, 0)], f"milp-(69)-2-{i}"
        MILP += Delta_ac_F_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_F[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[(i, 0)], f"milp-(69)-3-{i}"
        MILP += Delta_ac_F_plus[i] + Delta_ac_F_minus[i] <= 2 * MAX_CODE_ex * (1 - e_F[i]), f"milp-(69)-4-{i}"
    for nu in Gamma_tilde_ac_F:
        MILP += pulp.lpSum(delta_ac_F[(i, nu)] 
                for i in range(2, t_F + 1)) == ac_F[nu], f"milp-(69)-5-{nu}"

    # -------- Constraint (70) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += alpha_T[(i, 0)] + MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i]) >= \
                    alpha_CT[k], f"milp-(70)-1-{k}-{i}"
            MILP += alpha_CT[k] >= alpha_T[(i, 0)] - \
                    MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i]), f"milp-(70)-2-{k}-{i}"
        MILP += pulp.lpSum(m * delta_ac_CT[(k, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_CT) == beta_plus[k], f"milp-(70)-3-{k}"
        MILP += Delta_ac_CT_plus[k] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_CT[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CT) == alpha_C[(tail_C[k], 0)], f"milp-(70)-4-{k}"
        MILP += Delta_ac_CT_minus[k] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_CT[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CT) == alpha_CT[k], f"milp-(70)-5-{k}"
        MILP += Delta_ac_CT_plus[k] + Delta_ac_CT_minus[k] <= 2 * MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-(70)-6-{k}"
    for nu in Gamma_tilde_ac_CT:
        MILP += pulp.lpSum(delta_ac_CT[(k, nu)] for k in range(1, k_C + 1)) == ac_CT[nu], f"milp-(70)-7-{nu}"

    # -------- Constraint (71) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += alpha_T[(i, 0)] + MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i + 1]) >= \
                    alpha_TC[k], f"milp-(71)-1-{k}-{i}"
            MILP += alpha_TC[k] >= alpha_T[(i, 0)] - \
                    MAX_CODE_int * (1 - chi_T[(i, k)] + e_T[i + 1]), f"milp-(71)-2-{k}-{i}"
        MILP += pulp.lpSum(m * delta_ac_TC[(k, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_TC) == beta_minus[k], f"milp-(71)-3-{k}"
        MILP += Delta_ac_TC_plus[k] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_TC[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TC) == alpha_TC[k], f"milp-(71)-4-{k}"
        MILP += Delta_ac_TC_minus[k] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_TC[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TC) == alpha_C[(head_C[k], 0)], f"milp-(71)-5-{k}"
        MILP += Delta_ac_TC_plus[k] + Delta_ac_TC_minus[k] <= 2 * MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-(71)-6-{k}"
    for nu in Gamma_tilde_ac_TC:
        MILP += pulp.lpSum(delta_ac_TC[(k, nu)] for k in range(1, k_C + 1)) == ac_TC[nu], f"milp-(71)-7-{nu}"

    # -------- Constraint (72) --------
    for c in range(1, t_C_tilde + 1):
        for i in range(1, t_F + 1):
            MILP += alpha_F[(i, 0)] + MAX_CODE_int * (1 - chi_F[(i, c)] + e_F[i]) >= \
                    alpha_CF[c], f"milp-(72)-1-{c}-{i}"
            MILP += alpha_CF[c] >= alpha_F[(i, 0)] - \
                    MAX_CODE_int * (1 - chi_F[(i, c)] + e_F[i]), f"milp-(72)-2-{c}-{i}"
        MILP += pulp.lpSum(m * delta_ac_CF[(c, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_CF) == beta_in[c], f"milp-(72)-3-{c}"
        MILP += Delta_ac_CF_plus[c] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_CF[(c, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CF) == alpha_C[(c, 0)], f"milp-(72)-4-{c}"
        MILP += Delta_ac_CF_minus[c] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_CF[(c, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CF) == alpha_CF[c], f"milp-(72)-5-{c}"
        MILP += Delta_ac_CF_plus[c] + Delta_ac_CF_minus[c] <= \
                2 * max(MAX_CODE_int, MAX_CODE_int) * (1 - delta_chi_F[c]), f"milp-(72)-6-{c}"
    for nu in Gamma_tilde_ac_CF:
        MILP += pulp.lpSum(delta_ac_CF[(c, nu)] for c in range(1, t_C_tilde + 1)) == ac_CF[nu], f"milp-(72)-7-{nu}"

    # -------- Constraint (73) --------
    for i in range(1, t_T + 1):
        for j in range(1, t_F + 1):
            MILP += alpha_F[(j, 0)] + MAX_CODE_int * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]) >= \
                    alpha_TF[i], f"milp-(73)-1-{i}-{j}"
            MILP += alpha_TF[i] >= alpha_F[(j, 0)] - \
                    MAX_CODE_int * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]), f"milp-(73)-2-{i}-{j}"
        MILP += pulp.lpSum(m * delta_ac_TF[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_TF) == beta_in[i + t_C_tilde], f"milp-(73)-3-{i}"
        MILP += Delta_ac_TF_plus[i] + pulp.lpSum(Code_Lambda_int[a] * delta_ac_TF[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TF) == alpha_T[(i, 0)], f"milp-(73)-4-{i}"
        MILP += Delta_ac_TF_minus[i] + pulp.lpSum(Code_Lambda_int[b] * delta_ac_TF[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TF) == alpha_TF[i], f"milp-(73)-5-{i}"
        MILP += Delta_ac_TF_plus[i] + Delta_ac_TF_minus[i] <= \
                2 * max(MAX_CODE_int, MAX_CODE_int) * (1 - delta_chi_F[i + t_C_tilde]), f"milp-(73)-6-{i}"
    for nu in Gamma_tilde_ac_TF:
        MILP += pulp.lpSum(delta_ac_TF[(i, nu)] for i in range(1, t_T + 1)) == ac_TF[nu], f"milp-(73)-7-{nu}"

    # -------- Constraint (74) --------
    for (a, b, m) in Gamma_int_ac_less:
        MILP += ac_C[(a, b, m)] + ac_C[(b, a, m)] + \
                ac_T[(a, b, m)] + ac_T[(b, a, m)] + \
                ac_F[(a, b, m)] + ac_F[(b, a, m)] + \
                ac_CT[(a, b, m)] + ac_CT[(b, a, m)] + \
                ac_TC[(a, b, m)] + ac_TC[(b, a, m)] + \
                ac_CF[(a, b, m)] + ac_CF[(b, a, m)] + \
                ac_TF[(a, b, m)] + ac_TF[(b, a, m)] == ac_int[(a, b, m)], f"milp-(74)-1-{(a, b, m)}"
    for nu in Gamma_int_ac_equal:
        MILP += ac_C[nu] + ac_T[nu] + ac_F[nu] + \
                ac_CT[nu] + ac_TC[nu] + ac_CF[nu] + ac_TF[nu] == ac_int[nu], f"milp-(74)-2-{nu}"

    return MILP

# -------- A.9 Descriptor for the Number of Chemical Symbols --------
def prepare_variables_chemical_symbols(
    t_C,
    t_T,
    t_F,
    n_C,
    n_T,
    n_F,
    delta_i,
    n_star,
    n_LB_int,
    n_UB_int,
    Lambda_dg_int,
    epsilon_dg
):
    # Define ns_int
    ns_int = {mu: pulp.LpVariable(f"ns_int({mu})", 0, n_UB_int, cat=pulp.LpInteger) for mu in Lambda_dg_int}

    # Define delta_ns_C, delta_ns_T, delta_ns_F
    delta_ns_C = {(i, 0, mu): pulp.LpVariable(f"delta_ns_C({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_C + 1) for mu in (Lambda_dg_int + [epsilon_dg])}
    delta_ns_T = {(i, 0, mu): pulp.LpVariable(f"delta_ns_T({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_T + 1) for mu in (Lambda_dg_int + [epsilon_dg])}
    delta_ns_F = {(i, 0, mu): pulp.LpVariable(f"delta_ns_F({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_F + 1) for mu in (Lambda_dg_int + [epsilon_dg])}
    
    return ns_int, delta_ns_C, delta_ns_T, delta_ns_F
    
def add_constraints_chemical_symbols(
    # Model
    MILP,
    # Constants
    t_C,
    t_T,
    t_F,
    n_C,
    n_T,
    n_F,
    Lambda_dg_int,
    Lambda_tilde_dg_C,
    Lambda_tilde_dg_T,
    Lambda_tilde_dg_F,
    Code_Lambda,
    Code_Lambda_int,
    Code_Lambda_ex,
    epsilon_dg,
    delta_i,
    rho,
    # Binary Variables
    delta_ns_C,
    delta_ns_T,
    delta_ns_F,
    # Integer Variables
    alpha_C,
    alpha_T,
    alpha_F,
    deg_C,
    deg_T,
    deg_F,
    ns_int
):
    # -------- Constraint (75) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_ns_C[(i, 0, mu)] 
                for mu in (Lambda_tilde_dg_C + [epsilon_dg])) == 1, f"milp-(75)-C-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_int[a] * delta_ns_C[(i, 0, (a, d))] 
                for (a, d) in Lambda_tilde_dg_C) == alpha_C[(i, 0)], f"milp-(75)-C-2-{i}"
        MILP += pulp.lpSum(d * delta_ns_C[(i, 0, (a, d))]
                for (a, d) in Lambda_tilde_dg_C) == deg_C[(i, 0)], f"milp-(75)-C-3-{i}"
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(delta_ns_T[(i, 0, mu)] 
                for mu in (Lambda_tilde_dg_T + [epsilon_dg])) == 1, f"milp-(75)-T-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_int[a] * delta_ns_T[(i, 0, (a, d))] 
                for (a, d) in Lambda_tilde_dg_T) == alpha_T[(i, 0)], f"milp-(75)-T-2-{i}"
        MILP += pulp.lpSum(d * delta_ns_T[(i, 0, (a, d))]
                for (a, d) in Lambda_tilde_dg_T) == deg_T[(i, 0)], f"milp-(75)-T-3-{i}"
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_ns_F[(i, 0, mu)] 
                for mu in (Lambda_tilde_dg_F + [epsilon_dg])) == 1, f"milp-(75)-F-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_int[a] * delta_ns_F[(i, 0, (a, d))] 
                for (a, d) in Lambda_tilde_dg_F) == alpha_F[(i, 0)], f"milp-(75)-F-2-{i}"
        MILP += pulp.lpSum(d * delta_ns_F[(i, 0, (a, d))]
                for (a, d) in Lambda_tilde_dg_F) == deg_F[(i, 0)], f"milp-(75)-F-3-{i}"
    
    # -------- Constraint (76) --------
    for mu in Lambda_dg_int:
        MILP += pulp.lpSum(delta_ns_C[(i, 0, mu)] for i in range(1, t_C + 1)) + \
                pulp.lpSum(delta_ns_T[(i, 0, mu)] for i in range(1, t_T + 1)) + \
                pulp.lpSum(delta_ns_F[(i, 0, mu)] for i in range(1, t_F + 1)) == ns_int[mu], f"milp-(76)-{mu}"

    return MILP

# -------- A.10 Descriptor for the Number of Edge-configurations --------
def prepare_variables_edge_configuration(
    t_C,
    t_C_tilde,
    t_T,
    t_F,
    m_C,
    k_C,
    k_C_tilde,
    c_F,
    n_T,
    n_C,
    n_F,
    delta_i,
    rho,
    Gamma_int,
    Gamma_tilde_ec_C,
    Gamma_tilde_ec_T,
    Gamma_tilde_ec_F,
    Gamma_tilde_ec_CT,
    Gamma_tilde_ec_TC,
    Gamma_tilde_ec_CF,
    Gamma_tilde_ec_TF,
    ec_LB_int,
    ec_UB_int
):
    # Define ec_int
    ec_int = {gamma: pulp.LpVariable(f"ec_int{gamma}", ec_LB_int[gamma], ec_UB_int[gamma], cat=pulp.LpInteger)
            for gamma in Gamma_int}

    # Define ec_C, ec_T, ec_F
    ec_C = {gamma: pulp.LpVariable(f"ec_C({gamma})", 0, m_C, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_C}
    ec_T = {gamma: pulp.LpVariable(f"ec_T({gamma})", 0, t_T, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_T}
    ec_F = {gamma: pulp.LpVariable(f"ec_F({gamma})", 0, t_F, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_F}

    # Define ec_CT, ec_TC
    ec_CT = {gamma: pulp.LpVariable(f"ec_CT({gamma})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for gamma in Gamma_tilde_ec_CT}
    ec_TC = {gamma: pulp.LpVariable(f"ec_TC({gamma})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for gamma in Gamma_tilde_ec_TC} # Gamma_tilde_ec_TC instead of Gamma_tilde_ec_CT

    # Define ec_CF, ec_TF
    ec_CF = {gamma: pulp.LpVariable(f"ec_CF({gamma})", 0, t_C_tilde, cat=pulp.LpInteger)
            for gamma in Gamma_tilde_ec_CF}
    ec_TF = {gamma: pulp.LpVariable(f"ec_TF({gamma})", 0, t_T, cat=pulp.LpInteger)
            for gamma in Gamma_tilde_ec_TF}

    # Define delta_ec_C
    delta_ec_C = {(i, gamma): pulp.LpVariable(f"delta_ec_C({i},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(k_C_tilde + 1, m_C + 1) for gamma in Gamma_tilde_ec_C}

    # Define delta_ec_T
    delta_ec_T = {(i, gamma): pulp.LpVariable(f"delta_ec_T({i},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_T + 1) for gamma in Gamma_tilde_ec_T}
    
    # Define delta_ec_F
    delta_ec_F = {(i, gamma): pulp.LpVariable(f"delta_ec_F({i},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_F + 1) for gamma in Gamma_tilde_ec_F}

    # Define delta_ec_CT, delta_ec_TC
    delta_ec_CT = {(k, gamma): pulp.LpVariable(f"delta_ec_CT({k},{gamma})", 0, 1, cat=pulp.LpBinary)
            for k in range(1, k_C + 1) for gamma in Gamma_tilde_ec_CT}
    delta_ec_TC = {(k, gamma): pulp.LpVariable(f"delta_ec_TC({k},{gamma})", 0, 1, cat=pulp.LpBinary)
            for k in range(1, k_C + 1) for gamma in Gamma_tilde_ec_TC} # Gamma_tilde_ec_TC instead of Gamma_tilde_ec_CT

    # Define delta_ec_CF, delta_ec_TF
    delta_ec_CF = {(c, gamma): pulp.LpVariable(f"delta_ec_CF({c},{gamma})", 0, 1, cat=pulp.LpBinary)
            for c in range(1, t_C_tilde + 1) for gamma in Gamma_tilde_ec_CF}
    delta_ec_TF = {(i, gamma): pulp.LpVariable(f"delta_ec_TF({i},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_T + 1) for gamma in Gamma_tilde_ec_TF}
    
    # Define deg_T_CT, deg_T_TC
    deg_T_CT = {k: pulp.LpVariable(f"deg_T_CT({k})", 0, MAX_VAL, cat=pulp.LpInteger)
                for k in range(1, k_C + 1)}
    deg_T_TC = {k: pulp.LpVariable(f"deg_T_TC({k})", 0, MAX_VAL, cat=pulp.LpInteger)
                for k in range(1, k_C + 1)}

    # Define deg_F_CF
    deg_F_CF = {c: pulp.LpVariable(f"deg_F_CF({c})", 0, MAX_VAL, cat=pulp.LpInteger)
                for c in range(1, t_C_tilde + 1)}

    # Define deg_F_TF
    deg_F_TF = {i: pulp.LpVariable(f"deg_F_TF({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}

    # Define Delta_ec_C_plus, Delta_ec_C_minus
    Delta_ec_C_plus = {i: pulp.LpVariable(f"Delta_ec_C_plus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(k_C_tilde + 1, m_C + 1)}
    Delta_ec_C_minus = {i: pulp.LpVariable(f"Delta_ec_C_minus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(k_C_tilde + 1, m_C + 1)}
    
    # Define Delta_ec_T_plus, Delta_ec_T_minus
    Delta_ec_T_plus = {i: pulp.LpVariable(f"Delta_ec_T_plus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(2, t_T + 1)}
    Delta_ec_T_minus = {i: pulp.LpVariable(f"Delta_ec_T_minus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(2, t_T + 1)}

    # Define Delta_ec_F_plus, Delta_ec_F_minus
    Delta_ec_F_plus = {i: pulp.LpVariable(f"Delta_ec_F_plus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(2, t_F + 1)}
    Delta_ec_F_minus = {i: pulp.LpVariable(f"Delta_ec_F_minus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(2, t_F + 1)}

    # # Define Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, Delta_ec_F_plus, Delta_ec_F_minus
    # Delta_ec_C_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_C_plus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
    #                 for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    # Delta_ec_C_plus.update(Delta_ec_C_plus_temp)
    # Delta_ec_C_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_C_minus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
    #                 for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    # Delta_ec_C_minus.update(Delta_ec_C_minus_temp)
    # Delta_ec_T_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_T_plus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
    #                 for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    # Delta_ec_T_plus.update(Delta_ec_T_plus_temp)
    # Delta_ec_T_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_T_minus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
    #                 for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    # Delta_ec_T_minus.update(Delta_ec_T_minus_temp)
    # Delta_ec_F_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_F_plus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
    #                 for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    # Delta_ec_F_plus.update(Delta_ec_F_plus_temp)
    # Delta_ec_F_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_F_minus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
    #                 for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    # Delta_ec_F_minus.update(Delta_ec_F_minus_temp)
    
    # Define Delta_ec_CT_plus, Delta_ec_CT_minus
    Delta_ec_CT_plus = {k: pulp.LpVariable(f"Delta_ec_CT_plus({k})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    Delta_ec_CT_minus = {k: pulp.LpVariable(f"Delta_ec_CT_minus({k})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}

    # Define Delta_ec_TC_plus, Delta_ec_TC_minus
    Delta_ec_TC_plus = {k: pulp.LpVariable(f"Delta_ec_TC_plus({k})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    Delta_ec_TC_minus = {k: pulp.LpVariable(f"Delta_ec_TC_minus({k})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    
    # Define Delta_ec_CF_plus, Delta_ec_CF_minus
    Delta_ec_CF_plus = {c: pulp.LpVariable(f"Delta_ec_CF_plus({c})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for c in range(1, t_C_tilde + 1)}
    Delta_ec_CF_minus = {c: pulp.LpVariable(f"Delta_ec_CF_minus({c})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for c in range(1, t_C_tilde + 1)}

    # Define Delta_ec_TF_plus, Delta_ec_TF_minus
    Delta_ec_TF_plus = {i: pulp.LpVariable(f"Delta_ec_TF_plus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1)}
    Delta_ec_TF_minus = {i: pulp.LpVariable(f"Delta_ec_TF_minus({i})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1)}
    
    return ec_int, ec_C, ec_T, ec_F, ec_CT, ec_TC, ec_CF, ec_TF, \
        delta_ec_C, delta_ec_T, delta_ec_F, \
        delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF, \
        deg_T_CT, deg_T_TC, deg_F_CF, deg_F_TF, \
        Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, Delta_ec_F_plus, Delta_ec_F_minus, \
        Delta_ec_CT_plus, Delta_ec_CT_minus, Delta_ec_TC_plus, Delta_ec_TC_minus, \
        Delta_ec_CF_plus, Delta_ec_CF_minus, Delta_ec_TF_plus, Delta_ec_TF_minus

def add_constraints_edge_configuration(
    # Model
    MILP,
    # Constants
    t_C,
    t_C_tilde,
    t_T,
    t_F,
    n_C,
    n_T,
    n_F,
    k_C,
    k_C_tilde,
    m_C,
    c_F,
    tail_C,
    head_C,
    Gamma_int,
    Gamma_int_less,
    Gamma_int_equal,
    Gamma_tilde_ec_C,
    Gamma_tilde_ec_T,
    Gamma_tilde_ec_F,
    Gamma_tilde_ec_CT,
    Gamma_tilde_ec_TC,
    Gamma_tilde_ec_CF,
    Gamma_tilde_ec_TF,
    Gamma_tilde_ac_C,
    Gamma_tilde_ac_T,
    Gamma_tilde_ac_F,
    Gamma_tilde_ac_CT,
    Gamma_tilde_ac_TC,
    Gamma_tilde_ac_CF,
    Gamma_tilde_ac_TF,
    rho,
    delta_i,
    Code_Gamma_ac_int,
    Code_Gamma_ec_int,
    # Binary Variables
    delta_dg_C,
    delta_dg_T,
    delta_dg_F,
    delta_chi_T,
    delta_chi_F,
    chi_T,
    chi_F,
    delta_beta_C,
    delta_beta_T,
    delta_beta_F,
    delta_beta_plus,
    delta_beta_minus,
    delta_beta_in,
    delta_ac_C,
    delta_ac_T,
    delta_ac_F,
    delta_ac_CT,
    delta_ac_TC,
    delta_ac_CF,
    delta_ac_TF,
    delta_ec_C,
    delta_ec_T,
    delta_ec_F,
    delta_ec_CT,
    delta_ec_TC,
    delta_ec_CF,
    delta_ec_TF,
    e_C,
    e_T,
    e_F,
    v_T,
    v_F,
    # Integer Variables
    deg_C,
    deg_T,
    deg_F,
    deg_T_CT,
    deg_T_TC,
    deg_F_CF,
    deg_F_TF,
    ec_C,
    ec_T,
    ec_F,
    ec_CT,
    ec_TC,
    ec_CF,
    ec_TF,
    ec_int,
    Delta_ec_C_plus,
    Delta_ec_C_minus,
    Delta_ec_T_plus,
    Delta_ec_T_minus,
    Delta_ec_F_plus,
    Delta_ec_F_minus,
    Delta_ec_CT_plus,
    Delta_ec_CT_minus,
    Delta_ec_TC_plus,
    Delta_ec_TC_minus,
    Delta_ec_CF_plus,
    Delta_ec_CF_minus,
    Delta_ec_TF_plus,
    Delta_ec_TF_minus
):
    # -------- Constraint (77) --------
    for gamma in Gamma_int:
        if gamma not in Gamma_tilde_ec_C:
            MILP += ec_C[gamma] == 0, f"milp-(77)-1-{gamma}"
        if gamma not in Gamma_tilde_ec_T:
            MILP += ec_T[gamma] == 0, f"milp-(77)-2-{gamma}"
        if gamma not in Gamma_tilde_ec_F:
            MILP += ec_F[gamma] == 0, f"milp-(77)-3-{gamma}"
        if gamma not in Gamma_tilde_ec_CT:
            MILP += ec_CT[gamma] == 0, f"milp-(77)-4-{gamma}"
        if gamma not in Gamma_tilde_ec_TC:
            MILP += ec_TC[gamma] == 0, f"milp-(77)-5-{gamma}"
        if gamma not in Gamma_tilde_ec_CF:
            MILP += ec_CF[gamma] == 0, f"milp-(77)-6-{gamma}"
        if gamma not in Gamma_tilde_ec_TF:
            MILP += ec_TF[gamma] == 0, f"milp-(77)-7-{gamma}"

    # -------- Constraint (78) --------
    for m in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum(ec_C[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
                pulp.lpSum(delta_beta_C[(i, m)] for i in range(k_C_tilde + 1, m_C + 1)), f"milp-(78)-1-{m}"
        MILP += pulp.lpSum(ec_T[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
                pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)), f"milp-(78)-2-{m}"
        MILP += pulp.lpSum(ec_F[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
                pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)), f"milp-(78)-3-{m}"
        MILP += pulp.lpSum(ec_CT[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
                pulp.lpSum(delta_beta_plus[(k, m)] for k in range(1, k_C + 1)), f"milp-(78)-4-{m}"
        MILP += pulp.lpSum(ec_TC[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
                pulp.lpSum(delta_beta_minus[(k, m)] for k in range(1, k_C + 1)), f"milp-(78)-5-{m}"
        MILP += pulp.lpSum(ec_CF[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)), f"milp-(78)-6-{m}"
        MILP += pulp.lpSum(ec_TF[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_int if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-(78)-8-{m}"

    # -------- Constraint (79) --------
    for i in range(k_C_tilde + 1, m_C + 1):
        MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_C[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_C) == \
                pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_C[(i, nu)] 
                    for nu in Gamma_tilde_ac_C), f"milp-(79)-1-{i}"
        MILP += Delta_ec_C_plus[i] + pulp.lpSum(d * delta_ec_C[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_C) == deg_C[(tail_C[i], 0)], f"milp-(79)-2-{i}"
        MILP += Delta_ec_C_minus[i] + pulp.lpSum(d * delta_ec_C[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_C) == deg_C[(head_C[i], 0)], f"milp-(79)-3-{i}"
        MILP += Delta_ec_C_plus[i] + Delta_ec_C_minus[i] <= 8 * (1 - e_C[i]), f"milp-(79)-4-{i}"
    for gamma in Gamma_tilde_ec_C:
        MILP += pulp.lpSum(delta_ec_C[(i, gamma)] 
                for i in range(k_C_tilde + 1, m_C + 1)) == ec_C[gamma], f"milp-(79)-5-{gamma}"
    
    # -------- Constraint (80) --------
    for i in range(2, t_T + 1):
        MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_T[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_T) == \
                pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_T[(i, nu)] 
                    for nu in Gamma_tilde_ac_T), f"milp-(80)-1-{i}"
        MILP += Delta_ec_T_plus[i] + pulp.lpSum(d * delta_ec_T[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_T) == deg_T[(i - 1, 0)], f"milp-(80)-2-{i}"
        MILP += Delta_ec_T_minus[i] + pulp.lpSum(d * delta_ec_T[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_T) == deg_T[(i, 0)], f"milp-(80)-3-{i}"
        MILP += Delta_ec_T_plus[i] + Delta_ec_T_minus[i] <= 8 * (1 - e_T[i]), f"milp-(80)-4-{i}"
    for gamma in Gamma_tilde_ec_T:
        MILP += pulp.lpSum(delta_ec_T[(i, gamma)] 
                for i in range(2, t_T + 1)) == ec_T[gamma], f"milp-(80)-5-{gamma}"
  
    # -------- Constraint (81) --------
    for i in range(2, t_F + 1):
        MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_F[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_F) == \
                pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_F[(i, nu)] 
                    for nu in Gamma_tilde_ac_F), f"milp-(81)-1-{i}"
        MILP += Delta_ec_F_plus[i] + pulp.lpSum(d * delta_ec_F[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_F) == deg_F[(i - 1, 0)], f"milp-(81)-2-{i}"
        MILP += Delta_ec_F_minus[i] + pulp.lpSum(d * delta_ec_F[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_F) == deg_F[(i, 0)], f"milp-(81)-3-{i}"
        MILP += Delta_ec_F_plus[i] + Delta_ec_F_minus[i] <= 8 * (1 - e_F[i]), f"milp-(81)-4-{i}"
    for gamma in Gamma_tilde_ec_F:
        MILP += pulp.lpSum(delta_ec_F[(i, gamma)] 
                for i in range(2, t_F + 1)) == ec_F[gamma], f"milp-(81)-5-{gamma}"

    # -------- Constraint (82) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += deg_T[(i, 0)] + MAX_VAL * (1 - chi_T[(i, k)] + e_T[i]) >= \
                    deg_T_CT[k], f"milp-(82)-1-{k}-{i}"
            MILP += deg_T_CT[k] >= deg_T[(i, 0)] - \
                    MAX_VAL * (1 - chi_T[(i, k)] + e_T[i]), f"milp-(82)-2-{k}-{i}"
        MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_CT[(k, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_CT) == \
                pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_CT[(k, nu)] 
                    for nu in Gamma_tilde_ac_CT), f"milp-(82)-3-{k}"
        MILP += Delta_ec_CT_plus[k] + pulp.lpSum(d * delta_ec_CT[(k, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_CT) == deg_C[(tail_C[k], 0)], f"milp-(82)-4-{k}"
        MILP += Delta_ec_CT_minus[k] + pulp.lpSum(d * delta_ec_CT[(k, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_CT) == deg_T_CT[k], f"milp-(82)-5-{k}"
        MILP += Delta_ec_CT_plus[k] + Delta_ec_CT_minus[k] <= 8 * (1 - delta_chi_T[k]), f"milp-(82)-6-{k}"
    for gamma in Gamma_tilde_ec_CT:
        MILP += pulp.lpSum(delta_ec_CT[(k, gamma)] for k in range(1, k_C + 1)) == ec_CT[gamma], f"milp-(82)-7-{gamma}"

    # -------- Constraint (83) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += deg_T[(i, 0)] + MAX_VAL * (1 - chi_T[(i, k)] + e_T[i + 1]) >= \
                    deg_T_TC[k], f"milp-(83)-1-{k}-{i}"
            MILP += deg_T_TC[k] >= deg_T[(i, 0)] - \
                    MAX_VAL * (1 - chi_T[(i, k)] + e_T[i + 1]), f"milp-(83)-2-{k}-{i}"
        MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_TC[(k, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_TC) == \
                pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_TC[(k, nu)] 
                    for nu in Gamma_tilde_ac_TC), f"milp-(83)-3-{k}"
        MILP += Delta_ec_TC_plus[k] + pulp.lpSum(d * delta_ec_TC[(k, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_TC) == deg_T_TC[k], f"milp-(83)-4-{k}"
        MILP += Delta_ec_TC_minus[k] + pulp.lpSum(d * delta_ec_TC[(k, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_TC) == deg_C[(head_C[k], 0)], f"milp-(83)-5-{k}"
        MILP += Delta_ec_TC_plus[k] + Delta_ec_TC_minus[k] <= 8 * (1 - delta_chi_T[k]), f"milp-(83)-6-{k}"
    for gamma in Gamma_tilde_ec_TC:
        MILP += pulp.lpSum(delta_ec_TC[(k, gamma)] for k in range(1, k_C + 1)) == ec_TC[gamma], f"milp-(83)-7-{gamma}"

    # -------- Constraint (84) --------
    for c in range(1, t_C_tilde + 1):
        for i in range(1, t_F + 1):
            MILP += deg_F[(i, 0)] + MAX_VAL * (1 - chi_F[(i, c)] + e_F[i]) >= \
                    deg_F_CF[c], f"milp-(84)-1-{c}-{i}"
            MILP += deg_F_CF[c] >= deg_F[(i, 0)] - \
                    MAX_VAL * (1 - chi_F[(i, c)] + e_F[i]), f"milp-(84)-2-{c}-{i}"
        MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_CF[(c, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_CF) == \
                pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_CF[(c, nu)] 
                    for nu in Gamma_tilde_ac_CF), f"milp-(84)-3-{c}"
        MILP += Delta_ec_CF_plus[c] + pulp.lpSum(d * delta_ec_CF[(c, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_CF) == deg_C[(c, 0)], f"milp-(84)-4-{c}"
        MILP += Delta_ec_CF_minus[c] + pulp.lpSum(d * delta_ec_CF[(c, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_CF) == deg_F_CF[c], f"milp-(84)-5-{c}"
        MILP += Delta_ec_CF_plus[c] + Delta_ec_CF_minus[c] <= 8 * (1 - delta_chi_F[c]), f"milp-(84)-6-{c}"
    for gamma in Gamma_tilde_ec_CF:
        MILP += pulp.lpSum(delta_ec_CF[(c, gamma)] for c in range(1, t_C_tilde + 1)) == ec_CF[gamma], f"milp-(84)-7-{gamma}"

    # -------- Constraint (85) --------
    for i in range(1, t_T + 1):
        for j in range(1, t_F + 1):
            MILP += deg_F[(j, 0)] + MAX_VAL * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]) >= \
                    deg_F_TF[i], f"milp-(85)-1-{i}-{j}"
            MILP += deg_F_TF[i] >= deg_F[(j, 0)] - \
                    MAX_VAL * (1 - chi_F[(j, i + t_C_tilde)] + e_F[j]), f"milp-(85)-2-{i}-{j}"
        MILP += pulp.lpSum(Code_Gamma_ac_int[(a, b, m)] * delta_ec_TF[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_TF) == \
                pulp.lpSum(Code_Gamma_ac_int[nu] * delta_ac_TF[(i, nu)] 
                    for nu in Gamma_tilde_ac_TF), f"milp-(85)-3-{i}"
        MILP += Delta_ec_TF_plus[i] + pulp.lpSum(d * delta_ec_TF[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_TF) == deg_T[(i, 0)], f"milp-(85)-4-{i}"
        MILP += Delta_ec_TF_minus[i] + pulp.lpSum(d * delta_ec_TF[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_TF) == deg_F_TF[i], f"milp-(85)-5-{i}"
        MILP += Delta_ec_TF_plus[i] + Delta_ec_TF_minus[i] <= 8 * (1 - delta_chi_F[i + t_C_tilde]), f"milp-(85)-6-{i}"
    for gamma in Gamma_tilde_ec_TF:
        MILP += pulp.lpSum(delta_ec_TF[(i, gamma)] for i in range(1, t_T + 1)) == ec_TF[gamma], f"milp-(85)-7-{gamma}"

    # -------- Constraint (86) --------
    for (mu1, mu2, m) in Gamma_int_less:
        MILP += ec_C[(mu1, mu2, m)] + ec_C[(mu2, mu1, m)] + \
                ec_T[(mu1, mu2, m)] + ec_T[(mu2, mu1, m)] + \
                ec_F[(mu1, mu2, m)] + ec_F[(mu2, mu1, m)] + \
                ec_CT[(mu1, mu2, m)] + ec_CT[(mu2, mu1, m)] + \
                ec_TC[(mu1, mu2, m)] + ec_TC[(mu2, mu1, m)] + \
                ec_CF[(mu1, mu2, m)] + ec_CF[(mu2, mu1, m)] + \
                ec_TF[(mu1, mu2, m)] + ec_TF[(mu2, mu1, m)] == ec_int[(mu1, mu2, m)], f"milp-(86)-1-{(mu1, mu2, m)}"
    for gamma in Gamma_int_equal:
        MILP += ec_C[gamma] + ec_T[gamma] + ec_F[gamma] + \
                ec_CT[gamma] + ec_TC[gamma] + \
                ec_CF[gamma] + ec_TF[gamma] == ec_int[gamma], f"milp-(86)-2-{gamma}"

    return MILP

# -------- fv Descriptor for the Number of Fringe-configurations --------
def prepare_variables_fringe_configuration(
    t_C,
    t_T,
    t_F,
    set_F,
    Code_F
):
    #Define fc
    fc = {Code_F[psi]: pulp.LpVariable(f"fc({Code_F[psi]})", 0, t_C + t_T + t_F, cat=pulp.LpInteger)
        for psi in set_F}

    return fc

def add_constraints_fringe_configuration(
    # Model
    MILP,
    # Constants
    t_C,
    t_T,
    t_F,
    set_F,
    Code_F,
    # Binary Variables
    delta_fr_C,
    delta_fr_T,
    delta_fr_F,
    # Integer Variables
    fc
):
    # -------- Constraint (87) --------
    for psi in set_F:
        MILP += pulp.lpSum(delta_fr_C[(i, Code_F[psi])] for i in range(1, t_C + 1)) + \
                pulp.lpSum(delta_fr_T[(i, Code_F[psi])] for i in range(1, t_T + 1)) + \
                pulp.lpSum(delta_fr_F[(i, Code_F[psi])] for i in range(1, t_F + 1)) == fc[Code_F[psi]], f"milp-(87)-{Code_F[psi]}"

    return MILP

# -------- A.13 Constraints for the Graph Convolution --------
def prepare():

    return alpha_CT_C, alpha_CT_T, alpha_TC_C, alpha_TC_T, alpha_CF_C, alpha_CF_F, alpha_TF_F, alpha_TF_T

def add_constraints_graph_convolution(
    # Model
    MILP,
    # Constants
    t_T,
    t_F,
    t_C_tilde,
    k_C,
    I_ge_one,
    I_zero_one,
    I_equal_one,
    Lambda_int,
    MAX_CODE_int,
    head_C,
    tail_C,
    # Binary Variables
    e_C,
    e_T,
    chi_T,
    chi_F,
    delta_chi_T,
    delta_chi_F,
    # Integer Variables
    alpha_C,
    alpha_T,
    alpha_F,
    alpha_CT,
    alpha_TC,
    alpha_CF,
    alpha_TF,
):
    # -------- Constraint (88) --------
    for j in (I_ge_one + I_zero_one + I_equal_one):
        MILP += alpha_C[(j, 1)] >= alpha_C[(tail_C[j], 0)] - MAX_CODE_int * (1 - e_C[j]), f"milp-(88)-1-{j}-1"
        MILP += alpha_C[(j, 1)] <= alpha_C[(tail_C[j], 0)] + MAX_CODE_int * (1 - e_C[j]), f"milp-(88)-1-{j}-2"
        MILP += alpha_C[(j, 2)] >= alpha_C[(head_C[j], 0)] - MAX_CODE_int * (1 - e_C[j]), f"milp-(88)-2-{j}-1"
        MILP += alpha_C[(j, 2)] <= alpha_C[(head_C[j], 0)] + MAX_CODE_int * (1 - e_C[j]), f"milp-(88)-2-{j}-2"
        MILP += alpha_C[(j, 1)] + alpha_C[(j, 2)] <= 2 * MAX_CODE_int * e_C[j], f"milp-(88)-3-{j}"

    # -------- Constraint (89) --------
    for k in range(1, k_C + 1):
        MILP += alpha_CT_T[k] >= alpha_CT[k] - MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-(89)-1-{k}-1"
        MILP += alpha_CT_T[k] <= alpha_CT[k] + MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-(89)-1-{k}-2"
        MILP += alpha_CT_T[k] <= MAX_CODE_int * delta_chi_T[k], f"milp-(89)-2-{k}"

    # -------- Constraint (90) --------
    for k in range(1, k_C + 1):
        MILP += alpha_TC_T[k] >= alpha_TC[k] - MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-(90)-1-{k}-1"
        MILP += alpha_TC_T[k] <= alpha_TC[k] + MAX_CODE_int * (1 - delta_chi_T[k]), f"milp-(90)-1-{k}-2"
        MILP += alpha_TC_T[k] <= MAX_CODE_int * delta_chi_T[k], f"milp-(90)-2-{k}"

    # -------- Constraint (91) --------
    for i in range(1, t_T + 1):
        for k in range(1, k_C + 1):
            MILP += alpha_CT_C[i] >= \
                alpha_C[(tail_C[k], 0)] - MAX_CODE_int * (e_T[i] - chi_T[(i, k)] + 1), f"milp-(91)-1-{i}-{k}-1"
            MILP += alpha_CT_C[i] <= \
                alpha_C[(tail_C[k], 0)] + MAX_CODE_int * (e_T[i] - chi_T[(i, k)] + 1), f"milp-(91)-1-{i}-{k}-2"
            MILP += alpha_TC_C[i] >= \
                alpha_C[(head_C[k], 0)] - MAX_CODE_int * (e_T[i + 1] - chi_T[(i, k)] + 1), f"milp-(91)-2-{i}-{k}-1"
            MILP += alpha_TC_C[i] <= \
                alpha_C[(head_C[k], 0)] + MAX_CODE_int * (e_T[i + 1] - chi_T[(i, k)] + 1), f"milp-(91)-2-{i}-{k}-2"

        MILP += alpha_CT_C[i] + alpha_TC_C[i] <= 2 * MAX_CODE_int * pulp.lpSum(chi_T[(i, k)] 
            for k in range(1, k_C + 1)), f"milp-91-3-{i}"

    # -------- Constraint (92) --------
    for c in range(1, t_C_tilde + 1):
        MILP += alpha_CF_F[c] >= alpha_CF[c] - MAX_CODE_int * (chi_F[c] + 1), f"milp-(92)-1-{c}-1"
        MILP += alpha_CF_F[c] <= alpha_CF[c] + MAX_CODE_int * (chi_F[c] + 1), f"milp-(92)-1-{c}-2"
        MILP += alpha_CF_F[c] <= MAX_CODE_int + delta_chi_F[c], f"milp-(92)-2-{c}"

    # -------- Constraint (93) --------
    for i in range(1, t_F + 1):
        for c in range(1, t_C_tilde + 1):
            MILP += alpha_CF_C[i] >= alpha_C[(c, 0)] - MAX_CODE_int * (1 - chi_F[(i, c)]), f"milp-(93)-1-{i}-{c}-1"
            MILP += alpha_CF_C[i] <= alpha_C[(c, 0)] + MAX_CODE_int * (1 - chi_F[(i, c)]), f"milp-(93)-1-{i}-{c}-2"

        MILP += alpha_CF_C[i] <= MAX_CODE_int * pulp.lpSum(chi_F[(i, c)] for c in range(1, t_C_tilde + 1)), f"milp-93-2-{i}"

    # -------- Constraint (94) --------
    for i in range(1, t_T + 1):
        MILP += alpha_TF_F[i] >= alpha_TF[i] - MAX_CODE_int * (1 - delta_chi_F[t_C_tilde + i]), f"milp-94-1-{i}-1"
        MILP += alpha_TF_F[i] <= alpha_TF[i] + MAX_CODE_int * (1 - delta_chi_F[t_C_tilde + i]), f"milp-94-1-{i}-2"
        MILP += alpha_TF_F[i] <= MAX_CODE_int * delta_chi_F[t_C_tilde + i], f"milp-94-2-{i}"

    # -------- Constraint (95) --------
    for i in range(1, t_F + 1):
        for j in range(1, t_T + 1):
            MILP += alpha_TF_T[j] >= alpha_T[(j, 0)] - MAX_CODE_int * (1 - chi_F[(t_C_tilde + i, j)]), f"milp-95-1-{i}-{j}-1"
            MILP += alpha_TF_T[j] <= alpha_T[(j, 0)] + MAX_CODE_int * (1 - chi_F[(t_C_tilde + i, j)]), f"milp-95-1-{i}-{j}-2"

        MILP += alpha_TF_T[i] <= MAX_CODE_int * pulp.lpSum(delta_chi_F[(t_C_tilde + i, j)] 
            for j in range(1, t_T + 1)), f"milp-95-2-{i}"

    return MILP

def print_gstar_file(
    n_G,
    chi_T,

    t_C,
    t_T,
    index_C,
    index_T,
    graph_adj,
    E_C,
    ch_LB,
    ch_UB,
    I_ge_one,
    I_ge_two,
    I_zero_one,
    I_equal_one,

    set_F_v,
    set_F_E,

    outputfilename
):
    n = round(n_G.value())

    base_vertices = list()
    base_edges = list()
    base_edges_color = list()
    visited = set()
    core_set = set()
    for i in range(1, t_C + 1):
        base_vertices.append(index_C[(i, 0)])
        # print(index_C[(i, 0)])
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
                    # print(tmp)
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
                        # print(c)

    # print(base_edges)
    # print(base_edges_color)

    with open(outputfilename, "w") as f:
        f.write(f"{t_C}\n")
        # print(t_C)
        for i in range(t_C):
            f.write(f"{base_vertices[i]}\n")
            # print(base_vertices[i])
            f.write(f"{ch_LB[i + 1]} {ch_UB[i + 1]}\n")
            # print(f"{ch_LB[i + 1]} {ch_UB[i + 1]}\n")
            ###########
            # index of possible fringe tree
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
            ###########
            # index of possible fringe tree
            for j in set_F_E:
                f.write(f"{j.index} ")
            f.write("\n")


    return
    
def print_sdf_file(
    t_C,
    t_T,
    t_F,
    t_C_tilde,
    n_C,
    n_T,
    n_F,
    c_F,
    Lambda,
    Lambda_int,
    Lambda_ex,
    head_C,
    tail_C,
    I_equal_one,
    I_zero_one,
    I_ge_one,
    I_ge_two,
    set_F,
    Code_F,

    # Variables
    n_G,
    v_T,
    v_F,
    alpha_C,
    alpha_T,
    alpha_F,
    beta_C,
    beta_T,
    beta_F,
    beta_CT,
    beta_TC,
    beta_CF,
    beta_TF,
    chi_T,
    chi_F,
    e_T,
    e_F,
    delta_fr_C,
    delta_fr_T,
    delta_fr_F,

    outputfilename
    ):

    ''' A function to output sdf file.  '''
    n = round(n_G.value())

    graph_col = {i : " " for i in range(1, n + 1)}
    graph_adj = {(i, j) : 0 for i in range(1, n + 1) for j in range(1, n + 1)}

    graph_index = 0
    m = 0

    index_T = {(i, j) : 0 for i in range(1, t_T + 1) for j in range(n_T + 1)}
    index_C = {(i, j) : 0 for i in range(1, t_C + 1) for j in range(n_C + 1)}
    index_F = {(i, j) : 0 for i in range(1, t_F + 1) for j in range(n_F + 1)}


    for i in range(1, t_T + 1):
        if round(v_T[(i, 0)].value()) != 0:
            graph_index += 1
            index_T[(i, 0)] = graph_index
            graph_col[graph_index] = Lambda_int[round(alpha_T[(i, 0)].value()) - 1]

            psi = chemicalRootedTree()
            for psi_tmp in set_F:
                if round(delta_fr_T[i, Code_F[psi_tmp]].value()) != 0:
                    psi = psi_tmp

            for j in range(1, len(psi.vertex)):
                a_tmp, h_tmp = psi.vertex[j]
                graph_index += 1
                index_T[(i, j)] = graph_index
                graph_col[graph_index] = a_tmp

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
        graph_col[graph_index] = Lambda_int[round(alpha_C[(i, 0)].value()) - 1]

        psi = chemicalRootedTree()
        for psi_tmp in set_F:
            if round(delta_fr_C[i, Code_F[psi_tmp]].value()) != 0:
                psi = psi_tmp

        for j in range(1, len(psi.vertex)):
            a_tmp, h_tmp = psi.vertex[j]
            graph_index += 1
            index_C[(i, j)] = graph_index
            graph_col[graph_index] = a_tmp

        for j1 in range(len(psi.vertex)):
            for j2 in range(j1 + 1, len(psi.vertex)):
                if psi.beta[j1][j2] != 0:
                    m += 1
                    ind1 = index_C[(i, j1)]
                    ind2 = index_C[(i, j2)]
                    graph_adj[(ind1, ind2)] = psi.beta[j1][j2]
                    graph_adj[(ind2, ind1)] = psi.beta[j1][j2]

    for i in range(1, t_F + 1):
        if round(v_F[(i, 0)].value()) != 0:
            graph_index += 1
            index_F[(i, 0)] = graph_index
            graph_col[graph_index] = Lambda_int[round(alpha_F[(i, 0)].value()) - 1]

            psi = chemicalRootedTree()
            for psi_tmp in set_F:
                if round(delta_fr_F[i, Code_F[psi_tmp]].value()) != 0:
                    psi = psi_tmp

            for j in range(1, len(psi.vertex)):
                a_tmp, h_tmp = psi.vertex[j]
                graph_index += 1
                index_F[(i, j)] = graph_index
                graph_col[graph_index] = a_tmp

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
            # print("beta_T", ind1, ind2)

    for i in range(2, t_F + 1):
        mul = round(beta_F[i].value())
        if mul != 0:
            m += 1
            ind1 = index_F[(i - 1, 0)]
            ind2 = index_F[(i, 0)]
            graph_adj[(ind1, ind2)] = mul
            graph_adj[(ind2, ind1)] = mul
            # print("beta_F", ind1, ind2)

    for i in (I_equal_one + I_zero_one + I_ge_one):
        mul = round(beta_C[i].value())
        if mul != 0:
            m += 1
            ind1 = index_C[(head_C[i], 0)]
            ind2 = index_C[(tail_C[i], 0)]
            graph_adj[(ind1, ind2)] = mul
            graph_adj[(ind2, ind1)] = mul
            # print("beta_C", ind1, ind2)

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
                    # print("beta_CT", ind1, ind2)
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
                    # print("beta_TC", ind1, ind2)
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
                    # print("beta_CF", ind1, ind2)
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
                    # print("beta_TF", ind1, ind2)
                    break

    # for i in range(1, t_T + 1):
    #     for j in range(1, n_T + 1):
    #         mul = round(beta_T[(i, j)].value())
    #         if mul != 0:
    #             m += 1
    #             ind1 = index_T[(i, prt_T[j])]
    #             ind2 = index_T[(i, j)]
    #             graph_adj[(ind1, ind2)] = mul
    #             graph_adj[(ind2, ind1)] = mul
    #             # print("beta_T_FT", ind1, ind2)
    # for i in range(1, t_C + 1):
    #     for j in range(1, n_C[delta_i[i]] + 1):
    #         mul = round(beta_C[(i, j)].value())
    #         if mul != 0:
    #             m += 1
    #             ind1 = index_C[(i, prt_C[delta_i[i]][j])]
    #             ind2 = index_C[(i, j)]
    #             graph_adj[(ind1, ind2)] = mul
    #             graph_adj[(ind2, ind1)] = mul
    #             # print("beta_C_FT", ind1, ind2)
    # for i in range(1, t_F + 1):
    #     for j in range(1, n_F + 1):
    #         mul = round(beta_F[(i, j)].value())
    #         if mul != 0:
    #             m += 1
    #             ind1 = index_F[(i, prt_F[j])]
    #             ind2 = index_F[(i, j)]
    #             graph_adj[(ind1, ind2)] = mul
    #             graph_adj[(ind2, ind1)] = mul
    #             # print("beta_F_FT", ind1, ind2)

    with open(outputfilename, "w") as f:
        f.write("1\n")
        f.write("MILP_2L\n")
        f.write("fc_vector\n")
        f.write("{:3}{:3}  0  0  0  0  0  0  0  0999 V2000 \n".format(n, m))
        
        for i in range(1, n + 1):
            f.write("    0.0000    0.0000    0.0000 {:2}  0  0  0  0  0  0  0  0  0  0  0  0\n".format(graph_col[i]))

        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                if graph_adj[(i, j)] > 0:
                    f.write("{:3}{:3}{:3}  0  0  0  0\n".format(i, j, graph_adj[(i, j)]))

        f.write("M  END\n")
        f.write("$$$$\n")
        f.close()

    return index_C, index_T, graph_adj
    

