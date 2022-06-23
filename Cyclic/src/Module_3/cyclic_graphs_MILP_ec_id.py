import pulp
import time
from collections import namedtuple
import sys
import pandas as pd

CG_element = namedtuple("CG_element", "symbol, valence, mass")
MAX_BOND = 3
MAX_VAL = 4

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
    h_T : h^{\rm T}
    h_C : h^{\rm C}
    h_F : h^{\rm F}
    sigma : \sigma
# -------- A.4 Descriptor for the Number of Specified Degree --------
    deg_C : {\rm deg}^{\rm C}
    deg_T : {\rm deg}^{\rm T}
    deg_F : {\rm deg}^{\rm F}
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
    alpha_CT : \alpha^{rm CT}
    alpha_TC : \alpha^{rm TC}
    alpha_CF : \alpha^{rm CF}
    alpha_TF : \alpha^{rm TF}
    Delta_ac_C_plus : \Delta_{\rm ac}^{\rm C +}
    Delta_ac_C_minus : \Delta_{\rm ac}^{\rm C -}
    Delta_ac_T_plus : \Delta_{\rm ac}^{\rm T +}
    Delta_ac_T_minus : \Delta_{\rm ac}^{\rm T -}
    Delta_ac_F_plus : \Delta_{\rm ac}^{\rm F +}
    Delta_ac_F_minus : \Delta_{\rm ac}^{\rm F -}
    Delta_ac_CT_plus : \Delta_{\rm ac}^{\rm CT +}
    Delta_ac_CT_minus : \Delta_{\rm ac}^{\rm CT -}
    Delta_ac_TC_plus : \Delta_{\rm ac}^{\rm TC +}
    Delta_ac_TC_minus : \Delta_{\rm ac}^{\rm TC -}
    Delta_ac_CF_plus : \Delta_{\rm ac}^{\rm CF +}
    Delta_ac_CF_minus : \Delta_{\rm ac}^{\rm CF -}
    Delta_ac_TF_plus : \Delta_{\rm ac}^{\rm TF +}
    Delta_ac_TF_minus : \Delta_{\rm ac}^{\rm TF -}
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
    deg_T_CT : {\rm deg}_{\rm T}^{\rm CT}
    deg_T_TC : {\rm deg}_{\rm T}^{\rm TC}
    deg_F_CF : {\rm deg}_{\rm F}^{\rm CF}
    deg_F_TF : {\rm deg}_{\rm F}^{\rm TF}
    Delta_ec_C_plus : \Delta_{\rm ec}^{\rm C +}
    Delta_ec_C_minus : \Delta_{\rm ec}^{\rm C -}
    Delta_ec_T_plus : \Delta_{\rm ec}^{\rm T +}
    Delta_ec_T_minus : \Delta_{\rm ec}^{\rm T -}
    Delta_ec_F_plus : \Delta_{\rm ec}^{\rm F +}
    Delta_ec_F_minus : \Delta_{\rm ec}^{\rm F -}
    Delta_ec_CT_plus : \Delta_{\rm ec}^{\rm CT +}
    Delta_ec_CT_minus : \Delta_{\rm ec}^{\rm CT -}
    Delta_ec_TC_plus : \Delta_{\rm ec}^{\rm TC +}
    Delta_ec_TC_minus : \Delta_{\rm ec}^{\rm TC -}
    Delta_ec_CF_plus : \Delta_{\rm ec}^{\rm CF +}
    Delta_ec_CF_minus : \Delta_{\rm ec}^{\rm CF -}
    Delta_ec_TF_plus : \Delta_{\rm ec}^{\rm TF +}
    Delta_ec_TF_minus : \Delta_{\rm ec}^{\rm TF -}
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
    Lambda_co,
    Lambda_nc,
    Gamma_co,
    Gamma_in,
    Gamma_ex,
    Gamma_co_ac,
    Gamma_in_ac,
    Gamma_ex_ac
):
    #   A function to prepare the set Lambda_dg and Gamma
    set_Lambda_dg = list()
    for atom in set_Lambda:
        for i in range(1, atom.valence + 1):
            set_Lambda_dg.append((atom, i))

    # Lambda = [ele.symbol for ele in set_Lambda]
    # Extra_Lambda = Lambda[:]
    # Extra_Lambda.append("e")
    Lambda_dg = [(ele.symbol, i) for (ele, i) in set_Lambda_dg]  # The set Lambda_dg
    # Extra_Lambda_dg = Lambda_dg[:]
    # Extra_Lambda_dg.append(("e", 0))   # The set Lambda_dg with null element "epsilon"
    epsilon = "e"
    epsilon_dg = ("e", 0)
    Code_Lambda_dg = {ele_dg : i + 1 for i, ele_dg in enumerate(Lambda_dg)}
    Code_Lambda_dg[epsilon_dg] = 0
    Code_Lambda = {ele : i + 1 for i, ele in enumerate(Lambda)}
    Code_Lambda[epsilon] = 0

    Code_Lambda_co = {ele : i + 1 for i, ele in enumerate(Lambda_co)}
    Code_Lambda_co[epsilon] = 0

    Code_Lambda_nc = {ele : i + 1 for i, ele in enumerate(Lambda_nc)}
    Code_Lambda_nc[epsilon] = 0

    MAX_CODE = len(Lambda)
    MAX_CODE_dg = len(Lambda_dg)
    MAX_CODE_co = len(Lambda_co)
    MAX_CODE_nc = len(Lambda_nc)

    Code_Gamma_ec_co = {gamma : i + 1 for i, gamma in enumerate(Gamma_co)}
    Code_Gamma_ec_in = {gamma : i + 1 for i, gamma in enumerate(Gamma_in)}
    Code_Gamma_ec_ex = {gamma : i + 1 for i, gamma in enumerate(Gamma_ex)}
    Code_Gamma_ac_co = {gamma : i + 1 for i, gamma in enumerate(Gamma_co_ac)}
    Code_Gamma_ac_in = {gamma : i + 1 for i, gamma in enumerate(Gamma_in_ac)}
    Code_Gamma_ac_ex = {gamma : i + 1 for i, gamma in enumerate(Gamma_ex_ac)}

    val = {ele.symbol: ele.valence for ele in set_Lambda}
    val_temp = {(ele.symbol, i): ele.valence for (ele, i) in set_Lambda_dg}
    val.update(val_temp)
    mass = {ele.symbol: ele.mass for ele in set_Lambda}

    return Code_Lambda, MAX_CODE, \
        set_Lambda_dg, Lambda_dg, Code_Lambda_dg, MAX_CODE_dg, \
        val, mass, epsilon, epsilon_dg, \
        Code_Lambda_co, Code_Lambda_nc, MAX_CODE_co, MAX_CODE_nc, \
        Code_Gamma_ec_co, Code_Gamma_ec_in, Code_Gamma_ec_ex, Code_Gamma_ac_co, Code_Gamma_ac_in, Code_Gamma_ac_ex

# A function to prepare sets used in A.8 and A.10
def prepare_Gamma_ac(
    Gamma_co,
    Gamma_in,
    Gamma_ex,
    Gamma_co_ac,
    Gamma_in_ac,
    Gamma_ex_ac,
    Code_Lambda,
    Code_Lambda_dg
):
    # Here Gamma_co = Gamma_co_less \cup Gamma_co_equal
    Gamma_co_less = list()
    Gamma_co_equal = list()
    Gamma_co_great = list()

    for ((a1, d1), (a2, d2), m) in Gamma_co:
        if Code_Lambda_dg[(a1, d1)] < Code_Lambda_dg[(a2, d2)]:
            Gamma_co_less.append(((a1, d1), (a2, d2), m))
            # Gamma_co_great.append(((a2, d2), (a1, d1), m))
        elif Code_Lambda_dg[(a1, d1)] == Code_Lambda_dg[(a2, d2)]:
            Gamma_co_equal.append(((a1, d1), (a2, d2), m))
        elif Code_Lambda_dg[(a1, d1)] > Code_Lambda_dg[(a2, d2)]:
            # Gamma_co_less.append(((a2, d2), (a1, d1), m))
            Gamma_co_great.append(((a1, d1), (a2, d2), m))

    # Gamma_co_ac = list()
    # Gamma_in_ac = list()
    # Gamma_ex_ac = list()
    Gamma_co_ac_less = list()
    Gamma_co_ac_equal = list()
    Gamma_co_ac_great = list()

    # for ((a1, d1), (a2, d2), m) in Gamma_co:
    #     if (a1, a2, m) not in Gamma_co_ac:
    #         Gamma_co_ac.append((a1, a2, m))
    #     if (a2, a1, m) not in Gamma_co_ac:
    #         Gamma_co_ac.append((a2, a1, m))
    # for ((a1, d1), (a2, d2), m) in Gamma_in:
    #     if (a1, a2, m) not in Gamma_in_ac:
    #         Gamma_in_ac.append((a1, a2, m))
    # for ((a1, d1), (a2, d2), m) in Gamma_ex:
    #     if (a1, a2, m) not in Gamma_ex_ac:
    #         Gamma_ex_ac.append((a1, a2, m))  
    # for ((a1, d1), (a2, d2), m) in Gamma_co_less:
    #     if Code_Lambda[a1] < Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_less:
    #             Gamma_co_ac_less.append((a1, a2, m))
    #     elif Code_Lambda[a1] ==  Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_equal:
    #             Gamma_co_ac_equal.append((a1, a2, m))
    #     elif Code_Lambda[a1] > Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_great:
    #             Gamma_co_ac_great.append((a1, a2, m)) 
    # for ((a1, d1), (a2, d2), m) in Gamma_co_equal:
    #     if Code_Lambda[a1] < Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_less:
    #             Gamma_co_ac_less.append((a1, a2, m))
    #     elif Code_Lambda[a1] ==  Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_equal:
    #             Gamma_co_ac_equal.append((a1, a2, m))
    #     elif Code_Lambda[a1] > Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_great:
    #             Gamma_co_ac_great.append((a1, a2, m)) 
    # for ((a1, d1), (a2, d2), m) in Gamma_co_great:
    #     if Code_Lambda[a1] < Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_less:
    #             Gamma_co_ac_less.append((a1, a2, m))
    #     elif Code_Lambda[a1] ==  Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_equal:
    #             Gamma_co_ac_equal.append((a1, a2, m))
    #     elif Code_Lambda[a1] > Code_Lambda[a2]:
    #         if (a1, a2, m) not in Gamma_co_ac_great:
    #             Gamma_co_ac_great.append((a1, a2, m)) 
    for ((a1, d1), (a2, d2), m) in Gamma_co:
        if Code_Lambda[a1] < Code_Lambda[a2]:
            if (a1, a2, m) not in Gamma_co_ac_less:
                Gamma_co_ac_less.append((a1, a2, m))
            if (a2, a1, m) not in Gamma_co_ac_great:
                Gamma_co_ac_great.append((a2, a1, m))
        elif Code_Lambda[a1] ==  Code_Lambda[a2]:
            if (a1, a2, m) not in Gamma_co_ac_equal:
                Gamma_co_ac_equal.append((a1, a2, m))
        elif Code_Lambda[a1] > Code_Lambda[a2]:
            if (a1, a2, m) not in Gamma_co_ac_great:
                Gamma_co_ac_great.append((a1, a2, m))  
            if (a2, a1, m) not in Gamma_co_ac_less:
                Gamma_co_ac_less.append((a2, a1, m))

    Gamma_tilde_ac_C = Gamma_co_ac_equal + Gamma_co_ac_great + Gamma_co_ac_less
    Gamma_tilde_ac_T = Gamma_co_ac_equal + Gamma_co_ac_great + Gamma_co_ac_less
    Gamma_tilde_ac_CT = Gamma_co_ac_equal + Gamma_co_ac_great + Gamma_co_ac_less
    Gamma_tilde_ac_TC = Gamma_co_ac_equal + Gamma_co_ac_great + Gamma_co_ac_less
    Gamma_tilde_ac_F = Gamma_in_ac
    Gamma_tilde_ac_CF = Gamma_in_ac
    Gamma_tilde_ac_TF = Gamma_in_ac
    Gamma_tilde_ac_C_ex = Gamma_ex_ac
    Gamma_tilde_ac_T_ex = Gamma_ex_ac
    Gamma_tilde_ac_F_ex = Gamma_ex_ac

    Gamma_tilde_ec_C = Gamma_co_equal + Gamma_co_less + Gamma_co_great
    Gamma_tilde_ec_T = Gamma_co_equal + Gamma_co_less + Gamma_co_great
    Gamma_tilde_ec_CT = Gamma_co_equal + Gamma_co_less + Gamma_co_great
    Gamma_tilde_ec_TC = Gamma_co_equal + Gamma_co_less + Gamma_co_great
    Gamma_tilde_ec_F = Gamma_in
    Gamma_tilde_ec_CF = Gamma_in
    Gamma_tilde_ec_TF = Gamma_in
    Gamma_tilde_ec_C_ex = Gamma_ex
    Gamma_tilde_ec_T_ex = Gamma_ex
    Gamma_tilde_ec_F_ex = Gamma_ex

    return Gamma_co_less, Gamma_co_equal, Gamma_co_great, \
        Gamma_co_ac_less, Gamma_co_ac_equal, Gamma_co_ac_great, \
        Gamma_tilde_ac_C, Gamma_tilde_ac_T, Gamma_tilde_ac_CT, Gamma_tilde_ac_TC, \
        Gamma_tilde_ac_F, Gamma_tilde_ac_CF, Gamma_tilde_ac_TF, Gamma_tilde_ac_C_ex, Gamma_tilde_ac_T_ex, Gamma_tilde_ac_F_ex, \
        Gamma_tilde_ec_C, Gamma_tilde_ec_T, Gamma_tilde_ec_CT, Gamma_tilde_ec_TC, \
        Gamma_tilde_ec_F, Gamma_tilde_ec_CF, Gamma_tilde_ec_TF, Gamma_tilde_ec_C_ex, Gamma_tilde_ec_T_ex, Gamma_tilde_ec_F_ex

# A function to prepare sets used in A.9
def prepare_Lambda_tilde(
    Lambda_dg_co,
    Lambda_dg_nc
):

    Lambda_tilde_dg_C = Lambda_dg_co[:]
    Lambda_tilde_dg_T = Lambda_dg_co[:]
    Lambda_tilde_dg_F = Lambda_dg_nc[:]
    Lambda_tilde_dg_C_nc = Lambda_dg_nc[:]
    Lambda_tilde_dg_T_nc = Lambda_dg_nc[:]
    Lambda_tilde_dg_F_nc = Lambda_dg_nc[:]

    return Lambda_tilde_dg_C, Lambda_tilde_dg_T, Lambda_tilde_dg_F, Lambda_tilde_dg_C_nc, Lambda_tilde_dg_T_nc, Lambda_tilde_dg_F_nc

# Prepare Constants of Seed Graph
def prepare_dataset_for_scheme_graph(
    n_star,
    rho,
    V_C,
    E_C,
    E_ge_two,
    E_ge_one,
    cs_LB,
    cs_UB,
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
    Lambda_star,
    d_max,
    r_Gc
):
    t_C = len(V_C)
    t_C_tilde = sum(1 for i in V_C if bl_UB[i] == 1)
    t_T = cs_UB - len(V_C)
    ch_LB_tilde = sum(ch_LB[v] for v in V_C) + sum(ch_LB[e] for e in (E_ge_two + E_ge_one))
    bl_LB_star = sum(bl_LB[v] for v in V_C) + sum(bl_LB[e] for e in (E_ge_two + E_ge_one))
    bl_UB_star = sum(bl_UB[v] for v in V_C) + sum(bl_UB[e] for e in (E_ge_two + E_ge_one))
    ell_inl_star = sum(max(ch_UB[v] - rho, 0) for v in V_C) + \
                   sum(bl_UB[e] * max(ch_UB[e] - rho, 0) for e in (E_ge_one + E_ge_two))
    #Using cs_LB instead of cs_LB_star because cs_LB_star is not defined in paper
    t_F = min(n_star - cs_LB - max(ch_LB_tilde, rho * bl_LB_star), ell_inl_star)
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

    if d_max == 3:
        n_T = 5            # T(2, 2, 2)
        n_C = {0: 0, 2: 5, 1: 3} # 2: T(2, 2, 2), 1: T(1, 2, 2)
        n_F = 5            # T(2, 2, 2)
    elif d_max == 4:
        n_T = 6            # T(2, 3, 2)
        n_C = {0: 0, 2: 6, 1: 3} # 2: T(2, 3, 2), 1: T(1, 3, 2)
        n_F = 9            # T(3, 3, 2)

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

    Cld_T = [set() for _ in range(n_T + 1)]
    Cld_C = {i: [set() for _ in range(n_C[i] + 1)] for i in range(3)}
    Cld_F = [set() for _ in range(n_F + 1)]

    prt_T = list()
    prt_C = {i: list() for i in range(3)}
    prt_F = list()

    P_prc_T = list()
    P_prc_C = {i: list() for i in range(3)}
    P_prc_F = list()

    Dsn_T = [set() for _ in range(rho + 1)]
    Dsn_C = {i: [set() for _ in range(rho + 1)] for i in range(3)}
    Dsn_F = [set() for _ in range(rho + 1)]

    j_T = [0 for _ in range(rho + 1)]
    j_C = {i: [0 for _ in range(rho + 1)] for i in range(3)}
    j_F = [0 for _ in range(rho + 1)]

    tie_breaker_T = list()
    tie_breaker_C = {i: list() for i in range(3)}
    tie_breaker_F = list()

    if d_max == 3:
        Cld_T[0] = {1, 2}
        Cld_T[1] = {3, 4}
        Cld_T[2] = {5}
        Cld_C[0][0] = {}
        Cld_C[1][0] = {1}
        Cld_C[1][1] = {2, 3}
        Cld_C[2][0] = {1, 2}
        Cld_C[2][1] = {3, 4}
        Cld_C[2][2] = {5}
        Cld_F[0] = {1, 2}
        Cld_F[1] = {3, 4}
        Cld_F[2] = {5}

        prt_T = [0, 0, 0, 1, 1, 2]
        prt_C[0] = [0]
        prt_C[1] = [0, 0, 1, 1]
        prt_C[2] = [0, 0, 0, 1, 1, 2]
        prt_F = [0, 0, 0, 1, 1, 2]

        P_prc_T = [(0, 1), (1, 2), (1, 3), (2, 5), (3, 4), (3, 5)]
        P_prc_C[0] = []
        P_prc_C[1] = [(0, 1), (1, 2), (2, 3)]
        P_prc_C[2] = [(0, 1), (1, 2), (1, 3), (2, 5), (3, 4), (3, 5)]
        P_prc_F = [(0, 1), (1, 2), (1, 3), (2, 5), (3, 4), (3, 5)]
        
        Dsn_T[0] = {0}
        Dsn_T[1] = {1, 2}
        Dsn_T[2] = {3, 4, 5}
        Dsn_C[0][0] = {0}
        Dsn_C[0][1] = {}
        Dsn_C[0][2] = {}
        Dsn_C[1][0] = {0}
        Dsn_C[1][1] = {1}
        Dsn_C[1][2] = {2, 3}
        Dsn_C[2][0] = {0}
        Dsn_C[2][1] = {1, 2}
        Dsn_C[2][2] = {3, 4, 5}
        Dsn_F[0] = {0}
        Dsn_F[1] = {1, 2}
        Dsn_F[2] = {3, 4, 5}
        j_T = [0, 1, 3]
        j_C[0] = [0]
        j_C[1] = [0, 1, 2]
        j_C[2] = [0, 1, 3]
        j_F = [0, 1, 3]

        tie_breaker_T = {(3, 4)}
        tie_breaker_C[0] = {}
        tie_breaker_C[1] = {(2, 3)}
        tie_breaker_C[2] = {(3, 4)}
        tie_breaker_F = {(3, 4)}

    elif d_max == 4:
        Cld_T[0] = {1, 2}
        Cld_T[1] = {3, 4, 5}
        Cld_T[2] = {6}
        Cld_C[0][0] = {}
        Cld_C[1][0] = {1}
        Cld_C[1][1] = {2, 3}
        Cld_C[2][0] = {1, 2}
        Cld_C[2][1] = {3, 4, 5}
        Cld_C[2][2] = {6}
        Cld_F[0] = {1, 2, 3}
        Cld_F[1] = {4, 5, 6}
        Cld_F[2] = {7, 8}
        Cld_F[3] = {9}

        prt_T = [0, 0, 0, 1, 1, 1, 2]
        prt_C[0] = [0]
        prt_C[1] = [0, 0, 1, 1]
        prt_C[2] = [0, 0, 0, 1, 1, 1, 2]
        prt_F = [0, 0, 0, 0, 1, 1, 1, 2, 2, 3]

        P_prc_T = [(0, 1), (1, 2), (1, 3), (2, 6), (3, 4), (3, 6), (4, 5)]
        P_prc_C[0] = []
        P_prc_C[1] = [(0, 1), (1, 2), (2, 3)]
        P_prc_C[2] = [(0, 1), (1, 2), (1, 3), (2, 6), (3, 4), (3, 6), (4, 5)]
        P_prc_F = [(0, 1), (1, 2), (1, 4), (2, 3), (2, 7), (3, 9), (4, 5), (4, 7), (5, 6), (5, 8), (7, 8)]
        
        Dsn_T[0] = {0}
        Dsn_T[1] = {1, 2}
        Dsn_T[2] = {3, 4, 5, 6}
        Dsn_C[0][0] = {0}
        Dsn_C[0][1] = {}
        Dsn_C[0][2] = {}
        Dsn_C[1][0] = {0}
        Dsn_C[1][1] = {1}
        Dsn_C[1][2] = {2, 3}
        Dsn_C[2][0] = {0}
        Dsn_C[2][1] = {1, 2}
        Dsn_C[2][2] = {3, 4, 5, 6}
        Dsn_F[0] = {0}
        Dsn_F[1] = {1, 2, 3}
        Dsn_F[2] = {4, 5, 6, 7, 8, 9}
        j_T = [0, 1, 3]
        j_C[0] = [0]
        j_C[1] = [0, 1, 2]
        j_C[2] = [0, 1, 3]
        j_F = [0, 1, 4]

        tie_breaker_T = {(3, 4), (4, 5)}
        tie_breaker_C[0] = {}
        tie_breaker_C[1] = {(2, 3)}
        tie_breaker_C[2] = {(3, 4), (4, 5)}
        tie_breaker_F = {(4, 5), (5, 6), (7, 8)}

    # Here r(Gc) should be the rank of the graph, by Zhu
    m_UB_co = min(n_star + r_Gc, m_C + sum(ell_UB[e] for e in (I_ge_one + I_ge_two)))
    m_UB_nc = min(n_star, t_F + sum(n_C[delta_i[i]] for i in range(1, t_C + 1)) + n_T * t_T + n_F * t_F) - 1
    m_UB = min(n_star + r_Gc, m_UB_co + m_UB_nc)
        
    return t_C, t_T, t_F, k_C, k_C_tilde, c_F, m_C, head_C, tail_C, tail_F, E_C_plus, E_C_minus, \
        n_C, n_T, n_F, delta_i, Cld_T, Cld_C, Cld_F, prt_T, prt_C, prt_F, P_prc_T, P_prc_C, P_prc_F, Dsn_T, Dsn_C, Dsn_F, \
        j_T, j_C, j_F, t_C_tilde, I_ge_one_plus, I_ge_one_minus, I_ge_two_plus, I_ge_two_minus, I_zero_one_plus, I_zero_one_minus, \
        I_equal_one_plus, I_equal_one_minus, m_UB_co, m_UB_nc, m_UB, \
        tie_breaker_T, tie_breaker_C, tie_breaker_F

def prepare_fv( 
    ann_training_data_filename,
    Lambda_dg_co,
    Lambda_dg_nc,
    Gamma_co,
    Gamma_in,
    Gamma_ex,
    n_G,
    cs,
    ch_G,
    bl_G,
    MASS,
    bd_co,
    bd_in,
    bd_ex,
    n_H,
    dg_co,
    dg_nc,
    ns_co,
    ns_nc,
    ec_co,
    ec_in,
    ec_ex
):
    ''' A function to prepare some variables about the descriptors 
        from the feature vector '''

    descriptors = dict()      # get variable from fv file
    # stringoutput = dict()     # string used for output
    fv = pd.read_csv(ann_training_data_filename, sep=",")
    num_fv = len(list(fv.columns))

    # some variables used for AD constraint
    ceset = set()
    acset = set()
    desc_index_v = set()
    desc_index_e = set()
    mass_ind = -1

    for i, fv_name in enumerate(list(fv.columns)):
        if fv_name == "CID":
            pass
        elif fv_name == "n":
            descriptors[i] = n_G
            # stringoutput[i] = "               n : "
        elif fv_name == "cs":
            descriptors[i] = cs
        elif fv_name == "ch":
            descriptors[i] = ch_G
        elif fv_name == "bl_2":
            descriptors[i] = bl_G
        elif fv_name == "ms":
            descriptors[i] = MASS# / (10 * n_G)
            mass_ind = i
            # stringoutput[i] = "             M/n : "
        elif fv_name == "bd_co_2":
            descriptors[i] = bd_co[2]
            # stringoutput[i] = "    #double_bond_in : "
            desc_index_e.add(i)
        elif fv_name == "bd_co_3":
            descriptors[i] = bd_co[3]
            # stringoutput[i] = "    #double_bond_in : "
        elif fv_name == "bd_in_2":
            descriptors[i] = bd_in[2]
            # stringoutput[i] = "    #double_bond_in : "
        elif fv_name == "bd_in_3":
            descriptors[i] = bd_in[3]
            # stringoutput[i] = "    #double_bond_in : "
        elif fv_name == "bd_ex_2":
            descriptors[i] = bd_ex[2]
            # stringoutput[i] = "    #double_bond_in : "
        elif fv_name == "bd_ex_3":
            descriptors[i] = bd_ex[3]
            # stringoutput[i] = "    #double_bond_in : "

        elif fv_name == "nsH":
            descriptors[i] = n_H
            # stringoutput[i] = "               H : "
        elif fv_name[0] == "d" and fv_name[4] == "o":
            j = int(fv_name[6])
            descriptors[i] = dg_co[j]
            desc_index_v.add(i)
            # stringoutput[i] = "        #degree{}_in : ".format(j)
        elif fv_name[0] == "d" and fv_name[4] == "c":
            j = int(fv_name[6])
            descriptors[i] = dg_nc[j]
            desc_index_v.add(i)
            # stringoutput[i] = "        #degree{}_ex : ".format(j)
        elif fv_name[0] == "n" and fv_name[4] == "o":
            j = int(fv_name[-1])
            pos = len(fv_name) - 1
            while fv_name[pos] != '_':
                pos -= 1
            ele = fv_name[pos + 1: len(fv_name) - 1]
            if (ele, j) in Lambda_dg_co:
                descriptors[i] = ns_co[(ele, j)]
            else:
                descriptors[i] = 0
            desc_index_v.add(i)
            # stringoutput[i] = "        #degree{}_ex : ".format(j)
        elif fv_name[0] == "n" and fv_name[4] == "c":
            j = int(fv_name[-1])
            pos = len(fv_name) - 1
            while fv_name[pos] != '_':
                pos -= 1
            ele = fv_name[pos + 1 : len(fv_name) - 1]
            if (ele, j) in Lambda_dg_nc:
                descriptors[i] = ns_nc[(ele, j)]
            else:
                descriptors[i] = 0;
            desc_index_v.add(i)
            # stringoutput[i] = "        #degree{}_ex : ".format(j)
        else:
            str_tmp = fv_name
            pos = list()
            posnum = -1
            for k in range(len(str_tmp)):
                if str_tmp[k] == '_' and k >= 3:
                    pos.append(k)
            ele1 = str_tmp[pos[0] + 1: pos[1] - 1]
            d1 = int(str_tmp[pos[1] - 1])
            ele2 = str_tmp[pos[1] + 1: pos[2] - 1]
            d2 = int(str_tmp[pos[2] - 1])
            mul = int(str_tmp[-1])
            ec_tmp = ((ele1, d1), (ele2, d2), mul)
            if str_tmp[pos[0] - 1] == "o":
                # print(str_tmp, ec_tmp)
                if ec_tmp in Gamma_co:
                    descriptors[i] = ec_co[ec_tmp]
                else:
                    descriptors[i] = 0
            elif str_tmp[pos[0] - 1] == "n":
                # print(str_tmp, ec_tmp)
                if ec_tmp in Gamma_in:
                    descriptors[i] = ec_in[ec_tmp]
                else:
                    descriptors[i] = 0
            elif str_tmp[pos[0] - 1] == "x":
                # print(str_tmp, ec_tmp)
                if ec_tmp in Gamma_ex:
                    descriptors[i] = ec_ex[ec_tmp]
                else:
                    descriptors[i] = 0
            desc_index_e.add(i)

    return descriptors, num_fv, desc_index_v, desc_index_e, mass_ind

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
                        descriptors,
                        num_fv,
                        y, mass_ind,
                        mass_n, 
                        prop = "def"):
    ''' A function to add constraints used in ANN '''



    for i in range(1, num_fv):
        if i == mass_ind:  # Mass
            MILP += y[(1, i)] == mass_n,\
            "ann_input_{}_{}".format(i, prop)
        else:
            MILP += y[(1, i)] == descriptors[i], \
            "ann_input_{}_{}".format(i, prop)

    return MILP

# -------- A.1 Selecting Core-vertices and Core-edges   -------- 
def prepare_variables_selecting_core(
    # Constants
    t_C,
    t_T,
    k_C,
    m_C,
    ell_LB,
    ell_UB,
    bl_LB,
    bl_UB,
    cs_LB,
    cs_UB,
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
    clr_T = {c : pulp.LpVariable(f"clr_T({c})", ell_LB[c] - 1, ell_UB[c] - 1, cat = pulp.LpInteger)
            for c in range(1, k_C + 1)}
    clr_T[0] = pulp.LpVariable("clr_T(0)", 0, t_T, cat = pulp.LpInteger)

    # Define delta_chi_T
    delta_chi_T = {c : pulp.LpVariable(f"delta_chi_T({c})", 0, 1, cat = pulp.LpBinary)
            for c in range(0, k_C + 1)} 
    delta_chi_T_tmp = {(i, c) : pulp.LpVariable(f"delta_chi_T({i},{c})", 0, 1, cat = pulp.LpBinary)
            for i in range(1, t_T + 1) for c in range(0, k_C + 1)}
    delta_chi_T.update(delta_chi_T_tmp)

    # Define deg_tilde_C_plus
    deg_tilde_C_plus = {i : pulp.LpVariable(f"deg_tilde_C_plus({i})", 0, MAX_VAL, cat = pulp.LpInteger)
            for i in range(1, t_C + 1)}

    # Define deg_tilde_C_minus
    deg_tilde_C_minus = {i : pulp.LpVariable(f"deg_tilde_C_minus({i})", 0, MAX_VAL, cat = pulp.LpInteger)
            for i in range(1, t_C + 1)}

    # Define cs
    cs = pulp.LpVariable("cs", cs_LB, cs_UB, cat = pulp.LpInteger)

    return e_C, v_T, e_T, chi_T, clr_T, delta_chi_T, deg_tilde_C_plus, deg_tilde_C_minus, cs

def add_constraints_selecting_core(
    # Model
    MILP,

    # Constants
    t_C,
    t_T,
    k_C,
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
    deg_tilde_C_minus,
    cs
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
        MILP += pulp.lpSum(delta_chi_T[(i, c)]
                for c in range(0, k_C + 1)) == 1, f"milp-(7)-2-{i}"        
        MILP += pulp.lpSum(c * delta_chi_T[(i, c)]
                for c in range(0, k_C + 1)) == chi_T[i], f"milp-(7)-3-{i}"
        MILP += delta_chi_T[(i, 0)] == 1 - v_T[(i, 0)], f"milp-(7)-1-{i}"

    # -------- Constraint (8) --------
    # Constraint of delta_chi_T, clr_T, for c in [0, k_C]
    for c in range(0, k_C + 1):
        MILP += pulp.lpSum(delta_chi_T[(i, c)]
                for i in range(1, t_T + 1)) == clr_T[c], f"milp-(8)-1-{c}"
        MILP += pulp.lpSum(delta_chi_T[(i, c)]
                for i in range(1, t_T + 1)) <= t_T * delta_chi_T[c], f"milp-(8)-2-{c}-1"
        MILP += pulp.lpSum(delta_chi_T[(i, c)]
                for i in range(1, t_T + 1)) >= delta_chi_T[c], f"milp-(8)-2-{c}-2"

    # -------- Constraint (9) --------
    # Constraint of e_T, v_T, chi_T, for i in [2, t_T]
    for i in range(2, t_T + 1):
        MILP += v_T[(i - 1, 0)] >= v_T[(i, 0)], f"milp-(9)-1-{i}"
        MILP += chi_T[i - 1] - chi_T[i] <= k_C * (v_T[(i - 1, 0)] - e_T[i]), f"milp-(9)-2-{i}-1"
        MILP += chi_T[i - 1] - chi_T[i] >= v_T[(i - 1, 0)] - e_T[i], f"milp-(9)-2-{i}-2"

    # -------- Constraint (10) --------
    # Constraint of v_T (for specify sigma_co)
    MILP += t_C + pulp.lpSum(v_T[(i, 0)] for i in range(1, t_T + 1)) == cs, f"milp-(10)"

    return MILP

# -------- A.2 Constraints for Including Internal Vertices and Edges --------
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
    E_C
):
    # Define bl_G
    # Here max function should be sum, by Zhu.
    bl_G_UB = sum(bl_UB[v] for v in V_C) + sum(bl_UB[e] for e in (E_ge_one + E_ge_two))
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
    # modified range of c from paper([1,c_F] -> [0,c_F])
    clr_F = {i : pulp.LpVariable(f"clr_F({i})", 0, t_F, cat = pulp.LpInteger)
                for i in range(0, c_F + 1)}
    # Define delta_chi_F
    # modified range of c from paper([1,c_F] -> [0,c_F])
    delta_chi_F = {c : pulp.LpVariable(f"delta_chi_F({c})", bl_LB[c], 1, cat = pulp.LpInteger)
                for c in range(1, t_C_tilde + 1)}
    delta_chi_F[0] = pulp.LpVariable(f"delta_chi_F(0)", 0, 1, cat=pulp.LpBinary)
    delta_chi_F_tmp = {c: pulp.LpVariable(f"delta_chi_F({c})", 0, 1, cat=pulp.LpBinary)
                for c in range(t_C_tilde + 1, c_F + 1)}
    delta_chi_F.update(delta_chi_F_tmp)
    delta_chi_F_tmp = {(i, c) : pulp.LpVariable(f"delta_chi_F({i},{c})", 0, 1, cat = pulp.LpBinary)
                for i in range(1, t_F + 1) for c in range(0, c_F + 1)}
    delta_chi_F.update(delta_chi_F_tmp)
    # Define sigma
    sigma = {c : pulp.LpVariable(f"sigma({c})", 0, 1, cat = pulp.LpBinary)
                for c in range(0, c_F + 1)}
    # Define sigma_C, sigma_T, sigma_F
    sigma_C = {i : pulp.LpVariable(f"sigma_C({i})", 0, 1, cat = pulp.LpBinary)
                for i in range(1, t_C + 1)}
    sigma_T = {i : pulp.LpVariable(f"sigma_T({i})", 0, 1, cat = pulp.LpBinary)
                for i in range(1, t_T + 1)}
    sigma_F = {i : pulp.LpVariable(f"sigma_F({i})", 0, 1, cat = pulp.LpBinary)
                for i in range(1, t_F + 1)}
    # Define bl
    bl = {(c, i) : pulp.LpVariable(f"bl({c},{i})", 0, 1, cat = pulp.LpBinary)
                for c in (I_ge_two + I_ge_one) for i in range(1, t_T + 1)}
    return bl_G, v_F, e_F, chi_F, clr_F, delta_chi_F, sigma, sigma_C, sigma_T, sigma_F, bl

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
    delta_chi_T,
    e_F,
    v_F,
    v_T,
    sigma,
    sigma_C,
    sigma_T,
    sigma_F,
    # Integer Variables
    chi_F,
    clr_F,
    bl,
    bl_G
)->pulp.LpProblem:
    # -------- Constraint (11) --------
    # Constraint of delta_chi_F, chi_F, for i in [1, t_F]
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_chi_F[(i, c)]
                for c in range(0, c_F + 1)) == 1, f"milp-(11)-2-{i}"
        MILP += pulp.lpSum(c * delta_chi_F[(i, c)]
                for c in range(1, c_F + 1)) == chi_F[i], f"milp-(11)-3-{i}"
        MILP += delta_chi_F[(i, 0)] == 1 - v_F[(i, 0)], f"milp-(11)-1-{i}"

    # -------- Constraint (12) --------
    # Constraint of delta_chi_F, clr_F, for c in [0, c_F]
    for c in range(0, c_F + 1):
        MILP += pulp.lpSum(delta_chi_F[(i, c)]
                for i in range(1, t_F + 1)) == clr_F[c], f"milp-(12)-1-{c}"
        MILP += pulp.lpSum(delta_chi_F[(i, c)]
                for i in range(1, t_F + 1)) <= t_F * delta_chi_F[c], f"milp-(12)-2-{c}-1"
        MILP += pulp.lpSum(delta_chi_F[(i, c)]
                for i in range(1, t_F + 1)) >= delta_chi_F[c], f"milp-(12)-2-{c}-2"

    # -------- Constraint (13) --------
    # Constraint of e_F
    MILP += e_F[1] == 0, f"milp-(13)-1"
    MILP += e_F[t_F + 1] == 0, f"milp-(13)-2"

    # -------- Constraint (14) --------
    # Constraint of e_F, chi_F, v_F, for i in [2, t_F]
    for i in range(2, t_F + 1):
        MILP += v_F[(i - 1, 0)] >= v_F[(i, 0)], f"milp-(14)-1-{i}"
        MILP += chi_F[i - 1] - chi_F[i] <= c_F * (v_F[(i - 1, 0)] - e_F[i]), f"milp-(14)-2-{i}-1"
        MILP += chi_F[i - 1] - chi_F[i] >= v_F[(i - 1, 0)] - e_F[i], f"milp-(14)-2-{i}-2"

    #########################
    # Such constraint is necessary, or a vertex in v_F can attach to a vertex t_T which is not selected.
    # for i in range(1, t_T + 1):
    #     MILP += delta_chi_F[t_C_tilde + i] <= v_T[(i, 0)], f"milp-(test)-{i}"
    #########################  

    # -------- Constraint (15) --------
    # Constraint of delta_chi_F, sigma, sigma_C, sigma_T, sigma_F, bl_G
    MILP += pulp.lpSum(delta_chi_F[c] for c in range(1, c_F + 1)) == bl_G, f"milp-(15)-1"
    MILP += pulp.lpSum(sigma[c] for c in range(1, c_F + 1)) + \
            pulp.lpSum(sigma_C[i] for i in range(1, t_C + 1)) + \
            pulp.lpSum(sigma_T[i] for i in range(1, t_T + 1)) + \
            pulp.lpSum(sigma_F[i] for i in range(1, t_F + 1)) == 1, f"milp-(15)-2"

    # -------- Constraint (16) --------
    for c in (I_ge_two + I_ge_one):
        for i in range(1, t_T + 1):
            MILP += bl[(c, i)] >= delta_chi_F[t_C_tilde + i] + delta_chi_T[(i, c)] - 1, f"milp-(16)-{c}-{i}"

    # -------- Constraint (17) --------
    MILP += pulp.lpSum(bl[(c, i)] for c in (I_ge_one + I_ge_two) for i in range(1, t_T + 1)) <= \
            pulp.lpSum(delta_chi_F[t_C_tilde + i] for i in range(1, t_T + 1)), f"milp-(17)"

    # -------- Constraint (18) --------
    for c in (I_ge_two + I_ge_one):
        MILP += pulp.lpSum(bl[(c, i)] for i in range(1, t_T + 1)) >= bl_LB[E_C[c]], f"milp-(18)-{c}-1"
        MILP += pulp.lpSum(bl[(c, i)] for i in range(1, t_T + 1)) <= bl_UB[E_C[c]], f"milp-(18)-{c}-2"
    
    return MILP

# -------- A.3 Constraints for Including Fringe-trees --------
def prepare_variables_fringe_trees(
    n_LB, 
    n_star,
    rho,
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
    sigma
):
    # Define n_G
    n_G = pulp.LpVariable(f"n_G", n_LB, n_star, cat=pulp.LpInteger)

    # Define ch_G
    ch_G = pulp.LpVariable(f"ch_G", 0, max(max([ch_UB[v] for v in V_C]), 
                max([ch_UB[e] for e in (E_ge_one + E_ge_two)], default=0)), cat=pulp.LpInteger)

    # Define v_T
    v_T_temp = {(i, j): pulp.LpVariable(f"v_T({i},{j})", 0, 1, cat=pulp.LpBinary)
                for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    v_T.update(v_T_temp)

    # Define v_C, V_F
    v_C = {(i, j): pulp.LpVariable(f"v_C({i},{j})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_C + 1) for j in range(n_C[delta_i[i]] + 1)}
    v_F_temp = {(i, j): pulp.LpVariable(f"v_F({i},{j})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    v_F.update(v_F_temp)

    # Define h_T, h_C, h_F
    h_T = {i: pulp.LpVariable(f"h_T({i})", 0, rho, cat=pulp.LpInteger)
            for i in range(1, t_T + 1)}
    h_C = {i: pulp.LpVariable(f"h_C({i})", 0, rho, cat=pulp.LpInteger)
            for i in range(1, t_C + 1)}
    h_F = {i: pulp.LpVariable(f"h_F({i})", 0, rho, cat=pulp.LpInteger)
            for i in range(1, t_F + 1)}

    # Define sigma
    sigma_temp = {(c, i): pulp.LpVariable(f"sigma({c},{i})", 0, 1, cat=pulp.LpBinary)
                    for c in (I_ge_two + I_ge_one) for i in range(1, t_T + 1)} 
    sigma.update(sigma_temp)

    return n_G, ch_G, v_T, v_C, v_F, h_T, h_C, h_F, sigma

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
    j_T,
    j_F,
    j_C,
    n_T,
    n_C,
    n_F,
    Cld_T,
    Cld_C,
    Cld_F,
    P_prc_T,
    P_prc_C,
    P_prc_F,
    delta_i,
    ch_LB,
    ch_UB,
    E_C,
    # Binary Variables
    delta_chi_F,
    delta_chi_T,
    e_F,
    v_F,
    v_C,
    v_T,
    h_T,
    h_C,
    h_F,
    sigma,
    sigma_C,
    sigma_T,
    sigma_F,
    # Integer Variables
    chi_F,
    clr_F,
    n_G,
    ch_G 
)->pulp.LpProblem:
    # -------- Constraint (19) --------
    for i in range(1, t_C + 1):
        MILP += v_C[(i, 0)] == 1, f"milp-(19)-{i}"

    # -------- Constraint (20) --------
    for i in range(1, t_F + 1):
        MILP += v_F[(i, j_F[rho])] >= v_F[(i, 0)] - e_F[i + 1], f"milp-(20)-{i}"

    # -------- Constraint (21) --------
    for i in range(1, t_T + 1):
        for (j, h) in P_prc_T:
            MILP += v_T[(i, j)] >= v_T[(i, h)], f"milp-(21)-T-{i}-({j}-{h})"
    for i in range(1, t_C + 1):
        for (j, h) in P_prc_C[delta_i[i]]:
            MILP += v_C[(i, j)] >= v_C[(i, h)], f"milp-(21)-C-{i}-({j}-{h})"
    for i in range(1, t_F + 1):
        for (j, h) in P_prc_F:
            MILP += v_F[(i, j)] >= v_F[(i, h)], f"milp-(21)-F-{i}-({j}-{h})"

    # -------- Constraint (22) --------
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(v_T[(i, j_T[p])]
                for p in range(1, rho + 1)) == h_T[i], f"milp-(22)-T-{i}"
    for i in range(1, t_C + 1):
        if delta_i[i] != 0:
            MILP += pulp.lpSum(v_C[(i, j_C[delta_i[i]][p])]
                    for p in range(1, rho + 1)) == h_C[i], f"milp-(22)-C-{i}"
        else:
            MILP += h_C[i] == 0, f"milp-(22)-C-{i}"
        # just in case that delta_i[i] = 0
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(v_F[(i, j_F[p])]
                for p in range(1, rho + 1)) == h_F[i], f"milp-(22)-F-{i}"

    # -------- Constraint (23) --------
    # v_X in left side should be v_C, by Zhu
    # n_C(i) =  n_C[delta_i[i]]
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(v_C[(i, j)] for j in range(n_C[delta_i[i]] + 1)) <= \
                2 + 2 * pulp.lpSum(v_C[(i, h)] for h in Cld_C[delta_i[i]][0]), f"milp-(23)-{i}"

    # -------- Constraint (24) --------
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(v_T[(i, j)] for j in range(n_T + 1)) <= \
                2 + 2 * pulp.lpSum(v_T[(i, h)] for h in Cld_T[0]), f"milp-(24)-T-{i}"
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(v_F[(i, j)] for j in range(n_F + 1)) <= \
                2 + 2 * pulp.lpSum(v_F[(i, h)] for h in Cld_F[0]), f"milp-(24)-F-{i}"

    # -------- Constraint (25) --------
    for c in range(1, c_F + 1):
        MILP += clr_F[c] + rho <= ch_G, f"milp-(25)-{c}-1"
        MILP += clr_F[c] + rho >= ch_G - n_star * (1 - sigma[c]), f"milp-(25)-{c}-2"

    # -------- Constraint (26) --------
    # ht_X[i] is not defined, it should be h_X[i], by Zhu
    for i in range(1, t_C + 1):
        MILP += h_C[i] <= ch_G, f"milp-(26)-C-{i}-1"
        MILP += h_C[i] >= ch_G - n_star * (1 - sigma_C[i]), f"milp-(26)-C-{i}-2"
    for i in range(1, t_T + 1):
        MILP += h_T[i] <= ch_G, f"milp-(26)-T-{i}-1"
        MILP += h_T[i] >= ch_G - n_star * (1 - sigma_T[i]), f"milp-(26)-T-{i}-2"
    for i in range(1, t_F + 1):
        MILP += h_F[i] <= ch_G, f"milp-(26)-F-{i}-1"
        MILP += h_F[i] >= ch_G - n_star * (1 - sigma_F[i]), f"milp-(26)-F-{i}-2"

    # -------- Constraint (27) --------
    MILP += pulp.lpSum(v_C[(i, j)] for i in range(1, t_C + 1) for j in range(n_C[delta_i[i]] + 1)) + \
            pulp.lpSum(v_T[(i, j)] for i in range(1, t_T + 1) for j in range(n_T + 1)) + \
            pulp.lpSum(v_F[(i, j)] for i in range(1, t_F + 1) for j in range(n_F + 1)) == n_G, f"milp-(27)"

    # -------- Constraint (28) --------
    for i in range(1, t_C_tilde + 1):
        MILP += h_C[i] >= ch_LB[i] - n_star * delta_chi_F[i], f"milp-(28)-1-{i}"
        MILP += clr_F[i] + rho >= ch_LB[i], f"milp-(28)-2-{i}"
        MILP += h_C[i] <= ch_UB[i], f"milp-(28)-3-{i}"
        MILP += clr_F[i] + rho <= ch_UB[i] + n_star * (1 - delta_chi_F[i]), f"milp-(28)-4-{i}"

    # -------- Constraint (29) --------
    for i in range(t_C_tilde + 1, t_C + 1):
        MILP += h_C[i] >= ch_LB[i], f"milp-(29)-{i}-1"
        MILP += h_C[i] <= ch_UB[i], f"milp-(29)-{i}-2"

    # -------- Constraint (30) --------
    for i in range(1, t_T + 1):
        for c in (I_ge_one + I_ge_two):
            MILP += h_T[i] <= ch_UB[E_C[c]] + \
                n_star * (delta_chi_F[t_C_tilde + i] + 1 - delta_chi_T[(i, c)]), f"milp-(30)-1-{i}-{c}"
            MILP += clr_F[t_C_tilde + i] + rho <= ch_UB[E_C[c]] + \
                n_star * (2 - delta_chi_F[t_C_tilde + i] - delta_chi_T[(i, c)]), f"milp-(30)-2-{i}-{c}"

    # -------- Constraint (31) --------
    for c in (I_ge_one + I_ge_two):
        MILP += pulp.lpSum(sigma[(c, i)] for i in range(1, t_T + 1)) == delta_chi_T[c], f"milp-(31)-{c}"

    # -------- Constraint (32) --------
    for i in range(1, t_T + 1):
        for c in (I_ge_one + I_ge_two):
            MILP += delta_chi_T[(i, c)] >= sigma[(c, i)], f"milp-(32)-1-{i}-{c}"
            MILP += h_T[i] >= ch_LB[E_C[c]] - \
                n_star * (delta_chi_F[t_C_tilde + i] + 1 - sigma[(c, i)]), f"milp-(32)-2-{i}-{c}"
            MILP += clr_F[t_C_tilde + i] + rho >= ch_LB[E_C[c]] - \
                n_star * (2 - delta_chi_F[t_C_tilde + i] - sigma[(c, i)]), f"milp-(32)-3-{i}-{c}"

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
    cs_LB,
    cs_UB,
    dg_LB,
    dg_UB,
    rho
):
    # Define deg_C, deg_T, deg_F
    deg_C = {(i, j): pulp.LpVariable(f"deg_C({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_C + 1) for j in range(n_C[delta_i[i]] + 1)}
    deg_T = {(i, j): pulp.LpVariable(f"deg_T({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_T + 1) for j in range(n_T + 1)}
    deg_F = {(i, j): pulp.LpVariable(f"deg_F({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
            for i in range(1, t_F + 1) for j in range(n_F + 1)}

    # Define delta_dg_C, delta_dg,T, delta_dg_F
    # Here the range of d should be [0, 4] instead of [1, 4]
    delta_dg_C = {(i, j, d): pulp.LpVariable(f"delta_dg_C({i},{j},{d})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_C + 1) for j in range(n_C[delta_i[i]] + 1) for d in range(MAX_VAL + 1)}
    delta_dg_T = {(i, j, d): pulp.LpVariable(f"delta_dg_T({i},{j},{d})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_T + 1) for j in range(n_T + 1) for d in range(MAX_VAL + 1)}
    delta_dg_F = {(i, j, d): pulp.LpVariable(f"delta_dg_F({i},{j},{d})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_F + 1) for j in range(n_F + 1) for d in range(MAX_VAL + 1)}

    # Define dg
    dg = {d: pulp.LpVariable(f"dg({d})", dg_LB[d], dg_UB[d], cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}

    # Define dg_co, dg_C. dg_T
    dg_co = {d: pulp.LpVariable(f"dg_co({d})", 0, cs_UB, cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}
    dg_C = {d: pulp.LpVariable(f"dg_C({d})", 0, cs_UB, cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}
    dg_T = {d: pulp.LpVariable(f"dg_T({d})", 0, cs_UB, cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}

    # Define dg_nc
    dg_nc = {d: pulp.LpVariable(f"dg_nc({d})", 0, n_star - cs_LB, cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}
    
    # Define dg_in
    dg_in = {d: pulp.LpVariable(f"dg_in({d})", 0, n_star - cs_LB, cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}

    # Define dg_ex
    dg_ex = {d: pulp.LpVariable(f"dg_ex({d})", 0, n_star - cs_LB, cat=pulp.LpInteger) for d in range(1, MAX_VAL + 1)}

    # Define dg_Cp, dg_Tp, dg_Fp
    dg_Cp = {(p, d): pulp.LpVariable(f"dg_C({p})({d})", 0, n_star - cs_LB, cat=pulp.LpInteger)
            for d in range(1, MAX_VAL + 1) for p in range(1, rho + 1)}
    dg_Tp = {(p, d): pulp.LpVariable(f"dg_T({p})({d})", 0, n_star - cs_LB, cat=pulp.LpInteger)
            for d in range(1, MAX_VAL + 1) for p in range(1, rho + 1)}
    dg_Fp = {(p, d): pulp.LpVariable(f"dg_F({p})({d})", 0, n_star - cs_LB, cat=pulp.LpInteger)
            for d in range(1, MAX_VAL + 1) for p in range(1, rho + 1)}


    return deg_C, deg_T, deg_F, delta_dg_C, delta_dg_T, delta_dg_F, dg, dg_co, dg_C, dg_T, dg_nc, dg_in, dg_ex, dg_Cp, dg_Tp, dg_Fp

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
    Cld_T,
    Cld_C,
    Cld_F,
    Dsn_T,
    Dsn_C,
    Dsn_F,
    delta_i,
    dg4_UB_nc,
    I_ge_one_plus,
    I_ge_one_minus,
    I_ge_two_plus,
    I_ge_two_minus,
    rho,
    # Binary Variables
    delta_dg_C,
    delta_dg_T,
    delta_dg_F,
    e_T,
    e_F,
    v_F,
    v_C,
    v_T,
    delta_chi_T,
    delta_chi_F,
    # Integer Variables
    deg_C,
    deg_T,
    deg_F,
    dg,
    dg_co,
    dg_C,
    dg_T,
    dg_nc,
    dg_in,
    dg_ex,
    dg_Cp,
    dg_Tp,
    dg_Fp,
    deg_tilde_C_minus,
    deg_tilde_C_plus
):
    # -------- Constraint (33) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += v_C[(i, j)] + pulp.lpSum(v_C[(i, h)] 
                        for h in Cld_C[delta_i[i]][j]) == deg_C[(i, j)], f"milp-(33)-C-{i}-{j}"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += v_T[(i, j)] + pulp.lpSum(v_T[(i, h)] 
                        for h in Cld_T[j]) == deg_T[(i, j)], f"milp-(33)-T-{i}-{j}"
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            MILP += v_F[(i, j)] + pulp.lpSum(v_F[(i, h)] 
                        for h in Cld_F[j]) == deg_F[(i, j)], f"milp-(33)-F-{i}-{j}"

    # -------- Constraint (34) --------
    for i in range(1, t_C_tilde + 1):
        MILP += deg_tilde_C_minus[i] + deg_tilde_C_plus[i] + \
                pulp.lpSum(delta_chi_T[c] for c in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
                pulp.lpSum(delta_chi_T[c] for c in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
                delta_chi_F[i] + \
                pulp.lpSum(v_C[(i, h)] for h in Cld_C[delta_i[i]][0]) == deg_C[(i, 0)], f"milp-(34)-{i}"

    # -------- Constraint (35) --------
    for i in range(t_C_tilde + 1, t_C + 1):
        MILP += deg_tilde_C_minus[i] + deg_tilde_C_plus[i] + \
                pulp.lpSum(delta_chi_T[c] for c in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
                pulp.lpSum(delta_chi_T[c] for c in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
                pulp.lpSum(v_C[(i, h)] for h in Cld_C[delta_i[i]][0]) == deg_C[(i, 0)], f"milp-(35)-{i}"

    # -------- Constraint (36) --------
    MILP += e_T[1] == 0, f"milp-(36)-T1"
    MILP += e_T[t_T + 1] == 0, f"milp-(36)-T2"
    for i in range(1, t_T + 1):
        MILP += 2 * v_T[(i, 0)] + pulp.lpSum(v_T[(i, h)] for h in Cld_T[0]) + \
                delta_chi_F[t_C_tilde + i] == deg_T[(i, 0)], f"milp-(36)-{i}"

    # -------- Constraint (37) --------
    for i in range(1, t_F + 1):
        MILP += v_F[(i, 0)] + e_F[i + 1] + \
                pulp.lpSum(v_F[(i, h)] for h in Cld_F[0]) == deg_F[(i, 0)], f"milp-(37)-{i}"

    # -------- Constraint (38) --------
    for i in range(1, t_T + 1):
        for j in range(n_T + 1):
            MILP += pulp.lpSum(delta_dg_T[(i, j, d)] 
                    for d in range(MAX_VAL + 1)) == 1, f"milp-(38)-T-1-{i}-{j}"
            MILP += pulp.lpSum(d * delta_dg_T[(i, j, d)] 
                    for d in range(1, MAX_VAL + 1)) == deg_T[(i, j)], f"milp-(38)-T-2-{i}-{j}"
    for i in range(1, t_C + 1):
        for j in range(n_C[delta_i[i]] + 1):
            MILP += pulp.lpSum(delta_dg_C[(i, j, d)] 
                    for d in range(MAX_VAL + 1)) == 1, f"milp-(38)-C-1-{i}-{j}"
            MILP += pulp.lpSum(d * delta_dg_C[(i, j, d)] 
                    for d in range(1, MAX_VAL + 1)) == deg_C[(i, j)], f"milp-(38)-C-2-{i}-{j}"
    for i in range(1, t_F + 1):
        for j in range(n_F + 1):
            MILP += pulp.lpSum(delta_dg_F[(i, j, d)] 
                    for d in range(MAX_VAL + 1)) == 1, f"milp-(38)-F-1-{i}-{j}"
            MILP += pulp.lpSum(d * delta_dg_F[(i, j, d)] 
                    for d in range(1, MAX_VAL + 1)) == deg_F[(i, j)], f"milp-(38)-F-2-{i}-{j}"

    # -------- Constraint (39) --------
    # It should be i \in [1, t_X] instead of [1, t_C]
    for d in range(1, MAX_VAL + 1):
        for p in range(1, rho + 1):
            MILP += pulp.lpSum(delta_dg_C[(i, j, d)] for i in range(1, t_C + 1) 
                    for j in Dsn_C[delta_i[i]][p]) == dg_Cp[(p, d)], f"milp-(39)-C-{d}-{p}"
            MILP += pulp.lpSum(delta_dg_T[(i, j, d)] for i in range(1, t_T + 1) 
                    for j in Dsn_T[p]) == dg_Tp[(p, d)], f"milp-(39)-T-{d}-{p}"
            MILP += pulp.lpSum(delta_dg_F[(i, j, d)] for i in range(1, t_F + 1) 
                    for j in Dsn_F[p]) == dg_Fp[(p, d)], f"milp-(39)-F-{d}-{p}"
            
    # -------- Constraint (40) --------
    for d in range(1, MAX_VAL + 1):
        MILP += pulp.lpSum(delta_dg_C[(i, 0, d)] for i in range(1, t_C + 1)) == dg_C[d], f"milp-(40)-1-{d}"
        MILP += pulp.lpSum(delta_dg_T[(i, 0, d)] for i in range(1, t_T + 1)) == dg_T[d], f"milp-(40)-2-{d}"
        MILP += pulp.lpSum(delta_dg_F[(i, 0, d)] for i in range(1, t_F + 1)) == dg_in[d], f"milp-(40)-3-{d}"
        MILP += dg_in[d] + pulp.lpSum(dg_Cp[(p, d)] + dg_Tp[(p, d)] +  dg_Fp[(p, d)]
                for p in range(1, rho + 1)) == dg_nc[d], f"milp-(40)-4-{d}"

        MILP += dg_C[d] + dg_T[d] == dg_co[d], f"milp-(40)-5-{d}"

    # -------- Constraint (41) --------
    MILP += dg_nc[4] <= dg4_UB_nc, f"milp-(41)"

    return MILP

# -------- A.5 Assigning Multiplicity --------
def prepare_variables_multiplicity(
    t_C,
    t_C_tilde,
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
    m_max, 
    m_max_co,
    m_max_nc,
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

    # Define beta_C, beta_T, beta_F
    beta_C_temp = {(i, j): pulp.LpVariable(f"beta_C({i},{j})", 0, MAX_BOND, cat=pulp.LpInteger)
                for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    beta_C.update(beta_C_temp)
    beta_T_temp = {(i, j): pulp.LpVariable(f"beta_T({i},{j})", 0, MAX_BOND, cat=pulp.LpInteger)
                for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    beta_T.update(beta_T_temp)
    beta_F_temp = {(i, j): pulp.LpVariable(f"beta_F({i},{j})", 0, MAX_BOND, cat=pulp.LpInteger)
                for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    beta_F.update(beta_F_temp)

    # Define beta_plus, beta_minus
    beta_plus = {c: pulp.LpVariable(f"beta_plus({c})", 0, MAX_BOND, cat=pulp.LpInteger)
                for c in (I_ge_two + I_ge_one)}
    beta_minus = {c: pulp.LpVariable(f"beta_minus({c})", 0, MAX_BOND, cat=pulp.LpInteger)
                for c in (I_ge_two + I_ge_one)}

    # Define beta_in
    beta_in = {c: pulp.LpVariable(f"beta_in({c})", 0, MAX_BOND, cat=pulp.LpInteger)
                for c in range(1, c_F + 1)}

    # Define delta_beta_T, delta_beta_F
    delta_beta_T = {(i, m): pulp.LpVariable(f"delta_beta_T({i},{m})", 0, 1, cat=pulp.LpBinary)
                    for i in range(2, t_T + 1) for m in range(MAX_BOND + 1)}
    delta_beta_F = {(i, m): pulp.LpVariable(f"delta_beta_F({i},{m})", 0, 1, cat=pulp.LpBinary)
                    for i in range(2, t_F + 1) for m in range(MAX_BOND + 1)}

    # Define delta_beta_C
    delta_beta_C = {(i, m): pulp.LpVariable(f"delta_beta_C({i},{m})", 0, 1, cat=pulp.LpBinary)
                    for i in (I_equal_one + I_zero_one + I_ge_one) for m in range(MAX_BOND + 1)}

    # Define delta_beta_C, delta_beta_T, delta_beta_F
    delta_beta_C_temp = {(i, j, m): pulp.LpVariable(f"delta_beta_C({i},{j},{m})", 0, 1, cat=pulp.LpBinary)
                        for i in range(1, t_C + 1) for j in range(n_C[delta_i[i]] + 1) for m in range(MAX_BOND + 1)}
    delta_beta_C.update(delta_beta_C_temp)
    delta_beta_T_temp = {(i, j, m): pulp.LpVariable(f"delta_beta_T({i},{j},{m})", 0, 1, cat=pulp.LpBinary)
                        for i in range(1, t_T + 1) for j in range(n_T + 1) for m in range(MAX_BOND + 1)}
    delta_beta_T.update(delta_beta_T_temp)
    delta_beta_F_temp = {(i, j, m): pulp.LpVariable(f"delta_beta_F({i},{j},{m})", 0, 1, cat=pulp.LpBinary)
                        for i in range(1, t_F + 1) for j in range(n_F + 1) for m in range(MAX_BOND + 1)}
    delta_beta_F.update(delta_beta_F_temp)

    # Define delta_beta_plus, delta_beta_in
    # The range should be [0, 1] instead of [0, 3], by Zhu
    delta_beta_plus = {(c, m): pulp.LpVariable(f"delta_beta_plus({c},{m})", 0, 1, cat=pulp.LpBinary)
                        for c in (I_ge_two + I_ge_one) for m in range(MAX_BOND + 1)}
    delta_beta_minus = {(c, m): pulp.LpVariable(f"delta_beta_minus({c},{m})", 0, 1, cat=pulp.LpBinary)
                        for c in (I_ge_two + I_ge_one) for m in range(MAX_BOND + 1)}
    
    # Define delta_beta_in
    # The range should be [0, 1] instead of [0, 3], by Zhu
    delta_beta_in = {(c, m): pulp.LpVariable(f"delta_beta_in({c},{m})", 0, 1, cat=pulp.LpBinary)
                        for c in range(1, c_F + 1) for m in range(MAX_BOND + 1)}
    
    # Define bd
    bd = {m: pulp.LpVariable(f"bd({m})", 0, m_max, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}

    # Define bd_co
    bd_co = {m: pulp.LpVariable(f"bd_co({m})", 0, m_max_co, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}

    # Define bd_in
    bd_in = {m: pulp.LpVariable(f"bd_in({m})", 0, m_max_nc, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}

    # Define bd_ex
    bd_ex = {m: pulp.LpVariable(f"bd_ex({m})", 0, n_star, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}

    # Define bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF
    bd_T = {m: pulp.LpVariable(f"bd_T({m})", 0, m_max_co, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_C = {m: pulp.LpVariable(f"bd_C({m})", 0, m_max_co, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_CT = {m: pulp.LpVariable(f"bd_CT({m})", 0, m_max_co, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_TC = {m: pulp.LpVariable(f"bd_TC({m})", 0, m_max_co, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_F = {m: pulp.LpVariable(f"bd_F({m})", 0, m_max_nc, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_CF = {m: pulp.LpVariable(f"bd_CF({m})", 0, m_max_nc, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}
    bd_TF = {m: pulp.LpVariable(f"bd_TF({m})", 0, m_max_nc, cat=pulp.LpInteger) for m in range(1, MAX_BOND + 1)}

    # Define bd_Cp, bd_Tp, bd_Fp
    bd_Cp = {(p, m): pulp.LpVariable(f"bd_C({p})({m})", 0, n_star, cat=pulp.LpInteger)
            for p in range(1, rho + 1) for m in range(1, MAX_BOND + 1)}
    bd_Tp = {(p, m): pulp.LpVariable(f"bd_T({p})({m})", 0, n_star, cat=pulp.LpInteger) 
            for p in range(1, rho + 1) for m in range(1, MAX_BOND + 1)}
    bd_Fp = {(p, m): pulp.LpVariable(f"bd_F({p})({m})", 0, n_star, cat=pulp.LpInteger)
            for p in range(1, rho + 1) for m in range(1, MAX_BOND + 1)}
    
    return beta_T, beta_F, beta_C, beta_plus, beta_minus, beta_in, delta_beta_T, delta_beta_F, delta_beta_C, \
        delta_beta_plus, delta_beta_minus, delta_beta_in, bd, bd_co, bd_in, bd_ex, \
        bd_C, bd_T, bd_CT, bd_TC, bd_F, bd_CF, bd_TF, bd_Cp, bd_Tp, bd_Fp

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
    Dsn_C,
    Dsn_T,
    Dsn_F,
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
    # Binary Variables
    e_C,
    e_T,
    e_F,
    v_F,
    v_C,
    v_T,
    delta_chi_T,
    delta_chi_F,
    delta_beta_C,
    delta_beta_T,
    delta_beta_F,
    delta_beta_plus,
    delta_beta_minus,
    delta_beta_in,
    # Integer Variables
    beta_C,
    beta_T,
    beta_F,
    beta_plus, 
    beta_minus, 
    beta_in,
    bd,
    bd_co,
    bd_in,
    bd_ex,
    bd_C,
    bd_T,
    bd_CT,
    bd_TC,
    bd_F,
    bd_CF,
    bd_TF,
    bd_Cp,
    bd_Tp,
    bd_Fp
):
    # -------- Constraint (42) --------
    for i in (I_equal_one + I_zero_one + I_ge_one):
        MILP += beta_C[i] >= e_C[i], f"milp-(42)-{i}-1"
        MILP += beta_C[i] <= MAX_BOND * e_C[i], f"milp-(42)-{i}-2"

    # -------- Constraint (43) --------
    for i in range(2, t_T + 1):
        MILP += beta_T[i] >= e_T[i], f"milp-(43)-T-{i}-1"
        MILP += beta_T[i] <= MAX_BOND * e_T[i], f"milp-(43)-T-{i}-2"
    for i in range(2, t_F + 1):
        MILP += beta_F[i] >= e_F[i], f"milp-(43)-F-{i}-1"
        MILP += beta_F[i] <= MAX_BOND * e_F[i], f"milp-(43)-F-{i}-2"

    # -------- Constraint (44) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += beta_C[(i, j)] >= v_C[(i, j)], f"milp-(44)-C-{i}-{j}-1"
            MILP += beta_C[(i, j)] <= MAX_BOND * v_C[(i, j)], f"milp-(44)-C-{i}-{j}-2"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += beta_T[(i, j)] >= v_T[(i, j)], f"milp-(44)-T-{i}-{j}-1"
            MILP += beta_T[(i, j)] <= MAX_BOND * v_T[(i, j)], f"milp-(44)-T-{i}-{j}-2"
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            MILP += beta_F[(i, j)] >= v_F[(i, j)], f"milp-(44)-F-{i}-{j}-1"
            MILP += beta_F[(i, j)] <= MAX_BOND * v_F[(i, j)], f"milp-(44)-F-{i}-{j}-2"

    # -------- Constraint (45) --------
    for c in (I_ge_two + I_ge_one):
        MILP += beta_plus[c] >= delta_chi_T[c], f"milp-(45)-1-{c}-1"
        MILP += beta_plus[c] <= MAX_BOND * delta_chi_T[c], f"milp-(45)-1-{c}-2"
        MILP += beta_minus[c] >= delta_chi_T[c], f"milp-(45)-2-{c}-1"
        MILP += beta_minus[c] <= MAX_BOND * delta_chi_T[c], f"milp-(45)-2-{c}-2"

    # -------- Constraint (46) --------
    for c in range(1, c_F + 1):
        MILP += beta_in[c] >= delta_chi_F[c], f"milp-(46)-{c}-1"
        MILP += beta_in[c] <= MAX_BOND * delta_chi_F[c], f"milp-(46)-{c}-2"

    # -------- Constraint (47) --------
    for i in range(2, t_T + 1):
        MILP += pulp.lpSum(delta_beta_T[(i, m)] 
            for m in range(MAX_BOND + 1)) == 1, f"milp-(47)-T-1-{i}"
        MILP += pulp.lpSum(m * delta_beta_T[(i, m)] 
            for m in range(MAX_BOND + 1)) == beta_T[i], f"milp-(47)-T-2-{i}"
    for i in range(2, t_F + 1):
        MILP += pulp.lpSum(delta_beta_F[(i, m)] 
            for m in range(MAX_BOND + 1)) == 1, f"milp-(47)-F-1-{i}"
        MILP += pulp.lpSum(m * delta_beta_F[(i, m)] 
            for m in range(MAX_BOND + 1)) == beta_F[i], f"milp-(47)-F-2-{i}"

    # -------- Constraint (48) --------
    # Should be (I_equal_one + I_zero_one + I_ge_one) instead of [1, m_C], by Zhu
    for i in (I_equal_one + I_zero_one + I_ge_one):
        MILP += pulp.lpSum(delta_beta_C[(i, m)] 
            for m in range(MAX_BOND + 1)) == 1, f"milp-(48)-1-{i}"
        MILP += pulp.lpSum(m * delta_beta_C[(i, m)] 
            for m in range(MAX_BOND + 1)) == beta_C[i], f"milp-(48)-2-{i}"

    # -------- Constraint (49) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += pulp.lpSum(delta_beta_C[(i, j, m)] 
                for m in range(MAX_BOND + 1)) == 1, f"milp-(49)-C-1-{i}-{j}"
            MILP += pulp.lpSum(m * delta_beta_C[(i, j, m)] 
                for m in range(MAX_BOND + 1)) == beta_C[(i, j)], f"milp-(49)-C-2-{i}-{j}"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += pulp.lpSum(delta_beta_T[(i, j, m)] 
                for m in range(MAX_BOND + 1)) == 1, f"milp-(49)-T-1-{i}-{j}"
            MILP += pulp.lpSum(m * delta_beta_T[(i, j, m)] 
                for m in range(MAX_BOND + 1)) == beta_T[(i, j)], f"milp-(49)-T-2-{i}-{j}"
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            MILP += pulp.lpSum(delta_beta_F[(i, j, m)] 
                for m in range(MAX_BOND + 1)) == 1, f"milp-(49)-F-1-{i}-{j}"
            MILP += pulp.lpSum(m * delta_beta_F[(i, j, m)] 
                for m in range(MAX_BOND + 1)) == beta_F[(i, j)], f"milp-(49)-F-2-{i}-{j}"

    # -------- Constraint (50) --------
    for c in (I_ge_two + I_ge_one):
        MILP += pulp.lpSum(delta_beta_plus[(c, m)]
            for m in range(MAX_BOND + 1)) == 1, f"milp-(50)-1-{c}"
        MILP += pulp.lpSum(m * delta_beta_plus[(c, m)]
            for m in range(MAX_BOND + 1)) == beta_plus[c], f"milp-(50)-2-{c}"
    for c in (I_ge_two + I_ge_one):
        MILP += pulp.lpSum(delta_beta_minus[(c, m)]
            for m in range(MAX_BOND + 1)) == 1, f"milp-(50)-3-{c}"
        MILP += pulp.lpSum(m * delta_beta_minus[(c, m)]
            for m in range(MAX_BOND + 1)) == beta_minus[c], f"milp-(50)-4-{c}"
    for c in range(1, c_F + 1):
        MILP += pulp.lpSum(delta_beta_in[(c, m)]
            for m in range(MAX_BOND + 1)) == 1, f"milp-(50)-5-{c}"
        MILP += pulp.lpSum(m * delta_beta_in[(c, m)]
            for m in range(MAX_BOND + 1)) == beta_in[c], f"milp-(50)-6-{c}"

    # -------- Constraint (51) --------
    for p in range(1, rho + 1):
        for m in range(1, MAX_BOND + 1):
            MILP += pulp.lpSum(delta_beta_C[(i, j, m)] for i in range(1, t_C + 1) 
                    for j in Dsn_C[delta_i[i]][p]) == bd_Cp[(p, m)], f"milp-(51)-C-{p}-{m}"
            MILP += pulp.lpSum(delta_beta_T[(i, j, m)] for i in range(1, t_T + 1) 
                    for j in Dsn_T[p]) == bd_Tp[(p, m)], f"milp-(51)-T-{p}-{m}"
            MILP += pulp.lpSum(delta_beta_F[(i, j, m)] for i in range(1, t_F + 1) 
                    for j in Dsn_F[p]) == bd_Fp[(p, m)], f"milp-(51)-F-{p}-{m}"

    # -------- Constraint (52) --------
    for m in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum(delta_beta_C[(i, m)] for i in (I_equal_one + I_zero_one + I_ge_one)) == bd_C[m], f"milp-(52)-1-{m}"
        MILP += pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)) == bd_T[m], f"milp-(52)-2-{m}"
        MILP += pulp.lpSum(delta_beta_plus[(c, m)] for c in (I_ge_two + I_ge_one)) == bd_CT[m], f"milp-(52)-3-{m}"
        MILP += pulp.lpSum(delta_beta_minus[(c, m)] for c in (I_ge_two + I_ge_one)) == bd_TC[m], f"milp-(52)-4-{m}"
        MILP += bd_C[m] + bd_T[m] + bd_CT[m] + bd_TC[m] == bd_co[m], f"milp-(52)-5-{m}"
        MILP += pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)) == bd_F[m], f"milp-(52)-6-{m}"
        MILP += pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)) == bd_CF[m], f"milp-(52)-7-{m}"
        MILP += pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)) == bd_TF[m], f"milp-(52)-8-{m}"
        MILP += bd_F[m] + bd_TF[m] + bd_CF[m] == bd_in[m], f"milp-(52)-9-{m}"
        MILP += pulp.lpSum(bd_Cp[(p, m)] + bd_Tp[(p, m)] + bd_Fp[(p, m)] for p in range(1, rho + 1)) == bd_ex[m], f"milp-(52)-10-{m}"
        MILP += bd_co[m] + bd_in[m] + bd_ex[m] == bd[m], f"milp-(52)-11-{m}"

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
    MAX_CODE,
    MAX_CODE_co,
    MAX_CODE_nc,
    Lambda,
    Lambda_co,
    Lambda_nc,
    epsilon,
    na_LB,
    na_UB,
    na_LB_co,
    na_UB_co,
    na_LB_nc,
    na_UB_nc,
    rho
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
    alpha_C = {(i, 0): pulp.LpVariable(f"alpha_C({i},{0})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                for i in range(1, t_C + 1)}
    alpha_C_temp = {(i, j): pulp.LpVariable(f"alpha_C({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    alpha_C.update(alpha_C_temp)
    alpha_T = {(i, 0): pulp.LpVariable(f"alpha_T({i},{0})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    alpha_T_temp = {(i, j): pulp.LpVariable(f"alpha_T({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    alpha_T.update(alpha_T_temp)
    alpha_F = {(i, j): pulp.LpVariable(f"alpha_F({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                for i in range(1, t_F + 1) for j in range(n_F + 1)}

    # Define delta_alpha_C, delta_alpha_T, delta_alpha_F
    delta_alpha_C = {(i, j, mu): pulp.LpVariable(f"delta_alpha_C({i},{j},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_C + 1) for j in range(n_C[delta_i[i]] + 1) for mu in (Lambda + [epsilon])}
    delta_alpha_T = {(i, j, mu): pulp.LpVariable(f"delta_alpha_T({i},{j},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_T + 1) for j in range(n_T + 1) for mu in (Lambda + [epsilon])}
    delta_alpha_F = {(i, j, mu): pulp.LpVariable(f"delta_alpha_F({i},{j},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_F + 1) for j in range(n_F + 1) for mu in (Lambda + [epsilon])}

    # Define MASS
    MASS = pulp.LpVariable(f"Mass", cat=pulp.LpInteger)

    # Define n_H
    n_H = pulp.LpVariable(f"n_H", 0, 4 * n_star, cat=pulp.LpInteger)

    # Define na
    na = {atom: pulp.LpVariable(f"na({atom})", na_LB[atom], na_UB[atom], cat=pulp.LpInteger)
            for atom in Lambda}

    # Define na_co, na_C, na_T
    na_co = {atom: pulp.LpVariable(f"na_co({atom})", na_LB_co[atom], na_UB_co[atom], cat=pulp.LpInteger) 
            for atom in Lambda_co}
    na_C = {atom: pulp.LpVariable(f"na_C({atom})", 0, na_UB_co[atom], cat=pulp.LpInteger) 
            for atom in Lambda_co}
    na_T = {atom: pulp.LpVariable(f"na_T({atom})", 0, na_UB_co[atom], cat=pulp.LpInteger)
            for atom in Lambda_co}

    # Define na_nc, na_in, na_Cp, na_Tp, na_Fp
    na_nc = {atom: pulp.LpVariable(f"na_nc({atom})", na_LB_nc[atom], na_UB_nc[atom], cat=pulp.LpInteger)
            for atom in Lambda_nc}
    na_in = {atom: pulp.LpVariable(f"na_in({atom})", 0, na_UB_nc[atom], cat=pulp.LpInteger)
            for atom in Lambda_nc}
    na_ex = {atom: pulp.LpVariable(f"na_ex({atom})", 0, na_UB_nc[atom], cat=pulp.LpInteger)
            for atom in Lambda_nc} # Not defined in the draft but used in constraint (65)
    na_Cp = {(p, atom): pulp.LpVariable(f"na_C({p})({atom})", 0, na_UB_nc[atom], cat=pulp.LpInteger)
            for p in range(1, rho + 1) for atom in Lambda_nc}
    na_Tp = {(p, atom): pulp.LpVariable(f"na_T({p})({atom})", 0, na_UB_nc[atom], cat=pulp.LpInteger)
            for p in range(1, rho + 1) for atom in Lambda_nc}
    na_Fp = {(p, atom): pulp.LpVariable(f"na_F({p})({atom})", 0, na_UB_nc[atom], cat=pulp.LpInteger)
            for p in range(1, rho + 1) for atom in Lambda_nc}

    return beta_CT, beta_TC, beta_CF, beta_TF, alpha_C, alpha_T, alpha_F, delta_alpha_C, delta_alpha_T, delta_alpha_F, MASS, n_H, \
        na, na_co, na_C, na_T, na_nc, na_in, na_ex, na_Cp, na_Tp, na_Fp

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
    Cld_C,
    Cld_T,
    Cld_F,
    Dsn_C,
    Dsn_T,
    Dsn_F,
    val,
    mass,
    Lambda,
    Lambda_co,
    Lambda_nc,
    Code_Lambda,
    Code_Lambda_co,
    Code_Lambda_nc,
    Lambda_star,
    epsilon,
    rho,
    # Binary Variables
    e_T,
    e_F,
    delta_chi_T,
    delta_chi_F,
    delta_alpha_C,
    delta_alpha_T,
    delta_alpha_F,
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
    n_H,
    na,
    na_co,
    na_C,
    na_T,
    na_nc,
    na_in,
    na_ex,
    na_Cp,
    na_Tp,
    na_Fp,
    bd_co,
    bd_in,
    bd_ex
):
    # -------- Constraint (53) --------
    # for c in (I_ge_one + I_ge_two):
    #     MILP += pulp.lpSum(beta_CT[i] for i in range(1, t_T + 1)) == beta_plus[c], f"milp-(53)-3-{c}"
    #     MILP += pulp.lpSum(beta_TC[i] for i in range(1, t_T + 1)) == beta_minus[c], f"milp-(53)-4-{c}"
    # ##############################################################################################################
    for c in (I_ge_one + I_ge_two):
        for i in range(1, t_T + 1):
            MILP += beta_CT[i] >= beta_plus[c] - MAX_BOND * (e_T[i] - delta_chi_T[(i, c)] + 1), f"milp-(53)-1-{c}-{i}-1"
            MILP += beta_CT[i] <= beta_plus[c] + MAX_BOND * (e_T[i] - delta_chi_T[(i, c)] + 1), f"milp-(53)-1-{c}-{i}-2"
        for i in range(1, t_T + 1):
            MILP += beta_TC[i] >= beta_minus[c] - MAX_BOND * (e_T[i + 1] - delta_chi_T[(i, c)] + 1), f"milp-(53)-2-{c}-{i}-1"
            MILP += beta_TC[i] <= beta_minus[c] + MAX_BOND * (e_T[i + 1] - delta_chi_T[(i, c)] + 1), f"milp-(53)-2-{c}-{i}-2"
    # The constraint sum(beta_CT(i)) = beta_plus(c) looks very strange for me, it means beta_plus(c) has the same value for all c.
    # ##############################################################################################################
    # MILP += pulp.lpSum(beta_CT[i] for i in range(1, t_T + 1)) == \
    #         pulp.lpSum(beta_plus[c] for c in (I_ge_one + I_ge_two)), f"milp-(53)-3"
    # MILP += pulp.lpSum(beta_TC[i] for i in range(1, t_T + 1)) == \
    #         pulp.lpSum(beta_minus[c] for c in (I_ge_one + I_ge_two)), f"milp-(53)-4"

    # -------- Constraint (54) --------
    # for c in range(1, t_C_tilde + 1):
    #     MILP += pulp.lpSum(beta_CF[i] for i in range(1, t_F+ 1)) == beta_in[c], f"milp-(54)-2-{c}"
    # for c in range(t_C_tilde + 1, c_F + 1):
    #     MILP += pulp.lpSum(beta_TF[i] for i in range(1, t_F + 1)) == beta_in[c], f"milp-(54)-4-{c}"
    ##############################################################################################################
    for c in range(1, t_C_tilde + 1):
        for i in range(1, t_F + 1):
            MILP += beta_CF[i] >= beta_in[c] - MAX_BOND * (e_F[i] - delta_chi_F[(i, c)] + 1), f"milp-(54)-1-{c}-{i}-1"
            MILP += beta_CF[i] <= beta_in[c] + MAX_BOND * (e_F[i] - delta_chi_F[(i, c)] + 1), f"milp-(54)-1-{c}-{i}-2"
    for c in range(t_C_tilde + 1, c_F + 1):
        for i in range(1, t_F + 1):
            MILP += beta_TF[i] >= beta_in[c] - MAX_BOND * (e_F[i] - delta_chi_F[(i, c)] + 1), f"milp-(54)-3-{c}-{i}-1"
            MILP += beta_TF[i] <= beta_in[c] + MAX_BOND * (e_F[i] - delta_chi_F[(i, c)] + 1), f"milp-(54)-3-{c}-{i}-2"
    ##############################################################################################################
    # MILP += pulp.lpSum(beta_CF[i] for i in range(1, t_F + 1)) == \
    #         pulp.lpSum(beta_in[c] for c in range(1, t_C_tilde + 1)), f"milp-(54)-2"
    # MILP += pulp.lpSum(beta_TF[i] for i in range(1, t_F + 1)) == \
    #         pulp.lpSum(beta_in[c] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-(54)-4"

    # -------- Constraint (55) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_alpha_C[(i, 0, atom)] for atom in (Lambda_co + [epsilon])) == 1, \
                f"milp-(55)-C-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_co[atom] * delta_alpha_C[(i, 0, atom)] 
                for atom in (Lambda_co + [epsilon])) == alpha_C[(i, 0)], f"milp-(55)-C-2-{i}"
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(delta_alpha_T[(i, 0, atom)] for atom in (Lambda_co + [epsilon])) == 1, \
                f"milp-(55)-T-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_co[atom] * delta_alpha_T[(i, 0, atom)] 
                for atom in (Lambda_co + [epsilon])) == alpha_T[(i, 0)], f"milp-(55)-T-2-{i}"

    # -------- Constraint (56) --------
    for i in range(1, t_F + 1):
        for j in range(n_F + 1):
            MILP += pulp.lpSum(delta_alpha_F[(i, j, atom)] for atom in (Lambda_nc + [epsilon])) == 1, \
                    f"milp-(56)-F-1-{i}-{j}"
            MILP += pulp.lpSum(Code_Lambda_nc[atom] * delta_alpha_F[(i, j, atom)] 
                    for atom in (Lambda_nc + [epsilon])) == alpha_F[(i, j)], f"milp-(56)-2-{i}-{j}"

    # -------- Constraint (57) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += pulp.lpSum(delta_alpha_C[(i, j, atom)] for atom in (Lambda_nc + [epsilon])) == 1, \
                    f"milp-(57)-C-1-{i}-{j}"
            MILP += pulp.lpSum(Code_Lambda_nc[atom] * delta_alpha_C[(i, j, atom)] 
                    for atom in (Lambda_nc + [epsilon])) == alpha_C[(i, j)], f"milp-(57)-C-2-{i}-{j}"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += pulp.lpSum(delta_alpha_T[(i, j, atom)] for atom in (Lambda_nc + [epsilon])) == 1, \
                    f"milp-(57)-T-1-{i}-{j}"
            MILP += pulp.lpSum(Code_Lambda_nc[atom] * delta_alpha_T[(i, j, atom)] 
                    for atom in (Lambda_nc + [epsilon])) == alpha_T[(i, j)], f"milp-(57)-T-2-{i}-{j}"
    # For X=F, this part is same as (56)
    # for i in range(1, t_F + 1):
    #     for j in range(1, n_F + 1):
    #         MILP += pulp.lpSum(delta_alpha_F[(i, j, atom)] for atom in (Lambda_nc + [epsilon])) == 1, \
    #                 f"milp-(57)-F-1-{i}-{j}"
    #         MILP += pulp.lpSum(Code_Lambda[atom] * delta_alpha_F[(i, j, atom)] 
    #                 for atom in (Lambda_nc + [epsilon])) == alpha_F[(i, j)], f"milp-(57)-2-{i}-{j}"

    # -------- Constraint (58) --------
    for i in range(1, t_C_tilde + 1):
        MILP += pulp.lpSum(beta_C[j] for j in range(1, len(E_C) + 1) 
                    if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
                pulp.lpSum(beta_plus[c] for c in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
                pulp.lpSum(beta_minus[c] for c in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
                beta_in[i] + \
                pulp.lpSum(beta_C[(i, h)] for h in Cld_C[delta_i[i]][0]) <= \
                pulp.lpSum(val[atom] * delta_alpha_C[(i, 0, atom)] for atom in Lambda_co), f"milp-(58)-{i}"
    # For the first beta_C[i] in the left side, it should be beta_C[j], by Zhu

    # -------- Constraint (59) --------
    for i in range(t_C_tilde + 1, t_C + 1):
        MILP += pulp.lpSum(beta_C[j] for j in range(1, len(E_C) + 1) 
                    if j in (I_equal_one + I_zero_one + I_ge_one) and (E_C[j][1] == i or E_C[j][2] == i)) + \
                pulp.lpSum(beta_plus[c] for c in (I_ge_two_plus[i] + I_ge_one_plus[i])) + \
                pulp.lpSum(beta_minus[c] for c in (I_ge_two_minus[i] + I_ge_one_minus[i])) + \
                pulp.lpSum(beta_C[(i, h)] for h in Cld_C[delta_i[i]][0]) <= \
                pulp.lpSum(val[atom] * delta_alpha_C[(i, 0, atom)] for atom in Lambda_co), f"milp-(59)-{i}"

    # -------- Constraint (60) --------
    MILP += beta_T[1] == 0, f"milp-(60)-T1"
    MILP += beta_T[t_T + 1] == 0, f"milp-(60)-T2"
    for i in range(1, t_T + 1):
        MILP += beta_T[i] + beta_T[i + 1] + pulp.lpSum(beta_T[(i, h)] for h in Cld_T[0]) + \
                beta_CT[i] + beta_TC[i] + beta_in[t_C_tilde + i] <= \
                pulp.lpSum(val[atom] * delta_alpha_T[(i, 0, atom)] for atom in Lambda_co), f"milp-(60)-{i}"

    # -------- Constraint (61) --------
    MILP += beta_F[1] == 0, f"milp-(61)-F1"
    MILP += beta_F[t_F + 1] == 0, f"milp-(61)-F2"
    for i in range(1, t_F + 1):
        MILP += beta_F[i] + beta_F[i + 1] + pulp.lpSum(beta_F[(i, h)] for h in Cld_F[0]) + \
                beta_CF[i] + beta_TF[i] <= \
                pulp.lpSum(val[atom] * delta_alpha_F[(i, 0, atom)] for atom in Lambda_nc), f"milp-(61)-{i}"

    # -------- Constraint (62) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += beta_C[(i, j)] + pulp.lpSum(beta_C[(i, h)] for h in Cld_C[delta_i[i]][j]) <= \
                    pulp.lpSum(val[atom] * delta_alpha_C[(i, j, atom)] for atom in Lambda_nc), f"milp-(62)-C-{i}-{j}"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += beta_T[(i, j)] + pulp.lpSum(beta_T[(i, h)] for h in Cld_T[j]) <= \
                    pulp.lpSum(val[atom] * delta_alpha_T[(i, j, atom)] for atom in Lambda_nc), f"milp-(62)-T-{i}-{j}"
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            MILP += beta_F[(i, j)] + pulp.lpSum(beta_F[(i, h)] for h in Cld_F[j]) <= \
                    pulp.lpSum(val[atom] * delta_alpha_F[(i, j, atom)] for atom in Lambda_nc), f"milp-(62)-F-{i}-{j}"

    # -------- Constraint (63) --------
    for atom in Lambda_co:
        MILP += pulp.lpSum(delta_alpha_C[(i, 0, atom)] for i in range(1, t_C + 1)) == na_C[atom], f"milp-(63)-1-{atom}"
        MILP += pulp.lpSum(delta_alpha_T[(i, 0, atom)] for i in range(1, t_T + 1)) == na_T[atom], f"milp-(63)-2-{atom}"
    for atom in Lambda_nc:
        MILP += pulp.lpSum(delta_alpha_F[(i, 0, atom)] for i in range(1, t_F + 1)) == na_in[atom], f"milp-(63)-3-{atom}"

    # -------- Constraint (64) --------
    for p in range(1, rho + 1):
        for atom in Lambda_nc:
            MILP += pulp.lpSum(delta_alpha_C[(i, j, atom)] for i in range(1, t_C + 1) 
                    for j in Dsn_C[delta_i[i]][p]) == na_Cp[(p, atom)], f"milp-(64)-C-{p}-{atom}"
            MILP += pulp.lpSum(delta_alpha_T[(i, j, atom)] for i in range(1, t_T + 1) 
                    for j in Dsn_T[p]) == na_Tp[(p, atom)], f"milp-(64)-T-{p}-{atom}"
            MILP += pulp.lpSum(delta_alpha_F[(i, j, atom)] for i in range(1, t_F + 1) 
                    for j in Dsn_F[p]) == na_Fp[(p, atom)], f"milp-(64)-F-{p}-{atom}"

    # -------- Constraint (65) --------
    for atom in Lambda_co:
        MILP += na_C[atom] + na_T[atom] == na_co[atom], f"milp-(65)-1-{atom}"
        MILP += pulp.lpSum(na_Cp[(p, atom)] + na_Tp[(p, atom)] + na_Fp[(p, atom)] 
                for p in range(1, rho + 1)) == na_ex[atom], f"milp-(65)-2-{atom}"
        MILP += na_in[atom] + na_ex[atom] == na_nc[atom], f"milp-(65)-3-{atom}"
    for atom in Lambda:
        if atom in Lambda_co and atom in Lambda_nc:
            MILP += na_co[atom] + na_nc[atom] == na[atom], f"milp-(65)-4-{atom}"
        elif atom in Lambda_co:
            MILP += na_co[atom] == na[atom], f"milp-(65)-4-{atom}"
        elif atom in Lambda_nc:
            MILP += na_nc[atom] == na[atom], f"milp-(65)-4-{atom}"

    # -------- Constraint (66) --------
    MILP += pulp.lpSum(mass[atom] * na[atom] for atom in Lambda) == MASS, f"milp-(66)"

    # -------- Constraint (67) --------
    MILP += pulp.lpSum(val[atom] * na[atom] for atom in Lambda) - \
            2 * pulp.lpSum(m * (bd_co[m] + bd_in[m] + bd_ex[m]) for m in range(1, MAX_BOND + 1)) == \
            n_H, f"milp-(67)"

    # -------- Constraint (68) --------
    # for i in V_C_star:
    #     MILP += alpha_C[(i, 0)] == Code_Lambda[alpha_star[i]], f"milp-(68)-{i}"
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_alpha_C[(i, 0, a)] for a in Lambda_star[i]) == 1, f"milp-(68)-{i}"

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
    delta_chi_T,
    delta_beta_C,
    delta_beta_T,
    delta_beta_plus,
    delta_beta_minus,
    bd_T
    # Integer Variables
):
    # -------- Constraint (69) --------
    for i in (I_equal_one + I_zero_one):
        MILP += delta_beta_C[(i, 2)] >= bd2_LB[E_C[i]], f"milp-(69)-{i}-2-1"
        MILP += delta_beta_C[(i, 2)] <= bd2_UB[E_C[i]], f"milp-(69)-{i}-2-2"
        MILP += delta_beta_C[(i, 3)] >= bd3_LB[E_C[i]], f"milp-(69)-{i}-3-1"
        MILP += delta_beta_C[(i, 3)] <= bd3_UB[E_C[i]], f"milp-(69)-{i}-3-2"

    # -------- Constraint (70) --------
    for k in range(1, k_C + 1):
        for i in range(2, t_T + 1):
            for m in {2, 3}:
                MILP += bd_T[(k, i, m)] >= \
                    delta_beta_T[(i, m)] + delta_chi_T[(i, k)] - 1, f"milp-(70)-{k}-{i}-{m}"

    # -------- Constraint (71) --------
    for m in {2, 3}:
        MILP += pulp.lpSum(delta_beta_T[(j, m)] for j in range(2, t_T + 1)) >= \
                pulp.lpSum(bd_T[(k, i, m)] for k in range(1, k_C + 1) 
                                for i in range(2, t_T + 1)), f"milp-(71)-{m}"

    # -------- Constraint (72) --------
    for k in range(1, k_C + 1):
        MILP += pulp.lpSum(bd_T[(k, i, 2)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 2)] + delta_beta_minus[(k, 2)] >= bd2_LB[E_C[k]], f"milp-(72)-{k}-2-1"
        MILP += pulp.lpSum(bd_T[(k, i, 2)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 2)] + delta_beta_minus[(k, 2)] <= bd2_UB[E_C[k]], f"milp-(72)-{k}-2-2"
        MILP += pulp.lpSum(bd_T[(k, i, 3)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 3)] + delta_beta_minus[(k, 3)] >= bd3_LB[E_C[k]], f"milp-(72)-{k}-3-1"
        MILP += pulp.lpSum(bd_T[(k, i, 3)] for i in range(2, t_T + 1)) + \
                delta_beta_plus[(k, 3)] + delta_beta_minus[(k, 3)] <= bd3_UB[E_C[k]], f"milp-(72)-{k}-3-2"

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
    Dsn_C,
    Dsn_T,
    Dsn_F,
    rho,
    m_UB,
    Gamma_co_ac,
    Gamma_in_ac,
    Gamma_ex_ac,
    Gamma_tilde_ac_C,
    Gamma_tilde_ac_T,
    Gamma_tilde_ac_F,
    Gamma_tilde_ac_CT,
    Gamma_tilde_ac_TC,
    Gamma_tilde_ac_CF,
    Gamma_tilde_ac_TF,
    Gamma_tilde_ac_C_ex,
    Gamma_tilde_ac_T_ex,
    Gamma_tilde_ac_F_ex,
    ac_LB_co,
    ac_UB_co,
    ac_LB_in,
    ac_UB_in,
    ac_LB_ex,
    ac_UB_ex,
    Lambda_co,
    Lambda_nc,
    MAX_CODE_co,
    MAX_CODE_nc
):
    # Define delta_alpha_CT, delta_alpha_TC, delta_alpha_CF, delta_alpha_TF
    delta_alpha_CT = {(k, a): pulp.LpVariable(f"delta_alpha_CT({k},{a})", 0, 1, cat=pulp.LpBinary)
                    for k in range(1, k_C + 1) for a in Lambda_co}
    delta_alpha_TC = {(k, a): pulp.LpVariable(f"delta_alpha_TC({k},{a})", 0, 1, cat=pulp.LpBinary)
                    for k in range(1, k_C + 1) for a in Lambda_co}
    delta_alpha_CF = {(k, a): pulp.LpVariable(f"delta_alpha_CF({k},{a})", 0, 1, cat=pulp.LpBinary)
                    for k in range(1, t_C_tilde + 1) for a in Lambda_nc}
    delta_alpha_TF = {(k, a): pulp.LpVariable(f"delta_alpha_TF({k},{a})", 0, 1, cat=pulp.LpBinary)
                    for k in range(t_C_tilde + 1, c_F + 1) for a in Lambda_nc}
    # There is no \Lambda_in in this draft, maybe it should be \Lamnbda_nc, by Zhu

    # The variables ac_C, ac_T, ac_F, ac_Cp, ac_Tp, ac_Fp, ac_CT, ac_TC, ac_CF, ac_TF are defined twice in the draft.

    # Define ac_C, ac_T, ac_F
    # ac_C = {nu: pulp.LpVariable(f"ac_C({nu})", 0, m_C, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_C}
    # ac_T = {nu: pulp.LpVariable(f"ac_T({nu})", 0, t_T, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_T}
    # ac_F = {nu: pulp.LpVariable(f"ac_F({nu})", 0, t_F, cat=pulp.LpInteger) for nu in Gamma_tilde_ac_F}
    ac_C = {nu: pulp.LpVariable(f"ac_C({nu})", 0, m_C, cat=pulp.LpInteger) for nu in Gamma_co_ac}
    ac_T = {nu: pulp.LpVariable(f"ac_T({nu})", 0, t_T, cat=pulp.LpInteger) for nu in Gamma_co_ac}
    ac_F = {nu: pulp.LpVariable(f"ac_F({nu})", 0, t_F, cat=pulp.LpInteger) for nu in Gamma_in_ac}

    # Define ac_Cp, ac_Tp, ac_Fp
    # Here [1, n_X] should be [0, n_X], by Zhu
    # ac_Cp = {(p, nu): pulp.LpVariable(f"ac_C({p})({nu})", 0, n_C[2], cat=pulp.LpInteger)
    #         for p in range(1, rho + 1) for nu in Gamma_tilde_ac_C}
    # ac_Tp = {(p, nu): pulp.LpVariable(f"ac_T({p})({nu})", 0, n_T, cat=pulp.LpInteger)
    #         for p in range(1, rho + 1) for nu in Gamma_tilde_ac_T}
    # ac_Fp = {(p, nu): pulp.LpVariable(f"ac_F({p})({nu})", 0, n_F, cat=pulp.LpInteger)
    #         for p in range(1, rho + 1) for nu in Gamma_tilde_ac_F}
    ac_Cp = {(p, nu): pulp.LpVariable(f"ac_C({p})({nu})", 0, n_C[2], cat=pulp.LpInteger)
             for p in range(1, rho + 1) for nu in Gamma_ex_ac}
    ac_Tp = {(p, nu): pulp.LpVariable(f"ac_T({p})({nu})", 0, n_T, cat=pulp.LpInteger)
             for p in range(1, rho + 1) for nu in Gamma_ex_ac}
    ac_Fp = {(p, nu): pulp.LpVariable(f"ac_F({p})({nu})", 0, n_F, cat=pulp.LpInteger)
             for p in range(1, rho + 1) for nu in Gamma_ex_ac}

    # Define ac_CT, ac_TC
    # ac_CT = {nu: pulp.LpVariable(f"ac_CT({nu})", 0, min(k_C, t_T), cat=pulp.LpInteger)
    #         for nu in Gamma_tilde_ac_CT}
    # ac_TC = {nu: pulp.LpVariable(f"ac_TC({nu})", 0, min(k_C, t_T), cat=pulp.LpInteger)
    #         for nu in Gamma_tilde_ac_TC} # Gamma_tilde_ac_TC instead of Gamma_tilde_ac_CT
    ac_CT = {nu: pulp.LpVariable(f"ac_CT({nu})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for nu in Gamma_co_ac}
    ac_TC = {nu: pulp.LpVariable(f"ac_TC({nu})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for nu in Gamma_co_ac}

    # Define ac_CF, ac_TF
    # ac_CF = {nu: pulp.LpVariable(f"ac_CF({nu})", 0, t_C_tilde, cat=pulp.LpInteger)
    #         for nu in Gamma_tilde_ac_CF}
    # ac_TF = {nu: pulp.LpVariable(f"ac_TF({nu})", 0, t_T, cat=pulp.LpInteger)
    #         for nu in Gamma_tilde_ac_TF}
    ac_CF = {nu: pulp.LpVariable(f"ac_CF({nu})", 0, t_C_tilde, cat=pulp.LpInteger)
            for nu in Gamma_in_ac}
    ac_TF = {nu: pulp.LpVariable(f"ac_TF({nu})", 0, t_T, cat=pulp.LpInteger)
            for nu in Gamma_in_ac}

    # Define delta_ac_C
    delta_ac_C = {(i, nu): pulp.LpVariable(f"delta_ac_C({i},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(k_C_tilde + 1, m_C + 1) for nu in Gamma_tilde_ac_C}

    # Define delta_ac_T
    delta_ac_T = {(i, nu): pulp.LpVariable(f"delta_ac_T({i},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_T + 1) for nu in Gamma_tilde_ac_T}
    
    # Define delta_ac_F
    delta_ac_F = {(i, nu): pulp.LpVariable(f"delta_ac_F({i},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_F + 1) for nu in Gamma_tilde_ac_F}

    # Define delta_ac_Cp, delta_ac_Tp, delta_ac_Fp
    delta_ac_Cp = {(i, j, nu): pulp.LpVariable(f"delta_ac_C({i},{j},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_C + 1) 
            for j in range(1, n_C[delta_i[i]] + 1) for nu in Gamma_tilde_ac_C_ex}
    delta_ac_Tp = {(i, j, nu): pulp.LpVariable(f"delta_ac_T({i},{j},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_T + 1) 
            for j in range(1, n_T + 1) for nu in Gamma_tilde_ac_T_ex}
    delta_ac_Fp = {(i, j, nu): pulp.LpVariable(f"delta_ac_F({i},{j},{nu})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_F + 1) 
            for j in range(1, n_F + 1) for nu in Gamma_tilde_ac_F_ex}

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
    alpha_CT = {k: pulp.LpVariable(f"alpha_CT({k})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                for k in range(1, k_C + 1)}
    alpha_TC = {k: pulp.LpVariable(f"alpha_TC({k})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                for k in range(1, k_C + 1)}
    
    # Define alpha_CF
    alpha_CF = {c: pulp.LpVariable(f"alpha_CF({c})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                for c in range(1, t_C_tilde + 1)}

    # Define alpha_TF
    alpha_TF = {i: pulp.LpVariable(f"alpha_TF({i})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                for i in range(1, t_T + 1)}
    
    # Define Delta_ac_C_plus, Delta_ac_C_minus
    Delta_ac_C_plus = {i: pulp.LpVariable(f"Delta_ac_C_plus({i})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for i in range(k_C_tilde + 1, m_C + 1)}
    Delta_ac_C_minus = {i: pulp.LpVariable(f"Delta_ac_C_minus({i})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for i in range(k_C_tilde + 1, m_C + 1)}
    
    # Define Delta_ac_T_plus, Delta_ac_T_minus
    Delta_ac_T_plus = {i: pulp.LpVariable(f"Delta_ac_T_plus({i})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for i in range(2, t_T + 1)}
    Delta_ac_T_minus = {i: pulp.LpVariable(f"Delta_ac_T_minus({i})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for i in range(2, t_T + 1)}

    # Define Delta_ac_F_plus, Delta_ac_F_minus
    Delta_ac_F_plus = {i: pulp.LpVariable(f"Delta_ac_F_plus({i})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(2, t_F + 1)}
    Delta_ac_F_minus = {i: pulp.LpVariable(f"Delta_ac_F_minus({i})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(2, t_F + 1)}

    # Define Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus
    Delta_ac_C_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_C_plus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    Delta_ac_C_plus.update(Delta_ac_C_plus_temp)
    Delta_ac_C_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_C_minus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    Delta_ac_C_minus.update(Delta_ac_C_minus_temp)
    Delta_ac_T_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_T_plus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    Delta_ac_T_plus.update(Delta_ac_T_plus_temp)
    Delta_ac_T_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_T_minus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    Delta_ac_T_minus.update(Delta_ac_T_minus_temp)
    Delta_ac_F_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_F_plus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    Delta_ac_F_plus.update(Delta_ac_F_plus_temp)
    Delta_ac_F_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ac_F_minus({i},{j})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    Delta_ac_F_minus.update(Delta_ac_F_minus_temp)
    
    # Define Delta_ac_CT_plus, Delta_ac_CT_minus
    Delta_ac_CT_plus = {k: pulp.LpVariable(f"Delta_ac_CT_plus({k})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    Delta_ac_CT_minus = {k: pulp.LpVariable(f"Delta_ac_CT_minus({k})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}

    # Define Delta_ac_TC_plus, Delta_ac_TC_minus
    Delta_ac_TC_plus = {k: pulp.LpVariable(f"Delta_ac_TC_plus({k})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    Delta_ac_TC_minus = {k: pulp.LpVariable(f"Delta_ac_TC_minus({k})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for k in range(1, k_C + 1)}
    
    # Define Delta_ac_CF_plus, Delta_ac_CF_minus
    Delta_ac_CF_plus = {c: pulp.LpVariable(f"Delta_ac_CF_plus({c})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for c in range(1, t_C_tilde + 1)}
    Delta_ac_CF_minus = {c: pulp.LpVariable(f"Delta_ac_CF_minus({c})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for c in range(1, t_C_tilde + 1)}

    # Define Delta_ac_TF_plus, Delta_ac_TF_minus
    Delta_ac_TF_plus = {i: pulp.LpVariable(f"Delta_ac_TF_plus({i})", 0, MAX_CODE_co, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1)}
    Delta_ac_TF_minus = {i: pulp.LpVariable(f"Delta_ac_TF_minus({i})", 0, MAX_CODE_nc, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1)}

    # Define ac_co, ac_in, ac_ex
    ac_co = {nu: pulp.LpVariable(f"ac_co{nu}", ac_LB_co[nu], ac_UB_co[nu], cat=pulp.LpInteger)
            for nu in Gamma_co_ac}
    ac_in = {nu: pulp.LpVariable(f"ac_in{nu}", ac_LB_in[nu], ac_UB_in[nu], cat=pulp.LpInteger)
            for nu in Gamma_in_ac}
    ac_ex = {nu: pulp.LpVariable(f"ac_ex{nu}", ac_LB_ex[nu], ac_UB_ex[nu], cat=pulp.LpInteger)
            for nu in Gamma_ex_ac}

    # Define ac_nc
    ac_nc = {nu: pulp.LpVariable(f"ac_nc{nu}", 0, m_UB, cat=pulp.LpBinary) for nu in (set(Gamma_in_ac) | set(Gamma_ex_ac))}
    

    return ac_C, ac_T, ac_F, ac_Cp, ac_Tp, ac_Fp, ac_CT, ac_TC, ac_CF, ac_TF, \
        delta_ac_C, delta_ac_T, delta_ac_F, delta_ac_Cp, delta_ac_Tp, delta_ac_Fp, \
        ac_co, ac_in, ac_ex, ac_nc, \
        delta_alpha_CT, delta_alpha_TC, delta_alpha_CF, delta_alpha_TF, \
        alpha_CT, alpha_TC, alpha_CF, alpha_TF, \
        Delta_ac_C_plus, Delta_ac_C_minus, Delta_ac_T_plus, Delta_ac_T_minus, Delta_ac_F_plus, Delta_ac_F_minus, \
        Delta_ac_CT_plus, Delta_ac_CT_minus, Delta_ac_TC_plus, Delta_ac_TC_minus, \
        Delta_ac_CF_plus, Delta_ac_CF_minus, Delta_ac_TF_plus, Delta_ac_TF_minus, \
        delta_ac_CT, delta_ac_TC, delta_ac_CF, delta_ac_TF


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
    Dsn_C,
    Dsn_T,
    Dsn_F,
    prt_C,
    prt_T,
    prt_F,
    Gamma_co_ac,
    Gamma_in_ac,
    Gamma_ex_ac,
    Gamma_co_ac_less,
    Gamma_co_ac_equal,
    Gamma_tilde_ac_C,
    Gamma_tilde_ac_T,
    Gamma_tilde_ac_F,
    Gamma_tilde_ac_CT,
    Gamma_tilde_ac_TC,
    Gamma_tilde_ac_CF,
    Gamma_tilde_ac_TF,
    Gamma_tilde_ac_C_ex,
    Gamma_tilde_ac_T_ex,
    Gamma_tilde_ac_F_ex,
    rho,
    delta_i,
    Lambda_co,
    Lambda_nc,
    Code_Lambda_co,
    Code_Lambda_nc,
    MAX_CODE_co,
    MAX_CODE_nc,
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
    delta_ac_C,
    delta_ac_T,
    delta_ac_F,
    delta_ac_Cp,
    delta_ac_Tp,
    delta_ac_Fp,
    delta_ac_CT,
    delta_ac_TC,
    delta_ac_CF,
    delta_ac_TF,
    e_C,
    e_T,
    e_F,
    v_C,
    v_T,
    v_F,
    delta_alpha_CT,
    delta_alpha_TC,
    delta_alpha_CF,
    delta_alpha_TF,
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
    ac_Cp,
    ac_Tp,
    ac_Fp,
    ac_CT,
    ac_TC,
    ac_CF,
    ac_TF,
    ac_co,
    ac_in,
    ac_ex,
    ac_nc
):
    # -------- Constraint (73) --------
    for nu in Gamma_co_ac:
        if nu not in Gamma_tilde_ac_C:
            MILP += ac_C[nu] == 0, f"milp-(73)-1-{nu}"
        if nu not in Gamma_tilde_ac_T:
            MILP += ac_T[nu] == 0, f"milp-(73)-2-{nu}"
        if nu not in Gamma_tilde_ac_CT:
            MILP += ac_CT[nu] == 0, f"milp-(73)-5-{nu}"
        if nu not in Gamma_tilde_ac_TC:
            MILP += ac_TC[nu] == 0, f"milp-(73)-6-{nu}"
    for nu in Gamma_in_ac:
        if nu not in Gamma_tilde_ac_F:
            MILP += ac_F[nu] == 0, f"milp-(73)-3-{nu}"
        if nu not in Gamma_tilde_ac_CF:
            MILP += ac_CF[nu] == 0, f"milp-(73)-7-{nu}"
        if nu not in Gamma_tilde_ac_TF:
            MILP += ac_TF[nu] == 0, f"milp-(73)-8-{nu}"
    for nu in Gamma_ex_ac:
        if nu not in Gamma_tilde_ac_C_ex:
            for p in range(1, rho + 1):
                MILP += ac_Cp[(p, nu)] == 0, f"milp-(73)-4-C-{p}-{nu}"
        if nu not in Gamma_tilde_ac_T_ex:
            for p in range(1, rho + 1):
                MILP += ac_Tp[(p, nu)] == 0, f"milp-(73)-4-T-{p}-{nu}"
        if nu not in Gamma_tilde_ac_F_ex:
            for p in range(1, rho + 1):
                MILP += ac_Fp[(p, nu)] == 0, f"milp-(73)-4-F-{p}-{nu}"

    # -------- Constraint (74) --------
    for m in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum(ac_C[(a, b, m1)] for (a, b, m1) in Gamma_co_ac if m1 == m) == \
                pulp.lpSum(delta_beta_C[(i, m)] for i in range(k_C_tilde + 1, m_C + 1)), f"milp-(74)-1-{m}"
        MILP += pulp.lpSum(ac_T[(a, b, m1)] for (a, b, m1) in Gamma_co_ac if m1 == m) == \
                pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)), f"milp-(74)-2-{m}"
        MILP += pulp.lpSum(ac_F[(a, b, m1)] for (a, b, m1) in Gamma_in_ac if m1 == m) == \
                pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)), f"milp-(74)-3-{m}"
        for p in range(1, rho + 1):
            MILP += pulp.lpSum(ac_Cp[(p, (a, b, m1))] for (a, b, m1) in Gamma_ex_ac if m1 == m) == \
                    pulp.lpSum(delta_beta_C[(i, j, m)] for i in range(1, t_C + 1) for j in Dsn_C[delta_i[i]][p]), f"milp-(74)-4-C-{m}-{p}"
            MILP += pulp.lpSum(ac_Tp[(p, (a, b, m1))] for (a, b, m1) in Gamma_ex_ac if m1 == m) == \
                    pulp.lpSum(delta_beta_T[(i, j, m)] for i in range(1, t_T + 1) for j in Dsn_T[p]), f"milp-(74)-4-T-{m}-{p}"
            MILP += pulp.lpSum(ac_Fp[(p, (a, b, m1))] for (a, b, m1) in Gamma_ex_ac if m1 == m) == \
                    pulp.lpSum(delta_beta_F[(i, j, m)] for i in range(1, t_F + 1) for j in Dsn_F[p]), f"milp-(74)-4-F-{m}-{p}"
        MILP += pulp.lpSum(ac_CT[(a, b, m1)] for (a, b, m1) in Gamma_co_ac if m1 == m) == \
                pulp.lpSum(delta_beta_plus[(k, m)] for k in range(1, k_C + 1)), f"milp-(74)-5-{m}"
        MILP += pulp.lpSum(ac_TC[(a, b, m1)] for (a, b, m1) in Gamma_co_ac if m1 == m) == \
                pulp.lpSum(delta_beta_minus[(k, m)] for k in range(1, k_C + 1)), f"milp-(74)-6-{m}"
        MILP += pulp.lpSum(ac_CF[(a, b, m1)] for (a, b, m1) in Gamma_in_ac if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)), f"milp-(74)-7-{m}"
        MILP += pulp.lpSum(ac_TF[(a, b, m1)] for (a, b, m1) in Gamma_in_ac if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-(74)-8-{m}"

    # -------- Constraint (75) --------
    for i in range(k_C_tilde + 1, m_C + 1):
        MILP += pulp.lpSum(m * delta_ac_C[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_C) == beta_C[i], f"milp-(75)-1-{i}"
        MILP += Delta_ac_C_plus[i] + pulp.lpSum(Code_Lambda_co[a] * delta_ac_C[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_C) == alpha_C[(tail_C[i], 0)], f"milp-(75)-2-{i}"
        MILP += Delta_ac_C_minus[i] + pulp.lpSum(Code_Lambda_co[b] * delta_ac_C[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_C) == alpha_C[(head_C[i], 0)], f"milp-(75)-3-{i}"
        MILP += Delta_ac_C_plus[i] + Delta_ac_C_minus[i] <= 2 * MAX_CODE_co * (1 - e_C[i]), f"milp-(75)-4-{i}"
    for nu in Gamma_tilde_ac_C:
        MILP += pulp.lpSum(delta_ac_C[(i, nu)] 
                for i in range(k_C_tilde + 1, m_C + 1)) == ac_C[nu], f"milp-(75)-5-{nu}"

    # -------- Constraint (76) --------
    for i in range(2, t_T + 1):
        MILP += pulp.lpSum(m * delta_ac_T[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_T) == beta_T[i], f"milp-(76)-1-{i}"
        MILP += Delta_ac_T_plus[i] + pulp.lpSum(Code_Lambda_co[a] * delta_ac_T[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_T) == alpha_T[(i - 1, 0)], f"milp-(76)-2-{i}"
        MILP += Delta_ac_T_minus[i] + pulp.lpSum(Code_Lambda_co[b] * delta_ac_T[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_T) == alpha_T[(i, 0)], f"milp-(76)-3-{i}"
        MILP += Delta_ac_T_plus[i] + Delta_ac_T_minus[i] <= 2 * MAX_CODE_co * (1 - e_T[i]), f"milp-(76)-4-{i}"
    for nu in Gamma_tilde_ac_T:
        MILP += pulp.lpSum(delta_ac_T[(i, nu)] 
                for i in range(2, t_T + 1)) == ac_T[nu], f"milp-(76)-5-{nu}"

    # -------- Constraint (77) --------
    for i in range(2, t_F + 1):
        MILP += pulp.lpSum(m * delta_ac_F[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_F) == beta_F[i], f"milp-(77)-1-{i}"
        MILP += Delta_ac_F_plus[i] + pulp.lpSum(Code_Lambda_nc[a] * delta_ac_F[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[(i - 1, 0)], f"milp-(77)-2-{i}"
        MILP += Delta_ac_F_minus[i] + pulp.lpSum(Code_Lambda_nc[b] * delta_ac_F[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_F) == alpha_F[(i, 0)], f"milp-(77)-3-{i}"
        MILP += Delta_ac_F_plus[i] + Delta_ac_F_minus[i] <= 2 * MAX_CODE_nc * (1 - e_F[i]), f"milp-(77)-4-{i}"
    for nu in Gamma_tilde_ac_F:
        MILP += pulp.lpSum(delta_ac_F[(i, nu)] 
                for i in range(2, t_F + 1)) == ac_F[nu], f"milp-(77)-5-{nu}"

    # # -------- Constraint (78) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += pulp.lpSum(m * delta_ac_Cp[(i, j, (a, b, m))] 
                for (a, b, m) in Gamma_tilde_ac_C_ex) == beta_C[(i, j)], f"milp-(78)-C-1-{i}-{j}"
            MILP += Delta_ac_C_plus[(i, j)] + pulp.lpSum(Code_Lambda_nc[a] * delta_ac_Cp[(i, j, (a, b, m))]
                        for (a, b, m) in Gamma_tilde_ac_C_ex) == alpha_C[(i, prt_C[delta_i[i]][j])], f"milp-(78)-C-2-{i}-{j}"
            MILP += Delta_ac_C_minus[(i, j)] + pulp.lpSum(Code_Lambda_nc[b] * delta_ac_Cp[(i, j, (a, b, m))]
                        for (a, b, m) in Gamma_tilde_ac_C_ex) == alpha_C[(i, j)], f"milp-(77)-C-3-{i}-{j}"
            MILP += Delta_ac_C_plus[(i, j)] + Delta_ac_C_minus[(i, j)] <= 2 * MAX_CODE_nc * (1 - v_C[(i, j)]), f"milp-(77)-C-4-{i}-{j}"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += pulp.lpSum(m * delta_ac_Tp[(i, j, (a, b, m))] 
                for (a, b, m) in Gamma_tilde_ac_T_ex) == beta_T[(i, j)], f"milp-(78)-T-1-{i}-{j}"
            MILP += Delta_ac_T_plus[(i, j)] + pulp.lpSum(Code_Lambda_nc[a] * delta_ac_Tp[(i, j, (a, b, m))]
                        for (a, b, m) in Gamma_tilde_ac_T_ex) == alpha_T[(i, prt_T[j])], f"milp-(78)-T-2-{i}-{j}"
            MILP += Delta_ac_T_minus[(i, j)] + pulp.lpSum(Code_Lambda_nc[b] * delta_ac_Tp[(i, j, (a, b, m))]
                        for (a, b, m) in Gamma_tilde_ac_T_ex) == alpha_T[(i, j)], f"milp-(77)-T-3-{i}-{j}"
            MILP += Delta_ac_T_plus[(i, j)] + Delta_ac_T_minus[(i, j)] <= 2 * MAX_CODE_nc * (1 - v_T[(i, j)]), f"milp-(77)-T-4-{i}-{j}"
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            MILP += pulp.lpSum(m * delta_ac_Fp[(i, j, (a, b, m))] 
                for (a, b, m) in Gamma_tilde_ac_F_ex) == beta_F[(i, j)], f"milp-(78)-F-1-{i}-{j}"
            MILP += Delta_ac_F_plus[(i, j)] + pulp.lpSum(Code_Lambda_nc[a] * delta_ac_Fp[(i, j, (a, b, m))]
                        for (a, b, m) in Gamma_tilde_ac_F_ex) == alpha_F[(i, prt_F[j])], f"milp-(78)-F-2-{i}-{j}"
            MILP += Delta_ac_F_minus[(i, j)] + pulp.lpSum(Code_Lambda_nc[b] * delta_ac_Fp[(i, j, (a, b, m))]
                        for (a, b, m) in Gamma_tilde_ac_F_ex) == alpha_F[(i, j)], f"milp-(77)-F-3-{i}-{j}"
            MILP += Delta_ac_F_plus[(i, j)] + Delta_ac_F_minus[(i, j)] <= 2 * MAX_CODE_nc * (1 - v_F[(i, j)]), f"milp-(77)-F-4-{i}-{j}"      
    for p in range(1, rho + 1):
        for nu in Gamma_tilde_ac_C_ex:
            MILP += pulp.lpSum(delta_ac_Cp[(i, j, nu)] 
                    for i in range(1, t_C + 1) for j in Dsn_C[delta_i[i]][p]) == ac_Cp[(p, nu)], f"milp-(78)-C-5-{p}-{nu}"
        for nu in Gamma_tilde_ac_T_ex:
            MILP += pulp.lpSum(delta_ac_Tp[(i, j, nu)] 
                    for i in range(1, t_T + 1) for j in Dsn_T[p]) == ac_Tp[(p, nu)], f"milp-(78)-T-5-{p}-{nu}"
        for nu in Gamma_tilde_ac_F_ex:
            MILP += pulp.lpSum(delta_ac_Fp[(i, j, nu)] 
                    for i in range(1, t_F + 1) for j in Dsn_F[p]) == ac_Fp[(p, nu)], f"milp-(78)-F-5-{p}-{nu}"
        
    # -------- Constraint (79) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += alpha_T[(i, 0)] + MAX_CODE_co * (1 - delta_chi_T[(i, k)] + e_T[i]) >= \
                    alpha_CT[k], f"milp-(79)-1-{k}-{i}"
            MILP += alpha_CT[k] >= alpha_T[(i, 0)] - \
                    MAX_CODE_co * (1 - delta_chi_T[(i, k)] + e_T[i]), f"milp-(79)-2-{k}-{i}"
        MILP += pulp.lpSum(m * delta_ac_CT[(k, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_CT) == beta_plus[k], f"milp-(79)-3-{k}"
        MILP += Delta_ac_CT_plus[k] + pulp.lpSum(Code_Lambda_co[a] * delta_ac_CT[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CT) == alpha_C[(tail_C[k], 0)], f"milp-(79)-4-{k}"
        MILP += Delta_ac_CT_minus[k] + pulp.lpSum(Code_Lambda_co[b] * delta_ac_CT[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CT) == alpha_CT[k], f"milp-(79)-5-{k}"
        MILP += Delta_ac_CT_plus[k] + Delta_ac_CT_minus[k] <= 2 * MAX_CODE_co * (1 - delta_chi_T[k]), f"milp-(79)-6-{k}"
    for nu in Gamma_tilde_ac_CT:
        MILP += pulp.lpSum(delta_ac_CT[(k, nu)] for k in range(1, k_C + 1)) == ac_CT[nu], f"milp-(79)-7-{nu}"

    # -------- Constraint (80) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += alpha_T[(i, 0)] + MAX_CODE_co * (1 - delta_chi_T[(i, k)] + e_T[i + 1]) >= \
                    alpha_TC[k], f"milp-(80)-1-{k}-{i}"
            MILP += alpha_TC[k] >= alpha_T[(i, 0)] - \
                    MAX_CODE_co * (1 - delta_chi_T[(i, k)] + e_T[i + 1]), f"milp-(80)-2-{k}-{i}"
        MILP += pulp.lpSum(m * delta_ac_TC[(k, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_TC) == beta_minus[k], f"milp-(80)-3-{k}"
        MILP += Delta_ac_TC_plus[k] + pulp.lpSum(Code_Lambda_co[a] * delta_ac_TC[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TC) == alpha_TC[k], f"milp-(80)-4-{k}"
        MILP += Delta_ac_TC_minus[k] + pulp.lpSum(Code_Lambda_co[b] * delta_ac_TC[(k, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TC) == alpha_C[(head_C[k], 0)], f"milp-(80)-5-{k}"
        MILP += Delta_ac_TC_plus[k] + Delta_ac_TC_minus[k] <= 2 * MAX_CODE_co * (1 - delta_chi_T[k]), f"milp-(80)-6-{k}"
    for nu in Gamma_tilde_ac_TC:
        MILP += pulp.lpSum(delta_ac_TC[(k, nu)] for k in range(1, k_C + 1)) == ac_TC[nu], f"milp-(80)-7-{nu}"

    # -------- Constraint (81) --------
    for c in range(1, t_C_tilde + 1):
        for i in range(1, t_F + 1):
            MILP += alpha_F[(i, 0)] + MAX_CODE_nc * (1 - delta_chi_F[(i, c)] + e_F[i]) >= \
                    alpha_CF[c], f"milp-(81)-1-{c}-{i}"
            MILP += alpha_CF[c] >= alpha_F[(i, 0)] - \
                    MAX_CODE_nc * (1 - delta_chi_F[(i, c)] + e_F[i]), f"milp-(81)-2-{c}-{i}"
        MILP += pulp.lpSum(m * delta_ac_CF[(c, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_CF) == beta_in[c], f"milp-(81)-3-{c}"
        MILP += Delta_ac_CF_plus[c] + pulp.lpSum(Code_Lambda_co[a] * delta_ac_CF[(c, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CF) == alpha_C[(head_C[c], 0)], f"milp-(81)-4-{c}"
        MILP += Delta_ac_CF_minus[c] + pulp.lpSum(Code_Lambda_nc[b] * delta_ac_CF[(c, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_CF) == alpha_CF[c], f"milp-(81)-5-{c}"
        MILP += Delta_ac_CF_plus[c] + Delta_ac_CF_minus[c] <= \
                2 * max(MAX_CODE_co, MAX_CODE_nc) * (1 - delta_chi_F[c]), f"milp-(81)-6-{c}"
    for nu in Gamma_tilde_ac_CF:
        MILP += pulp.lpSum(delta_ac_CF[(c, nu)] for c in range(1, t_C_tilde + 1)) == ac_CF[nu], f"milp-(81)-7-{nu}"

    # -------- Constraint (82) --------
    for i in range(1, t_T + 1):
        for j in range(1, t_F + 1):
            MILP += alpha_F[(j, 0)] + MAX_CODE_nc * (1 - delta_chi_F[(j, i + t_C_tilde)] + e_F[j]) >= \
                    alpha_TF[i], f"milp-(82)-1-{i}-{j}"
            MILP += alpha_TF[i] >= alpha_F[(j, 0)] - \
                    MAX_CODE_nc * (1 - delta_chi_F[(j, i + t_C_tilde)] + e_F[j]), f"milp-(82)-2-{i}-{j}"
        MILP += pulp.lpSum(m * delta_ac_TF[(i, (a, b, m))] 
            for (a, b, m) in Gamma_tilde_ac_TF) == beta_in[i + t_C_tilde], f"milp-(82)-3-{i}"
        MILP += Delta_ac_TF_plus[i] + pulp.lpSum(Code_Lambda_co[a] * delta_ac_TF[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TF) == alpha_T[(i, 0)], f"milp-(82)-4-{i}"
        MILP += Delta_ac_TF_minus[i] + pulp.lpSum(Code_Lambda_nc[b] * delta_ac_TF[(i, (a, b, m))]
                    for (a, b, m) in Gamma_tilde_ac_TF) == alpha_TF[i], f"milp-(82)-5-{i}"
        MILP += Delta_ac_TF_plus[i] + Delta_ac_TF_minus[i] <= \
                2 * max(MAX_CODE_co, MAX_CODE_nc) * (1 - delta_chi_F[i + t_C_tilde]), f"milp-(82)-6-{i}"
    for nu in Gamma_tilde_ac_TF:
        MILP += pulp.lpSum(delta_ac_TF[(i, nu)] for i in range(1, t_T + 1)) == ac_TF[nu], f"milp-(82)-7-{nu}"

    # -------- Constraint (83) --------
    for (a, b, m) in Gamma_co_ac_less:
        MILP += ac_C[(a, b, m)] + ac_C[(b, a, m)] + \
                ac_T[(a, b, m)] + ac_T[(b, a, m)] + \
                ac_CT[(a, b, m)] + ac_CT[(b, a, m)] + \
                ac_TC[(a, b, m)] + ac_TC[(b, a, m)] == ac_co[(a, b, m)], f"milp-(83)-1-{(a, b, m)}"
    for nu in Gamma_co_ac_equal:
        MILP += ac_C[nu] + ac_T[nu] + ac_CT[nu] + ac_TC[nu] == ac_co[nu], f"milp-(83)-2-{nu}"
    for nu in Gamma_in_ac:
        MILP += ac_F[nu] + ac_CF[nu] + ac_TF[nu] == ac_in[nu], f"milp-(83)-3-{nu}"
    for nu in Gamma_ex_ac:
        MILP += pulp.lpSum(ac_Cp[(p, nu)] + ac_Tp[(p, nu)] + ac_Fp[(p, nu)] 
                for p in range(1, rho + 1)) == ac_ex[nu], f"milp-(83)-4-{nu}"

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
    cs_LB,
    cs_UB,
    Lambda_dg_co,
    Lambda_dg_nc,
    epsilon_dg
):
    # Define ns_co
    ns_co = {mu: pulp.LpVariable(f"ns_co({mu})", 0, cs_UB, cat=pulp.LpInteger) for mu in Lambda_dg_co}

    # Define ns_nc
    ns_nc = {mu: pulp.LpVariable(f"ns_nc({mu})", 0, n_star - cs_LB, cat=pulp.LpInteger) for mu in Lambda_dg_nc}

    # Define delta_ns_C, delta_ns_T
    delta_ns_C = {(i, 0, mu): pulp.LpVariable(f"delta_ns_C({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_C + 1) for mu in (Lambda_dg_co + [epsilon_dg])} # no j is used here
    delta_ns_T = {(i, 0, mu): pulp.LpVariable(f"delta_ns_T({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_T + 1) for mu in (Lambda_dg_co + [epsilon_dg])} # no j is used here
    
    # Define delta_ns_F
    delta_ns_F = {(i, 0, mu): pulp.LpVariable(f"delta_ns_F({i},{0},{mu})", 0, 1, cat=pulp.LpBinary)
                    for i in range(1, t_F + 1) for mu in (Lambda_dg_nc + [epsilon_dg])}

    # Define delta_ns_C, delta_ns_T, delta_ns_F
    delta_ns_C_temp = {(i, j, mu): pulp.LpVariable(f"delta_ns_C({i},{j},{mu})", 0, 1, cat=pulp.LpBinary)
                        for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1) 
                        for mu in (Lambda_dg_nc + [epsilon_dg])}
    delta_ns_C.update(delta_ns_C_temp)
    delta_ns_T_temp = {(i, j, mu): pulp.LpVariable(f"delta_ns_T({i},{j},{mu})", 0, 1, cat=pulp.LpBinary)
                        for i in range(1, t_T + 1) for j in range(1, n_T + 1) 
                        for mu in (Lambda_dg_nc + [epsilon_dg])}
    delta_ns_T.update(delta_ns_T_temp)
    delta_ns_F_temp = {(i, j, mu): pulp.LpVariable(f"delta_ns_F({i},{j},{mu})", 0, 1, cat=pulp.LpBinary)
                        for i in range(1, t_F + 1) for j in range(1, n_F + 1) 
                        for mu in (Lambda_dg_nc + [epsilon_dg])}
    delta_ns_F.update(delta_ns_F_temp)
    
    return ns_co, ns_nc, delta_ns_C, delta_ns_T, delta_ns_F
    
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
    Lambda_dg_co,
    Lambda_dg_nc,
    Lambda_tilde_dg_C,
    Lambda_tilde_dg_T,
    Lambda_tilde_dg_F,
    Lambda_tilde_dg_C_nc,
    Lambda_tilde_dg_T_nc,
    Lambda_tilde_dg_F_nc,
    Code_Lambda,
    Code_Lambda_co,
    Code_Lambda_nc,
    epsilon_dg,
    delta_i,
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
    ns_co,
    ns_nc
):
    # -------- Constraint (84) --------
    for i in range(1, t_C + 1):
        MILP += pulp.lpSum(delta_ns_C[(i, 0, mu)] 
                for mu in (Lambda_tilde_dg_C + [epsilon_dg])) == 1, f"milp-(84)-C-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_co[a] * delta_ns_C[(i, 0, (a, d))] 
                for (a, d) in Lambda_tilde_dg_C) == alpha_C[(i, 0)], f"milp-(84)-C-2-{i}"
        MILP += pulp.lpSum(d * delta_ns_C[(i, 0, (a, d))]
                for (a, d) in Lambda_tilde_dg_C) == deg_C[(i, 0)], f"milp-(84)-C-3-{i}"
    for i in range(1, t_T + 1):
        MILP += pulp.lpSum(delta_ns_T[(i, 0, mu)] 
                for mu in (Lambda_tilde_dg_T + [epsilon_dg])) == 1, f"milp-(84)-T-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_co[a] * delta_ns_T[(i, 0, (a, d))] 
                for (a, d) in Lambda_tilde_dg_T) == alpha_T[(i, 0)], f"milp-(84)-T-2-{i}"
        MILP += pulp.lpSum(d * delta_ns_T[(i, 0, (a, d))]
                for (a, d) in Lambda_tilde_dg_T) == deg_T[(i, 0)], f"milp-(84)-T-3-{i}"
    
    # -------- Constraint (85) --------
    for i in range(1, t_F + 1):
        MILP += pulp.lpSum(delta_ns_F[(i, 0, mu)] 
                for mu in (Lambda_tilde_dg_F + [epsilon_dg])) == 1, f"milp-(85)-1-{i}"
        MILP += pulp.lpSum(Code_Lambda_nc[a] * delta_ns_F[(i, 0, (a, d))] 
                for (a, d) in Lambda_tilde_dg_F) == alpha_F[(i, 0)], f"milp-(85)-2-{i}"
        MILP += pulp.lpSum(d * delta_ns_F[(i, 0, (a, d))]
                for (a, d) in Lambda_tilde_dg_F) == deg_F[(i, 0)], f"milp-(85)-3-{i}"
    
    # -------- Constraint (86) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += pulp.lpSum(delta_ns_C[(i, j, mu)] 
                    for mu in (Lambda_tilde_dg_C_nc + [epsilon_dg])) == 1, f"milp-(86)-C-1-{i}-{j}"
            MILP += pulp.lpSum(Code_Lambda_nc[a] * delta_ns_C[(i, j, (a, d))] 
                    for (a, d) in Lambda_tilde_dg_C_nc) == alpha_C[(i, j)], f"milp-(86)-C-2-{i}-{j}"
            MILP += pulp.lpSum(d * delta_ns_C[(i, j, (a, d))]
                    for (a, d) in Lambda_tilde_dg_C_nc) == deg_C[(i, j)], f"milp-(86)-C-3-{i}-{j}"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += pulp.lpSum(delta_ns_T[(i, j, mu)] 
                    for mu in (Lambda_tilde_dg_T_nc + [epsilon_dg])) == 1, f"milp-(86)-T-1-{i}-{j}"
            MILP += pulp.lpSum(Code_Lambda_nc[a] * delta_ns_T[(i, j, (a, d))] 
                    for (a, d) in Lambda_tilde_dg_T_nc) == alpha_T[(i, j)], f"milp-(86)-T-2-{i}-{j}"
            MILP += pulp.lpSum(d * delta_ns_T[(i, j, (a, d))]
                    for (a, d) in Lambda_tilde_dg_T_nc) == deg_T[(i, j)], f"milp-(86)-T-3-{i}-{j}"
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            MILP += pulp.lpSum(delta_ns_F[(i, j, mu)] 
                    for mu in (Lambda_tilde_dg_F_nc + [epsilon_dg])) == 1, f"milp-(86)-F-1-{i}-{j}"
            MILP += pulp.lpSum(Code_Lambda_nc[a] * delta_ns_F[(i, j, (a, d))]
                    for (a, d) in Lambda_tilde_dg_F_nc) == alpha_F[(i, j)], f"milp-(86)-F-2-{i}-{j}"
            MILP += pulp.lpSum(d * delta_ns_F[(i, j, (a, d))]
                    for (a, d) in Lambda_tilde_dg_F_nc) == deg_F[(i, j)], f"milp-(86)-F-3-{i}-{j}"
    
    # -------- Constraint (87) --------
    for mu in Lambda_dg_co:
        MILP += pulp.lpSum(delta_ns_C[(i, 0, mu)] for i in range(1, t_C + 1)) + \
                pulp.lpSum(delta_ns_T[(i, 0, mu)] for i in range(1, t_T + 1)) == ns_co[mu], f"milp-(87)-1-{mu}"
    for mu in Lambda_dg_nc:
        MILP += pulp.lpSum(delta_ns_F[(i, 0, mu)] for i in range(1, t_F + 1)) + \
                pulp.lpSum(delta_ns_C[(i, j, mu)] for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)) + \
                pulp.lpSum(delta_ns_T[(i, j, mu)] for i in range(1, t_T + 1) for j in range(1, n_T + 1)) + \
                pulp.lpSum(delta_ns_F[(i, j, mu)] for i in range(1, t_F + 1) for j in range(1, n_F + 1)) == \
                ns_nc[mu], f"milp-(87)-2-{mu}"

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
    Dsn_C,
    Dsn_T,
    Dsn_F,
    rho,
    m_UB,
    Gamma_co,
    Gamma_in,
    Gamma_ex,
    Gamma_tilde_ec_C,
    Gamma_tilde_ec_T,
    Gamma_tilde_ec_F,
    Gamma_tilde_ec_CT,
    Gamma_tilde_ec_TC,
    Gamma_tilde_ec_CF,
    Gamma_tilde_ec_TF,
    Gamma_tilde_ec_C_ex,
    Gamma_tilde_ec_T_ex,
    Gamma_tilde_ec_F_ex,
    ec_LB_co,
    ec_UB_co,
    ec_LB_in,
    ec_UB_in,
    ec_LB_ex,
    ec_UB_ex
):
    # Define delta_dg_CT, delta_dg_TC, delta_dg_CF, delta_dg_TF
    delta_dg_CT = {(k, d): pulp.LpVariable(f"delta_dg_CT({k},{d})", 0, 1, cat=pulp.LpBinary)
                    for d in range(1, MAX_VAL + 1) for k in range(1, k_C + 1)}
    delta_dg_TC = {(k, d): pulp.LpVariable(f"delta_dg_TC({k},{d})", 0, 1, cat=pulp.LpBinary)
                    for d in range(1, MAX_VAL + 1) for k in range(1, k_C + 1)}
    delta_dg_CF = {(c, d): pulp.LpVariable(f"delta_dg_CF({c},{d})", 0, 1, cat=pulp.LpBinary)
                    for d in range(1, MAX_VAL + 1) for c in range(1, t_C_tilde + 1)}
    delta_dg_TF = {(c, d): pulp.LpVariable(f"delta_dg_TF({c},{d})", 0, 1, cat=pulp.LpBinary)
                    for d in range(1, MAX_VAL + 1) for c in range(t_C_tilde, c_F + 1)}    

    # The variables ec_C, ec_T, ec_F, ec_Cp, ec_Tp, ec_Fp, ec_CT, ec_TC, ec_CF, ec_TF are defined twice in the draft.

    # Define ec_C, ec_T, ec_F
    # ec_C = {gamma: pulp.LpVariable(f"ec_C({gamma})", 0, m_C, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_C}
    # ec_T = {gamma: pulp.LpVariable(f"ec_T({gamma})", 0, t_T, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_T}
    # ec_F = {gamma: pulp.LpVariable(f"ec_F({gamma})", 0, t_F, cat=pulp.LpInteger) for gamma in Gamma_tilde_ec_F}
    ec_C = {gamma: pulp.LpVariable(f"ec_C({gamma})", 0, m_C, cat=pulp.LpInteger) for gamma in Gamma_co}
    ec_T = {gamma: pulp.LpVariable(f"ec_T({gamma})", 0, t_T, cat=pulp.LpInteger) for gamma in Gamma_co}
    ec_F = {gamma: pulp.LpVariable(f"ec_F({gamma})", 0, t_F, cat=pulp.LpInteger) for gamma in Gamma_in}

    # Define ec_Cp, ec_Tp, ec_Fp
    # Here [1, n_X] should be [0, n_X], by Zhu
    # ec_Cp = {(p, gamma): pulp.LpVariable(f"ec_C({p})({gamma})", 0, n_C[2], cat=pulp.LpInteger)
    #         for p in range(1, rho + 1) for gamma in Gamma_tilde_ec_C}
    # ec_Tp = {(p, gamma): pulp.LpVariable(f"ec_T({p})({gamma})", 0, n_T, cat=pulp.LpInteger)
    #         for p in range(1, rho + 1) for gamma in Gamma_tilde_ec_T}
    # ec_Fp = {(p, gamma): pulp.LpVariable(f"ec_F({p})({gamma})", 0, n_F, cat=pulp.LpInteger)
    #         for p in range(1, rho + 1) for gamma in Gamma_tilde_ec_F}
    ec_Cp = {(p, gamma): pulp.LpVariable(f"ec_C({p})({gamma})", 0, n_C[2], cat=pulp.LpInteger)
             for p in range(1, rho + 1) for gamma in Gamma_ex}
    ec_Tp = {(p, gamma): pulp.LpVariable(f"ec_T({p})({gamma})", 0, n_T, cat=pulp.LpInteger)
             for p in range(1, rho + 1) for gamma in Gamma_ex}
    ec_Fp = {(p, gamma): pulp.LpVariable(f"ec_F({p})({gamma})", 0, n_F, cat=pulp.LpInteger)
             for p in range(1, rho + 1) for gamma in Gamma_ex}

    # Define ec_CT, ec_TC
    # ec_CT = {gamma: pulp.LpVariable(f"ec_CT({gamma})", 0, min(k_C, t_T), cat=pulp.LpInteger)
    #         for gamma in Gamma_tilde_ec_CT}
    # ec_TC = {gamma: pulp.LpVariable(f"ec_TC({gamma})", 0, min(k_C, t_T), cat=pulp.LpInteger)
    #         for gamma in Gamma_tilde_ec_TC} # Gamma_tilde_ec_TC instead of Gamma_tilde_ec_CT
    ec_CT = {gamma: pulp.LpVariable(f"ec_CT({gamma})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for gamma in Gamma_co}
    ec_TC = {gamma: pulp.LpVariable(f"ec_TC({gamma})", 0, min(k_C, t_T), cat=pulp.LpInteger)
            for gamma in Gamma_co}

    # Define ec_CF, ec_TF
    # ec_CF = {gamma: pulp.LpVariable(f"ec_CF({gamma})", 0, t_C_tilde, cat=pulp.LpInteger)
    #         for gamma in Gamma_tilde_ec_CF}
    # ec_TF = {gamma: pulp.LpVariable(f"ec_TF({gamma})", 0, t_T, cat=pulp.LpInteger)
    #         for gamma in Gamma_tilde_ec_TF}
    ec_CF = {gamma: pulp.LpVariable(f"ec_CF({gamma})", 0, t_C_tilde, cat=pulp.LpInteger)
            for gamma in Gamma_in}
    ec_TF = {gamma: pulp.LpVariable(f"ec_TF({gamma})", 0, t_T, cat=pulp.LpInteger)
            for gamma in Gamma_in}

    # Define delta_ec_C
    delta_ec_C = {(i, gamma): pulp.LpVariable(f"delta_ec_C({i},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(k_C_tilde + 1, m_C + 1) for gamma in Gamma_tilde_ec_C}

    # Define delta_ec_T
    delta_ec_T = {(i, gamma): pulp.LpVariable(f"delta_ec_T({i},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_T + 1) for gamma in Gamma_tilde_ec_T}
    
    # Define delta_ec_F
    delta_ec_F = {(i, gamma): pulp.LpVariable(f"delta_ec_F({i},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(2, t_F + 1) for gamma in Gamma_tilde_ec_F}

    # Define delta_ec_Cp, delta_ec_Tp, delta_ec_Fp
    delta_ec_Cp = {(i, j, gamma): pulp.LpVariable(f"delta_ec_C({i},{j},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_C + 1) 
            for j in range(1, n_C[delta_i[i]] + 1) for gamma in Gamma_tilde_ec_C_ex}
    delta_ec_Tp = {(i, j, gamma): pulp.LpVariable(f"delta_ec_T({i},{j},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_T + 1) 
            for j in range(1, n_T + 1) for gamma in Gamma_tilde_ec_T_ex}
    delta_ec_Fp = {(i, j, gamma): pulp.LpVariable(f"delta_ec_F({i},{j},{gamma})", 0, 1, cat=pulp.LpBinary)
            for i in range(1, t_F + 1) 
            for j in range(1, n_F + 1) for gamma in Gamma_tilde_ec_F_ex}

    # Define delta_ec_CT, delta_ex_TC
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

    # Define Delta_ec_C_plus, Delta_ec_C_minus, Delta_ec_T_plus, Delta_ec_T_minus, Delta_ec_F_plus, Delta_ec_F_minus
    Delta_ec_C_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_C_plus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    Delta_ec_C_plus.update(Delta_ec_C_plus_temp)
    Delta_ec_C_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_C_minus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_C + 1) for j in range(1, n_C[delta_i[i]] + 1)}
    Delta_ec_C_minus.update(Delta_ec_C_minus_temp)
    Delta_ec_T_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_T_plus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    Delta_ec_T_plus.update(Delta_ec_T_plus_temp)
    Delta_ec_T_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_T_minus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_T + 1) for j in range(1, n_T + 1)}
    Delta_ec_T_minus.update(Delta_ec_T_minus_temp)
    Delta_ec_F_plus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_F_plus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    Delta_ec_F_plus.update(Delta_ec_F_plus_temp)
    Delta_ec_F_minus_temp = {(i, j): pulp.LpVariable(f"Delta_ec_F_minus({i},{j})", 0, MAX_VAL, cat=pulp.LpInteger)
                    for i in range(1, t_F + 1) for j in range(1, n_F + 1)}
    Delta_ec_F_minus.update(Delta_ec_F_minus_temp)
    
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

    # Define ec_co, ec_in, ec_ex
    ec_co = {gamma: pulp.LpVariable(f"ec_co{gamma}", ec_LB_co[gamma], ec_UB_co[gamma], cat=pulp.LpInteger)
            for gamma in Gamma_co}
    ec_in = {gamma: pulp.LpVariable(f"ec_in{gamma}", ec_LB_in[gamma], ec_UB_in[gamma], cat=pulp.LpInteger)
            for gamma in Gamma_in}
    ec_ex = {gamma: pulp.LpVariable(f"ec_ex{gamma}", ec_LB_ex[gamma], ec_UB_ex[gamma], cat=pulp.LpInteger)
            for gamma in Gamma_ex}

    # Define ec_nc
    ec_nc = {gamma: pulp.LpVariable(f"ec_nc{gamma}", 0, m_UB, cat=pulp.LpBinary) for gamma in (set(Gamma_in) | set(Gamma_ex))}
    
    return ec_C, ec_T, ec_F, ec_Cp, ec_Tp, ec_Fp, ec_CT, ec_TC, ec_CF, ec_TF, \
        delta_ec_C, delta_ec_T, delta_ec_F, delta_ec_Cp, delta_ec_Tp, delta_ec_Fp, \
        delta_ec_CT, delta_ec_TC, delta_ec_CF, delta_ec_TF, ec_co, ec_in, ec_ex, ec_nc, \
        delta_dg_CT, delta_dg_TC, delta_dg_CF, delta_dg_TF, \
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
    Dsn_C,
    Dsn_T,
    Dsn_F,
    prt_C,
    prt_T,
    prt_F,
    Gamma_co,
    Gamma_in,
    Gamma_ex,
    Gamma_co_less,
    Gamma_co_equal,
    Gamma_tilde_ec_C,
    Gamma_tilde_ec_T,
    Gamma_tilde_ec_F,
    Gamma_tilde_ec_CT,
    Gamma_tilde_ec_TC,
    Gamma_tilde_ec_CF,
    Gamma_tilde_ec_TF,
    Gamma_tilde_ec_C_ex,
    Gamma_tilde_ec_T_ex,
    Gamma_tilde_ec_F_ex,
    Gamma_tilde_ac_C,
    Gamma_tilde_ac_T,
    Gamma_tilde_ac_F,
    Gamma_tilde_ac_CT,
    Gamma_tilde_ac_TC,
    Gamma_tilde_ac_CF,
    Gamma_tilde_ac_TF,
    Gamma_tilde_ac_C_ex,
    Gamma_tilde_ac_T_ex,
    Gamma_tilde_ac_F_ex,
    rho,
    delta_i,
    Code_Gamma_ac_co,
    Code_Gamma_ac_in,
    Code_Gamma_ac_ex,
    Code_Gamma_ec_co,
    Code_Gamma_ec_in,
    Code_Gamma_ec_ex,
    # Binary Variables
    delta_dg_C,
    delta_dg_T,
    delta_dg_F,
    delta_chi_T,
    delta_chi_F,
    delta_beta_C,
    delta_beta_T,
    delta_beta_F,
    delta_beta_plus,
    delta_beta_minus,
    delta_beta_in,
    delta_ac_C,
    delta_ac_T,
    delta_ac_F,
    delta_ac_Cp,
    delta_ac_Tp,
    delta_ac_Fp,
    delta_ac_CT,
    delta_ac_TC,
    delta_ac_CF,
    delta_ac_TF,
    delta_ec_C,
    delta_ec_T,
    delta_ec_F,
    delta_ec_Cp,
    delta_ec_Tp,
    delta_ec_Fp,
    delta_ec_CT,
    delta_ec_TC,
    delta_ec_CF,
    delta_ec_TF,
    e_C,
    e_T,
    e_F,
    v_C,
    v_T,
    v_F,
    delta_dg_CT,
    delta_dg_TC,
    delta_dg_CF,
    delta_dg_TF,
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
    ec_Cp,
    ec_Tp,
    ec_Fp,
    ec_CT,
    ec_TC,
    ec_CF,
    ec_TF,
    ec_co,
    ec_in,
    ec_ex,
    ec_nc,
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
    # -------- Constraint (88) --------
    for gamma in Gamma_co:
        if gamma not in Gamma_tilde_ec_C:
            MILP += ec_C[gamma] == 0, f"milp-(88)-1-{gamma}"
        if gamma not in Gamma_tilde_ec_T:
            MILP += ec_T[gamma] == 0, f"milp-(88)-2-{gamma}"
        if gamma not in Gamma_tilde_ec_CT:
            MILP += ec_CT[gamma] == 0, f"milp-(88)-5-{gamma}"
        if gamma not in Gamma_tilde_ec_TC:
            MILP += ec_TC[gamma] == 0, f"milp-(88)-6-{gamma}"
    for gamma in Gamma_in:
        if gamma not in Gamma_tilde_ec_F:
            MILP += ec_F[gamma] == 0, f"milp-(88)-3-{gamma}"
        if gamma not in Gamma_tilde_ec_CF:
            MILP += ec_CF[gamma] == 0, f"milp-(88)-7-{gamma}"
        if gamma not in Gamma_tilde_ec_TF:
            MILP += ec_TF[gamma] == 0, f"milp-(88)-8-{gamma}"
    for gamma in Gamma_ex:
        if gamma not in Gamma_tilde_ec_C_ex:
            for p in range(1, rho + 1):
                MILP += ec_Cp[(p, gamma)] == 0, f"milp-(88)-4-C-{p}-{gamma}"
        if gamma not in Gamma_tilde_ec_T_ex:
            for p in range(1, rho + 1):
                MILP += ec_Tp[(p, gamma)] == 0, f"milp-(88)-4-T-{p}-{gamma}"
        if gamma not in Gamma_tilde_ec_F_ex:
            for p in range(1, rho + 1):
                MILP += ec_Fp[(p, gamma)] == 0, f"milp-(88)-4-F-{p}-{gamma}"

    # -------- Constraint (89) --------
    for m in range(1, MAX_BOND + 1):
        MILP += pulp.lpSum(ec_C[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_co if m1 == m) == \
                pulp.lpSum(delta_beta_C[(i, m)] for i in range(k_C_tilde + 1, m_C + 1)), f"milp-(89)-1-{m}"
        MILP += pulp.lpSum(ec_T[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_co if m1 == m) == \
                pulp.lpSum(delta_beta_T[(i, m)] for i in range(2, t_T + 1)), f"milp-(89)-2-{m}"
        MILP += pulp.lpSum(ec_F[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_in if m1 == m) == \
                pulp.lpSum(delta_beta_F[(i, m)] for i in range(2, t_F + 1)), f"milp-(89)-3-{m}"
        for p in range(1, rho + 1):
            MILP += pulp.lpSum(ec_Cp[(p, (mu1, mu2, m1))] for (mu1, mu2, m1) in Gamma_ex if m1 == m) == \
                    pulp.lpSum(delta_beta_C[(i, j, m)] for i in range(1, t_C + 1) for j in Dsn_C[delta_i[i]][p]), f"milp-(89)-4-C-{m}-{p}"
            MILP += pulp.lpSum(ec_Tp[(p, (mu1, mu2, m1))] for (mu1, mu2, m1) in Gamma_ex if m1 == m) == \
                    pulp.lpSum(delta_beta_T[(i, j, m)] for i in range(1, t_T + 1) for j in Dsn_T[p]), f"milp-(89)-4-T-{m}-{p}"
            MILP += pulp.lpSum(ec_Fp[(p, (mu1, mu2, m1))] for (mu1, mu2, m1) in Gamma_ex if m1 == m) == \
                    pulp.lpSum(delta_beta_F[(i, j, m)] for i in range(1, t_F + 1) for j in Dsn_F[p]), f"milp-(89)-4-F-{m}-{p}"
        MILP += pulp.lpSum(ec_CT[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_co if m1 == m) == \
                pulp.lpSum(delta_beta_plus[(k, m)] for k in range(1, k_C + 1)), f"milp-(89)-5-{m}"
        MILP += pulp.lpSum(ec_TC[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_co if m1 == m) == \
                pulp.lpSum(delta_beta_minus[(k, m)] for k in range(1, k_C + 1)), f"milp-(89)-6-{m}"
        MILP += pulp.lpSum(ec_CF[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_in if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(1, t_C_tilde + 1)), f"milp-(89)-7-{m}"
        MILP += pulp.lpSum(ec_TF[(mu1, mu2, m1)] for (mu1, mu2, m1) in Gamma_in if m1 == m) == \
                pulp.lpSum(delta_beta_in[(c, m)] for c in range(t_C_tilde + 1, c_F + 1)), f"milp-(89)-8-{m}"

    # -------- Constraint (90) --------
    for i in range(k_C_tilde + 1, m_C + 1):
        MILP += pulp.lpSum(Code_Gamma_ac_co[(a, b, m)] * delta_ec_C[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_C) == \
                pulp.lpSum(Code_Gamma_ac_co[nu] * delta_ac_C[(i, nu)] 
                    for nu in Gamma_tilde_ac_C), f"milp-(90)-1-{i}"
        MILP += Delta_ec_C_plus[i] + pulp.lpSum(d * delta_ec_C[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_C) == deg_C[(tail_C[i], 0)], f"milp-(90)-2-{i}"
        MILP += Delta_ec_C_minus[i] + pulp.lpSum(d * delta_ec_C[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_C) == deg_C[(head_C[i], 0)], f"milp-(90)-3-{i}"
        MILP += Delta_ec_C_plus[i] + Delta_ec_C_minus[i] <= 8 * (1 - e_C[i]), f"milp-(90)-4-{i}"
    for gamma in Gamma_tilde_ec_C:
        MILP += pulp.lpSum(delta_ec_C[(i, gamma)] 
                for i in range(k_C_tilde + 1, m_C + 1)) == ec_C[gamma], f"milp-(90)-5-{gamma}"
    
    # -------- Constraint (91) --------
    for i in range(2, t_T + 1):
        MILP += pulp.lpSum(Code_Gamma_ac_co[(a, b, m)] * delta_ec_T[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_T) == \
                pulp.lpSum(Code_Gamma_ac_co[nu] * delta_ac_T[(i, nu)] 
                    for nu in Gamma_tilde_ac_T), f"milp-(91)-1-{i}"
        MILP += Delta_ec_T_plus[i] + pulp.lpSum(d * delta_ec_T[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_T) == deg_T[(i - 1, 0)], f"milp-(91)-2-{i}"
        MILP += Delta_ec_T_minus[i] + pulp.lpSum(d * delta_ec_T[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_T) == deg_T[(i, 0)], f"milp-(91)-3-{i}"
        MILP += Delta_ec_T_plus[i] + Delta_ec_T_minus[i] <= 8 * (1 - e_T[i]), f"milp-(91)-4-{i}"
    for gamma in Gamma_tilde_ec_T:
        MILP += pulp.lpSum(delta_ec_T[(i, gamma)] 
                for i in range(2, t_T + 1)) == ec_T[gamma], f"milp-(91)-5-{gamma}"
  
    # -------- Constraint (92) --------
    for i in range(2, t_F + 1):
        MILP += pulp.lpSum(Code_Gamma_ac_in[(a, b, m)] * delta_ec_F[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_F) == \
                pulp.lpSum(Code_Gamma_ac_in[nu] * delta_ac_F[(i, nu)] 
                    for nu in Gamma_tilde_ac_F), f"milp-(92)-1-{i}"
        MILP += Delta_ec_F_plus[i] + pulp.lpSum(d * delta_ec_F[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_F) == deg_F[(i - 1, 0)], f"milp-(92)-2-{i}"
        MILP += Delta_ec_F_minus[i] + pulp.lpSum(d * delta_ec_F[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_F) == deg_F[(i, 0)], f"milp-(92)-3-{i}"
        MILP += Delta_ec_F_plus[i] + Delta_ec_F_minus[i] <= 8 * (1 - e_F[i]), f"milp-(92)-4-{i}"
    for gamma in Gamma_tilde_ec_F:
        MILP += pulp.lpSum(delta_ec_F[(i, gamma)] 
                for i in range(2, t_F + 1)) == ec_F[gamma], f"milp-(92)-5-{gamma}"

    # -------- Constraint (93) --------
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            MILP += pulp.lpSum(Code_Gamma_ac_ex[(a, b, m)] * delta_ec_Cp[(i, j, ((a, d1), (b, d2), m))] 
                        for ((a, d1), (b, d2), m) in Gamma_tilde_ec_C_ex) == \
                    pulp.lpSum(Code_Gamma_ac_ex[nu] * delta_ac_Cp[(i, j, nu)] 
                        for nu in Gamma_tilde_ac_C_ex), f"milp-(93)-C-1-{i}-{j}"
            MILP += Delta_ec_C_plus[(i, j)] + pulp.lpSum(d * delta_ec_Cp[(i, j, ((a, d), b, m))]
                        for ((a, d), b, m) in Gamma_tilde_ec_C_ex) == deg_C[(i, prt_C[delta_i[i]][j])], f"milp-(93)-C-2-{i}-{j}"
            MILP += Delta_ec_C_minus[(i, j)] + pulp.lpSum(d * delta_ec_Cp[(i, j, (a, (b, d), m))]
                        for (a, (b, d), m) in Gamma_tilde_ec_C_ex) == deg_C[(i, j)], f"milp-(93)-C-3-{i}-{j}"
            MILP += Delta_ec_C_plus[(i, j)] + Delta_ec_C_minus[(i, j)] <= 8 * (1 - v_C[(i, j)]), f"milp-(93)-C-4-{i}-{j}"
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            MILP += pulp.lpSum(Code_Gamma_ac_ex[(a, b, m)] * delta_ec_Tp[(i, j, ((a, d1), (b, d2), m))] 
                        for ((a, d1), (b, d2), m) in Gamma_tilde_ec_T_ex) == \
                    pulp.lpSum(Code_Gamma_ac_ex[nu] * delta_ac_Tp[(i, j, nu)] 
                        for nu in Gamma_tilde_ac_T_ex), f"milp-(93)-T-1-{i}-{j}"
            MILP += Delta_ec_T_plus[(i, j)] + pulp.lpSum(d * delta_ec_Tp[(i, j, ((a, d), b, m))]
                        for ((a, d), b, m) in Gamma_tilde_ec_T_ex) == deg_T[(i, prt_T[j])], f"milp-(93)-T-2-{i}-{j}"
            MILP += Delta_ec_T_minus[(i, j)] + pulp.lpSum(d * delta_ec_Tp[(i, j, (a, (b, d), m))]
                        for (a, (b, d), m) in Gamma_tilde_ec_T_ex) == deg_T[(i, j)], f"milp-(93)-T-3-{i}-{j}"
            MILP += Delta_ec_T_plus[(i, j)] + Delta_ec_T_minus[(i, j)] <= 8 * (1 - v_T[(i, j)]), f"milp-(93)-T-4-{i}-{j}"
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            MILP += pulp.lpSum(Code_Gamma_ac_ex[(a, b, m)] * delta_ec_Fp[(i, j, ((a, d1), (b, d2), m))] 
                        for ((a, d1), (b, d2), m) in Gamma_tilde_ec_F_ex) == \
                    pulp.lpSum(Code_Gamma_ac_ex[nu] * delta_ac_Fp[(i, j, nu)] 
                        for nu in Gamma_tilde_ac_F_ex), f"milp-(93)-F-1-{i}-{j}"
            MILP += Delta_ec_F_plus[(i, j)] + pulp.lpSum(d * delta_ec_Fp[(i, j, ((a, d), b, m))]
                        for ((a, d), b, m) in Gamma_tilde_ec_F_ex) == deg_F[(i, prt_F[j])], f"milp-(93)-F-2-{i}-{j}"
            MILP += Delta_ec_F_minus[(i, j)] + pulp.lpSum(d * delta_ec_Fp[(i, j, (a, (b, d), m))]
                        for (a, (b, d), m) in Gamma_tilde_ec_F_ex) == deg_F[(i, j)], f"milp-(93)-F-3-{i}-{j}"
            MILP += Delta_ec_F_plus[(i, j)] + Delta_ec_F_minus[(i, j)] <= 8 * (1 - v_F[(i, j)]), f"milp-(93)-F-4-{i}-{j}"      
    for p in range(1, rho + 1):
        for gamma in Gamma_tilde_ec_C_ex:
            MILP += pulp.lpSum(delta_ec_Cp[(i, j, gamma)] 
                    for i in range(1, t_C + 1) for j in Dsn_C[delta_i[i]][p]) == ec_Cp[(p, gamma)], f"milp-(93)-C-5-{p}-{gamma}"
        for gamma in Gamma_tilde_ec_T_ex:
            MILP += pulp.lpSum(delta_ec_Tp[(i, j, gamma)] 
                    for i in range(1, t_T + 1) for j in Dsn_T[p]) == ec_Tp[(p, gamma)], f"milp-(93)-T-5-{p}-{gamma}"
        for gamma in Gamma_tilde_ec_F_ex:
            MILP += pulp.lpSum(delta_ec_Fp[(i, j, gamma)] 
                    for i in range(1, t_F + 1) for j in Dsn_F[p]) == ec_Fp[(p, gamma)], f"milp-(93)-F-5-{p}-{gamma}"
        
    # -------- Constraint (94) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += deg_T[(i, 0)] + MAX_VAL * (1 - delta_chi_T[(i, k)] + e_T[i]) >= \
                    deg_T_CT[k], f"milp-(94)-1-{k}-{i}"
            MILP += deg_T_CT[k] >= deg_T[(i, 0)] - \
                    MAX_VAL * (1 - delta_chi_T[(i, k)] + e_T[i]), f"milp-(94)-2-{k}-{i}"
        MILP += pulp.lpSum(Code_Gamma_ac_co[(a, b, m)] * delta_ec_CT[(k, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_CT) == \
                pulp.lpSum(Code_Gamma_ac_co[nu] * delta_ac_CT[(k, nu)] 
                    for nu in Gamma_tilde_ac_CT), f"milp-(94)-3-{k}"
        MILP += Delta_ec_CT_plus[k] + pulp.lpSum(d * delta_ec_CT[(k, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_CT) == deg_C[(tail_C[k], 0)], f"milp-(94)-4-{k}"
        MILP += Delta_ec_CT_minus[k] + pulp.lpSum(d * delta_ec_CT[(k, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_CT) == deg_T_CT[k], f"milp-(94)-5-{k}"
        MILP += Delta_ec_CT_plus[k] + Delta_ec_CT_minus[k] <= 8 * (1 - delta_chi_T[k]), f"milp-(94)-6-{k}"
    for gamma in Gamma_tilde_ec_CT:
        MILP += pulp.lpSum(delta_ec_CT[(k, gamma)] for k in range(1, k_C + 1)) == ec_CT[gamma], f"milp-(94)-7-{gamma}"

    # -------- Constraint (95) --------
    for k in range(1, k_C + 1):
        for i in range(1, t_T + 1):
            MILP += deg_T[(i, 0)] + MAX_VAL * (1 - delta_chi_T[(i, k)] + e_T[i + 1]) >= \
                    deg_T_TC[k], f"milp-(95)-1-{k}-{i}"
            MILP += deg_T_TC[k] >= deg_T[(i, 0)] - \
                    MAX_VAL * (1 - delta_chi_T[(i, k)] + e_T[i + 1]), f"milp-(95)-2-{k}-{i}"
        MILP += pulp.lpSum(Code_Gamma_ac_co[(a, b, m)] * delta_ec_TC[(k, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_TC) == \
                pulp.lpSum(Code_Gamma_ac_co[nu] * delta_ac_TC[(k, nu)] 
                    for nu in Gamma_tilde_ac_TC), f"milp-(95)-3-{k}"
        MILP += Delta_ec_TC_plus[k] + pulp.lpSum(d * delta_ec_TC[(k, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_TC) == deg_T_TC[k], f"milp-(95)-4-{k}"
        MILP += Delta_ec_TC_minus[k] + pulp.lpSum(d * delta_ec_TC[(k, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_TC) == deg_C[(head_C[k], 0)], f"milp-(95)-5-{k}"
        MILP += Delta_ec_TC_plus[k] + Delta_ec_TC_minus[k] <= 8 * (1 - delta_chi_T[k]), f"milp-(95)-6-{k}"
    for gamma in Gamma_tilde_ec_TC:
        MILP += pulp.lpSum(delta_ec_TC[(k, gamma)] for k in range(1, k_C + 1)) == ec_TC[gamma], f"milp-(95)-7-{gamma}"

    # -------- Constraint (96) --------
    for c in range(1, t_C_tilde + 1):
        for i in range(1, t_F + 1):
            MILP += deg_F[(i, 0)] + MAX_VAL * (1 - delta_chi_F[(i, c)] + e_F[i]) >= \
                    deg_F_CF[c], f"milp-(96)-1-{c}-{i}"
            MILP += deg_F_CF[c] >= deg_F[(i, 0)] - \
                    MAX_VAL * (1 - delta_chi_F[(i, c)] + e_F[i]), f"milp-(96)-2-{c}-{i}"
        MILP += pulp.lpSum(Code_Gamma_ac_in[(a, b, m)] * delta_ec_CF[(c, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_CF) == \
                pulp.lpSum(Code_Gamma_ac_in[nu] * delta_ac_CF[(c, nu)] 
                    for nu in Gamma_tilde_ac_CF), f"milp-(96)-3-{c}"
        MILP += Delta_ec_CF_plus[c] + pulp.lpSum(d * delta_ec_CF[(c, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_CF) == deg_C[(c, 0)], f"milp-(96)-4-{c}"
        MILP += Delta_ec_CF_minus[c] + pulp.lpSum(d * delta_ec_CF[(c, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_CF) == deg_F_CF[c], f"milp-(96)-5-{c}"
        MILP += Delta_ec_CF_plus[c] + Delta_ec_CF_minus[c] <= 8 * (1 - delta_chi_F[c]), f"milp-(96)-6-{c}"
    for gamma in Gamma_tilde_ec_CF:
        MILP += pulp.lpSum(delta_ec_CF[(c, gamma)] for c in range(1, t_C_tilde + 1)) == ec_CF[gamma], f"milp-(96)-7-{gamma}"

    # -------- Constraint (97) --------
    for i in range(1, t_T + 1):
        for j in range(1, t_F + 1):
            MILP += deg_F[(j, 0)] + MAX_VAL * (1 - delta_chi_F[(j, i + t_C_tilde)] + e_F[j]) >= \
                    deg_F_TF[i], f"milp-(97)-1-{i}-{j}"
            MILP += deg_F_TF[i] >= deg_F[(j, 0)] - \
                    MAX_VAL * (1 - delta_chi_F[(j, i + t_C_tilde)] + e_F[j]), f"milp-(97)-2-{i}-{j}"
        MILP += pulp.lpSum(Code_Gamma_ac_in[(a, b, m)] * delta_ec_TF[(i, ((a, d1), (b, d2), m))] 
                    for ((a, d1), (b, d2), m) in Gamma_tilde_ec_TF) == \
                pulp.lpSum(Code_Gamma_ac_in[nu] * delta_ac_TF[(i, nu)] 
                    for nu in Gamma_tilde_ac_TF), f"milp-(97)-3-{i}"
        MILP += Delta_ec_TF_plus[i] + pulp.lpSum(d * delta_ec_TF[(i, ((a, d), b, m))]
                    for ((a, d), b, m) in Gamma_tilde_ec_TF) == deg_T[(i, 0)], f"milp-(97)-4-{i}"
        MILP += Delta_ec_TF_minus[i] + pulp.lpSum(d * delta_ec_TF[(i, (a, (b, d), m))]
                    for (a, (b, d), m) in Gamma_tilde_ec_TF) == deg_F_TF[i], f"milp-(97)-5-{i}"
        MILP += Delta_ec_TF_plus[i] + Delta_ec_TF_minus[i] <= 8 * (1 - delta_chi_F[i + t_C_tilde]), f"milp-(97)-6-{i}"
    for gamma in Gamma_tilde_ec_TF:
        MILP += pulp.lpSum(delta_ec_TF[(i, gamma)] for i in range(1, t_T + 1)) == ec_TF[gamma], f"milp-(97)-7-{gamma}"

    # -------- Constraint (98) --------
    for (mu1, mu2, m) in Gamma_co_less:
        MILP += ec_C[(mu1, mu2, m)] + ec_C[(mu2, mu1, m)] + \
                ec_T[(mu1, mu2, m)] + ec_T[(mu2, mu1, m)] + \
                ec_CT[(mu1, mu2, m)] + ec_CT[(mu2, mu1, m)] + \
                ec_TC[(mu1, mu2, m)] + ec_TC[(mu2, mu1, m)] == ec_co[(mu1, mu2, m)], f"milp-(98)-1-{(mu1, mu2, m)}"
    for gamma in Gamma_co_equal:
        MILP += ec_C[gamma] + ec_T[gamma] + ec_CT[gamma] + ec_TC[gamma] == ec_co[gamma], f"milp-(98)-2-{gamma}"
    for gamma in Gamma_in:
        MILP += ec_F[gamma] + ec_CF[gamma] + ec_TF[gamma] == ec_in[gamma], f"milp-(98)-3-{gamma}"
    for gamma in Gamma_ex:
        MILP += pulp.lpSum(ec_Cp[(p, gamma)] + ec_Tp[(p, gamma)] + ec_Fp[(p, gamma)] 
                for p in range(1, rho + 1)) == ec_ex[gamma], f"milp-(98)-4-{gamma}"

    return MILP

# Add constraints of fringe tree tie breaker
def add_constraints_tie_breaker(
    # Model
    MILP,
    # Constants
    t_C,
    t_T,
    t_F,
    tie_breaker_C,
    tie_breaker_T,
    tie_breaker_F,
    delta_i,
    # Binary Variables
    # Integer Variables
    alpha_C,
    alpha_T,
    alpha_F,
    beta_C,
    beta_T,
    beta_F
):
    # -------- Constraint (99) --------
    # This m in this constraint should not exist, by Zhu
    for i in range(1, t_C + 1):
        for (j, h) in tie_breaker_C[delta_i[i]]:
            MILP += 4 * alpha_C[(i, j)] + beta_C[(i, j)] >= \
                    4 * alpha_C[(i, h)] + beta_C[(i, h)], f"milp-(99)-C-{i}-{j}-{h}"
    for i in range(1, t_T + 1):
        for (j, h) in tie_breaker_T:
            MILP += 4 * alpha_T[(i, j)] + beta_T[(i, j)] >= \
                    4 * alpha_T[(i, h)] + beta_T[(i, h)], f"milp-(99)-T-{i}-{j}-{h}"
    for i in range(1, t_F + 1):
        for (j, h) in tie_breaker_F:
            MILP += 4 * alpha_F[(i, j)] + beta_F[(i, j)] >= \
                    4 * alpha_F[(i, h)] + beta_F[(i, h)], f"milp-(99)-F-{i}-{j}-{h}"

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

    with open(outputfilename, "w", newline='\n') as f:
        f.write(f"{t_C}\n")
        # print(t_C)
        for i in range(t_C):
            f.write(f"{base_vertices[i]}\n")
            # print(base_vertices[i])
            f.write(f"{ch_LB[i + 1]} {ch_UB[i + 1]}\n")
            # print(f"{ch_LB[i + 1]} {ch_UB[i + 1]}\n")
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
    delta_i,
    Lambda,
    Lambda_co,
    Lambda_nc,
    prt_T,
    prt_C,
    prt_F,
    head_C,
    tail_C,
    I_equal_one,
    I_zero_one,
    I_ge_one,
    I_ge_two,

    # Variables
    n_G,
    v_C,
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
    delta_chi_T,
    delta_chi_F,
    e_T,
    e_F,

    outputfilename
    ):

    ''' A function to output sdf file.  '''
    n = round(n_G.value())

    graph_col = {i : " " for i in range(1, n + 1)}
    graph_adj = {(i, j) : 0 for i in range(1, n + 1) for j in range(1, n + 1)}

    graph_index = 0

    index_T = {(i, j) : 0 for i in range(1, t_T + 1) for j in range(n_T + 1)}
    index_C = {(i, j) : 0 for i in range(1, t_C + 1) for j in range(n_C[delta_i[i]] + 1)}
    index_F = {(i, j) : 0 for i in range(1, t_F + 1) for j in range(n_F + 1)}

    for i in range(1, t_T + 1):
        for j in range(n_T + 1):
            if round(v_T[(i, j)].value()) != 0:
                graph_index += 1
                index_T[(i, j)] = graph_index
                if j == 0:
                    graph_col[graph_index] = Lambda_co[round(alpha_T[(i, j)].value()) - 1]
                else:
                    graph_col[graph_index] = Lambda_nc[round(alpha_T[(i, j)].value()) - 1]
    for i in range(1, t_C + 1):
        for j in range(n_C[delta_i[i]] + 1):
            if round(v_C[(i, j)].value()) != 0:
                graph_index += 1
                index_C[(i, j)] = graph_index
                if j == 0:
                    graph_col[graph_index] = Lambda_co[round(alpha_C[(i, j)].value()) - 1]
                else:
                    graph_col[graph_index] = Lambda_nc[round(alpha_C[(i, j)].value()) - 1]
    for i in range(1, t_F + 1):
        for j in range(n_F + 1):
            if round(v_F[(i, j)].value()) != 0:
                graph_index += 1
                index_F[(i, j)] = graph_index
                graph_col[graph_index] = Lambda_nc[round(alpha_F[(i, j)].value()) - 1]

    m = 0

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
                if round(delta_chi_T[(i, c)].value()) == 1:
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
                if round(delta_chi_T[(i, c)].value()) == 1:
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
                if round(delta_chi_F[(i, c)].value()) == 1:
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
                if round(delta_chi_F[(i, c)].value()) == 1:
                    m += 1
                    ind1 = index_T[(c - t_C_tilde, 0)]
                    ind2 = index_F[(i, 0)]
                    graph_adj[(ind1, ind2)] = mul
                    graph_adj[(ind2, ind1)] = mul 
                    # print("beta_TF", ind1, ind2)
                    break
    
    for i in range(1, t_T + 1):
        for j in range(1, n_T + 1):
            mul = round(beta_T[(i, j)].value())
            if mul != 0:
                m += 1
                ind1 = index_T[(i, prt_T[j])]
                ind2 = index_T[(i, j)]
                graph_adj[(ind1, ind2)] = mul
                graph_adj[(ind2, ind1)] = mul
                # print("beta_T_FT", ind1, ind2)
    for i in range(1, t_C + 1):
        for j in range(1, n_C[delta_i[i]] + 1):
            mul = round(beta_C[(i, j)].value())
            if mul != 0:
                m += 1
                ind1 = index_C[(i, prt_C[delta_i[i]][j])]
                ind2 = index_C[(i, j)]
                graph_adj[(ind1, ind2)] = mul
                graph_adj[(ind2, ind1)] = mul
                # print("beta_C_FT", ind1, ind2)
    for i in range(1, t_F + 1):
        for j in range(1, n_F + 1):
            mul = round(beta_F[(i, j)].value())
            if mul != 0:
                m += 1
                ind1 = index_F[(i, prt_F[j])]
                ind2 = index_F[(i, j)]
                graph_adj[(ind1, ind2)] = mul
                graph_adj[(ind2, ind1)] = mul
                # print("beta_F_FT", ind1, ind2)

    with open(outputfilename, "w", newline='\n') as f:
        f.write("1\n")
        f.write("MILP_cyclic\n")
        f.write("\n")
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
    

