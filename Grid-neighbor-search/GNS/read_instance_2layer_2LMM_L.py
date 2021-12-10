"""
read_instance_BH-cyclic.py
"""

'''
[seed graph]
   V_C        : "V_C"
   E_C        : "E_C"

[core specification]
   ell_LB     : "\ell_{\rm LB}"
   ell_UB     : "\ell_{\rm UB}"
   cs_LB      : "\textsc{cs}_{\rm LB}"
   cs_UB      : "\textsc{cs}_{\rm UB}"
'''

import sys

def read_pmax_file(filename):
    with open(filename,'r') as f:
        F = [line.rstrip('\n') for line in f if line[0]!='#']
    
    p_max = int(F.pop(0))
    s = F.pop(0)
    delta = list(map(float, s.split(' ')))
    
    s = F.pop(0)
    r = list(map(int, s.split(' ')))
    return p_max, delta, r
    
    


def read_seed_graph(filename):
    with open(filename,'r') as f:
        F = [line.rstrip('\n') for line in f if line[0]!='#']

    ### read V_C ###
    num_V_C = int(F.pop(0))
    V_C = tuple(range(1,num_V_C+1))
    
    ### read E_C ###
    num_E_C = int(F.pop(0))
    E_C = {}
    for e in range(num_E_C):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        E_C[arr[0]] = (arr[0], arr[1], arr[2])   # Add arr[0] to distinguish two edges with same starting and ending vertices, by Zhu
        
    ### read ell_LB and ell_UB ###
    ell_LB = {}
    ell_UB = {}
    for e in range(num_E_C):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        ell_LB[arr[0]] = arr[1]
        ell_UB[arr[0]] = arr[2]

    ### compute E_ge_two, E_ge_one, E_zero_one, E_equal_one ###
    E_ge_two = []
    E_ge_one = []
    E_zero_one = []
    E_equal_one = []
    I_ge_two = []
    I_ge_one = []
    I_zero_one = []
    I_equal_one = []
    for e in E_C:
        if ell_LB[e] >= 2:
            E_ge_two.append(E_C[e])
            I_ge_two.append(e)
        elif ell_LB[e] == 1 and ell_UB[e] >= 2:
            E_ge_one.append(E_C[e])
            I_ge_one.append(e)
        elif ell_LB[e] == 0 and ell_UB[e] == 1:
            E_zero_one.append(E_C[e])
            I_zero_one.append(e)
        elif ell_LB[e] == 1 and ell_UB[e] == 1:
            E_equal_one.append(E_C[e])
            I_equal_one.append(e)
        else:
            sys.stderr.write('error: a strange edge is found.\n')
            sys.exit(1)
            
    ### read n_LB_int and n_UB_int ###
    n_LB_int = int(F.pop(0))
    n_UB_int = int(F.pop(0))

    # read n_LB and n_star
    n_LB = int(F.pop(0))
    n_star = int(F.pop(0))

    # read rho
    rho = int(F.pop(0))

    ### read ch_LB and ch_UB ###
    ch_LB = {}
    ch_UB = {}
    for v in range(num_V_C):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        ch_LB[arr[0]] = arr[1]
        ch_UB[arr[0]] = arr[2]
    for e in range(len(E_ge_two + E_ge_one)):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        ch_LB[E_C[arr[0]]] = arr[1]
        ch_UB[E_C[arr[0]]] = arr[2]

    ### read bl_LB and bl_UB ###
    bl_LB = {}
    bl_UB = {}
    for v in range(num_V_C):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        bl_LB[arr[0]] = arr[1]
        bl_UB[arr[0]] = arr[2]
    for e in range(len(E_ge_two + E_ge_one)):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        bl_LB[E_C[arr[0]]] = arr[1]
        bl_UB[E_C[arr[0]]] = arr[2]

    # read Lambda
    s = F.pop(0)
    Lambda = list(s.split(' '))

    # read Lambda_dg_int
    s = F.pop(0)
    num = int(s)
    Lambda_dg_int = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        Lambda_dg_int.append((arr[0], int(arr[1])))

    # read Gamma_int_ac
    s = F.pop(0)
    num = int(s)
    Gamma_int_ac = list()
    nu_int = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        tmp_1 = (arr[0], arr[1], int(arr[2]))
        tmp_2 = (arr[1], arr[0], int(arr[2]))
        nu_int.append(tmp_1)
        if tmp_1 not in Gamma_int_ac:
            Gamma_int_ac.append(tmp_1)
        if tmp_2 not in Gamma_int_ac:
            Gamma_int_ac.append(tmp_2)

    # read Gamma_int
    s = F.pop(0)
    num = int(s)
    Gamma_int = list()
    gam_int = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        tmp_1 = ((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4]))
        tmp_2 = ((arr[2], int(arr[3])), (arr[0], int(arr[1])), int(arr[4]))
        gam_int.append(tmp_1)
        if tmp_1 not in Gamma_int:
            Gamma_int.append(tmp_1)
        if tmp_2 not in Gamma_int:
            Gamma_int.append(tmp_2)

    # read Lambda_star
    Lambda_star = {i: set() for i in range(1, num_V_C + 1)}
    for i in range(1, num_V_C + 1):
        s = F.pop(0)
        arr = list(s.split(' '))
        ind = int(arr[0])
        arr.pop(0)
        for a in arr:
            Lambda_star[ind].add(a)

    Lambda_int = list()
    # read na_LB and na_UB
    s = F.pop(0)
    num = int(s)
    na_LB = {}
    na_UB = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        na_LB[arr[0]] = int(arr[1])
        na_UB[arr[0]] = int(arr[2])

    # read na_LB_int and na_UB_int
    s = F.pop(0)
    num = int(s)
    na_LB_int = {}
    na_UB_int = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        na_LB_int[arr[0]] = int(arr[1])
        na_UB_int[arr[0]] = int(arr[2])
        Lambda_int.append(arr[0])

    # read ns_LB_int and ns_UB_int
    s = F.pop(0)
    num = int(s)
    ns_LB_int = {}
    ns_UB_int = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        ns_LB_int[(arr[0], int(arr[1]))] = int(arr[2])
        ns_UB_int[(arr[0], int(arr[1]))] = int(arr[3])

    # read ac_LB_int and ac_UB_int
    s = F.pop(0)
    num = int(s)
    ac_LB_int = {}
    ac_UB_int = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = nu_int[int(arr[0]) - 1]
        ac_LB_int[(a1, a2, m)] = int(arr[1])
        ac_LB_int[(a2, a1, m)] = int(arr[1])
        ac_UB_int[(a1, a2, m)] = int(arr[2])
        ac_UB_int[(a2, a1, m)] = int(arr[2])

    # read ec_LB_int and ec_UB_int
    s = F.pop(0)
    num = int(s)
    ec_LB_int = {}
    ec_UB_int = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = gam_int[int(arr[0]) - 1]
        ec_LB_int[(a1, a2, m)] = int(arr[1])
        ec_LB_int[(a2, a1, m)] = int(arr[1])
        ec_UB_int[(a1, a2, m)] = int(arr[2])
        ec_UB_int[(a2, a1, m)] = int(arr[2])

    # read bd2_LB and bd2_UB
    bd2_LB = {}
    bd2_UB = {}
    for e in range(len(E_C)):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        bd2_LB[E_C[arr[0]]] = arr[1]
        bd2_UB[E_C[arr[0]]] = arr[2]

    # read bd3_LB and bd3_UB
    bd3_LB = {}
    bd3_UB = {}
    for e in range(len(E_C)):
        s = F.pop(0)
        arr = list(map(int, s.split(' ')))
        bd3_LB[E_C[arr[0]]] = arr[1]
        bd3_UB[E_C[arr[0]]] = arr[2]

    # read ac_LB_lf and ac_UB_lf
    s = F.pop(0)
    num = int(s)
    ac_LB_lf = dict()
    ac_UB_lf = dict()
    for e in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        ac_LB_lf[(arr[0], arr[1], int(arr[2]))] = int(arr[3])
        ac_UB_lf[(arr[0], arr[1], int(arr[2]))] = int(arr[4])
    s = F.pop(0)
    arr = list(map(int, s.split(' ')))
    ac_LB_lf_common = arr[0]
    ac_UB_lf_common = arr[1]    
        
    ####################################
    # # Undefined constants for instances but used in MILP
    r_GC = num_E_C - (num_V_C - 1)
    dg_LB = [0,0,0,0,0]
    dg_UB = [n_star,n_star,n_star,n_star,n_star]
    
    return V_C, E_C, \
        E_ge_two, E_ge_one, E_zero_one, E_equal_one, \
        I_ge_two, I_ge_one, I_zero_one, I_equal_one, \
        ell_LB, ell_UB, n_LB_int, n_UB_int, \
        n_LB, n_star, rho, \
        ch_LB, ch_UB, bl_LB, bl_UB, \
        Lambda, Lambda_dg_int, Gamma_int_ac, Gamma_int, \
        Lambda_star, na_LB, na_UB, Lambda_int, \
        na_LB_int, na_UB_int, ns_LB_int, ns_UB_int, \
        ac_LB_int, ac_UB_int, ec_LB_int, ec_UB_int, \
        bd2_LB, bd2_UB, bd3_LB, bd3_UB, \
        dg_LB, dg_UB, ac_LB_lf, ac_UB_lf, ac_LB_lf_common, ac_UB_lf_common, r_GC

def get_value(filename):
    
    y_min = 0
    y_max = 0

    ind = 0

    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if len(line.split(",")) < 2:
                continue
            if line.split(",")[0] == "CID":
                continue
            if ind == 0:
                y_min = float(line.split(",")[1])
                y_max = float(line.split(",")[1])
                ind = 1
            else:
                y_tmp = float(line.split(",")[1])
                if y_tmp > y_max:
                    y_max = y_tmp
                if y_tmp < y_min:
                    y_min = y_tmp

    return y_min, y_max 


# prepare a set of chemical rooted tree
class chemicalRootedTree():
    def __init__(self):
        self.root = ("e", 0)
        self.index = 0
        self.vertex = []
        self.adj = []
        self.alpha = []
        self.beta = []
        self.height = 0
        self.chg = []

def prepare_fringe_trees(fringe_filename, Lambda):

# modified for 2LMM, 0527
    set_F = list()
    strF = dict()

    fc_LB = dict()
    fc_UB = dict()
    
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
                fc_LB[ind] = int(LB_tmp)
                UB_tmp = line.split(",")[5].replace('\n', '')
                UB_tmp = UB_tmp.replace(' ', '')
                fc_UB[ind] = int(UB_tmp)
            else:
                fc_LB[ind] = 0
                fc_UB[ind] = 10
            
            psi = chemicalRootedTree()
            seq1 = str1.split()
            seq2 = [int(mul) for mul in line.split(",")[2].split()]
            seq3 = [int(chg) for chg in line.split(",")[3].split()]
            psi.index = ind
            psi.vertex = [(seq1[j], int(seq1[j + 1])) for j in range(0, len(seq1), 2)]
            psi.root = psi.vertex[0]
            psi.height = max(psi.vertex[v][1] for v in range(len(psi.vertex)) if psi.vertex[v][0] != "H1")
            psi.adj = [set() for _ in range(len(psi.vertex))]
            psi.beta = [[0 for _ in range(len(psi.vertex))] for _ in range(len(psi.vertex))]
            psi.chg = [chg for chg in seq3]
            for j in range(len(seq2)):
                cld = j + 1
                prt = max(v for v in range(j + 1) if psi.vertex[v][1] == psi.vertex[cld][1] - 1)   
                psi.adj[prt].add(cld)
                psi.adj[cld].add(prt)
                psi.beta[prt][cld] = seq2[j]
                psi.beta[cld][prt] = seq2[j]
                # print(str(prt) + " " + str(cld) + " " + str(j) + " " + str(seq2[j]))
            flag = True
            for (a, d) in psi.vertex:
                if a not in Lambda:
                    flag = False
                    break
            if flag:
                strF[ind] = (str1, str2, str3)
                set_F.append(psi)

    Lambda_ex = list()
    for psi in set_F:
        for (a, d) in psi.vertex[1:]:
            if a not in Lambda_ex and a in Lambda:
                Lambda_ex.append(a)

    return set_F, Lambda_ex, strF, fc_LB, fc_UB
    
if __name__=="__main__":
    
    V_C, E_C, \
    E_ge_two, E_ge_one, E_zero_one, E_equal_one, \
    I_ge_two, I_ge_one, I_zero_one, I_equal_one, \
    ell_LB, ell_UB, n_LB_int, n_UB_int, \
    n_LB, n_star, rho, \
    ch_LB, ch_UB, bl_LB, bl_UB, \
    Lambda, Lambda_dg_int, Gamma_int_ac, Gamma_int, \
    Lambda_star, na_LB, na_UB, Lambda_int, \
    na_LB_int, na_UB_int, ns_LB_int, ns_UB_int, \
    ac_LB_int, ac_UB_int, ec_LB_int, ec_UB_int, \
    bd2_LB, bd2_UB, bd3_LB, bd3_UB, dg_LB, dg_UB = read_seed_graph(sys.argv[1])

    set_F, psi_epsilon, Code_F, n_psi, deg_r, \
    beta_r, atom_r, ht, Lambda_ex = prepare_fringe_trees(sys.argv[2])
    
    # print(V_C)
    # print(E_C)
    # print(E_ge_two)
    # print(E_ge_one)
    # print(E_zero_one)
    # print(E_equal_one)
    # print(ell_LB)
    # print(ell_UB)
    # print(bl_UB)

    for psi in (set_F + [psi_epsilon]):
        print(str(Code_F[psi]) + " " + str(n_psi[Code_F[psi]]) + " " + \
                str(ht[Code_F[psi]]) + " " + str(atom_r[Code_F[psi]]) + " " + \
                str(deg_r[Code_F[psi]]) + " " + str(beta_r[Code_F[psi]]))

    # print(Lambda_ex)

    # set_F_v = {v : set_F for v in V_C}
    # set_F_E = set_F
    # n_C = max(psi.numVertex - 1 for v in V_C for psi in set_F_v[v])
    # n_T = max(psi.numVertex - 1 for psi in set_F_E)
    # n_F = max(psi.numVertex - 1 for psi in set_F_E)
    # print(str(n_C) + " " + str(n_T) + " " + str(n_F))

    MAX_VAL = 4
    val = {"C": 4, "O": 2, "N": 3}

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

    print(n_H)
    print(na_alpha_ex)


