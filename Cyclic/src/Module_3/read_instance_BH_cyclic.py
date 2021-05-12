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
        s= F.pop(0)
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
            
    ### read cs_LB and cs_UB ###
    cs_LB = int(F.pop(0))
    cs_UB = int(F.pop(0))

    # read n_LB and n_star
    n_LB = int(F.pop(0))
    n_star = int(F.pop(0))

    # read rho
    rho = int(F.pop(0))

    # read dg4_UB_nc
    dg4_UB_nc = int(F.pop(0))

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

    # read Lambda
    s = F.pop(0)
    Lambda = list(s.split(' '))

    # read Lambda_dg_co
    s = F.pop(0)
    num = int(s)
    Lambda_dg_co = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        Lambda_dg_co.append((arr[0], int(arr[1])))

    # read Lambda_dg_nc
    s = F.pop(0)
    num = int(s)
    Lambda_dg_nc = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        Lambda_dg_nc.append((arr[0], int(arr[1])))

    # read Gamma_co_ac
    s = F.pop(0)
    num = int(s)
    Gamma_co_ac = list()
    nu_co = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        tmp_1 = (arr[0], arr[1], int(arr[2]))
        tmp_2 = (arr[1], arr[0], int(arr[2]))
        nu_co.append(tmp_1)
        if tmp_1 not in Gamma_co_ac:
            Gamma_co_ac.append(tmp_1)
        if tmp_2 not in Gamma_co_ac:
            Gamma_co_ac.append(tmp_2)
        # Gamma_co.append(((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4])))

    # read Gamma_in_ac
    s = F.pop(0)
    num = int(s)
    Gamma_in_ac = list()
    nu_in = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        tmp_1 = (arr[0], arr[1], int(arr[2]))
        # tmp_2 = (arr[1], arr[0], int(arr[2]))
        nu_in.append(tmp_1)
        if tmp_1 not in Gamma_in_ac:
            Gamma_in_ac.append(tmp_1)
        # if tmp_2 not in Gamma_in_ac:
        #     Gamma_in_ac.append(tmp_2)

    # read Gamma_ex_ac
    s = F.pop(0)
    num = int(s)
    Gamma_ex_ac = list()
    nu_ex = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        tmp_1 = (arr[0], arr[1], int(arr[2]))
        # tmp_2 = (arr[1], arr[0], int(arr[2]))
        nu_ex.append(tmp_1)
        if tmp_1 not in Gamma_ex_ac:
            Gamma_ex_ac.append(tmp_1)
        # if tmp_2 not in Gamma_ex_ac:
        #     Gamma_ex_ac.append(tmp_2)

    # read Gamma_co
    s = F.pop(0)
    num = int(s)
    Gamma_co = list()
    gam_co = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        tmp_1 = ((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4]))
        tmp_2 = ((arr[2], int(arr[3])), (arr[0], int(arr[1])), int(arr[4]))
        gam_co.append(tmp_1)
        if tmp_1 not in Gamma_co:
            Gamma_co.append(tmp_1)
        if tmp_2 not in Gamma_co:
            Gamma_co.append(tmp_2)
        # Gamma_co.append(((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4])))

    # read Gamma_in
    s = F.pop(0)
    num = int(s)
    Gamma_in = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        Gamma_in.append(((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4])))

    # read Gamma_ex
    s = F.pop(0)
    num = int(s)
    Gamma_ex = list()
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        Gamma_ex.append(((arr[0], int(arr[1])), (arr[2], int(arr[3])), int(arr[4])))

    Lambda_co = list()
    Lambda_nc = list()
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

    # read na_LB_co and na_UB_co
    s = F.pop(0)
    num = int(s)
    na_LB_co = {}
    na_UB_co = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        na_LB_co[arr[0]] = int(arr[1])
        na_UB_co[arr[0]] = int(arr[2])
        Lambda_co.append(arr[0])

    # read na_LB_nc and na_UB_nc
    s = F.pop(0)
    num = int(s)
    na_LB_nc = {}
    na_UB_nc = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        na_LB_nc[arr[0]] = int(arr[1])
        na_UB_nc[arr[0]] = int(arr[2])
        Lambda_nc.append(arr[0])

    # read ns_LB and ns_UB
    s = F.pop(0)
    num = int(s)
    ns_LB = {}
    ns_UB = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        ns_LB[(arr[0], int(arr[1]))] = int(arr[2])
        ns_UB[(arr[0], int(arr[1]))] = int(arr[3])

    # read ns_LB_co and ns_UB_co
    s = F.pop(0)
    num = int(s)
    ns_LB_co = {}
    ns_UB_co = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        ns_LB_co[(arr[0], int(arr[1]))] = int(arr[2])
        ns_UB_co[(arr[0], int(arr[1]))] = int(arr[3])

    # read ns_LB_nc and ns_UB_nc
    s = F.pop(0)
    num = int(s)
    ns_LB_nc = {}
    ns_UB_nc = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        ns_LB_nc[(arr[0], int(arr[1]))] = int(arr[2])
        ns_UB_nc[(arr[0], int(arr[1]))] = int(arr[3])

    # read ac_LB_co and ac_UB_co
    s = F.pop(0)
    num = int(s)
    ac_LB_co = {}
    ac_UB_co = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = nu_co[int(arr[0]) - 1]
        ac_LB_co[(a1, a2, m)] = int(arr[1])
        ac_LB_co[(a2, a1, m)] = int(arr[1])
        ac_UB_co[(a1, a2, m)] = int(arr[2])
        ac_UB_co[(a2, a1, m)] = int(arr[2])

    # read ac_LB_in and ac_UB_in
    s = F.pop(0)
    num = int(s)
    ac_LB_in = {}
    ac_UB_in = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = nu_in[int(arr[0]) - 1]
        ac_LB_in[(a1, a2, m)] = int(arr[1])
        # ac_LB_in[(a2, a1, m)] = int(arr[1])
        ac_UB_in[(a1, a2, m)] = int(arr[2])
        # ac_UB_in[(a2, a1, m)] = int(arr[2])

    # read ac_LB_ex and ac_UB_ex
    s = F.pop(0)
    num = int(s)
    ac_LB_ex = {}
    ac_UB_ex = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = nu_ex[int(arr[0]) - 1]
        ac_LB_ex[(a1, a2, m)] = int(arr[1])
        # ac_LB_ex[(a2, a1, m)] = int(arr[1])
        ac_UB_ex[(a1, a2, m)] = int(arr[2])
        # ac_UB_ex[(a2, a1, m)] = int(arr[2])

    # read ec_LB_co and ec_UB_co
    s = F.pop(0)
    num = int(s)
    ec_LB_co = {}
    ec_UB_co = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = gam_co[int(arr[0]) - 1]
        ec_LB_co[(a1, a2, m)] = int(arr[1])
        ec_LB_co[(a2, a1, m)] = int(arr[1])
        ec_UB_co[(a1, a2, m)] = int(arr[2])
        ec_UB_co[(a2, a1, m)] = int(arr[2])

    # read ec_LB_in and ec_UB_in
    s = F.pop(0)
    num = int(s)
    ec_LB_in = {}
    ec_UB_in = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = Gamma_in[int(arr[0]) - 1]
        ec_LB_in[(a1, a2, m)] = int(arr[1])
        ec_UB_in[(a1, a2, m)] = int(arr[2])

    # read ec_LB_ex and ec_UB_ex
    s = F.pop(0)
    num = int(s)
    ec_LB_ex = {}
    ec_UB_ex = {}
    for i in range(num):
        s = F.pop(0)
        arr = list(s.split(' '))
        a1, a2, m = Gamma_ex[int(arr[0]) - 1]
        ec_LB_ex[(a1, a2, m)] = int(arr[1])
        ec_UB_ex[(a1, a2, m)] = int(arr[2])

    # # read V_C_star
    # s = F.pop(0)
    # V_C_star = list(map(int, s.split(' ')))

    # # read alpha_star
    # alpha_star = {}
    # for i in range(len(V_C_star)):
    #     s = F.pop(0)
    #     arr = list(s.split(' '))
    #     alpha_star[int(arr[0])] = arr[1]

    # read Lambda_star
    Lambda_star = {i: set() for i in range(1, num_V_C + 1)}
    for i in range(1, num_V_C + 1):
        s = F.pop(0)
        arr = list(s.split(' '))
        ind = int(arr[0])
        arr.pop(0)
        for a in arr:
            Lambda_star[ind].add(a)

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
        
    ####################################
    # Undefined constants for instances but used in MILP
    r_GC = 5
    dg_LB = [0,0,0,0,0]
    dg_UB = [n_star,n_star,n_star,n_star,n_star]
    
    return V_C, E_C, E_ge_two, E_ge_one, E_zero_one, E_equal_one, \
        I_ge_two, I_ge_one, I_zero_one, I_equal_one, \
        ell_LB, ell_UB, cs_LB, cs_UB, bl_LB, bl_UB, ch_LB, ch_UB, \
        Lambda, Lambda_dg_co, Lambda_dg_nc, Gamma_co, Gamma_in, Gamma_ex, \
        ns_LB, ns_UB, ns_LB_co, ns_UB_co, ns_LB_nc, ns_UB_nc, Lambda_star, \
        bd2_LB, bd2_UB, bd3_LB, bd3_UB, dg4_UB_nc, n_LB, n_star, rho, r_GC, dg_LB, dg_UB, \
        na_LB, na_UB, na_LB_co, na_UB_co, na_LB_nc, na_UB_nc, Lambda_co, Lambda_nc, \
        ec_LB_co, ec_UB_co, ec_LB_in, ec_UB_in, ec_LB_ex, ec_UB_ex, \
        Gamma_co_ac, Gamma_in_ac, Gamma_ex_ac, \
        ac_LB_co, ac_UB_co, ac_LB_in, ac_UB_in, ac_LB_ex, ac_UB_ex

    
if __name__=="__main__":
    V_C, E_C, E_ge_two, E_ge_one, E_zero_one, E_equal_one,\
        ell_LB, ell_UB, cs_LB, cs_UB = read_seed_graph(sys.argv[1])
    print(V_C)
    print(E_C)
    print(E_ge_two)
    print(E_ge_one)
    print(E_zero_one)
    print(E_equal_one)
    print(ell_LB)
    print(ell_UB)
    print(cs_LB)
    print(cs_UB)
