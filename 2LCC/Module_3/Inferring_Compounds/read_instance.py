#global
Cmax=6
Cmin=4
max_degree = 4
small_val = 5
large_val = 1e5
eps = 1e-5

E_Cu_mu = dict() # E_Cu_mu = {1: [1, 6, 7, 8], 2: [1, 2], 3: [2, 3], 4: [3, 4, 7], 5: [4, 5, 8], 6: [5, 6]}
E_Cu = dict()   # E_Cu=[1:(1,2),2:(2,3),3: (3,4),4: (4,5),5: (5,6),6: (1,6),7: (1,4),8: (1,5)]

# E_cu and E_cu_mu
for i in range(1,1+2*Cmax-Cmin):
    if i <= Cmax-1:
        E_Cu[i] = (i,i+1)
    elif i == Cmax:
        E_Cu[i] = (1,i)
    else:
        E_Cu[i] = (1,i-Cmax+Cmin-1)
for mu in range(1,1+Cmax):
    E_Cu_mu[mu] = [ei for ei in E_Cu if mu in E_Cu[ei]]

def element_info():
    ''' A function to prepare information of chemical elements. '''
    set_element = dict() # element symbol : (val, mass)
    set_element["B3"] = (3, 108)
    set_element["C4"] = (4,120)
    set_element["N3"] = (3,140)
    set_element["O2"] = (2,160)
    set_element["F1"] = (1, 190)
    set_element["Si4"] = (4, 280)
    set_element["P5"] = (5,309)
    set_element["S2"] = (2,320)
    set_element["S6"] = (6,320)
    set_element["Cl1"] = (1, 354)
    set_element["V3"] = (3, 509)
    set_element["Br1"]= ( 1, 799)
    set_element["Cd2"] =( 2, 1124)
    set_element["I1"]= (1, 1269)
    set_element["Hg2"] = (2, 2006)
    set_element["K1"] = (1,391)
    set_element["Pb2"] = (2, 2072)
    set_element["Al3"] = (3, 269)
    set_element["H1"] = (1,10)
   # set_element["e*"] = (2, 0)
    return set_element
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

def prepare_gamma_lf_ac(set_Lambda):
    gamma_lf_ac = list()
    for a1 in set_Lambda:
        for a2 in set_Lambda:
            for m in range(1, min(set_Lambda[a1][0], set_Lambda[a2][0]) + 1):
                if m > 3:
                    continue
                ac_lf_tmp = (a1, a2, m)
                gamma_lf_ac.append(ac_lf_tmp)
    return gamma_lf_ac

def read_seedgraph(seed_graph_path, gamma_lf_ac, set_Lambda):
    E_T = list()
    V_T = list()
    Vo = list()
    Vo_bar = list()
    Eo = list() #list Eo[Ring edge]
    Eo_bar = list()
    Eo_u = dict() #dict,Eo_u{u:[RingEdge(u)]}
    Eo_bar_u = dict() #dict, EoLine_u{u:[NonRingEdge(u)]}
    Xi_T = list()
    Xi_u = dict()
    Lambda = dict() # { element: code }
    Lambda_int = list()
    gammaInt_ac = dict() # {ac: code}
    gammaInt_ec = dict() # {ec: code}
    leaf_of_non_ring_node = list()
    N_T_v = dict()#dict,{v:N_T(v)}
    val_a = dict()  # a:valence of a
    mass_star_a = dict() # a: mass of a
    naInt_LB_a = dict()# dict,{a:xx}
    naInt_UB_a = dict() # dict,{a:xx}
    naEx_LB_a = dict() # dict
    naEx_UB_a = dict() # dict
    na_LB_a = dict()
    na_UB_a = dict()
    acLf_LB = dict()
    acLf_UB = dict()
    acInt_LB = dict()
    acInt_UB = dict()
    ecInt_LB = dict()
    ecInt_UB = dict()
    #n_LB = None
    #n_UB = None
    M1 = large_val # large number
    epsilon1 = small_val # small number
    #n_int_LB = None
    #n_int_UB = None
    M_ms = large_val  # large number

    # -------- read seed graph instance --------
    with open(seed_graph_path, 'r') as file:
        file.readline()
        file.readline()
        num_node = int(file.readline())
        for i in range(1, num_node+1):
            V_T.append(i)
        file.readline()
        num_Vo = int(file.readline())
        file.readline()
        for i in range(num_Vo):
            ring_node = int(file.readline())
            Vo.append(ring_node)
        Vo_bar.extend([x for x in V_T if x not in Vo])
        file.readline()
        num_edge = int (file.readline()) #the number of edge
        adjacency_list = [[] for _ in range(num_node)] #creat adjacency list of seed graph

        file.readline()   # ring edge
        num_ring_edge = int (file.readline())
        for i in range(num_ring_edge):
            edge = tuple(map(int,file.readline().strip().split(" ")))
            adjacency_list[edge[0]-1].append(edge[1]-1)
            adjacency_list[edge[1]-1].append(edge[0]-1)
            E_T.append(edge)
            Eo.append(edge)
        for u in V_T:
            ring_edge_u =[]
            Eo_u[u] = ring_edge_u
            for e in Eo:
                if u in e:
                    Eo_u[u].append(e)

        file.readline()   # non-ring edge
        num_non_ring_edge = int(file.readline())
        for i in range(num_non_ring_edge):
            edge = tuple(map(int, file.readline().strip().split(" ")))
            adjacency_list[edge[0] - 1].append(edge[1] - 1)
            adjacency_list[edge[1] - 1].append(edge[0] - 1)
            E_T.append(edge)
            Eo_bar.append(edge)
        for u in V_T:
            non_ring_edge_u =[]
            Eo_bar_u[u] = non_ring_edge_u
            for eP in Eo_bar:
                if u in eP:
                    Eo_bar_u[u].append(eP)
        for v in Vo_bar:
            if len(adjacency_list[v-1]) == 1:
                leaf_of_non_ring_node.append(v)
            N_T_v[v] = [x+1 for x in adjacency_list[v - 1]]

        file.readline()
        n_int_LB = int(file.readline())
        n_int_UB = int(file.readline())
        file.readline()
        file.readline()
        n_LB = int(file.readline())
        n_UB = int(file.readline())
        file.readline()
        rho = int(file.readline())
        file.readline()
        file.readline()
        line = file.readline()
        clean_line = line.strip()
        A = clean_line.split(" ")
        for ele in A:
            try:
                mass_star_a[ele] = set_Lambda[ele][1]
            except:
                exit(f"Do not find element '{ele}'")
        for i in range(len(A)):
            val_a[A[i]] = int(A[i][-1])
            Lambda[A[i]] = i+1

        file.readline()  # Gamma_int_ac
        Int_ac =[] #just for acInt_LB_UB
        Int_ac_bar = [] #just for acInt_LB_UB
        tmp_gammaInt_ac =[]
        num_gammaInt_ac =int(file.readline())
        for i in range(num_gammaInt_ac):
            clean_line = file.readline().strip().split(" ") #remove \n
            ac = (clean_line[0],clean_line[1],int(clean_line[2]))
            ac_bar = (ac[1],ac[0],ac[2])
            tmp_gammaInt_ac.append(ac)
            if ac_bar != ac:
                tmp_gammaInt_ac.append(ac_bar)
            Int_ac.append(ac)
            Int_ac_bar.append(ac_bar)
        for i in range(len(tmp_gammaInt_ac)):
            gammaInt_ac[tmp_gammaInt_ac[i]] = i+1

        file.readline() # Gamma_int_ec
        Int_ec = []  #just for ecInt_LB_UB
        Int_ec_bar = []  #just for ecInt_LB_UB
        tmp_gammaInt_ec = []
        num_gammaInt_ec = int(file.readline())
        for i in range(num_gammaInt_ec):
            line = file.readline()
            clean_line = line.strip().split(" ")
            tmp_ec = tuple(int(item) if item.isdigit() else item for item in clean_line)
            ec = ((tmp_ec[0],tmp_ec[1]),(tmp_ec[2],tmp_ec[3]),tmp_ec[4])
            ec_bar = (ec[1],ec[0],ec[2])
            tmp_gammaInt_ec.append(ec)
            if ec_bar != ec:
                tmp_gammaInt_ec.append(ec_bar)
            Int_ec.append(ec)
            Int_ec_bar.append(ec_bar)
        for i in range(len(tmp_gammaInt_ec)):
            gammaInt_ec[tmp_gammaInt_ec[i]] = i+1

        file.readline()  # na_LB and na_UB
        num_na_LB_UB = int(file.readline())
        for i in range(num_na_LB_UB):
            clean_line =file.readline().strip().split(" ")
            na_LB_a[clean_line[0]] = int(clean_line[1]) #na_LB
            na_UB_a[clean_line[0]] = int(clean_line[2]) #na_UB

        file.readline()# naInt_LB and naInt_UB
        num_naInt_LB_UB = int(file.readline())
        for i in range(num_naInt_LB_UB):
            clean_line = file.readline().strip().split(" ")
            naInt_LB_a[clean_line[0]] = int(clean_line[1]) # naInt_LB
            naInt_UB_a[clean_line[0]] = int(clean_line[2]) # naInt_UB
            Lambda_int.append(clean_line[0])
        for ele in Lambda:
            if ele not in naInt_LB_a:
                naInt_LB_a[ele] = 0
                naInt_UB_a[ele] = 0

        file.readline() # na_LB_ex and na_UB_ex
        num_naEx_LB_UB = int(file.readline())
        for i in range(num_naEx_LB_UB):
            clean_line = file.readline().strip().split(" ")
            naEx_LB_a[clean_line[0]] = int(clean_line[1])#  naEx_LB
            naEx_UB_a[clean_line[0]] = int(clean_line[2])  # naEx_UB
        for ele in Lambda:
            if ele not in naEx_LB_a:
                naEx_LB_a[ele] = 0
                naEx_UB_a[ele] = 0

        file.readline() # ac_LB_int and ac_UB_int
        num_acInt_LB_UB = int(file.readline())
        for i in range(num_acInt_LB_UB):
            clean_line = file.readline().strip().split(" ")
            acInt_LB[Int_ac[i]] = int(clean_line[1])
            acInt_UB[Int_ac[i]] = int(clean_line[2])
            acInt_LB[Int_ac_bar[i]] = int(clean_line[1])
            acInt_UB[Int_ac_bar[i]] = int(clean_line[2])

        file.readline() # ec_LB_int and ec_UB_int
        num_ecInt_LB_UB = int(file.readline())
        for i in range(num_ecInt_LB_UB):
            clean_line = file.readline().strip().split(" ")
            ecInt_LB[Int_ec[i]] = int(clean_line[1])
            ecInt_UB[Int_ec[i]] = int(clean_line[2])
            ecInt_LB[Int_ec_bar[i]] = int(clean_line[1])
            ecInt_UB[Int_ec_bar[i]] = int(clean_line[2])

        file.readline() # ac_LB_lf and ac_UB_lf
        num_ac_lf_LB_UB = int(file.readline())
        for i in range(num_ac_lf_LB_UB):
            clean_line = file.readline().strip().split(" ")
            acLf_LB[(clean_line[0],clean_line[1],int(clean_line[2]))] = int(clean_line[3])
            acLf_UB[(clean_line[0],clean_line[1],int(clean_line[2]))] = int(clean_line[4])

        clean_line = file.readline().strip().split(" ")#0 8
        for gamma in gamma_lf_ac:
            if gamma not in acLf_LB:
                acLf_UB[gamma] = int(clean_line[1])
                acLf_LB[gamma] = int(clean_line[0])

        file.readline()    ### cycle-configuration ###
        for u in Vo:
            Xi_u[u] = []
        while True:
            try:
                ID = list(map(int,file.readline().strip().split()))
                num_cc = int(file.readline())
                for i in range(num_cc):
                    cc = tuple(map(int,file.readline().strip().split()))
                    Xi_T.append(cc)
                    for id in ID:
                        Xi_u[id].append(cc)
            except:
                file.close()
                break

    return (V_T, E_T, Vo, Vo_bar, Eo, Eo_bar, Eo_u, Eo_bar_u, Xi_T, Xi_u,
            Lambda, gammaInt_ac, gammaInt_ec, leaf_of_non_ring_node,
            N_T_v, val_a, mass_star_a, naInt_LB_a, naInt_UB_a,
            naEx_LB_a, naEx_UB_a, na_LB_a, na_UB_a, acLf_LB, acLf_UB,
            acInt_LB, acInt_UB, ecInt_LB, ecInt_UB,
            n_LB, n_UB, M1, epsilon1, n_int_LB, n_int_UB, M_ms, rho, Lambda_int)

def prepare_fringe_trees(fringe_path, Lambda):
    fringe_set = list()
    strF = dict()

    fc_LB = dict()
    fc_UB = dict()

    with open(fringe_path, 'r') as f:
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
            flag = True
            for (a, d) in psi.vertex:
                if a not in Lambda:
                    flag = False
                    break
            if flag:
                strF[ind] = (str1, str2, str3)
                fringe_set.append(psi)

    Lambda_ex = list()
    for psi in fringe_set:
        for (a, d) in psi.vertex[1:]:
            if a not in Lambda_ex and a in Lambda:
                Lambda_ex.append(a)

    return fringe_set, Lambda_ex, strF, fc_LB, fc_UB

def read_fringe_tree(fv_fringe_tree_filename, strF):

    index_set = dict()
    with open(fv_fringe_tree_filename,'r') as f:
        lines = f.readlines()
        if lines[0][0] == 'F':
            lines = lines[1:]
        for line in lines:
            if len(line.split(",")) < 4:
                continue
            ind = int(line.split(",")[0])
            str1 = line.split(",")[1]
            str2 = line.split(",")[2]
            str3 = line.split(",")[3].replace('\n', '')
            for key, f in strF.items():
                if (str1, str2, str3) == f:
                    index_set[ind] = key
    return index_set

def function_of_fringe(fringe_set,V_T,gamma_lf_ac,Lambda,set_Lambda):
    F_u = dict()  # dict,  F_u{u:available FT(u)}
    F_t = dict()  # dict, code : fringe tree
    alpha_r_phi = dict()  # dict
    valex_F_phi = dict()  # dict
    eledeg_F_phi = dict()  # dict
    n_a_phi = dict()  # dict
    deg_H_bar_phi = dict()  # dict {phi:deg of root of phi}
    ms_F_phi = dict()  # dict,bijection: phi -> R
    ht_F_phi = dict() # dict,bijection: phi -> R
    n_H_bar_phi = dict()  # dict,bijection: phi -> R
    ac_lf_gamma_phi = dict() # dict,bijection: (gamma,phi) -> R

    Code_F = {phi: phi.index for phi in fringe_set}
    ac_phi_lf = {gamma: {Code_F[psi]: 0} for psi in fringe_set for gamma in gamma_lf_ac} # it equal to ac_lf_gamma_phi
    for phi in fringe_set:
        ac_psi_lf_tmp = {gamma: 0 for gamma in gamma_lf_ac}
        for i in range(1, len(phi.vertex)):
            if phi.vertex[i][0] == "H1":
                continue
            tmp = len([1 for j in phi.adj[i] if phi.vertex[j][0] != "H1"])
            if tmp == 1:
                symbol1 = phi.vertex[i][0]
                symbol2 = ""
                m = 0
                for j in phi.adj[i]:
                    if phi.vertex[j][0] != "H1":
                        symbol2 = phi.vertex[j][0]
                        m = phi.beta[i][j]
                        break
                ac_psi_lf_tmp[(symbol1, symbol2, m)] += 1
        for gamma, d in  ac_phi_lf.items():
            d[Code_F[phi]] = ac_psi_lf_tmp[gamma]
        for a in Lambda:
            n_a_phi[(a,phi.index)] = 0
        for u ,(ele ,dep) in enumerate(phi.vertex[1:]):
            if ele in Lambda:
                n_a_phi[(ele,phi.index)] += 1

        F_t[phi.index] = phi
        ht_F_phi[phi.index] = phi.height
        n_H_bar_phi[phi.index] = sum([1 for v in phi.vertex if v[0] != "H1"]) - 1
        alpha_r_phi[phi.index] = phi.root[0]
        valex_F_phi[phi.index] = sum(phi.beta[0][v] for v in phi.adj[0])
        eledeg_F_phi[phi.index] = phi.chg[0]
        deg_H_bar_phi[phi.index] = sum([1 for v in phi.adj[0] if phi.vertex[v][0] != "H1"]) #without "H"
        # deg_H_bar_phi[phi.index] = len(phi.adj[0]) #with "H"
        mass_of_phi = 0
        for chem_ele in phi.vertex:
            try:
                mass_of_phi += set_Lambda[chem_ele[0]][1]
            except:
                exit(f"Do not find element '{chem_ele[0]}'")
        ms_F_phi[phi.index] = mass_of_phi
    for u in V_T:
        F_u[u] = [phi.index for phi in fringe_set]

    for gamma, d in ac_phi_lf.items():
        for phi in d:
            ac_lf_gamma_phi[(gamma,phi)] = ac_phi_lf[gamma][phi]

    return  (F_u, F_t, alpha_r_phi, valex_F_phi,
             eledeg_F_phi, n_a_phi,
             deg_H_bar_phi, ms_F_phi, ht_F_phi,
             n_H_bar_phi, ac_lf_gamma_phi)


def prepare_ksi_delta_r_mu0(Xi_T):
    ksi_delta_r_mu = {(ksi, delta, r, mu) :set() for ksi in Xi_T for delta in {"P","M"} for r in range(1, 1 + Cmax) for mu in range(1, Cmax + 1)}
    ksi_P_r = {(ksi,r) : [] for ksi in Xi_T for r in range(1, 1 + max(ksi))}
    ksi_M_r = {(ksi,r) : [] for ksi in Xi_T for r in range(1, 1 + max(ksi))}
    for ksi in Xi_T:
        for l, r in enumerate(ksi):
            for mu0 in range(1 , 1 + len(ksi)):
                if mu0 + l == len(ksi):
                    ksi_P_r[(ksi, r)].append((mu0, mu0 + l ))
                    ksi_delta_r_mu[(ksi,"P",r,mu0)].add(mu0 + l)
                if mu0 - l == 0 or mu0 - l == len(ksi):
                    ksi_M_r[(ksi,r)].append((mu0, len(ksi)))
                    ksi_delta_r_mu[(ksi, "M", r, mu0)].add(len(ksi))
                if mu0 - l != 0 and mu0 - l != len(ksi):
                    ksi_M_r[(ksi,r)].append((mu0,(mu0 - l) % len(ksi)))
                    ksi_delta_r_mu[(ksi, "M", r, mu0)].add((mu0 - l) % len(ksi))
                if mu0 + l != len(ksi):
                    ksi_P_r[(ksi,r)].append((mu0,(mu0+l)%len(ksi)))
                    ksi_delta_r_mu[(ksi, "P", r, mu0)].add((mu0+l)%len(ksi))
        # for mu0 in range(1, len(ksi) + 1):
        #     for mu in range(1, len(ksi) + 1):
        #         l_P = (mu - mu0) % len(ksi)
        #         l_M = (mu0 - mu) % len(ksi)
        #         ksi_P_r[(ksi, ksi[l_P])].append((mu0, mu))
        #         ksi_M_r[(ksi, ksi[l_M])].append((mu0, mu))
        #         ksi_delta_r_mu[(ksi,"P",ksi[l_P],mu0)].append(mu)
        #         ksi_delta_r_mu[(ksi,"M",ksi[l_M],mu0)].append(mu)
    return ksi_delta_r_mu