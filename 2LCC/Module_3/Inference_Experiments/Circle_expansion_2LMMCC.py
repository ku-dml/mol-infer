import pulp,re
from read_instance import *
import pandas as pd

# -------- 3.1 Assigning Cycles to Ring Nodes --------
def assigning_cycles_to_ring_nodes(
        Vo,
        M1,
        epsilon1,
        Xi_T, #list, all the CC of T
        Xi_u, #dict,{u:list[available cc of u in Vo]} [(1,2,3,4,5),(1,2,3,4,5)]
        ksi_delta_r_mu,
        #model
        MILP):

    #-----prepare_variables-----#
    #Define y_u_[mu]
    y_u_mu = {(u,mu):pulp.LpVariable(f"y_{u}_{mu}",0,cat=pulp.LpContinuous)
            for u in Vo for mu in range(1,1+Cmax)}

    #Define z_u_r
    z_u_r = {(u,r):pulp.LpVariable(f"z_{u}_{r}",0,cat=pulp.LpContinuous)
            for u in Vo for r in range(1,Cmax+1)}

    #Define x_u_[ksi]_[mu0]_[delta]
    x_u_ksi_mu0_delta={(u,ksi,mu0,delta):pulp.LpVariable(f"x_{u}_{ksi}_{mu0}_{delta}",cat=pulp.LpBinary)
                       for u in Vo for ksi in Xi_u[u]
                       for mu0 in range(1,len(ksi)+1) for delta in ["P","M"]}

    #Define x_u_[ksi]
    x_u_ksi={(u,ksi):pulp.LpVariable(f"x_{u}_{ksi}",0,1,cat=pulp.LpBinary)
              for u in Vo for ksi in Xi_u[u]}

    #Define cc([ksi])
    cc_ksi={ksi:pulp.LpVariable(f"cc({ksi})",0,cat=pulp.LpInteger)
            for ksi in Xi_T}

    #-----add_constraints-----#
    # -------- Constraint (1) --------
    for u in Vo:
        for r in range(1, Cmax + 1):
            MILP+=z_u_r[(u,r)] >= 0 ,f"Milp-1_{u}_{r}-1"
            MILP+=z_u_r[(u,r)] <= M1 , f"Milp-1_{u}_{r}-2"

    # -------- Constraint (2) --------
    for u in Vo:
        for r in range(1,Cmax):
            MILP += z_u_r[(u,r)]+epsilon1 <= z_u_r[(u,r+1)],f"Milp-2_{u}_{r}"

    # -------- Constraint (3) --------
    for u in Vo:
        for ksi in Xi_u[u]:
            MILP += pulp.lpSum(x_u_ksi_mu0_delta[(u,ksi,mu0,delta)]
                              for mu0 in range(1,len(ksi)+1) for delta in ["P","M"]) == x_u_ksi[(u,ksi)],f"Milp-3_{u}_{ksi}"

    # -------- Constraint (4) --------
    for u in Vo:
        MILP += pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u]) == 1 ,f"MILP-4_{u}"
    # -------- Constraint (5) --------
    for u in Vo:
        for mu in range(1, Cmax + 1 ):
            for r in range(1,Cmax + 1):
                MILP += y_u_mu[(u,mu)] <= z_u_r[(u, r)] + M1 * (1 -
                pulp.lpSum(x_u_ksi_mu0_delta[(u,ksi,mu0,delta)] for ksi in Xi_u[u] for delta in ["P","M"] for mu0 in ksi_delta_r_mu[(ksi,delta,r,mu)])),f"MILP-5_{u}_{mu}_{r}"
    # -------- Constraint (6) --------
                MILP += y_u_mu[(u, mu)] >= z_u_r[(u, r)] - M1 * (1 -
                pulp.lpSum(x_u_ksi_mu0_delta[(u, ksi, mu0, delta)] for ksi in Xi_u[u] for delta in ["P", "M"] for mu0 in ksi_delta_r_mu[(ksi, delta, r, mu)])),f"MILP-6_{u}_{mu}_{r}"
    # -------- Constraint (7) --------
    for ksi in Xi_T:
        MILP += cc_ksi[ksi] == pulp.lpSum(x_u_ksi[(u,ksi)] for u in Vo if ksi in Xi_u[u]) # how many times CC(ksi) was used.

    return y_u_mu, z_u_r, x_u_ksi_mu0_delta, x_u_ksi, cc_ksi, MILP

# -------- 3.2 Attaching Two Cycles --------
def attaching_two_cycles(
        Xi_u,
        x_u_ksi,
        Vo, #list, Vo[RingNode]
        Eo_u, #dict,Eo_u{u:RingEdge(u)}
        Eo_bar_u, #dict, EoLine_u{u:NonRingEdge(u)}
        #model
        MILP):

    # -----prepare_variables-----#
    #Define e_u_i
    e_u_i={(u,i):pulp.LpVariable(f"e_{u}_{i}",0,cat=pulp.LpBinary)
           for u in Vo for i in range(1,2*Cmax-Cmin+1)}

    #Define xEdge_u_[e]_[v]
    xEdge_u_e_v = {(u,e,v):pulp.LpVariable(f"xEdge_{u}_{e}_{v}",0,cat=pulp.LpBinary)
                   for u in Vo for e in Eo_u[u] for v in range(1,2*Cmax-Cmin+1)}

    #Define xNode_u_[e']_[mu]
    xNode_u_ePrime_mu={(u,ep,mu):pulp.LpVariable(f"xNode_{u}_{ep}_{mu}",0,cat=pulp.LpBinary)
                       for u in Vo for ep in Eo_bar_u[u] for mu in range(1, Cmax + 1)}

    # -----add_constraints-----#

    # -------- Constraint (8) --------
    for u in Vo:
        for e in Eo_u[u]:
            MILP += pulp.lpSum(xEdge_u_e_v[(u,e,v)] for v in range(1,2*Cmax-Cmin+1)) == 1,f"MILP-8_{u}_{e}"

    # -------- Constraint (9) --------
    for u in Vo:
        for ePrime in Eo_bar_u[u]:
            MILP += pulp.lpSum(xNode_u_ePrime_mu[(u,ePrime,mu)] for mu in range(1,Cmax+1)) == 1 ,f"MILP-9_{u}_{ePrime}"
    # -------- Constraint (10) --------
    for u in Vo:
        for mu in range(1,Cmax+1):
            MILP += pulp.lpSum(xEdge_u_e_v[(u,e,v)] for e in Eo_u[u] for v in E_Cu_mu[mu]) <= 1,f"MILP-10_{u}_{mu}"
    # -------- Constraint (11) --------
    for u in Vo:
        for v in range(1,2*Cmax-Cmin+1):
            MILP += pulp.lpSum(xEdge_u_e_v[(u,e,v)] for e in Eo_u[u]) <= 1, f"MILP-11_{u}_{v}"
    # -------- Constraint (12) --------
    for u in Vo:
        for i in range(1,Cmin):
            MILP += e_u_i[(u,i)] == 1,f"MILP-12_{u}_{i}-1"
        for i in range(Cmin,Cmax-1):
            MILP+=pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) > i) == e_u_i[(u,i)] ,f"MILP-12_{u}_{i}-2"#e4 was used if 5,6-cycle
        for i in [Cmax-1,Cmax]:
            MILP += pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) == Cmax) == e_u_i[(u,i)],f"MILP-12_{u}_{i}-3"#e5,e6 was used in 6-cycle
        for i in range(Cmax+1,2*Cmax-Cmin+1):
            MILP += pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) == i-Cmax+Cmin-1) == e_u_i[(u,i)],f"MILP-12_{u}_{i}-4"#e7 should be used in 4-cycle,e8 should be used in 5-cycle
    # -------- Constraint (13) --------
    for u in Vo:
        for i in range(1, 2*Cmax-Cmin):
            for e in Eo_u[u]:
                MILP += xEdge_u_e_v[(u,e,i)] <=e_u_i[(u,i)],f"MILP-13_{u}_{i}_{e}" #就算用了某个circle edge，也不一定指派ring Edge
    # -------- Constraint (14) --------
    for u in Vo:
        for ePrime in Eo_bar_u[u]:
            for i in range(Cmin+1,Cmax+1):
                MILP +=  pulp.lpSum([x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) >= i]) >= xNode_u_ePrime_mu[(u,ePrime,i)],f"MILP-14_{u}_{ePrime}_{i}" #除非五元环以上，才会在5,6号circle顶点指派Non-ring-Edge

    return e_u_i, xEdge_u_e_v, xNode_u_ePrime_mu, MILP
# -------- 3.3 Constraints for Including Fringe-Trees --------
def constraint_for_including_fringe_tree(
        #Constants
        rho,
        leaf_of_non_ring_node,
        ms_F_phi,#dict,bijection: phi -> R
        ht_F_phi,#dict,bijection: phi -> R
        n_H_bar_phi,#dict,bijection: phi -> R
        ac_lf_gamma_phi,#dict,bijection: phi -> R
        Vo,# Ring node
        Vo_bar, # V(T)-Vo
        x_u_ksi,
        y_u_mu,
        xEdge_u_e_v,
        Xi_T, #list
        Xi_u,# dict
        Eo, #list Eo[Ring edge]
        Eo_u,#dict, {u:ring edge(u)}
        F_u, #dict,  F_u{u:available FT(u)}
        F_t, #list, all the fringe-tree
        n_LB,
        n_UB,
        n_int_LB,
        n_int_UB,
        gamma_lf_ac,
        acLf_LB,
        acLf_UB,
        #model
        MILP
):
    # -----prepare_variables-----#
    #Define delta_F(u,[mu];[fi])
    delta_F_u_mu_phi={(u,mu,phi):pulp.LpVariable(f"delta_F({u},{mu},{phi})",0,cat=pulp.LpBinary)
                    for u in Vo for mu in range(1,1+Cmax) for phi in F_u[u]}

    #Define delta_F(v;[phi])
    delta_F_v_phi={(v,phi):pulp.LpVariable(f"delta_F({v},{phi})",0,cat=pulp.LpBinary)
                   for v in Vo_bar for phi in F_u[v]}

    #Define fc([e];[phi])
    fc_e_phi={(e,phi):pulp.LpVariable(f"fc({e},{phi})",0,2,cat=pulp.LpInteger)
             for e in Eo for phi in F_t}

    #Define rank
    rank=pulp.LpVariable("rank",cat=pulp.LpInteger)

    #Define n_G,n_int
    n_G=pulp.LpVariable("n_G",lowBound=n_LB,upBound=n_UB,cat=pulp.LpInteger)
    n_int=pulp.LpVariable("n_int",lowBound=n_int_LB,upBound=n_int_UB,cat=pulp.LpInteger)

    #Define fc([phi])
    fc_phi={phi:pulp.LpVariable(f"fc({phi})",cat=pulp.LpInteger)
            for phi in F_t}

    #Define ac_lf([gamma])
    ac_lf_gamma={gamma: pulp.LpVariable(f"ac_lf({gamma})",lowBound=acLf_LB[gamma],upBound=acLf_UB[gamma],cat=pulp.LpInteger)
                 for gamma in gamma_lf_ac} #acLF_LB里面可能没有gamma

    # -----add_constraints-----#

    # -------- Constraint (15) --------
    for u in Vo:
        for mu in range(1,Cmax+1):
            MILP += (pulp.lpSum(delta_F_u_mu_phi[(u,mu,phi)] for phi in F_u[u])  ==
                     pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) >= mu)),f"MILP-15_{u}_{mu}"
    # -------- Constraint (16) --------
    for u in Vo:
        for mu in range(1,Cmax+1):
            MILP += pulp.lpSum(ms_F_phi[phi] * delta_F_u_mu_phi[(u, mu, phi)] for phi in F_u[u]) == y_u_mu[(u, mu)],f"MILP-16_{u}_{mu}"
    # -------- Constraint (17) --------
    for v in Vo_bar:
        MILP += pulp.lpSum(delta_F_v_phi[(v,phi)] for phi in F_u[v]) == 1,f"MILP-17_{v}-1"
        if v in leaf_of_non_ring_node:
            MILP += pulp.lpSum(delta_F_v_phi[(v,phi)] for phi in F_u[v] if ht_F_phi[phi] == rho) == 1,f"MILP-17_{v}-2"

    # -------- Constraint (18) --------
    for e in Eo: #u = e[0] u'=e[1],E_Cu[v][0]=i1,  E_Cu[v][1]=i2   E_Cu[vPrime][0]=j1 E_Cu[vPrime][1] =j2
        for v in range(1, 2 * Cmax - Cmin + 1):
            for vPrime in range(1, 2 * Cmax - Cmin + 1):
                MILP += (pulp.lpSum(phi * delta_F_u_mu_phi[(e[0],E_Cu[v][0],phi)] for phi in F_u[e[0]]) -
                         pulp.lpSum(phi * delta_F_u_mu_phi[(e[1],E_Cu[vPrime][1],phi)] for phi in F_u[e[1]]) <=
                         len(F_t)*(2-xEdge_u_e_v[(e[0],e,v)]-xEdge_u_e_v[(e[1],e,vPrime)])),f"MILP-18_{e}_{v}_{vPrime}-1"

                MILP += (pulp.lpSum(phi * delta_F_u_mu_phi[(e[0], E_Cu[v][0], phi)] for phi in F_u[e[0]]) -
                         pulp.lpSum(phi * delta_F_u_mu_phi[(e[1], E_Cu[vPrime][1], phi)] for phi in F_u[e[1]]) >=
                         len(F_t) * ( xEdge_u_e_v[(e[0], e, v)] + xEdge_u_e_v[(e[1], e, vPrime)]-2)),f"MILP-18_{e}_{v}_{vPrime}-2"

                MILP += (pulp.lpSum(phi * delta_F_u_mu_phi[(e[0], E_Cu[v][1], phi)] for phi in F_u[e[0]]) -
                         pulp.lpSum(phi * delta_F_u_mu_phi[(e[1], E_Cu[vPrime][0], phi)] for phi in F_u[e[1]]) <=
                         len(F_t) * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vPrime)])),f"MILP-18_{e}_{v}_{vPrime}-3"

                MILP += (pulp.lpSum(phi * delta_F_u_mu_phi[(e[0], E_Cu[v][1], phi)] for phi in F_u[e[0]]) -
                         pulp.lpSum(phi * delta_F_u_mu_phi[(e[1], E_Cu[vPrime][0], phi)] for phi in F_u[e[1]]) >=
                         len(F_t) * (xEdge_u_e_v[(e[0], e, v)] + xEdge_u_e_v[(e[1], e, vPrime)] - 2)),f"MILP-18_{e}_{v}_{vPrime}-4"
    # -------- Constraint (19) --------
    for e in Eo:
        for v in range(1,2*Cmax-Cmin+1):
            for phi in F_t:
                MILP += fc_e_phi[e,phi] - delta_F_u_mu_phi[(e[0],E_Cu[v][0],phi)] -delta_F_u_mu_phi[(e[0],E_Cu[v][1],phi)] <= 2*(1-xEdge_u_e_v[(e[0],e,v)]),f"MILP-19_{e}_{v}_{phi}-1"

                MILP += fc_e_phi[e,phi] - delta_F_u_mu_phi[(e[0],E_Cu[v][0],phi)] -delta_F_u_mu_phi[(e[0],E_Cu[v][1],phi)] >= 2*(xEdge_u_e_v[(e[0],e,v)]-1),f"MILP-19_{e}_{v}_{phi}-2"
    # -------- Constraint (20) --------
    for phi in F_t:
        MILP += fc_phi[phi] == pulp.lpSum([delta_F_u_mu_phi[(u,mu,phi)] for u in Vo for mu in range(1,1+Cmax) if (u,mu,phi) in delta_F_u_mu_phi]) + \
                pulp.lpSum(delta_F_v_phi[(v,phi)] for v in Vo_bar if (v, phi) in delta_F_v_phi) - pulp.lpSum(fc_e_phi[(e, phi)] for e in Eo),f"MILP-20_{phi}"
    # -------- Constraint (21) --------
    for gamma in gamma_lf_ac:
        MILP += ac_lf_gamma[gamma] == pulp.lpSum(ac_lf_gamma_phi[(gamma,phi)]*fc_phi[phi] for phi in F_t),f"MILP-21_{gamma}-1"
    # -------- Constraint (22) --------
    MILP += rank == len(Vo),f"MILP-22"
    # -------- Constraint (23) --------
    MILP += pulp.lpSum(len(ksi)*x_u_ksi[(u,ksi)] for u in Vo for ksi in Xi_u[u]) + len(Vo_bar) - 2 * len(Eo) == n_int,f"MILP-23"
    # -------- Constraint (24) --------
    MILP += n_G==n_int + pulp.lpSum(n_H_bar_phi[phi] * fc_phi[phi]
                                    for phi in F_t) ,f"MILP-24"

    return delta_F_u_mu_phi, delta_F_v_phi, fc_e_phi , rank, n_G, n_int, fc_phi, ac_lf_gamma, MILP

# -------- 3.4 Descriptors for the Number of Specified Degree --------
def descriptprs_for_the_number_of_specified_degree(
        #Constants
        deg_H_bar_phi, #dict {phi:deg of root of phi}
        Vo,
        Vo_bar,
        Eo,#list,Eo[ring-edge]
        Eo_u,#dict
        Eo_bar_u,#dict {u: non ring edge[u]}
        Xi_u,
        x_u_ksi,
        E_Cu_mu,
        N_T_v,#dict,{v:N_T(v)}
        F_u,
        xEdge_u_e_v,
        xNode_u_ePrime_mu,
        delta_F_v_phi,
        delta_F_u_mu_phi,
        #model
        MILP
):
    # -----prepare_variables-----#

    #Define delta_Deg(u,[mu];d)
    delta_Deg_u_mu_d={(u,mu,d):pulp.LpVariable(f"delta_Deg({u},{mu},{d})",cat=pulp.LpBinary)
                     for u in Vo for mu in range(1,1+Cmax) for d in range(1,1+max_degree)}

    #Define delta_Deg(v;d)
    delta_Deg_v_d={(v,d): pulp.LpVariable(f"delata_Deg({v},{d})",cat=pulp.LpBinary)
                   for v in Vo_bar for d in range(1, 1 + max_degree)}

    #Define deltaInt_Deg(u,[mu];d)
    deltaInt_Deg_u_mu_d={(u,mu,d):pulp.LpVariable(f"deltaInt_Deg({u},{mu},{d})",cat=pulp.LpBinary)
                        for u in Vo for mu in range(1,1+Cmax) for d in range(1 , 1+max_degree)}

    #Define deltaIntDeg(v;d)
    deltaInt_Deg_v_d={(v,d):pulp.LpVariable(f"deltaInt_Deg({v},{d}",cat=pulp.LpBinary)
                      for v in Vo_bar for d in range(1, 1 + max_degree)}

    #Define deg(d)
    deg_d={d:pulp.LpVariable(f"deg({d})",cat=pulp.LpInteger)
           for d in range(1,1+max_degree)}

    #Define degInt(d)
    degInt_d={d:pulp.LpVariable(f"degInt({d})",cat=pulp.LpInteger)
           for d in range(1,1+max_degree)}

    #Define deg([e];d)
    deg_e_d={(e,d):pulp.LpVariable(f"deg({e},{d})",0,2,cat=pulp.LpInteger)
             for e in Eo for d in range(1,max_degree+1)}

    #Define degInt([e];d)
    degInt_e_d={(e,d):pulp.LpVariable(f"degInt({e},{d})",0,2,cat=pulp.LpInteger)
             for e in Eo for d in range(1,max_degree+1)}
    #Define degEdge_u_e_mu_plus
    degEdge_u_e_mu_plus = {(u, e, mu) : pulp.LpVariable(f"degEdge_{u}_{e}_{mu}_puls",cat=pulp.LpInteger)
                           for u in Vo for mu in range(1, 1 + Cmax) for e in Eo_u[u]}

    # Define degEdge_u_e_mu_minus
    degEdge_u_e_mu_minus = {(u, e, mu): pulp.LpVariable(f"degEdge_{u}_{e}_{mu}_minus", cat=pulp.LpInteger)
                           for u in Vo for mu in range(1, 1 + Cmax) for e in Eo_u[u]}

    # -----add_constraints-----#

    # -------- Constraint (25) --------
    for u in Vo:
        for mu in range(1, 1 + Cmax):
            for e in Eo_u[u]:
                MILP += 0 <= degEdge_u_e_mu_plus[(u,e,mu)] ,f"MILP-25_{u}_{e}_{mu}-1"
                MILP += degEdge_u_e_mu_plus[(u,e,mu)] <= 4 * pulp.lpSum(xEdge_u_e_v[(u,e,v)] for v in E_Cu_mu[mu]),f"MILP-25_{u}_{e}_{mu}-2"
                MILP += 0 <= degEdge_u_e_mu_minus[(u, e, mu)], f"MILP-25_{u}_{e}_{mu}-3"
                MILP += degEdge_u_e_mu_minus[(u, e, mu)] <= 4 * pulp.lpSum(xEdge_u_e_v[(u, e, v)] for v in E_Cu_mu[mu]), f"MILP-25_{u}_{e}_{mu}-4"

    # -------- Constraint (26) --------
                for v in E_Cu_mu[mu]:
                    MILP += 4 * (xEdge_u_e_v[(u,e,v)] - 1) <= degEdge_u_e_mu_plus[(u,e,mu)] - \
                        1 - pulp.lpSum(xNode_u_ePrime_mu[(u,ep,mu)] for ep in Eo_bar_u[u]),f"MILP-26_{u}_{mu}_{e}_{v}-1"
                    MILP += degEdge_u_e_mu_plus[(u,e,mu)] - \
                        1 - pulp.lpSum(xNode_u_ePrime_mu[(u,ep,mu)] for ep in Eo_bar_u[u]) <= 4 * (1 - xEdge_u_e_v[(u,e,v)]),f"MILP-26_{u}_{mu}_{e}_{v}-2"

    # -------- Constraint (27) --------
    # u = e[0] ,up = e[1] ,i1 = E_cu[v][0] , i2 = E_cu[v][1] , j1 = E_cu[vp][0], j2 = E_cu[vp][1]
    for e in Eo:
        for v in range(1, 2 * Cmax - Cmin + 1):
            for vp in range(1, 2 * Cmax - Cmin + 1):
                MILP += 3 * (xEdge_u_e_v[(e[0],e,v)] + xEdge_u_e_v[(e[1],e,vp)] - 2) <= \
                        degEdge_u_e_mu_minus[(e[0], e, E_Cu[v][0])] - degEdge_u_e_mu_plus[(e[1], e, E_Cu[vp][1])],f"MILP-27_{e}_{v}_{vp}-1"
                MILP += degEdge_u_e_mu_minus[(e[0], e, E_Cu[v][0])] - degEdge_u_e_mu_plus[(e[1], e, E_Cu[vp][1])] <=\
                        3 * (2 - xEdge_u_e_v[(e[0],e,v)] - xEdge_u_e_v[(e[1],e,vp)]),f"MILP-27_{e}_{v}_{vp}-2"

                MILP += 3 * (xEdge_u_e_v[(e[0],e,v)] + xEdge_u_e_v[(e[1],e,vp)] - 2) <= \
                        degEdge_u_e_mu_minus[(e[0], e, E_Cu[v][1])] - degEdge_u_e_mu_plus[(e[1], e, E_Cu[vp][0])],f"MILP-27_{e}_{v}_{vp}-3"
                MILP += degEdge_u_e_mu_minus[(e[0], e, E_Cu[v][1])] - degEdge_u_e_mu_plus[(e[1], e, E_Cu[vp][0])] <= \
                        3 * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vp)]), f"MILP-27_{e}_{v}_{vp}-4"

                MILP += 3 * (xEdge_u_e_v[(e[0], e, v)] + xEdge_u_e_v[(e[1], e, vp)] - 2) <= \
                        degEdge_u_e_mu_minus[(e[1], e, E_Cu[vp][1])] - degEdge_u_e_mu_plus[(e[0], e, E_Cu[v][0])], f"MILP-27_{e}_{v}_{vp}-5"
                MILP += degEdge_u_e_mu_minus[(e[1], e, E_Cu[vp][1])] - degEdge_u_e_mu_plus[(e[0], e, E_Cu[v][0])] <= \
                        3 * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vp)]), f"MILP-27_{e}_{v}_{vp}-6"

                MILP += 3 * (xEdge_u_e_v[(e[0],e,v)] + xEdge_u_e_v[(e[1],e,vp)] - 2) <= \
                        degEdge_u_e_mu_minus[(e[1], e, E_Cu[vp][0])] - degEdge_u_e_mu_plus[(e[0], e, E_Cu[v][1])],f"MILP-27_{e}_{v}_{vp}-7"
                MILP += degEdge_u_e_mu_minus[(e[1], e, E_Cu[vp][0])] - degEdge_u_e_mu_plus[(e[0], e, E_Cu[v][1])] <= \
                        3 * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vp)]),f"MILP-27_{e}_{v}_{vp}-8"
    # -------- Constraint (28) --------
    for u in Vo:
        for mu in range(1,1+Cmax):
            MILP += pulp.lpSum(delta_Deg_u_mu_d[(u,mu,d)] for d in range(1,1+max_degree)) == pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) >= mu),f"MILP-28_{u}_{mu}"
    # -------- Constraint (29) --------
    for u in Vo:
        for mu in range(1,1+Cmax):
            MILP += 2 * (pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) >= mu) - 1) <= \
                    pulp.lpSum(d*delta_Deg_u_mu_d[(u,mu,d)] for d in range(1,1 + max_degree)) - \
                    (2 + pulp.lpSum(deg_H_bar_phi[phi] * delta_F_u_mu_phi[(u, mu, phi)] for phi in F_u[u]) +
                    pulp.lpSum(xNode_u_ePrime_mu[(u,ePrime,mu)] for ePrime in Eo_bar_u[u]) +
                    pulp.lpSum(degEdge_u_e_mu_minus[(u,e,mu)] for e in Eo_u[u]) ),f"MILP-29_{u}_{mu}-1"

            MILP += pulp.lpSum(d*delta_Deg_u_mu_d[(u,mu,d)] for d in range(1,1 + max_degree)) <= \
                    2 + pulp.lpSum(deg_H_bar_phi[phi] * delta_F_u_mu_phi[(u, mu, phi)] for phi in F_u[u]) + \
                    pulp.lpSum(xNode_u_ePrime_mu[(u, ePrime, mu)] for ePrime in Eo_bar_u[u]) + \
                    pulp.lpSum(degEdge_u_e_mu_minus[(u,e,mu)] for e in Eo_u[u]), f"MILP-29_{u}_{mu}-2"
    # -------- Constraint (30) --------
    for u in Vo:
        for mu in range(1, 1 + Cmax):
            MILP += pulp.lpSum(deltaInt_Deg_u_mu_d[(u,mu,d)] for d in range(1,1+max_degree)) == \
                    pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) >= mu),f"MILP-30_{u}_{mu}"
    # -------- Constraint (31) --------
    for u in Vo:
        for mu in range(1, 1 + Cmax):
            MILP += 2 * ( pulp.lpSum( x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) >= mu ) - 1) <= \
                    pulp.lpSum(d * deltaInt_Deg_u_mu_d[(u,mu,d)] for d in range(1,1 + max_degree)) - \
                    (2 + pulp.lpSum(xNode_u_ePrime_mu[(u,ePrime,mu)] for ePrime in Eo_bar_u[u]) +
                     pulp.lpSum(degEdge_u_e_mu_minus[(u,e,mu)] for e in Eo_u[u] ) ) , f"MILP-31_{u}_{mu}-1"

            MILP += pulp.lpSum(d * deltaInt_Deg_u_mu_d[(u, mu, d)] for d in range(1, 1 + max_degree)) <= \
                2 + pulp.lpSum(xNode_u_ePrime_mu[(u, ePrime, mu)] for ePrime in Eo_bar_u[u]) + \
                 pulp.lpSum( degEdge_u_e_mu_minus[(u,e,mu)] for e in Eo_u[u] ) , f"MILP-31_{u}_{mu}-2"

    # -------- Constraint (32) --------
    for v in Vo_bar:
        MILP += pulp.lpSum(delta_Deg_v_d[(v,d)] for d in range(1,1+max_degree)) == 1,f"MILP-32_{v}"
    # -------- Constraint (33) --------
    for v in Vo_bar:
        MILP += (pulp.lpSum(d * delta_Deg_v_d[(v, d)] for d in range(1, 1 + max_degree)) ==
                 len(N_T_v[v]) + pulp.lpSum(deg_H_bar_phi[phi] * delta_F_v_phi[(v, phi)] for phi in F_u[v])),f"MILP-33_{v}"
    # -------- Constraint (34) --------
    for v in Vo_bar:
        MILP += pulp.lpSum(deltaInt_Deg_v_d[(v,d)] for d in range(1,1+max_degree)) ==1,f"MILP-34_{v}"
    # -------- Constraint (35) --------
    for v in Vo_bar:
        MILP += pulp.lpSum(d*deltaInt_Deg_v_d[(v,d)] for d in range(1,1+max_degree)) == len(N_T_v[v]),f"MILP-35_{v}"
    # -------- Constraint (36) --------
    for e in Eo:
        for v in range(1,1+2*Cmax-Cmin):
            for d in range(1,max_degree+1):
                MILP += deg_e_d[(e,d)] -delta_Deg_u_mu_d[(e[0],E_Cu[v][0],d)] - delta_Deg_u_mu_d[(e[0],E_Cu[v][1],d)] <=\
                        2*(1-xEdge_u_e_v[(e[0],e,v)]),f"MILP-36_{e}_{v}_{d}-1"

                MILP += deg_e_d[(e,d)] -delta_Deg_u_mu_d[(e[0],E_Cu[v][0],d)] - delta_Deg_u_mu_d[(e[0],E_Cu[v][1],d)] >=\
                        2*(xEdge_u_e_v[(e[0],e,v)]-1),f"MILP-36_{e}_{v}_{d}-2"
    # -------- Constraint (37) --------
    for e in Eo:
        for v in range(1,1+2*Cmax-Cmin):
            for d in range(1,max_degree+1):
                MILP += degInt_e_d[(e,d)] - deltaInt_Deg_u_mu_d[(e[0], E_Cu[v][0], d)] - deltaInt_Deg_u_mu_d[(e[0], E_Cu[v][1],d)] <= \
                        2*(1-xEdge_u_e_v[(e[0],e,v)]),f"MILP-37_{e}_{v}_{d}-1"
                MILP += degInt_e_d[(e,d)] - deltaInt_Deg_u_mu_d[(e[0], E_Cu[v][0], d)] - deltaInt_Deg_u_mu_d[(e[0], E_Cu[v][1],d)] >= \
                        2 * (xEdge_u_e_v[(e[0], e, v)] - 1),f"MILP-37_{e}_{v}_{d}-2"
    # -------- Constraint (38) --------
    for d in range(1, 1+max_degree):
        MILP += (deg_d[d] == pulp.lpSum(delta_Deg_u_mu_d[(u,mu,d)] for u in Vo for mu in range(1,1+Cmax)) +
                 pulp.lpSum(delta_Deg_v_d[(v,d)] for v in Vo_bar) - pulp.lpSum(deg_e_d[(e, d)] for e in Eo)),f"MILP-38_{d}"
    # -------- Constraint (39) --------
    for d in range(1, 1 + max_degree):
        MILP += (degInt_d[d] == pulp.lpSum(deltaInt_Deg_u_mu_d[(u,mu,d)] for u in Vo for mu in range(1,1+Cmax)) +
                 pulp.lpSum(deltaInt_Deg_v_d[(v,d)] for v in Vo_bar) - pulp.lpSum(degInt_e_d[(e, d)] for e in Eo)),f"MILP-39_{d}"

    return  (delta_Deg_u_mu_d, delta_Deg_v_d, deltaInt_Deg_u_mu_d,
             deltaInt_Deg_v_d, deg_d, degInt_d, deg_e_d, degInt_e_d,
             degEdge_u_e_mu_minus, degEdge_u_e_mu_plus, MILP)
    # -------- 3.5 Assigning Bond-Multiplicity --------
def assigning_bond_multiplicity(
        Vo,
        E_T,
        Eo,
        Eo_u,
        Eo_bar, #list, Non ring edge
        xEdge_u_e_v,
        e_u_i,
        MILP

):
    # -----prepare_variables-----#

    #Define beta_u_i
    beta_u_i={(u,i):pulp.LpVariable(f"beta_{u}_{i}",0,3,cat=pulp.LpInteger)
              for u in Vo for i in range(1,2*Cmax-Cmin+1) }

    #Define beta_[e]
    beta_e={e:pulp.LpVariable(f"beta_{e}",1,3,cat=pulp.LpInteger)
            for e in E_T}

    #Define delta_beta(u,i;m)
    delta_beta_u_i_m={(u,i,m):pulp.LpVariable(f"delta_beta({u},{i},{m})",cat=pulp.LpBinary)
                     for u in Vo for i in range(1,2*Cmax-Cmin+1) for m in range(1,4)}

    #Define delta_beta([e];m)
    delta_beta_e_m={(e,m):pulp.LpVariable(f"delta_beta({e},{m})",cat=pulp.LpBinary)
                   for e in E_T for m in range(1,4)}

    #Define bd(m)
    bd_m={m:pulp.LpVariable(f"bd({m})",cat=pulp.LpInteger)
          for m in range(1,4)}

    # -----add_constraints-----#

    # -------- Constraint (40) --------
    for u in Vo:
        for i in range(1,1+2*Cmax-Cmin):
            MILP += beta_u_i[(u,i)] >=e_u_i[(u,i)],f"MILP-40_{u}_{i}-1"
            MILP += beta_u_i[(u,i)] <= 3*e_u_i[(u,i)],f"MILP-40_{u}_{i}-2"
    # -------- Constraint (41) --------
    for e in E_T:
        MILP += beta_e[e] <= 3,f"MILP-41_{e}-1"
        MILP += beta_e[e] >= 1,f"MILP-41_{e}-2"
    # -------- Constraint (42) --------
    for u in Vo:
        for i in range(1,1+2*Cmax-Cmin):
            MILP += pulp.lpSum(delta_beta_u_i_m[(u,i,m)] for m in range(1,4)) == e_u_i[(u,i)],f"MILP-42_{u}_{i}-1"
            MILP += pulp.lpSum(m*delta_beta_u_i_m[(u,i,m)] for m in range(1,4)) == beta_u_i[(u,i)],f"MILP-42_{u}_{i}-2"
    # -------- Constraint (43) --------
    for e in E_T:
        MILP += pulp.lpSum(delta_beta_e_m[(e,m)] for m in range(1,4)) == 1,f"MILP-43_{e}-1"
        MILP += pulp.lpSum(m*delta_beta_e_m[(e,m)] for m in range(1,4)) == beta_e[e],f"MILP-43_{e}-2"
    # -------- Constraint (44) --------
    for u in Vo:
        for e in Eo_u[u]:
            for i in range(1,1+2*Cmax-Cmin):
                MILP += 3*(xEdge_u_e_v[(u,e,i)]-1) <= beta_u_i[(u,i)] - beta_e[e],f"MILP-44_{u}_{e}_{i}-1"
                MILP += 3*(1-xEdge_u_e_v[(u,e,i)]) >= beta_u_i[(u,i)] - beta_e[e],f"MILP-44_{u}_{e}_{i}-2"
    # -------- Constraint (45) --------
    for m in range(1,4):
        MILP += bd_m[m] == pulp.lpSum(delta_beta_u_i_m[(u,i,m)] for u in Vo for i in range(1,1+2*Cmax-Cmin)) +\
                 pulp.lpSum(delta_beta_e_m[(eP,m)] for eP in Eo_bar) -\
                 pulp.lpSum(delta_beta_e_m[(e,m)] for e in Eo),f"MILP-45_{m}-1"

    return beta_u_i, beta_e, delta_beta_u_i_m, delta_beta_e_m, bd_m, MILP

   # -------- 3.6 Assigning Chemical Elements and Valence Condition --------
def assigning_chemical_elements_and_valence_condition(
        #Constants
        Lambda,  #dict,A[element:code]
        alpha_r_phi, #dict
        valEx_F_phi,#dict
        eledeg_F_phi,#dict
        n_a_phi,# dict
        val_a,#dict
        mass_star_a,#dict
        M_ms, #large number
        naInt_LB_a,#dict,{a:xx}
        naInt_UB_a,#dict,{a:xx}
        naEx_LB_a,#dict
        naEx_UB_a,#dict
        n_LB,
        n_UB,
        na_LB_a,
        na_UB_a,
        F_T,
        F_u,
        Xi_u,
        Vo,
        Vo_bar,
        Eo,#list,Eo[ring edge]
        Eo_u,#dict,Eo_u{u:ring edge incident to u}
        Eo_bar_u,#dict,Eo_line_u{u:non ring edge incident to u}
        delta_F_u_mu_phi,
        delta_F_v_phi,
        xNode_u_ePrime_mu,
        xEdge_u_e_v,
        x_u_ksi,
        beta_e,
        beta_u_i,
        fc_e_phi,
        fc_phi,
        n_G,
        MILP
        ):
    # -----prepare_variables-----#

    #Define alpha(u,[mu])
    alpha_u_mu={(u,mu):pulp.LpVariable(f"alpha({u},{mu})",cat=pulp.LpInteger)
                for u in Vo for mu in range(1,Cmax+1)}

    #Define alpha(v)
    alpha_v={v:pulp.LpVariable(f"alpha({v})",cat=pulp.LpInteger)
             for v in Vo_bar}

    #Define delta_alpha(u,[mu];[a])
    delta_alpha_u_mu_a={(u,mu,a):pulp.LpVariable(f"delta_alpha({u},{mu},{a})",cat=pulp.LpBinary)
                       for u in Vo for mu in range(1,Cmax+1) for a in Lambda}

    #Define delta_alpha(v;[a])
    delta_alpha_v_a={(v,a):pulp.LpVariable(f"delta_alpha({v},{a})",cat=pulp.LpBinary)
                     for v in Vo_bar for a in Lambda}

    #Define betaNode_u_ePrime_mu
    betaNode_u_ePrime_mu={(u,ePrime,mu):pulp.LpVariable(f"betaNode_{u}_{ePrime}_{mu}",cat=pulp.LpInteger)
                          for u in Vo for ePrime in Eo_bar_u[u] for mu in range(1, Cmax + 1)}

    #Define betaEdge_u_e_mu_+
    betaEdge_u_e_mu_plus={(u,e,mu):pulp.LpVariable(f"betaEdge_{u}_{e}_{mu}_plus",cat=pulp.LpInteger)
                          for u in Vo for e in Eo_u[u] for mu in range(1,Cmax+1)}

    #Define betaEdge_u_e_mu_-
    betaEdge_u_e_mu_mins={(u,e,mu):pulp.LpVariable(f"betaEdge_{u}_{e}_{mu}_mins",cat=pulp.LpInteger)
                           for u in Vo for e in Eo_u[u] for mu in range(1,Cmax+1)}

    #Define na([e];[a])
    na_e_a={(e,a):pulp.LpVariable(f"na({e},{a})",0,2,cat=pulp.LpInteger)
            for e in Eo for a in Lambda}

    #Define naInt([a]) 不会有H吗？如果没有H还需要改后面的约束部分
    naInt_a={a:pulp.LpVariable(f"naInt({a})",lowBound=naInt_LB_a[a],upBound=naInt_UB_a[a],cat=pulp.LpInteger)
             for a in Lambda}
    #Define naEx([a])
    naEx_a = {a: pulp.LpVariable(f"naEx({a})", lowBound=naEx_LB_a[a], upBound=naEx_UB_a[a], cat=pulp.LpInteger)
               for a in Lambda}
    #Define na([a])
    na_a = {a:pulp.LpVariable(f"na({a})",lowBound=na_LB_a[a],upBound=na_UB_a[a])
                for a in Lambda}
    #Define delta_atm(i)
    delta_atm_i={i:pulp.LpVariable(f"delta_atm({i})",cat=pulp.LpBinary)
                 for i in range(n_LB+na_LB_a["H1"],1+n_UB+na_UB_a["H1"])}
    #Define Mass
    Mass=pulp.LpVariable("Mass",0,cat=pulp.LpInteger)

    #Define msLine
    msLine=pulp.LpVariable("msLine",0,cat=pulp.LpContinuous)

    # -----add_constraints-----#

    # -------- Constraint (46) --------
    for u in Vo:
        for mu in range(1,Cmax+1):
            MILP += alpha_u_mu[(u,mu)] == pulp.lpSum(Lambda[alpha_r_phi[phi]] * delta_F_u_mu_phi[(u,mu,phi)] for phi in F_u[u]),f"MILP-46_{u}_{mu}"
    # -------- Constraint (47) --------
    for u in Vo:
        for mu in range(1, Cmax + 1):
            MILP += pulp.lpSum(delta_alpha_u_mu_a[(u,mu,a)] for a in Lambda) == pulp.lpSum(x_u_ksi[(u,ksi)] for ksi in Xi_u[u] if len(ksi) >= mu),f"MILP-47_{u}_{mu}"
    # -------- Constraint (48) --------
    for u in Vo:
        for mu in range(1, Cmax + 1):
            MILP += pulp.lpSum(Lambda[a]*delta_alpha_u_mu_a[(u,mu,a)] for a in Lambda) == alpha_u_mu[(u,mu)],f"MILP-48_{u}_{mu}"
    # -------- Constraint (49) --------
    for v in Vo_bar:
        MILP += alpha_v[v] == pulp.lpSum(Lambda[alpha_r_phi[phi]]*delta_F_v_phi[(v,phi)] for phi in F_u[v]),f"MILP-49_{v}"
    # -------- Constraint (50) --------
    for v in Vo_bar:
        MILP += pulp.lpSum(delta_alpha_v_a[(v,a)] for a in Lambda) == 1,f"MILP-50_{v}"
    # -------- Constraint (51) --------
    for v in Vo_bar:
        MILP += pulp.lpSum(Lambda[a]* delta_alpha_v_a[(v,a)] for a in Lambda) == alpha_v[v],f"MILP-51_{v}"
    # -------- Constraint (52) --------
    for u in Vo:
        for mu in range(1,Cmax+1):
            for ep in Eo_bar_u[u]:
                MILP += betaNode_u_ePrime_mu[(u,ep,mu)] <= 3 * xNode_u_ePrime_mu[(u,ep,mu)],f"MILP-52_{u}_{mu}_{ep}-1"
                MILP += betaNode_u_ePrime_mu[(u,ep,mu)] >= 0,f"MULP-49_{u}_{mu}_{ep}-2"
    # -------- Constraint (53) --------
    for u in Vo:
        for mu in range(1, Cmax + 1):
            for ep in Eo_bar_u[u]:
                MILP += 3*(xNode_u_ePrime_mu[(u,ep,mu)]-1) <= beta_e[ep]-betaNode_u_ePrime_mu[(u,ep,mu)],f"MILP-53_{u}_{mu}_{ep}-1"
                MILP += beta_e[ep]-betaNode_u_ePrime_mu[(u,ep,mu)] <= 3 * (1-xNode_u_ePrime_mu[(u,ep,mu)]),f"MILP-53_{u}_{mu}_{ep}-2"
    # -------- Constraint (54) --------
    for u in Vo:
        for mu in range(1, Cmax + 1):
            for e in Eo_u[u]:
                MILP += betaEdge_u_e_mu_plus[(u, e, mu)] <= 6 * pulp.lpSum(xEdge_u_e_v[(u,e,v)] for v in E_Cu_mu[mu]),f"MILP-54_{u}_{mu}_{e}-1"
                MILP += betaEdge_u_e_mu_plus[(u, e, mu)] >= 0, f"MILP-54_{u}_{mu}_{e}-2"
                MILP += betaEdge_u_e_mu_mins[(u, e, mu)] <= 6 * pulp.lpSum(xEdge_u_e_v[(u, e, v)] for v in E_Cu_mu[mu]), f"MILP-54_{u}_{mu}_{e}-3"
                MILP += betaEdge_u_e_mu_mins[(u, e, mu)] >= 0, f"MILP-54_{u}_{mu}_{e}-4"
    # -------- Constraint (55) --------
    for u in Vo:
        for mu in range(1, Cmax + 1):
            for e in Eo_u[u]:
                for v in E_Cu_mu[mu]:
                    MILP += 6 * (xEdge_u_e_v[(u,e,v)] - 1) <= \
                            betaEdge_u_e_mu_plus[(u,e,mu)] + beta_e[e] - pulp.lpSum(beta_u_i[(u,vP)] for vP in E_Cu_mu[mu]) - pulp.lpSum(betaNode_u_ePrime_mu[(u,ep,mu)] for ep in Eo_bar_u[u]),f"MILP-55_{u}_{mu}_{e}_{v}-1"
                    MILP += betaEdge_u_e_mu_plus[(u,e,mu)] + beta_e[e] - pulp.lpSum(beta_u_i[(u,vP)] for vP in E_Cu_mu[mu]) - pulp.lpSum(betaNode_u_ePrime_mu[(u,ep,mu)] for ep in Eo_bar_u[u]) <= \
                            6 * (1 - xEdge_u_e_v[(u,e,v)]),f"MILP-55_{u}_{mu}_{e}_{v}-2"
    # -------- Constraint (56) --------
    for e in Eo:
        for v in range(1,1+2*Cmax-Cmin): #u = e[0] u'=e[1],E_Cu[v][0]=i1,  E_Cu[v][1]=i2   E_Cu[vPrime][0]=j1 E_Cu[vPrime][1] =j2
            for vPrime in range(1, 1 + 2 * Cmax - Cmin):
                MILP += 3 * (xEdge_u_e_v[(e[0],e,v)] + xEdge_u_e_v[(e[1],e,vPrime)] - 2) <= betaEdge_u_e_mu_mins[(e[0],e,E_Cu[v][0])] -betaEdge_u_e_mu_plus[(e[1],e,E_Cu[vPrime][1])],f"MILP-56_{e}_{v}_{vPrime}-1"
                MILP += betaEdge_u_e_mu_mins[(e[0], e, E_Cu[v][0])] - betaEdge_u_e_mu_plus[(e[1], e, E_Cu[vPrime][1])] <= 3 * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vPrime)]),f"MILP-56_{e}_{v}_{vPrime}-2"

                MILP += 3 * (xEdge_u_e_v[(e[0], e, v)] + xEdge_u_e_v[(e[1], e, vPrime)] - 2) <= betaEdge_u_e_mu_mins[(e[0], e, E_Cu[v][1])] - betaEdge_u_e_mu_plus[(e[1], e, E_Cu[vPrime][0])],f"MILP-56_{e}_{v}_{vPrime}-3"
                MILP += betaEdge_u_e_mu_mins[(e[0], e, E_Cu[v][1])] - betaEdge_u_e_mu_plus[(e[1], e, E_Cu[vPrime][0])] <= 3 * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vPrime)]),f"MILP-56_{e}_{v}_{vPrime}-4"

                MILP += 3 * (xEdge_u_e_v[(e[0], e, v)] + xEdge_u_e_v[(e[1], e, vPrime)] - 2) <= betaEdge_u_e_mu_mins[(e[1], e, E_Cu[vPrime][1])] - betaEdge_u_e_mu_plus[(e[0], e, E_Cu[v][0])],f"MILP-56_{e}_{v}_{vPrime}-5"
                MILP += betaEdge_u_e_mu_mins[(e[1], e, E_Cu[vPrime][1])] - betaEdge_u_e_mu_plus[(e[0], e, E_Cu[v][0])] <= 3 * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vPrime)]),f"MILP-56_{e}_{v}_{vPrime}-6"

                MILP += 3 * (xEdge_u_e_v[(e[0], e, v)] + xEdge_u_e_v[(e[1], e, vPrime)] - 2) <= betaEdge_u_e_mu_mins[(e[1], e, E_Cu[vPrime][0])] - betaEdge_u_e_mu_plus[(e[0], e, E_Cu[v][1])],f"MILP-56_{e}_{v}_{vPrime}-7"
                MILP += betaEdge_u_e_mu_mins[(e[1], e, E_Cu[vPrime][0])] - betaEdge_u_e_mu_plus[(e[0], e, E_Cu[v][1])] <= 3 * (2 - xEdge_u_e_v[(e[0], e, v)] - xEdge_u_e_v[(e[1], e, vPrime)]),f"MILP-56_{e}_{v}_{vPrime}-8"
    # -------- Constraint (57) --------
    for u in Vo:
        for mu in range(1,Cmax+1):
            MILP += (pulp.lpSum(val_a[a]*delta_alpha_u_mu_a[(u,mu,a)] for a in Lambda) ==
                     pulp.lpSum(beta_u_i[u,v] for v in E_Cu_mu[mu]) +
                     pulp.lpSum(betaNode_u_ePrime_mu[(u,ep,mu)] for ep in Eo_bar_u[u]) +
                     pulp.lpSum(betaEdge_u_e_mu_mins[(u,e,mu)] for e in Eo_u[u]) +
                     pulp.lpSum((valEx_F_phi[phi]-eledeg_F_phi[phi])* delta_F_u_mu_phi[(u,mu,phi)] for phi in F_u[u])),f"MILP-57_{u}_{mu}"
    # -------- Constraint (58) --------
    for v in Vo_bar:
        MILP += (pulp.lpSum(val_a[a]*delta_alpha_v_a[(v,a)] for a in Lambda) == pulp.lpSum(beta_e[eP] for eP in Eo_bar_u[v]) +
                 pulp.lpSum((valEx_F_phi[phi]-eledeg_F_phi[phi])*delta_F_v_phi[(v,phi)] for phi in F_u[v] )),f"MILP-58_{v}"
    # -------- Constraint (59) --------
    for a in Lambda:
        for e in Eo:
            MILP += na_e_a[(e,a)] == pulp.lpSum(fc_e_phi[(e,phi)] for phi in F_T if alpha_r_phi[phi] == a),f"MILP-59_{a}_{e}"
    # -------- Constraint (60) --------
    for a in Lambda:
        MILP += (naInt_a[a] == pulp.lpSum(delta_alpha_u_mu_a[(u, mu, a)] for u in Vo for mu in range(1,1+Cmax)) +
                 pulp.lpSum(delta_alpha_v_a[(v, a)] for v in Vo_bar) -
                 pulp.lpSum(na_e_a[(e, a)] for e in Eo)),f"MILP-60_{a}"
    # -------- Constraint (61) --------
    for a in Lambda:
        MILP += naEx_a[a] == pulp.lpSum(n_a_phi[(a,phi)] * fc_phi[phi] for phi in F_T),f"MILP-61_{a}"
    # -------- Constraint (62) --------
    for a in Lambda:
        MILP += na_a[a] == naInt_a[a] + naEx_a[a],f"MILP-62_{a}"
    # -------- Constraint (63) --------
    MILP += Mass == pulp.lpSum(mass_star_a[a] * na_a[a] for a in Lambda),f"MILP-63"
    # -------- Constraint (64) --------
    MILP += pulp.lpSum(delta_atm_i[i] for i in range(n_LB+na_LB_a["H1"],n_UB+na_UB_a["H1"]+1)) == 1,f"MILP-64"
    # -------- Constraint (65) --------
    MILP += pulp.lpSum(i*delta_atm_i[i] for i in range(n_LB+na_LB_a["H1"],n_UB+na_UB_a["H1"]+1)) == n_G+naEx_a["H1"],f"MILP-65"
    # -------- Constraint (66) --------
    for i in range(n_LB + na_LB_a["H1"], n_UB + na_UB_a["H1"] + 1):
        MILP += M_ms*(delta_atm_i[i] - 1) <= msLine - Mass * 1/i, f"MILP-66_{i}-1"
        MILP += M_ms*(1 - delta_atm_i[i]) >= msLine - Mass * 1/i, f"MILP-66_{i}-2"

    return (alpha_u_mu, alpha_v, delta_alpha_u_mu_a, delta_alpha_v_a, betaNode_u_ePrime_mu, betaEdge_u_e_mu_plus,
              betaEdge_u_e_mu_mins, na_e_a, naInt_a, naEx_a, na_a, delta_atm_i, Mass, msLine, MILP)
   # -------- 3.7 Descriptors for the Number of Adjacency-configurations --------
def desriptors_for_the_number_of_ac(
        Lambda,
        gammaInt_ac,#dict {ac : code }
        acInt_LB, #dict {ac : acInt_LB}
        acInt_UB, #dict {ac : acInt_UB}
        Vo,
        E_T,
        xNode_u_ePrime_mu,
        xEdge_u_e_v,
        e_u_i,
        beta_u_i,
        alpha_u_mu,
        alpha_v,
        beta_e,
        Eo_bar,
        Eo,
        MILP
):
    #Define delta_ac(u,[v];[gamma])
    delta_ac_u_v_gamma={(u,v,gamma):pulp.LpVariable(f"delta_ac({u},{v},{gamma})",cat=pulp.LpBinary)
                       for u in Vo for v in range(1,2*Cmax-Cmin+1) for gamma in gammaInt_ac}

    # Define delta_ac([e];[gamma])
    delta_ac_e_gamma = {(e, gamma): pulp.LpVariable(f"delta_ac({e},{gamma})", cat=pulp.LpBinary)
                         for e in E_T for gamma in gammaInt_ac}

    #Define acInt([gamma])
    acInt_gamma={gamma:pulp.LpVariable(f"acInt({gamma})",lowBound=acInt_LB[gamma],upBound=acInt_UB[gamma],cat=pulp.LpInteger)
                 for gamma in gammaInt_ac}

    # -----add_constraints-----#

    # -------- Constraint (67) --------
    for u in Vo:
        for v in range(1, 2 * Cmax - Cmin + 1):
            MILP += pulp.lpSum(delta_ac_u_v_gamma[(u,v,gamma)] for gamma in gammaInt_ac) == e_u_i[(u,v)],f"MILP-67_{u}_{v}"
    # -------- Constraint (68) --------
    for u in Vo:
        for v in range(1, 2 * Cmax - Cmin + 1):
            MILP += pulp.lpSum(gamma[2]*delta_ac_u_v_gamma[(u,v,gamma)] for gamma in gammaInt_ac) - beta_u_i[(u,v)] >= 3 * (e_u_i[(u,v)] - 1),f"MILP-68_{u}_{v}-1"

            MILP += pulp.lpSum(gamma[2] * delta_ac_u_v_gamma[(u, v, gamma)] for gamma in gammaInt_ac) <= beta_u_i[(u, v)],f"MILP-68_{u}_{v}-2"
    # -------- Constraint (69) --------
    for u in Vo:
        for v in range(1, 2 * Cmax - Cmin + 1):
            MILP += pulp.lpSum(Lambda[gamma[0]] * delta_ac_u_v_gamma[(u, v, gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(u, E_Cu[v][0])] >= len(Lambda) * (e_u_i[(u, v)] - 1),f"MILP-69_{u}_{v}-1"
            MILP += pulp.lpSum(Lambda[gamma[0]] * delta_ac_u_v_gamma[(u, v, gamma)] for gamma in gammaInt_ac) <= alpha_u_mu[(u, E_Cu[v][0])],f"MILP-69_{u}_{v}-2"
    # -------- Constraint (70) --------
            MILP += pulp.lpSum(Lambda[gamma[1]] * delta_ac_u_v_gamma[(u, v, gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(u, E_Cu[v][1])] >= len(Lambda) * (e_u_i[(u, v)] - 1), f"MILP-70_{u}_{v}-1"
            MILP += pulp.lpSum(Lambda[gamma[1]] * delta_ac_u_v_gamma[(u, v, gamma)] for gamma in gammaInt_ac) <= alpha_u_mu[(u, E_Cu[v][1])],f"MILP-70_{u}_{v}-2"
    # -------- Constraint (71) --------
    for e in E_T:
        MILP += pulp.lpSum(delta_ac_e_gamma[(e,gamma)] for gamma in gammaInt_ac) == 1,f"MILP-71_{e}"
    # -------- Constraint (72) --------
    for e in E_T:
        MILP += pulp.lpSum(gamma[2]* delta_ac_e_gamma[(e,gamma)] for gamma in gammaInt_ac) == beta_e[e],f"MILP-72_{e}"
    # -------- Constraint (73) --------
    for ePrime in Eo_bar: #normally recognize that e'=uv in E(T) satisfied u < v
        if ePrime[0] not in Vo: #if ePrime[0]==u in VoLine: VoLine \cup Vo = V(T)
            MILP += pulp.lpSum(Lambda[gamma[0]] * delta_ac_e_gamma[(ePrime,gamma)] for gamma in gammaInt_ac) == alpha_v[ePrime[0]],f"MILP-73_{ePrime}-1"
        elif ePrime[0] in Vo:
           for mu in range(1,Cmax+1):
               MILP += len(Lambda) * (xNode_u_ePrime_mu[(ePrime[0],ePrime,mu)]-1) <= \
                   pulp.lpSum(Lambda[gamma[0]] * delta_ac_e_gamma[(ePrime,gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(ePrime[0],mu)],f"MILP-73_{ePrime}_{mu}-2"
               MILP += pulp.lpSum(Lambda[gamma[0]] * delta_ac_e_gamma[(ePrime,gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(ePrime[0],mu)] <= \
                       len(Lambda) * (1-xNode_u_ePrime_mu[(ePrime[0], ePrime, mu)]),f"MILP-73_{ePrime}_{mu}-3"
    # -------- Constraint (74) --------
        if ePrime[1] not in Vo:
            MILP += pulp.lpSum(Lambda[gamma[1]]*delta_ac_e_gamma[(ePrime,gamma)] for gamma in gammaInt_ac) == alpha_v[ePrime[1]],f"MILP-74_{ePrime}-1"
        elif ePrime[1] in Vo:
           for mu in range(1,Cmax+1):
               MILP += len(Lambda) * (xNode_u_ePrime_mu[(ePrime[1],ePrime,mu)]-1) <= \
                   pulp.lpSum(Lambda[gamma[1]] * delta_ac_e_gamma[(ePrime,gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(ePrime[1],mu)],f"MILP-74_{ePrime}_{mu}-2"
               MILP += pulp.lpSum(Lambda[gamma[1]] * delta_ac_e_gamma[(ePrime,gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(ePrime[1],mu)] <= \
                       len(Lambda) * (1-xNode_u_ePrime_mu[(ePrime[1], ePrime, mu)]),f"MILP-74_{ePrime}_{mu}-3"
    # -------- Constraint (75) --------
    for e in Eo:#normolly recognize that e=uu' satisfied u < u'
        for v in range(1,2*Cmax-Cmin+1):
            MILP += (len(Lambda) *(xEdge_u_e_v[(e[0],e,v)]-1) <=
                     pulp.lpSum(Lambda[gamma[0]]*delta_ac_e_gamma[(e,gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(e[0],E_Cu[v][0])]),f"MILP-75_{e}_{v}-1"
            MILP += (pulp.lpSum(Lambda[gamma[0]]*delta_ac_e_gamma[(e,gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(e[0],E_Cu[v][0])] <=
                     len(Lambda) *(1 - xEdge_u_e_v[(e[0],e,v)])),f"MILP-75_{e}_{v}-2"
            MILP += (len(Lambda) * (xEdge_u_e_v[(e[0],e,v)] - 1) <=
                     pulp.lpSum(Lambda[gamma[1]] * delta_ac_e_gamma[(e, gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(e[0], E_Cu[v][1])]),f"MILP-75_{e}_{v}-3"
            MILP += (pulp.lpSum(Lambda[gamma[1]] * delta_ac_e_gamma[(e, gamma)] for gamma in gammaInt_ac) - alpha_u_mu[(e[0], E_Cu[v][1])] <=
                     len(Lambda) *(1 - xEdge_u_e_v[(e[0],e,v)])),f"MILP-75_{e}_{v}-4"
    # -------- Constraint (76) --------
    for gamma in gammaInt_ac:
        gamma_bar = (gamma[1],gamma[0],gamma[2])
        if gamma[0] != gamma[1]:
            MILP += (acInt_gamma[gamma] == pulp.lpSum((delta_ac_u_v_gamma[(u,v,gamma)]+delta_ac_u_v_gamma[(u,v,gamma_bar)]) for u in Vo for v in range(1,2*Cmax-Cmin+1)) +
                     pulp.lpSum((delta_ac_e_gamma[(ePrime, gamma)] + delta_ac_e_gamma[(ePrime, gamma_bar)]) for ePrime in Eo_bar) -
                     pulp.lpSum((delta_ac_e_gamma[(e, gamma)] + delta_ac_e_gamma[(e, gamma_bar)]) for e in Eo)),f"MILP-76_{gamma}-1"

        elif gamma[0] == gamma[1]:
            MILP += (acInt_gamma[gamma] ==
                     pulp.lpSum(delta_ac_u_v_gamma[(u,v,gamma)] for u in Vo for v in range(1,2*Cmax-Cmin+1)) +
                     pulp.lpSum(delta_ac_e_gamma[(ePrime,gamma)] for ePrime in Eo_bar) -
                     pulp.lpSum(delta_ac_e_gamma[(e,gamma)] for e in Eo)),f"MILP-76_{gamma}-1"

    return  delta_ac_u_v_gamma, delta_ac_e_gamma , acInt_gamma ,MILP


    # -------- 3.8 Descriptors for the Number of edge-configurations --------

def descriptor_for_the_number_of_ec(
        gammaInt_ec,
        gammaInt_ac,
        ecInt_LB,
        ecInt_UB,
        delta_ac_u_v_gamma,
        delta_ac_e_gamma,
        delta_deg_u_mu_d,
        delta_deg_v_d, #v in VoLine
        xNode_u_e_mu,
        xEdge_u_e_v,
        Vo,
        E_T,
        Eo_bar,
        Eo,
        e_u_i,
        MILP
):
    # Define delta_ec(u,[v];[tau])
    delta_ec_u_v_tau={(u,v,tau):pulp.LpVariable(f"delta_ec({u},{v},{tau})",cat=pulp.LpBinary)
                      for u in Vo for v in range(1,2*Cmax-Cmin+1) for tau in gammaInt_ec}

    # Define delta_ec([e];[tau])
    delta_ec_e_tau = {(e, tau): pulp.LpVariable(f"delta_ec({e},{tau})", cat=pulp.LpBinary)
                        for e in E_T for tau in gammaInt_ec}

    # Define ecInt([tau])
    ecInt_tau = {tau: pulp.LpVariable(f"ecInt({tau})",lowBound=ecInt_LB[tau],upBound=ecInt_UB[tau], cat=pulp.LpInteger)
                   for tau in gammaInt_ec}

    # -----add_constraints-----#

    # -------- Constraint (77) --------
    for u in Vo:
        for v in range(1,2*Cmax-Cmin+1):
            MILP += pulp.lpSum(delta_ec_u_v_tau[(u,v,tau)] for tau in gammaInt_ec) == e_u_i[(u,v)],f"MILP-77_{u}_{v}"
    # -------- Constraint (78) --------
            MILP += (pulp.lpSum(gammaInt_ac[(tau[0][0],tau[1][0],tau[2])]*delta_ec_u_v_tau[(u,v,tau)] for tau in gammaInt_ec) ==
                     pulp.lpSum(gammaInt_ac[gamma]*delta_ac_u_v_gamma[(u,v,gamma)] for gamma in gammaInt_ac)),f"MILP-78_{u}_{v}"
    # -------- Constraint (79) --------
    for u in Vo: #d = tau[0][1]
        for v in range(1, 2 * Cmax - Cmin + 1):
            MILP += pulp.lpSum(tau[0][1] * delta_ec_u_v_tau[(u,v,tau)] for tau in gammaInt_ec) - \
                    pulp.lpSum(d * delta_deg_u_mu_d[(u, E_Cu[v][0], d)] for d in range(1,5)) >= 4 * (e_u_i[(u,v)] - 1),f"MILP-79_{u}_{v}-1"
            MILP += pulp.lpSum(tau[0][1] * delta_ec_u_v_tau[(u, v, tau)] for tau in gammaInt_ec) <= \
                    pulp.lpSum(d * delta_deg_u_mu_d[(u, E_Cu[v][0], d)] for d in range(1,5)), f"MILP-79_{u}_{v}-2"
    # -------- Constraint (80) --------
    for u in Vo:  # d' = tau[1][1]
        for v in range(1, 2 * Cmax - Cmin + 1):
            MILP += pulp.lpSum(tau[1][1] * delta_ec_u_v_tau[(u, v, tau)] for tau in gammaInt_ec) - \
                    pulp.lpSum(d * delta_deg_u_mu_d[(u, E_Cu[v][1], d)] for d in range(1, 5)) >= 4 * (e_u_i[(u, v)] - 1),f"MILP-80_{u}_{v}-1"
            MILP += pulp.lpSum(tau[1][1] * delta_ec_u_v_tau[(u, v, tau)] for tau in gammaInt_ec) <= \
                    pulp.lpSum(d * delta_deg_u_mu_d[(u, E_Cu[v][1], d)] for d in range(1, 5)), f"MILP-80_{u}_{v}-2"
    # -------- Constraint (81) --------
    for e in E_T:
        MILP += pulp.lpSum(delta_ec_e_tau[(e,tau)] for tau in gammaInt_ec) == 1,f"MILP-81_{e}"
    # -------- Constraint (82) --------
        MILP += (pulp.lpSum(gammaInt_ac[(tau[0][0],tau[1][0],tau[2])]*delta_ec_e_tau[(e,tau)] for tau in gammaInt_ec ) ==
                 pulp.lpSum(gammaInt_ac[gamma]*delta_ac_e_gamma[(e,gamma)] for gamma in gammaInt_ac)),f"MILP-82_{e}"
    # -------- Constraint (83) --------
    for ePrime in Eo_bar:
        if ePrime[0] not in Vo:
            MILP += (pulp.lpSum(tau[0][1]*delta_ec_e_tau[(ePrime,tau)] for tau in gammaInt_ec) ==
                     pulp.lpSum(d*delta_deg_v_d[(ePrime[0],d)] for d in range(1,max_degree+1))),f"MILP-83_{ePrime}-1"
        if ePrime[0] in Vo:
            for mu in range(1,Cmax+1):
                MILP += (pulp.lpSum(tau[0][1] * delta_ec_e_tau[(ePrime,tau)] for tau in gammaInt_ec) - pulp.lpSum(d * delta_deg_u_mu_d[(ePrime[0],mu,d)] for d in range(1,max_degree+1)) >=
                         4 * (xNode_u_e_mu[(ePrime[0],ePrime,mu)]-1)),f"MILP-83_{ePrime}_{mu}-2"

                MILP += (pulp.lpSum(tau[0][1] * delta_ec_e_tau[(ePrime, tau)] for tau in gammaInt_ec) - pulp.lpSum(d * delta_deg_u_mu_d[(ePrime[0], mu, d)] for d in range(1, max_degree + 1)) <=
                         4 * (1 - xNode_u_e_mu[(ePrime[0], ePrime, mu)])),f"MILP-83_{ePrime}_{mu}-3"
    # -------- Constraint (84) --------
        if ePrime[1] not in Vo:
            MILP += (pulp.lpSum(tau[1][1] * delta_ec_e_tau[(ePrime,tau)] for tau in gammaInt_ec) ==
                     pulp.lpSum(dP*delta_deg_v_d[(ePrime[1],dP)] for dP in range(1,max_degree+1))),f"MILP-84_{ePrime}-1"
        if ePrime[1] in Vo:
            for mu in range(1,Cmax+1):
                MILP += (pulp.lpSum(tau[1][1] * delta_ec_e_tau[(ePrime,tau)] for tau in gammaInt_ec) - pulp.lpSum(dP * delta_deg_u_mu_d[(ePrime[1],mu,dP)] for dP in range(1,max_degree + 1)) >=
                         4 * (xNode_u_e_mu[(ePrime[1],ePrime,mu)]-1)),f"MILP-84_{ePrime}_{mu}-2"

                MILP += (pulp.lpSum(tau[1][1] * delta_ec_e_tau[(ePrime, tau)] for tau in gammaInt_ec) - pulp.lpSum(dP * delta_deg_u_mu_d[(ePrime[1], mu, dP)] for dP in range(1, max_degree + 1)) <=
                         4 * (1 - xNode_u_e_mu[(ePrime[1], ePrime, mu)])),f"MILP-84_{ePrime}_{mu}-3"
    # -------- Constraint (85) --------
    for e in Eo:
        for v in range(1,1+2*Cmax-Cmin):
            MILP += (4 * (xEdge_u_e_v[(e[0],e,v)]-1) <=
                     pulp.lpSum(tau[0][1] * delta_ec_e_tau[(e,tau)] for tau in gammaInt_ec) - pulp.lpSum(d*delta_deg_u_mu_d[(e[0],E_Cu[v][0],d)] for d in range(1,1+max_degree))),f"MILP-83_{e}_{v}-1"
            MILP += (pulp.lpSum(tau[0][1] * delta_ec_e_tau[(e,tau)] for tau in gammaInt_ec) - pulp.lpSum(d*delta_deg_u_mu_d[(e[0],E_Cu[v][0],d)] for d in range(1,1+max_degree)) <=
                     4 * (1-xEdge_u_e_v[(e[0],e,v)])),f"MILP-85_{e}_{v}-2"

            MILP += (4 * (xEdge_u_e_v[(e[0],e,v)]-1) <=
                     pulp.lpSum(tau[1][1] * delta_ec_e_tau[(e,tau)] for tau in gammaInt_ec) - pulp.lpSum(dP * delta_deg_u_mu_d[(e[0],E_Cu[v][1],dP)] for dP in range(1,1+max_degree))),f"MILP-83_{e}_{v}-3"
            MILP += (pulp.lpSum(tau[1][1] * delta_ec_e_tau[(e,tau)] for tau in gammaInt_ec) - pulp.lpSum(dP * delta_deg_u_mu_d[(e[0],E_Cu[v][1],dP)] for dP in range(1,1+max_degree)) <=
                     4 * (1-xEdge_u_e_v[(e[0],e,v)])),f"MILP-85_{e}_{v}-4"
    # -------- Constraint (86) --------
    for tau in gammaInt_ec:
        if tau[0] != tau[1]:
            tau_bar=(tau[1],tau[0],tau[2])
            MILP += (ecInt_tau[tau] ==
                     pulp.lpSum( (delta_ec_u_v_tau[(u,v,tau)]+delta_ec_u_v_tau[(u,v,tau_bar)]) for u in Vo for v in range(1,2*Cmax-Cmin+1)) +
                     pulp.lpSum((delta_ec_e_tau[(ePrime,tau)]+ delta_ec_e_tau[(ePrime,tau_bar)]) for ePrime in Eo_bar) -
                     pulp.lpSum((delta_ec_e_tau[(e,tau)]+delta_ec_e_tau[(e,tau_bar)]) for e in Eo)),f"MILP-86_{tau}-1"
        elif tau[0] == tau[1]:
            MILP += (ecInt_tau[tau] ==
                     pulp.lpSum(delta_ec_u_v_tau[(u,v,tau)] for u in Vo for v in range(1,2*Cmax-Cmin+1)) +
                     pulp.lpSum(delta_ec_e_tau[(ePrime,tau)] for ePrime in Eo_bar) -
                     pulp.lpSum(delta_ec_e_tau[(e,tau)] for e in Eo)),f"MILP-86_{tau}-2"
    return delta_ec_u_v_tau, delta_ec_e_tau, ecInt_tau, MILP

def prepare_fv(
    training_data_filename,
    Xi_T,
    cc_ksi,
    Lambda_int,
    Lambda_ex,
    gammaInt_ec,
    gamma_lf_ac,
    index_set,
    rank,
    n_G,
    n_int,
    msLine,
    deg_d,degInt_d,
    bd_m,
    naInt_a, naEx_a,
    ecInt_tau,
    fc_phi,
    ac_lf_gamma
):
    ''' A function to prepare some variables about the descriptors
        from the feature vector '''

    descriptors = dict()      # get variable from fv file
    fv = pd.read_csv(training_data_filename, sep=",")
    num_fv = len(list(fv.columns))

    # prepare for normalization and standardization
    max_dcp = dict()
    min_dcp = dict()
    for i in range(1, num_fv):
        max_dcp[i] = fv.iloc[:, i].max()
        min_dcp[i] = fv.iloc[:, i].min()

    for i, fv_name in enumerate(list(fv.columns)):
        if fv_name == "CID":
            pass
        elif fv_name == "rank":
            descriptors[i] = rank

        elif fv_name == "n":
            descriptors[i] = n_G

        elif fv_name == "n_in":
            descriptors[i] = n_int

        elif fv_name == "ms":
            descriptors[i] = msLine

        elif fv_name[0] == "d":
            if fv_name[3] == "i":
                j = int(fv_name[6])
                descriptors[i] = degInt_d[j]

            else:
                j = int(fv_name[3])
                descriptors[i] = deg_d[j]

        elif fv_name == "bd_in_2":
            descriptors[i] = bd_m[2]

        elif fv_name == "bd_in_3":
            descriptors[i] = bd_m[3]

        elif fv_name[0] == "n" and fv_name[3] == "i":
            ele = fv_name[6:] # may need change, Zhu, 0109
            if ele in Lambda_int:
                descriptors[i] = naInt_a[ele]

            else:
                descriptors[i] = 0

        elif fv_name[0] == "n" and fv_name[3] == "e":
            ele = fv_name[6:]
            if ele in Lambda_ex:
                descriptors[i] = naEx_a[ele]
            else:
                descriptors[i] = 0

        elif fv_name[0] == "F" and fv_name[1] == "C":
            index = int(fv_name.split("_")[1])
            if index in index_set.keys():
                descriptors[i] = fc_phi[index_set[index]]

            else:
                descriptors[i] = 0

        elif fv_name[0] == "L":
            lines = fv_name.split("_")
            nu = (lines[1], lines[2], int(lines[3]))
            if nu in gamma_lf_ac:
                descriptors[i] = ac_lf_gamma[nu]
            else:
                descriptors[i] = 0

        elif fv_name[0] == "C" and fv_name[1] == "C":
            lines = fv_name.split("_")
            ksi =tuple(int(ele) for ele in lines[2:])
            if ksi in Xi_T:
                descriptors[i] = cc_ksi[ksi]
            else:
                descriptors[i] = 0

        else:
            str_tmp = fv_name.split("_")
            ec_tmp1 = ((str_tmp[0], int(str_tmp[1])), (str_tmp[2], int(str_tmp[3])), int(str_tmp[4]))
            ec_tmp2 = ((str_tmp[2], int(str_tmp[3])), (str_tmp[0], int(str_tmp[1])), int(str_tmp[4]))
            if (ec_tmp1 in gammaInt_ec) or (ec_tmp2 in gammaInt_ec):
                descriptors[i] = ecInt_tau[ec_tmp1]
            else:
                descriptors[i] = 0

    return descriptors, num_fv, max_dcp, min_dcp

def add_constraints_nor_std_fv(
    MILP,
    num_fv,
    descriptors,
    max_dcp,
    min_dcp,
    eps
):
    x_hat = {i: pulp.LpVariable(f"x_hat({i})") for i in range(1, num_fv)} #标准化以后的相对应的值

    for i in range(1, num_fv):
        if max_dcp[i] == min_dcp[i]:
            MILP += x_hat[i] == 0, f"milp-2LMH-std-{i}-1(nor)"
        else:
            MILP += x_hat[i] >= (descriptors[i] - eps - min_dcp[i]) / (max_dcp[i] - min_dcp[i]), f"milp-2LMH-std-{i}-1(nor)"
            MILP += x_hat[i] <= (descriptors[i] + eps - min_dcp[i]) / (max_dcp[i] - min_dcp[i]), f"milp-2LMH-std-{i}-2(nor)"

    return MILP, x_hat


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

def add_constraints__LR( MILP,
                        descriptors,
                        num_fv,
                        LR,
                        prop = "def"):

    for i in range(1, num_fv):
        MILP += LR.weight_var[i] == descriptors[i], \
            "LR_connect_x_{}_{}".format(i, prop)

    return MILP

def output_sdf_file(
        Lambda,
        F_t,
        F_u,
        Vo,
        Vo_bar,
        Eo_u,
        Eo_bar_u,
        Eo,
        Eo_bar,

        #variable
        xEdge_u_e_v, #ring edge 给了哪个边
        xNode_u_ep_mu,#non ring node 给了哪个点
        e_u_i,#indicating whether the edge ei in Cu is used;
        beta_e,#that stores the bond-multiplicty of the edge e in E_T
        beta_u_i,#that stores the bond-multiplicty of the edge ei in Cu;
        delta_F_u_mu_phi,# indicating whether the fringe- tree phi is attached to vertex μ ∈ V (Cu);
        delta_F_v_phi, #indicating whether the fringe-tree phi is attached to node v ∈ V (T) \ V ◦;
        alpha_u_mu, #that represents the chemical element assigned to the vertex μ in Cu;
        alpha_v, #that represents the chemical element assigned to the vertex v ∈ V (T) \ V ◦;

        output_path,
        instance_name
):
    # n_tmp = 1000
    Lambda_bar = {v: k for k, v in Lambda.items()}  #exchange value and key

    node_tmp = {}#store the node in CG, use the name :node in Cu -> (u,mu), non ring node -> v ,node in fringe on Cu -> (u,mu,i), node in fringe on non ring node v-> ("F",v,i)
    edge_tmp = {}#store the edge in CG.
    graph_node = {}
    graph_adj = {}
    node_index= {}

    chg_index = list()
    chg = list()
    chg_ele_index = dict()
    chg_ele = dict()

    two_in_one = dict()

    for e in Eo:
        l_small, l_big, r_small, r_big = 0, 0, 0, 0
        for v in range(1, 1 + 2 * Cmax - Cmin):
            if round(xEdge_u_e_v[(e[0], e, v)].value()) != 0:
                l_small = E_Cu[v][0]
                l_big = E_Cu[v][1]
            if round(xEdge_u_e_v[(e[1], e, v)].value()) != 0:
                r_small = E_Cu[v][0]
                r_big = E_Cu[v][1]
        two_in_one[(e[0], l_small)] = (e[1], r_big)
        two_in_one[(e[0], l_big)] = (e[1], r_small)

    #以ring node展开计算
    for u in Vo:
        for i in range(1, 1 + 2 * Cmax - Cmin):
            if round(e_u_i[(u, i)].value()) != 0:
                node_tmp[(u, E_Cu[i][0])] = Lambda_bar[round(alpha_u_mu[(u, E_Cu[i][0])].value())]
                node_tmp[(u, E_Cu[i][1])] = Lambda_bar[round(alpha_u_mu[(u, E_Cu[i][1])].value())]
                edge_tmp[( (u, E_Cu[i][0]), (u, E_Cu[i][1]) )] = round(beta_u_i[(u,i)].value())

        #attach fringe tree on a circle
        for mu in range(1,1+Cmax):
            if (u,mu) not in two_in_one:
                for phi_tmp in F_u[u]:
                    if round(delta_F_u_mu_phi[(u,mu,phi_tmp)].value()) != 0:
                        phi = F_t[phi_tmp]
                        for j in range(1, len(phi.vertex)):
                            a_tmp, h_tmp = phi.vertex[j]
                            node_tmp[(u,mu,j)] = a_tmp #name is (u,mu,j)

                        for j1 in range(len(phi.vertex)):
                            for j2 in range(j1 + 1, len(phi.vertex)):
                                if j1 == 0 and phi.beta[j1][j2] != 0:
                                    edge_tmp[((u,mu),(u,mu,j2))] = phi.beta[0][j2]
                                elif j1 != 0 and phi.beta[j1][j2] != 0:
                                    edge_tmp[((u,mu,j1),(u,mu,j2))] = phi.beta[j1][j2]

                        if phi.chg[0] != 0:
                            chg.append((u,mu))
                            chg.append(phi.chg[0])
                            chg_ele[(u,mu)] = 4-phi.chg[0]
                        else:
                            chg_ele[(u,mu)] = 0

                        for j in range(1, len(phi.vertex)):
                            if phi.chg[j] != 0:
                                chg.append((u, mu, j))
                                chg.append(phi.chg[j])
                                chg_ele[(u, mu, j)] = 4 - phi.chg[j]
                            else:
                                chg_ele[(u, mu, j)] = 0

    #添加non-ring node
    for v in Vo_bar:
        node_tmp[v] = Lambda_bar[round(alpha_v[v].value())]

    # attach fringe tree on non ring node
    for v in Vo_bar:
        for phi_tmp in F_u[v]:
            if round(delta_F_v_phi[(v, phi_tmp)].value()) != 0:
                phi = F_t[phi_tmp]
                for j in range(1, len(phi.vertex)):
                    a_tmp, h_tmp = phi.vertex[j]
                    node_tmp[("F", v, j)] = a_tmp

                for j1 in range(len(phi.vertex)):
                    for j2 in range(j1 + 1, len(phi.vertex)):
                        if j1 == 0 and phi.beta[j1][j2] != 0:
                            edge_tmp[((v), ("F", v, j2))] = phi.beta[j1][j2]
                        elif j1 != 0 and phi.beta[j1][j2] != 0:
                            edge_tmp[(("F", v, j1), ("F", v, j2))] = phi.beta[j1][j2]

                if phi.chg[0] != 0:
                    chg.append(v)
                    chg.append(phi.chg[0])
                    chg_ele[v] = 4 - phi.chg[0]
                else:
                    chg_ele[v] = 0
                for j in range(1, len(phi.vertex)):
                    if phi.chg[j] != 0:
                        chg.append(("F", v, j))
                        chg.append(phi.chg[j])
                        chg_ele[("F", v, j)] = 4 - phi.chg[j]
                    else:
                        chg_ele[("F", v, j)] = 0

    # two_in_one = dict()
    # for e in Eo:
    #     l_small,l_big,r_small,r_big = 0,0,0,0
    #     for v in range(1, 1 + 2 * Cmax - Cmin):
    #         if round(xEdge_u_e_v[(e[0],e,v)].value()) != 0:
    #             l_small = E_Cu[v][0]
    #             l_big = E_Cu[v][1]
    #         if round(xEdge_u_e_v[(e[1],e,v)].value()) != 0:
    #             r_small = E_Cu[v][0]
    #             r_big = E_Cu[v][1]
    #     two_in_one[(e[0], l_small)] = (e[1], r_big)
    #     two_in_one[(e[0],l_big)] = (e[1],r_small)

    #计算non ring edge，1，ring to ring 2.non to ring 3. non to non
    for ep in Eo_bar:
        if ep[0] in Vo and ep[1] in Vo:
            mu1 = 0
            mu2 = 0
            for mu in range (1, 1 + Cmax):
                if(round(xNode_u_ep_mu[(ep[0],ep,mu)].value() ) ) != 0:
                    mu1 = mu
                if(round(xNode_u_ep_mu[ep[1],ep,mu].value() ) ) != 0:
                    mu2 = mu
            edge_tmp[( (ep[0],mu1), (ep[1],mu2) )] = round(beta_e[ep].value())
        elif ep[0] in Vo and ep[1] not in Vo:
            for mu in range(1, 1 + Cmax):
                if(round(xNode_u_ep_mu[(ep[0],ep,mu)].value() )) != 0:
                    edge_tmp[( (ep[0],mu), ep[1] )] = round(beta_e[ep].value())
        elif ep[0] not in Vo and ep[1] in Vo:
            for mu in range(1, 1 + Cmax):
                if(round(xNode_u_ep_mu[(ep[1],ep,mu)].value() ) ) != 0:
                    edge_tmp[( (ep[1],mu), ep[0] )] = round(beta_e[ep].value())
        elif ep[0] not in Vo and ep[1] not in Vo:
            edge_tmp[( ep[0], ep[1] )] = round(beta_e[ep].value())

    #连接两个圆环
    for del_node in two_in_one:
        new_node_tmp = dict()
        new_edge_tmp = dict()
        for node in node_tmp:
            if del_node != node:
                new_node_tmp[node] = node_tmp[node]
        node_tmp = new_node_tmp

        for edge in edge_tmp:
            if del_node not in edge:
                new_edge_tmp[edge] = edge_tmp[edge]
            else:
                if edge[0] == del_node:
                    new_edge_tmp[ (two_in_one[del_node], edge[1]) ] = edge_tmp[edge]
                elif edge[1] == del_node:
                    new_edge_tmp[ (two_in_one[del_node], edge[0]) ] = edge_tmp[edge]
        edge_tmp = new_edge_tmp

        # chg_ele.pop(del_node)

    x = 0
    for node in node_tmp:
        x += 1
        node_index[node] = x #index of node is x
        graph_node[x] = node_tmp[node]

    for ele in chg:
        if chg.index(ele) % 2 == 0:
            chg_index.append(node_index[ele])
        else:
            chg_index.append(ele)
    for ele in chg_ele:
        chg_ele_index[node_index[ele]] = chg_ele[ele]
    #index all edge
    for edge in edge_tmp:
        graph_adj[(node_index[edge[0]],node_index[edge[1]])] = edge_tmp[edge]

    new_graph_adj = dict()
    for edge in graph_adj:
        new_graph_adj[tuple(sorted(edge))] = graph_adj[edge]
    graph_adj = new_graph_adj

    n = len(graph_node)
    m = len(graph_adj)

    with open(output_path, "w") as f:
        f.write(f"{instance_name}\n")
        f.write("MILP_2LMMCC\n")
        f.write("MILP_2LMM_INSTANCE_1\n")
        f.write("{:3}{:3}  0  0  0  0  0  0  0  0999 V2000 \n".format(n, m))

        for node in graph_node:
            atom_tmp = re.sub(r'[0-9]', '', graph_node[node])
            f.write("    0.0000    0.0000    0.0000 {:2}  0  {:1}  0  0  0  0  0  0  0  0  0  0\n".format(atom_tmp,chg_ele_index[node]))
        for edge in graph_adj:
            f.write("{:3}{:3}{:3}  0  0  0  0\n".format(edge[0], edge[1],graph_adj[edge]))

        if len(chg_index) != 0:
            f.write("M  CHG{:3}".format(int(len(chg_index) / 2)))
            for tmp in chg_index:
                f.write("{:4}".format(tmp))
            f.write("\n")

        f.write("M  END\n")
        f.write("$$$$\n")
        f.close()


