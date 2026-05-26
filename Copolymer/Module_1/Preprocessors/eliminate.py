#!/usr/bin/python

# references:
#   - http://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx

import os,sys

print("=== {} ===".format(sys.argv[0]))

#
# This script ADMITS charged/multivalence molecules.
#
# v062: you can specify the output filename.
#       the number of atoms is not considered for eliminated molecules
#       this is essential difference from v061


try:
    import networkx as nx
except:
    sys.stderr.write('error: you need networkx.\n')
    sys.stderr.write('please install networkx as follows.\n')
    sys.stderr.write('pip3 install networkx.\n')
    exit(1)

ChgExpr = {7:-3, 6:-2, 5:-1, 0:0, 3:1, 2:2, 1:3, 4:"radical"}

############################################
def read_mol(fp):
    SDF = []
    ### initialize flags ###
    is_charged = False
    is_radical = False
    is_aromar = False
    
    ### CID ###
    s = fp.readline()
    SDF.append(s)

    if s == "":
        return (None, None, None, None, None, None, None, None, None, None, None)

    cid = s.replace('\n','')
    tmp = fp.readline()
    SDF.append(tmp)
    tmp = fp.readline()
    SDF.append(tmp)

    ### n & m ###
    s = fp.readline()
    SDF.append(s)
    n = int(s[:3])
    m = int(s[3:6])

    ### nodes ###
    Alpha = {}
    Chg = {}
    for i in range(1,n+1):
        s = fp.readline()
        SDF.append(s)
        elem = s[31:34].replace(' ','')
        Alpha[i] = elem
        monoisotope = int(s[34:36])  ### currently not used
        Chg[i] = ChgExpr[int(s[36:39])]
        if Chg[i] == "radical":
            is_radical = True
        elif Chg[i] != 0:
            is_charged = True

    ### edges ###
    Edge = {}
    Deg = {}
    for i in range(1,n+1):
        Deg[i] = 0
    for _ in range(1,m+1):
        s = fp.readline()
        SDF.append(s)
        i = int(s[:3])
        j = int(s[3:6])
        bond = int(s[6:9])
        Edge[(i,j)] = bond
        if bond == 4:
            is_aromar = True
        if Alpha[j] != 'H':
            Deg[i] += 1
        if Alpha[i] != 'H':
            Deg[j] += 1
    
    ### detect "M CHG" or "M RAD" in the remaining part ###
    while 1:
        s = fp.readline()
        SDF.append(s)
        if s.find("$$$$") != -1:
            break
        if s.find("M  CHG") != -1:
            num = int(s[6:9])
            idx = 9
            for _ in range(num):
                i = int(s[idx:idx+4])
                chg = int(s[idx+4:idx+8])
                Chg[i] = chg
                idx += 8
            is_charged = True
        if s.find("M  RAD") != -1:
            is_radical = True

    return (cid, n, m, Alpha, Chg, Edge, Deg, SDF, is_charged, is_radical, is_aromar)

############################################
def decide_connectivity(n, m, E):
    if n<=1:
        return False
    import networkx as nx
    G = nx.Graph()
    for i in range(1,n+1):
        G.add_node(i)
    for (i,j) in E:
        G.add_edge(i,j)
    return nx.is_connected(G)

############################################
def write_mol(fp, data):
    for line in data:
        fp.write(line)    

############################################
def isEliminated(cid, n, m, A, C, E, D, SDF, is_charged, is_radical, is_aromar):
    ##### test whether the current molecule should be excluded #####
    # decide the graph has radical
    # is_radical
        
    # decide the connectivity
    is_connected = decide_connectivity(n, m, E)
            
    # count the carbon num
    carbon_num = len([i for i in A if A[i]=='C'])
        
    # decide if the graph has node with degree>=5
    is_large_deg = False
    if len(D)>0 and max(list(D.values())) >= 5:
        is_large_deg = True
    
    if is_radical == True or is_connected == False or carbon_num < 4 or is_large_deg == True:
        sys.stdout.write("warning: cid {} is eliminated:".format(cid))
        if is_radical == True:
            sys.stdout.write(' radical')
        if is_connected == False:
            sys.stdout.write(' disconnected')
        if carbon_num < 4:
            sys.stdout.write(' |C|={}<=3'.format(carbon_num))
        if is_large_deg == True:
            sys.stdout.write(' degree={}'.format(max(list(D.values()))))
        sys.stdout.write("\n")
        return True
    else:
        return False

############################################
def main(argv):
    if len(argv) < 3: 
        print("usage: {} (INPUT.sdf) (OUTPUT.sdf)".format(argv[0]))
        sys.exit()

    SDFfile = argv[1]
    SDFoutputfile = argv[2]
        
    inputFile = open(SDFfile, 'r')
    outputFile = open(SDFoutputfile, 'w')

    lambdaDict = dict()

    Num = 0
    Num_radical = 0
    Num_disconnected = 0
    Num_few_carbons = 0
    Num_feasible = 0
    Num_charged = 0
    Num_aromar = 0
    Num_large_deg = 0

    while 1:
        ##### read SDF file #####
        cid, n, m, A, C, E, D, SDF, is_charged, is_radical, is_aromar = read_mol(inputFile)
        if cid == None:
            break

        Num += 1

        ##### count the number of compounds have certain atom #####
        '''
        atomSet = set()
        for a in A.values():
            atomSet.add(a)
        for a in atomSet:
            if a in lambdaDict:
                lambdaDict[a] += 1
            else:
                lambdaDict[a] = 1
        '''
        
        ##### test whether the current molecule should be excluded #####
        # the molecule should be eliminated if:
        # - the number of carbon is less than 4,
        # - radical is contained,
        # - the graph is disconnected, and
        # - there is an atom whose degree in H-suppressed model is >= 5.

        warning_message = ''

        # decide the graph contains a radical
        if is_radical == True:
            Num_radical += 1
            warning_message += ' radical'
            
        # decide the connectivity
        is_connected = decide_connectivity(n, m, E)
        if is_connected == False:
            Num_disconnected += 1
            warning_message += ' disconnected'
                
        # count the carbon num
        carbon_num = len([i for i in A if A[i]=='C'])
        if carbon_num < 4:
            Num_few_carbons += 1
            warning_message += ' |C|={}<=3'.format(carbon_num)
            
        # decide if the graph has node with degree>=5
        is_large_deg = False
        if len(D)>0 and max(list(D.values())) >= 5:
            is_large_deg = True
            Num_large_deg += 1
            warning_message += ' degree={}'.format(max(list(D.values())))
        
        if is_radical == True or is_connected == False or carbon_num < 4 or is_large_deg == True:
            sys.stdout.write('warning: cid ' + str(cid) + ' is eliminated:'  +warning_message + '\n')
            continue

        ##### write feasible compounds #####
        write_mol(outputFile, SDF)
        Num_feasible += 1
        if is_charged == True:
            Num_charged += 1
        if is_aromar == True:
            Num_aromar += 1

        ##### count the number of compounds have certain atom #####
        atomSet = set()
        for a in A.values():
            atomSet.add(a)
        for a in atomSet:
            if a in lambdaDict:
                lambdaDict[a] += 1
            else:
                lambdaDict[a] = 1

    inputFile.close()
    outputFile.close()

    lambdaSorted = sorted(lambdaDict.items(), key=lambda l:l[1], reverse=True)

    print("\n---------- reading {} is finished ----------".format((SDFfile.split('/')[-1])))
    print("NUM_ALL_GRAPHS:", Num)
    print("NUM_RADICAL:", Num_radical)
    print("NUM_DISCONNECTED:", Num_disconnected)
    print("NUM_FEW_CARBONS:", Num_few_carbons)
    print("NUM_LARGE_DEG:", Num_large_deg)
    print("NUM_FEASIBLE:", Num_feasible)
    print("NUM_CHARGED:", Num_charged)
    print("NUM_AROMAR:", Num_aromar)

    print("Finished. " + str(Num_feasible) + " compounds are extracted.")
    for (elem, value) in lambdaSorted:
        print(elem + ":" + str(value))

main(sys.argv)
