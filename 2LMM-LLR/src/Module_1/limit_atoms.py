#!/Library/Frameworks/Python.framework/Versions/3.7/bin/python3

import sys, os

# Fix Python 2.x.
try: input = raw_input
except NameError: pass


def getSdfString(inFile):
    ans = str()
    ll = "somestring"
    while ll != "" and ll.find("$$$$") == -1:
        ll = inFile.readline()
        ans += ll
    return ans

def parseSdfString(sdfString):
    lines = sdfString.splitlines()
    #print(lines)
    graphName = lines[0].strip()
    try:
        #print(lines[3][:3])
        numVerts = int(lines[3][:3])
        #print(numVerts)
        #print(lines[3][3:6])
        numEdges = int(lines[3][3:6])
        #print(numEdges)
        
    except:
        print(
            '''
            There is something strange with the sdf format.
            Could not read number of vertices and edges
            ''')
        quit()
    vrtxLineBegin = 4
    vertexColors = list()
    for line in lines[vrtxLineBegin : vrtxLineBegin + numVerts]:
        try:
            info = line.split()
            vLabel = info[3] # the atom type in the SDF format
            vertexColors.append(vLabel)
        except:
            print(
                '''
                There was an error in processing
                the connection table
                {}
                '''.format(line))
            quit()
    if numVerts != len(vertexColors):
        print(
            '''
            Error, the number of extracted vertices does not match at {}
            '''.format(graphName)
        )
        quit()

    edgeLineBegin = numVerts+4
    #print(edgeLineBegin)
    edgeMult = dict()
    #print(lines[edgeLineBegin : edgeLineBegin + numEdges])
    for line in lines[edgeLineBegin : edgeLineBegin + numEdges]:
        try:
            v1 = int(line[:3]) - 1
            v2 = int(line[3:6]) - 1
            mult = int(line[6:9])
        except:
            print(
                '''
                There was an error in processing
                the connection table
                {}
                '''.format(line))
            quit()
        edgeMult[(v1, v2)] = mult
    if numEdges != len(edgeMult):
        print(
            '''
            Warning: the number of extracted edges does not match at {} ({}!={})
                     this guy has {}
            '''.format(graphName,numEdges,len(edgeMult),set(vertexColors))
        )
        #quit()
    return set(vertexColors)

# small test
def main(argv):
    if len(argv) < 3: 
        print("Please supply an input file name and set of symbols")
        sys.exit()

    ### set is changed to list on 20201201 ###
    #input_set = set(argv[2:])
    input_set = list(argv[2:])
    ##########################################
    inf = argv[1]
    outf = inf[:-4]
    outf = "_".join([outf, *input_set])
    counter = 0
    print(outf)
    with open(outf, "w", newline='\n') as off:
        with open(inf) as iff:
            sdfString = getSdfString(iff)
            while (sdfString):
                cur_set = parseSdfString(sdfString)
                if not cur_set.difference(input_set):
                    off.write(sdfString)
                    counter += 1
                sdfString = getSdfString(iff)
    #new_outf = outf + f"_{counter}.sdf"
    new_outf = outf + f".sdf"
    os.rename(outf, new_outf)
    # comment out on 20201201
    #with open(new_outf+"memo", "w", newline='\n') as ff:
    #    pass
    #    # just keep a memo of the output result...
   
                
main(sys.argv)

#def test():
#    main(["1", "../sdf/C9_N1_O3.sdf", "O", "N", "C"])
                
# test()

