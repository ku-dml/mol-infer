# -*- coding: utf-8 -*-
import sys
import os
import re

class toEliminate():

    def valDict(self):
        """
        Give each element an assumpted valence.
        Return a dictionary.
        """
        val = dict()
        val["C"] = 4
        val["S"] = 2
        val["O"] = 2
        val["N"] = 3
        val["Cl"] = 1
        val["Br"] = 1
        val["I"] = 1
        val["F"] = 1
        val["Si"] = 4
        val["P"] = 3
        val["Pb"] = 2
        val["Al"] = 3

        return val    
    
    def toEli(self, inFile, valDict):
        """
        Given an SDfile open for reading.
        read the information of molecules in sdf format.
        the molecule should be eliminated if:
        - the number of carbon is less than 4,
        - has atoms with charge,
        - valence of some elements is not the assumpted one.
        Return a list of CIDs of eliminated molecules,
               the string describing the molecules.
        """
        sdfString = str()
        cidSet = set()
        ll = inFile.readline()
        lambdaDict = dict()
        while ll != "":
            ans = str()
            charged = 0
            while ll.find("$$$$") == -1:
                ans += ll
                # check charge
                if re.match("M(.*)CHG", ll):
                    charged = 1
                ll = inFile.readline()
            ans += ll
            lines = ans.splitlines()
            #cid = int(re.findall(r"\d+",lines[0].strip())[0])
            cid = lines[0].replace("\n","")

            atom_num = int(lines[3][0:3])
            bond_num = int(lines[3][3:6])
            # print(atom_num)
            # print(bond_num)
            carbon_num = 0
            ver_dict = dict()
            deg_dict = dict()
            # get number of carbon
            # atomset
            atomSet = set()
            for i in range(1, atom_num + 1):
                elem = lines[i + 3].split()[3]
                if elem == 'C':
                    carbon_num += 1
                ver_dict[i] = elem
                deg_dict[i] = 0
                if elem not in atomSet:
                    atomSet.add(elem)
            # print(carbon_num)
            # print(ver_dict)
            for a in atomSet:
                if a in lambdaDict.keys():
                    lambdaDict[a] += 1
                else:
                    lambdaDict[a] = 1

            for j in range(1, bond_num + 1):
                # print(len(lines[j + atom_num + 3].split()))
                v1 = int(lines[j + atom_num + 3][0:3])
                v2 = int(lines[j + atom_num + 3][3:6])
                mul = int(lines[j + atom_num + 3][6:9])
                # print("{}_{}_{}".format(v1,v2,mul))
                deg_dict[v1] += mul
                deg_dict[v2] += mul
            #print(deg_dict)
            # check degree of each vertex 
            un_val = 0
            for k in range(1, atom_num + 1):
                if ver_dict[k] in valDict.keys():
                    if deg_dict[k] > valDict[ver_dict[k]]:
                        un_val += 1

            if charged == 0 and carbon_num > 3 and un_val == 0:
                sdfString += ans
                if cid in cidSet:
                    sys.stderr.write("warning: duplication of cid {} is found.\n".format(cid))        
                cidSet.add(cid)
            else:
                sys.stderr.write("warning: cid {} is eliminated.\n".format(cid))        

            ll = inFile.readline()


        return cidSet, sdfString, lambdaDict

    def writeSdf(self, sdfString, cidSet, inFile, filename=None):
        #filename = inFile.split('.')[0] + '_eli_' + str(len(cidSet)) + '.sdf' 
        if not filename:
            filename = inFile.split('.')[0] + '_eli.sdf' 
        with open(filename, 'w') as f:
            f.write(sdfString)

def main(argv):
    if len(argv) < 1: 
        print("Please input the SDFfile")
        sys.exit()
    SDFfile = argv[1]
    c = toEliminate()
    valDict = c.valDict()
    with open(SDFfile,'r') as inFile:
        cidSet, sdfString, lambdaDict = c.toEli(inFile, valDict)
    c.writeSdf(sdfString, cidSet, SDFfile, argv[2] if 2 < len(argv) else None)
    print("Done.")
                
main(sys.argv)
