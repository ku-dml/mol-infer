#!/usr/bin/env python3

############################################################
NOTE ='''
usage: norm_csv2_csv.py (TrRawCSV)(TrNormCSV)[(TsRawCSV)(TsNormCSV)]

The script normalizes feature vectors in TrRawCSV (in our FV format)
into [0,1]-vectors by linear transformation (i.e., (x-x_{min})/(x_{max}-x_{min})
and then writes them into TrNormCSV.

If TsRawCSV (optional) is also supplied, then the feature vectors
are also transformed using min/max values in TrRawCSV
and then are written in TsNormCSV.
The values in TsNormCSV are not necessarily in [0,1].

'''
############################################################

import numpy as np
import pandas as pd
import sys

try:
    TrRawCSV, TrNormCSV = sys.argv[1], sys.argv[2]
    if len(sys.argv)>3:
        TsRawCSV, TsNormCSV = sys.argv[3], sys.argv[4]
    else:
        TsRawCSV, TsNormCSV = "", ""
except:
    sys.stderr.write(NOTE)
    exit(1)


def outputNormalizedCSV(Data, Min, Max, output_csv):
    # open output_csv
    f = open(output_csv, "w")

    # write header to output_csv
    f.write('CID')
    for desc in Data.columns:
        f.write(',{}'.format(desc))
    f.write('\n')

    # write normalized vectors to output_csv
    for i,cid in enumerate(Data.index):
        fv = Data.values[i]
        s = str(cid)
        for j,desc in enumerate(Data.columns):
            if Min[desc]==Max[desc]:
                s += ",0"
            else:
                s += ",{:.6f}".format((fv[j]-Min[desc])/(Max[desc]-Min[desc]))
        f.write(s+"\n")
    f.close()
    print('# {} vectors are normalized and written to {}'.format(i+1,output_csv))

# read TrRawCSV    
Tr = pd.read_csv(TrRawCSV, header=0, index_col=0)

# compute min/max values for each descriptor
Min = {}
Max = {}
print("#name\tmin\tmax")
for desc in Tr.columns:
    Min[desc] = Tr[desc].min()
    Max[desc] = Tr[desc].max()
    print('{}\t{}\t{}'.format(desc,Min[desc],Max[desc]))

# output normalized  feature vectors to TrNormCSV
outputNormalizedCSV(Tr, Min, Max, TrNormCSV)

# read TsRawCSV if specified
if TsRawCSV == '':
    sys.exit(0)
Ts = pd.read_csv(TsRawCSV, header=0, index_col=0)

# check compatibility
if len(Tr.columns) != len(Ts.columns):
    sys.stderr.write('error: {} and {} are not compatible:\n# of descriptors are different: {} and {}\n'.format(TrRawCSV,TsRawCSV,len(Tr.columns),len(Ts.columns)))
    exit(1)
    
for j,desc in enumerate(Tr.columns):
    if desc != Ts.columns[j]:
        sys.stderr.write('error: {} and {} are not compatible.\nthe names for descriptor {} are different: {} and {}\n'.format(TrRawCSV,TsRawCSV,j,desc,Ts.columns[j]))
        exit(1)

outputNormalizedCSV(Ts, Min, Max, TsNormCSV)
