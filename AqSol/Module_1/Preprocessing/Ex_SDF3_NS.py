# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 10:47:46 2023

@author: ITSolutions
"""

import csv
import sys
import numpy as np
from pubchempy import get_compounds, get_sdf
from openbabel import openbabel

def main(argv):
    csv_file = argv[1]
    input_file = open(csv_file, 'r')
    reader = csv.DictReader(input_file)

    print("Input File:", csv_file)

    sdf_file = (csv_file.split('/')[-1])[:-4] + '.SDF'
    output_sdf_file = open(sdf_file, 'w')

    print("Output SDF File:", sdf_file)

    csv_output_file = (csv_file.split('/')[-1])[:-4] + '_values.csv'
    output_file = open(csv_output_file, 'w', newline='')
    csv_writer = csv.writer(output_file)
    csv_writer.writerow(['CID', 'logS'])

    print("Output CSV File:", csv_output_file)

    norm_Target_file = (csv_file.split('/')[-1])[:-4] + '_norm_values.txt'
    norm_file = open(norm_Target_file, 'w', newline='')
    norm_writer = csv.writer(norm_file)
    norm_writer.writerow(['CID', 'a'])

    print("Normalized CSV File:", norm_Target_file)

    sdfs = []
    CID = []
    Target = []
    norm_Target = []
    num = 0  # initialize to calculate number of compounds sdf file

    for row in reader:
        input_smiles = row.get('Smiles')
        input_name = row.get('Name')
        input_target = row.get('logS')

        if input_smiles:
            try:
                mol = openbabel.OBMol()

                # Read the SMILES string
                conv = openbabel.OBConversion()
                conv.SetInFormat("smi")
                conv.ReadString(mol, input_smiles)

                # Convert to canonical SMILES
                conv.SetOutFormat("can")
                canonical_smiles = conv.WriteString(mol).strip()

                cids = get_compounds(canonical_smiles, 'smiles')
                cids = cids[0]  # Get the first compound

                sdf = get_sdf(cids.cid, 'cid')
                sdfs.append(sdf)
                num += 1

                tar = float(input_target)  # converting logS value to floating points
                csv_writer.writerow([cids.cid, tar])

                # Write CID and logS to the CSV file
                CID.append(cids.cid)
                Target.append(tar)

            except:
                print("Compound not found:", input_smiles)
                #print("Error:", str(e))
                continue

        elif input_name:
            try:
                cids = get_compounds(input_name, 'name')
                cids = cids[0]  # Get the first compound

                canonical_smiles = cids.canonical_smiles

                sdf = get_sdf(cids.cid, 'cid')
                sdfs.append(sdf)
                num += 1

                tar = float(input_target)  # converting logS value to floating points
                csv_writer.writerow([cids.cid, tar])

                # Write CID and logS to the CSV file
                CID.append(cids.cid)
                Target.append(tar)

            except:
                print("Compound not found:", input_name)
                #print("Error:", str(e))
                continue

    combined_sdf = ''.join(sdfs)
    output_sdf_file.write(combined_sdf)

    max_value = max(Target)
    min_value = min(Target)
    for value in Target:
        norm_Targets = (value - min_value) / (max_value - min_value)
        norm_Target.append(norm_Targets)

    for i in range(len(CID)):
        val1 = CID[i]
        val2 = norm_Target[i]
        norm_writer.writerow([val1, val2])
    
    # calculating the number of compounds in the input CSV file
    input_file.seek(0)
    num_compounds = sum(1 for _ in input_file) - 1

    print("Number of compounds in input:", num_compounds)
    print("Number of compounds in output:", num)
    print("=================================")

    input_file.close()
    output_sdf_file.close()
    output_file.close()
    norm_file.close()

main(sys.argv)
