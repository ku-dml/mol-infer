# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 03:24:11 2025

@author: Dell Vostro
"""

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

if len(sys.argv) != 3:
    print("Usage: python RDKit-SDF-extraction.py input.csv output.sdf")
    sys.exit(1)

df = pd.read_csv(sys.argv[1])
w = Chem.SDWriter(sys.argv[2])

ok = 0
fail = 0

for i, row in df.iterrows():
    smiles = row["CID"]
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        fail += 1
        continue

    mol = Chem.AddHs(mol)

    # Method 1: ETKDG
    ps = AllChem.ETKDGv3()
    res = AllChem.EmbedMolecule(mol, ps)

    # Method 2: Random coords retry
    if res != 0:
        ps.useRandomCoords = True
        res = AllChem.EmbedMolecule(mol, ps)

    # Method 3: 2D → 3D fallback
    if res != 0:
        print(f"Using fallback 2D→3D embedding for: {smiles}")
        Chem.AllChem.Compute2DCoords(mol)
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol, randomSeed=1)
        except:
            pass

    # Optimize (if possible)
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except:
        print(f"UFF failed for {smiles}, writing unoptimized structure.")

    # Store SDF
    mol.SetProp("_Name", smiles)
    w.write(mol)   # <-- FIXED HERE
    ok += 1

    if i % 1000 == 0:
        print(f"{i} molecules processed...")

w.close()   # <-- FIXED HERE

print("SDF generation complete.")
print("Total:", len(df), "Success:", ok)
