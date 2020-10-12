#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import numpy as np

# exec(
#         f"""smiles = subprocess.run(
#             'obabel Tetramethylamonium.mol2 -osmi',
#             shell=True,
#             capture_output=True,
#             text=True
#         ).stdout.split()[0]"""
#     )
# print(smiles)
# RDLogger.DisableLog('rdApp.*')
# mol = Chem.MolFromSmiles(smiles)
# print(mol.GetNumAtoms())
# mol.SetProp('_Name', 'Samu')
# print([a.GetAtomicNum() for a in mol.GetAtoms()])
mol = Chem.MolFromSmiles('[N+](C)(C)(C)C')

mol2 = Chem.AddHs(mol)
# Chem.rdPartialCharges.ComputeGasteigerCharges(mol2)
# n = mol.GetNumAtoms()
# n2 = mol2.GetNumAtoms()
# print(n2)
# partialCharges = np.zeros((n2, 1), dtype=np.float32)
# for i in range(n2):
#     a = mol2.GetAtomWithIdx(i)
#     partialCharges[i, 0] = a.GetProp("_GasteigerCharge")
# print(partialCharges)
# netCharge = np.sum(partialCharges, axis=0, dtype=np.float32)[0]
# print(netCharge)
radicalElectrons = Descriptors.NumRadicalElectrons(mol2)
print(radicalElectrons)
valenceElectrons = Descriptors.NumValenceElectrons(mol2)
print(valenceElectrons)
