#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import SDWriter

converter = ob.OBConversion()
converter.SetInAndOutFormats('mol2', 'pdbqt')
converter.AddOption('h')  # mantener los H
mol_ob = ob.OBMol()

with open('TCE.mol2', 'r') as file:
    mol2_str = file.read()

print(mol2_str)
# mol_rdkit = Chem.MolFromMol2Block(mol2_str)
# mol_rdkit = Chem.AddHs(mol_rdkit)
# AllChem.EmbedMolecule(mol_rdkit, randomSeed=0xf00d)
# AllChem.MMFFOptimizeMolecule(mol_rdkit)
# Chem.rdPartialCharges.ComputeGasteigerCharges(mol_rdkit, nIter=50)
# writer = SDWriter('out.sdf')
# sdf_str = writer.GetText(mol_rdkit)
# print(sdf_str)

converter.ReadString(mol_ob, mol2_str)

output = converter.WriteString(mol_ob)
print(output)

