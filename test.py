#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Geometry import Point3D

output_file = 'molecules/NLDYACGHTUPAQU-UHFFFAOYSA-N/Optimization_0.log'
with open(output_file, 'r') as file:
    log = file.read()
match_begin = [m for m in re.finditer('Standard orientation:', log)]
beginning = match_begin[-1].span()[1]
match_end = [m for m in re.finditer('Rotational constants', log)]
end = match_end[-1].span()[0]
pattern = re.compile('\s+\d+\s+\d+\s+\d\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+')
lines = [i.group() for i in pattern.finditer(log, beginning, end)]
coords = np.array([line.split() for line in lines], dtype=float)
coords = np.delete(coords, (0, 1, 2), 1)

def create_conformers(mol, num_rot_bonds):
    if num_rot_bonds <= 7:
        n_conf = 50
    elif 8 <= num_rot_bonds <= 12:
        n_conf = 200
    else:
        n_conf = 300
    return AllChem.EmbedMultipleConfs(mol, n_conf)


with open('TCE.mol2', 'r') as file:
    mol2_block = file.read()
mol = Chem.MolFromMol2Block(mol2_block)
num_rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
num_atoms = mol.GetNumAtoms()
conf_ids = create_conformers(mol, num_rot_bonds)
# for i in range(mol.GetNumConformers()):
    # print(mol.GetConformer(i).GetPositions())
conf = mol.GetConformer(0)
for i in range(num_atoms):
    x, y, z = coords[i]
    conf.SetAtomPosition(i, Point3D(x, y, z))
# conf_id = mol.AddConformer(conf, assignId=True)
# for i in range(mol.GetNumConformers()):
    # print(i, mol.GetConformer(i).GetPositions())
# print(Chem.MolToMolBlock(mol, confId=49))
for i in range(mol.GetNumConformers()):
    with open(f'mol_{i}.pdb', 'w') as file:
        file.write(Chem.MolToPDBBlock(mol, confId=i))
# with open('molecules/NLDYACGHTUPAQU-UHFFFAOYSA-N/Optimization_0.log', 'r') as file:
    # log = file.read()
    # matches_begin = [m for m in re.finditer('Standard orientation:', log)]
    # beginning = matches_begin[-1].span()[1]
    # matches_end = [m for m in re.finditer('Rotational constants', log)]
    # end = matches_end[-1].span()[0]
    # pattern = re.compile('\s+\d+\s+\d+\s+\d\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+')
    # lines = [i.group() for i in pattern.finditer(log, beginning, end)]
    # coords = np.array([line.split() for line in lines])
    # print(np.delete(coords, (0, 1, 2), 1))
