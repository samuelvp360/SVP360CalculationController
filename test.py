#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D
from rdkit import RDLogger
from openbabel import openbabel as ob


path_1 = 'ogk.pdb'
# path_1 = 'molecules/FMGBNISRFNDECK-OQXSLNJDSA-N/FMGBNISRFNDECK-OQXSLNJDSA-N_COI1_JAZ_degron_receptor.pdbqt'
path_2 = 'molecules/FMGBNISRFNDECK-OQXSLNJDSA-N/FMGBNISRFNDECK-OQXSLNJDSA-N_COI1_JAZ_degron_receptor_{}.pdbqt'


def mol_from_file(path):
    file_format = path.split('.')[-1]
    with open(path, 'r') as file:
        converter = ob.OBConversion()
        converter.SetInAndOutFormats(file_format, 'mol2')
        block = file.read()
        mol_ob = ob.OBMol()
        converter.ReadString(mol_ob, block)
        new_block = converter.WriteString(mol_ob)
        mol = Chem.MolFromMol2Block(new_block)
        # print(new_block)
        return mol
        # conf = mol.GetConformer(0)
        # num_atoms = mol.GetNumAtoms()
        # atoms_coords = np.empty((num_atoms, 3))
        # for i, atom in enumerate(mol.GetAtoms()):
            # position = conf.GetAtomPosition(i)
            # atoms_coords[i] = np.array((
                # position.x, position.y, position.z
            # ))
        # print(atoms_coords)

mol_ref = mol_from_file(path_1)
mol_ref_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol_ref))])))[1]
mol_ref = Chem.RenumberAtoms(mol_ref, mol_ref_neworder)
rmsd = np.empty(20)
for i in range(1, 21):
    mol_prob = mol_from_file(path_2.format(i))
    mol_prob_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol_prob))])))[1]
    mol_prob = Chem.RenumberAtoms(mol_prob, mol_prob_neworder)
    rmsd[i - 1] = rdMolAlign.CalcRMS(mol_prob, mol_ref)
# print(Chem.MolToSmiles(mol_ref), Chem.MolToSmiles(mol_prob))
# try:
# except RuntimeError:
    # print('No hay subestructura común entre la prueba y la referencia')
print(rmsd)
