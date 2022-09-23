#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.Chem import rdMolAlign, PyMol
from rdkit.Geometry import Point3D
from rdkit import RDLogger
from openbabel import openbabel as ob
import subprocess


smiles = 'CC[C@H](C)[C@H](NC(=O)c1cccc2c1CC[C@@H]2O)C(=O)O'

class Molecule():

    def __init__(self, smiles):
        self.mol = Chem.MolFromSmiles(smiles)

    def __set_conformer(self, conf):
        self.mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(self.mol, randomSeed=0xf00d)
        AllChem.EmbedMultipleConfs(self.mol, 1)
        self.mol.AddConformer(conf, assignId=0)
        self._p_changed = True

    def __create_conformers(self):
        new_mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(new_mol, randomSeed=0xf00d)
        num_rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(new_mol)
        if num_rot_bonds <= 7:
            n_conf = 50
        elif 8 <= num_rot_bonds <= 12:
            n_conf = 200
        else:
            n_conf = 300
        AllChem.EmbedMultipleConfs(new_mol, n_conf)
        return new_mol

    def __minimize(self):
        new_mol = self.__create_conformers()
        energies = AllChem.MMFFOptimizeMoleculeConfs(
            new_mol, maxIters=2000, nonBondedThresh=100.
        )
        energies_list = [e[1] for e in energies]
        min_e_index = energies_list.index(min(energies_list))
        conf = new_mol.GetConformer(min_e_index)
        self.__set_conformer(conf)

    def dock_prep(self):
        self.__minimize()
        converter = ob.OBConversion()
        mol_ob = ob.OBMol()
        converter.SetInAndOutFormats('pdb', 'mol2')
        file_name = 'mol_prueba.mol2'
        converter.ReadString(mol_ob, Chem.MolToPDBBlock(self.mol, confId=0))
        converter.WriteFile(mol_ob, file_name)
        new_name = f'{file_name.split(".")[0]}.pdbqt'
        subprocess.run(
            f'obabel {file_name} -O {new_name} -p 7.4 -xp', # cambio aquí
            shell=True,
            capture_output=True,
            text=True
        )

    def depict_3D(self):
        self.__minimize()
        self.viewer = PyMol.MolViewer()
        self.viewer.ShowMol(self.mol, name='Samu')

if __name__ == '__main__':
    a_molecule = Molecule(smiles)
# a_molecule.dock_prep()
    a_molecule.depict_3D()
