#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import rdMolAlign, PyMol
from rdkit.Geometry import Point3D
from rdkit import RDLogger
from openbabel import openbabel as ob
import subprocess


smiles = 'CC[C@H](C)[C@H](NC(=O)c1cccc2c1CC[C@@H]2O)C(=O)O'

class Molecule():

    def __init__(self, smiles):
        self.mol = Chem.MolFromSmiles(smiles)
        self.Rg = 0.
        self.descriptors = self.set_descriptors()
        self.morgan_fp = self.get_morgan_fp

    def set_descriptors(self):
        self.__minimize()
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        header = calc.GetDescriptorNames()
        des = calc.CalcDescriptors(self.mol)
        df = pd.DataFrame([des], columns=header)
        for c in df.columns:
            if c.startswith('fr_'):
                df.drop(c, inplace=True, axis=1)
        df['Rg'] = self.Rg if self.Rg else self.calculate_Rg()
        return df

    @property
    def get_morgan_fp(self):
        bit_info = {}
        fp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(
            self.mol, radius=2, nBits=2048, bitInfo=bit_info
        )
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(self.mol))
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        # tpls = [(mol, x, bit_info) for x in fp.GetOnBits()]
        if not os.path.exists('morgan_fp/'):
            os.makedirs('morgan_fp/')
        for bit in bit_info.keys():
            svg_img = Draw.DrawMorganBit(
                mol, bit, bit_info,
                # whichExample=1
                # , molsPerRow=5, legends=[str(x) for x in fp.GetOnBits()],
                # useSVG=True
            )
            with open(f'morgan_fp/{bit}.svg', 'w') as file:
                file.write(svg_img)
        fp = np.array(fp)
        df = pd.DataFrame([fp], columns=[f'{i}' for i in range(2048)])
        return df

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

    def calculate_Rg(self):
        conf = self.mol.GetConformer(0)
        centroid = Chem.rdMolTransforms.ComputeCentroid(conf)
        centroid_coords = centroid.x, centroid.y, centroid.z
        # new_mol = Chem.AddHs(self.mol)
        # AllChem.EmbedMolecule(new_mol, randomSeed=0xf00d)
        # AllChem.EmbedMultipleConfs(new_mol, 1)
        # new_mol.AddConformer(conf, assignId=0)
        Rg = 0.
        num_atoms = 0
        for i, atom in enumerate(self.mol.GetAtoms()):
            if not atom.GetAtomicNum() == 1:
                position = conf.GetAtomPosition(i)
                Rg += (position.x - centroid[0]) ** 2 + (position.y - centroid[1]) ** 2 + (position.z - centroid[2]) ** 2
                num_atoms += 1
        return np.sqrt(Rg / num_atoms)

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
            f'obabel {file_name} -O {new_name} -p 7.4', # cambio aquí
            shell=True,
            capture_output=True,
            text=True
        )


if __name__ == '__main__':
    a_molecule = Molecule(smiles)
    # a_molecule.set_descriptors()
    # a_molecule.dock_prep()
