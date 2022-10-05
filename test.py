#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import rdMolAlign, PyMol
from rdkit.Geometry import Point3D
from rdkit import RDLogger
from openbabel import openbabel as ob
from loguru import logger
import subprocess


smiles = 'CC/C=C\CN1C(=O)CCN1CC(=O)N[C@H](C(=O)OC)[C@@H](C)CC'

class Molecule():

    def __init__(self, smiles):
        self.conf_description = {}
        self.mol = Chem.MolFromSmiles(smiles)
        self.Rg = 0.
        self.output_file = '/home/samuelvip/Documentos/Python_projects/SVP360CalculationManager/molecules/KHMVBZYLCKUTCC-PFOXVCDASA-N/Optimization_1.log'
        self.__create_conformers()
        self.__minimize()
        # self.descriptors = self.set_descriptors()
        # self.morgan_fp = self.get_morgan_fp

    @logger.catch
    def __str__(self):
        for k, v in self.conf_description.items():
            mol_block = Chem.MolToMolBlock(
                self.mol, confId=v,
            )
            with open(f'a_molecule_{k}.mol', 'w') as file:
                file.write(mol_block)
        return 'done'

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

    def __set_conformer(self, conf, method):
        conf_id = self.mol.AddConformer(conf, assignId=True)
        self.conf_description[method] = conf_id
        print(conf_id)
        # self._p_changed = True

    def __create_conformers(self):
        self.mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(self.mol, randomSeed=0xf00d)
        num_rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(self.mol)
        if num_rot_bonds <= 7:
            n_conf = 50
        elif 8 <= num_rot_bonds <= 12:
            n_conf = 200
        else:
            n_conf = 300
        AllChem.EmbedMultipleConfs(self.mol, n_conf)
        mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(self.mol))])))[1]
        self.mol = Chem.RenumberAtoms(self.mol, mol_neworder)

    def __minimize(self):
        method = 'MMFF94s'
        energies = AllChem.MMFFOptimizeMoleculeConfs(
            self.mol, maxIters=1000, nonBondedThresh=100.,
            mmffVariant=method
        )
        energies_list = [e[1] for e in energies]
        min_e_index = energies_list.index(min(energies_list))
        conf = self.mol.GetConformer(min_e_index)
        self.__set_conformer(conf, method)

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

    def extract_coords(self):
        # try:
        cmd = f'obabel {self.output_file} -omol2'
        out = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            # with open(self.output_file, 'r') as file:
                # log = file.read()
        # except FileNotFoundError:
            # return False
        opt_mol = Chem.MolFromMol2Block(
            out.stdout, removeHs=False,
            # sanitize=False,
        )
        mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(opt_mol))])))[1]
        opt_mol = Chem.RenumberAtoms(opt_mol, mol_neworder)
        opt_conf = opt_mol.GetConformer()
        # match_begin = [m for m in re.finditer('Standard orientation:', log)]
        # beginning = match_begin[-1].span()[1]
        # match_end = [m for m in re.finditer('Rotational constants', log)]
        # end = match_end[-1].span()[0]
        # pattern = re.compile('\s+\d+\s+\d+\s+\d\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+')
        # lines = [i.group() for i in pattern.finditer(log, beginning, end)]
        # coords = np.array([line.split() for line in lines])
        # coords = np.delete(coords, (0, 1, 2), 1)
        # coords = coords.astype(float)
        # conf_id = len(self.conf_description)
        # if not conf_id:
            # self.mol = Chem.AddHs(self.mol)
            # AllChem.EmbedMolecule(self.mol, randomSeed=0xf00d)
            # AllChem.EmbedMultipleConfs(self.mol, 1)
        # conf = self.mol.GetConformer(conf_id)
        # num_atoms = self.mol.GetNumAtoms()
        # for i, atom in enumerate(self.mol.GetAtoms()):
            # conf.SetAtomPosition(i, Point3D(*coords[i]))
        self.__set_conformer(opt_conf, 'PM6')
        # return True


if __name__ == '__main__':
    a_molecule = Molecule(smiles)
    # cmd = f'obabel {a_molecule.output_file} -o{a_molecule.output_file.replace("log", "mol2")}'
    # print(a_molecule)
    a_molecule.extract_coords()
    print(a_molecule)
    # # a_molecule.dock_prep()
