# -*- coding: utf-8 -*-

import os
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger
from persistent import Persistent


class Molecule(Persistent):

    def __init__(self, path, file_format):
        RDLogger.DisableLog('rdApp.*')  # evita los mensajes adicionales no cruciales
        self.mol = self.__create_mol(path, file_format)
        if self.mol:
            self.set_name(path, init=True)
            self.inchi_key = Chem.MolToInchiKey(self.mol)
            self.smiles = Chem.MolToSmiles(self.mol)
            self.MW = Chem.rdMolDescriptors.CalcExactMolWt(self.mol)
            self.formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
            self.__set_mol_picture()

    def __create_mol(self, path, file_format):
        block = ''
        with open(path, 'r') as file:
            block = file.read()
        if file_format == 'mol':
            return Chem.MolFromMolBlock(block)
        elif file_format == 'mol2':
            return Chem.MolFromMol2Block(block)
        elif file_format == 'pdb':
            return Chem.MolFromPDBBlock(block)

    def __set_mol_picture(self):
        AllChem.Compute2DCoords(self.mol)
        self.mol_pic = 'images/' + self.inchi_key + '.png'
        if not os.path.exists('images/'):
            os.makedirs('images/')
        Draw.MolToFile(self.mol, self.mol_pic)

    def set_name(self, name, init=False):
        if init:
            self.name = name.split('/')[-1].split('.')[0]
        else:
            self.name = name
            self._p_changed = True
        self.mol.SetProp('_Name', self.name)

