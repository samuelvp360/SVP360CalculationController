# -*- coding: utf-8 -*-

import os
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger
from openbabel import openbabel as ob
from persistent import Persistent


class Molecule(Persistent):

    def __init__(self, path, file_format):
        RDLogger.DisableLog('rdApp.*')  # evita los mensajes adicionales no cruciales
        self.mol = self.__create_mol(path, file_format)
        if self.mol:
            AllChem.EmbedMolecule(self.mol, randomSeed=0xf00d)
            self.set_name(path, init=True)
            self.inchi_key = Chem.MolToInchiKey(self.mol)
            self.smiles = Chem.MolToSmiles(self.mol)
            self.MW = Chem.rdMolDescriptors.CalcExactMolWt(self.mol)
            self.formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
            self.__set_mol_picture()
            self.mol = Chem.AddHs(self.mol)

    def __create_mol(self, path, file_format):
        with open(path, 'r') as file:
            block = file.read()
        if file_format == 'mol':
            return Chem.MolFromMolBlock(block)
        elif file_format == 'mol2':
            mol = Chem.MolFromMol2Block(block)
            if mol:
                return mol
            else:
                converter = ob.OBConversion()
                converter.SetInAndOutFormats('mol2', 'pdb')
                mol_ob = ob.OBMol()
                converter.ReadString(mol_ob, block)
                new_block = converter.WriteString(mol_ob)
                return Chem.MolFromPDBBlock(new_block)
        elif file_format == 'pdb':
            return Chem.MolFromPDBBlock(block)

    def __set_mol_picture(self):
        AllChem.Compute2DCoords(self.mol)
        self.mol_pic = f'molecules/{self.inchi_key}/{self.inchi_key}.png'
        if not os.path.exists(f'molecules/{self.inchi_key}/'):
            os.makedirs(f'molecules/{self.inchi_key}/')
        Draw.MolToFile(self.mol, self.mol_pic, size=(200, 200))

    def set_name(self, name, init=False):
        if init:
            self.name = name.split('/')[-1].split('.')[0]
        else:
            self.name = name
            self._p_changed = True
        self.mol.SetProp('_Name', self.name)

    @property
    def get_name(self):
        return self.name

    @property
    def get_coordinates(self):
        # hay que agregar manualmente la carga y multiplicidad
        coords = Chem.MolToXYZBlock(self.mol).splitlines()
        new_coords = []
        for line in coords[2:]:
            line = line.split()
            line = ' {0:<20}{1:>14}{2:>14}{3:>14}'.format(*line)
            new_coords.append(line)
        return '\n'.join(new_coords)

    @property
    def get_zmatrix(self):
        converter = ob.OBConversion()
        converter.SetInAndOutFormats('mol', 'gzmat')
        converter.AddOption('h')  # mantener los H
        mol_ob = ob.OBMol()
        mol_str = Chem.MolToMolBlock(self.mol)
        converter.ReadString(mol_ob, mol_str)
        zmat = converter.WriteString(mol_ob).splitlines()
        new_zmat = [zmat[5]]
        for line in zmat[6:]:
            if line.startswith('Var'):
                new_zmat.append(' ')
            else:
                new_zmat.append(' ' + line.replace('=', '\t'))
        return '\n'.join(new_zmat)

