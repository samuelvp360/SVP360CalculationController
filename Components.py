#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import shutil
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit import RDLogger
from openbabel import openbabel as ob
from persistent import Persistent
from datetime import datetime


class Molecule(Persistent):

    def __init__(self, path, file_format):
        RDLogger.DisableLog('rdApp.*')  # evita los mensajes adicionales no cruciales
        self.mol = self.__create_mol(path, file_format)
        self.calculations = []
        if self.mol:
            self.set_name(path, init=True)
            self.inchi_key = Chem.MolToInchiKey(self.mol)
            self.smiles = Chem.MolToSmiles(self.mol)
            self.MW = Chem.rdMolDescriptors.CalcExactMolWt(self.mol)
            self.formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
            self.__set_mol_picture()
            self.__create_conformers()

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

    def __set_mol_picture(self):
        AllChem.Compute2DCoords(self.mol)
        self.mol_pic = f'molecules/{self.inchi_key}/{self.inchi_key}.png'
        if not os.path.exists(f'molecules/{self.inchi_key}/'):
            os.makedirs(f'molecules/{self.inchi_key}/')
        Draw.MolToFile(self.mol, self.mol_pic, size=(200, 200))

    def set_name(self, name, init=False):
        if init:
            self.__name = name.split('/')[-1].split('.')[0]
        else:
            self.__name = name
            self._p_changed = True
        self.mol.SetProp('_Name', self.__name)

    def add_calculation(self, calculation_id):
        self.calculations.append(calculation_id)
        self._p_changed = True

    @property
    def get_name(self):
        return self.__name

    @property
    def get_coordinates(self):
        # hay que agregar manualmente la carga y multiplicidad
        coords = Chem.MolToXYZBlock(self.mol, confId=0).splitlines()
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
        mol_str = Chem.MolToMolBlock(self.mol, confId=0)
        converter.ReadString(mol_ob, mol_str)
        zmat = converter.WriteString(mol_ob).splitlines()
        new_zmat = []
        for line in zmat[6:]:
            if line.startswith('Var'):
                new_zmat.append(' ')
            else:
                new_zmat.append(' ' + line.replace('=', '\t'))
        return '\n'.join(new_zmat)

    @property
    def get_valence_electrons(self):
        return Descriptors.NumValenceElectrons(self.mol)

    @property
    def get_formal_charge(self):
        return Chem.rdmolops.GetFormalCharge(self.mol)

    @property
    def get_heavy_atoms(self):
        return tuple(
            set([a.GetSymbol() for a in self.mol.GetAtoms() if a.GetAtomicNum() > 22.9])
        )

    @property
    def get_ligth_atoms(self):
        return tuple(
            set([a.GetSymbol() for a in self.mol.GetAtoms() if a.GetAtomicNum() <= 22.9])
        )


class Optimization(Persistent):

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.programmed = datetime.now()
        self.__status = 'Programmed'

    def set_status(self, status):
        self.__status = status
        if status == 'Running':
            self.started = datetime.now()
        elif status in 'Finished':
            self.opt_coords = self.__get_coordinates()
            self.mulliken_charges = self.__get_mulliken()
            self.homo = self.__get_homo()
            self.lumo = self.__get_lumo()
            self.finished = datetime.now()
            self.elapsed = self.finished - self.started
            print(self.opt_coords)
        else:
            self.finished = datetime.now()
            self.elapsed = self.finished - self.started
        self._p_changed = True

    @property
    def get_status(self):
        return self.__status

    def __get_coordinates(self):
        with open(self.output_file, 'r') as file:
            log = file.read()
        match_begin = [m for m in re.finditer('Standard orientation:', log)]
        beginning = match_begin[-1].span()[1]
        match_end = [m for m in re.finditer('Rotational constants', log)]
        end = match_end[-1].span()[0]
        pattern = re.compile('\s+\d+\s+\d+\s+\d\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+')
        lines = [i.group() for i in pattern.finditer(log, beginning, end)]
        coords = np.array([line.split() for line in lines])
        return np.delete(coords, (0, 1, 2), 1)

    def __get_mulliken(self):
        pass

    def __get_homo(self):
        pass

    def __get_lumo(self):
        pass

