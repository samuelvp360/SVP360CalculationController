#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import shutil
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D
from rdkit import RDLogger
from openbabel import openbabel as ob
from persistent import Persistent
from vina import Vina
from datetime import datetime
ob.obErrorLog.StopLogging()  # se puede cambiar por SetOutputLevel para logging


class Molecule(Persistent):

    def __init__(self, path=None, file_format=None, smiles=None):
        RDLogger.DisableLog('rdApp.*')  # evita los mensajes adicionales no cruciales
        if smiles:
            if len(smiles.split(' ')) == 2:
                smiles, name = smiles.split(' ')
            else:
                name = 'NN'
            self.mol = self.__from_smiles(smiles)
            self.set_name(name)
        else:
            self.mol = self.__from_path(path, file_format)
            self.set_name(path, init=True)
        self.calculations = []
        if self.mol:
            self.inchi_key = Chem.MolToInchiKey(self.mol)
            self.smiles = Chem.MolToSmiles(self.mol)
            self.MW = Chem.rdMolDescriptors.CalcExactMolWt(self.mol)
            self.formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
            self.Rg = 0.
            self.__set_mol_picture()

    def __from_smiles(self, smiles):
        converter = ob.OBConversion()
        converter.SetInAndOutFormats('smi', 'can')
        mol_ob = ob.OBMol()
        converter.ReadString(mol_ob, smiles)
        new_block = converter.WriteString(mol_ob)
        return Chem.MolFromSmiles(new_block)

    def __from_path(self, path, file_format):
        with open(path, 'r') as file:
            block = file.read()
            converter = ob.OBConversion()
            converter.SetInAndOutFormats(file_format, 'can')
            mol_ob = ob.OBMol()
            converter.ReadString(mol_ob, block)
            new_block = converter.WriteString(mol_ob)
            return Chem.MolFromSmiles(new_block)

    def __set_mol_picture(self):
        AllChem.Compute2DCoords(self.mol)
        self.mol_pic = f'molecules/{self.inchi_key}/{self.inchi_key}.png'
        if not os.path.exists(f'molecules/{self.inchi_key}/'):
            os.makedirs(f'molecules/{self.inchi_key}/')
        Draw.MolToFile(self.mol, self.mol_pic, size=(300, 300))

    def dock_prep(self, center, receptor_name):
        half = self.mol.GetNumAtoms() // 2
        mol_ref = Chem.MolFromSmiles('C')
        AllChem.EmbedMultipleConfs(mol_ref, 1)
        conf = mol_ref.GetConformer(0)
        conf.SetAtomPosition(0, Point3D(*center))
        rdMolAlign.AlignMol(self.mol, mol_ref, atomMap=[(half, 0)])
        converter = ob.OBConversion()
        mol_ob = ob.OBMol()
        converter.SetInAndOutFormats('pdb', 'mol2')
        file_name = f'molecules/{self.inchi_key}/{self.inchi_key}.mol2'
        converter.ReadString(mol_ob, Chem.MolToPDBBlock(self.mol, confId=0))
        converter.WriteFile(mol_ob, file_name)
        new_name = f'{file_name.split(".")[0]}_{receptor_name}.pdbqt'
        subprocess.run(
            f'obabel {file_name} -O {new_name} -p 7.4 -xp',
            shell=True,
            capture_output=True,
            text=True
        )
        # os.remove(file_name)
        return new_name

    def set_conformer(self, conf):
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
        return conf

    def calculate_Rg(self):
        conf = self.__minimize()
        centroid = Chem.rdMolTransforms.ComputeCentroid(conf)
        centroid_coords = centroid.x, centroid.y, centroid.z
        new_mol = Chem.AddHs(self.mol)
        AllChem.EmbedMolecule(new_mol, randomSeed=0xf00d)
        AllChem.EmbedMultipleConfs(new_mol, 1)
        new_mol.AddConformer(conf, assignId=0)
        Rg = 0.
        num_atoms = 0
        for i, atom in enumerate(new_mol.GetAtoms()):
            if not atom.GetAtomicNum() == 1:
                position = conf.GetAtomPosition(i)
                Rg += (position.x - centroid[0]) ** 2 + (position.y - centroid[1]) ** 2 + (position.z - centroid[2]) ** 2
                num_atoms += 1
        return np.sqrt(Rg / num_atoms), conf

    def set_Rg(self, Rg):
        self.Rg = Rg
        self._p_changed = True

    def set_name(self, name, init=False):
        if init:
            self.__name = name.split('/')[-1].split('.')[0]
        else:
            self.__name = name.strip()
            self._p_changed = True
        if self.mol:
            self.mol.SetProp('_Name', self.__name)
        else:
            return False

    def add_calculation(self, job_id):
        self.calculations.append(job_id)
        self._p_changed = True

    def remove_calculation(self, job_id):
        self.calculations.remove(job_id)
        self._p_changed = True

    @property
    def get_name(self):
        return self.__name

    @property
    def get_coordinates(self): # arreglar
        new_mol = self.__minimize()
        conf = new_mol.GetConformer(0)
        num_atoms = new_mol.GetNumAtoms()
        atoms_coords = np.empty((num_atoms, 3))
        for i, atom in enumerate(new_mol.GetAtoms()):
            position = conf.GetAtomPosition(i)
            atoms_coords[i] = np.array((position.x, position.y, position.z))
        return atoms_coords

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
        elif status == 'Finished':
            self.opt_coords = self.__get_coordinates()
            self.mulliken_charges = self.__get_mulliken()
            self.homo = self.__get_homo()
            self.lumo = self.__get_lumo()
            self.finished = datetime.now()
            self.elapsed = self.finished - self.started
        else:
            self.finished = datetime.now()
            self.elapsed = self.finished - self.started
        # self._p_changed = True

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


class Project(Persistent):

    def __init__(self, name):
        self.name = name
        self.molecules = []
        self.calculations = []
        self.job_ids = []
        self.status = 'Pending'
        self.grid_img = f'projects/{self.name}/{self.name}.png'
        self.__create_grid_img()

    def __create_grid_img(self):
        if not os.path.exists(f'projects/{self.name}/'):
            os.makedirs(f'projects/{self.name}/')
        if self.molecules:
            img = Draw.MolsToGridImage(
                [Chem.MolFromSmiles(m.smiles) for m in self.molecules],
                molsPerRow=3,
                legends=[m.get_name for m in self.molecules]
            )
            img.save(self.grid_img)

    # def set_status(self, status):
        # self.status = status
        # self._p_changed = True

    def add_molecule(self, molecule):
        self.molecules.append(molecule)
        self.__create_grid_img()
        self._p_changed = True

    def remove_molecule(self):
        self.molecules.pop()
        self.__create_grid_img()
        self._p_changed = True

    def add_calculation(self, calc_data):
        self.calculations.append(calc_data)
        self._p_changed = True

    def add_job_id(self, job_id):
        self.job_ids.append(job_id)
        self._p_changed = True

    def remove_calculation(self):
        self.calculations.pop()
        self._p_changed = True

    def run_protocol(self):
        pass


class Docking(Persistent):

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.energies = []
        self.programmed = datetime.now()
        self.__status = 'Programmed'

    def set_status(self, status):
        self.__status = status
        if status == 'Running':
            self.started = datetime.now()
        elif status in 'Finished':
            self.energies = self.get_binding_energies(self.output_file)
            if self.redocking:
                self.rmsd = self.get_rmsd(self.nat_lig_path)
                print(self.rmsd)
            print(self.energies)
            self.finished = datetime.now()
            self.elapsed = self.finished - self.started
            self._p_changed = True
        else:
            self.finished = datetime.now()
            self.elapsed = self.finished - self.started

    @property
    def get_status(self):
        return self.__status

    def get_rmsd(self, nat_lig_path):
        lig_format = nat_lig_path.split('.')[-1]
        converter = ob.OBConversion()
        rmsd = np.empty(len(self.output_file))
        with open(nat_lig_path, 'r') as file:
            nat_lig_block = file.read()
        # if lig_format == 'pdb':
        converter.SetInAndOutFormats(lig_format, 'mol2')
            # mol_ref = Chem.rdmolfiles.MolFromPDBBlock(nat_lig_block)
        # elif lig_format == 'pdbqt':
            # converter.SetInAndOutFormats('pdbqt', 'mol2')
        mol_ob = ob.OBMol()
        converter.ReadString(mol_ob, nat_lig_block)
        # mol_ob.AddHydrogens()
        new_block = converter.WriteString(mol_ob)
        # new_name = nat_lig_path.split(".")[0] + '.mol2'
        # subprocess.run(
            # f'obabel {nat_lig_path} -O {new_name}',
            # shell=True,
            # capture_output=True,
            # text=True
        # )
        # with open(new_name, 'r') as file:
            # block = file.read()
        # mol_ref = Chem.MolFromMol2Block(block)
        # os.remove(new_name)
        mol_ref = Chem.MolFromMol2Block(new_block)
        mol_ref_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol_ref))])))[1]
        mol_ref = Chem.RenumberAtoms(mol_ref, mol_ref_neworder)
        # print(Chem.MolToMolBlock(mol_ref))
        if mol_ref:
            output_files = self.output_file
            for i, out in enumerate(output_files):
                with open(out, 'r') as file:
                    block = file.read()
                    mol_ob = ob.OBMol()
                    converter.SetInAndOutFormats('pdbqt', 'mol2')
                    converter.ReadString(mol_ob, block)
                    # mol_ob.AddHydrogens()
                    new_block = converter.WriteString(mol_ob)
                    mol = Chem.MolFromMol2Block(new_block)
                    mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]
                    mol = Chem.RenumberAtoms(mol, mol_neworder)
                    # print(Chem.MolToMolBlock(mol))
                    try:
                        rmsd[i] = rdMolAlign.GetBestRMS(mol, mol_ref)
                    except RuntimeError:
                        print('An error occurred. Try again.')
                        return rmsd
            return rmsd

    def get_binding_energies(self, outputs):
        energies = np.zeros(len(outputs))
        for index, out in enumerate(outputs):
            with open(out, 'r') as file:
                for line in file.readlines():
                    if line.startswith('REMARK VINA RESULT:'):
                        energies[index] = float(line.split()[3])
                        break
        return energies

    def create_config(self):
        conf = ''
        with open(f'{self.config.get("ligand").split(".")[0]}_conf.txt', 'w') as file:
            for k, v in self.config.items():
                conf += f'{k} = {v}\n'
            file.write(conf)
        # input_file = []
        # template = 'molecules/{0}/{0}_{1}_out_{2}.conf'
        # receptor_name = self.config.get('receptor').split('/')[-1].split('_receptor')[0]
        # previous = len(
            # glob(
                # template.format(
                    # self.molecule.inchi_key,
                    # receptor_name,
                    # '*'
                # )
            # )
        # )
        # for i in range(previous, previous + self.times):
            # output_file.append(
                # template.format(
                    # self.molecule.inchi_key,
                    # receptor_name,
                    # i
                # )
            # )
        # self.config.update({'out': tuple(output_file)})




