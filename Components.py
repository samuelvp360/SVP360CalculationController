#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import os
import re
import subprocess
import numpy as np
import pandas as pd
from PIL import Image, ImageQt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolAlign
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Geometry import Point3D
from rdkit import RDLogger
from openbabel import openbabel as ob
from map4 import MAP4Calculator
from persistent import Persistent
from datetime import datetime
import MDAnalysis as mda
import prolif as plf
from requests.exceptions import ConnectionError
try:
    from chembl_webresource_client.new_client import new_client
    conected = True
except ConnectionError:
    conected = False
ob.obErrorLog.StopLogging()  # se puede cambiar por SetOutputLevel para logging


class Molecule(Persistent):

    def __init__(self, mol_input, file_format='smiles'):
        RDLogger.DisableLog('rdApp.*')  # avoid non crucial messages
        self.conf_dict = {}
        self.calculations = []
        self.activities = []
        self.mol = self.get_rdkit_mol(mol_input, file_format)
        if self.mol:
            input_conf = self.mol.GetConformer()
            self.set_conformer(input_conf, 'input')
            self.inchi_key = Chem.MolToInchiKey(self.mol)
            self.smiles = Chem.MolToSmiles(Chem.RemoveHs(self.mol))
            self.MW = Chem.rdMolDescriptors.CalcExactMolWt(self.mol)
            self.formula = Chem.rdMolDescriptors.CalcMolFormula(self.mol)
            self.__set_mol_img()
            # self.descriptors = self.get_descriptors
            # self.fps_dict = self.get_fps_dict

    def __eq__(self, other):
        if isinstance(other, Molecule):
            return self.smiles == other.smiles or \
                    self.inchi_key == other.inchi_key
        return False

    def get_rdkit_mol(self, mol_input, file_format):
        if file_format == 'smiles':
            return self.__from_smiles(mol_input)
        elif file_format in 'mol2':
            mol = self.__from_path(mol_input, file_format)
            conf = self.__from_opt_file(mol_input)
            conf = conf.GetConformer()
            mol.AddConformer(conf)
            return mol
        elif file_format == 'log':
            return self.__from_opt_file(mol_input)
        else:
            return self.__from_path(mol_input, file_format)

    def __canonize(self, mol):
        mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]
        return Chem.RenumberAtoms(mol, mol_neworder) if mol else False

    def __from_smiles(self, smiles):
        if len(smiles.split(' ')) == 2:
            smiles, name = smiles.split(' ')
        else:
            name = 'NN'
        self.set_name(name)
        converter = ob.OBConversion()
        converter.SetInAndOutFormats('smi', 'can')
        mol_ob = ob.OBMol()
        converter.ReadString(mol_ob, smiles)
        new_block = converter.WriteString(mol_ob)
        mol = Chem.MolFromSmiles(new_block)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        AllChem.EmbedMultipleConfs(mol, 1)
        mol = self.__canonize(mol)
        return mol

    def __from_path(self, path, file_format):
        self.set_name(path.split('/')[-1], init=True)
        converter = ob.OBConversion()
        converter.SetInAndOutFormats(file_format, 'can')
        mol_ob = ob.OBMol()
        with open(path, 'r') as file:
            block = file.read()
            converter.ReadString(mol_ob, block)
        new_block = converter.WriteString(mol_ob)
        mol = Chem.MolFromSmiles(new_block)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        AllChem.EmbedMultipleConfs(mol, 1)
        mol = self.__canonize(mol)
        return mol

    def __from_opt_file(self, path, set_name=True):
        if set_name:
            self.set_name(path.split('/')[-1], init=True)
        cmd = f'obabel {path} -omol2 -h'
        out = subprocess.run(
            cmd, shell=True, capture_output=True, text=True
        )
        mol = Chem.MolFromMol2Block(
            out.stdout, removeHs=False,
            sanitize=False,
        )
        # mol = self.__canonize(mol)
        return mol

    @classmethod
    def get_mol_img(self, mol):
        drawer = Draw.MolDraw2DCairo(350, 350)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        data = drawer.GetDrawingText()
        bio = io.BytesIO(data)
        img = Image.open(bio)
        return img

    def __set_mol_img(self):
        self.mol_img = f'molecules/{self.inchi_key}/{self.inchi_key}.png'
        if not os.path.exists(f'molecules/{self.inchi_key}/'):
            os.makedirs(f'molecules/{self.inchi_key}/')
        else:
            self.mol_img = f'molecules/{self.inchi_key}_2/{self.inchi_key}.png'
            os.makedirs(f'molecules/{self.inchi_key}_2/')
            # return  # to avoid redo the process for existing molecules
        Draw.MolToFile(
            Chem.MolFromSmiles(self.smiles), self.mol_img, size=(350, 350)
        )

    def ligand_prep(self, center, receptor_name, conformer):
        half = self.mol.GetNumAtoms() // 2
        mol_ref = Chem.MolFromSmiles('C')
        AllChem.EmbedMultipleConfs(mol_ref, 1)
        conf = mol_ref.GetConformer(0)
        conf.SetAtomPosition(0, Point3D(*center))
        rdMolAlign.AlignMol(
            self.mol, mol_ref, atomMap=[(half, 0)],
            prbCid=self.conf_dict[conformer]
        )
        converter = ob.OBConversion()
        mol_ob = ob.OBMol()
        converter.SetInAndOutFormats('pdb', 'mol2')
        try:
            file_name = f'molecules/{self.inchi_key}/{self.inchi_key}.mol2'
        except FileNotFoundError:
            file_name = f'molecules/{self.inchi_key}_2/{self.inchi_key}.mol2'
        converter.ReadString(
            mol_ob, Chem.MolToPDBBlock(
                self.mol, confId=self.conf_dict[conformer]
            )
        )
        converter.WriteFile(mol_ob, file_name)
        output_filename = f'{file_name.split(".")[0]}_{receptor_name}.pdbqt'
        subprocess.run(
            f'obabel {file_name} -O {output_filename} -p 7.4',
            shell=True, capture_output=True, text=True
        )
        # os.remove(file_name)
        # return output_filename

    def set_conformer(self, conf, method):
        if method in self.conf_dict:
            conf_id = self.conf_dict[method]
            self.mol.RemoveConformer(conf_id)
            conf_id = self.mol.AddConformer(conf, assignId=True)
            self.conf_dict.update({method: conf_id})
        else:
            conf_id = self.mol.AddConformer(conf, assignId=True)
            self.conf_dict[method] = conf_id
        self._p_changed = True

# GETTERS -------------------------------------------------------
    @property
    def get_descriptors(self):
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(
            [x[0] for x in Descriptors._descList]
        )
        header = calc.GetDescriptorNames()
        Chem.SanitizeMol(self.mol)
        des = calc.CalcDescriptors(self.mol)
        descriptors = pd.DataFrame([des], columns=header)
        for c in descriptors.columns:  # to get rid of the fingerprints
            if c.startswith('fr_'):
                descriptors.drop(c, inplace=True, axis=1)
        descriptors.loc[0, 'Rg'] = self.get_Rg
        return descriptors

    @property
    def get_fps_dict(self):
        fps_dict = {}
        mol = Chem.AddHs(self.get_2d_mol, addCoords=True)
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        for i in range(1,4):
            bit_info = {}
            fps_dict[f'Morgan{i}'] = np.array(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol, radius=i, bitInfo=bit_info
            ))
        for i in ('chiral', 'achiral'):
            chirality = True if i == 'chiral' else False
            fps_dict[f'Hashed Atom Pairs ({i})'] = np.array(Chem.rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                mol, includeChirality=chirality
            ))
        rdkit_fp = Chem.RDKFingerprint(mol, maxPath=5)  #, bitInfo=rdkit_bit_info)
        fps_dict['RDKit'] = np.array(rdkit_fp)
        fps_dict['Layered'] = np.array(Chem.rdmolops.LayeredFingerprint(mol))
        fps_dict['Pattern'] = np.array(Chem.rdmolops.PatternFingerprint(mol))
        fps_dict['Hasehd Topological Torsions'] = np.array(Chem.rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol))
        calculator = MAP4Calculator(dimensions=2048, is_folded=True)
        fps_dict['MAP4 (folded)'] = calculator.calculate(mol)
        return fps_dict

    @property
    def get_2d_mol(self):
        return Chem.MolFromSmiles(self.smiles)

    def __create_conformers(self):
        new_mol = Chem.AddHs(self.mol, addCoords=True)
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

    def minimize(self):
        method = 'MMFF94s'
        new_mol = self.__create_conformers()
        energies = AllChem.MMFFOptimizeMoleculeConfs(
            new_mol, maxIters=2000, nonBondedThresh=100.,
            mmffVariant=method
        )
        energies_list = [e[1] for e in energies]
        min_e_index = energies_list.index(min(energies_list))
        conf = new_mol.GetConformer(min_e_index)
        return conf, method

    @property
    def get_Rg(self):
        conf = self.mol.GetConformer(self.conf_dict['input'])
        centroid = Chem.rdMolTransforms.ComputeCentroid(conf)
        Rg = 0.
        num_atoms = 0
        for i, atom in enumerate(self.mol.GetAtoms()):
            if not atom.GetAtomicNum() == 1:
                position = conf.GetAtomPosition(i)
                Rg += (position.x - centroid.x) ** 2 + \
                    (position.y - centroid.y) ** 2 + \
                    (position.z - centroid.z) ** 2
                num_atoms += 1
        return np.sqrt(Rg / num_atoms)

    def set_name(self, name, init=False):
        if init:
            self.__name = name.split('/')[-1].split('.')[0]
        else:
            self.__name = name.strip()
            self._p_changed = True
        # if self.mol:
            # self.mol.SetProp('_Name', self.__name)
        # else:
            # return False

    def add_calculation(self, job_id):
        self.calculations.append(job_id)
        self._p_changed = True

    def add_activity_value(self, activity):
        self.activities.append(activity)
        self._p_changed = True

    def remove_calculation(self, job_id):
        self.calculations.remove(job_id)
        self._p_changed = True

    @property
    def get_name(self):
        return self.__name

    def add_conf_from_opt_file(self, output_file, method):
        # primero hay que checkear si convegió, de lo contrario saltar un aviso
        # recomendando cambiar el sistema de coordenadas
        mol = self.__from_opt_file(output_file, set_name=False)
        if mol:
            init_num_atoms = self.mol.GetNumAtoms()
            new_num_atoms = mol.GetNumAtoms()
            if not new_num_atoms == init_num_atoms:
                print(f'Number of atoms mismatch old: {init_num_atoms} vs new: {new_num_atoms}')
                return
            conf = mol.GetConformer()
            self.set_conformer(conf, method)
        else:
            return  # arrojar un error

    @property
    def get_coordinates(self):
        converter = ob.OBConversion()
        converter.SetInAndOutFormats('mol', 'xyz')
        mol_ob = ob.OBMol()
        mol_str = Chem.MolToMolBlock(self.mol, confId=0)
        converter.ReadString(mol_ob, mol_str)
        converter.AddOption('h')  # mantener los H
        coordinates = converter.WriteString(mol_ob).splitlines()[2:]
        coordinates = [' ' + c for c in coordinates]
        return '\n'.join(coordinates)

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


class Protein(Persistent):

    def __init__(self, path):
        self.prot = self.get_prot(path)
        print(self.prot.n_residues)

    def get_prot(self, path):
        file_format = path.split('.')[-1]
        univ = mda.Universe(path)
        try:
            prot = plf.Molecule.from_mda(univ, force=True)
            return prot
        except:
            file_format = path.split('.')[-1]
            univ = self.fix_kekule(univ)
            if file_format == 'mol2':
                prot = self.fix_aromatic(univ, path)
            else:
                prot = plf.Molecule.from_mda(univ, force=True)
            return prot

    def fix_kekule(self, univ):
        # replace aromatic bonds with single bonds
        for i, bond_order in enumerate(univ._topology.bonds.order):
            # you may need to replace double bonds ("2") as well
            if bond_order == "ar":
                univ._topology.bonds.order[i] = 1
        # clear the bond cache, just in case
        univ._topology.bonds._cache.pop("bd", None)
        # infer bond orders again
        return univ

    def fix_aromatic(self, univ, path):
        mol = Chem.MolFromMol2File(path, removeHs=False)
        for atom, resname in zip(mol.GetAtoms(), univ.atoms.resnames):
            resid = plf.ResidueId.from_string(resname)
            mi = Chem.AtomPDBResidueInfo()
            mi.SetResidueNumber(resid.number)
            mi.SetResidueName(resid.name)
            atom.SetMonomerInfo(mi)
        prot = plf.Molecule(mol)
        return prot


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
            # self.opt_coords = self.__extract_coords()
            self.mulliken_charges = self.__get_mulliken()
            self.homo = self.__get_homo()
            self.lumo = self.__get_lumo()
            self.finished = datetime.now()
            self.elapsed = self.finished - self.started
        else:
            self.finished = datetime.now()
            if hasattr(self, 'started'):
                self.elapsed = self.finished - self.started
            else:
                self.elapsed = ''
        self._p_changed = True

    @property
    def get_status(self):
        return self.__status

    # def __extract_coords(self):
        # try:
            # with open(self.output_file, 'r') as file:
                # log = file.read()
        # except FileNotFoundError:
            # return False
        # match_begin = [m for m in re.finditer('Standard orientation:', log)]
        # beginning = match_begin[-1].span()[1]
        # match_end = [m for m in re.finditer('Rotational constants', log)]
        # end = match_end[-1].span()[0]
        # pattern = re.compile('\s+\d+\s+\d+\s+\d\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+')
        # lines = [i.group() for i in pattern.finditer(log, beginning, end)]
        # coords = np.array([line.split() for line in lines])
        # return np.delete(coords, (0, 1, 2), 1)

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
        self.descriptors = pd.DataFrame([])
        self.fps = []
        self.clusters = {}
        self.grid_img = f'projects/{self.name}/{self.name}.png'
        self.grid_img_svg = f'projects/{self.name}/{self.name}.svg'
        self.__create_grid_img()

    def __eq__(self, other):
        if isinstance(other, Project):
            return self.name == other.name
        return False

    def __create_grid_img(self):
        if not os.path.exists(f'projects/{self.name}/'):
            os.makedirs(f'projects/{self.name}/')
        if self.molecules:
            img = Draw.MolsToGridImage(
                [Chem.MolFromSmiles(m.smiles) for m in self.molecules],
                molsPerRow=5,
                legends=[m.get_name for m in self.molecules]
            )
            img.save(self.grid_img)
            img_svg = Draw.MolsToGridImage(
                [Chem.MolFromSmiles(m.smiles) for m in self.molecules],
                molsPerRow=5, useSVG=True,
                legends=[m.get_name for m in self.molecules]
            )
            with open(self.grid_img_svg, 'w') as file:
                file.write(img_svg)

    def add_molecule(self, molecule):
        self.molecules.append(molecule)
        self.__create_grid_img()
        self._p_changed = True

    def pop_molecule(self):
        self.molecules.pop()
        self.fps.pop()
        self.descriptors.drop(index=self.descriptors.index[-1], inplace=True)
        self.__create_grid_img()
        self._p_changed = True

    def add_calculation(self, calc_data):
        self.calculations.append(calc_data)
        self._p_changed = True

    def add_job_id(self, job_id):
        self.job_ids.append(job_id)
        self._p_changed = True

    def add_fps(self, fps_dict):
        self.fps.append(fps_dict)
        self._p_changed = True

    def add_descriptors(self, descriptors, inchi_key):
        descriptors = descriptors.copy()
        descriptors['inchi_key'] = inchi_key
        descriptors.set_index('inchi_key', inplace=True)
        self.descriptors = pd.concat([self.descriptors, descriptors])
        self.descriptors.replace(0, np.nan, inplace=True)
        self.descriptors.dropna(how='all', axis=1, inplace=True)
        self.descriptors.replace(np.nan, 0, inplace=True)
        self._p_changed = True

    def set_clusters(
        self, cluster_type, labels,
        total, silhouette_avg,
        threshold
    ):
        self.clusters = {
            'type': cluster_type,
            'labels': labels,
            'total': total,
            'silhouette_avg': silhouette_avg,
            'threshold': threshold
        }
        # self._p_changed = True

    def get_cluster(self, cluster_num):
        labels = self.clusters.get('labels')
        mols_in_cluster = [
            m for l, m in zip(labels, self.molecules) if l == cluster_num
        ]
        return mols_in_cluster

    def get_docking_poses(self):
        # primero pasar del pdbqt a pdb y luego asignar los órdenes de enlace
        # según la plantilla, que es la molécula rdkit
        newMol = AllChem.AssignBondOrdersFromTemplate(template, docked_pose)


    def remove_calculation(self):
        if self.calculations:
            self.calculations.pop()
            self._p_changed = True


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


class Similarity():

    def __init__(self, fps, names):
        self.fps = fps
        self.names = names

    def similarity_matrix(self, fp_type, simil_metric):
        df = pd.DataFrame()
        fps = [fp_dict.get(fp_type) for fp_dict in self.fps]
        if simil_metric == 'Tanimoto':
            simil_function = self.tanimoto
        elif simil_metric == 'Dice':
            simil_function = self.dice
        elif simil_metric == 'Cosine':
            simil_function = self.cosine
        elif simil_metric == 'Sokal':
            simil_function = self.sokal
        elif simil_metric == 'Russel':
            simil_function = self.russel
        elif simil_metric == 'Hamann':
            simil_function = self.hamann
        elif simil_metric == 'Kulczynski':
            simil_function = self.kulczynski
        for i_index, i in enumerate(fps):
            for j_index, j in enumerate(fps):
                similarity = simil_function(i, j)
                df.loc[i_index, j_index] = round(similarity, 2)
        df.set_axis(self.names, axis='columns', inplace=True)
        df.set_axis(self.names, axis='index', inplace=True)
        return df

    @classmethod
    def tanimoto(self, i, j):
        a = np.sum(i & j)
        b = np.sum(i) - a
        c = np.sum(j) - a
        d = np.sum((1 - i) & (1 - j))
        return a / (a + b + c)

    def dice(self, i, j):
        a = np.sum(i & j)
        b = np.sum(i) - a
        c = np.sum(j) - a
        d = np.sum((1 - i) & (1 - j))
        return 2 * a / (2 * a + b + c)

    def sokal(self, i, j):
        a = np.sum(i & j)
        b = np.sum(i) - a
        c = np.sum(j) - a
        d = np.sum((1 - i) & (1 - j))
        return a / (a + (2 * (b + c)))

    def russel(self, i, j):
        a = np.sum(i & j)
        b = np.sum(i) - a
        c = np.sum(j) - a
        d = np.sum((1 - i) & (1 - j))
        return a / (a + b + c + d)

    def hamann(self, i, j):
        a = np.sum(i & j)
        b = np.sum(i) - a
        c = np.sum(j) - a
        d = np.sum((1 - i) & (1 - j))
        return (a + d - b - c) / (a + b + c + d)

    def kulczynski(self, i, j):
        a = np.sum(i & j)
        b = np.sum(i) - a
        c = np.sum(j) - a
        d = np.sum((1 - i) & (1 - j))
        return (a / 2) * ((1 / (a + b)) + (1 / (a + c)))

    def cosine(self, i, j):
        return np.dot(i, j) / (np.linalg.norm(i) * np.linalg.norm(j))


class MyChembl():

    def __init__(self):
        self.target = new_client.target if conected else None
        self.activity = new_client.activity if conected else None
        self.molecule = new_client.molecule if conected else None

    def get_target_id(self, target_name):
        target_id = self.target.filter(pref_name__icontains=target_name).only(
            ['organism', 'pref_name', 'target_chembl_id']
        )
        return pd.DataFrame(target_id) if target_id else pd.DataFrame()

    def get_activities_for_a_target(self, target_id):
        columns = [
            'molecule_chembl_id', 'standard_type', 'standard_relation',
            'standard_units', 'value'
        ]
        query = self.activity.filter(target_chembl_id=target_id).only(columns).order_by('type')
        return query

    def get_query_df(self, page):
        df = pd.DataFrame.from_records(page)
        df.dropna(inplace=True, subset=[
            'standard_type', 'standard_units',
            'standard_relation', 'value'
        ])
        df.drop(columns=['relation', 'type', 'units'], inplace=True)
        return df

    def get_similarities(self, chembl_id, fp_method, mol_ref):
        structures = self.molecule.get(chembl_id).get('molecule_structures')
        if structures:
            smiles = structures.get('canonical_smiles')
            mol_prob = Chem.MolFromSmiles(smiles)
            if 'Morgan' in fp_method:
                radius = int(fp_method[-1])
                mol_prob_fp = np.array(
                    Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(
                        mol_prob, radius=radius
                ))
                mol_ref_fp = np.array(
                    Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(
                        mol_ref, radius=radius
                ))
            elif 'Hashed Atom Pairs' in fp_method:
                chirality = True if 'chiral' in fp_method else False
                mol_prob_fp = np.array(
                    Chem.rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                        mol_prob, includeChirality=chirality
                    ))
                mol_ref_fp = np.array(
                    Chem.rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
                        mol_ref, includeChirality=chirality
                    ))
            elif fp_method == 'RDKit':
                mol_prob_fp = np.array(Chem.RDKFingerprint(mol_prob, maxPath=5))
                mol_ref_fp = np.array(Chem.RDKFingerprint(mol_ref, maxPath=5))
            similarity = Similarity.tanimoto(mol_ref_fp, mol_prob_fp)
            result = {
                'similarity': round(similarity, 3),
                'smiles': smiles,
            }
        return result


class Activity(Persistent):

    def __init__(self, act_dict):
        self.id = act_dict['id']
        self.name = act_dict['name']
        self.species = act_dict['species']
        self.unit = act_dict['unit']
        self.values = {}

    def check(self, mol_id):
        return mol_id in self.values

    def add_activity_value(self, mol_id, value, update=False):
        if not update:
            self.values[mol_id] = value
        else:
            self.values.update({mol_id: value})
        self._p_changed = True

