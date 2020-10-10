#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
# import multiprocessing
# from psutil import virtual_memory
import re
import os
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
# from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
import numpy as np
from persistent import Persistent
import transaction


class Gaussian(Persistent):

    def __init__(self, moleculePath):

        RDLogger.DisableLog('rdApp.*')  # evita los mensajes adicionales no cruciales
        self.path = moleculePath.replace(' ', '\ ')
        self.__properties = {}
        exec(
            f"""self.smiles = subprocess.run(
                'obabel {self.path} -osmi',
                shell=True,
                capture_output=True,
                text=True
            ).stdout.split()[0]"""
        )
        self.__properties['SMILES'] = self.smiles
        if Chem.MolFromMol2File(moleculePath) is not None:
            self.mol = Chem.MolFromMol2File(moleculePath)
        else:
            mol = Chem.MolFromSmiles(self.__properties['SMILES'])
            self.mol = Chem.AddHs(mol)
        
        self.__calculations = {
            'OPT': {}, 'SPEN': {}, 'SPERA': {}, 'SPERC': {}, 'REACTIVITY INDEX': {}
        }

        self.stored = False
        self._p_changed = False
        self.lightAtoms = []
        self.heavyAtoms = []
        self.totalAtoms = []
        self.SetName()
        self.SetCharges()
        self.SetLightAndHeavyAtoms()
        self.__properties['NET CHARGE'] = int(round(np.sum(self.__properties['PARTIAL CHARGES'], axis=0, dtype=np.float32)[0]))
        print(self.__properties['NET CHARGE'])
        self.__properties['FORMULA'] = ''.join(set((i+str(self.totalAtoms.count(i)) for i in self.totalAtoms)))
#         self.__properties['INIT COORD'] = {}
        self.__properties['MOLAR MASS'] = round(Descriptors.ExactMolWt(self.mol), 2)
        self.__properties['INCHIKEY'] = re.sub('\n', '', re.sub('-', '_', Chem.inchi.MolToInchiKey(self.mol)))
        self.Set2DImage()

    def SetName(self, name=None):

        if name is not None and len(str(name).split(' ')) < 2:
            self.__properties['NAME'] = name
            self.mol.SetProp('_Name', self.__properties['NAME'])
            if self.stored:
                self._p_changed = True
                transaction.commit()
        else:
            self.__properties['NAME'] = self.path.split('/')[-1].split('.')[0]
            self.mol.SetProp('_Name', self.__properties['NAME'])

    def SetLightAndHeavyAtoms(self):

        for atom in self.mol.GetAtoms():
            self.totalAtoms.append(atom.GetSymbol())
            if atom.GetAtomicNum() <= 22.9:
                self.lightAtoms.append(atom.GetSymbol())
            else:
                self.heavyAtoms.append(atom.GetSymbol())

        if len(self.heavyAtoms) > 0:
            self.heavyAtomsPresent = True
        else:
            self.heavyAtomsPresent = False

    def SetCharges(self):

        Chem.rdPartialCharges.ComputeGasteigerCharges(self.mol)
        n = self.mol.GetNumAtoms()
        self.__properties['PARTIAL CHARGES'] = np.zeros((n, 1), dtype=np.float32)
        for i in range(n):
            a = self.mol.GetAtomWithIdx(i)
            self.__properties['PARTIAL CHARGES'][i, 0] = a.GetProp("_GasteigerCharge")
        print(self.__properties['PARTIAL CHARGES'])

    def Set2DImage(self):

        self.mol2D = Chem.MolFromSmiles(self.__properties['SMILES'])
        self.__properties['2D IMAGE'] = Draw.MolToImage(self.mol2D, size=(500, 500))

    def SetCalculations(
            self, kind, method, functional, basis, basis2, saveFile, status, memory, cpu
    ):

        if kind == 'opt':
            position = len(self.__calculations['OPT']) + 1
            self.__calculations['OPT'][position] = {}
            self.__calculations['OPT'][position]['METHOD'] = method
            self.__calculations['OPT'][position]['FUNCTIONAL'] = functional
            self.__calculations['OPT'][position]['BASIS'] = basis
            self.__calculations['OPT'][position]['BASIS2'] = basis2
            results = self.__Optimization(method, functional, basis, basis2, memory, cpu)
            self.__calculations['OPT'][position]['RESULT'] = results
            if results is not None:  # estoy hay que optimizarlo
                self.__calculations['OPT'][position]['STATUS'] = 'Finished'
            else:
                self.__calculations['OPT'][position]['STATUS'] = 'Failed'

            if saveFile:
                self.__ReplaceCoordinates(self.__properties['NAME'], results)

        if self.stored:
            print('se da el cambio')  # ojo
            self._p_changed = True
            transaction.commit()

    @property
    def GetProperties(self):
        return self.__properties

    @property
    def GetName(self):
        return self.__properties['NAME']

    @property
    def GetPartialCharges(self):
        return self.__properties['PARTIAL CHARGES']

    @property
    def GetNetCharge(self):
        return self.__properties['NET CHARGE']

    @property
    def GetMolarMass(self):
        return self.__properties['MOLAR MASS']

    @property
    def GetForm(self):
        return self.__properties['FORMULA']

    @property
    def GetSmiles(self):
        return self.__properties['SMILES']

    @property
    def GetInchikey(self):
        return self.__properties['INCHIKEY']

    @property
    def Get2DImage(self):
        return self.__properties['2D IMAGE']

    def GetCalculations(self, kind):

        if kind == 'opt':
            return self.__calculations['OPT']

    def __str__(self):
        return f"""NAME: {self.__properties['NAME']}
FORMULA: {self.__properties['FORMULA']}
MOLAR MASS: {str(self.__properties['MOLAR MASS'])} g/mol
SMILES: {self.__properties['SMILES']}
INCHIKEY: {self.__properties['INCHIKEY']}
STORED: {self.stored}
{self.__calculations}"""

    def __Inputs(self, kind, method, functional, basis, basis2, memory, cpu):
        """
        Parameters
        ----------
        kind : opt, spen, spera or sperc.
        method : DFT or SEMIEMPIRICAL.
        functional : a valid functional from the list.
        basis : a valid basis from the list.
        basis2 : a second basis set in case of heavy atoms.

        Returns
        -------
        an input file with .com extension to run Gaussian
        """

        # mem = virtual_memory()
        # memory = mem.total//1000000 - 2000
        # cpu = multiprocessing.cpu_count() - 1

        self.method = method
        self.functional = functional
        self.basis = basis
        self.basis2 = basis2

        if kind == 'opt':

            exec(f"""self.run = subprocess.run(
                'obabel {self.path} -ogzmat',
                shell=True,
                capture_output=True,
                text=True).stdout""")
            self.gzmat = self.run.splitlines()
            self.gzmat.insert(0, '%Mem='+str(memory)+'MB')
            self.gzmat[1] = '%NProcShared='+str(cpu)

            if self.functional == 'MPWB1K':
                if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
                    self.gzmat[2] = ('#P mpwb95/GENECP IOp(3/76=0560004400) Opt')
                else:
                    self.gzmat[2] = (f'#P mpwb95/{basis} IOp(3/76=0560004400) Opt')
            elif self.functional == 'M06-2X':
                if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
                    self.gzmat[2] = ('#P mpwb95/GENECP IOp(3/76=0560004400) Opt')
                else:
                    self.gzmat[2] = (f'#P mpwb95/{basis} integral=ultrafine Opt')
            elif self.functional == 'AM1' or self.functional == 'PM3' or self.functional == 'PM6':
                self.gzmat[2] = (f'#P {self.functional} Opt')
            else:
                if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
                    self.gzmat[2] = (f'#P {self.functional}/GENECP Opt')
                else:
                    self.gzmat[2] = (f'#P {self.functional}/{self.basis} Opt')

            self.gzmat[4] = ' '+self.__properties['NAME']

            if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
                self.gzmat.append(
                    ' '.join(self.heavyAtoms) + ' 0\n' + basis2 + '\n****\n' +
                    ' '.join(self.lightAtoms) + ' 0\n' + basis + '\n****\n\n' +
                    ' '.join(self.heavyAtoms) + ' 0\n'+basis2
                )

            with open('Opt_'+self.__properties['NAME']+'.com', 'w') as optFile:
                optFile.write('\n'.join(self.gzmat))

    def __Run(self, inputFile):
        """
        Parameters
        ----------
        inputFile : the .com file as Gaussian input

        Returns
        -------
        log file corresponding to the output of required calculation
        """

        homeDir = '/home/'+subprocess.run(
            'ls /home/', shell=True, text=True, capture_output=True
            ).stdout.replace('\n', '')

        with open(f'{homeDir}/.bashrc', 'r') as f:
            bashrcInfo = f.read()
            root = re.search(r'g[0-9]+root=[/a-zA-Z0-9]+', bashrcInfo)[0]
            scrdir = re.search(r'GAUSS_SCRDIR=[/a-zA-Z0-9]+', bashrcInfo)[0]
            source = '$'+re.search(r'g[0-9]+root[/a-zA-Z0-9]+\.profile', bashrcInfo)[0]
            gaussianExecutable = source.replace('source ', '').replace('bsd/', '').replace('.profile', '')

        cmd = (f"export {root}; export {scrdir}; source {source}; {gaussianExecutable} < {inputFile} > Opt_{self.__properties['NAME']}.log")

        subprocess.run(cmd, shell=True, executable='/bin/bash')

        os.remove(inputFile)

    def __ExtractCoordinates(self, name):
        """
        Parameters
        ----------
        name : the name of the molecule

        Returns
        -------
        A numpy object containing the matrix of optimized coordinates
        """
        with open(f"Opt_{name}.log", 'r') as logFile:
            opt = []
            log = logFile.read()
            maxForce = re.findall(r' Maximum Force\s+\d+\.\d+\s+\d+\.\d+\s+(YES|NO)', log)
            rmsForce = re.findall(r' RMS     Force\s+\d+\.\d+\s+\d+\.\d+\s+(YES|NO)', log)
            maxDispl = re.findall(r' Maximum Displacement\s+\d+\.\d+\s+\d+\.\d+\s+(YES|NO)', log)
            rmsDispl = re.findall(r' RMS     Displacement\s+\d+\.\d+\s+\d+\.\d+\s+(YES|NO)', log)

            if maxForce[-1] == 'YES' and rmsForce[-1] == 'YES' and maxDispl[-1] == 'YES' and rmsDispl[-1] == 'YES':

                minSpan = re.search(r'\s+!\s+Optimized Parameters\s+!', log).end()
                for match in re.finditer(r'\s+Distance matrix\s\(angstroms\):', log):
                    maxSpan = match.start()

                optCoord = re.finditer(r'\s+\d+\s+\d+\s+\d\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+', log)
                for i in optCoord:
                    if i.start() > minSpan and i.end() < maxSpan:
                        opt.append(i.group(0).split()[3:6])
        # os.remove(f'Opt_{name}.log')

        return np.array(opt)

    def __ReplaceCoordinates(self, name, coordinatesMatrix):
        """

        Parameters
        ----------
        name : the name of the molecule
        coordinatesMatrix : numpy object containing the matrix of optimized coordinates

        Returns
        -------
        A new {name}.mol2 file with the optimized coordinates

        """

        with open(f'{name}.mol2', 'r') as inputFile:
            lines = inputFile.readlines()
            optimizedMol2File = []
            i = 0
            for line in lines:
                if len(line.split()) == 9:
                    fields = line.split()
                    fields[2] = str(coordinatesMatrix[i, 0])
                    fields[3] = str(coordinatesMatrix[i, 1])
                    fields[4] = str(coordinatesMatrix[i, 2])
                    i += 1
                    optimizedMol2File.append(
                        '      {} {}          {}   {}   {} {}     {}  {}        {}\n'.format(*fields)
                    )
                else:
                    optimizedMol2File.append(line)

        with open(f'Opt_{name}.mol2', 'w') as outputFile:
            # preguntar si quiero guardar este fichero
            outputFile.writelines(optimizedMol2File)

    def __Optimization(self, method, functional, basis, basis2, memory, cpu):
        """
        Parameters
        ----------
        functional : main functional to use.
        basis : main basis for light atoms.
        basis2 : secondary basis if heavy atoms are present. The default is None.

        Returns
        -------
        Opt_{name}.mol2 file

        """
        self.__Inputs('opt', method, functional, basis, basis2, memory, cpu)

        self.__Run(f"Opt_{self.__properties['NAME']}.com")

        result = self.__ExtractCoordinates(self.__properties['NAME'])

        return result

    def SPEN(self):
        pass

    def SPERA(self):
        pass

    def SPERC(self):
        pass

    def ReactivityIndices(self):
        pass
