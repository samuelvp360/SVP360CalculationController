#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import re
import os
import datetime
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
# from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
import pandas as pd
from persistent import Persistent
import transaction


class Gaussian(Persistent):

    def __init__(self, moleculePath):

        RDLogger.DisableLog('rdApp.*')  # evita los mensajes adicionales no cruciales
        self.moleculePath = moleculePath
        self.path = self.moleculePath.replace(' ', '\ ')
        exec(
            f"""self.smiles = subprocess.run(
                'obabel {self.path} -ocan',
                shell=True,
                capture_output=True,
                text=True
            ).stdout.split()[0]"""
        )
        self.__properties = pd.Series({'SMILES': self.smiles}, dtype='string')
        if Chem.MolFromMol2File(self.moleculePath) is not None:
            self.mol = Chem.MolFromMol2File(self.moleculePath)
        else:
            mol = Chem.MolFromSmiles(self.__properties['SMILES'])
            self.mol = Chem.AddHs(mol)

        self.stored = False
        self._p_changed = False
        self.lightAtoms = []
        self.heavyAtoms = []
        self.totalAtoms = []
        self.SetName()
        self.SetInitMatrix()
        self.SetZMatrix()
        self.SetLightAndHeavyAtoms()
        self.__properties['FORMULA'] = ''.join(set((i + str(self.totalAtoms.count(i)) for i in self.totalAtoms)))
        self.__properties['MOLAR MASS'] = round(Descriptors.ExactMolWt(self.mol), 2)
        self.__properties['VALENCE ELECTRONS'] = Descriptors.NumValenceElectrons(self.mol)
        self.__properties['RADICAL ELECTRONS'] = Descriptors.NumRadicalElectrons(self.mol)  # arreglar
        self.__properties['INCHIKEY'] = re.sub('\n', '', re.sub('-', '_', Chem.inchi.MolToInchiKey(self.mol)))
        self.Set2DImage()
        self.__calculations = {}
        # self.__calculations = pd.DataFrame(columns=['JOB TYPE', 'DESCRIPTION', 'RESULTS', 'STATUS'])

    def SetName(self, nameMol=None):

        if nameMol is not None and len(str(nameMol).split(' ')) < 2:
            self.__properties['NAME'] = nameMol
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

    def Set2DImage(self):

        self.mol2D = Chem.MolFromSmiles(self.__properties['SMILES'])
        if self.mol2D is None:
            self.__properties['SMILES'] = Chem.MolToSmiles(self.mol)
            self.mol2D = Chem.MolFromSmiles(self.__properties['SMILES'])
        self.__properties['2D IMAGE'] = Draw.MolToImage(self.mol2D, size=(500, 500))

    def SetInitMatrix(self):
        """
        docstring
        """
        with open(self.moleculePath, 'r') as f:
            lines = f.readlines()
            initialCoords = [i.split()[:] for i in lines if len(i.split()) == 9]
            col = ['ID', 'ATOM', 'X', 'Y', 'Z', 'ATOM_TYPE', 'SUBS_ID', 'SUBS_NAME', 'CHARGE']
            self.__initialMatrix = pd.DataFrame(initialCoords, columns=col)
            self.__initialMatrix[['ID', 'SUBS_ID']] = self.__initialMatrix[['ID', 'SUBS_ID']].astype(int)
            self.__initialMatrix[['X', 'Y', 'Z', 'CHARGE']] = self.__initialMatrix[['X', 'Y', 'Z', 'CHARGE']].astype(float)
            self.__initialMatrix['ATOM'] = self.__initialMatrix['ATOM'].astype(str)
            doubleLetterSymbol = self.__initialMatrix['ATOM'].str.len() == 2
            self.__initialMatrix['ATOM'].loc[doubleLetterSymbol] = [''.join([i[0], i[1].lower()]) for i in self.__initialMatrix['ATOM'] if len(i) == 2]
            self.__initialMatrix.set_index('ID', inplace=True)
            self.atomsList = [i for i in self.__initialMatrix['ATOM']]

    def SetZMatrix(self):
        """
        docstring
        """
        exec(f"""self.run = subprocess.run(
                'obabel {self.path} -ogzmat',
                shell=True,
                capture_output=True,
                text=True).stdout""")
        self.__zMatrix = self.run.splitlines()[6:]

    def SetCalculations(self, jobType, keywords, status, inputFile=None):

        position = len(self.__calculations)
        self.__calculations[position] = {}
        self.__calculations[position]['JOB TYPE'] = jobType
        self.__calculations[position]['DESCRIPTION'] = keywords
        start = datetime.datetime.now()
        if jobType == 'Optimization':
            results = self.__Optimization(inputFile)
            self.__calculations[position]['RESULTS'] = results
        if results is not None:
            self.__calculations[position]['STATUS'] = 'Finished'
        else:
            self.__calculations[position]['STATUS'] = 'Failed'
        end = datetime.datetime.now()
        elapsedTime = end - start
        self.__calculations[position]['DATE OF RUN'] = str(datetime.datetime.now())
        self.__calculations[position]['ELAPSED TIME'] = str(elapsedTime)
        if self.stored:
            self._p_changed = True
            transaction.commit()
#  ----------------------------------------------------GETTERS---------------------------------------
    @property
    def GetProperties(self):
        return self.__properties

    @property
    def GetName(self):
        return self.__properties['NAME']

    @property
    def GetInitMatrix(self):
        return self.__initialMatrix

    @property
    def GetNetCharge(self):
        return int(round(self.__initialMatrix['CHARGE'].sum(), 1))

    @property
    def GetGasteigerCharges(self):
        return self.__initialMatrix['CHARGE']

    @property
    def GetValenceElectrons(self):
        return self.__properties['VALENCE ELECTRONS']

    @property
    def GetRadicalElectrons(self):
        return self.__properties['RADICAL ELECTRONS']

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

    @property
    def GetInitCoords(self):
        return self.__initialMatrix.iloc[:, 0:4]

    @property
    def GetZMatrix(self):
        return '\n'.join(self.__zMatrix)

    @property
    def GetCalculations(self):
        return self.__calculations

    def __str__(self):
        return f"""NAME: {self.__properties['NAME']}
FORMULA: {self.__properties['FORMULA']}
MOLAR MASS: {str(self.__properties['MOLAR MASS'])} g/mol
SMILES: {self.__properties['SMILES']}
INCHIKEY: {self.__properties['INCHIKEY']}
STORED: {self.stored}
{self.__calculations}"""

    # def __Inputs(self, kind, method, functional, basis, basis2, memory, cpu):
    #     """
    #     Parameters
    #     ----------
    #     kind : opt, spen, spera or sperc.
    #     method : DFT or SEMIEMPIRICAL.
    #     functional : a valid functional from the list.
    #     basis : a valid basis from the list.
    #     basis2 : a second basis set in case of heavy atoms.

    #     Returns
    #     -------
    #     an input file with .com extension to run Gaussian
    #     """

    #     self.method = method
    #     self.functional = functional
    #     self.basis = basis
    #     self.basis2 = basis2

    #     if kind == 'opt':

    #         exec(f"""self.run = subprocess.run(
    #             'obabel {self.path} -ogzmat',
    #             shell=True,
    #             capture_output=True,
    #             text=True).stdout""")
    #         self.gzmat = self.run.splitlines()
    #         self.gzmat.insert(0, '%Mem='+str(memory)+'MB')
    #         self.gzmat[1] = '%NProcShared='+str(cpu)

    #         if self.functional == 'MPWB1K':
    #             if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
    #                 self.gzmat[2] = ('#P mpwb95/GENECP IOp(3/76=0560004400) Opt')
    #             else:
    #                 self.gzmat[2] = (f'#P mpwb95/{basis} IOp(3/76=0560004400) Opt')
    #         elif self.functional == 'M06-2X':
    #             if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
    #                 self.gzmat[2] = ('#P mpwb95/GENECP IOp(3/76=0560004400) Opt')
    #             else:
    #                 self.gzmat[2] = (f'#P mpwb95/{basis} integral=ultrafine Opt')
    #         elif self.functional == 'AM1' or self.functional == 'PM3' or self.functional == 'PM6':
    #             self.gzmat[2] = (f'#P {self.functional} Opt')
    #         else:
    #             if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
    #                 self.gzmat[2] = (f'#P {self.functional}/GENECP Opt')
    #             else:
    #                 self.gzmat[2] = (f'#P {self.functional}/{self.basis} Opt')

    #         self.gzmat[4] = ' '+self.__properties['NAME']

    #         if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
    #             self.gzmat.append(
    #                 ' '.join(self.heavyAtoms) + ' 0\n' + basis2 + '\n****\n' +
    #                 ' '.join(self.lightAtoms) + ' 0\n' + basis + '\n****\n\n' +
    #                 ' '.join(self.heavyAtoms) + ' 0\n'+basis2
    #             )

    #         with open('Opt_'+self.__properties['NAME']+'.com', 'w') as optFile:
    #             optFile.write('\n'.join(self.gzmat))
# ----------------------------------------------------METHODS---------------------------------------
    def __Run(self, inputFile, outputFile):
        '''
        Parameters
        ----------
        inputFile : the .com file as Gaussian input

        Returns
        -------
        log file corresponding to the output of required calculation
        '''
        myEnv = os.environ.copy()
        gaussExecutable = myEnv['GAUSS_EXEDIR'].split('/')[-1]
        
        # homeDir = '/home/'+subprocess.run(
        #     'ls /home/', shell=True, text=True, capture_output=True
        # ).stdout.replace('\n', '')

        # with open(f'{homeDir}/.bashrc', 'r') as f:
        #     bashrcInfo = f.read()
        #     root = re.search(r'g[0-9]+root=[/a-zA-Z0-9]+', bashrcInfo)[0]
        #     scrdir = re.search(r'GAUSS_SCRDIR=[/a-zA-Z0-9]+', bashrcInfo)[0]
        #     source = '$' + re.search(r'g[0-9]+root[/a-zA-Z0-9]+\.profile', bashrcInfo)[0]
        #     gaussianExecutable = source.replace('source ', '').replace('bsd/', '').replace('.profile', '')

        # cmd = (f"export {root}; export {scrdir}; source {source}; {gaussianExecutable} < {inputFile} > Opt_{self.__properties['NAME']}.log")
        cmd = f'{gaussExecutable} < {inputFile} > {outputFile}'
        
        subprocess.run(cmd, shell=True)

        # os.remove(inputFile)

    def __ExtractCoordinates(self, outputFile):
        '''
        Parameters
        ----------
        outputFile

        Returns
        -------
        A Pandas DataFrame containing the matrix of optimized coordinates
        '''
        with open(f'{outputFile}', 'r') as logFile:
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
                        opt.append([i.group(0).split()[a] for a in range(6) if a != 2])

        col = ['ID', 'ATOM', 'X', 'Y', 'Z']
        optMatrix = pd.DataFrame(opt, columns=col)
        optMatrix.astype({'ID': int, 'ATOM': str, 'X': float, 'Y': float, 'Z': float})
        optMatrix.set_index('ID', inplace=True)
        optMatrix['ATOM'] = self.atomsList
        return optMatrix

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

    def __Optimization(self, inputFile):
        '''
        Parameters
        ----------
        input: .com file as inputFile

        Returns
        -------
        A Pandas DataFrame with the results
        '''
        outputFile = inputFile.split('.')[0] + '.log'

        self.__Run(inputFile, outputFile)

        coordinates = self.__ExtractCoordinates(outputFile)

        return coordinates

    def SPEN(self):
        pass

    def SPERA(self):
        pass

    def SPERC(self):
        pass

    def ReactivityIndices(self):
        pass
