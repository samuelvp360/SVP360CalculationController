import openbabel as ob
import subprocess
import re
import os
from PIL import Image, ImageQt
import cairosvg
from io import BytesIO
import numpy as np


class Gaussian():
		
	def __init__(self,molecule):
		ob.obErrorLog.SetOutputLevel(0) # evita los mensajes adicionales no cruciales mientras se usa openbabel
		self._path = molecule
		self.__properties = {}
		self.__calculations = {}
		self.__properties['NAME'] = ''
		self.mol = ob.OBMol()
		self.conversion = ob.OBConversion()
		self.conversion.ReadFile(self.mol,molecule)

		self.lightAtoms = []
		self.heavyAtoms = []
		self.__properties['INIT COORD'] = {}
		self.__properties['MOLAR MASS'] = round(self.mol.GetMolWt(),2)
		self.__properties['FORMULA'] = self.mol.GetFormula()
		self.__properties['SMILES'] = ''

		self.conversion.SetOutFormat('inchikey')
		self.__properties['INCHIKEY'] = re.sub('\n','',re.sub('-','_',self.conversion.WriteString(self.mol)))
		self.__properties['2D IMAGE'] = ''
		self.SetName(molecule)
		self.SetSmiles(molecule)
		self.SetLightAndHeavyAtoms()
		self.Set2DImage()
		
	"""SETTER'S"""
	def SetName(self, name=None):
		if name != None and len(str(name).split(' ')) < 2:
			self.__properties['NAME'] = name
		else:
			self.__properties['NAME'] = self._path.split('/')[-1].split('.')[0]

	def SetSmiles(self, name=None):
		if name != None and len(name.split(' ')) == 1:
			self.conversion.SetOutFormat('smi')
			self.__properties['SMILES'] = self.conversion.WriteString(self.mol).split()[0]+'\t'+name
		else:
			self.conversion.SetOutFormat('smi')
			self.__properties['SMILES'] = self.conversion.WriteString(self.mol).split()[0]+'\t'+self.__properties['NAME']

	def SetLightAndHeavyAtoms(self):		
		for atom in ob.OBMolAtomIter(self.mol):
			at=ob.OBMol()
			a=at.NewAtom()
			a.SetAtomicNum(atom.GetAtomicNum())
			if at.GetMolWt() < 22.9 and self.lightAtoms.count(at.GetFormula()) == 0: # change here the criteria of a heavy atom
				self.lightAtoms.append(at.GetFormula())
			if at.GetMolWt() > 22.9 and self.heavyAtoms.count(at.GetFormula()) == 0:
				self.heavyAtoms.append(at.GetFormula())
			del at

	def Set2DImage(self):
		self.mol.SetTitle('')
		self.conversion.SetOutFormat('svg')
		self.conversion.WriteFile(self.mol,'image.svg')
		out = BytesIO()

		cairosvg.svg2png(url='image.svg', write_to=out, output_width=500, output_height=500)
		os.remove('image.svg')
		PILImage = Image.open(out)
		self.__properties['2D IMAGE'] = ImageQt.ImageQt(PILImage)
		del out, PILImage

	"""GETTER'S"""
	def GetProperties(self):
		return self.__properties

	def GetName(self):
		return self.__properties['NAME']
			
	def GetMolarMass(self):
		return self.__properties['MOLAR MASS']
			
	def GetForm(self):
		return self.__properties['FORMULA']
			
	def GetSmiles(self):
		return self.__properties['SMILES']
			
	def GetInchikey(self):
		return self.__properties['INCHIKEY']

	def Get2DImage(self):
		return self.__properties['2D IMAGE']

	def __str__(self):
		return '{}\n{}\n{}\n{}\n{}\n{}'.format('NAME: '+self.__properties['NAME'],'FORMULA: '+self.__properties['FORMULA'],
			'MOLAR MASS: '+str(self.__properties['MOLAR MASS'])+' g/mol','SMILES: '+self.__properties['SMILES'],'INCHIKEY: '+
			self.__properties['INCHIKEY'], self.__properties['OPT COORDS'])
		
	def __Inputs(self, kind, functional, basis, basis2=None):
		"""Generate the Gaussian inputs for different procedures"""

		with open('/proc/meminfo','r') as f:
			meminfo=f.read()
			matched=re.search(r'^MemTotal:\s+(\d+)',meminfo)
			if matched:
				memory=int(int(matched.groups()[0])/1e3-2e3)

		with open('/proc/cpuinfo','r') as f:
			cpuinfo=f.read()
			matched=re.search(r'siblings\t:\s\d',cpuinfo)
			cpu=str(int(matched[0].split()[2])-1)

		self.functional = functional
		self.basis = basis
		self.basis2 = basis2
		
		if kind=='OPT':
			self.conversion.SetInAndOutFormats('mol2','gzmat')		
			self.mol.SetTotalCharge(0)
			self.mol.SetTotalSpinMultiplicity(1)
			self.gzmat=self.conversion.WriteString(self.mol).splitlines()
			self.gzmat.insert(0,'%Mem='+str(memory)+'MB')
			self.gzmat[1]='%NProcShared='+str(cpu)
			
				
			if self.functional=='MPWB1K':
				if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
					self.gzmat[2]=('#P mpwb95/GENECP IOp(3/76=0560004400) Opt')
				else:
					self.gzmat[2]=(f'#P mpwb95/{basis} IOp(3/76=0560004400) Opt')

			elif self.functional=='M06-2X':
				if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
					self.gzmat[2]=('#P mpwb95/GENECP IOp(3/76=0560004400) Opt')
				else:
					self.gzmat[2]=(f'#P mpwb95/{basis} integral=ultrafine Opt')
			else:
				if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
					self.gzmat[2]=(f'#P {self.functional}/GENECP Opt')
				else:
					self.gzmat[2]=(f'#P {self.functional}/{self.basis} Opt')
			
			self.gzmat[4]=' '+self.__properties['NAME']+'\n'

			if self.basis2 == 'LANL2DZ' or self.basis2 == 'SDD':
				self.gzmat.append(' '.join(self.heavyAtoms)+' 0\n'+basis2+'\n****\n'+' '.join(self.lightAtoms)+\
					' 0\n'+basis+'\n****\n\n'+' '.join(self.heavyAtoms)+' 0\n'+basis2)
			
			with open('Opt_'+self.__properties['NAME']+'.com','w') as self.optFile:
				self.optFile.write('\n'.join(self.gzmat))
	
	def __Run(self, inputFile):
		"""
		

		Parameters
		----------
		inputFile : the .com file as Gaussian input

		Returns
		-------
		.log file corresponding to the output of required calculation

		"""

		homeDir = '/home/'+subprocess.run('ls /home/', shell=True, text=True, capture_output=True).stdout.replace('\n', '')

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
		self.__properties['OPT COORDS']
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
				self.__properties['OPT COORDS'] = np.array(opt)
				del(opt)
		os.remove(f"Opt_{name}.log")

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
					fields[2] = str(coordinatesMatrix[i,0])
					fields[3] = str(coordinatesMatrix[i,1])
					fields[4] = str(coordinatesMatrix[i,2])
					i += 1
					optimizedMol2File.append(
						'      {} {}          {}   {}   {} {}     {}  {}        {}\n'.format(
							fields[0], fields[1], fields[2], fields[3], fields[4],
							fields[5], fields[6], fields[7], fields[8]
							)
						)
				else:
					optimizedMol2File.append(line)

		with open(f'{name}.mol2', 'w') as outputFile: # preguntar si quiero
		# guardar este fichero
			outputFile.writelines(optimizedMol2File)




	def Optimization(self, functional, basis, basis2=None):
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
		if basis2 is not None:
			self.__Inputs('OPT', functional, basis, basis2)
		else:
			self.__Inputs('OPT', functional, basis)

		self.__Run(f"Opt_{self.__properties['NAME']}.com")

		self.__ExtractCoordinates(self.__properties['NAME'])

		self.__ReplaceCoordinates(self.__properties['NAME'], self.__properties['OPT COORDS'])



	def SPEN(self):
		pass

	def SPERA(self):
		pass

	def SPERC(self):
		pass

	def ReactivityIndices(self):
		pass

