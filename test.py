#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import re
# import numpy as np
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from openbabel import openbabel as ob

# output_file = 'molecules/NLDYACGHTUPAQU-UHFFFAOYSA-N/Optimization_0.log'
# with open(output_file, 'r') as file:
    # log = file.read()
# match_begin = [m for m in re.finditer('Standard orientation:', log)]
# beginning = match_begin[-1].span()[1]
# match_end = [m for m in re.finditer('Rotational constants', log)]
# end = match_end[-1].span()[0]
# pattern = re.compile('\s+\d+\s+\d+\s+\d\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+\s+(\d|-\d)+\.\d+')
# lines = [i.group() for i in pattern.finditer(log, beginning, end)]
# coords = np.array([line.split() for line in lines], dtype=float)
# coords = np.delete(coords, (0, 1, 2), 1)

# Eliminar aguas y otros ligandos
subprocess.run(
    f'obabel -ipdb COI1_MODELLER_JAZ_degron.pdb -OCOI1_MODELLER_JAZ_degron.mol2 -h',
    shell=True, capture_output=True, text=True
)
subprocess.run(
    f'obabel -imol2 COI1_MODELLER_JAZ_degron.mol2 -OCOI1_MODELLER_JAZ_degron.pdbqt -p7.4 -xr',
    shell=True, capture_output=True, text=True
)
