#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# import os
# import subprocess
# myEnv = os.environ.copy()
# gaussExecutable = myEnv['GAUSS_EXEDIR'].split('/')[-1]
# # print(myEnv)
# subprocess.run(f'{gaussExecutable} < Opt_TCE.com > Opt_TCE.log', shell=True, env=myEnv)
from Computational_calculations.Gaussian import Gaussian
TCE = Gaussian('SVP-POV_IMINE2-InCl3.mol2')
print(
    f'''{' '.join(set(TCE.heavyAtoms))} 0
LANL2DZ
****
{' '.join(set(TCE.lightAtoms))} 0
6-31G(d)
****

{' '.join(TCE.heavyAtoms)} 0
LANL2DZ
'''
)
