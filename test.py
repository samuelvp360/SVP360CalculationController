#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import  numpy as np

with open('ogk.pdb', 'r') as file:
    lines = file.readlines()
    coords = np.empty((0, 3), dtype=float)
    for line in lines:
        if line.startswith('HETATM'):
            x, y, z = line.split()[-6:-3]
            if not x.count('.'):
                x, y, z = line.split()[-5:-2]
            coords = np.append(coords, np.array([[float(x), float(y), float(z)]]),
                      axis=0)
    center = tuple(np.average(coords, axis=0))
    print(center)
