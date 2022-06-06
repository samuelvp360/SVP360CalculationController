#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import  pandas as pd
from io import StringIO


path = '/home/samuelvip/Documentos/Bayreuth/IR/JAN2IM.csv'
with open(path, 'r') as file:
    lines = file.readlines()
    new_lines = []
    for line in lines:
        if line[0].isdigit():
            new_lines.append(line)

    csv_str = ''.join(new_lines)

try:
    df = pd.read_csv(StringIO(csv_str), dtype=float, header=None)
except ValueError:
    df = pd.read_csv(
        StringIO(csv_str), delimiter=';',
        decimal=',', dtype=float,
        names=('Wavenumber', 'Raw Data')
    )

print(df)
print(df.dtypes)
