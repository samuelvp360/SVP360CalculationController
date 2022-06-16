#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import  pandas as pd
from io import StringIO

arr = np.array([
    1., 2., 3., 4., 5., 6., 7., 8., 9., 10.
])

print(np.searchsorted(arr, 6.23123))
