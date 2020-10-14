#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re

sss = 'AUTO=3'

sss = re.sub(r'\d+$', '4', sss)

print(sss)
# print(aaa)
