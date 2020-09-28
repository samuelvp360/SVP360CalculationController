#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 15:16:41 2020

@author: samuelvip
"""


myDict = {
    'casa': 123, 'calle': 321, 'apto': 401, 'cuarto': 'el de al lado'
    }


class test():

    def __init__(self):

        self.name = 'Samu'

    def reviewAttr(self, parameters):

#         for key in parameters:
#             setattr(self, key, parameters[key])
#         print(f'''esta en una prueba en la que determinaré si es correcto
# el desempaquetamiento del diccionario,
# el apartamento es el {self.apto},
# la casa es la {self.casa},
# la calle la {self.calle},
# y el cuarto es {self.cuarto}.''')
        parameters.keys() = parameters.items()

prueba = test()
prueba.reviewAttr(myDict)
