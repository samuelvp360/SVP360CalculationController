#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd


class PeriodicTable():
    """
    docstring
    """
    def __init__(self):
        """
        docstring
        """
        self.periodicTable = pd.DataFrame(
            {
                'Number':
                [
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                    22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                    42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
                    62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,
                    82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
                    102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118
                ],
                'Name':
                [
                    'Neutron', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron', 'Carbon', 'Nitrogen',
                    'Oxygen', 'Fluorine', 'Neon', 'Sodium', 'Magnesium', 'Aluminum', 'Silicon', 'Phosphorus',
                    'Sulfur', 'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium',
                    'Chromium', 'Manganese', 'Iron', 'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
                    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium', 'Strontium', 'Yttrium', 'Zirconium',
                    'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver', 'Cadmium',
                    'Indium', 'Tin', 'Antimony', 'Tellurium', 'Iodine', 'Xenon', 'Cesium', 'Barium', 'Lanthanum',
                    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium', 'Samarium', 'Europium', 'Gadolinium',
                    'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium', 'Lutetium', 'Hafnium',
                    'Tantalum', 'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury', 'Thallium',
                    'Lead', 'Bismuth', 'Polonium', 'Astatine', 'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
                    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium', 'Americium', 'Curium', 'Berkelium',
                    'Californium', 'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium', 'Lawrencium', 'Rutherfordium',
                    'Dubnium', 'Seaborgium', 'Bohrium', 'Hassium', 'Meitnerium', 'Darmstadtium', 'Roentgenium',
                    'Copernicium', 'Nihonium', 'Flerovium', 'Moscovium', 'Livermorium', 'Tennessine', 'Oganesson'
                ],
                'Symbol':
                [
                    'n', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',
                    'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
                    'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
                    'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
                    'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
                    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
                ],
                'Common Ions':
                [
                    [], [-1, 1], [], [1], [2], [3], [-4, -3, -2, -1, 1, 2, 3, 4], [-3, 3, 5], [-2], [-1],
                    [], [1], [2], [3], [-4, 4], [-3, 3, 5], [-2, 2, 4, 6], [-1, 1, 3, 5, 7], [], [1], [2],
                    [3], [4], [5], [3, 6], [2, 4, 7], [2, 3, 6], [2, 3], [2], [2], [2], [3], [-4, 2, 4],
                    [-3, 3, 5], [-2, 2, 4, 6], [-1, 1, 3, 5], [2], [1], [2], [3], [4], [5], [4, 6], [4, 7],
                    [3, 4], [3], [2, 4], [1], [2], [3], [-4, 2, 4], [-3, 3, 5], [-2, 2, 4, 6], [-1, 1, 3, 5, 7],
                    [2, 4, 6], [1], [2], [3], [3, 4], [3], [3], [3], [3], [2, 3], [3], [3], [3], [3], [3], [3],
                    [3], [3], [4], [5], [4, 6], [4], [4], [3, 4], [2, 4], [3], [1, 2], [1, 3], [2, 4], [3],
                    [-2, 2, 4], [-1, 1], [2], [1], [2], [3], [4], [5], [6], [5], [4], [3], [3], [3], [3], [3],
                    [3], [3], [2], [3], [4], [5], [6], [7], [8], [], [], [], [2], [], [], [], [], [], []
                ],
                'Standard Mass':
                [
                    None, 1.0079, 4.0026, 6.941, 9.0122, 10.811, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797, 22.9897,
                    24.305, 26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.0983, 39.948, 40.078, 44.9559, 47.867,
                    50.9415, 51.9961, 54.938, 55.845, 58.6934, 58.9332, 63.546, 65.39, 69.723, 72.64, 74.9216,
                    78.96, 79.904, 83.8, 85.4678, 87.62, 88.9059, 91.224, 92.9064, 95.94, 98, 101.07, 102.9055,
                    106.42, 107.8682, 112.411, 114.818, 118.71, 121.76, 126.9045, 127.6, 131.293, 132.9055,
                    137.327, 138.9055, 140.116, 140.9077, 144.24, 145, 150.36, 151.964, 157.25, 158.9253, 162.5,
                    164.9303, 167.259, 168.9342, 173.04, 174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217,
                    195.078, 196.9665, 200.59, 204.3833, 207.2, 208.9804, 209, 210, 222, 223, 226, 227, 231.0359,
                    232.0381, 237, 238.0289, 243, 244, 247, 247, 251, 252, 257, 258, 259, 262, 261, 262, 266, 264, 277,
                    268, None, 272, None, None, None, None, None, None, None
                ]
            }
        )

        self.periodicTable.set_index('Number', inplace=True)

    def GetSymbol(self, number):
        """
        docstring
        """
        return self.periodicTable.loc[number, 'Symbol']

    def GetName(self, number):
        """
        docstring
        """
        return self.periodicTable.loc[number, 'Name']

    def GetCommonIons(self, number):
        """
        docstring
        """
        return self.periodicTable.loc[number, 'Common Ions']

    def GetNumber(self, symbolOrName):
        """
        docstring
        """
        if len(symbolOrName) > 2:
            mask = (self.periodicTable.loc[:, 'Name'] == symbolOrName)
            return self.periodicTable.loc[mask].index[0]
        else:
            mask = (self.periodicTable.loc[:, 'Symbol'] == symbolOrName)
            return self.periodicTable.loc[mask].index[0]

    def GetMass(self, numberOrSymbol):
        """
        docstring
        """
        if type(numberOrSymbol) == int:
            return self.periodicTable.loc[numberOrSymbol, 'Standard Mass']
        else:
            number = self.GetNumber(numberOrSymbol)
            return self.periodicTable.loc[number, 'Standard Mass']
