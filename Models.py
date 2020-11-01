#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from Views import resources


class MoleculesModel(qtc.QAbstractListModel):
    """Model to populate the list of the uploaded molecules"""

    def __init__(self, molecules):
        super(MoleculesModel, self).__init__()
        self.moleculesList = molecules

    def data(self, index, role):
        if role == qtc.Qt.DisplayRole:
            return self.moleculesList[index.row()].GetName

        if role == qtc.Qt.DecorationRole:
            if self.moleculesList[index.row()].stored:
                return qtg.QIcon(':/icons/okDB.png')

    def rowCount(self, index):
        return len(self.moleculesList)

    def headerData(self, section, orientation, role):
        if role == qtc.Qt.DisplayRole:
            if section == 0:
                return 'Molecules'

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


# class PropertiesModel(qtc.QAbstractListModel):
#     """Model to populate the Properties Widgets"""

#     def __init__(self, molecule):
#         super(PropertiesModel, self).__init__()
#         self.molecule = molecule

#     def data(self, index, role):
#         if role == qtc.Qt.DisplayRole or role == qtc.Qt.EditRole:
#             if index.row() == 0:
#                 return self.molecule.GetName
#             elif index.row() == 1:
#                 return self.molecule.GetForm
#             elif index.row() == 2:
#                 return self.molecule.GetMolarMass
#             elif index.row() == 3:
#                 return self.molecule.GetInchikey
#             elif index.row() == 4:
#                 return self.molecule.GetSmiles

#     def setData(self, index, value, role=qtc.Qt.EditRole):
#         if role == qtc.Qt.EditRole:
#             if index.row() == 0:
#                 self.molecule.SetName(value)
#                 self.dataChanged.emit(index, index)
#                 return True
#             else:
#                 return False

#     def rowCount(self, index):
#         return 6

#     def headerData(self, section, orientation, role):
#         if role == qtc.Qt.DisplayRole:
#             if section == 0:
#                 return 'Properties'

#     def flags(self, index):
#         return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable | qtc.Qt.ItemIsEditable


class StatusModel(qtc.QAbstractTableModel):
    """Model to populate the list of calculations to do"""

    def __init__(self, calcToDoList):

        super(StatusModel, self).__init__()
        self.calcToDoList = calcToDoList
        self._headers = (
            'Molecule', 'Job Type', 'Keywords', 'Status'
        )

    def data(self, index, role):

        if role == qtc.Qt.DisplayRole:
            if index.column() == 0:
                return self.calcToDoList[index.row()][index.column()].GetName
            return self.calcToDoList[index.row()][index.column()]

        if role == qtc.Qt.ForegroundRole:
            if self.calcToDoList[index.row()][3] == 'Pending':
                return qtg.QColor('orange')
            elif self.calcToDoList[index.row()][3] == 'Running':
                return qtg.QColor('blue')
            elif self.calcToDoList[index.row()][3] == 'Finished':
                return qtg.QColor('green')
            elif self.calcToDoList[index.row()][3] == 'Failed':
                return qtg.QColor('red')

    def rowCount(self, index):
        return len(self.calcToDoList)

    def columnCount(self, index):
        return 4

    def headerData(self, section, orientation, role):

        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Horizontal:
                return self._headers[section]
            if orientation == qtc.Qt.Vertical:
                return section

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable
