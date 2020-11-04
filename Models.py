#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from Views import resources


class MoleculesModel(qtc.QAbstractListModel):
    '''Model to populate the list of the uploaded molecules'''

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


class AvailableCalcModel(qtc.QAbstractTableModel):
    '''Model to populate the Available Calculations Table'''

    def __init__(self, molecule):
        super(AvailableCalcModel, self).__init__()
        self._calculations = molecule.GetCalculations
        self._headers = (
            'JOB TYPE', 'DESCRIPTION', 'RESULTS', 'STATUS', 'DATE OF RUN', 'ELAPSED TIME'
        )

    def data(self, index, role):
        if role == qtc.Qt.DisplayRole:
            return self._calculations[index.row()][self._headers[index.column()]]

    def rowCount(self, index):
        return len(self._calculations)

    def columnCount(self, index):
        return 6
        # return len(self._calculations[index.row()])

    def headerData(self, section, orientation, role):
        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Horizontal:
                return self._headers[section]
            if orientation == qtc.Qt.Vertical:
                return section

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


class StatusModel(qtc.QAbstractTableModel):
    '''Model to populate the list of calculations to do'''

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
