# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from Views import resources


class StandardItem(qtg.QStandardItem):
    def __init__(
        self, txt='', img='', font_size=12, set_bold=False,
        color='black'
    ):
        super().__init__()

        fnt = qtg.QFont('Open Sans', font_size)
        fnt.setBold(set_bold)
        self.setEditable(False)
        self.setFont(fnt)
        self.setText(txt)
        self.setForeground(qtg.QColor(color))
        if img:
            image = qtg.QPixmap(img)
            image = image.scaledToWidth(150)
            self.setData(image, qtc.Qt.DecorationRole)


class MoleculesModel(qtg.QStandardItemModel):
    def __init__(self, molecules_list):
        super().__init__()
        self.populated_tree = self.populate_tree(molecules_list)

    def populate_tree(self, molecules_list):
        if molecules_list:
            std_item_list = []
            for mol in molecules_list:
                m_std_item_2 = StandardItem('', img=mol.mol_pic)
                m_std_item_1 = StandardItem(
                    f'{mol.name:120}', font_size=14, set_bold=True,
                    color='#236e96'
                )
                summary = StandardItem(
                    'Propiedades', font_size=12, set_bold=True
                )
                formula_1 = StandardItem('Formula', font_size=10)
                formula_2 = StandardItem(f'{mol.formula}', font_size=10)
                FW_1 = StandardItem('FW', font_size=10)
                FW_2 = StandardItem(f'{mol.MW:.2f} uma', font_size=10)
                inchi_1 = StandardItem('Inchi Key', font_size=10)
                inchi_2 = StandardItem(f'{mol.inchi_key}', font_size=10)
                smiles_1 = StandardItem('Smiles', font_size=10)
                smiles_2 = StandardItem(f'{mol.smiles}', font_size=10)
                std_item_list.append((m_std_item_1, m_std_item_2))
                summary.appendRow([formula_1, formula_2])
                summary.appendRow([FW_1, FW_2])
                summary.appendRow([inchi_1, inchi_2])
                summary.appendRow([smiles_1, smiles_2])
                m_std_item_1.appendRows([summary])
            return std_item_list
        else:
            return []

    def create_model(self):
        tree_model = qtg.QStandardItemModel()
        tree_model.setColumnCount(2)
        rood_node = tree_model.invisibleRootItem()
        if self.populated_tree:
            for item_1, item_2 in self.populated_tree:
                rood_node.appendRow([item_1, item_2])
        return tree_model


class AvailableCalcModel(qtc.QAbstractTableModel):
    '''Model to populate the Available Calculations Table'''

    def __init__(self, molecule):
        super(AvailableCalcModel, self).__init__()
        self._calculations = molecule.GetCalculations
        self._headers = (
            'JOB TYPE', 'DESCRIPTION', 'CHARGE/MULT', 'RESULTS',
            'STATUS', 'DATE OF RUN', 'ELAPSED TIME'
        )

    def data(self, index, role):
        if role == qtc.Qt.DisplayRole:
            if index.column() != 3:
                return self._calculations[index.row()][self._headers[index.column()]]
            else:
                return 'Optimized coordinates'

    def rowCount(self, index):
        return len(self._calculations)

    def columnCount(self, index):
        return len(self._headers)

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
                return qtg.QColor('#FD971F')
            elif self.calcToDoList[index.row()][3] == 'Running':
                return qtg.QColor('#66D9EF')
            elif self.calcToDoList[index.row()][3] == 'Finished':
                return qtg.QColor('#A6E22E')
            elif self.calcToDoList[index.row()][3] == 'Failed':
                return qtg.QColor('#F92672')

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


class ResultsModel(qtc.QAbstractTableModel):
    '''Model to populate the Results Table'''

    def __init__(self, calculation):
        super(ResultsModel, self).__init__()
        self._calculation[0] = calculation
        self._headers = ('ATOM', 'X', 'Y', 'Z')

    def data(self, index, role):
        if index.isValid():
            if role == qtc.Qt.DisplayRole:
                return str(self._calculation.iloc[index.row(), index.column()])
        return None

    def rowCount(self, index):
        return self._calculation.shape[0]

    def columnCount(self, index):
        return self._calculation.shape[1]

    def headerData(self, section, orientation, role):
        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Horizontal:
                return self._calculation.columns[section]
            if orientation == qtc.Qt.Vertical:
                return section + 1
        return None

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


class IRDataModel(qtc.QAbstractTableModel):
    '''Model to populate the IR Bands Table'''

    def __init__(self, bandsDict, yAxis):
        super(IRDataModel, self).__init__()
        self._bands = bandsDict
        self._axis = yAxis

    def data(self, index, role):
        if role == qtc.Qt.DisplayRole:
            if index.column() == 0:
                return str(round(self._bands.iloc[index.row(), index.column()], 2))
            elif index.column() == 1:
                return str(round(self._bands.loc[index.row(), self._axis], 6))
            elif index.column() == 2:
                return str(round(self._bands.loc[index.row(), 'HWHM'], 2))
            elif index.column() == 3:
                if self._bands.get('Beta') is not None:
                    return str(round(self._bands.loc[index.row(), 'Beta'], 4))
                else:
                    return ''
        if role == qtc.Qt.TextAlignmentRole:
            return qtc.Qt.AlignCenter

    def rowCount(self, index):
        return self._bands.shape[0]

    def columnCount(self, index):
        return 4

    def headerData(self, section, orientation, role):
        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Horizontal:
                if section == 0:
                    return 'Wavenumber'
                elif section == 1:
                    return self._axis
                elif section == 2:
                    return 'HWHM'
                else:
                    return 'Beta'
            if orientation == qtc.Qt.Vertical:
                return section + 1

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


class AvailableSpectraModel(qtc.QAbstractTableModel):
    '''Model to populate the IR available data'''

    def __init__(self, availableIRData):
        super(AvailableSpectraModel, self).__init__()
        self._availableIRData = availableIRData

    def data(self, index, role):
        keys = [i for i in self._availableIRData.keys()]
        data = self._availableIRData.get(keys[index.row()])
        if role == qtc.Qt.DisplayRole:
            if index.column() == 0:
                return data['TYPE']
            elif index.column() == 1:
                return data['SOLVENT']
            elif index.column() == 2:
                return data['DATE']

    def rowCount(self, index):
        return len(self._availableIRData.items())

    def columnCount(self, index):
        return 3

    def headerData(self, section, orientation, role):
        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Horizontal:
                if section == 0:
                    return 'Type'
                elif section == 1:
                    return 'Solvent'
                elif section == 2:
                    return 'Date'

            if orientation == qtc.Qt.Vertical:
                return section + 1

        if role == qtc.Qt.ForegroundRole:
            return qtg.QColor('#66D9EF')

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable
