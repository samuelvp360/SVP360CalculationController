# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from Views import resources
from rdkit.Chem import Draw
from loguru import logger


class PandasModel(qtc.QAbstractTableModel):
    '''Generic model to populate a pandas containing model'''

    def __init__(self, data, chunk=1):
        super().__init__()
        if chunk == 1:
            self.data = data if data.shape[0] <= 101 else data[:100]
        else:
            self.data = data[(chunk - 1) * 101:] \
                    if data.shape[0] <= (chunk * 100) + 1 \
                    else data[(chunk - 1) * 101:chunk * 101]

    def data(self, index, role):
        value = self.data.iloc[index.row(), index.column()]
        if role == qtc.Qt.DisplayRole:
            return str(value)
        elif role == qtc.Qt.TextAlignmentRole:
            if index.column() == 0:
                return qtc.Qt.AlignLeft | qtc.Qt.AlignVCenter
            else:
                return qtc.Qt.AlignRight | qtc.Qt.AlignVCenter

    def rowCount(self, index):
        return self.data.shape[0]

    def columnCount(self, index):
        return self.data.shape[1]

    def headerData(self, section, orientation, role):
        headers_v = self.data.index.values.tolist()
        headers_h = self.data.columns.values.tolist()
        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Vertical:
                return str(headers_v[section])
            if orientation == qtc.Qt.Horizontal:
                return str(headers_h[section])
        elif role == qtc.Qt.TextAlignmentRole:
            if orientation == qtc.Qt.Horizontal:
                return qtc.Qt.AlignCenter | qtc.Qt.AlignVCenter
            if orientation == qtc.Qt.Vertical:
                return qtc.Qt.AlignLeft | qtc.Qt.AlignVCenter

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


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
            self.setData(image, qtc.Qt.DecorationRole)


class MoleculesModel(qtg.QStandardItemModel):

    def __init__(self, molecules_list, calculations_list):
        super().__init__()
        self.populated_tree = self.populate_tree(
            molecules_list, calculations_list
        )

    def populate_tree(self, molecules_list, calculations_list):
        if molecules_list:
            std_item_list = []
            for mol in molecules_list:
                m_std_item_2 = StandardItem(img=mol.mol_pic)
                m_std_item_1 = StandardItem(
                    f'{mol.get_name:30}', font_size=14,
                    set_bold=True, color='#236e96'
                )
                summary = StandardItem(
                    'Propiedades', font_size=12, set_bold=True
                )
                formula_1 = StandardItem('Formula', font_size=10)
                formula_2 = StandardItem(
                    self.translate_formula(mol.formula), font_size=10
                )
                FW_1 = StandardItem('FW', font_size=10)
                FW_2 = StandardItem(f'{mol.MW:.2f} uma', font_size=10)
                inchi_1 = StandardItem('Inchi Key', font_size=10)
                inchi_2 = StandardItem(f'{mol.inchi_key}', font_size=10)
                smiles_1 = StandardItem('Smiles', font_size=10)
                smiles_2 = StandardItem(f'{mol.smiles}', font_size=10)
                Rg_1 = StandardItem('Rg', font_size=10)
                Rg_2 = StandardItem(f'{mol.descriptors.loc[0, "Rg"]:.2f}', font_size=10)
                std_item_list.append((m_std_item_1, m_std_item_2))
                summary.appendRow([formula_1, formula_2])
                summary.appendRow([FW_1, FW_2])
                summary.appendRow([inchi_1, inchi_2])
                summary.appendRow([smiles_1, smiles_2])
                summary.appendRow([Rg_1, Rg_2])
                calculations = StandardItem(
                    'Cálculos disponibles', font_size=12, set_bold=True
                )
                id_calc = StandardItem('ID', font_size=10, set_bold=True)
                calc_type = StandardItem(
                    'Tipo de cálculo', font_size=10, set_bold=True
                )
                calculations.appendRow([id_calc, calc_type])
                for c in mol.calculations:
                    calc_type = [
                        calc.type for calc in calculations_list if calc.id == c
                    ]
                    if not calc_type:
                        continue
                    else:
                        calc_type = calc_type[0]
                    calc_type = StandardItem(calc_type, font_size=10)
                    calc_id = StandardItem(str(c), font_size=10)
                    calculations.appendRow([calc_id, calc_type])
                m_std_item_1.appendRows([summary, calculations])
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

    def translate_formula(self, formula):
        trans = str.maketrans({
            '1': '\u2081', '2': '\u2082', '3': '\u2083', '4': '\u2084',
            '5': '\u2085', '6': '\u2086', '7': '\u2087', '8': '\u2088',
            '9': '\u2089', '0': '\u2080', '+': '\u207A', '-': '\u207B'
        })
        return formula.translate(trans)


class ProjectsModel(qtg.QStandardItemModel):

    def __init__(self, projects_list):
        super().__init__()
        self.populated_tree = self.populate_tree(projects_list)

    def populate_tree(self, projects_list):
        if projects_list:
            std_item_list = []
            for project in projects_list:
                p_std_item = StandardItem(
                    f'{project.name}', font_size=12,
                    set_bold=True, color='#236e96'
                )
                molecules = StandardItem(
                    'Molecules', font_size=10, set_bold=True
                )
                if project.molecules:
                    mol_std_item = StandardItem(img=project.grid_img)
                else:
                    mol_std_item = StandardItem(img='')
                molecules.appendRow(mol_std_item)
                calculations = StandardItem(
                    'Calculations', font_size=10, set_bold=True
                )
                for index, calc in enumerate(project.calculations):
                    calc_type = f'\nType:\t\t{calc["type"]}'
                    keywords = f'\nKeywords:\t\t{calc.get("keywords")}'
                    multipl = f'\nMultipl.:\t\t{calc["charge_mult"].split("/")[1]}'
                    status = f'\nStatus:\t\t{calc.get("status")}'
                    string = f'{index + 1}.{calc_type}{keywords}{multipl}{status}'
                    calc_std_item = StandardItem(string, font_size=8)
                    calculations.appendRow(calc_std_item)
                p_std_item.appendRows([molecules, calculations])
                std_item_list.append(p_std_item)
            return std_item_list
        else:
            return []

    def create_model(self):
        tree_model = qtg.QStandardItemModel()
        tree_model.setColumnCount(2)
        rood_node = tree_model.invisibleRootItem()
        if self.populated_tree:
            for item_1 in self.populated_tree:
                rood_node.appendRow([item_1])
        return tree_model


class JobsModel(qtc.QAbstractTableModel):

    def __init__(self, jobs_list):
        super().__init__()
        self.jobs_list = jobs_list[::-1]
        self.headers = (
            'ID', 'Molécula', 'Tipo de trabajo', 'Carga/Mult.',
            'Keywords', 'Output', 'Estatus', 'Programado', 'Comienzo',
            'Final', 'Transcurrido'
        )

    def data(self, index, role):
        value = self.jobs_list[index.row()]
        time_format = '%H:%M:%S; %d.%m.%Y'
        if role == qtc.Qt.DisplayRole:
            if index.column() == 0:
                return str(value.id)
            if index.column() == 1:
                return f"{value.molecule} ({value.molecule_id})"
            elif index.column() == 2:
                return str(value.type)
            elif index.column() == 3:
                return str(value.charge_mult)
            elif index.column() == 4:
                return str(value.keywords)
            elif index.column() == 5:
                if isinstance(value.output_file, (list, tuple)):
                    out = (i.split('/')[2] for i in value.output_file)
                    return '\n'.join(out)
                else:
                    return str(value.output_file.split('/')[2])
            elif index.column() == 6:
                return str(value.get_status)
            elif index.column() == 7:
                time = getattr(value, 'programmed', '')
                return time.strftime(time_format) if time else ''
            elif index.column() == 8:
                time = getattr(value, 'started', '')
                return time.strftime(time_format) if time else ''
            elif index.column() == 9:
                time = getattr(value, 'finished', '')
                return time.strftime(time_format) if time else ''
            elif index.column() == 10:
                time = getattr(value, 'elapsed', '')
                return str(time) if time else ''

        if role == qtc.Qt.ForegroundRole:
            status = index.column() == 6
            if status and value.get_status == 'Programmed':
                return qtg.QColor('orange')
            elif status and value.get_status == 'Running':
                return qtg.QColor('blue')
            elif status and value.get_status == 'Finished':
                return qtg.QColor('green')
            elif status and value.get_status == 'Failed':
                return qtg.QColor('red')

    def rowCount(self, index):
        return len(self.jobs_list)

    def columnCount(self, index):
        return len(self.headers)

    def headerData(self, section, orientation, role):
        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Horizontal:
                return self.headers[section]


# class AvailableCalcModel(qtc.QAbstractTableModel):
    # '''Model to populate the Available Calculations Table'''

    # def __init__(self, molecule):
        # super(AvailableCalcModel, self).__init__()
        # self._calculations = molecule.GetCalculations
        # self._headers = (
            # 'JOB TYPE', 'DESCRIPTION', 'CHARGE/MULT', 'RESULTS',
            # 'STATUS', 'DATE OF RUN', 'ELAPSED TIME'
        # )

    # def data(self, index, role):
        # if role == qtc.Qt.DisplayRole:
            # if index.column() != 3:
                # return self._calculations[index.row()][self._headers[index.column()]]
            # else:
                # return 'Optimized coordinates'

    # def rowCount(self, index):
        # return len(self._calculations)

    # def columnCount(self, index):
        # return len(self._headers)

    # def headerData(self, section, orientation, role):
        # if role == qtc.Qt.DisplayRole:
            # if orientation == qtc.Qt.Horizontal:
                # return self._headers[section]
            # if orientation == qtc.Qt.Vertical:
                # return section

    # def flags(self, index):
        # return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


# class StatusModel(qtc.QAbstractTableModel):
    # '''Model to populate the list of calculations to do'''

    # def __init__(self, calcToDoList):

        # super(StatusModel, self).__init__()
        # self.calcToDoList = calcToDoList
        # self._headers = (
            # 'Molecule', 'Job Type', 'Keywords', 'Status'
        # )

    # def data(self, index, role):

        # if role == qtc.Qt.DisplayRole:
            # if index.column() == 0:
                # return self.calcToDoList[index.row()][index.column()].GetName
            # return self.calcToDoList[index.row()][index.column()]

        # if role == qtc.Qt.ForegroundRole:
            # if self.calcToDoList[index.row()][3] == 'Pending':
                # return qtg.QColor('#FD971F')
            # elif self.calcToDoList[index.row()][3] == 'Running':
                # return qtg.QColor('#66D9EF')
            # elif self.calcToDoList[index.row()][3] == 'Finished':
                # return qtg.QColor('#A6E22E')
            # elif self.calcToDoList[index.row()][3] == 'Failed':
                # return qtg.QColor('#F92672')

    # def rowCount(self, index):
        # return len(self.calcToDoList)

    # def columnCount(self, index):
        # return 4

    # def headerData(self, section, orientation, role):

        # if role == qtc.Qt.DisplayRole:
            # if orientation == qtc.Qt.Horizontal:
                # return self._headers[section]
            # if orientation == qtc.Qt.Vertical:
                # return section

    # def flags(self, index):
        # return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


# class ResultsModel(qtc.QAbstractTableModel):
    # '''Model to populate the Results Table'''

    # def __init__(self, calculation):
        # super(ResultsModel, self).__init__()
        # self._calculation[0] = calculation
        # self._headers = ('ATOM', 'X', 'Y', 'Z')

    # def data(self, index, role):
        # if index.isValid():
            # if role == qtc.Qt.DisplayRole:
                # return str(self._calculation.iloc[index.row(), index.column()])
        # return None

    # def rowCount(self, index):
        # return self._calculation.shape[0]

    # def columnCount(self, index):
        # return self._calculation.shape[1]

    # def headerData(self, section, orientation, role):
        # if role == qtc.Qt.DisplayRole:
            # if orientation == qtc.Qt.Horizontal:
                # return self._calculation.columns[section]
            # if orientation == qtc.Qt.Vertical:
                # return section + 1
        # return None

    # def flags(self, index):
        # return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


class IRDataModel(qtc.QAbstractTableModel):
    '''Model to populate the IR Bands Table'''

    def __init__(self, bands_df, y_axis):
        super().__init__()
        self.bands = bands_df
        self.y_axis = y_axis

    def data(self, index, role):
        if role == qtc.Qt.DisplayRole:
            if index.column() == 0:
                return str(round(self.bands.iloc[index.row(), index.column()], 2))
            elif index.column() == 1:
                return str(round(self.bands.loc[index.row(), self.y_axis], 6))
            elif index.column() == 2:
                return str(round(self.bands.loc[index.row(), 'HWHM'], 2))
            elif index.column() == 3:
                if self.bands.get('Beta') is not None:
                    return str(round(self.bands.loc[index.row(), 'Beta'], 4))
                else:
                    return ''
        if role == qtc.Qt.TextAlignmentRole:
            return qtc.Qt.AlignCenter

    def rowCount(self, index):
        return self.bands.shape[0]

    def columnCount(self, index):
        return 4

    def headerData(self, section, orientation, role):
        if role == qtc.Qt.DisplayRole:
            if orientation == qtc.Qt.Horizontal:
                if section == 0:
                    return 'Wavenumber'
                elif section == 1:
                    return self.y_axis
                elif section == 2:
                    return 'HWHM'
                else:
                    return '\u03b2'
            if orientation == qtc.Qt.Vertical:
                return section + 1

    def flags(self, index):
        return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable


class SimilarityModel(PandasModel):

    def data(self, index, role):
        value = self.data.iloc[index.row(), index.column()]
        if role == qtc.Qt.DisplayRole:
            return str(value)
        elif role == qtc.Qt.TextAlignmentRole:
            return qtc.Qt.AlignCenter | qtc.Qt.AlignVCenter
        # Background
        if role == qtc.Qt.BackgroundRole:
            bg_color = qtg.QColor('red')
            value = value if value > 0 else 0
            bg_color.setAlphaF(value)
            return bg_color


class DisplayMFPModel(PandasModel):

    def data(self, index, role):
        value = self.data.iloc[index.row(), index.column()]
        if role == qtc.Qt.DisplayRole:
            return str(round(value, 3))
        elif role == qtc.Qt.TextAlignmentRole:
            return qtc.Qt.AlignCenter | qtc.Qt.AlignVCenter
        # Foreground
        if role == qtc.Qt.ForegroundRole:
            return qtg.QColor('white')
        # Background
        if role == qtc.Qt.BackgroundRole:
            if value:
                return qtg.QColor('green')
            else:
                return qtg.QColor('red')



# class AvailableSpectraModel(qtc.QAbstractTableModel):
    # '''Model to populate the IR available data'''

    # def __init__(self, availableIRData):
        # super(AvailableSpectraModel, self).__init__()
        # self._availableIRData = availableIRData

    # def data(self, index, role):
        # keys = [i for i in self._availableIRData.keys()]
        # data = self._availableIRData.get(keys[index.row()])
        # if role == qtc.Qt.DisplayRole:
            # if index.column() == 0:
                # return data['TYPE']
            # elif index.column() == 1:
                # return data['SOLVENT']
            # elif index.column() == 2:
                # return data['DATE']

    # def rowCount(self, index):
        # return len(self._availableIRData.items())

    # def columnCount(self, index):
        # return 3

    # def headerData(self, section, orientation, role):
        # if role == qtc.Qt.DisplayRole:
            # if orientation == qtc.Qt.Horizontal:
                # if section == 0:
                    # return 'Type'
                # elif section == 1:
                    # return 'Solvent'
                # elif section == 2:
                    # return 'Date'

            # if orientation == qtc.Qt.Vertical:
                # return section + 1

        # if role == qtc.Qt.ForegroundRole:
            # return qtg.QColor('#66D9EF')

    # def flags(self, index):
        # return qtc.Qt.ItemIsEnabled | qtc.Qt.ItemIsSelectable
