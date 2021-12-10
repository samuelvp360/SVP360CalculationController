# -*- coding: utf-8 -*-

import sys
from Components import Molecule
from Models import MoleculesModel, SelectionModel
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from DB.molecules_db import MyZODB
# import pdb
# from PyQt5.QtCore import pyqtRemoveInputHook


class MainWindow(qtw.QMainWindow):

    closed = qtc.pyqtSignal()

    def __init__(self):
        super().__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
        self.database = MyZODB()
        self.selection_list = []
        self.selected_mol = None
        self.set_models()
        # pyqtRemoveInputHook()

    def set_models(self):
        self.molecules_list = list(self.database.get_db)
        tree_model = MoleculesModel(self.molecules_list)
        self.uiMoleculesTree.setModel(tree_model.create_model())
        self.uiMoleculesTree.header().setSectionResizeMode(qtw.QHeaderView.ResizeToContents)
        list_model = SelectionModel(self.selection_list)
        self.uiSelectionList.setModel(list_model)
        self.selected_mol = None

    def set_selected_mol(self, value):
        if self.molecules_list:
            has_parent = value.parent().data()
            while has_parent:
                value = value.parent()
                has_parent = value.parent().data()
            if value.data():
                self.selected_mol = [
                    m for m in self.molecules_list if m.get_name == value.data().strip()
                ][0]
        else:
            self.selected_mol = None

    def select(self):
        selected = self.selected_mol is not None
        already_in = self.selected_mol in self.selection_list
        if selected and not already_in:
            self.selection_list.append(self.selected_mol)
            self.set_models()

    def unselect(self):
        if self.selected_mol is not None and self.selected_mol in self.selection_list:
            index = self.selection_list.index(self.selected_mol)
            self.selection_list.pop(index)
            self.set_models()

    def add_molecule(self):
        mol_path, _ = qtw.QFileDialog.getOpenFileNames(
            self, 'Selecciona la molécula a cargar',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Archivos de moléculas (*.mol *.mol2 *pdb)'
        )
        if mol_path:
            for m in mol_path:
                file_format = m.split('.')[-1]
                molecule = Molecule(m, file_format)
                if molecule.mol:
                    exists = self.database.check(molecule.inchi_key)
                    if not exists:
                        self.database.set(molecule)
                        self.set_models()
                    else:
                        qtw.QMessageBox.critical(
                            self, 'Molécula existente', f'La molécula con Inchy key: {molecule.inchi_key} ya existe en la base de datos. No será agregada.'
                        )

    def remove_molecule(self):
        if self.selected_mol:
            to_remove = self.selected_mol.inchi_key
            self.database.remove(to_remove)
            self.set_models()

    def closeEvent(self, event):
        self.database.close()
        self.closed.emit()
        event.accept()

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

