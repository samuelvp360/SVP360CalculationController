# -*- coding: utf-8 -*-

import sys
from Components import Molecule
from Models import MoleculesModel
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from DB.molecules_db import MyZODB


class MainWindow(qtw.QMainWindow):

    closed = qtc.pyqtSignal()

    def __init__(self):
        super().__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
        self.database = MyZODB()
        self.set_models()

    def set_models(self):
        self.molecules_list = list(self.database.get_db)
        tree_model = MoleculesModel(self.molecules_list)
        self.uiMoleculesTree.setModel(tree_model.create_model())
        self.uiMoleculesTree.header().setSectionResizeMode(qtw.QHeaderView.ResizeToContents)

    def set_selected_mol(self, value):
        if self.molecules_list:
            has_parent = value.parent().data()
            while has_parent:
                value = value.parent()
                has_parent = value.parent().data()
            if value.data():
                self.selected_mol = [
                    i for i, m in enumerate(self.molecules_list) if m.name == value.data().strip()
                ][0]
            else:
                self.selected_mol = None

    def add_molecule(self):
        mol_path, _ = qtw.QFileDialog.getOpenFileName(
            self, 'Selecciona la molécula a cargar',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Archivos de moléculas (*.mol)'
        )
        if mol_path:
            molecule = Molecule(mol_path, 'mol')
            self.database.set(molecule)
            self.set_models()

    def remove_molecule(self):
        if self.selected_mol:
            to_remove = self.molecules_list[self.selected_mol].inchi_key
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
