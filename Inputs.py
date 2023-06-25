# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import uic


class ActivityInput(qtw.QWidget):

    submitted = qtc.pyqtSignal(dict, str)

    def __init__(self, molecule, activities_db):
        super().__init__()
        uic.loadUi('Views/uiActivity.ui', self)
        self.molecule = molecule
        self.activities_db = activities_db
        self.id = ''
        self.unit = ''
        self.species = ''
        self.act_names = [act.name for act in activities_db]
        if self.act_names:
            completer = qtw.QCompleter(self.act_names)
            completer.setCaseSensitivity(qtc.Qt.CaseInsensitive)
            self.uiNameLine.setCompleter(completer)
        self.activity = {}
        self.show()

    def set_activity(self):
        act_name = self.uiNameLine.text().strip().lower()
        if act_name in self.act_names:
            act = [a for a in self.activities_db if a.name == act_name][0]
            self.id = act.id
            self.unit = act.unit
            self.species = act.species
            self.uiUnitLine.setText(self.unit)
            self.uiSpeciesLine.setText(self.species)
            self.uiUnitLine.setEnabled(False)
            self.uiSpeciesLine.setEnabled(False)
        else:
            self.name = act_name
            current_ids = [a.id for a in self.activities_db]
            self.id = max(current_ids) + 1 if current_ids else 0
        self.uiIDLabel.setText(str(self.id))
        self.uiSubmitButton.setEnabled(True)

    def submit_activity(self):
        self.activity = {
            'id': self.id,
            'name': self.uiNameLine.text().strip().lower(),
            'value': self.uiValueSpin.value(),
            'unit': self.uiUnitLine.text(),
            'species': self.uiSpeciesLine.text()
        }
        self.submitted.emit(self.activity, self.molecule.smiles)
        self.close()
