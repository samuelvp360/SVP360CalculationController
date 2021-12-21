#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Components import Molecule, Optimization
from Worker import Worker
from Models import MoleculesModel, SelectionModel, JobsModel
from Calculations import Gaussian
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from DB.molecules_db import MyZODB
from datetime import datetime
# import pdb
# from PyQt5.QtCore import pyqtRemoveInputHook


class MainWindow(qtw.QMainWindow):

    closed = qtc.pyqtSignal()

    def __init__(self):
        super().__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
        self.database = MyZODB()
        self.mols_in_project_list = []
        self.selected_mol = None
        self.set_models()
        self.resume_queue()
        # pyqtRemoveInputHook()

    def set_models(self):
        self.jobs_list = list(self.database.get_jobs_db)
        self.master_queue = [
            j for j in self.jobs_list if j.get_status in ('Programmed', 'Running')
        ]
        self.molecules_list = list(self.database.get_molecules_db)
        tree_model = MoleculesModel(self.molecules_list, self.jobs_list)
        self.uiMoleculesTree.setModel(tree_model.create_model())
        self.uiMoleculesTree.header().setSectionResizeMode(qtw.QHeaderView.ResizeToContents)
        list_model = SelectionModel(self.mols_in_project_list)
        self.uiProjectsList.setModel(list_model)
        self.selected_mol = None
        jobs_model = JobsModel(self.jobs_list)
        self.uiJobsTableView.setModel(jobs_model)
        self.uiJobsTableView.resizeColumnsToContents()

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
        already_in = self.selected_mol in self.mols_in_project_list
        if selected and not already_in:
            self.mols_in_project_list.append(self.selected_mol)
            self.set_models()

    def unselect(self):
        if self.selected_mol is not None and self.selected_mol in self.mols_in_project_list:
            index = self.mols_in_project_list.index(self.selected_mol)
            self.mols_in_project_list.pop(index)
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
                    exists = self.database.check('molecules', molecule.inchi_key)
                    if not exists:
                        self.database.set(
                            'molecules', molecule.inchi_key, molecule
                        )
                        self.set_models()
                    else:
                        qtw.QMessageBox.critical(
                            self, 'Molécula existente', f'La molécula con Inchy key: {molecule.inchi_key} ya existe en la base de datos. No será agregada.'
                        )

    def remove_molecule(self):
        # hay que colorear de rojo los trabajos de moléculas que ya no estén en
        # la bd. Si están pendientes, no correr tampoco
        if self.selected_mol:
            to_remove = self.selected_mol.inchi_key
            self.database.remove(to_remove)
            self.set_models()

    def gaussian_setup(self):
        self.gauss_controller = Gaussian(self.selected_mol)
        self.gauss_controller.submitted.connect(self.queue_manager)
        self.gauss_controller.show()

    @qtc.pyqtSlot(dict)
    def queue_manager(self, calculation):
        calculation['id'] = self.database.get_job_id
        if calculation.get('type') == 'Optimization':
            job = Optimization(**calculation)
        self.selected_mol.add_calculation(job.id)
        self.database.set('jobs', job.id, job)
        self.set_models()

    def resume_queue(self):
        if self.master_queue:
            reply = qtw.QMessageBox.question(
                self, 'Trabajo pendiente en cola',
                'Aún hay trabajo en cola. ¿Desea reanudar?',
                qtw.QMessageBox.Yes | qtw.QMessageBox.No
            )
            if reply == qtw.QMessageBox.Yes:
                self.start_master_queue()

    def pause_master_queue(self):
        if hasattr(self, 'worker'):
            self.worker.pause()

    @qtc.pyqtSlot(int, str)
    def workflow(self, job_id, status):
        job = self.database.get('jobs', job_id)
        molecule = self.database.get('molecules', job.molecule_id)
        job.set_status(status)
        self.database.commit()
        self.set_models()

    def start_master_queue(self):
        self.uiStartQueueButton.setEnabled(False)
        self.uiPauseQueueButton.setEnabled(True)
        self.worker = Worker(self.master_queue)
        self.thread = qtc.QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.start_queue)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(lambda: self.uiStartQueueButton.setEnabled(True))
        self.worker.workflow.connect(self.workflow)
        self.uiPauseQueueButton.clicked.connect(self.worker.pause)
        self.thread.start()

    def right_click(self, position):
        if self.selected_mol:
            menu = qtw.QMenu()
            gauss = qtw.QMenu('Gaussian')
            vina = qtw.QMenu('Vina')
            menu.addMenu(gauss)
            menu.addMenu(vina)
            calculate = gauss.addAction('Calcular')
            prepare = vina.addAction('Preparar ligando')
            prepare_2 = vina.addAction('Preparar receptor')
            docking = vina.addAction('Hacer Docking')
            action = menu.exec_(self.uiMoleculesTree.mapToGlobal(position))
            if action == calculate:
                self.gaussian_setup()
            elif action == prepare:
                print('Preparando ligando')
            elif action == prepare_2:
                print('Preparando receptor')
            elif action == docking:
                print('Docking')

    def closeEvent(self, event):
        if hasattr(self, 'worker'):
            self.worker.pause()
        # permite que se termine el último trabajo y luego detiene la lista,
        # para evitar que los demás trabajos queden sin tiempo de ejecución.
        self.database.close()
        self.closed.emit()
        event.accept()

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())

