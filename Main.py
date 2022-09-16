#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Components import Molecule, Optimization, Project, Docking
from Worker import Worker, MolWorker
from Models import MoleculesModel, ProjectsModel, JobsModel
from Calculations import Gaussian, MyVina, DockingPlotter
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from DB.molecules_db import MyZODB
from datetime import datetime
from loguru import logger
# import threading
# from PyQt5.QtCore import pyqtRemoveInputHook


class MainWindow(qtw.QMainWindow):

    closed = qtc.pyqtSignal()

    # @logger.catch
    def __init__(self):
        super().__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
        self.database = MyZODB()
        self.selected_mol = None
        self.selected_project = None
        self.set_models()
        self.resume_queue()
        calc_items = (
            'Gaussian', 'Vina', 'Reactividad local'
        )
        # to_calculate_Rg = [mol for mol in self.molecules_list if mol.Rg == 0]
        # if to_calculate_Rg:
            # self.Rg_calculation(to_calculate_Rg)
        self.uiCalculationsComboBox.addItems(calc_items)
        # pyqtRemoveInputHook()

    def set_models(self):
        # Jobs
        self.jobs_list = list(self.database.get_jobs_db)
        self.master_queue = [
            j for j in self.jobs_list if j.get_status in ('Programmed', 'Running')
        ]
        jobs_model = JobsModel(self.jobs_list)
        self.uiJobsTableView.setModel(jobs_model)
        self.uiJobsTableView.resizeColumnsToContents()
        self.uiJobsTableView.resizeRowsToContents()
        # Molecules
        self.molecules_list = list(self.database.get_molecules_db)
        mol_tree_model = MoleculesModel(self.molecules_list, self.jobs_list)
        self.uiMoleculesTree.setModel(mol_tree_model.create_model())
        self.uiMoleculesTree.header().setSectionResizeMode(qtw.QHeaderView.ResizeToContents)
        self.selected_mol = None
        to_calculate_Rg = [mol for mol in self.molecules_list if mol.Rg == 0]
        if to_calculate_Rg and not hasattr(self, 'mol_thread'):
            self.Rg_calculation(to_calculate_Rg)
        elif hasattr(self, 'mol_thread') and self.mol_thread.isFinished():
            self.Rg_calculation(to_calculate_Rg)
        # Projects
        self.projects_list = list(self.database.get_projects_db)
        print([p.__dict__ for p in self.projects_list])
        proj_tree_model = ProjectsModel(self.projects_list)
        self.uiProjectsTreeView.setModel(proj_tree_model.create_model())
        self.uiProjectsTreeView.header().setSectionResizeMode(qtw.QHeaderView.ResizeToContents)
        self.uiProjectsTreeView.expandAll()

    @logger.catch
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
                self.statusBar().showMessage(
                    f'Molécula seleccionada: {self.selected_mol.get_name} ({self.selected_mol.inchi_key})'
                )
        else:
            self.selected_mol = None

    def set_selected_project(self, value):
        if self.projects_list:
            has_parent = value.parent().data()
            while has_parent:
                value = value.parent()
                has_parent = value.parent().data()
            if value.data():
                self.selected_project = [
                    p for p in self.projects_list if p.name == value.data().strip()
                ][0]
        else:
            self.selected_project = None

    def select_mol(self):
        selected_mol = self.selected_mol is not None
        if self.selected_project is not None:
            already_in = self.selected_mol in self.selected_project.molecules
            if selected_mol and not already_in:
                self.selected_project.add_molecule(self.selected_mol)
                self.database.commit()
                self.set_models()

    def unselect_mol(self):
        if self.selected_project is not None:
            self.selected_project.remove_molecule()
            self.database.commit()
            self.set_models()

    @qtc.pyqtSlot(float, str, object)
    def Rg_workflow(self, Rg, inchi_key, conf):
        molecule = self.database.get('molecules', inchi_key)
        molecule.set_Rg(Rg)
        molecule.set_conformer(conf)
        self.database.commit()
        message = f'Rg value for {molecule.get_name} -> {Rg:.2f}'
        self.statusBar().showMessage(message)
        self.set_models()

    @logger.catch
    def Rg_calculation(self, mol_list):
        self.mol_worker = MolWorker(tuple(mol_list))
        self.mol_thread = qtc.QThread()
        self.mol_worker.moveToThread(self.mol_thread)
        self.mol_thread.started.connect(self.mol_worker.start)
        self.mol_worker.finished.connect(self.mol_thread.quit)
        self.mol_worker.workflow.connect(self.Rg_workflow)
        self.mol_thread.start()

    def add_molecule(self):
        mol_path, _ = qtw.QFileDialog.getOpenFileNames(
            self, 'Selecciona la molécula a cargar',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Archivos de moléculas (*.mol *.mol2 *.pdb *.txt *.smi)'
        )
        if mol_path:
            # mol_list = []
            for m in mol_path:
                file_format = m.split('.')[-1]
                if file_format in ('txt', 'smi'):
                    with open(m, 'r') as file:
                        smi_lines = file.readlines()
                        for smi in smi_lines:
                            molecule = Molecule(smiles=smi)
                            if molecule.mol:
                                self.store_molecule(molecule)
                            else:
                                del molecule # a message can be displayed
                            # mol_list.append(molecule)
                        # self.Rg_calculation(mol_list)
                        return
                molecule = Molecule(path=m, file_format=file_format)
                if molecule.mol:
                    self.store_molecule(molecule)
                else:
                    del molecule
                # mol_list.append(molecule)
                # self.Rg_calculation(mol_list)

    def store_molecule(self, molecule):
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
            for project in self.projects_list:
                if self.selected_mol in project.molecules:
                    qtw.QMessageBox.critical(
                        self, 'Molécula en proyecto',
                        f'La molécula {self.selected_mol.get_name} hace parte del proyecto {project.name}. No será borrada.'
                    )
                    return False
            self.database.remove(to_remove)
            self.set_models()

    def create_project(self):
        name = self.uiProjectNameLine.text()
        exists = self.database.check('projects', name)
        if exists:
            qtw.QMessageBox.critical(
                self, 'Proyecto existente', f'El proyecto con nombre: {name} ya existe en la base de datos. No será creado.'
            )
        else:
            if name:
                project = Project(name)
                self.database.set('projects', project.name, project)
                self.uiProjectNameLine.clear()
                self.set_models()

    # revisar
    def add_group_calculation(self):
        calc = self.uiCalculationsComboBox.currentText()
        if calc == 'Gaussian':
            qtw.QMessageBox.critical(
                self, 'Cálculo grupal', 'Los valores (keywords) utilizados para la primera molécula serán aplicados el resto.'
            )
            self.gaussian_setup(group=True)

    def remove_group_calc(self):
        if self.selected_project is not None:
            self.selected_project.remove_calculation()
            self.database.commit()
            self.set_models()

    @qtc.pyqtSlot(dict, str)
    def queue_manager(self, calculation, project_name):
        calculation['id'] = self.database.get_job_id
        if calculation.get('type') == 'Optimization':
            job = Optimization(**calculation)
        elif calculation.get('type') == 'Docking':
            job = Docking(**calculation)
            job.create_config()
        mol = self.database.get('molecules', calculation['molecule_id'])
        project = self.database.get('projects', project_name)
        mol.add_calculation(job.id)
        project.add_job_id(job.id)
        self.database.set('jobs', job.id, job)
        self.set_models()

    @logger.catch
    def start_project(self, log):
        # print(log)
        # escogre si ejecutar en serie o en paralelo
        project = self.selected_project
        has_calculations = project.calculations
        # already_programmed = project.status == 'Programmed'
        if project is not None and has_calculations:
            for calc in self.selected_project.calculations:
                if calc['type'] in ('Optimization', 'Energy', 'Frequency'):
                    # esto sería para hacer en paralelo; invertir el if anterior y
                    # el siguiente for para hacer en serie (según moléculas)
                    for molecule in project.molecules:
                        self.gauss_controller = Gaussian(
                            molecule,
                            keywords=calc['keywords'],
                            job_type=calc['type']
                        )
                        self.gauss_controller.submitted.connect(self.queue_manager)
                        self.gauss_controller.queue_calculation()
                elif calc['type'] == 'Docking': # buscar otra forma con menos args
                    for molecule in project.molecules:
                        self.vina_controller = MyVina(
                            molecule,
                            config=calc['config'],
                            project=project,
                            times=calc['times'],
                            auto_box_size=calc['auto_box_size'],
                            nat_lig_path=calc['nat_lig_path'],
                            redocking=calc['redocking']
                        )
                        self.vina_controller.submitted.connect(self.queue_manager)
                        self.vina_controller.queue_calculation()
            # project.set_status('Programmed')
            self.database.commit()

    def gaussian_setup(self, group=False):
        """gaussian_setup.
        Opens the Gaussian controller to get the parameters for a calculation
        with g11

        Parameters
        ----------
        group :
            bool, if True, the parameters for the first molecule will be used
            for all the others.
        """
        if not group:
            self.gauss_controller = Gaussian(self.selected_mol)
            self.gauss_controller.submitted.connect(self.queue_manager)
        else:
            self.gauss_controller = Gaussian(
                self.selected_project.molecules[0], group=True
            )
            self.gauss_controller.submitted.connect(self.set_group_calc)

    def vina_setup(self, group=False):
        """vina_setup.
        Opens the Vina controller to get the parameters for a calculation
        with vina

        Parameters
        ----------
        group :
            bool, if True, the parameters for the first molecule will be used
            for all the others.
        """
        if not group:
            self.vina_controller = MyVina(self.selected_mol)
            self.vina_controller.submitted.connect(self.queue_manager)
        else:
            Rg_pending = [mol for mol in self.selected_project.molecules if mol.Rg == 0]
            if Rg_pending:
                message = 'Some molecules in this project are pending for the calculation of Rg values. wait for a moment until they are done, then try again'
                qtw.QMessageBox.critical(self, 'Rg calculation pending', message)
                return
            self.vina_controller = MyVina(
                self.selected_project.molecules[0],
                project=self.selected_project
            )
            self.vina_controller.submitted.connect(self.set_group_calc)

    @qtc.pyqtSlot(dict)
    def set_group_calc(self, calculation):
        self.selected_project.add_calculation(calculation)
        self.database.commit()
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
        # molecule = self.database.get('molecules', job.molecule_id)
        job.set_status(status)
        self.database.commit()
        message = f'{job.type}: {job.molecule} -> {status}'
        self.statusBar().showMessage(message)
        self.set_models()

    def start_master_queue(self):
        self.uiStartQueueButton.setEnabled(False)
        self.uiPauseQueueButton.setEnabled(True)
        self.worker = Worker(tuple(self.master_queue))
        self.thread = qtc.QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.start_queue)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(
            lambda: self.uiStartQueueButton.setEnabled(True)
        )
        self.worker.workflow.connect(self.workflow)
        self.uiPauseQueueButton.clicked.connect(self.worker.pause)
        self.thread.start()

    def plot_docking_results(self, project):
        if self.selected_project.calculations:
            for index, calculation in enumerate(self.selected_project.calculations):
                if calculation.get('type') == 'Docking':
                    self.docking_plotter = DockingPlotter(
                        project, [j for j in self.jobs_list if j.id in project.job_ids]
                    )

    def projects_right_click(self, position):
        if self.selected_project:
            menu = qtw.QMenu()
            calculation = qtw.QMenu('Cálculo grupal')
            results = qtw.QMenu('Visualizar resultados')
            menu.addMenu(calculation)
            menu.addMenu(results)
            gauss = calculation.addAction('Gaussian')
            vina = calculation.addAction('Vina')
            local_react = calculation.addAction('Reactividad local')
            docking_results = results.addAction('Docking')
            action = menu.exec_(self.uiProjectsTreeView.mapToGlobal(position))
            if action == gauss and self.selected_project is not None:
                self.gaussian_setup(group=True)
            elif action == vina and self.selected_project is not None:
                self.vina_setup(group=True)
            elif action == local_react:
                pass
            elif action == docking_results and self.selected_project is not None:
                self.plot_docking_results(self.selected_project)

    def molecules_right_click(self, position):
        if self.selected_mol:
            menu = qtw.QMenu()
            gauss = qtw.QMenu('Gaussian')
            vina = qtw.QMenu('Vina')
            menu.addMenu(gauss)
            menu.addMenu(vina)
            calculate = gauss.addAction('Calcular')
            docking = vina.addAction('Hacer Docking')
            action = menu.exec_(self.uiMoleculesTree.mapToGlobal(position))
            if action == calculate:
                self.gaussian_setup()
            elif action == docking:
                self.vina_setup()

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

# TODO
# Docking para proyecto,
# Arreglar lo de los hidrógenos en la representación de grilla en los proyectos
# Gráfica de energías vs moléculas en el proyecto
# poder descartar un cálculo que fue programado, borrando los archivos que haya
# generado

