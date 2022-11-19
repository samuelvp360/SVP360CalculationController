#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shutil
import sys
from Components import Molecule, Optimization, Project, Docking
from Worker import Worker, MolWorker, GenericWorker
from Models import MoleculesModel, ProjectsModel, JobsModel
from Calculations import Gaussian, MyVina
from Plotters import DockingPlotter, RedockingPlotter
from Explorers import FPExplorer, SimilarityExplorer, \
        DescriptorsExplorer, ChEMBLExplorer
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import uic
from DB.molecules_db import MyZODB
from loguru import logger
import pandas as pd
# import numpy as np
# from PyQt5.QtCore import pyqtRemoveInputHook


class MainWindow(qtw.QMainWindow):

    closed = qtc.pyqtSignal()

    # @logger.catch
    def __init__(self):
        super().__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
        self.database = MyZODB()
        self.selected_mol = []
        self.selected_project = None
        self.queue_status = 'Paused'
        self.set_models()
        self.resume_queue()
        calc_items = (
            'Gaussian', 'Vina', 'Reactividad local'
        )
        self.uiCalculationsComboBox.addItems(calc_items)
        search_items = (
            'Name', 'InchiKey', 'Smiles', 'Molar mass', 'Formula'
        )
        self.uiSearchCombo.addItems(search_items)
        # pyqtRemoveInputHook()

    def set_models(self):
        # Jobs
        self.jobs_list = list(self.database.get_jobs_db)
        # for job in self.jobs_list:
            # self.database.remove('jobs', job.id)
        # self.database.commit()
        jobs_model = JobsModel(self.jobs_list)
        self.uiJobsTableView.setModel(jobs_model)
        self.uiJobsTableView.resizeColumnsToContents()
        self.uiJobsTableView.resizeRowsToContents()
        # Molecules
        if not self.uiSearchLine.text():
            self.molecules_list = list(self.database.get_molecules_db)
        mol_tree_model = MoleculesModel(self.molecules_list, self.jobs_list)
        self.uiMoleculesTree.setModel(mol_tree_model.create_model())
        self.uiMoleculesTree.header().setSectionResizeMode(qtw.QHeaderView.ResizeToContents)
        # Projects
        self.projects_list = list(self.database.get_projects_db)
        proj_tree_model = ProjectsModel(self.projects_list)
        self.uiProjectsTreeView.setModel(proj_tree_model.create_model())
        self.uiProjectsTreeView.header().setSectionResizeMode(qtw.QHeaderView.ResizeToContents)
        self.uiProjectsTreeView.expandAll()
        # Master queue
        self.master_queue = [
            j for j in self.jobs_list if j.get_status in ('Programmed', 'Running')
        ]
        if self.master_queue and self.queue_status != 'Running':
            self.uiStartQueueButton.setEnabled(True)
        else:
            self.uiStartQueueButton.setEnabled(False)
        if self.master_queue and not hasattr(self, 'queue_thread') \
           and self.queue_status != 'Paused':
            self.start_master_queue()
        elif hasattr(self, 'queue_thread') and self.queue_thread.isFinished() \
           and self.queue_status != 'Paused':
            self.start_master_queue()
        else:
            self.queue_status = 'Paused'
        self.uiTotalMoleculesLabel.setText(
            f'Total molecules: {len(self.molecules_list)}; Selected: {len(self.selected_mol)}'
        )

    def search_molecule(self):
        value = self.uiSearchLine.text().lower()
        criteriom = self.uiSearchCombo.currentText()
        self.molecules_list = list(self.database.get_molecules_db)
        if criteriom == 'Name':
            self.molecules_list = [
                m for m in self.molecules_list if value in m.get_name.lower()
            ]
        elif criteriom == 'InchiKey':
            self.molecules_list = [
                m for m in self.molecules_list if value in m.inchi_key.lower()
            ]
        elif criteriom == 'Smiles':
            self.molecules_list = [
                m for m in self.molecules_list if value in m.smiles.lower()
            ]
        elif criteriom == 'Molar mass':
            pass
        elif criteriom == 'Formula':
            self.molecules_list = [
                m for m in self.molecules_list if value in m.formula.lower()
            ]
        self.set_models()

    @logger.catch
    def set_selected_mol(self, value):
        indexes = self.uiMoleculesTree.selectedIndexes()
        if indexes:
            rows = [
                j.row() for i, j in enumerate(indexes) if i % 2 \
                and j.parent().row() == -1
            ]
            self.selected_mol = [
                m for i, m in enumerate(self.molecules_list) if i in rows
            ]
        if self.molecules_list and not rows:
            has_parent = value.parent().data()
            while has_parent:
                value = value.parent()
                has_parent = value.parent().data()
            if value.data():
                self.selected_mol = [
                    m for m in self.molecules_list if m.get_name == value.data().strip()
                ]
        if self.selected_mol:
            mol_names = [m.get_name for m in self.selected_mol]
            self.statusBar().showMessage(
                f'Selected Molecules: {", ".join(mol_names)}'
            )
        self.uiTotalMoleculesLabel.setText(
            f'Total molecules: {len(self.molecules_list)}; Selected: {len(self.selected_mol)}'
        )

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

    def include_mol(self):
        if self.selected_project is not None and self.selected_mol:
            already_in = [
                m for m in self.selected_mol if m \
                in self.selected_project.molecules
            ]
            for m in self.selected_mol:
                if m not in already_in:
                    self.selected_project.add_molecule(m)
                    self.selected_project.add_fps(m.get_fps_dict)
                    self.selected_project.add_descriptors(
                        m.get_descriptors, m.inchi_key
                    )
                    self.database.commit()
                    self.set_models()

    def drop_mol(self):
        if self.selected_project is not None:
            self.selected_project.pop_molecule()
            self.database.commit()
            self.set_models()

    # @logger.catch
    # def descriptors_calc(self, mol_list):
        # self.mol_worker = MolWorker(tuple(mol_list))
        # self.mol_thread = qtc.QThread()
        # self.mol_worker.moveToThread(self.mol_thread)
        # self.mol_thread.started.connect(self.mol_worker.start)
        # self.mol_worker.finished.connect(self.mol_thread.quit)
        # self.mol_worker.workflow.connect(self.descriptors_workflow)
        # self.mol_thread.start()

    def open_file(self):
        mol_path, _ = qtw.QFileDialog.getOpenFileNames(
            self, 'Select the molecule',
            filter='*.mol;;*.mol2;;*.pdb;;*.txt;;*.smi;;*.log;;*.csv',
            options=qtw.QFileDialog.DontUseNativeDialog,
        )
        if mol_path:
            for path in mol_path:
                file_format = path.split('.')[-1]
                if file_format in ('txt', 'smi', 'csv'):
                    if file_format == 'csv':
                        df = pd.read_csv(path, sep=';')
                        columns = [c.lower() for c in df.columns.values.tolist()]
                        df.dropna(subset=['Smiles'], inplace=True)
                        self.create_molecule(
                            df['Smiles'].values.tolist(),
                            df['ChEMBL ID'].values.tolist()
                        )
                    else:
                        with open(path, 'r') as file:
                            smi_lines = file.readlines()
                            self.create_molecule(smi_lines)
                else:
                    self.create_molecule(path=path)

    def send_to_worker(
        self, function, sequence, workflow, *,
        mapping=False, kwargs={}
    ):
        self.worker = GenericWorker(
            function, sequence=sequence,
            mapping=mapping, **kwargs
        )
        self.thread = qtc.QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.workflow.connect(workflow)
        self.thread.start()

    def create_molecule(self, smiles=[], names=[], path=False):
        kwargs = {}
        mapping = True
        smi_lines = smiles
        if path:
            smi_lines = [path]
            kwargs = {'file_format': path.split('.')[-1]}
            mapping = False
        elif names:
            smi_lines = [
                f'{smiles} {name}' for smiles, name in \
                zip(
                    df['Smiles'].values.tolist(),
                    df['ChEMBL ID'].values.tolist()
                )
            ]
        self.send_to_worker(
            Molecule, sequence=smi_lines,
            workflow=self.mol_workflow,
            mapping=mapping, kwargs=kwargs
        )

    def create_project(self):
        name = self.uiProjectNameLine.text()
        exists = self.database.check('projects', name)
        if exists:
            qtw.QMessageBox.critical(
                self, 'Existing project', f'The project with name: {name} already exists in the database. It shall not be added.'
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

    def store_molecule(self, molecule):
        if molecule.mol:
            exists = self.database.check('molecules', molecule.inchi_key)
            if not exists:
                self.database.set(
                    'molecules', molecule.inchi_key, molecule
                )
            else:
                qtw.QMessageBox.critical(
                    self, 'Molécula existente', f'La molécula con Inchy key: {molecule.inchi_key} ya existe en la base de datos. No será agregada.'
                )

    def remove_molecule(self):
        # hay que colorear de rojo los trabajos de moléculas que ya no estén en
        # la bd. Si están pendientes, no correr tampoco
        if self.selected_mol:
            for m in self.selected_mol:
                in_project = [p for p in self.projects_list if m in p.molecules]
                if in_project:
                    names = ', '.join([p.name for p in in_project])
                    qtw.QMessageBox.critical(
                        self, 'Molecule in project',
                        f'The molecule {m.get_name} is part of the project(s) {names}. It will not be removed'
                    )
                    continue
                self.database.remove('molecules', m.inchi_key)
                self.database.commit()
                shutil.rmtree(f'molecules/{m.inchi_key}/', ignore_errors=True)
                self.set_models()

    def remove_group_calc(self):
        if self.selected_project is not None:
            self.selected_project.remove_calculation()
            self.database.commit()
            self.set_models()

    def remove_project(self):
        if self.selected_project is not None:
            reply = qtw.QMessageBox.question(
                self, 'Remove a project',
                f'All the calculations made in this project will be removed as well. Do you want to continue removing {self.selected_project.name} project',
                qtw.QMessageBox.Yes | qtw.QMessageBox.No
            )
            if reply == qtw.QMessageBox.Yes:
                for job_id in self.selected_project.job_ids:
                    job = self.database.get('jobs', job_id)
                    if not job:
                        continue
                    molecule = self.database.get('molecules', job.molecule_id)
                    molecule.remove_calculation(job_id)
                    self.database.remove('jobs', job_id)
                self.database.remove('projects', self.selected_project.name)
                self.database.commit()
                self.set_models()
                shutil.rmtree(
                    f'projects/{self.selected_project.name}/',
                    ignore_errors=True
                )

    def remove_job(self):
        index = self.uiJobsTableView.currentIndex()
        if not index:
            return
        reply = qtw.QMessageBox.question(
            self, 'Remove a run',
            'Are you sure you want to remove this run from the queue?',
            qtw.QMessageBox.Yes | qtw.QMessageBox.No
        )
        if reply == qtw.QMessageBox.Yes:
            new_index = self.uiJobsTableView.model().index(index.row(), 0)
            job_id = int(new_index.data())
            self.database.remove('jobs', job_id)
            self.database.commit()
            self.set_models()

    @qtc.pyqtSlot(dict, str)
    def queue_manager(self, calculation, project_name):
        calculation['id'] = self.database.get_job_id
        mol = self.database.get('molecules', calculation['molecule_id'])
        project = self.database.get('projects', project_name)
        if calculation.get('type') == 'Optimization':
            job = Optimization(**calculation)
            # coords = job.opt_coords
            # mol.set_conformer(coords=coords, method=job.keywords)
        elif calculation.get('type') == 'Docking':
            job = Docking(**calculation)
            job.create_config()
        mol.add_calculation(job.id)
        project.add_job_id(job.id)
        self.database.set('jobs', job.id, job)
        self.set_models()

    @qtc.pyqtSlot(object, float)
    def mol_workflow(self, results, progess):
        if not isinstance(results, list):
            results = [results]
        for molecule in results:
            if not molecule:
                return  # mensaje de error
            exists = self.database.check('molecules', molecule.inchi_key)
            if not exists and molecule.mol:
                self.database.set('molecules', molecule.inchi_key, molecule)
            else:
                qtw.QMessageBox.critical(
                    self, 'Existing molecule', f'The molecule with inchi_key: {molecule.inchi_key} already exists in the database. It shall not be added.'
                )
        self.set_models()
        reply = qtw.QMessageBox.question(
            self, 'Link to a project',
            f'Do you want to link these molecules to a project',
            qtw.QMessageBox.Yes | qtw.QMessageBox.No
        )
        if reply == qtw.QMessageBox.Yes:
            project_name, done = qtw.QInputDialog.getText(
                self, 'Project name',
                'Enter the name of the project:'
            )
            if done:
                project_name = project_name.strip()
                exist = self.database.check('projects', project_name)
                if not exist:
                    self.uiProjectNameLine.setText(project_name)
                    self.create_project()
                self.selected_mol = results
                self.selected_project = [
                    p for p in self.projects_list \
                    if p.name == project_name
                ][0]
                self.include_mol()
        self.set_models()

    @logger.catch
    def start_project(self, log):
        # print(log)
        # escogre si ejecutar en serie o en paralelo
        project = self.selected_project
        has_calculations = project.calculations
        # already_programmed = project.status == 'Programmed'
        if project is not None and has_calculations:
            for index, calc in enumerate(project.calculations):
                programmed = calc.get('status') == 'Programmed'
                if calc['type'] in ('Optimization', 'Energy', 'Frequency')\
                   and not programmed:
                    # esto sería para hacer en paralelo; invertir el if anterior y
                    # el siguiente for para hacer en serie (según moléculas)
                    for molecule in project.molecules:
                        self.gauss_controller = Gaussian(
                            molecule,
                            # keywords=calc['keywords'],
                            # job_type=calc['type'],
                            # coordinates=calc.get('coordinates'),
                            project=project,
                            calc=calc
                        )
                        self.gauss_controller.submitted.connect(self.queue_manager)
                        self.gauss_controller.queue_calculation()
                elif calc['type'] == 'Docking' and not programmed:  # buscar otra forma con menos args
                    for molecule in project.molecules:
                        self.vina_controller = MyVina(
                            molecule,
                            project=project,
                            calc=calc
                        )
                        self.vina_controller.submitted.connect(self.queue_manager)
                        self.vina_controller.queue_calculation()
                project.calculations[index]['status'] = 'Programmed'
                project._p_changed = True
            self.database.commit()
            self.set_models()

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
            self.gauss_controller = Gaussian(self.selected_mol[0])
            self.gauss_controller.submitted.connect(self.queue_manager)
        else:
            self.gauss_controller = Gaussian(
                self.selected_project.molecules[0],  #group=True,
                project=self.selected_project
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
        if not group and self.selected_mol:
            self.vina_controller = MyVina(self.selected_mol[0])
            self.vina_controller.submitted.connect(self.queue_manager)
        elif group and self.selected_project:
            self.vina_controller = MyVina(
                self.selected_project.molecules[0],
                project=self.selected_project
            )
            self.vina_controller.submitted.connect(self.set_group_calc)

    @qtc.pyqtSlot(dict, str)
    def set_group_calc(self, calculation, project_name):
        self.selected_project.add_calculation(calculation)
        self.database.commit()
        self.set_models()

    def resume_queue(self):
        if self.master_queue:
            reply = qtw.QMessageBox.question(
                self, 'Work pending in the queue',
                'There are some works pending in the queue, do you want to resume?',
                qtw.QMessageBox.Yes | qtw.QMessageBox.No
            )
            if reply == qtw.QMessageBox.Yes:
                self.start_master_queue()
            else:
                self.queue_status = 'Paused'

    def pause_master_queue(self):
        self.queue_status = 'Paused'

    @qtc.pyqtSlot(int, str)
    def workflow(self, job_id, status):
        job = self.database.get('jobs', job_id)
        job.set_status(status)
        if status == 'Finished' and job.type == 'Optimization':
            mol = self.database.get('molecules', job.molecule_id)
            mol.add_conf_from_opt_file(job.output_file, job.keywords)
        self.database.commit()
        message = f'{job.type}: {job.molecule} -> {status}'
        self.statusBar().showMessage(message)
        self.set_models()

    def start_master_queue(self):
        next_job = self.master_queue[0]
        if next_job.type == 'Docking':
            mol = self.database.get('molecule', job.molecule_id)
            center = (
                job.config['center_x'],
                job.config['center_y'],
                job.config['center_z']
            )
            mol.ligand_prep(
                center, job.receptor_name, job.conformer
            )
        self.uiStartQueueButton.setEnabled(False)
        self.uiPauseQueueButton.setEnabled(True)
        self.worker = Worker(next_job)
        self.queue_thread = qtc.QThread()
        self.worker.moveToThread(self.queue_thread)
        self.queue_thread.started.connect(self.worker.start_queue)
        self.worker.finished.connect(self.queue_thread.quit)
        # self.worker.finished.connect(
            # lambda: self.uiStartQueueButton.setEnabled(True)
        # )
        self.worker.workflow.connect(self.workflow)
        self.queue_thread.start()
        self.queue_status = 'Running'

    def plot_docking_results(self, dock_type):
        if dock_type == 'docking':
            projects = [
                p for p in self.projects_list if p.calculations \
            ]
            self.docking_plotter = DockingPlotter(projects, self.jobs_list)
        elif dock_type == 'redocking':
            projects = [
                p for p in self.projects_list if p.calculations \
                and p.calculations[0].get('redocking')
            ]
            self.docking_plotter = RedockingPlotter(projects, self.jobs_list)

    def ChEMBL_explorer(self):
        self.chembl_explorer = ChEMBLExplorer(self)
        self.chembl_explorer.show()

    def display_data(self, kind):
        if kind == 'fingerprints':
            self.display = FPExplorer(self.selected_project)
        elif kind == 'similarities':
            self.display = SimilarityExplorer(self.selected_project)
        elif kind == 'descriptors':
            self.display = DescriptorsExplorer(self.selected_project)
        self.display.show()

    def projects_right_click(self, position):
        if self.selected_project:
            menu = qtw.QMenu()
            calculation = qtw.QMenu('Group calculation')
            visualize = qtw.QMenu('Visualize')
            menu.addMenu(calculation)
            menu.addMenu(visualize)
            gauss = calculation.addAction('Gaussian')
            vina = calculation.addAction('Vina')
            local_react = calculation.addAction('Local reactivity')
            docking_results = visualize.addAction('Docking results')
            redocking_results = visualize.addAction('Redocking results')
            descriptors = visualize.addAction('Descriptors')
            fingerprints = visualize.addAction('Fingerprints')
            similarities = visualize.addAction('Similarities')
            action = menu.exec_(self.uiProjectsTreeView.mapToGlobal(position))
            if action == gauss and self.selected_project is not None:
                self.gaussian_setup(group=True)
            elif action == vina and self.selected_project is not None:
                self.vina_setup(group=True)
            elif action == local_react:
                pass
            elif action == docking_results:
                self.plot_docking_results('docking')
            elif action == redocking_results:
                self.plot_docking_results('redocking')
            elif action == fingerprints:
                self.display_data('fingerprints')
            elif action == descriptors:
                self.display_data('descriptors')
            elif action == similarities:
                self.display_data('similarities')

    def molecules_right_click(self, position):
        if self.selected_mol[-1]:
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
        # if hasattr(self, 'worker'):
            # self.worker.pause()
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

