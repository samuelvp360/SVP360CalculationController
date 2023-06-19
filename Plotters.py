#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import matplotlib
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from PIL import Image, ImageQt
import pandas as pd
import numpy as np
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from Models import PandasModel
from sklearn.metrics import r2_score
from loguru import logger
matplotlib.use('Qt5Agg')


class PlotCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, parent=None):
        self.fig = Figure(
            figsize=(12, 8), dpi=100,
            tight_layout=True
        )
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)


class DockingPlotter(qtw.QWidget):

    # @logger.catch
    def __init__(self, projects, jobs):
        super().__init__()
        uic.loadUi('Views/uiDockingResults.ui', self)
        self.projects = projects
        self.jobs = jobs
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        self.selected_jobs = []
        self.selected_projects = []
        self.create_checks()
        self.create_df()
        self.show()

    def create_checks(self):
        self.checks = [
            qtw.QCheckBox(f'{project.name}') for project in self.projects
        ]
        for i, check in enumerate(self.checks):
            check.stateChanged.connect(self.select_project)
            self.uiProjectsLayout.addWidget(check, i)
        spacer = qtw.QSpacerItem(20, 40, qtw.QSizePolicy.Minimum, qtw.QSizePolicy.Expanding)
        self.uiProjectsLayout.addItem(spacer)

    def set_model(self):
        model = PandasModel(self.df)
        self.uiResultsTableView.setModel(model)
        self.uiResultsTableView.resizeColumnsToContents()
        self.uiResultsTableView.resizeRowsToContents()
        self.plot()

    def select_project(self):
        indexes = [p.isChecked() for p in self.checks]
        self.selected_projects = [
            p for check, p in zip(indexes, self.projects) if check
        ]
        self.create_df()

    def select_jobs(self):
        indexes = self.uiResultsTableView.selectedIndexes()
        if indexes:
            self.selected_jobs = [i.row() for i in indexes]
            self.plot()
        else:
            self.selected_jobs = []

    def create_df(self):
        dfs = []
        for project in self.selected_projects:
            if not project.job_ids:
                continue
            names = pd.Series(
                [mol.get_name for mol in project.molecules], name='Names'
            )
            energies = []
            names = []
            receptor = []
            box_size = []
            project_name = []
            for job in self.jobs:
                if job.id in project.job_ids and job.type == 'Docking':
                    if len(job.energies) == 0:
                        continue
                    energies.append(job.energies)
                    names.append(job.molecule)
                    receptor.append(job.receptor_name)
                    size_x = job.config.get('size_x')
                    size_y = job.config.get('size_y')
                    size_z = job.config.get('size_z')
                    box_size.append(f'{size_x}x{size_y}x{size_z}')
                    project_name.append(project.name)
            values = {
                'Molecule': names,
                'Receptor': receptor,
                'Box size (\u212B\u00B3)': box_size,
                'Binding energies': list(energies),
                'Project': project_name
            }
            dfs.append(pd.DataFrame(values))
        self.df = pd.concat(dfs) if dfs else pd.DataFrame([])
        self.set_model()

    def plot(self):
        self.canvas.ax.clear()
        if not self.selected_jobs or not self.df.shape[0]:
            self.uiSortedCheck.setChecked(False)
            self.uiCompareProjectButton.setEnabled(False)
            self.uiCompareReceptorButton.setEnabled(False)
            return
        self.sorted = self.uiSortedCheck.isChecked()
        self.uiSortedCheck.setEnabled(True)
        self.canvas.ax.set_title('Docking')
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.axhline(-7., linewidth=.5, color='red')
        to_plot_both = self.df.iloc[self.selected_jobs]
        receptors = pd.unique(to_plot_both['Receptor'])
        if len(receptors) > 1:
            self.uiCompareReceptorButton.setEnabled(True)
        else:
            self.uiCompareReceptorButton.setEnabled(False)
        projects = pd.unique(to_plot_both['Project'])
        if len(projects) > 1:
            self.uiCompareProjectButton.setEnabled(True)
        else:
            self.uiCompareProjectButton.setEnabled(False)
        total = to_plot_both.shape[0]
        if self.sorted:
            to_plot_both['Medians'] = to_plot_both['Binding energies'].apply(np.median)
            to_plot_both.sort_values(by='Medians', inplace=True)
        props={'color': 'blue', 'linewidth': 1.0}
        positions = np.arange(total)
        self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
        self.canvas.ax.boxplot(
            positions=positions,
            x=to_plot_both['Binding energies'],
            labels=to_plot_both['Molecule'],
            boxprops=props, medianprops=props,
            whiskerprops=props, capprops=props,
        )
        self.canvas.draw()

    def compare_project(self):
        self.canvas.ax.clear()
        self.canvas.ax.set_title('Docking')
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.axhline(-7., linewidth=.5, color='red')
        colors = ['blue', 'red', 'green', 'magenta', 'orange']
        to_plot_both = self.df.iloc[self.selected_jobs]
        projects = pd.unique(to_plot_both['Project'])
        if len(projects) != 2:
            return # here a critical warning
        total = to_plot_both.shape[0]
        both_df = []
        for i, project in enumerate(projects):
            mask = to_plot_both['Project'] == project
            to_plot = to_plot_both.loc[mask, ['Binding energies', 'Molecule']]
            both_df.append(to_plot)
        molecules_p1 = both_df[0]['Molecule'].values.tolist()
        molecules_p2 = both_df[1]['Molecule'].values.tolist()
        if both_df[0].shape[0] != both_df[1].shape[0]:
            return # avisar que no tienen el mismo número de compuestos
        elif molecules_p1 != molecules_p2:
            return # avisar que no están en el mismo orden
        for i, project in enumerate(projects):
            props = {'color': colors[i], 'linewidth': 1.0}
            total = to_plot.shape[0]
            positions = np.arange(total) * 2 - 0.4 * (-1) ** (i + 1)
            self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
            self.canvas.ax.plot([], c=colors[i], label=projects[i])
            self.canvas.ax.boxplot(
                positions=positions,
                x=both_df[i].loc[:, 'Binding energies'],
                labels=molecules_p1,
                boxprops=props, medianprops=props,
                whiskerprops=props, capprops=props,
            )
            self.canvas.ax.set_xticks(np.arange(total) * 2)
        self.canvas.ax.legend()
        self.canvas.draw()
        x1 = both_df[0]['Binding energies'].apply(np.median).values.tolist()
        x2 = both_df[1]['Binding energies'].apply(np.median).values.tolist()
        plt.figure(figsize=(10,10))
        plt.scatter(x1, x2, alpha=0.3)
        plt.title('Linear correlation between projects')
        plt.xlabel(projects[0])
        plt.ylabel(projects[1])
        fit = np.polyfit(x1, x2, 1)
        f = np.poly1d(fit)
        m, b = fit
        r2 = r2_score(x2, f(x1))
        x = np.linspace(min(x1), max(x1), num=100)
        plt.plot(
            x, f(x), c='orange', label='Trend line', linewidth=0.5
        )
        plt.annotate(
            f'r\u00B2 = {r2:.5f}\nm = {m:.5f}\nb = {b:.5f}', (min(x1), max(x2) - 0.5)
        )
        plt.legend(loc='lower right')
        plt.show()

    def compare_receptor(self):
        self.canvas.ax.clear()
        self.canvas.ax.set_title('Docking')
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.axhline(-7., linewidth=.5, color='red')
        colors = ['blue', 'red', 'green', 'magenta', 'orange']
        to_plot_both = self.df.iloc[self.selected_jobs]
        receptors = pd.unique(to_plot_both['Receptor'])
        if len(receptors) != 2:
            return # here a critical warning
        total = to_plot_both.shape[0]
        both_df = []
        for i, receptor in enumerate(receptors):
            mask = to_plot_both['Receptor'] == receptor
            to_plot = to_plot_both.loc[mask, ['Binding energies', 'Molecule']]
            to_plot.reset_index(drop=True, inplace=True)
            both_df.append(to_plot)
        molecules_p1 = both_df[0]['Molecule'].values.tolist()
        molecules_p2 = both_df[1]['Molecule'].values.tolist()
        if both_df[0].shape[0] != both_df[1].shape[0]:
            return # avisar que no tienen el mismo número de compuestos
        elif molecules_p1 != molecules_p2:
            return # avisar que no están en el mismo orden
        for i, receptor in enumerate(receptors):
            props = {'color': colors[i], 'linewidth': 1.0}
            total = both_df[i].shape[0]
            positions = np.arange(total) * 2 - 0.4 * (-1) ** (i + 1)
            self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
            self.canvas.ax.plot([], c=colors[i], label=receptors[i])
            self.canvas.ax.boxplot(
                positions=positions,
                x=both_df[i].loc[:, 'Binding energies'],
                labels=molecules_p1,
                boxprops=props, medianprops=props,
                whiskerprops=props, capprops=props,
            )
            self.canvas.ax.set_xticks(np.arange(total) * 2)
        self.canvas.ax.legend()
        self.canvas.draw()
        # if working with all data is wanted
        # both_df = [df['Binding energies'].apply(np.sort) for df in both_df]
        # x1 = np.concatenate(both_df[0].values.tolist())
        # x2 = np.concatenate(both_df[1].values.tolist())
        # if working with the medians of each group of data is wanted
        x1 = both_df[0]['Binding energies'].apply(np.median).values.tolist()
        x2 = both_df[1]['Binding energies'].apply(np.median).values.tolist()
        plt.figure(figsize=(5,5))
        plt.scatter(x1, x2, alpha=0.3, label='Medians by molecule')
        plt.title('Linear correlation between receptors')
        plt.xlabel(receptors[0])
        plt.ylabel(receptors[1])
        fit = np.polyfit(x1, x2, 1)
        f = np.poly1d(fit)
        m, b = fit
        r2 = r2_score(x2, f(x1))
        x = np.linspace(min(x1), max(x1), num=100)
        plt.plot(
            x, f(x), c='orange', label='Trend line', linewidth=0.5
        )
        plt.annotate(
            f'r\u00B2 = {r2:.5f}\nm = {m:.5f}\nb = {b:.5f}', (min(x1), max(x2) - 0.5)
        )
        plt.legend(loc='lower right')
        plt.show()

    def export_to_excel(self):
        filename, _ = qtw.QFileDialog.getSaveFileName(
            self, 'Save results as excel data sheet',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Excel file (*.xlsx *.xls)'
        )
        if filename:
            if '.xlsx' not in filename or '.xls' not in filename:
                self.df.to_excel(filename + '.xlsx')
            else:
                self.df.to_excel(filename)

    def import_from_excel(self):
        filename, _ = qtw.QFileDialog.getOpenFileName(
            self, 'Save results as excel data sheet',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Excel file (*.xlsx *xls)'
        )
        if filename:
            new_df = pd.read_excel(filename)
            new_df['Binding energies'] = new_df['Binding energies'].apply(
                lambda x: np.array(
                    [float(i) for i in x.replace('[', '').replace(']', '').split()]
            ))
            new_df = new_df.iloc[:, 1:]
            self.df = pd.concat([self.df, new_df], ignore_index=True)
            self.set_model()


class RedockingPlotter(qtw.QWidget):

    def __init__(self, projects, jobs):
        super().__init__()
        uic.loadUi('Views/uiRedockingResults.ui', self)
        self.projects = projects
        self.jobs = jobs
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        self.selected_jobs = []
        self.selected_projects = []
        self.create_checks()
        self.create_df()
        self.show()

    def create_checks(self):
        self.checks = [
            qtw.QCheckBox(f'{project.name}') for project in self.projects
        ]
        for i, check in enumerate(self.checks):
            check.stateChanged.connect(self.select_project)
            self.uiProjectsLayout.addWidget(check, i)
        spacer = qtw.QSpacerItem(20, 40, qtw.QSizePolicy.Minimum, qtw.QSizePolicy.Expanding)
        self.uiProjectsLayout.addItem(spacer)

    def create_df(self):
        dfs = []
        for project in self.selected_projects:
            if not project.job_ids:
                continue
            names = pd.Series(
                [mol.get_name for mol in project.molecules], name='Names'
            )
            energies = []
            names = []
            receptor = []
            box_size = []
            rmsd = []
            project_name = []
            for job in self.jobs:
                if job.id in project.job_ids:
                    if len(job.energies) == 0:
                        continue
                    energies.append(job.energies)
                    names.append(job.molecule)
                    receptor.append(job.receptor_name)
                    size_x = job.config.get('size_x')
                    size_y = job.config.get('size_y')
                    size_z = job.config.get('size_z')
                    box_size.append(f'{size_x}x{size_y}x{size_z}')
                    project_name.append(project.name)
                    rmsd.append(job.rmsd)
            values = {
                'Molecule': names,
                'Receptor': receptor,
                'Box size (\u212B\u00B3)': box_size,
                'Binding energies': list(energies),
                'Project': project_name,
                'RMSD': rmsd
            }
            dfs.append(pd.DataFrame(values))
        self.df = pd.concat(dfs) if dfs else pd.DataFrame([])
        self.set_model()

    def set_model(self):
        model = PandasModel(self.df)
        self.uiResultsTableView.setModel(model)
        self.uiResultsTableView.resizeColumnsToContents()
        self.uiResultsTableView.resizeRowsToContents()
        self.plot()

    def plot(self):
        self.canvas.ax.clear()
        if not self.selected_jobs or not self.df.shape[0]:
            return
        self.ax2 = self.canvas.ax.twinx()
        self.ax2.clear()
        self.canvas.ax.set_ylabel('RMSD')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.set_ylim(bottom=0., top=2.5)
        props={'color': 'blue', 'linewidth': 1.0}
        to_plot = self.df.iloc[self.selected_jobs]
        self.canvas.ax.set_xlabel(to_plot['Molecule'])
        receptor_name = to_plot.loc[0, 'Receptor']
        self.canvas.ax.set_title(f'Docking on {receptor_name}')
        self.canvas.ax.boxplot(
            x=to_plot['RMSD'], positions=[-0.5],
            # labels=to_plot['Molecule'],
            boxprops=props, medianprops=props,
            whiskerprops=props, capprops=props,
        )
        self.ax2.set_ylabel('Binding Energy (kcal/mol)')
        props={'color': 'red', 'linewidth': 1.0}
        self.ax2.boxplot(
            x=to_plot['Binding energies'], positions=[0.5],
            labels=to_plot['Molecule'],
            boxprops=props, medianprops=props,
            whiskerprops=props, capprops=props
        )
        self.canvas.draw()

    def select_jobs(self):
        indexes = self.uiResultsTableView.selectedIndexes()
        if indexes:
            self.selected_jobs = [i.row() for i in indexes]
            self.plot()
        else:
            self.selected_jobs = []

    def select_project(self):
        indexes = [p.isChecked() for p in self.checks]
        self.selected_projects = [
            p for check, p in zip(indexes, self.projects) if check
        ]
        self.create_df()


class Heatmap(qtw.QWidget):

    def __init__(self, data, fp_type, simil_metric):
        super().__init__()
        uic.loadUi('Views/uiHeatmap.ui', self)
        self.data = data
        self.fp_type = fp_type
        self.simil_metric = simil_metric
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)

    def plot(self):
        self.canvas.ax.clear()
        simil_matrix = self.data.to_numpy()
        labels = self.data.columns.values.tolist()
        self.canvas.ax.set_title(
            f'Similarity Matrix ({self.fp_type}; {self.simil_metric})'
        )
        im = self.canvas.ax.imshow(simil_matrix, cmap='magma_r')
        # Create colorbar
        cbar = self.canvas.ax.figure.colorbar(im, ax=self.canvas.ax)
        cbar.ax.set_ylabel('', rotation=-90, va="bottom")
        # Show all ticks and label them with the respective list entries.
        self.canvas.ax.set_xticks(
            np.arange(simil_matrix.shape[1]), labels=labels,
            rotation=90
        )
        self.canvas.ax.set_yticks(np.arange(simil_matrix.shape[0]), labels=labels)
        # Let the horizontal axes labeling appear on top.
        self.canvas.ax.tick_params(
            top=True, bottom=False,
            labeltop=True, labelbottom=False
        )
        self.canvas.ax.tick_params(axis='x', labelsize=8)
        self.canvas.ax.tick_params(axis='y', labelsize=8)
        # Turn spines off and create white grid.
        self.canvas.ax.spines[:].set_visible(False)
        self.canvas.ax.set_xticks(np.arange(simil_matrix.shape[1] + 1) - .5, minor=True)
        self.canvas.ax.set_yticks(np.arange(simil_matrix.shape[0] + 1) - .5, minor=True)
        self.canvas.ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
        self.canvas.ax.tick_params(which="minor", bottom=False, left=False)
        self.canvas.draw()


class SimilMap(qtw.QWidget):

    def __init__(self, ref, prob, fp_type, ref_name='NN', prob_name='NN'):
        super().__init__()
        uic.loadUi('Views/uiSimilMap.ui', self)
        self.img = self.get_simil_map(ref, prob, fp_type)
        qimage = qtg.QPixmap.fromImage(ImageQt.ImageQt(self.img))
        self.uiSimilMapLabel.setPixmap(qimage)
        self.uiLabel.setText(
            f'Ref: {ref_name}\tProb: {prob_name}\nFP type: {fp_type}'
        )
        self.show()

    def save_image(self):
        filename, extension = qtw.QFileDialog.getSaveFileName(
            self, 'Save image', 'Image.png', 'PNG (*.png)',
            options=qtw.QFileDialog.DontUseNativeDialog,
        )
        if filename:
            has_extension = '.' in filename
            if not has_extension:
                filename += '.png'
            self.img.save(filename)

    @classmethod
    def get_simil_map(self, ref, prob, fp_type, size=(550, 550)):
        d = Draw.MolDraw2DCairo(*size)
        if 'Morgan' in fp_type:
            radius = int(fp_type[-1])
            simil_function = lambda m, i: SimilarityMaps.GetMorganFingerprint(
                m, i, radius=radius, fpType='bv'
            )
        elif fp_type == 'Hasehd Topological Torsions':
            simil_function = lambda m, i: SimilarityMaps.GetTTFingerprint(
                m, i, fpType='bv'
            )
        elif 'Hashed Atom Pairs' in fp_type:
            simil_function = lambda m, i: SimilarityMaps.GetAPFingerprint(
                m, i, fpType='bv'
            )
        elif fp_type == 'RDKit':
            simil_function = lambda m, i: SimilarityMaps.GetRDKFingerprint(
                m, i, fpType='bv'
            )
        d.DrawMolecule(ref)
        _, maxWeight = SimilarityMaps.GetSimilarityMapForFingerprint(
            ref, prob, simil_function, draw2d=d
        )
        d.FinishDrawing()
        data = d.GetDrawingText()
        bio = io.BytesIO(data)
        return Image.open(bio)


