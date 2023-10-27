#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import matplotlib
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import QtCore as qtc
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
# from loguru import logger
from scipy.stats import tukey_hsd
from string import ascii_letters as letters
from Worker import GenericWorker
from threading import active_count
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
    def __init__(self, parent, projects, jobs):
        super().__init__()
        uic.loadUi('Views/uiDockingResults.ui', self)
        self.parent = parent
        self.parent.closed.connect(self.close)
        self.projects = projects
        self.jobs = jobs
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        self.selected_jobs = []
        self.selected_projects = []
        self.clusters = {}
        self.to_plot = pd.DataFrame({})
        self.sig_diff_letters = False
        items = [
            '',
            'min to max',
            'max to min',
            'clusters'
        ]
        self.uiSortedCombo.addItems(items)
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
        model = PandasModel(self.df, chunk=0)
        self.uiResultsTableView.setModel(model)
        self.uiResultsTableView.resizeColumnsToContents()
        self.uiResultsTableView.resizeRowsToContents()
        self.sel_to_plot()

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
            self.sel_to_plot()
        else:
            self.selected_jobs = []

    def create_df(self):
        dfs = []
        clusters = {}
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
            clusters[project.name] = project.clusters.copy()
        self.df = pd.concat(dfs) if dfs else pd.DataFrame([])
        self.clusters = clusters
        self.set_model()

    def sel_to_plot(self):
        if not self.selected_jobs or not self.df.shape[0]:
            self.uiCompareProjectButton.setEnabled(False)
            self.uiCompareReceptorButton.setEnabled(False)
            return
        self.sorted = self.uiSortedCombo.currentIndex()
        self.to_plot = self.df.iloc[self.selected_jobs]
        receptors = pd.unique(self.to_plot['Receptor'])
        if len(receptors) > 1:
            self.uiCompareReceptorButton.setEnabled(True)
        else:
            self.uiCompareReceptorButton.setEnabled(False)
        projects = pd.unique(self.to_plot['Project'])
        if len(projects) > 1:
            self.uiCompareProjectButton.setEnabled(True)
        else:
            self.uiCompareProjectButton.setEnabled(False)
        if self.sorted == 1:
            self.to_plot.loc[:, 'Medians'] = self.to_plot.loc[:, 'Binding energies'].apply(np.median)
            self.to_plot.sort_values(by='Medians', inplace=True)
        elif self.sorted == 2:
            self.to_plot.loc[:, 'Medians'] = self.to_plot.loc[:, 'Binding energies'].apply(np.median)
            self.to_plot.sort_values(by='Medians', ascending=False, inplace=True)
        elif self.sorted == 3:
            self.to_plot.loc[:, 'Medians'] = self.to_plot.loc[:, 'Binding energies'].apply(np.median)
            if len(projects) == 1:
                clusters = self.clusters[projects[0]]
                sorted_indexes = clusters.get('sorted_indexes')
                if len(sorted_indexes) == len(self.to_plot):
                    self.to_plot = self.to_plot.iloc[sorted_indexes]
        if self.to_plot.shape[0] > 1:
            self.parent.send_to_worker(
                tukey_hsd,
                sequence=self.to_plot.loc[:, 'Binding energies'].tolist(),
                workflow=self.workflow,
                mapping=False,
                unpack=True
            )
        self.uiSigDiffCheck.setEnabled(False)
        self.plot()

    @qtc.pyqtSlot(object, float)
    def workflow(self, results, process):
        diff_array = [i for i in map(self.diff_tukey, results.pvalue)]
        diff_matrix = np.asmatrix(diff_array)
        diff_matrix_tril = np.tril(diff_matrix)
        self.sig_diff_letters = self.sig_diff(diff_matrix_tril)
        print(process)
        self.uiSigDiffCheck.setEnabled(True)

    def set_sig_diff(self, checked):
        if self.to_plot.shape[0] == len(self.sig_diff_letters) and checked:
            self.plot(self.sig_diff_letters)
        else:
            self.plot()

    def plot(self, sig_diff_letters=False):
        self.canvas.ax.clear()
        self.canvas.ax.set_title('Docking')
        total = self.to_plot.shape[0]
        min_value = min([min(i) for i in self.to_plot['Binding energies']])
        max_value = max([max(i) for i in self.to_plot['Binding energies']])
        avg_value = sum([sum(i) / len(i) for i in self.to_plot['Binding energies']]) / len(self.to_plot['Binding energies'])
        self.uiMinValueLabel.setText(f'Min: {min_value:.2f}')
        self.uiMaxValueLabel.setText(f'Max: {max_value:.2f}')
        self.uiAvgValueLabel.setText(f'Avg: {avg_value:.2f}')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.axhline(min_value, linewidth=.5, color='green')
        self.canvas.ax.axhline(avg_value, linewidth=.5, color='blue', linestyle='--')
        self.canvas.ax.axhline(max_value, linewidth=.5, color='red')
        positions = np.arange(total)
        box_patches = self.canvas.ax.boxplot(
            positions=positions,
            x=self.to_plot['Binding energies'],
            labels=self.to_plot['Molecule'],
            patch_artist=True,
            notch=True
        )
        self.customize_boxes(
            box_patches, self.canvas.ax,
            max_value, sig_diff_letters, color='blue'
        )
        self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
        self.canvas.ax.set_ylim(bottom=min_value - 2, top=max_value + 2)
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.draw()

    def customize_boxes(self, patches, ax, max_val, sig_diff_letters, color='black'):
        for box in patches['boxes']:
            box.set(
                facecolor='white',
                edgecolor=color,
                # alpha=0.5,
                linewidth=1.0,
            )
        for median in patches['medians']:
            (x_l, y), (x_r, _) = median.get_xydata()
            x_center = (x_r + x_l) / 2
            ax.text(
                x_center, max_val + 1, f'{y:.2f}',
                fontsize=10, color=color,
                verticalalignment='bottom',
                horizontalalignment='center',
                rotation=90
            )
        if sig_diff_letters:
            for i, whisker in enumerate(patches['whiskers']):
                if i % 2:
                    bottom, top = whisker.get_xydata()
                    offset = (top[1] - bottom[1])
                    ax.text(
                        top[0], top[1] + offset, sig_diff_letters[(i - 1) // 2],
                        fontsize=15, color=color,
                        verticalalignment='bottom',
                        horizontalalignment='center',
                    )

    def compare_project(self):
        self.canvas.ax.clear()
        self.canvas.ax.set_title('Docking')
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        colors = ['blue', 'red', 'green', 'magenta', 'orange']
        to_plot_both = self.df.iloc[self.selected_jobs]
        projects = pd.unique(to_plot_both['Project'])
        if len(projects) != 2:
            return  # here a critical warning
        total = to_plot_both.shape[0]
        both_df = []
        for i, project in enumerate(projects):
            mask = to_plot_both['Project'] == project
            to_plot = to_plot_both.loc[mask, ['Binding energies', 'Molecule']]
            both_df.append(to_plot)
        molecules_p1 = both_df[0]['Molecule'].values.tolist()
        molecules_p2 = both_df[1]['Molecule'].values.tolist()
        if both_df[0].shape[0] != both_df[1].shape[0]:
            return  # avisar que no tienen el mismo número de compuestos
        elif molecules_p1 != molecules_p2:
            return  # avisar que no están en el mismo orden
        for i, project in enumerate(projects):
            total = to_plot.shape[0]
            positions = np.arange(total) * 2 - 0.4 * (-1) ** (i + 1)
            self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
            self.canvas.ax.plot([], label=projects[i])
            box_patches = self.canvas.ax.boxplot(
                positions=positions,
                x=both_df[i].loc[:, 'Binding energies'],
                labels=molecules_p1,
                patch_artist=True,
                notch=True
            )
            self.canvas.ax.set_xticks(np.arange(total) * 2)
            self.customize_boxes(box_patches, self.canvas.ax, color=colors[i])
        self.canvas.ax.legend()
        self.canvas.draw()
        x1 = both_df[0]['Binding energies'].apply(np.median).values.tolist()
        x2 = both_df[1]['Binding energies'].apply(np.median).values.tolist()
        plt.figure(figsize=(10, 10))
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
        colors = ['blue', 'red', 'green', 'magenta', 'orange']
        to_plot_both = self.df.iloc[self.selected_jobs]
        receptors = pd.unique(to_plot_both['Receptor'])
        if len(receptors) != 2:
            return  # here a critical warning
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
            return  # avisar que no tienen el mismo número de compuestos
        elif molecules_p1 != molecules_p2:
            return  # avisar que no están en el mismo orden
        max_values = []
        min_values = []
        avg_values = []
        patches = []
        for i, receptor in enumerate(receptors):
            total = both_df[i].shape[0]
            positions = np.arange(total) * 2 - 0.4 * (-1) ** (i + 1)
            self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
            self.canvas.ax.plot([], label=receptors[i])
            box_patches = self.canvas.ax.boxplot(
                positions=positions,
                x=both_df[i].loc[:, 'Binding energies'],
                labels=molecules_p1,
                patch_artist=True,
                notch=True
            )
            patches.append(box_patches)
            self.canvas.ax.set_xticks(np.arange(total) * 2)
            max_values.append(max([max(i) for i in both_df[i]['Binding energies']]))
            min_values.append(min([min(i) for i in both_df[i]['Binding energies']]))
            avg_values.append(
                sum([sum(i) / len(i) for i in both_df[i]['Binding energies']]) / len(both_df[i]['Binding energies'])
            )
        min_value = min(min_values)
        max_value = max(max_values)
        avg_value = sum(avg_values) / len(avg_values)
        self.uiMinValueLabel.setText(f'Min: {min_value:.2f}')
        self.uiMaxValueLabel.setText(f'Max: {max_value:.2f}')
        self.uiAvgValueLabel.setText(f'Avg: {avg_value:.2f}')
        for i, receptor in enumerate(receptors):
            self.customize_boxes(
                patches[i], self.canvas.ax, max_value,
                self.sig_diff_letters, color=colors[i]
            )
        self.canvas.ax.set_ylim(bottom=min_value - 2, top=max_value + 2)
        self.canvas.ax.legend()
        self.canvas.draw()
        # if working with all data is wanted
        # both_df = [df['Binding energies'].apply(np.sort) for df in both_df]
        # x1 = np.concatenate(both_df[0].values.tolist())
        # x2 = np.concatenate(both_df[1].values.tolist())
        # if working with the medians of each group of data is wanted
        x1 = both_df[0]['Binding energies'].apply(np.median).values.tolist()
        x2 = both_df[1]['Binding energies'].apply(np.median).values.tolist()
        plt.figure(figsize=(5, 5))
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

    def diff_tukey(self, arr):
        '''
        arr:
            numpy range.
        return:
            a binary list. 1 if there are significative differences, 0 if not.
        '''
        return [i for i in map(lambda x: 1 if abs(x) < 0.05 else 0, arr)]

    def sig_diff(self, diff_matrix):
        sig_diff_letters = []
        last = ''
        for i, row in enumerate(diff_matrix):
            if i == 0:
                sig_diff_letters.append(letters[i])
                last = letters[i]
                continue
            this_row = []
            for j, col in enumerate(row):
                if j < i:
                    to_add = sig_diff_letters[j]
                    if col == 0 and len(to_add) == 1:
                        this_row.append(to_add)
                    else:
                        continue
                else:
                    break
            if not this_row:
                idx_last_letter = letters.index(last[-1])
                last = letters[idx_last_letter + 1]
                this_row.append(last)
            this_row = ''.join(sorted(set(this_row)))
            sig_diff_letters.append(this_row)
        return sig_diff_letters


class RedockingPlotter(qtw.QWidget):

    def __init__(self, parent, projects, jobs):
        super().__init__()
        uic.loadUi('Views/uiRedockingResults.ui', self)
        self.parent = parent
        self.parent.closed.connect(self.close)
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
        if not hasattr(self, 'ax2'):
            self.ax2 = self.canvas.ax.twinx()
        self.ax2.clear()
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.tick_params(axis='y', colors='red')
        self.ax2.tick_params(axis='y', colors='blue')
        # self.canvas.ax.set_title(f'Redocking')
        to_plot = self.df.iloc[self.selected_jobs]
        min_value_energy = min([min(i) for i in to_plot['Binding energies']])
        max_value_energy = max([max(i) for i in to_plot['Binding energies']])
        avg_value_energy = sum([sum(i) / len(i) for i in to_plot['Binding energies']]) / len(to_plot['Binding energies'])
        # avg_value_energy = (min_value_energy + max_value_energy) / 2
        min_value_rmsd = min([min(i) for i in to_plot['RMSD']])
        max_value_rmsd = max([max(i) for i in to_plot['RMSD']])
        total = to_plot.shape[0]
        positions = np.arange(total)
        positions_2 = np.arange(total) * -1
        positions_2 = positions_2 + -1
        positions_2 = positions_2[::-1]
        positions_3 = np.arange(total) * 1
        positions_3 = positions_3 + 1
        # self.ax2.axhline(min_value_energy, linewidth=.5, color='green')
        # self.ax2.axhline(avg_value_energy, linewidth=.5, color='blue', linestyle='--')
        # self.ax2.axhline(max_value_energy, linewidth=.5, color='red')
        receptor_name = to_plot.loc[0, 'Receptor']
        rmsd_patches = self.canvas.ax.boxplot(
            to_plot['RMSD'], positions=positions_2,
            patch_artist=True, # to fill with colors
            notch=True,
            labels=[
                f'{mol} ({proj})' for mol, proj in zip(to_plot['Molecule'], to_plot['Project'])
            ]
        )
        self.customize_boxes(rmsd_patches, self.canvas.ax, max_value_rmsd, color='red')
        energy_patches = self.ax2.boxplot(
            x=to_plot['Binding energies'], positions=positions_3,
            notch=True,
            patch_artist=True, # to fill with colors
            labels=[
                f'{mol} ({proj})' for mol, proj in zip(to_plot['Molecule'], to_plot['Project'])
            ]
        )
        self.customize_boxes(energy_patches, self.ax2, max_value_energy, color='blue')
        self.canvas.ax.set_ylabel('RMSD (\u212B)', color='red')
        self.canvas.ax.set_ylim(bottom=0., top=2.5)
        self.ax2.set_ylim(bottom=min_value_energy - 0.6, top=max_value_energy + 1)
        self.ax2.set_ylabel('Binding Energy (kcal/mol)', color='blue')
        self.canvas.draw()

    def customize_boxes(self, patches, ax, max_val, color='black'):
        for box in patches['boxes']:
            box.set(
                facecolor='white',
                edgecolor=color,
                # alpha=0.5,
                linewidth=1.0,
            )
        for median in patches['medians']:
            (x_l, y), (x_r, _) = median.get_xydata()
            x_center = (x_r + x_l) / 2
            ax.text(
                x_center, max_val + 0.3, f'{y:.2f}',
                fontsize=10, color=color,
                verticalalignment='bottom',
                horizontalalignment='center',
                rotation=90
            )

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
        cbar = self.canvas.ax.figure.colorbar(
            im, ax=self.canvas.ax, pad=0.15,
            aspect=50
        )
        # cbar.ax.set_ylabel('Similarity', rotation=-90, va="bottom")
        cbar.ax.set_ylabel('Similarity', rotation=90, va="center")
        # Show all ticks and label them with the respective list entries.
        self.canvas.ax.set_xticks(
            np.arange(simil_matrix.shape[1]), labels=labels,
            rotation=90
        )
        self.canvas.ax.set_yticks(np.arange(simil_matrix.shape[0]), labels=labels)
        # Let the horizontal axes labeling appear on top.
        self.canvas.ax.tick_params(
            top=True, bottom=False,
            labeltop=True, labelbottom=False,
            left=False, right=True,
        )
        self.canvas.ax.yaxis.set_label_position("right")
        self.canvas.ax.yaxis.tick_right()
        self.canvas.ax.tick_params(axis='x', labelsize=8)
        self.canvas.ax.tick_params(axis='y', labelsize=8)
        # Turn spines off and create white grid.
        self.canvas.ax.spines[:].set_visible(False)
        self.canvas.ax.set_xticks(np.arange(simil_matrix.shape[1] + 1) - .5, minor=True)
        self.canvas.ax.set_yticks(np.arange(simil_matrix.shape[0] + 1) - .5, minor=True)
        # self.canvas.ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
        self.canvas.ax.tick_params(
            which="minor", bottom=False, left=False, right=False, top=False
        )
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


