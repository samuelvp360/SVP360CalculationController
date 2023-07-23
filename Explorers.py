#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import pandas as pd
import numpy as np
from functools import reduce
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from PIL import Image, ImageQt
from Models import PandasModel, DisplayMFPModel, SimilarityModel
from Components import Similarity, MyChembl, Molecule
from Plotters import PlotCanvas, Heatmap, SimilMap
from Worker import GenericWorker
from django.core.paginator import Paginator
from rdkit import Chem
from rdkit.Chem import Draw
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, leaves_list
from sklearn.metrics import silhouette_samples, silhouette_score
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from loguru import logger


class FPExplorer(qtw.QWidget):

    def __init__(self, project):
        super().__init__()
        uic.loadUi('Views/uiDisplayFP.ui', self)
        self.project = project
        fp_types = self.project.fps[0].keys()
        self.uiFPTypeCombo.addItems(fp_types)
        self.uiFPTypeCombo.setCurrentText('Morgan2')

    def change_fp_type(self, fp_type):
        if fp_type:
            self.fp_type = fp_type
            self.df = self.create_df(fp_type)
            self.set_model()

    def create_df(self, fp_type):
        fps = pd.DataFrame([fp.get(fp_type) for fp in self.project.fps])
        fps.set_axis([str(i) for i in range(2048)], axis=1, inplace=True)
        fps.set_axis(
            [m.get_name for m in self.project.molecules], axis=0, inplace=True
        )
        fps.replace(0, np.nan, inplace=True)
        fps.dropna(how='all', axis=1, inplace=True)
        fps.replace(np.nan, 0, inplace=True)
        fps = fps.applymap(int)
        return fps

    def set_model(self):
        model = DisplayMFPModel(self.df)
        self.uiFPTable.setModel(model)
        self.uiFPTable.resizeColumnsToContents()
        self.uiFPTable.resizeRowsToContents()

    def draw_qimage(self, mol, mol_size=(350, 350), **kwargs):
        drawer = Draw.MolDraw2DCairo(mol_size[0], mol_size[1])
        mc = Chem.Mol(mol.ToBinary())
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
        if not mc.GetNumConformers():
            Chem.rdDepictor.Compute2DCoords(mc)
        drawer.DrawMolecule(mc, **kwargs)
        drawer.FinishDrawing()
        data = drawer.GetDrawingText()
        bio = io.BytesIO(data)
        self.img = Image.open(bio)
        qimage = qtg.QPixmap.fromImage(ImageQt.ImageQt(self.img))
        return qimage

    def depict_fp(self, index):
        has_bit = self.df.iloc[index.row(), index.column()]
        if not has_bit or not self.fp_type.startswith('Morgan'):
            self.uiFPLabel.setPixmap(qtg.QPixmap(
                self.project.molecules[index.row()].mol_img
            ))
            return
        all_bits = self.df.columns.values.tolist()
        bit_id = int(all_bits[index.column()])
        molecule = self.project.molecules[index.row()].get_2d_mol
        fp_radius = int(self.fp_type[-1])
        info = {}
        fp = Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(
            molecule, fp_radius, 2048, bitInfo=info
        )
        this_bit = info.get(bit_id)
        if not this_bit:
            self.uiFPLabel.setPixmap(qtg.QPixmap(
                self.project.molecules[index.row()].mol_img
            ))
            return
        atom_id, radius = this_bit[0]
        if radius > 0:
            env = Chem.FindAtomEnvironmentOfRadiusN(molecule, radius , atom_id)
            atoms_to_use = []
            for b in env:
                atoms_to_use.append(molecule.GetBondWithIdx(b).GetBeginAtomIdx())
                atoms_to_use.append(molecule.GetBondWithIdx(b).GetEndAtomIdx())
            atoms_to_use = list(set(atoms_to_use))
        else:
            atoms_to_use = [atom_id]
            env = None
        qimage = self.draw_qimage(
            molecule, highlightAtoms=atoms_to_use,
            highlightAtomColors={atom_id: (0.3, 0.3, 1)}
        )
        self.uiFPLabel.setPixmap(qimage)


class DescriptorsExplorer(qtw.QWidget):

    def __init__(self, project):
        super().__init__()
        uic.loadUi('Views/uiDisplayDescriptors.ui', self)
        self.project = project
        self.set_model()

    def set_model(self):
        self.df = self.project.descriptors
        model = PandasModel(self.df)
        self.uiDescriptorsTable.setModel(model)
        self.uiDescriptorsTable.resizeColumnsToContents()
        self.uiDescriptorsTable.resizeRowsToContents()


class SimilarityExplorer(qtw.QWidget):

    def __init__(self, parent, project):
        super().__init__()
        uic.loadUi('Views/uiSimilarity.ui', self)
        self.parent = parent
        self.parent.closed.connect(self.close)
        self.project = project
        self.names = [m.get_name for m in project.molecules]
        self.molecules = project.molecules
        self.similarity = Similarity(project.fps, self.names)
        self.simil_metric = 'Tanimoto'
        self.fp_type = 'Morgan2'
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        simil_metrics = (
            'Tanimoto', 'Dice', 'Cosine', 'Sokal',
            'Russel', 'Kulczynski', 'Hamann'
        )
        self.uiSimilMetricCombo.addItems(simil_metrics)
        # self.uiSimilMetricCombo.setCurrentText('Tanimoto')
        fp_types = project.fps[0].keys()
        self.uiFPTypeCombo.addItems(fp_types)
        # self.uiFPTypeCombo.setCurrentText('Morgan2')
        self.uiShowSimilMapButton.clicked.connect(self.show_simil_map)
        self.clusters = None

    def set_model(self):
        self.simil_df = self.similarity.similarity_matrix(
            self.fp_type, self.simil_metric
        )
        new_simil_df = self.get_HCL()
        simil_model = SimilarityModel(new_simil_df)
        self.uiSimilarityTable.setModel(simil_model)
        self.uiSimilarityTable.resizeColumnsToContents()
        self.uiSimilarityTable.resizeRowsToContents()

    def set_simil_metric(self, simil_metric):
        self.simil_metric = simil_metric
        self.set_model()

    def set_fp_type(self, fp_type):
        self.fp_type = fp_type
        if fp_type in (
            'Morgan1', 'Morgan2', 'Morgan3',
            'Hasehd Topological Torsions', 'RDKit',
            'Hashed Atom Pairs (chiral)',
            'Hashed Atom Pairs (achiral)'
        ):
            self.uiShowSimilMapButton.setEnabled(True)
        else:
            self.uiShowSimilMapButton.setEnabled(False)
        self.set_model()

    def show_heatmap(self):
        new_simil_df = self.get_HCL()
        self.heatmap = Heatmap(new_simil_df, self.fp_type, self.simil_metric)
        self.heatmap.plot()
        self.heatmap.show()

    def get_HCL(self, clus_num=0):
        if self.simil_df.shape[0] < 4:
            return
        simil_matrix = self.simil_df.to_numpy()
        linked = linkage(simil_matrix, 'single')
        n_cluster_to_test = range(2, 11)
        silhouette_avg = 0
        best = clus_num
        if not clus_num:
            for n_clusters in n_cluster_to_test:
                cluster_labels = fcluster(linked, n_clusters, criterion='maxclust')
                score = silhouette_score(
                    simil_matrix, cluster_labels, metric='euclidean'
                )
                if n_clusters > 2 and score > silhouette_avg:
                    best = n_clusters
                    silhouette_avg = score
        for i in np.linspace(0.10, 10, 10000):
            n_clusters = max(fcluster(
                linked, i, criterion='distance'
            ))
            if n_clusters == best:
                threshold = i
                break
            else:
                threshold = 0.01
        self.canvas.ax.clear()
        sorted_indexes = leaves_list(linked)
        new_df = self.simil_df.iloc[sorted_indexes, sorted_indexes] # ojo aquí
        names = self.simil_df.columns.tolist()
        dend = dendrogram(
            linked, orientation='left', labels=names,
            distance_sort='descending', show_leaf_counts=True,
            ax=self.canvas.ax, get_leaves=True, color_threshold=threshold
        )
        self.canvas.ax.axvline(threshold, linestyle='--', linewidth=.5, c='red')
        cluster_labels = fcluster(linked, best, criterion='maxclust')
        silhouette_avg = silhouette_score(
            simil_matrix, cluster_labels, metric='euclidean'
        )
        self.clusters = {
            'type': 'HCL',
            'labels': cluster_labels,
            'sorted_indexes': sorted_indexes,
            'total': best,
            'silhouette_avg': round(silhouette_avg, 3),
            'threshold': threshold
        }
        self.uiSilhouetteLabel.setText(
            f'Silhouette score (avg): {silhouette_avg:.3f}'
        )
        self.uiTLabel.setText(f't: {round(threshold, 2)}')
        self.uiNumClustersSpin.setValue(best)
        self.canvas.ax.set_title(
            f'Similarity HCL cluster ({self.fp_type}; {self.simil_metric})'
        )
        self.canvas.ax.spines['left'].set_visible(False)
        self.canvas.ax.spines['top'].set_visible(False)
        self.canvas.ax.spines['right'].set_visible(False)
        self.canvas.draw()
        return new_df

    def save_clusters(self):
        self.project.set_clusters(self.clusters)
        self.parent.database.commit()

    def show_simil_map(self):
        indexes = self.uiSimilarityTable.selectedIndexes()
        if not indexes:
            return
        i = indexes[0].row()
        j = indexes[0].column()
        ref = self.molecules[i].get_2d_mol
        ref_name = self.molecules[i].get_name
        prob = self.molecules[j].get_2d_mol
        prob_name = self.molecules[j].get_name
        setattr(
            self, f'simil_map_{self.fp_type}_{i}_{j}',
            SimilMap(ref, prob, self.fp_type, ref_name, prob_name)
        )


class ChEMBLExplorer(qtw.QWidget):

    def __init__(self, parent):
        super().__init__()
        uic.loadUi('Views/uiChembl.ui', self)
        self.parent = parent
        self.my_chembl = MyChembl()
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        self.available_targets = pd.DataFrame()
        self.std_type = pd.DataFrame()
        self.std_relation = pd.DataFrame()
        self.std_units = pd.DataFrame()
        self.activities = pd.DataFrame()
        self.filtered = self.activities.copy()
        self.similarities = []
        self.smiles = []
        self.target = ''
        self.uiProgressBar.setVisible(False)
        self.current_chunk = 1
        self.pending = 0
        self.mol_ref = None
        self.mol_ref_in_db = False
        self.fp_method = ''
        fp_methods = (
            'Morgan2', 'Morgan3', 'Morgan1',
            'Hashed Atom Pairs (achiral)',
            'Hashed Atom Pairs (chiral)',
            'RDKit'
        )
        self.uiFPMethodCombo.addItems(fp_methods)

    def search_target(self):
        target = self.uiTargetLine.text()
        if not target:
            return False
        self.available_targets = self.my_chembl.get_target_id(target)
        target_model = PandasModel(self.available_targets)
        self.uiTargetsTable.setModel(target_model)
        self.uiTargetsTable.resizeColumnsToContents()
        self.uiTargetsTable.resizeRowsToContents()

    @qtc.pyqtSlot(object, float)
    def query_workflow(self, df, process):
        if df is not False:
            self.activities = pd.concat(
                [self.activities, df], ignore_index=True
            )
            self.uiProgressBar.setVisible(True)
            self.uiProgressBar.setValue(process)
            self.pending = 100 - process
            if not self.pending:
                self.activities.loc[:, 'value'] = pd.to_numeric(
                    self.activities['value'], downcast='float'
                )
                self.homogenize_units()
                self.uiProgressBar.setVisible(False)
            self.display_available_types()
        else:
            self.pending = 0
            print('An error has occurred')

    @qtc.pyqtSlot(object, float)
    def simil_workflow(self, result, process):
        if result is not False:
            self.similarities.append(result.get('similarity'))
            self.smiles.append(result.get('smiles'))
            if len(self.similarities) == self.filtered.shape[0]:
                self.filtered.loc[:, 'similarity'] = self.similarities
                self.filtered.loc[:, 'smiles'] = self.smiles
                self.similarities = []
                self.smiles = []
            self.uiProgressBar.setVisible(True)
            self.uiProgressBar.setValue(process)
            self.pending = 100 - process
            if not self.pending:
                self.uiProgressBar.setVisible(False)
                self.uiFilterByMoleculeButton.setEnabled(True)
            self.plot_histogram()
            self.display_chunk()
        else:
            print('Hubo un error')  # intercambiar por un mensaje crítico

    def plot_histogram(self):
        if 'similarity' in self.filtered.columns.values:
            self.canvas.ax.clear()
            threshold = self.uiThresSpin.value() / 100
            self.canvas.ax.hist(self.filtered[f'similarity'], 10)
            self.canvas.ax.axvline(threshold, color='red', linestyle='--')
            self.canvas.ax.set_title('Distribution after filtering')
            self.canvas.ax.set_xlabel('Tanimoto similarity')
            self.canvas.ax.set_ylabel('Count')
            self.canvas.draw()
            self.uiFilterByMoleculeButton.setEnabled(True)
            self.current_chunk = 1

    # @logger.catch
    def select_target(self):
        index = self.uiTargetsTable.currentIndex()
        if index:
            self.target = self.available_targets.loc[index.row(), 'target_chembl_id']
            self.uiTargetLabel.setText(self.target)
            query = self.my_chembl.get_activities_for_a_target(self.target)
            if len(query) > 50:
                columns = [
                    'molecule_chembl_id', 'standard_type', 'standard_relation',
                    'standard_units', 'value'
                ]
                self.activities = pd.DataFrame(columns=columns)
                paginator = Paginator(query, 50)
                self.total_process = paginator.num_pages
                self.current_process = 0
                pages = [paginator.page(page_num) for page_num in paginator.page_range]
                self.send_to_worker(
                    function=self.my_chembl.get_query_df,
                    workflow=self.query_workflow,
                    sequence=pages
                )
            else:
                self.activities = pd.DataFrame.from_records(query)
                self.display_available_types()

    def send_to_worker(self, function, workflow, sequence=[], kwargs={}, mapping=False):
        self.worker = GenericWorker(
            function, sequence, **kwargs, mapping=mapping
        )
        self.thread = qtc.QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.workflow.connect(workflow)
        self.uiCancelProcessButton.clicked.connect(self.worker.pause)
        self.thread.start()

    def display_available_types(self):
        self.available_types = pd.unique(self.activities['standard_type'])
        items = pd.Series({
            i: sum(self.activities['standard_type'] == i) \
            for i in self.available_types
        }, dtype=str)
        std_types = pd.DataFrame(items, columns=['Available'])
        std_type_model = PandasModel(std_types)
        self.uiStdTypeTable.setModel(std_type_model)
        self.uiStdTypeTable.resizeColumnsToContents()
        self.uiStdTypeTable.resizeRowsToContents()
        self.total = self.activities.shape[0]
        self.uiTotalLabel.setText(str(self.total))

    def cancel_process(self):
        self.pending = 0
        self.uiProgressBar.setVisible(False)

    def homogenize_units(self):
        mask = self.activities['standard_units'] == 'nM'
        nM_to_uM = self.activities[mask].copy()
        nM_to_uM.loc[:, 'value'] = nM_to_uM.loc[:, 'value'] / 1000
        nM_to_uM.loc[:, 'standard_units'] = 'uM'
        self.activities.update(nM_to_uM)
        mask = self.activities['standard_units'] == 'mM'
        mM_to_uM = self.activities[mask].copy()
        mM_to_uM.loc[:, 'value'] = mM_to_uM.loc[:, 'value'] * 1000
        mM_to_uM.loc[:, 'standard_units'] = 'uM'
        self.activities.update(mM_to_uM)

    def filter_by_type(self):
        indexes = self.uiStdTypeTable.selectedIndexes()
        if indexes and not self.pending:
            selected_items = [self.available_types[i.row()] for i in indexes]
            filters = [
                self.activities['standard_type'] == i for i in selected_items
            ]
            self.type_filter = reduce(lambda x, y: x + y, filters)
            self.filtered = self.activities.loc[self.type_filter].copy()
            self.filtered.reset_index(inplace=True, drop=True)
            self.available_units = pd.unique(self.filtered['standard_units'])
            items = pd.Series({
                unit: sum(self.filtered['standard_units'] == unit) \
                for unit in self.available_units
            }, dtype=str)
            std_units = pd.DataFrame(items, columns=['Amounts'])
            std_unit_model = PandasModel(std_units)
            self.uiStdUnitsTable.setModel(std_unit_model)
            self.uiStdUnitsTable.resizeColumnsToContents()
            self.uiStdUnitsTable.resizeRowsToContents()
            self.current_chunk = 1
            self.display_chunk()

    def filter_by_unit(self):
        indexes = self.uiStdUnitsTable.selectedIndexes()
        if indexes and not self.pending:
            selected_items = [self.available_units[i.row()] for i in indexes]
            filters = [
                self.activities['standard_units'] == i for i in selected_items
            ]
            self.unit_filter = reduce(lambda x, y: x + y, filters)
            self.type_unit_filter = reduce(
                lambda x, y: x * y, [self.type_filter, self.unit_filter]
            )
            self.filtered = self.activities.loc[self.type_unit_filter].copy()
            self.filtered.reset_index(inplace=True, drop=True)
            self.available_relations = pd.unique(self.filtered['standard_relation'])
            items = pd.Series({
                relation: sum(self.filtered['standard_relation'] == relation) \
                for relation in self.available_relations
            }, dtype=str)
            std_relations = pd.DataFrame(items, columns=['Amount'])
            std_relation_model = PandasModel(std_relations)
            self.uiStdRelationTable.setModel(std_relation_model)
            self.uiStdRelationTable.resizeColumnsToContents()
            self.uiStdRelationTable.resizeRowsToContents()
            self.current_chunk = 1
            self.display_chunk()

    def filter_by_relation(self):
        indexes = self.uiStdRelationTable.selectedIndexes()
        if indexes and not self.pending:
            selected_items = [self.available_relations[i.row()] for i in indexes]
            filters = [
                self.activities['standard_relation'] == i for i in selected_items
            ]
            self.relation_filter = reduce(lambda x, y: x + y, filters)
            self.type_unit_relation_filter = reduce(
                lambda x, y: x * y, [
                    self.type_filter, self.unit_filter, self.relation_filter
                ]
            )
            self.filtered = self.activities.loc[self.type_unit_relation_filter].copy()
            self.filtered.reset_index(inplace=True, drop=True)
            self.current_chunk = 1
            self.display_chunk()
            self.uiUploadMoleculeButton.setEnabled(True)

    def filter_by_molecule(self):
        self.uiFilterByMoleculeButton.setEnabled(False)
        mol_ref = self.mol_ref.get_2d_mol
        kwargs = {
            'fp_method': self.fp_method,
            'mol_ref': mol_ref
        }
        self.send_to_worker(
            function=self.my_chembl.get_similarities,
            workflow=self.simil_workflow,
            sequence=self.filtered['molecule_chembl_id'],
            kwargs=kwargs
        )

    def display_chunk(self):
        self.total = self.filtered.shape[0]
        self.uiTotalLabel.setText(str(self.total))
        self.chunks = self.total // 100 + 1
        if self.chunks == 1 or self.current_chunk == self.chunks:
            self.uiNextChunkButton.setEnabled(False)
        elif self.current_chunk == 1:
            self.uiPrevChunkButton.setEnabled(False)
            self.uiNextChunkButton.setEnabled(True)
        else:
            self.uiNextChunkButton.setEnabled(True)
        self.uiPageLabel.setText(f'Page: {self.current_chunk}/{self.chunks}')
        filtered_model = PandasModel(self.filtered, chunk=self.current_chunk)
        self.uiFilteredTable.setModel(filtered_model)
        self.uiFilteredTable.resizeRowsToContents()
        self.uiFilteredTable.resizeColumnsToContents()

    def display_selected_mol(self, index):
        has_smiles = 'smiles' in self.filtered.columns.values.tolist()
        if index and has_smiles:
            smiles = self.filtered.loc[index.row(), 'smiles']
            mol_prob = Chem.MolFromSmiles(smiles)
            mol_ref = self.mol_ref.get_2d_mol
            self.img = SimilMap.get_simil_map(
                mol_ref, mol_prob, self.fp_method, size=(350, 350)
            )
            qimage = qtg.QPixmap.fromImage(ImageQt.ImageQt(self.img))
            self.uiSelectedMoleculeLabel.setPixmap(qimage)
            self.uiSaveImageButton.setEnabled(True)
        else:
            self.uiSaveImageButton.setEnabled(False)

    def save_img(self):
        filename, extension = qtw.QFileDialog.getSaveFileName(
            self, 'Save image', 'Image.png', 'PNG (*.png)',
            options=qtw.QFileDialog.DontUseNativeDialog,
        )
        if filename:
            has_extension = '.' in filename
            if not has_extension:
                filename += '.png'
            self.img.save(filename)

    def next_chunk(self):
        self.current_chunk += 1
        self.uiPrevChunkButton.setEnabled(True)
        self.uiPageLabel.setText(f'Page: {self.current_chunk}/{self.chunks}')
        self.display_chunk()

    def prev_chunk(self):
        self.current_chunk -= 1
        self.uiPageLabel.setText(f'Page: {self.current_chunk}/{self.chunks}')
        self.uiNextChunkButton.setEnabled(True)
        # if self.current_chunk == 1:
            # self.uiPrevChunkButton.setEnabled(False)
        self.display_chunk()

    def upload_filter_molecule(self):
        mol_path, _ = qtw.QFileDialog.getOpenFileName(
            self, 'Select the molecule',
            filter='*.mol;;*.mol2;;*.pdb;;*.txt;;*.smi;;*.log;;*.csv',
            options=qtw.QFileDialog.DontUseNativeDialog,
        )
        if mol_path:
            file_format = mol_path.split('.')[-1]
            mol = Molecule(mol_path, file_format=file_format)
            self.mol_ref_in_db = self.parent.database.check(
                'molecules', mol.smiles
            )
            if self.mol_ref_in_db:
                self.mol_ref = self.parent.database.get(
                    'molecules', mol.smiles
                )
            else:
                self.mol_ref = mol
            self.uiMoleculeLabel.setEnabled(True)
            self.uiMoleculeLabel.setPixmap(
                qtg.QPixmap(self.mol_ref.mol_img)
            )
            self.uiThresSpin.setEnabled(True)
            self.uiFPMethodCombo.setEnabled(True)
            self.uiRemoveMoleculeButton.setEnabled(True)
            self.uiFilterByMoleculeButton.setEnabled(True)

    def remove_filter_molecule(self):
        self.mol_ref = None
        self.uiMoleculeLabel.setEnabled(False)
        self.uiRemoveMoleculeButton.setEnabled(False)
        self.uiThresSpin.setEnabled(False)
        self.uiFPMethodCombo.setEnabled(False)
        self.uiFilterByMoleculeButton.setEnabled(False)

    def drop_molecule(self):
        index = self.uiFilteredTable.currentIndex()
        if index:
            self.filtered.drop(index=index.row(), inplace=True)
            self.filtered.reset_index(inplace=True, drop=True)
            self.display_chunk()

    def set_fp_method(self, value):
        self.fp_method = value
        if value and self.mol_ref is not None:
            self.uiFilterByMoleculeButton.setEnabled(True)
        else:
            self.uiFilterByMoleculeButton.setEnabled(False)

    def send_molecules(self):
        pass


# TODO
# 1. preguntar si quiero conservar la molécula de referencia creada,
# siempre y cuando no estuviera previamente en la BD
# 3. 
