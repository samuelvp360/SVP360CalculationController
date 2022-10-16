#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from PyQt5 import QtWidgets as qtw
from PyQt5 import uic
from Models import PandasModel, DisplayMFPModel, SimilarityModel
from Components import Similarity
from Plotters import PlotCanvas, Heatmap, SimilMap
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.metrics import silhouette_samples, silhouette_score
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar


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

    def __init__(self, project):
        super().__init__()
        uic.loadUi('Views/uiSimilarity.ui', self)
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
        self.uiShowSimilMapButton.clicked.connect(
            lambda: self.show_simil_map(self.fp_type)
        )

    def set_model(self):
        self.simil_df = self.similarity.similarity_matrix(
            self.fp_type, self.simil_metric
        )
        self.HCL = self.get_HCL()
        simil_model = SimilarityModel(self.simil_df)
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
        self.heatmap = Heatmap(self.simil_df, self.fp_type, self.simil_metric)
        self.heatmap.plot()
        self.heatmap.show()

    def get_HCL(self, clus_num=0):
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
        dend = dendrogram(
            linked, orientation='left', labels=self.names,
            distance_sort='descending', show_leaf_counts=True,
            ax=self.canvas.ax, get_leaves=True, color_threshold=threshold
        )
        self.canvas.ax.axvline(threshold, linestyle='--', linewidth=.5, c='red')
        cluster_labels = fcluster(linked, best, criterion='maxclust')
        silhouette_avg = silhouette_score(
            simil_matrix, cluster_labels, metric='euclidean'
        )
        # self.project.set_clusters(
            # 'HCL', cluster_labels, best,
            # round(silhouette_avg, 3), threshold
        # )
        # for i in range(1, best + 1):
            # print(
                # [m.get_name for m in self.project.get_cluster(i)]
            # )
        self.uiSilhouetteLabel.setText(
            f'Silhouette score (avg): {silhouette_avg}'
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

    def show_simil_map(self, fp_type):
        indexes = self.uiSimilarityTable.selectedIndexes()
        if not indexes:
            return
        i = indexes[0].row()
        j = indexes[0].column()
        ref = self.molecules[i].mol
        ref_name = self.molecules[i].get_name
        ref = Chem.MolFromSmiles(Chem.MolToSmiles(ref))
        prob = self.molecules[j].mol
        prob_name = self.molecules[j].get_name
        prob = Chem.MolFromSmiles(Chem.MolToSmiles(prob))
        d = Draw.MolDraw2DCairo(550, 550)
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
        setattr(
            self, f'simil_map_{fp_type}_{i}_{j}',
            SimilMap(data, ref_name, prob_name, fp_type)
        )


