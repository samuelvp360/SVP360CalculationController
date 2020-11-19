#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import uic
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from scipy.signal import argrelextrema
import math

matplotlib.use('Qt5Agg')


class IRCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(11, 4), dpi=100, facecolor='#f8b4a5', tight_layout=True)
        self.axes = self.fig.add_subplot(111)
        super(IRCanvas, self).__init__(self.fig)


class IRPlotter(qtw.QWidget):
    """
    docstring
    """
    def __init__(self, csvData1):
        super(IRPlotter, self).__init__()
        uic.loadUi('Views/uiIRPlotter.ui', self)
        self.df = pd.read_csv(csvData1, names=('Wavenumber', 'Raw Data'))
        self._title = csvData1.replace('.csv', '')
        self._baseThreshold = self.uiBaseThresholdSpinBox.value()
        self._peaksThreshold = self.uiPeaksThresholdSpinBox.value()
        self._orderBase = self.uiOrderBaseSpinBox.value()
        self._orderPeaks = self.uiOrderPeaksSpinBox.value()
        self._showGrid = self.uiShowGridButton.isChecked()
        self._peaksList = []
        self.canvas = IRCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiSpectraLayout.addWidget(self.canvas)
        self.cid = self.canvas.fig.canvas.mpl_connect('button_press_event', self.PeaksDetectionManual)
        self.SetYAxis()
        self.SelectYAxis()
# ------------------------------------SIGNALS---------------------------------------------------------------
        self.uiDetectBaselineButton.clicked.connect(self.BaselineDetection)
        self.uiSubsBaselineButton.clicked.connect(self.BaselineSubs)
        self.uiBaseThresholdSpinBox.valueChanged.connect(self.SetBaseThreshold)
        self.uiPeaksThresholdSpinBox.valueChanged.connect(self.SetPeaksThreshold)
        self.uiOrderBaseSpinBox.valueChanged.connect(self.SetOrderBase)
        self.uiOrderPeaksSpinBox.valueChanged.connect(self.SetOrderPeaks)
        self.uiShowGridButton.clicked.connect(self.SetGrid)
        self.uiTransButton.clicked.connect(self.SelectYAxis)
        self.uiPercentTransButton.clicked.connect(self.SelectYAxis)
        self.uiAbsorbanceButton.clicked.connect(self.SelectYAxis)
        self.uiRestoreButton.clicked.connect(self.Restore)
        self.uiPeakAutoButton.clicked.connect(self.PeaksDetectionAuto)
# ------------------------------------METHODS---------------------------------------------------------------

    def Plot(self, yList):
        """
        docstring
        """
        self.canvas.axes.clear()
        self.canvas.axes.plot(self.df['Wavenumber'], self.df[yList], linewidth=1.0)
        self.canvas.axes.set_title(self._title)
        self.canvas.axes.set_xlabel('Wavenumber [cm-1]')
        self.canvas.axes.set_ylabel(yList[0]) if type(yList) == list else self.canvas.axes.set_ylabel(yList)
        self.canvas.axes.set_xlim(4000, 500)
        self.canvas.axes.grid(True, color='gray', linestyle=':', linewidth=0.5) if self._showGrid else self.canvas.axes.grid(False)
        if len(self._peaksList) > 0:
            if yList == 'Transmittance':
                intensity = 1
            elif yList == '%Transmittance':
                intensity = 2
            elif yList == 'Absorbance':
                intensity = 3
            else:  # quitar
                intensity = 1
            for label in self._peaksList:
                self.canvas.axes.annotate(
                    str(round(label[0], 0)), xy=(label[0], label[intensity]), xytext=(0, -50),
                    textcoords='offset points', ha='center', va='bottom', rotation=90,
                    bbox=dict(boxstyle='round', fc='white', alpha=0.5),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                )
        self.canvas.draw()

    def SetBaseThreshold(self):
        self._baseThreshold = self.uiBaseThresholdSpinBox.value()
        if self.df.get('baseline') is not None:
            self.BaselineDetection()

    def SetPeaksThreshold(self):
        self._peaksThreshold = self.uiPeaksThresholdSpinBox.value()
        if self.df.get('peaks') is not None:
            del self.df['peaks']
            self.PeaksDetectionAuto()

    def SetOrderBase(self):
        self._orderBase = self.uiOrderBaseSpinBox.value()
        if self.df.get('baseline') is not None:
            self.BaselineDetection()

    def SetOrderPeaks(self):
        self._orderPeaks = self.uiOrderPeaksSpinBox.value()
        if self.df.get('peaks') is not None:
            self.PeaksDetectionAuto()

    def SetGrid(self):
        self._showGrid = self.uiShowGridButton.isChecked()
        self.SelectYAxis()

    def SetYAxis(self):
        """
        docstring
        """
        if 0.0 < self.df['Raw Data'].median() < 1.0:
            self.df['Transmittance'] = self.df['Raw Data']
        self.df['%Transmittance'] = self.df['Transmittance'] * 100
        self.df['Absorbance'] = pd.Series([2.0 - math.log10(i) for i in self.df['%Transmittance']])

    def SelectYAxis(self):
        if self.uiTransButton.isChecked():
            self.uiPercentTransButton.setChecked(False)
            self.uiAbsorbanceButton.setChecked(False)
            self.Plot('Transmittance')
        elif self.uiPercentTransButton.isChecked():
            self.uiAbsorbanceButton.setChecked(False)
            self.uiTransButton.setChecked(False)
            self.Plot('%Transmittance')
        elif self.uiAbsorbanceButton.isChecked():
            self.uiPercentTransButton.setChecked(False)
            self.uiTransButton.setChecked(False)
            self.Plot('Absorbance')

    def BaselineDetection(self):
        """
        docstring
        """
        self.df['baseline'] = self.df.loc[argrelextrema(self.df.Transmittance.values, np.greater_equal, order=self._orderBase)[0], 'Transmittance']
        mask1 = (self.df['baseline'] < self._baseThreshold)
        self.df.loc[mask1, 'baseline'] = np.nan

        mask2 = (self.df['baseline'].isnull() == False)
        baselineX = self.df.loc[mask2, 'Wavenumber']
        baselineY = self.df.loc[mask2, 'Transmittance']
        coef = np.polyfit(baselineX, baselineY, 1)
        self.baselineFit = np.poly1d(coef)
        self.df['baseline'] = pd.Series([self.baselineFit(i) for i in self.df['Wavenumber']])
        self.Plot(['Transmittance', 'baseline'])

    def BaselineSubs(self):
        """
        docstring
        """
        self.df['Transmittance'] = 1 - self.df['baseline'] + self.df['Transmittance']
        self.df['%Transmittance'] = self.df['Transmittance'] * 100
        self.df['Absorbance'] = pd.Series([2.0 - math.log10(i) for i in self.df['%Transmittance']])
        self.uiDetectBaselineButton.setEnabled(False)
        self.uiSubsBaselineButton.setEnabled(False)
        self.uiBaseThresholdSpinBox.setEnabled(False)
        self.uiOrderBaseSpinBox.setEnabled(False)
        self.SelectYAxis()

    def Restore(self):
        """
        docstring
        """
        self.uiDetectBaselineButton.setEnabled(True)
        self.uiSubsBaselineButton.setEnabled(True)
        self.uiBaseThresholdSpinBox.setEnabled(True)
        self.uiOrderBaseSpinBox.setEnabled(True)
        self.uiPeakManualButton.setChecked(False)
        if self.df.get('peaks') is not None:
            del self.df['peaks']
        self._peaksList = []
        self.SetYAxis()
        self.SelectYAxis()

    def PeaksDetectionAuto(self):
        """
        docstring
        """
        self.df['peaks'] = self.df.loc[argrelextrema(self.df.Transmittance.values, np.less_equal, order=self._orderPeaks)[0], 'Transmittance']
        mask3 = (self.df['peaks'] > self._peaksThreshold)
        self.df.loc[mask3, 'peaks'] = np.nan
        mask4 = (self.df['peaks'].isnull() == False)
        self._peaksList = [
            (i, j, k, l) for i, j, k, l in zip(
                self.df.loc[mask4, 'Wavenumber'], self.df.loc[mask4, 'Transmittance'],
                self.df.loc[mask4, '%Transmittance'], self.df.loc[mask4, 'Absorbance']
            )
        ]
        self.SelectYAxis()

    def PeaksDetectionManual(self, event):
        """
        docstring
        """
        if self.uiPeakManualButton.isChecked():
            ix, iy = event.xdata, event.ydata
            if self.uiTransButton.isChecked():
                self._peaksList.append((ix, iy, iy * 100, 2.0 - math.log10(iy * 100)))
            elif self.uiPercentTransButton.isChecked():
                self._peaksList.append((ix, iy / 100, iy, 2.0 - math.log10(iy)))
            elif self.uiAbsorbanceButton.isChecked():
                self._peaksList.append((ix, (10 ** (2 - iy)) / 100, 10 ** (2 - iy), iy))
        self.SelectYAxis()
