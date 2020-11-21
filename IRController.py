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
        self._bandThreshold = self.uiBandsThresholdSpinBox.value()
        self._orderBase = self.uiOrderBaseSpinBox.value()
        self._orderBands = self.uiOrderBandsSpinBox.value()
        self._showGrid = self.uiShowGridButton.isChecked()
        self._baseline = False
        self._bandFitting = False
        self._bandDict = pd.DataFrame({
            'Wavenumber': [],
            'Transmittance': [],
            '%Transmittance': [],
            'Absorbance': [],
            'Height (T)': [],
            'Height (%T)': [],
            'Height (A)': []
        })
        self._predictedBands = {}
        self.canvas = IRCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiSpectraLayout.addWidget(self.canvas)
        self.cid = self.canvas.fig.canvas.mpl_connect('button_press_event', self.BandsDetectionManual)
        self.SetYAxis()
        self.SelectYAxis('T')
# ------------------------------------SIGNALS---------------------------------------------------------------
        self.uiDetectBaselineButton.clicked.connect(self.BaselineDetection)
        self.uiSubsBaselineButton.clicked.connect(self.BaselineSubs)
        self.uiBaseThresholdSpinBox.valueChanged.connect(self.SetBaseThreshold)
        self.uiBandsThresholdSpinBox.valueChanged.connect(self.SetBandsThreshold)
        self.uiOrderBaseSpinBox.valueChanged.connect(self.SetOrderBase)
        self.uiOrderBandsSpinBox.valueChanged.connect(self.SetOrderBands)
        self.uiShowGridButton.clicked.connect(self.SetGrid)
        self.uiTransButton.clicked.connect(lambda: self.SelectYAxis('T'))
        self.uiPercentTransButton.clicked.connect(lambda: self.SelectYAxis('%T'))
        self.uiAbsorbanceButton.clicked.connect(lambda: self.SelectYAxis('A'))
        self.uiRestoreButton.clicked.connect(self.Restore)
        self.uiBandAutoButton.clicked.connect(self.BandsDetectionAuto)
        self.uiBandFittingAutoButton.clicked.connect(self.BandFitting)
# ------------------------------------METHODS---------------------------------------------------------------

    def Plot(self, yAxis):
        """
        docstring
        """
        self.canvas.axes.clear()
        self.canvas.axes.plot(self.df['Wavenumber'], self.df[yAxis], linewidth=1.0)
        self.canvas.axes.set_title(self._title)
        self.canvas.axes.set_xlabel('Wavenumber [cm-1]')
        self.canvas.axes.set_ylabel(yAxis)
        self.canvas.axes.set_xlim(4000, 500)
        self.canvas.axes.grid(True, color='gray', linestyle=':', linewidth=0.5) if self._showGrid else self.canvas.axes.grid(False)
        if self._baseline:
            self.canvas.axes.plot(self.df['Wavenumber'], self.df['baseline'], linewidth=1.0, color='red')
        if self._bandFitting:
            for _, predicted in self._predictedBands[yAxis].iteritems():
                self.canvas.axes.plot(self.df['Wavenumber'], predicted, linewidth=0.0)
                if yAxis == 'Transmittance':
                    self.canvas.axes.plot(
                        self.df['Wavenumber'], self._predictedBands['Transmittance']['Fitted Curve Lorentzian'],
                        color='red', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.axes.fill_between(self.df['Wavenumber'], predicted, 1, alpha=0.2)
                elif yAxis == '%Transmittance':
                    self.canvas.axes.plot(
                        self.df['Wavenumber'], self._predictedBands['%Transmittance']['Fitted Curve Lorentzian'],
                        color='red', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.axes.fill_between(self.df['Wavenumber'], predicted, 100, alpha=0.2)
                else:
                    self.canvas.axes.plot(
                        self.df['Wavenumber'], self._predictedBands['Absorbance']['Fitted Curve Lorentzian'],
                        color='red', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.axes.fill_between(self.df['Wavenumber'], predicted, alpha=0.2)
        if self._bandDict.shape[0] > 0:
            for index, row in self._bandDict.iterrows():
                self.canvas.axes.annotate(
                    str(round(row['Wavenumber'], 0)), xy=(row['Wavenumber'], row[yAxis]), xytext=(0, -50),
                    textcoords='offset points', ha='center', va='bottom', rotation=90,
                    bbox=dict(boxstyle='round', fc='white', alpha=0.5),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                )
        self.canvas.draw()

    def SetBaseThreshold(self):
        self._baseThreshold = self.uiBaseThresholdSpinBox.value()
        if self.df.get('baseline') is not None:
            self.BaselineDetection()

    def SetBandsThreshold(self):
        self._bandThreshold = self.uiBandsThresholdSpinBox.value()
        if self.df.get('Bands') is not None:
            del self.df['Bands']
            self.BandsDetectionAuto()

    def SetOrderBase(self):
        self._orderBase = self.uiOrderBaseSpinBox.value()
        if self.df.get('baseline') is not None:
            self.BaselineDetection()

    def SetOrderBands(self):
        self._orderBands = self.uiOrderBandsSpinBox.value()
        if self.df.get('Bands') is not None:
            self.BandsDetectionAuto()

    def SetGrid(self):
        self._showGrid = self.uiShowGridButton.isChecked()
        if self.uiTransButton.isChecked():
            self.SelectYAxis('T')
        elif self.uiPercentTransButton.isChecked():
            self.SelectYAxis('%T')
        else:
            self.SelectYAxis('A')

    def SetYAxis(self):
        """
        docstring
        """
        if 0.0 < self.df['Raw Data'].median() < 1.0:
            self.df['Transmittance'] = self.df['Raw Data']
        self.df['%Transmittance'] = round(self.df['Transmittance'] * 100, 6)
        self.df['Absorbance'] = pd.Series([round(2.0 - math.log10(i), 6) for i in self.df['%Transmittance']])

    def SelectYAxis(self, unit):

        if unit == 'T':
            self.uiTransButton.setChecked(True)
            self.uiPercentTransButton.setChecked(False)
            self.uiAbsorbanceButton.setChecked(False)
            self.Plot('Transmittance')
        elif unit == '%T':
            self.uiPercentTransButton.setChecked(True)
            self.uiAbsorbanceButton.setChecked(False)
            self.uiTransButton.setChecked(False)
            self.Plot('%Transmittance')
        elif unit == 'A':
            self.uiAbsorbanceButton.setChecked(True)
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
        self._baseline = True
        self.uiSubsBaselineButton.setEnabled(True)
        self.uiDetectBaselineButton.setEnabled(False)
        if self.uiTransButton.isChecked():
            self.SelectYAxis('T')
        elif self.uiPercentTransButton.isChecked():
            self.SelectYAxis('%T')
        else:
            self.SelectYAxis('A')

    def BaselineSubs(self):
        """
        docstring
        """
        self.df['Transmittance'] = 1 - self.df['baseline'] + self.df['Transmittance']
        self.df['%Transmittance'] = self.df['Transmittance'] * 100
        self.df['Absorbance'] = pd.Series([2.0 - math.log10(i) for i in self.df['%Transmittance']])
        self.uiDetectBaselineButton.setEnabled(True)
        self.uiSubsBaselineButton.setEnabled(False)
        self.uiBaseThresholdSpinBox.setEnabled(False)
        self.uiOrderBaseSpinBox.setEnabled(False)
        self._baseline = False
        if self._bandDict.shape[0] > 0:
            self.BandsDetectionAuto()
        else:
            if self.uiTransButton.isChecked():
                self.SelectYAxis('T')
            elif self.uiPercentTransButton.isChecked():
                self.SelectYAxis('%T')
            else:
                self.SelectYAxis('A')

    def Restore(self):
        """
        docstring
        """
        self.uiDetectBaselineButton.setEnabled(True)
        self.uiBaseThresholdSpinBox.setEnabled(True)
        self.uiOrderBaseSpinBox.setEnabled(True)
        self.uiSubsBaselineButton.setEnabled(False)
        self.uiBandManualButton.setChecked(False)
        self.uiBandFittingAutoButton.setEnabled(False)
        self._bandFitting = False
        if self.df.get('Bands') is not None:
            del self.df['Bands']
        self._bandDict = pd.DataFrame({
            'Wavenumber': [],
            'Transmittance': [],
            '%Transmittance': [],
            'Absorbance': [],
            'Height (T)': [],
            'Height (%T)': [],
            'Height (A)': []
        })
        self._predictedBands = {}
        self.SetYAxis()
        self.SelectYAxis('T')

    def BandsDetectionAuto(self):
        """
        docstring
        """
        self.df['Bands'] = self.df.loc[argrelextrema(self.df.Transmittance.values, np.less_equal, order=self._orderBands)[0], 'Transmittance']
        mask3 = (self.df['Bands'] > self._bandThreshold)
        self.df.loc[mask3, 'Bands'] = np.nan
        mask4 = (self.df['Bands'].isnull() == False)
        self._bandDict = self.df.loc[mask4, ['Wavenumber', 'Transmittance', '%Transmittance', 'Absorbance']].copy()
        self._bandDict = self._bandDict.reset_index(drop=True)
        self._bandDict['Height (T)'] = pd.Series([1 - i for i in self._bandDict['Transmittance']])
        self._bandDict['Height (%T)'] = pd.Series([100 - i for i in self._bandDict['%Transmittance']])
        self._bandDict['Height (A)'] = self._bandDict['Absorbance']
        self._bandDict['HWHM'] = pd.Series(map(self.HWHM, self._bandDict['Wavenumber'], self._bandDict['Absorbance']))
        self.uiBandFittingAutoButton.setEnabled(True)
        if self.uiTransButton.isChecked():
            self.SelectYAxis('T')
        elif self.uiPercentTransButton.isChecked():
            self.SelectYAxis('%T')
        else:
            self.SelectYAxis('A')

    def BandsDetectionManual(self, event):
        """
        docstring
        """
        if self.uiBandManualButton.isChecked():
            ix, iy = event.xdata, event.ydata
            if self.uiTransButton.isChecked():
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [iy],
                    '%Transmittance': [iy * 100],
                    'Absorbance': [2.0 - math.log10(iy * 100)],
                    'Height (T)': [1 - iy],
                    'Height (%T)': [100 - (iy * 100)],
                    'Height (A)': [2.0 - math.log10(iy * 100)],
                    'HWHM': [self.HWHM(ix, 2.0 - math.log10(iy * 100))]
                })
                self._bandDict = self._bandDict.append(newBand, ignore_index=True)
            elif self.uiPercentTransButton.isChecked():
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [iy / 100],
                    '%Transmittance': [iy],
                    'Absorbance': [2.0 - math.log10(iy)],
                    'Height (T)': [1 - (iy / 100)],
                    'Height (%T)': [100 - iy],
                    'Height (A)': [2.0 - math.log10(iy)],
                    'HWHM': [self.HWHM(ix, 2.0 - math.log10(iy))]
                })
                self._bandDict = self._bandDict.append(newBand, ignore_index=True)
            elif self.uiAbsorbanceButton.isChecked():
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [(10 ** (2 - iy)) / 100],
                    '%Transmittance': [10 ** (2 - iy)],
                    'Absorbance': [iy],
                    'Height (T)': [1 - (10 ** (2 - iy)) / 100],
                    'Height (%T)': [100 - (10 ** (2 - iy))],
                    'Height (A)': [iy],
                    'HWHM': [self.HWHM(ix, iy)]
                })
                self._bandDict = self._bandDict.append(newBand, ignore_index=True)
            self.uiBandFittingAutoButton.setEnabled(True)
            if self._bandFitting:
                self.BandFitting()
            else:
                if self.uiTransButton.isChecked():
                    self.SelectYAxis('T')
                elif self.uiPercentTransButton.isChecked():
                    self.SelectYAxis('%T')
                else:
                    self.SelectYAxis('A')

    def HWHM(self, x, y):
        '''
        docstring
        '''
        halfMaximum = y / 2.0
        candidates = [round(abs(j - x), 6) for i, j in zip(self.df['Absorbance'], self.df['Wavenumber']) if abs(i - halfMaximum) < 1e-3]
        hwhm = sorted(candidates)[0]
        if hwhm == 0.000000:
            hwhm = 1.0
        return hwhm

    def BandFitting(self):
        """
        docstring
        """
        self._predictedBands['Absorbance'] = pd.DataFrame(
            map(self.LorentzianFunction, self._bandDict['Wavenumber'], self._bandDict['Height (A)'], self._bandDict['HWHM'])
        ).transpose()
        self._predictedBands['Transmittance'] = self._predictedBands['Absorbance'].applymap(lambda x: (10 ** (2 - x)) / 100)
        self._predictedBands['%Transmittance'] = self._predictedBands['Transmittance'].applymap(lambda x: 100 * x)
        self._predictedBands['Absorbance']['Fitted Curve Lorentzian'] = self._predictedBands['Absorbance'].sum(axis=1)
        print(self._predictedBands['Absorbance'])
        self._predictedBands['Transmittance']['Fitted Curve Lorentzian'] = self._predictedBands['Absorbance']['Fitted Curve Lorentzian'].apply(lambda x: (10 ** (2 - x)) / 100)
        print(self._predictedBands['Transmittance'])
        self._predictedBands['%Transmittance']['Fitted Curve Lorentzian'] = self._predictedBands['Transmittance']['Fitted Curve Lorentzian'].apply(lambda x: 100 * x)
        print(self._predictedBands['%Transmittance'])
        self._bandFitting = True
        self.uiBandFittingAutoButton.setEnabled(False)
        if self.uiTransButton.isChecked():
            self.SelectYAxis('T')
        elif self.uiPercentTransButton.isChecked():
            self.SelectYAxis('%T')
        else:
            self.SelectYAxis('A')

    def LorentzianFunction(self, x0, amp, hwhm):
        """
        docstring
        """
        xRange = self.df['Wavenumber']
        return pd.Series([amp / (1 + ((x - x0) / hwhm) ** 2) for x in xRange])
