#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import uic
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from scipy.signal import argrelextrema
import math
from Models import IRDataModel
from Worker import IRWorkerThread

matplotlib.use('Qt5Agg')


class IRCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(22, 4), dpi=100, facecolor='#2d2a2e')
        self.grid = GridSpec(17, 1, left=0.1, bottom=0.15, right=0.94, top=0.94, wspace=0.3, hspace=0.3)
        self.axes = self.fig.add_subplot(self.grid[0:12, 0])
        self.axes2 = self.fig.add_subplot(self.grid[14:, 0])
        super(IRCanvas, self).__init__(self.fig)


class IRPlotter(qtw.QMainWindow):
    """
    docstring
    """
    def __init__(self, csvData1):
        super(IRPlotter, self).__init__()
        uic.loadUi('Views/uiIRPlotter.ui', self)
        self.df = pd.read_csv(csvData1, names=('Wavenumber', 'Raw Data'))
        self._title = csvData1.split('/')[-1].replace('.csv', '')
        self.statusBar().showMessage(f'Showing IR spectrum from {self._title}')
        self._baseThreshold = self.uiBaseThresholdSpinBox.value()
        self._bandThreshold = self.uiBandsThresholdSpinBox.value()
        self._orderBase = self.uiOrderBaseSpinBox.value()
        self._orderBands = self.uiOrderBandsSpinBox.value()
        self._showGrid = self.uiShowGridButton.isChecked()
        self._baseline = False
        self._bandFitting = False
        self._normalized = False
        self._bandsDict = pd.DataFrame({
            'Wavenumber': [],
            'Transmittance': [],
            '%Transmittance': [],
            'Absorbance': [],
            'Height (T)': [],
            'Height (%T)': [],
            'Height (A)': []
        })
        # self._bandsModel = IRDataModel(self._bandsDict)
        # self.uiIRBandsTableView.setModel(self._bandsModel)
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
        self.uiNormalizeButton.clicked.connect(self.SetNormalization)
# ------------------------------------METHODS---------------------------------------------------------------

    def Plot(self, yAxis):
        """
        docstring
        """
        self.canvas.axes.clear()
        self.canvas.axes2.clear()
        self.canvas.axes.plot(self.df['Wavenumber'], self.df[yAxis], linewidth=1.5, color='#272822')
        self.canvas.axes.set_title(self._title, color='#ae81ff')
        self.canvas.axes.set_xlabel('Wavenumber [cm-1]', color='#f92672')
        self.canvas.axes.set_ylabel(yAxis, color='#f92672')
        self.canvas.axes.tick_params(axis='x', colors='#66d9ef')
        self.canvas.axes.tick_params(axis='y', colors='#66d9ef')
        self.canvas.axes2.set_ylabel('Residuals', color='#f92672')
        self.canvas.axes2.tick_params(axis='x', colors='#66d9ef')
        self.canvas.axes2.tick_params(axis='y', colors='#66d9ef')
        self.canvas.axes2.axhline()
        self.canvas.axes.set_xlim(4000, 500)
        self.canvas.axes2.set_xlim(4000, 500)
        self.canvas.axes.grid(True, color='#2d2a2e', linestyle=':', linewidth=0.5) if self._showGrid else self.canvas.axes.grid(False)
        if self._baseline:
            self.canvas.axes.plot(self.df['Wavenumber'], self.df['baseline'], linewidth=1.0, color='#f92672')
        if self._bandFitting:
            print(self.df)
            print(self._bandsDict)
            print(self._predictedBands)
            for _, predicted in self._predictedBands[yAxis].iteritems():
                self.canvas.axes.plot(self.df['Wavenumber'], predicted, linewidth=0.0)
                self.canvas.axes2.plot(
                    self.df['Wavenumber'], self._predictedBands['Residuals'], 'o',
                    markersize=0.2, linewidth=0.1, color='#ae81ff'
                )
                if yAxis == 'Transmittance':
                    self.canvas.axes.plot(
                        self.df['Wavenumber'], self._predictedBands['Transmittance']['Fitted Curve Lorentzian'],
                        color='#ae81ff', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.axes.fill_between(self.df['Wavenumber'], predicted, 1, alpha=0.2)
                elif yAxis == '%Transmittance':
                    self.canvas.axes.plot(
                        self.df['Wavenumber'], self._predictedBands['%Transmittance']['Fitted Curve Lorentzian'],
                        color='#ae81ff', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.axes.fill_between(self.df['Wavenumber'], predicted, 100, alpha=0.2)
                else:
                    self.canvas.axes.plot(
                        self.df['Wavenumber'], self._predictedBands['Absorbance']['Fitted Curve Lorentzian'],
                        color='#ae81ff', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.axes.fill_between(
                        self.df['Wavenumber'], predicted, 0, alpha=0.2
                    )
        if self._bandsDict.shape[0] > 0:
            for index, row in self._bandsDict.iterrows():
                self.canvas.axes.annotate(
                    str(round(row['Wavenumber'], 0)), xy=(row['Wavenumber'], row[yAxis]), xytext=(0, -50),
                    textcoords='offset points', ha='center', va='bottom', rotation=90,
                    bbox=dict(boxstyle='round', color='white', alpha=0.5, ec="white"),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                )
        self.canvas.draw()
        self._bandsModel = IRDataModel(self._bandsDict, yAxis)
        self.uiIRBandsTableView.setModel(self._bandsModel)
        self._bandsModel.layoutChanged.emit()

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
        """SetGrid."""
        self._showGrid = self.uiShowGridButton.isChecked()
        if self.uiTransButton.isChecked():
            self.SelectYAxis('T')
        elif self.uiPercentTransButton.isChecked():
            self.SelectYAxis('%T')
        else:
            self.SelectYAxis('A')

    def SetNormalization(self):
        """SetNormalization.
        """
        if self.uiNormalizeButton.isChecked():
            self._normalized = True
        else:
            self._normalized = False
        self.SetYAxis()
        self.SelectYAxis('T')

    def SetYAxis(self):
        """
        docstring
        """
        if 0.0 < self.df['Raw Data'].median() < 1.0:
            self.df['Transmittance'] = self.df['Raw Data']
            self.df['%Transmittance'] = round(self.df['Transmittance'] * 100, 6)
            self.df['Absorbance'] = pd.Series([round(2.0 - math.log10(i), 6) for i in self.df['%Transmittance']])
        elif self.df['Raw Data'].min() < 1.0 and 2 > self.df['Raw Data'].max() > 0:
            self.df['Absorbance'] = self.df['Raw Data']
            self.df['Transmittance'] = pd.Series([round((10 ** (2 - i)) / 100, 6) for i in self.df['Absorbance']])
            self.df['%Transmittance'] = round(self.df['Transmittance'] * 100, 6)
        else:
            self.df['%Transmittance'] = self.df['Raw Data']
            self.df['Transmittance'] = round(self.df['%Transmittance'] / 100, 6)
            self.df['Absorbance'] = pd.Series([round(2.0 - math.log10(i), 6) for i in self.df['%Transmittance']])

        if self._normalized:
            minTrans = self.df['Transmittance'].min()
            maxTrans = self.df['Transmittance'].max()
            self.df['Transmittance'] = self.df['Transmittance'].apply(
                lambda x: 1.01 - (maxTrans - x) / (maxTrans - minTrans)
            )
            self.df['%Transmittance'] = self.df['Transmittance'] * 100
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
        self.uiSubsBaselineButton.setEnabled(False)
        self.uiBaseThresholdSpinBox.setEnabled(False)
        self.uiOrderBaseSpinBox.setEnabled(False)
        self._baseline = False
        if self._bandsDict.shape[0] > 0:
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
        self.uiNormalizeButton.setChecked(False)
        self._bandFitting = False
        self._normalized = False
        if self.df.get('Bands') is not None:
            del self.df['Bands']
        self._bandsDict = pd.DataFrame({
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
        self._bandsDict = self.df.loc[mask4, ['Wavenumber', 'Transmittance', '%Transmittance', 'Absorbance']].copy()
        self._bandsDict = self._bandsDict.reset_index(drop=True)
        self._bandsDict['Height (T)'] = pd.Series([1 - i for i in self._bandsDict['Transmittance']])
        self._bandsDict['Height (%T)'] = pd.Series([100 - i for i in self._bandsDict['%Transmittance']])
        self._bandsDict['Height (A)'] = self._bandsDict['Absorbance']
        self._bandsDict['HWHM'] = pd.Series(map(self.HWHM, self._bandsDict['Wavenumber'], self._bandsDict['Transmittance']))
        self.uiBandFittingAutoButton.setEnabled(True)
        if self.uiTransButton.isChecked():
            self.SelectYAxis('T')
        elif self.uiPercentTransButton.isChecked():
            self.SelectYAxis('%T')
        else:
            self.SelectYAxis('A')

    def BandsDetectionManual(self, event):  # y debe ser el valor correspondiente en la cura para el x seleccionado, no iy
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
                    'HWHM': [self.HWHM(ix, iy)]
                })
                self._bandsDict = self._bandsDict.append(newBand, ignore_index=True)
            elif self.uiPercentTransButton.isChecked():
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [iy / 100],
                    '%Transmittance': [iy],
                    'Absorbance': [2.0 - math.log10(iy)],
                    'Height (T)': [1 - (iy / 100)],
                    'Height (%T)': [100 - iy],
                    'Height (A)': [2.0 - math.log10(iy)],
                    'HWHM': [self.HWHM(ix, (iy / 100))]
                })
                self._bandsDict = self._bandsDict.append(newBand, ignore_index=True)
            elif self.uiAbsorbanceButton.isChecked():
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [(10 ** (2 - iy)) / 100],
                    '%Transmittance': [10 ** (2 - iy)],
                    'Absorbance': [iy],
                    'Height (T)': [1 - (10 ** (2 - iy)) / 100],
                    'Height (%T)': [100 - (10 ** (2 - iy))],
                    'Height (A)': [iy],
                    'HWHM': [self.HWHM(ix, (10 ** (2 - iy)) / 100)]
                })
                self._bandsDict = self._bandsDict.append(newBand, ignore_index=True)
            self.uiBandFittingAutoButton.setEnabled(True)
            if self.uiTransButton.isChecked():
                self.SelectYAxis('T')
            elif self.uiPercentTransButton.isChecked():
                self.SelectYAxis('%T')
            else:
                self.SelectYAxis('A')

    def HWHM(self, x, y):
        """HWHM.

        Parameters
        ----------
        x :
            x
        y :
            y
        """
        halfMax = round((1 + y) / 2, 6)

        mask = self.df['Wavenumber'].between(x, x + 10)
        mask2 = self.df['Wavenumber'].between(x - 10, x)

        bandsToFitX = self.df.loc[mask, 'Wavenumber'].values
        bandsToFitY = self.df.loc[mask, 'Transmittance'].values
        bandsToFitX2 = self.df.loc[mask2, 'Wavenumber'].values
        bandsToFitY2 = self.df.loc[mask2, 'Transmittance'].values

        coef = np.polyfit(bandsToFitX, bandsToFitY, 1)
        bandSteepFit = np.poly1d(coef)
        coef2 = np.polyfit(bandsToFitX2, bandsToFitY2, 1)
        bandSteepFit2 = np.poly1d(coef2)

        xRange = np.linspace(x - 500, x + 500, 100000)
        fitted = pd.DataFrame({'x': xRange, 'y1': bandSteepFit(xRange), 'y2': bandSteepFit2(xRange)})

        tolerance = 1e-5

        mask3 = fitted['y1'].between(halfMax - tolerance, halfMax + tolerance)
        mask4 = fitted['y2'].between(halfMax - tolerance, halfMax + tolerance)

        found = fitted.loc[mask3, 'x'].values
        found2 = fitted.loc[mask4, 'x'].values
        if len(found) > 0 and len(found2) > 0:
            hwhm1 = abs(x - found[0])
            hwhm2 = abs(x - found2[0])
        else:
            tolerance = 1e-4
            mask3 = fitted['y1'].between(
                halfMax - tolerance, halfMax + tolerance
            )
            mask4 = fitted['y2'].between(halfMax - tolerance, halfMax + tolerance)
            found = fitted.loc[mask3, 'x'].values
            found2 = fitted.loc[mask4, 'x'].values
            if len(found) > 0 and len(found2) > 0:
                hwhm1 = abs(x - found[0])
                hwhm2 = abs(x - found2[0])
            else:
                return 1.000

        mean = (hwhm1 + hwhm2) / 2
        print(abs((hwhm1 - hwhm2) * 100 / mean))
        return mean

    def BandFitting(self):
        """
        docstring
        """
        self.uiBandFittingAutoButton.setEnabled(False)
        self.worker = IRWorkerThread(self._bandsDict, self.df)
        self.thread = qtc.QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.Fit)
        self.worker.finished.connect(self.thread.quit)
        self.statusBar().showMessage(f'Now fitting {self._title} IR spectrum')
        self.worker.okFit.connect(self.OkFitManager)
        self.worker.error.connect(self.ErrorFitManager)
        self.thread.start()

    @qtc.pyqtSlot(dict, object, object)
    def OkFitManager(self, predictedBands, bands, df):

        self._predictedBands = predictedBands
        self.df = df
        self._bandsDict = bands
        self.statusBar().showMessage(f'Successful fitting for {self._title}')
        self._bandFitting = True

        if self.uiTransButton.isChecked():
            self.SelectYAxis('T')
        elif self.uiPercentTransButton.isChecked():
            self.SelectYAxis('%T')
        else:
            self.SelectYAxis('A')

    @qtc.pyqtSlot()
    def ErrorFitManager(self):
        self.statusBar().showMessage(f'Unsuccessful fitting for {self._title}')
        self.uiBandFittingAutoButton.setEnabled(True)

    def FirstDerivative(self):
        pass

    def SecondDerivative(self):
        pass
