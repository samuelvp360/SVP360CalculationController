#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import uic
import pandas as pd
import numpy as np
import nmrglue as ng
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
from scipy.signal import argrelextrema, savgol_filter
import math
from Models import IRDataModel, AvailableSpectraModel
from Worker import IRWorkerThread
from Views import resources
from datetime import datetime

matplotlib.use('Qt5Agg')


class IRCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(6, 4), dpi=100, facecolor='#2d2a2e')
        self.grid = GridSpec(17, 1, left=0.1, bottom=0.15, right=0.94, top=0.94, wspace=0.3, hspace=0.3)
        self.ax = self.fig.add_subplot(self.grid[0:, 0])
        self.ax2 = self.ax.twinx()
        super(IRCanvas, self).__init__(self.fig)


class NMRCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(6, 4), dpi=100, facecolor='#2d2a2e')
        self.grid = GridSpec(17, 1, left=0.1, bottom=0.15, right=0.94, top=0.94, wspace=0.3, hspace=0.3)
        self.ax = self.fig.add_subplot(self.grid[0:, 0])
        self.ax2 = self.ax.twinx()
        self.ax3 = self.ax.twinx()
        super(NMRCanvas, self).__init__(self.fig)


class IRPlotter(qtw.QMainWindow):
    """
    docstring
    """
    dataSent = qtc.pyqtSignal(dict)

    def __init__(
        self,
        title,
        solvent,
        csvData1=None,
        df=None,
        bandsDict=None,
        predictedBands=None
    ):
        super().__init__()
        uic.loadUi('Views/uiIRPlotter.ui', self)
        self.uiTransButton.setIcon(qtg.QIcon(':/icons/transIcon.png'))
        self.uiAbsorbanceButton.setIcon(qtg.QIcon(':/icons/absIcon.png'))
        self.uiPercentTransButton.setIcon(qtg.QIcon(':/icons/percentIcon.png'))
        self.uiDetectBaselineButton.setIcon(qtg.QIcon(':/icons/baselineIcon.png'))
        self.uiNormalizeButton.setIcon(qtg.QIcon(':/icons/normalizeIcon.png'))
        self.uiShowGridButton.setIcon(qtg.QIcon(':/icons/gridIcon.png'))

        self.canvas = IRCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiSpectraLayout.addWidget(self.canvas)
        self.cid = self.canvas.fig.canvas.mpl_connect('button_press_event', self.BandsDetectionManual)
        self._title = title
        self._solvent = solvent
        self.statusBar().showMessage(
            f'Showing IR spectrum from {self._title} ({self._solvent})'
        )
        if bandsDict is not None:
            self.uiBandFittingAutoButton.setEnabled(True)
            self.uiRemoveBandButton.setEnabled(True)

        self._baseThreshold = self.uiBaseThresholdSpinBox.value()
        self._bandThreshold = self.uiBandsThresholdSpinBox.value()
        self._orderBase = self.uiOrderBaseSpinBox.value()
        self._orderBands = self.uiOrderBandsSpinBox.value()
        self._showGrid = self.uiShowGridButton.isChecked()
        self._normalized = self.uiNormalizeButton.isChecked()
        self._baseline = False

        self.df = pd.read_csv(csvData1, names=('Wavenumber', 'Raw Data')) if df is None else df

        self._bandsDict = pd.DataFrame({
            'Wavenumber': [],
            'Transmittance': [],
            '%Transmittance': [],
            'Absorbance': [],
            'Transmittance(N)': [],
            '%Transmittance(N)': [],
            'Absorbance(N)': [],
            'HWHM': []
        }) if bandsDict is None else bandsDict
        self._predictedBands = {} if predictedBands is None else predictedBands
        self._bandFitting = False if predictedBands is None else True
        if csvData1 is not None:
            self.SetYAxis()
        self.SelectYAxis('Transmittance')
# ------------------------------------SIGNALS---------------------------------------------------------------
        self.uiDetectBaselineButton.clicked.connect(self.BaselineDetection)
        self.uiSubsBaselineButton.clicked.connect(self.BaselineSubs)
        self.uiBaseThresholdSpinBox.valueChanged.connect(self.SetBaseThreshold)
        self.uiBandsThresholdSpinBox.valueChanged.connect(self.SetBandsThreshold)
        self.uiOrderBaseSpinBox.valueChanged.connect(self.SetOrderBase)
        self.uiOrderBandsSpinBox.valueChanged.connect(self.SetOrderBands)
        self.uiShowGridButton.clicked.connect(self.SetGrid)
        self.uiTransButton.clicked.connect(lambda: self.SelectYAxis('Transmittance'))
        self.uiPercentTransButton.clicked.connect(lambda: self.SelectYAxis('%Transmittance'))
        self.uiAbsorbanceButton.clicked.connect(lambda: self.SelectYAxis('Absorbance'))
        self.uiRestoreButton.clicked.connect(self.Restore)
        self.uiBandAutoButton.clicked.connect(self.BandsDetectionAuto)
        self.uiBandFittingAutoButton.clicked.connect(self.BandFitting)
        self.uiNormalizeButton.clicked.connect(self.SetNormalization)
        self.uiSaveDataButton.clicked.connect(self.SaveData)
        self.uiRemoveBandButton.clicked.connect(self.RemoveBand)
# ------------------------------------METHODS---------------------------------------------------------------

    def Plot(self):
        """
        docstring
        """
        self.canvas.ax.clear()
        self.canvas.ax2.clear()
        self.canvas.ax.plot(self.df['Wavenumber'], self.df[self._yAxis], linewidth=1.5, color='#272822')
        self.canvas.ax.set_title(f'{self._title} ({self._solvent})', color='#ae81ff')
        self.canvas.ax.set_xlabel('Wavenumber [cm-1]', color='#f92672')
        self.canvas.ax.set_ylabel(self._yAxis, color='#f92672')
        self.canvas.ax.tick_params(axis='x', colors='#66d9ef')
        self.canvas.ax.tick_params(axis='y', colors='#66d9ef')
        self.canvas.ax2.set_ylabel('Residuals', color='#f92672')
        self.canvas.ax2.tick_params(axis='y', colors='#66d9ef')
        self.canvas.ax.set_xlim(4000, 500)
        self.canvas.ax.grid(True, color='#2d2a2e', linestyle=':', linewidth=0.5) if self._showGrid else self.canvas.ax.grid(False)
        if self._baseline:
            self.canvas.ax.plot(
                self.df['Wavenumber'], self._baselineDataFrame[self._yAxis],
                linewidth=1.0, color='#f92672'
            )
        if self._bandFitting:
            for _, predicted in self._predictedBands[self._yAxis].iteritems():
                self.canvas.ax.plot(self.df['Wavenumber'], predicted, linewidth=0.0)
                axis = self._yAxis.split('(')[0]
                if axis == 'Transmittance' or axis == '%Transmittance':
                    self.canvas.ax2.set_ylim(0.1, -0.005)
                else:
                    self.canvas.ax2.set_ylim(-0.005, 0.1)
                self.canvas.ax2.plot(
                    self.df['Wavenumber'], self._predictedBands['Residuals'],
                    '-r', markersize=0.2, linewidth=0.1
                )
                self.canvas.ax.plot(
                    self.df['Wavenumber'], self._predictedBands[self._yAxis]['Lorentzian'],
                    color='#ae81ff', linewidth=0.5, linestyle=':'
                )
                if self._yAxis == 'Transmittance' or self._yAxis == 'Transmittance(N)':
                    self.canvas.ax.fill_between(self.df['Wavenumber'], predicted, 1, alpha=0.2)
                elif self._yAxis == '%Transmittance' or self._yAxis == '%Transmittance(N)':
                    self.canvas.ax.fill_between(self.df['Wavenumber'], predicted, 100, alpha=0.2)
                else:
                    self.canvas.ax.fill_between(
                        self.df['Wavenumber'], predicted, 0, alpha=0.2
                    )
        if self._bandsDict.shape[0] > 0:
            for index, row in self._bandsDict.iterrows():
                axis = self._yAxis.split('(')[0]
                if axis == 'Transmittance' or axis == '%Transmittance':
                    offset = -50
                else:
                    offset = 15
                self.canvas.ax.annotate(
                    str(round(row['Wavenumber'], 0)), xy=(row['Wavenumber'], row[self._yAxis]),
                    xytext=(0, offset), textcoords='offset points', ha='center', va='bottom',
                    rotation=90, bbox=dict(boxstyle='round', color='white', alpha=0.5, ec="white"),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                )
        self.canvas.draw()
        self._bandsModel = IRDataModel(self._bandsDict, self._yAxis)
        self.uiIRBandsTableView.setModel(self._bandsModel)
        self.uiIRBandsTableView.resizeColumnsToContents()
        self.uiIRBandsTableView.resizeRowsToContents()
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
        self.Plot()

    def SetNormalization(self):
        """SetNormalization.
        """
        self._normalized = True if self.uiNormalizeButton.isChecked() else False
        self.SelectYAxis(self._yAxis) if self._normalized else self.SelectYAxis(self._yAxis.replace('(N)', ''))

    def SetYAxis(self):
        """
        docstring
        """
        if 0.0 < self.df['Raw Data'].median() < 1.0:
            self.df['Transmittance'] = self.df['Raw Data']
            self.df['Absorbance'] = self.df['Transmittance'].apply(self.TransToAbs)
            self.df['%Transmittance'] = self.df['Transmittance'].apply(self.TransToPercent)
        elif self.df['Raw Data'].min() < 1.0 and 2 > self.df['Raw Data'].max() > 0:
            self.df['Absorbance'] = self.df['Raw Data']
            self.df['Transmittance'] = self.df['Absorbance'].apply(self.AbsToTrans)
            self.df['%Transmittance'] = self.df['Transmittance'].apply(self.TransToPercent)
        else:
            self.df['%Transmittance'] = self.df['Raw Data']
            self.df['Transmittance'] = self.df['%Transmittance'].apply(self.TransToPercent)
            self.df['Absorbance'] = self.df['%Transmittance'].apply(self.PercentToAbs)
        self._minAbs = self.df['Absorbance'].min()
        self._maxAbs = self.df['Absorbance'].max()
        self._yRangeAbs = self._maxAbs - self._minAbs
        self._minTrans = self.df['Transmittance'].min()
        self._maxTrans = self.df['Transmittance'].max()
        self._yRangeTrans = self._maxTrans - self._minTrans
        self.df['Absorbance(N)'] = self.df['Absorbance'].apply(self.NormalizeA)
        self.df['Transmittance(N)'] = self.df['Transmittance'].apply(self.NormalizeT)
        self.df['%Transmittance(N)'] = self.df['Transmittance(N)'] * 100

    def SelectYAxis(self, unit):

        self.uiTransButton.setChecked(True) if unit == 'Transmittance' else self.uiTransButton.setChecked(False)
        self.uiPercentTransButton.setChecked(True) if unit == '%Transmittance' else self.uiPercentTransButton.setChecked(False)
        self.uiAbsorbanceButton.setChecked(True) if unit == 'Absorbance' else self.uiAbsorbanceButton.setChecked(False)
        self._yAxis = unit + '(N)' if self._normalized else unit
        self.Plot()

    def BaselineDetection(self):
        """
        docstring
        """
        self._baselineDataFrame = pd.DataFrame({'Wavenumber': self.df['Wavenumber'].to_list()})
        self._baselineDataFrame['Transmittance'] = self.df.loc[argrelextrema(self.df.Transmittance.values, np.greater_equal, order=self._orderBase)[0], 'Transmittance']
        mask1 = (self._baselineDataFrame['Transmittance'] < self._baseThreshold)
        self._baselineDataFrame.loc[mask1, 'Transmittance'] = np.nan
        mask2 = (self._baselineDataFrame['Transmittance'].isnull() == False)
        baselineX = self.df.loc[mask2, 'Wavenumber']
        baselineY = self.df.loc[mask2, 'Transmittance']
        coef = np.polyfit(baselineX, baselineY, 1)
        baselineFit = np.poly1d(coef)
        self._baselineDataFrame['Transmittance'] = pd.Series([baselineFit(i) for i in self.df['Wavenumber']])
        self._baselineDataFrame['Transmittance(N)'] = self._baselineDataFrame['Transmittance'].apply(self.NormalizeT)
        self._baselineDataFrame['%Transmittance'] = self._baselineDataFrame['Transmittance'] * 100
        self._baselineDataFrame['%Transmittance(N)'] = self._baselineDataFrame['Transmittance(N)'].apply(self.TransToPercent)
        self._baselineDataFrame['Absorbance'] = self._baselineDataFrame['Transmittance'].apply(self.TransToAbs)
        self._baselineDataFrame['Absorbance(N)'] = self._baselineDataFrame['Absorbance'].apply(self.NormalizeA)
        self._baseline = True
        self.uiSubsBaselineButton.setEnabled(True)
        self.uiDetectBaselineButton.setEnabled(False)
        self.Plot()

    def BaselineSubs(self):

        self.df['Transmittance'] = 1 - self._baselineDataFrame['Transmittance'] + self.df['Transmittance']
        self.df['%Transmittance'] = self.df['Transmittance'].apply(self.TransToPercent)
        self.df['Absorbance'] = self.df['Transmittance'].apply(self.TransToAbs)
        self._minAbs = self.df['Absorbance'].min()
        self._maxAbs = self.df['Absorbance'].max()
        self._yRangeAbs = self._maxAbs - self._minAbs
        self._minTrans = self.df['Transmittance'].min()
        self._maxTrans = self.df['Transmittance'].max()
        self._yRangeTrans = self._maxTrans - self._minTrans
        self.df['Transmittance(N)'] = self.df['Transmittance'].apply(self.NormalizeT)
        self.df['%Transmittance(N)'] = self.df['Transmittance(N)'].apply(self.TransToPercent)
        self.df['Absorbance(N)'] = self.df['Absorbance'].apply(self.NormalizeA)
        self.uiSubsBaselineButton.setEnabled(False)
        self.uiBaseThresholdSpinBox.setEnabled(False)
        self.uiOrderBaseSpinBox.setEnabled(False)
        self._baseline = False
        if self._bandsDict.shape[0] > 0:
            self.BandsDetectionAuto()
        else:
            self.Plot()

    def TransToPercent(self, y):
        return round(y * 100, 6)

    def TransToAbs(self, y):
        return round(2.0 - math.log10(y * 100), 6)

    def PercentToTrans(self, y):
        return round(y / 100, 6)

    def PercentToAbs(self, y):
        return round(2.0 - math.log10(y), 6)

    def AbsToTrans(self, y):
        return round((10 ** (2 - y)) / 100, 6)

    def AbsToPercent(self, y):
        return round((10 ** (2 - y)), 6)

    def NormalizeT(self, y):
        return round((1.0 - (self._maxTrans - y) / self._yRangeTrans), 6)

    def NormalizeA(self, y):
        return round(((y - self._minAbs) / self._yRangeAbs), 6)

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
        self.uiRemoveBandButton.setEnabled(False)
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
            'Transmittance(N)': [],
            '%Transmittance(N)': [],
            'Absorbance(N)': [],
        })
        self._predictedBands = {}
        self.SetYAxis()
        self.SelectYAxis('Transmittance')

    def BandsDetectionAuto(self):
        """
        docstring
        """
        self.df['Bands'] = self.df.loc[argrelextrema(self.df.Transmittance.values, np.less_equal, order=self._orderBands)[0], 'Transmittance']
        mask3 = (self.df['Bands'] > self._bandThreshold)
        self.df.loc[mask3, 'Bands'] = np.nan
        mask4 = (self.df['Bands'].isnull() == False)
        self._bandsDict = self.df.loc[
            mask4, [
                'Wavenumber', 'Transmittance', '%Transmittance', 'Absorbance',
                'Transmittance(N)', '%Transmittance(N)', 'Absorbance(N)'
            ]
        ].copy()
        self._bandsDict = self._bandsDict.reset_index(drop=True)
        self._bandsDict['HWHM'] = pd.Series(
            map(self.HWHM, self._bandsDict['Wavenumber'], self._bandsDict['Transmittance'])
        )
        self._bandsDict = self._bandsDict.round(6)
        self.uiBandFittingAutoButton.setEnabled(True)
        self.uiRemoveBandButton.setEnabled(True)
        self.Plot()

    def BandsDetectionManual(self, event):
        """
        docstring
        """
        if self.uiBandManualButton.isChecked():
            ix, _ = event.xdata, event.ydata
            tolerance = 1  # this value can be moved to avoid errors if the x
            # sample data changes.
            mask = self.df['Wavenumber'].between(ix - tolerance, ix + tolerance)
            try:
                transmittance = self.df.loc[mask, 'Transmittance'].values[0]
                transmittanceN = self.df.loc[mask, 'Transmittance(N)'].values[0]
                percent = self.df.loc[mask, '%Transmittance'].values[0]
                percentN = self.df.loc[mask, '%Transmittance(N)'].values[0]
                absorbance = self.df.loc[mask, 'Absorbance'].values[0]
                absorbanceN = self.df.loc[mask, 'Absorbance(N)'].values[0]
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [transmittance],
                    '%Transmittance': [percent],
                    'Absorbance': [absorbance],
                    'Transmittance(N)': [transmittanceN],
                    '%Transmittance(N)': [percentN],
                    'Absorbance(N)': [absorbanceN],
                    'HWHM': [self.HWHM(ix, transmittance)]
                })
                self._bandsDict = self._bandsDict.append(newBand, ignore_index=True)
                self._bandsDict = self._bandsDict.round(6)
            except IndexError:
                pass
            self.uiBandFittingAutoButton.setEnabled(True)
            self.uiRemoveBandButton.setEnabled(True)
            self.Plot()

    def RemoveBand(self):

        indexes = self.uiIRBandsTableView.selectedIndexes()
        if indexes:
            index = indexes[0]
            self._bandsDict.drop([index.row()], inplace=True)
            self._bandsDict.reset_index(drop=True)
            self._bandsModel.layoutChanged.emit()
            self.Plot()

        if self._bandFitting:
            reply = qtw.QMessageBox.question(
                self, 'Redo Fitting',
                'You have removed a band, do you want to redo the fitting?',
                qtw.QMessageBox.Yes | qtw.QMessageBox.No,
                qtw.QMessageBox.No
            )

            if reply == qtw.QMessageBox.Yes:
                self.BandFitting()

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
        self.uiRemoveBandButton.setEnabled(False)
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
        print(self._bandsDict)
        self._bandsDict = self._bandsDict.round(6)
        self._predictedBands['Absorbance(N)'] = self._predictedBands['Absorbance'].applymap(self.NormalizeA)
        self._predictedBands['Transmittance'] = self._predictedBands['Absorbance'].applymap(self.AbsToTrans)
        self._predictedBands['Transmittance(N)'] = self._predictedBands['Transmittance'].applymap(self.NormalizeT)
        self._predictedBands['%Transmittance'] = self._predictedBands['Transmittance'].applymap(self.TransToPercent)
        self._predictedBands['%Transmittance(N)'] = self._predictedBands['Transmittance(N)'].applymap(lambda x: x * 100)
        self._predictedBands['Absorbance']['Lorentzian'] = self._predictedBands['Absorbance'].sum(axis=1)
        self._predictedBands['Absorbance(N)']['Lorentzian'] = self._predictedBands['Absorbance']['Lorentzian'].apply(self.NormalizeA)
        self._predictedBands['Transmittance']['Lorentzian'] = self._predictedBands['Absorbance']['Lorentzian'].apply(self.AbsToTrans)
        self._predictedBands['Transmittance(N)']['Lorentzian'] = self._predictedBands['Transmittance']['Lorentzian'].apply(self.NormalizeT)
        self._predictedBands['%Transmittance']['Lorentzian'] = self._predictedBands['Transmittance']['Lorentzian'] * 100
        self._predictedBands['%Transmittance(N)']['Lorentzian'] = self._predictedBands['Transmittance(N)']['Lorentzian'] * 100
        self._predictedBands['Residuals'] = self.df['Absorbance'] - self._predictedBands['Absorbance']['Lorentzian']
        self.statusBar().showMessage(f'Successful fitting of {self._title}')
        self._bandFitting = True
        self.Plot()

    @qtc.pyqtSlot()
    def ErrorFitManager(self):
        self.statusBar().showMessage(f'Unsuccessful fitting for {self._title}')
        self.uiBandFittingAutoButton.setEnabled(True)
        self.uiRemoveBandButton.setEnabled(True)

    def FirstDerivative(self):
        pass

    def SecondDerivative(self):
        pass

    def SaveData(self):
        date = datetime.now()
        bandsDict = self._bandsDict if self._bandsDict.shape[0] > 0 else None
        predicted = self._predictedBands if len(self._predictedBands.items()) > 0 else None
        self.dataSent.emit({
            'TYPE': 'FTIR',
            'DATE': str(date),
            'SOLVENT': self._solvent,
            'DATA FRAME': self.df,
            'BANDS': bandsDict,
            'PREDICTED': predicted
        })
        reply = qtw.QMessageBox.question(
            self, 'Window Close',
            'The spectrum has been stored, do you want to close the window?',
            qtw.QMessageBox.Yes | qtw.QMessageBox.No,
            qtw.QMessageBox.No
        )

        if reply == qtw.QMessageBox.Yes:
            self.close()


class SpectrumSelector(qtw.QWidget):
    """
    docstring
    """
    def __init__(self, molecule):
        super(SpectrumSelector, self).__init__()
        uic.loadUi('Views/uiSpectrumSelector.ui', self)
        self._molecule = molecule
        self._model = AvailableSpectraModel(self._molecule.GetExperimental)
        self.uiAvailableSpectraTableView.setModel(self._model)
        self.uiAvailableSpectraTableView.resizeColumnsToContents()
        self.uiAvailableSpectraTableView.resizeRowsToContents()
        self.uiLoadIRButton.clicked.connect(self.LoadIR)
        self.uiLoadNMRButton.clicked.connect(self.LoadNMR)
        self.uiRemoveButton.clicked.connect(self.RemoveSpectrum)
        self.uiShowButton.clicked.connect(self.ShowSpectrum)
        if len(self._molecule.GetExperimental.items()) > 0:
            self.uiRemoveButton.setEnabled(True)
            self.uiShowButton.setEnabled(True)

    def LoadIR(self):

        filePath, ok1 = qtw.QFileDialog.getOpenFileName(
            self,
            'Select your IR csv data',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='CSV files(*.csv)'
        )

        slv = ['KBr', 'CHCl3', 'CH2Cl2', 'Neat', 'Other']
        if ok1:
            solvent, ok2 = qtw.QInputDialog.getItem(
                self, 'Solvent', 'Please select the used solvent', slv, 0, False
            )

        if ok1 and ok2 and solvent:
            self.IRPlotter = IRPlotter(
                self._molecule.GetName, solvent, csvData1=filePath
            )
            self.IRPlotter.show()
            self.IRPlotter.dataSent.connect(self.SaveSpectrum)

    def LoadNMR(self):

        dirPath = qtw.QFileDialog.getExistingDirectory(
            self,
            'Select your NMR data',
            options=qtw.QFileDialog.ShowDirsOnly,
        )
        self.NMRPlotter = NMRPlotter(
            self._molecule.GetName, dataDir=dirPath
        )
        self.NMRPlotter.show()
        self.NMRPlotter.dataSent.connect(self.SaveSpectrum)

    def ShowSpectrum(self):

        indexes = self.uiAvailableSpectraTableView.selectedIndexes()
        if indexes:
            index = indexes[0]
            keys = [i for i in self._molecule.GetExperimental.keys()]
            data = self._molecule.GetExperimental.get(keys[index.row()])
            if data['TYPE'] == 'FTIR':
                solvent = data['SOLVENT']
                dataFrame = data['DATA FRAME']
                bands = data['BANDS']
                predicted = data['PREDICTED']
                self.IRPlotter = IRPlotter(
                    self._molecule.GetName,
                    solvent,
                    df=dataFrame,
                    bandsDict=bands,
                    predictedBands=predicted
                )
                self.IRPlotter.show()
                self.IRPlotter.dataSent.connect(self.SaveSpectrum)
            elif data['TYPE'] == 'NMR 1H':
                pass

    @qtc.pyqtSlot(dict)
    def SaveSpectrum(self, data):

        self._molecule.SetIR(data)
        self._model.layoutChanged.emit()
        self.uiAvailableSpectraTableView.resizeColumnsToContents()
        self.uiAvailableSpectraTableView.resizeRowsToContents()
        self.uiShowButton.setEnabled(True)
        self.uiRemoveButton.setEnabled(True)

    def RemoveSpectrum(self):

        indexes = self.uiAvailableSpectraTableView.selectedIndexes()
        if indexes:
            index = indexes[0]
            self._molecule.RemoveExperimental(index.row())
            self._model.layoutChanged.emit()
            self.uiAvailableSpectraTableView.resizeColumnsToContents()
            self.uiAvailableSpectraTableView.resizeRowsToContents()

        if len(self._molecule.GetExperimental.items()) == 0:
            self.uiRemoveButton.setEnabled(False)
            self.uiShowButton.setEnabled(False)


class NMRPlotter(qtw.QMainWindow):
    """
    docstring
    """
    dataSent = qtc.pyqtSignal(dict)

    def __init__(
        self,
        title,
        dataDir=None,
        data=None,
        udic=None,
        ppm=None,
        bandsDict=None,
        predictedBands=None
    ):
        super(NMRPlotter, self).__init__()
        uic.loadUi('Views/uiNMRPlotter.ui', self)
        self.canvas = NMRCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiSpectraLayout.addWidget(self.canvas)
        self.cid = self.canvas.fig.canvas.mpl_connect(
            'button_press_event', self.PeaksDetectionManual
        )
        self.GSD = GSD()
        self._title = title
        self._showGrid = True
        self._smoothed = False
        self._xAxis = 'ppm'
        self._yAxis = 'Intensity'
        if dataDir is not None:
            dic, data = ng.bruker.read(dataDir)
            data = ng.bruker.remove_digital_filter(dic, data)
            self.udic = ng.bruker.guess_udic(dic, data)
            data = ng.proc_base.zf_size(data, self.udic[0]['size'])
            data = ng.proc_base.fft(data)
            data = ng.proc_base.ps(data, p0=-10.0)  # este valor debe
            # ir cambiado, y esto debe ser opcional
            data = ng.proc_base.di(data)
            data = ng.proc_base.rev(data)
            uc = ng.fileiobase.uc_from_udic(self.udic)
            ppm = uc.ppm_scale()
            hz = uc.hz_scale()
            self.df = pd.DataFrame({'Intensity': data})
            self.df['ppm'] = ppm
            self.df['hz'] = hz
            self.df['First'] = self.df['Intensity'].diff() / self.df['ppm'].diff()
            self.df['Second'] = self.df['First'].diff() / self.df['ppm'].diff()
        else:
            self.df = data
            self.udic = udic
        self._nmrType = self.udic[0]['label']
        self.statusBar().showMessage(f'Showing {self._nmrType} spectra for {self._title}')
        self.Plot()
# ------------------------------------SIGNALS---------------------------------------------------------------
        self.uiShowGridButton.clicked.connect(self.SetGrid)
        self.uiPpmButton.clicked.connect(lambda: self.SelectXAxis('ppm'))
        self.uiHzButton.clicked.connect(lambda: self.SelectXAxis('hz'))
        self.uiSmoothButton.clicked.connect(self.Smooth)
# ------------------------------------METHODS---------------------------------------------------------------

    def Plot(self):

        self.canvas.ax.clear()
        self.canvas.ax2.clear()
        self.canvas.ax3.clear()
        self.canvas.ax.plot(self.df[self._xAxis], self.df[self._yAxis], linewidth=0.5, color='#272822')
        self.canvas.ax.set_title(f'{self._title}', color='#ae81ff')
        if self._xAxis == 'ppm':
            self.canvas.ax.set_xlabel('Chemical Shift (ppm)', color='#f92672')
        else:
            self.canvas.ax.set_xlabel('Frequency (Hz)', color='#f92672')
        self.canvas.ax.set_ylabel('Intensity', color='#f92672')
        self.canvas.ax.tick_params(axis='x', colors='#66d9ef')
        self.canvas.ax.tick_params(axis='y', colors='#66d9ef')
        self.canvas.ax2.set_ylabel('First derivative', color='#f92672')
        self.canvas.ax2.tick_params(axis='y', colors='#66d9ef')
        self.canvas.ax2.plot(self.df[self._xAxis], self.df['First'], linewidth=0.5, color='red')
        self.canvas.ax3.set_ylabel('Second derivative', color='#f92672')
        self.canvas.ax3.tick_params(axis='y', colors='#66d9ef')
        self.canvas.ax3.plot(self.df[self._xAxis], self.df['Second'], linewidth=0.5, color='blue')
        if self._nmrType == '1H' and self._xAxis == 'ppm':
            self.canvas.ax.set_xlim(15, -1)
            self.canvas.ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
        elif self._nmrType == '13C' and self._xAxis == 'ppm':
            self.canvas.ax.set_xlim(200, -10)
            self.canvas.ax.xaxis.set_major_locator(ticker.MultipleLocator(10.0))
        else:
            self.canvas.ax.set_xlim(18000, -6000)
        self.canvas.ax.grid(True, color='#2d2a2e', linestyle=':', linewidth=0.2) if self._showGrid else self.canvas.ax.grid(False)
        self.canvas.draw()

    def PeaksDetectionManual(self, event):
        pass

    def PeaksDetectionAuto(self):
        pass

    def SetGrid(self):

        self._showGrid = self.uiShowGridButton.isChecked()
        self.Plot()

    def SelectXAxis(self, unit):

        self.uiPpmButton.setChecked(True) if unit == 'ppm' else self.uiPpmButton.setChecked(False)
        self.uiHzButton.setChecked(True) if unit == 'hz' else self.uiHzButton.setChecked(False)
        self._xAxis = unit
        self.Plot()

    def Smooth(self):
        if self.uiSmoothButton.isChecked():
            self.df['smoothed'] = self.GSD.SavGol(self.df['Intensity'])
            self._yAxis = 'smoothed'
        else:
            self._yAxis = 'Intensity'
        self.Plot()


class GSD():

    def __init__(self):
        pass
        # self.data = data
        # self.smoothedPpm = self.SavGol(data['ppm'])
        # self.smoothedHz = self.SavGol(data['hz'])

    def SavGol(self, yRange):
        return savgol_filter(yRange, 5, 3)

