#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import uic
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from scipy.signal import argrelextrema
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
        self.fig = Figure(figsize=(22, 4), dpi=100, facecolor='#2d2a2e')
        self.grid = GridSpec(17, 1, left=0.1, bottom=0.15, right=0.94, top=0.94, wspace=0.3, hspace=0.3)
        self.ax = self.fig.add_subplot(self.grid[0:11, 0])
        self.ax2 = self.fig.add_subplot(self.grid[14:, 0])
        super(IRCanvas, self).__init__(self.fig)


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
        super(IRPlotter, self).__init__()
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
        self.statusBar().showMessage(f'Showing IR spectrum from {self._title}')
        if bandsDict is not None:
            self.uiBandFittingAutoButton.setEnabled(True)
            self.uiRemoveBandButton.setEnabled(True)

        self._baseThreshold = self.uiBaseThresholdSpinBox.value()
        self._bandThreshold = self.uiBandsThresholdSpinBox.value()
        self._orderBase = self.uiOrderBaseSpinBox.value()
        self._orderBands = self.uiOrderBandsSpinBox.value()
        self._showGrid = self.uiShowGridButton.isChecked()
        self._baseline = False
        self._normalized = False

        self.df = pd.read_csv(csvData1, names=('Wavenumber', 'Raw Data')) if df is None else df

        self._bandsDict = pd.DataFrame({
            'Wavenumber': [],
            'Transmittance': [],
            '%Transmittance': [],
            'Absorbance': [],
            'Height (T)': [],
            'Height (%T)': [],
            'Height (A)': []
        }) if bandsDict is None else bandsDict
        self._predictedBands = {} if predictedBands is None else predictedBands
        self._bandFitting = False if predictedBands is None else True
        if csvData1 is not None:
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
        self.uiSaveDataButton.clicked.connect(self.SaveData)
        self.uiRemoveBandButton.clicked.connect(self.RemoveBand)
# ------------------------------------METHODS---------------------------------------------------------------

    def Plot(self, yAxis):
        """
        docstring
        """
        self.canvas.ax.clear()
        self.canvas.ax2.clear()
        self.canvas.ax.plot(self.df['Wavenumber'], self.df[yAxis], linewidth=1.5, color='#272822')
        self.canvas.ax.set_title(self._title, color='#ae81ff')
        self.canvas.ax.set_xlabel('Wavenumber [cm-1]', color='#f92672')
        self.canvas.ax.set_ylabel(yAxis, color='#f92672')
        self.canvas.ax.tick_params(axis='x', colors='#66d9ef')
        self.canvas.ax.tick_params(axis='y', colors='#66d9ef')
        self.canvas.ax2.set_ylabel('Residuals', color='#f92672')
        self.canvas.ax2.tick_params(axis='x', colors='#66d9ef')
        self.canvas.ax2.tick_params(axis='y', colors='#66d9ef')
        self.canvas.ax2.axhline()
        self.canvas.ax.set_xlim(4000, 500)
        self.canvas.ax2.set_xlim(4000, 500)
        self.canvas.ax.grid(True, color='#2d2a2e', linestyle=':', linewidth=0.5) if self._showGrid else self.canvas.ax.grid(False)
        if self._baseline:
            self.canvas.ax.plot(self.df['Wavenumber'], self.df['baseline'], linewidth=1.0, color='#f92672')
        if self._bandFitting:
            for _, predicted in self._predictedBands[yAxis].iteritems():
                self.canvas.ax.plot(self.df['Wavenumber'], predicted, linewidth=0.0)
                self.canvas.ax2.plot(
                    self.df['Wavenumber'], self._predictedBands['Residuals'], 'o',
                    markersize=0.2, linewidth=0.1, color='#ae81ff'
                )
                if yAxis == 'Transmittance':
                    self.canvas.ax.plot(
                        self.df['Wavenumber'], self._predictedBands['Transmittance']['Fitted Curve Lorentzian'],
                        color='#ae81ff', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.ax.fill_between(self.df['Wavenumber'], predicted, 1, alpha=0.2)
                elif yAxis == '%Transmittance':
                    self.canvas.ax.plot(
                        self.df['Wavenumber'], self._predictedBands['%Transmittance']['Fitted Curve Lorentzian'],
                        color='#ae81ff', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.ax.fill_between(self.df['Wavenumber'], predicted, 100, alpha=0.2)
                else:
                    self.canvas.ax.plot(
                        self.df['Wavenumber'], self._predictedBands['Absorbance']['Fitted Curve Lorentzian'],
                        color='#ae81ff', linewidth=0.5, linestyle=':'
                    )
                    self.canvas.ax.fill_between(
                        self.df['Wavenumber'], predicted, 0, alpha=0.2
                    )
        if self._bandsDict.shape[0] > 0:
            for index, row in self._bandsDict.iterrows():
                self.canvas.ax.annotate(
                    str(round(row['Wavenumber'], 0)), xy=(row['Wavenumber'], row[yAxis]), xytext=(0, -50),
                    textcoords='offset points', ha='center', va='bottom', rotation=90,
                    bbox=dict(boxstyle='round', color='white', alpha=0.5, ec="white"),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                )
        self.canvas.draw()
        self._bandsModel = IRDataModel(self._bandsDict, yAxis)
        self.uiIRBandsTableView.setModel(self._bandsModel)
        self.uiIRBandsTableView.resizeColumnsToContents()
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
        self._bandsDict = self._bandsDict.round(6)
        self.uiBandFittingAutoButton.setEnabled(True)
        self.uiRemoveBandButton.setEnabled(True)
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
            ix, _ = event.xdata, event.ydata
            tolerance = 1  # this value can be moved to avoid errors if the x
            # sample data changes.
            mask = self.df['Wavenumber'].between(ix - tolerance, ix + tolerance)
            if self.uiTransButton.isChecked():
                found = self.df.loc[mask, 'Transmittance'].values
                y = found[0]
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [y],
                    '%Transmittance': [y * 100],
                    'Absorbance': [2.0 - math.log10(y * 100)],
                    'Height (T)': [1 - y],
                    'Height (%T)': [100 - (y * 100)],
                    'Height (A)': [2.0 - math.log10(y * 100)],
                    'HWHM': [self.HWHM(ix, y)]
                })
                self._bandsDict = self._bandsDict.append(newBand, ignore_index=True)
            elif self.uiPercentTransButton.isChecked():
                found = self.df.loc[mask, '%Transmittance'].values
                y = found[0]
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [y / 100],
                    '%Transmittance': [y],
                    'Absorbance': [2.0 - math.log10(y)],
                    'Height (T)': [1 - (y / 100)],
                    'Height (%T)': [100 - y],
                    'Height (A)': [2.0 - math.log10(y)],
                    'HWHM': [self.HWHM(ix, (y / 100))]
                })
                self._bandsDict = self._bandsDict.append(newBand, ignore_index=True)
            elif self.uiAbsorbanceButton.isChecked():
                found = self.df.loc[mask, '%Transmittance'].values
                y = found[0]
                newBand = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [(10 ** (2 - y)) / 100],
                    '%Transmittance': [10 ** (2 - y)],
                    'Absorbance': [y],
                    'Height (T)': [1 - (10 ** (2 - y)) / 100],
                    'Height (%T)': [100 - (10 ** (2 - y))],
                    'Height (A)': [y],
                    'HWHM': [self.HWHM(ix, (10 ** (2 - y)) / 100)]
                })
                self._bandsDict = self._bandsDict.append(newBand, ignore_index=True)
            self._bandsDict = self._bandsDict.round(6)
            self.uiBandFittingAutoButton.setEnabled(True)
            self.uiRemoveBandButton.setEnabled(True)
            if self.uiTransButton.isChecked():
                self.SelectYAxis('T')
            elif self.uiPercentTransButton.isChecked():
                self.SelectYAxis('%T')
            else:
                self.SelectYAxis('A')

    def RemoveBand(self):

        indexes = self.uiIRBandsTableView.selectedIndexes()
        if indexes:
            index = indexes[0]
            self._bandsDict.drop([index.row()], inplace=True)
            self._bandsDict.reset_index(drop=True)
            self._bandsModel.layoutChanged.emit()
            if self.uiTransButton.isChecked():
                self.SelectYAxis('T')
            elif self.uiPercentTransButton.isChecked():
                self.SelectYAxis('%T')
            else:
                self.SelectYAxis('A')

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
        self.uiLoadButton.clicked.connect(self.LoadSpectrum)
        self.uiRemoveButton.clicked.connect(self.RemoveSpectrum)
        self.uiShowButton.clicked.connect(self.ShowSpectrum)
        if len(self._molecule.GetExperimental.items()) > 0:
            self.uiRemoveButton.setEnabled(True)
            self.uiShowButton.setEnabled(True)

    def LoadSpectrum(self):

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

    @qtc.pyqtSlot(dict)
    def SaveSpectrum(self, data):

        self._molecule.SetIR(data)
        self._model.layoutChanged.emit()
        self.uiAvailableSpectraTableView.resizeColumnsToContents()
        self.uiShowButton.setEnabled(True)
        self.uiRemoveButton.setEnabled(True)

    def RemoveSpectrum(self):

        indexes = self.uiAvailableSpectraTableView.selectedIndexes()
        if indexes:
            index = indexes[0]
            self._molecule.RemoveExperimental(index.row())
            self._model.layoutChanged.emit()
            self.uiAvailableSpectraTableView.resizeColumnsToContents()

        if len(self._molecule.GetExperimental.items()) == 0:
            self.uiRemoveButton.setEnabled(False)
            self.uiShowButton.setEnabled(False)
