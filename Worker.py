#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
import pandas as pd
from scipy.optimize import curve_fit


class WorkerThread(qtc.QObject):

    finished = qtc.pyqtSignal()
    startingCalc = qtc.pyqtSignal(object, int, str)
    okCalc = qtc.pyqtSignal(object, int, str)
    error = qtc.pyqtSignal(object, int, str)

    def __init__(self, masterQueue):
        super(WorkerThread, self).__init__(parent=None)
        self.masterQueue = masterQueue
        self.isStopped = False

    @qtc.pyqtSlot()
    def Calculate(self):

        try:
            for index, calc in enumerate(self.masterQueue):
                self.molecule, *self.parameters = calc
                if self.parameters[2] == 'Pending':
                    self.startingCalc.emit(self.molecule, index, 'Running')
                    self.molecule.SetCalculations(*self.parameters)
                    self.okCalc.emit(self.molecule, index, 'Finished')
                    print(self.molecule.GetCalculations)
                if self.isStopped:
                    break
        except:
            self.error.emit(self.molecule, index, 'Aborted')
        finally:
            self.finished.emit()
            self.deleteLater()

    def Stop(self):
        self.isStopped = True

    def Pause(self):
        self.isPaused = True


class IRWorkerThread(qtc.QObject):

    finished = qtc.pyqtSignal()
    okFit = qtc.pyqtSignal(dict, object, object)
    error = qtc.pyqtSignal()

    def __init__(self, bands, df):
        super(IRWorkerThread, self).__init__(parent=None)
        self.bands = bands
        self.predictedBands = {}
        self.df = df

    @qtc.pyqtSlot()
    def Fit(self):

        n = self.bands.shape[0]
        position = self.bands['Wavenumber'].to_list()
        amplitude = self.bands['Height (A)'].to_list()
        hwhm = self.bands['HWHM'].to_list()
        parameters = position + amplitude + hwhm
        lowerBounds = tuple(n * [500] + n * [0] + n * [1e-5])
        upperBounds = tuple(n * [4000] + n * [2] + n * [1000])

        try:
            fittedParams, cov = curve_fit(
                self.LorentzianModel, self.df['Wavenumber'],
                self.df['Absorbance'], p0=parameters, absolute_sigma=True,
                bounds=(lowerBounds, upperBounds)
            )
            fittedPosition = fittedParams[:n]
            fittedAmplitude = fittedParams[n:2 * n]
            fittedHwhm = fittedParams[2 * n:]
            self.bands['Wavenumber'] = pd.Series(fittedPosition)
            self.bands['Height (A)'] = pd.Series(fittedAmplitude)
            self.bands['HWHM'] = pd.Series(fittedHwhm)
            self.predictedBands['Absorbance'] = pd.DataFrame(
                map(self.LorentzianByBand, self.bands['Wavenumber'], self.bands['Height (A)'], self.bands['HWHM'])
            ).transpose()
            self.predictedBands['Transmittance'] = self.predictedBands['Absorbance'].applymap(lambda x: (10 ** (2 - x)) / 100)
            self.predictedBands['%Transmittance'] = self.predictedBands['Transmittance'].applymap(lambda x: 100 * x)
            self.predictedBands['Absorbance']['Fitted Curve Lorentzian'] = self.predictedBands['Absorbance'].sum(axis=1)
            self.predictedBands['Residuals'] = self.df['Absorbance'] - self.predictedBands['Absorbance']['Fitted Curve Lorentzian']
            self.predictedBands['Transmittance']['Fitted Curve Lorentzian'] = self.predictedBands['Absorbance']['Fitted Curve Lorentzian'].apply(lambda x: (10 ** (2 - x)) / 100)
            self.predictedBands['%Transmittance']['Fitted Curve Lorentzian'] = self.predictedBands['Transmittance']['Fitted Curve Lorentzian'].apply(lambda x: 100 * x)
            self.okFit.emit(self.predictedBands, self.bands, self.df)
        except RuntimeError:
            self.error.emit()
        finally:
            self.finished.emit()
            self.deleteLater()

    def LorentzianFunction(self, x, x0, amp, hwhm):
        return amp / (1 + ((x - x0) / hwhm) ** 2)

    def LorentzianByBand(self, x0, amp, hwhm):
        """
        docstring
        """
        xRange = self.df['Wavenumber']
        return pd.Series([(amp / (1 + ((x - x0) / hwhm) ** 2)) for x in xRange])

    def LorentzianModel(self, x, *args):
        """LorentzianModel.

        Args:
            x:
            args:
        """
        n = self.bands.shape[0]
        x0 = args[:n]
        amp = args[n:2 * n]
        hwhm = args[2 * n:]

        res = [self.LorentzianFunction(x, x0[i], amp[i], hwhm[i]) for i in range(n)]

        return sum(res)


