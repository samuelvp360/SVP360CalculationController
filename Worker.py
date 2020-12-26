#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
import pandas as pd
from scipy.optimize import curve_fit
from numpy import exp, sqrt


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
        amplitude = self.bands['Absorbance'].to_list()
        hwhm = self.bands['HWHM'].to_list()
        beta = n * [0.5]
        parameters = position + amplitude + hwhm + beta
        lowerBounds = tuple([i - 50 for i in position] + [i / 2.0 for i in amplitude] + n * [1e-5] + n * [0])
        upperBounds = tuple([i + 50 for i in position] + [i * 1.5 for i in amplitude] + n * [1000] + n * [1])

        try:
            fittedParams, cov = curve_fit(
                self.FitModel, self.df['Wavenumber'],
                self.df['Absorbance'], p0=parameters, absolute_sigma=True,
                bounds=(lowerBounds, upperBounds)
            )
            fittedPosition = fittedParams[:n]
            fittedAmplitude = fittedParams[n:2 * n]
            fittedHwhm = fittedParams[2 * n:3 * n]
            fittedBeta = fittedParams[3 * n:]
            self.bands['Wavenumber'] = pd.Series(fittedPosition)
            self.bands['Absorbance'] = pd.Series(fittedAmplitude)
            self.bands['HWHM'] = pd.Series(fittedHwhm)
            self.bands['Beta'] = pd.Series(fittedBeta)
            self.predictedBands['Absorbance'] = pd.DataFrame(
                map(
                    self.ModelByBand, self.bands['Wavenumber'],
                    self.bands['Absorbance'], self.bands['HWHM'],
                    self.bands['Beta']
                )
            ).transpose()
            self.okFit.emit(self.predictedBands, self.bands, self.df)
        except RuntimeError:
            self.error.emit()
        finally:
            self.finished.emit()
            self.deleteLater()

    def L(self, x, x0, amp, hwhm):
        """Lorentzian Function.

        Parameters
        ----------
        x :
            frequency variable
        x0 :
            position of the band
        amp :
            height of the band
        hwhm :
            half widht at half maximum of the band
        """
        print('Lorentzian')
        return amp / (1 + ((x - x0) / hwhm) ** 2)

    def G(self, x, x0, amp, hwhm):
        """Gaussian Function.

        Parameters
        ----------
        x :
            frequency variable
        x0 :
            position of the band
        amp :
            height of the band
        hwhm :
            half widht at half maximum of the band
        """
        print('Gaussian')
        return amp * exp(-1 * (x - x0) ** 2 / (2 * hwhm ** 2))

    def K(self, x, x0, amp, hwhm, beta):
        """Brown Function.

        Parameters
        ----------
        x :
            frequency variable
        x0 :
            position of the band
        amp :
            height of the band
        hwhm :
            half widht at half maximum of the band
        beta :
            Voigt-likeness parameter.
            If beta=0, K resembles a Gaussian shape.
            If beta=1, K resembles a Lorentzian shape.
        """
        print('Brown')
        return amp / (1 + (beta * (x - x0) / (sqrt(2) * hwhm)) ** 2) ** (1 / beta ** 2)

    def ModelByBand(self, x0, amp, hwhm, beta):
        """LorentzianByBand.

        Parameters
        ----------
        x0 :
            x0
        amp :
            amp
        hwhm :
            hwhm
        """
        return pd.Series([self.K(x, x0, amp, hwhm, beta) for x in self.df['Wavenumber']])

    def FitModel(self, x, *args):
        """FitModel.

        Args:
            x:
            args:
        """
        n = self.bands.shape[0]
        x0 = args[:n]
        amp = args[n:2 * n]
        hwhm = args[2 * n:3 * n]
        beta = args[3 * n:]

        res = [self.K(x, x0[i], amp[i], hwhm[i], beta[i]) for i in range(n)]

        return sum(res)

