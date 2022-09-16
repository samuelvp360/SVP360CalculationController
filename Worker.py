#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
import multiprocessing
from PyQt5 import QtCore as qtc
from vina import Vina
from loguru import logger
from scipy.optimize import curve_fit
from numpy import exp, sqrt
import pandas as pd
import copy


class Worker(qtc.QObject):

    finished = qtc.pyqtSignal()
    workflow = qtc.pyqtSignal(int, str)

    def __init__(self, queue):
        super().__init__(parent=None)
        self.master_queue = copy.deepcopy(queue)
        self.paused = False

    @qtc.pyqtSlot()
    def start_queue(self):
        for job in self.master_queue:
            self.job = job
            try:
                if self.paused:
                    break
                if job.type == 'Docking':
                    done = self.check_vina()
                elif job.type in ('Optimization', 'Energy', 'Frequency'):
                    done = self.check_gauss()
                if done:
                    self.workflow.emit(job.id, 'Finished')
                else:
                    self.workflow.emit(job.id, 'Running')
                    if job.type == 'Docking':
                        self.run_vina()
                    elif job.type in ('Optimization', 'Energy', 'Frequency'):
                        self.run_gauss()
                    self.workflow.emit(job.id, 'Finished')
            except:
                self.workflow.emit(job.id, 'Failed')
        self.finished.emit()
        self.deleteLater()

    @qtc.pyqtSlot()
    def pause(self):
        self.paused = True

    def run_gauss(self):
        '''
        Parameters
        ----------
        input_file: the .com file as Gaussian input
        output_file: the name of the .log file to be done

        Returns
        -------
        True if the input file is successfully run in the shell
        False if something went wrong with the gaussian executable
        '''
        env = os.environ.copy()
        try:
            gauss_exec = env['GAUSS_EXEDIR'].split('/')[-1]
        except KeyError:
            return False
            print('no leyó el ejecutable de Gaussian')
        cmd = f'{gauss_exec} < {self.job.input_file} > {self.job.output_file}'
        subprocess.run(cmd, shell=True)
        return True

    def check_gauss(self):
        if not os.path.exists(self.job.output_file):
            return False
        if 'Optimization' in self.job.output_file:
            with open(output_file, 'r') as f:
                file = f.read()
                finished = re.findall('Optimization completed.', file)
                if finished:
                    return True
                return False

    # @logger.catch
    def run_vina(self):
        cmd = 'vina --verbosity 0'
        for k, v in self.job.config.items():
            cmd += f' --{k} "{v}"'
        for out in self.job.output_file:
            new_cmd = cmd + f' --out "{out}"'
            subprocess.run(
                new_cmd,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )

    def check_vina(self):
        return all((os.path.exists(out) for out in self.job.output_file))


class IRWorkerThread(qtc.QObject):

    finished = qtc.pyqtSignal()
    ok_fit = qtc.pyqtSignal(dict, object, object)
    error = qtc.pyqtSignal()

    def __init__(self, bands, df):
        super().__init__()
        self.bands = bands
        self.predicted_bands = {}
        self.df = df

    @qtc.pyqtSlot()
    def fit(self):
        n = self.bands.shape[0]
        position = self.bands['Wavenumber'].to_list()
        amplitude = self.bands['Absorbance'].to_list()
        hwhm = self.bands['HWHM'].to_list()
        beta = n * [0.5]
        parameters = position + amplitude + hwhm + beta
        lower_bounds = tuple([i - 50 for i in position] + [i / 2.0 for i in amplitude] + n * [1e-5] + n * [0])
        upper_bounds = tuple([i + 50 for i in position] + [i * 1.5 for i in amplitude] + n * [1000] + n * [1])

        try:
            fitted_params, cov = curve_fit(
                self.fit_model, self.df['Wavenumber'],
                self.df['Absorbance'], p0=parameters, absolute_sigma=True,
                bounds=(lower_bounds, upper_bounds)
            )
            fitted_position = fitted_params[:n]
            fitted_amplitude = fitted_params[n:2 * n]
            fitted_hwhm = fitted_params[2 * n:3 * n]
            fitted_beta = fitted_params[3 * n:]
            self.bands['Wavenumber'] = pd.Series(fitted_position)
            self.bands['Absorbance'] = pd.Series(fitted_amplitude)
            self.bands['HWHM'] = pd.Series(fitted_hwhm)
            self.bands['Beta'] = pd.Series(fitted_beta)
            self.predicted_bands['Absorbance'] = pd.DataFrame(
                map(
                    self.model_by_band, self.bands['Wavenumber'],
                    self.bands['Absorbance'], self.bands['HWHM'],
                    self.bands['Beta']
                )
            ).transpose()
            self.ok_fit.emit(self.predicted_bands, self.bands, self.df)
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

    def model_by_band(self, x0, amp, hwhm, beta):
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

    def fit_model(self, x, *args):
        """fit_model.

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


class MolWorker(qtc.QObject):

    finished = qtc.pyqtSignal()
    workflow = qtc.pyqtSignal(float, str, object)

    def __init__(self, mol_list):
        super().__init__(parent=None)
        self.mol_list = copy.deepcopy(mol_list)

    @qtc.pyqtSlot()
    def start(self):
        for mol in self.mol_list:
            Rg, conf = mol.calculate_Rg()
            self.workflow.emit(Rg, mol.inchi_key, conf)
        self.finished.emit()
        self.deleteLater()


