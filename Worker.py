#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from PyQt5 import QtCore as qtc


class Worker(qtc.QObject):

    finished = qtc.pyqtSignal()
    workflow = qtc.pyqtSignal(int, str)

    def __init__(self, queue):
        super().__init__()
        self.master_queue = queue
        self.paused = False

    @qtc.pyqtSlot()
    def start_queue(self):
        for job in self.master_queue:
            try:
                if self.paused:
                    break
                if job.type == 'Docking':
                    done = self.check_vina(job.output_file)
                elif job.type in ('Optimization', 'Energy', 'Frequency'):
                    done = self.check_gauss(job.output_file)
                if done:
                    self.workflow.emit(job.id, 'Finished')
                else:
                    self.workflow.emit(job.id, 'Running')
                    if job.type == 'Docking':
                        self.receptor_prep(
                            job.keywords.get('receptor_in'),
                            job.keywords.get('receptor_out')
                        )
                        self.run_vina(job.keywords)
                    elif job.type in ('Optimization', 'Energy', 'Frequency'):
                        self.run_gauss(job.input_file, job.output_file)
                    self.workflow.emit(job.id, 'Finished')
            except:
                self.workflow.emit(job.id, 'Failed')
        self.finished.emit()
        self.deleteLater()

    @qtc.pyqtSlot()
    def pause(self):
        self.paused = True

    def run_gauss(self, input_file, output_file):
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
        cmd = f'{gauss_exec} < {input_file} > {output_file}'
        subprocess.run(cmd, shell=True)
        return True

    def check_gauss(self, output_file):
        if not os.path.exists(output_file):
            return False
        if 'Optimization' in output_file:
            with open(output_file, 'r') as f:
                file = f.read()
                finished = re.findall('Optimization completed.', file)
                if finished:
                    return True
                return False

    def run_vina(self, config):
        vina = Vina(sf_name='vina')
        # primero hay que preparar el receptor para convertirlo en pdbqt
        vina.set_receptor(config.get('receptor_out'))
        vina.set_ligand_from_file(config.get('ligand'))
        vina.compute_vina_maps(
            center=config.get('center'),
            box_size=config.get('size')
        )
        vina.optimize()
        vina.dock(
            exhaustiveness=config.get('exhaustiveness'),
            n_poses=20  # se evalúan 20 aunque solo se escriben las que diga el usuario
        )
        vina.write_poses(
            config.get('out'),
            n_poses=config.get('num_modes'),
            energy_range=config.get('energy_range'),
            overwrite=True
        )
        return True

    def check_vina(self, output_file):
        return os.path.exists(output_file)

