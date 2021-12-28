#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from Calculations import GaussianWorker

class Worker(qtc.QObject):

    finished = qtc.pyqtSignal()
    workflow = qtc.pyqtSignal(int, str)

    def __init__(self, queue):
        super().__init__()
        self.gaussian_worker = GaussianWorker()
        self.master_queue = queue
        self.paused = False

    @qtc.pyqtSlot()
    def start_queue(self):
        for job in self.master_queue:
            try:
                if self.paused:
                    break
                done = self.gaussian_worker.check(job.output_file)
                if done:
                    self.workflow.emit(job.id, 'Finished')
                else:
                    self.workflow.emit(job.id, 'Running')
                    self.gaussian_worker.run(job.input_file, job.output_file)
                    self.workflow.emit(job.id, 'Finished')
            except:
                self.workflow.emit(job.id, 'Failed')
        self.finished.emit()
        self.deleteLater()

    @qtc.pyqtSlot()
    def pause(self):
        self.paused = True

