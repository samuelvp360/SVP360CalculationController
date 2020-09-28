#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc


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
                self.index = index
                self.molecule, self.parameters = calc
                if self.parameters[6] == 'Pending':
                    self.startingCalc.emit(self.molecule, self.index, 'Running')
                    self.molecule.SetCalculations(*self.parameters)
                    self.okCalc.emit(self.molecule, self.index, 'Finished')
                if self.isStopped:
                    break
        except:
            self.error.emit(self.molecule, self.index, 'Aborted')
        finally:
            self.finished.emit()
            self.deleteLater()

    def Stop(self):
        self.isStopped = True

    def Pause(self):
        self.isPaused = True
