#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from Models import StatusModel
import multiprocessing
from psutil import virtual_memory

calcWindowBase, calcWindowForm = uic.loadUiType('Views/uiCalculationsWindow2.ui')
geomWidgetBase, geomWidgetForm = uic.loadUiType('Views/uiGeometryWidget.ui')
freqWidgetBase, freqWidgetForm = uic.loadUiType('Views/uiFreqWidget.ui')
ircWidgetBase, ircWidgetForm = uic.loadUiType('Views/uiIRCWidget.ui')


class CalculationsController(calcWindowBase, calcWindowForm):

    def __init__(self, parent=None):
        super(calcWindowBase, self).__init__(parent)
        self.setupUi(self)

        # mem = virtual_memory()
        self.memory = virtual_memory().total//1000000
        self.cpu = multiprocessing.cpu_count()
        self._jobTypes = ('Energy', 'Optimization', 'Frequency', 'Opt+Freq', 'IRC',
                          'Scan', 'Stability', 'NMR')
        self._optimizeToA = ('Minimum', 'TS (Berny)')
        self._forceConstants = (None, 'Calculate at First Point', 'Calculate at all Point',
                                'Read Internal', 'Read Cartesian')
        self._computeRaman = (None, 'No', 'Yes')
        self._saveNormalModesOptions = (None, 'Yes', 'No')
        self._computeROA = ('No', 'Yes')
        self._readIncidentLightFreqs = (None, 'Yes', 'No')
        self._state = ('Ground State', 'ZINDO', 'CIS', 'SAC-CI', 'TD-SCF', 'TDA', 'EOM-CCSD')
        self._methods = ('Mechanics...', 'Semi-empirical...', 'Hartree-Fock', 'DFT...',
                         'MP2', 'MP4', 'CCSD', 'BD', 'CASSCF')
        self._spin = ('Default Spin', 'Restricted', 'Unrestricted', 'Restricted-Open')
        self._semiempiricalFunctionals = ('PM6', 'PDDG', 'AM1', 'PM3', 'PM3MM', 'INDO', 'CNDO')
        self._dftFunctionals = ('LSDA', 'BPV86', 'B3LYP', 'CAM-B3LYP', 'B3PW91', 'MPW1PW91',
                                'PBEPBE', 'HSEH1PBE', 'HCTH', 'TPSSTPSS', 'WB97XD', 'APFD',
                                'MPWB1K', 'M06-2X')
        self._basisSet = ('STO-3G', '3-21G', '6-31G', '6-31G\'', '6-311G', 'cc-pVDZ', 'cc-pVTZ',
                          'cc-pVQZ', 'LanL2DZ', 'LanL2MB', 'SDD', 'DGDZVP', 'DGDZVP2', 'DGTZVP',
                          'GEN', 'GENECP')

        self.chosenParameters = {
            'METHOD': None, 'FUNCTIONAL': None, 'BASIS': None,
            'BASIS2': None, 'SAVE FILE': False, 'STATUS': 'Pending'
            }
        self._methodWidgets = [
            self.uiStateComboBox,
            self.uiMethodComboBox,
            self.uiShellTypeComboBox,
            self.uiFunctionalComboBox,
            self.uiBasisComboBox,
            self.uiDiffuseComboBox,
            self.uiPolarizedComboBox,
            self.uiSecondPolarizedComboBox,
            self.uiBasis2ComboBox,
            self.uiDiffuse2ComboBox,
            self.uiPolarized2ComboBox,
            self.uiSecondPolarizedComboBox,
            self.uiChargeSpinBox,
            self.uiSpinComboBox
            ]

        # por aquí
        self._widgets[0].addItems(self._methods)
        self._widgets[1].addItems(self._dftFunctionals)
        self._widgets[2].addItems(self._semiempiricalFunctionals)
        self._widgets[3].addItems(self._basisSetsLightAtoms)
        self._widgets[4].addItems(self._basisSetsHeavyAtoms)
        for w in self._widgets[0:5]:
            w.setCurrentIndex(0)
        if not self._heavy:
            self._widgets[4].setEnabled(False)
        self._widgets[5].setEnabled(False)
        self._widgets[6].setEnabled(False)

    # --------------------------------------------------SIGNALS--------------------------------------------------
        self._widgets[0].currentIndexChanged.connect(lambda: self.setChosenParameters(1))
        self._widgets[1].currentIndexChanged.connect(lambda: self.setChosenParameters(2))
        self._widgets[2].currentIndexChanged.connect(lambda: self.setChosenParameters(3))
        self._widgets[3].currentIndexChanged.connect(lambda: self.setChosenParameters(4))
        self._widgets[4].currentIndexChanged.connect(lambda: self.setChosenParameters(5))
        self._widgets[5].stateChanged.connect(lambda: self.setChosenParameters(6))
        self._widgets[6].pressed.connect(lambda: self.SetQueue('opt'))  # debe cambiar para poner en cola
        self._widgets[6].pressed.connect(window.QueueManager)
    # --------------------------------------------------METHODS--------------------------------------------------

    def setChosenParameters(self, number):

        if number == 1:
            self.chosenParameters['METHOD'] = self._methods[
                self._widgets[0].currentIndex()
                ]
            if self._widgets[0].currentIndex() == 1:
                self._widgets[1].setEnabled(True)
                self._widgets[2].setEnabled(False)
                self._widgets[3].setEnabled(True)
                if self._heavy:
                    self._widgets[4].setEnabled(True)
                self._widgets[2].setCurrentIndex(0)
            elif self._widgets[0].currentIndex() == 2:
                self._widgets[1].setEnabled(False)
                self._widgets[3].setEnabled(False)
                self._widgets[4].setEnabled(False)
                self._widgets[2].setEnabled(True)
                self._widgets[1].setCurrentIndex(0)
                self._widgets[3].setCurrentIndex(0)
                self._widgets[4].setCurrentIndex(0)
            else:
                for w in self._widgets[1:5]:
                    w.setEnabled(False)
                    w.setCurrentIndex(0)
        elif number == 2:
            self.chosenParameters['FUNCTIONAL'] = self._dftFunctionals[
                self._widgets[1].currentIndex()
                ]
        elif number == 3:
            self.chosenParameters['FUNCTIONAL'] = self._semiempiricalFunctionals[
                self._widgets[2].currentIndex()
                ]
            if self._widgets[2].currentIndex() > 0 and self._widgets[0].currentIndex() > 0:
                self._widgets[5].setEnabled(True)
                self._widgets[6].setEnabled(True)
            else:
                self._widgets[5].setEnabled(False)
                self._widgets[6].setEnabled(False)
        elif number == 4:
            self.chosenParameters['BASIS'] = self._basisSetsLightAtoms[
                self._widgets[3].currentIndex()
                ]
            if self._widgets[3].currentIndex() > 0 and self._widgets[0].currentIndex() > 0:
                self._widgets[5].setEnabled(True)
                self._widgets[6].setEnabled(True)
            else:
                self._widgets[5].setEnabled(False)
                self._widgets[6].setEnabled(False)
        elif number == 5:
            self.chosenParameters['BASIS2'] = self._basisSetsHeavyAtoms[
                self._widgets[4].currentIndex()
                ]
        elif number == 6:
            if self._widgets[5].isChecked():
                self.chosenParameters['SAVE FILE'] = True
            else:
                self.chosenParameters['SAVE FILE'] = False

    def SetQueue(self, kind):

        self.queue.append([
            kind,
            self.chosenParameters['METHOD'],
            self.chosenParameters['FUNCTIONAL'],
            self.chosenParameters['BASIS'],
            self.chosenParameters['BASIS2'],
            self.chosenParameters['SAVE FILE'],
            self.chosenParameters['STATUS'],
            self.memory,
            self.cpu
            ])
