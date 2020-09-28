#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from Models import StatusModel
import multiprocessing
from psutil import virtual_memory

class CalculationsController(qtw.QWidget):

    def __init__(self):
        super(CalculationsController, self).__init__()
        uic.loadUi('Views/uiCalculationsWindow2.ui', self)
        # uic.loadUi('Views/uiGeometryWidget.ui', self)
        # uic.loadUi('Views/uiFreqWidget.ui', self)
        # uic.loadUi('Views/uiIRCWidget.ui', self)

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
