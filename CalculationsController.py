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
        # ------------------------------------PROPERTIES------------------------------------------------------

        # mem = virtual_memory()
        self.memory = virtual_memory().total // 1000000
        self.cpu = multiprocessing.cpu_count()
        # self._parametersList = []
        self._link0MemoryLine = ''
        self._link0CPULine = ''
        self._link0ChkLine = ''
        self._methodInputLine = []
        self._nameLine = ''
        self._chargeMultiplicityLine = []
        self._addInputs = []

        # ------------------------------------WIDGETS------------------------------------------------------
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
        self._mechanics = ('UFF', 'Dreiding', 'Amber')
        self._augmented = (None, 'aug-')
        self._basisSet = ('STO-3G', '3-21G', '6-31G', '6-31G\'', '6-311G', 'cc-pVDZ', 'cc-pVTZ',
                          'cc-pVQZ', 'LanL2DZ', 'LanL2MB', 'SDD', 'DGDZVP', 'DGDZVP2', 'DGTZVP',
                          'GEN', 'GENECP')
        self._asterisk = (None, '*', '**')
        self._diffuse = (None, '+', '++')
        self._firstPolarizedOrbitals = (None, 'd', '2d', '3d', 'df', '2df', '3df', '3d2f')
        self._secondPolarizedOrbitals = (None, 'p', '2p', '3p', 'pd', '2pd', '3pd', '3p2d')
        self.chosenParameters = {
            'METHOD': None, 'FUNCTIONAL': None, 'BASIS': None,
            'BASIS2': None, 'SAVE FILE': False, 'STATUS': 'Pending'
        }
        self._jobsWidgets = [
            self.uiJobTypeComboBox
        ]
        self._methodWidgets = [
            self.uiStateComboBox,  # 0
            self.uiMethodComboBox,  # 1
            self.uiShellTypeComboBox,  # 2
            self.uiDftFunctionalComboBox,  # 3
            self.uiSemiFunctionalComboBox,  # 4
            self.uiMechanicsComboBox,  # 5
            self.uiBasisLabel,  # 6
            self.uiAugmentedComboBox,  # 7
            self.uiBasisComboBox,  # 8
            self.uiDiffuseComboBox,  # 9
            self.uiAsteriskComboBox,  # 10
            self.uiParentheses1,  # 11
            self.uiPolarizedComboBox,  # 12
            self.uiComma1,  # 13
            self.uiSecondPolarizedComboBox,  # 14
            self.uiParentheses2,  # 15
            self.uiBasis2Label,  # 16
            self.uiAugmented2ComboBox,  # 17
            self.uiBasis2ComboBox,  # 18
            self.uiDiffuse2ComboBox,  # 19
            self.uiAsterisk2ComboBox,  # 20
            self.uiParentheses3,  # 21
            self.uiPolarized2ComboBox,  # 22
            self.uiComma2,  # 23
            self.uiSecondPolarized2ComboBox,  # 24
            self.uiParentheses4,  # 25
            self.uiChargeSpinBox,  # 26
            self.uiEvenSpinComboBox,  # 27
            self.uiOddSpinComboBox,  # 28
        ]

        # ------------------------------------POPULATE WIDGETS------------------------------------------------------
        self._jobsWidgets[0].addItems(self._jobTypes)
        self._methodWidgets[0].addItems(self._state)
        self._methodWidgets[1].addItems(self._methods)
        self._methodWidgets[1].setCurrentIndex(2)
        self._methodWidgets[2].addItems(self._spin)
        self._methodWidgets[3].addItems(self._dftFunctionals)
        self._methodWidgets[4].addItems(self._semiempiricalFunctionals)
        self._methodWidgets[5].addItems(self._mechanics)
        self._methodWidgets[7].addItems(self._augmented)
        self._methodWidgets[8].addItems(self._basisSet)
        self._methodWidgets[9].addItems(self._diffuse)
        self._methodWidgets[10].addItems(self._asterisk)
        self._methodWidgets[12].addItems(self._firstPolarizedOrbitals)
        self._methodWidgets[14].addItems(self._secondPolarizedOrbitals)
        self._methodWidgets[17].addItems(self._augmented)
        self._methodWidgets[18].addItems(self._basisSet)
        self._methodWidgets[19].addItems(self._diffuse)
        self._methodWidgets[20].addItems(self._asterisk)
        self._methodWidgets[22].addItems(self._firstPolarizedOrbitals)
        self._methodWidgets[24].addItems(self._secondPolarizedOrbitals)

        [self._methodWidgets[i].setVisible(False) for i in range(3, 26) if i != 8 and i != 10]
        # ------------------------------------SIGNALS---------------------------------------------------------------
        self._methodWidgets[1].currentIndexChanged.connect(lambda: self.MethodController(1, self._methodWidgets[1].currentIndex()))
        # ------------------------------------METHODS---------------------------------------------------------------

    def MethodController(self, widgetNumber, selection):

        if widgetNumber == 1:
            if selection == 0:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 11) if i not in (5, 7, 9)]
                self._methodWidgets[5].setVisible(True)
            elif selection == 1:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 11) if i not in (4, 7, 9)]
                self._methodWidgets[4].setVisible(True)
            elif selection == 2:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 6)]
                [self._methodWidgets[i].setVisible(True) for i in range(6, 11) if i not in (7, 9)]
            elif selection == 3:
                [self._methodWidgets[i].setVisible(False) for i in range(4, 6)]
                [self._methodWidgets[i].setVisible(True) for i in range(3, 11) if i not in (4, 5, 7, 9)]

        self._parametersList = [self._methodWidgets[i].currentIndex() for i in range(0, 25) if i not in (6, 11, 13, 15, 16, 21, 23)]
        # print(self._parametersList)
        # por aquí
