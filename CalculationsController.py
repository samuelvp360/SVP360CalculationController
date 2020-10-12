#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
# from PyQt5 import QtGui as qtg
from PyQt5 import uic
# from Models import StatusModel
import multiprocessing
from psutil import virtual_memory


class CalculationsController(qtw.QWidget):

    def __init__(self):
        super(CalculationsController, self).__init__()
        uic.loadUi('Views/uiCalculationsWindow2.ui', self)
        self.uiOptGroupBox.setVisible(False)
        self.uiFreqGroupBox.setVisible(False)
        self.uiIRCGroupBox.setVisible(False)
        # ------------------------------------PROPERTIES------------------------------------------------------

        # mem = virtual_memory()
        self.memory = virtual_memory().total // 1000000
        self.cpu = multiprocessing.cpu_count()
        # self._parametersList = []
        self._keywordsLine = ['# ', '', '', '', '', '', '']
        self._link0MemoryLine = ''
        self._link0CPULine = ''
        self._link0ChkLine = ''
        self._chargeMultiplicityLine = ['', ' ', '']
        self._addInputs = []

        # ------------------------------------WIDGETS------------------------------------------------------
        self._jobTypes = ('Energy', 'Optimization', 'Frequency', 'Opt+Freq', 'IRC',
                          'Scan', 'Stability', 'NMR')
        self._optimizeToA = ('Minimum', 'TS (Berny)')
        self._forceConstants = (None, 'Calculate at First Point', 'Calculate at all Points',
                                'Read Internal', 'Read Cartesian')
        self._computeRaman = (None, 'No', 'Yes')
        self._saveNormalModesOptions = (None, 'Yes', 'No')
        self._computeROA = ('No', 'Yes')
        self._readIncidentLightFreqs = (None, 'Yes', 'No')
        self._forceConstants2 = ('Calculate Once', 'Calculate Always', 'Read from .CHK')
        self._recorrect = (None, 'Never', 'Yes', 'Always', 'Test')
        self._followIRC = ('Both directions', 'Forward only', 'Reverse only')
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
            self.uiJobTypeComboBox,  # 0
            self.uiTightCheckBox,  # 1
            self.uiOptimizeToAComboBox,  # 2
            self.uiForceConstComboBox,  # 3
            self.uiUseRFOStepsCheckBox,  # 4
            self.uiAnharmonicCheckBox,  # 5
            self.uiAnharmonicLine,  # 6
            self.uiSkipCheckBox,  # 7
            self.uiROAComboBox,  # 8
            self.uiProjectedFreqCheckBox,  # 9
            self.uiNormalModeComboBox,  # 10
            self.uiAnharmonicModesCheckBox,  # 11
            self.uiIncidentLightComboBox,  # 12
            self.uiRamanComboBox,  # 13
            self.uiVCDCheckBox,  # 14
            self.uiRecalculateCheckBox,  # 15
            self.uiForceConst2ComboBox,  # 16
            self.uiIRCMaxCheckBox,  # 17
            self.uiMorePointsLine,  # 18
            self.uiMorePointsCheckBox,  # 19
            self.uiRecorrectComboBox,  # 20
            self.uiRecalculateLine,  # 21
            self.uiFollowIRCComboBox  # 22
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
            self.uiEvenSpinSpinBox,  # 27
            self.uiOddSpinSpinBox,  # 28
        ]

        # ------------------------------------POPULATE WIDGETS------------------------------------------------------
        self._methodWidgets[0].addItems(self._state)
        self._methodWidgets[1].addItems(self._methods)
        self._methodWidgets[1].setCurrentIndex(2)
        self._methodWidgets[2].addItems(self._spin)
        self._methodWidgets[3].addItems(self._dftFunctionals)
        self._methodWidgets[3].setCurrentIndex(2)
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
        [self._methodWidgets[i].setVisible(False) for i in range(3, 26) if i not in (6, 8, 10)]

        self._jobsWidgets[0].addItems(self._jobTypes)
        self._jobsWidgets[2].addItems(self._optimizeToA)
        self._jobsWidgets[3].addItems(self._forceConstants)
        self._jobsWidgets[8].addItems(self._computeROA)
        self._jobsWidgets[10].addItems(self._saveNormalModesOptions)
        self._jobsWidgets[12].addItems(self._readIncidentLightFreqs)
        self._jobsWidgets[13].addItems(self._computeRaman)
        self._jobsWidgets[16].addItems(self._forceConstants2)
        self._jobsWidgets[20].addItems(self._recorrect)
        self._jobsWidgets[22].addItems(self._followIRC)

        self.MethodController(1, self._methodWidgets[1].currentIndex())

        # ------------------------------------SIGNALS---------------------------------------------------------------
        self._methodWidgets[1].currentIndexChanged.connect(lambda: self.MethodController(1, self._methodWidgets[1].currentIndex()))
        self._jobsWidgets[0].currentIndexChanged.connect(lambda: self.JobTypeController(0, self._jobsWidgets[0].currentIndex()))
        self.uiTitleLineEdit.textChanged.connect(self.SetTitle)
        self._methodWidgets[26].valueChanged.connect(self.SetCahrgeMult)
        self._methodWidgets[27].valueChanged.connect(self.SetCahrgeMult)
        self._methodWidgets[28].valueChanged.connect(self.SetCahrgeMult)
        # ------------------------------------METHODS---------------------------------------------------------------

    def DetectMolecule(self, molecule):
        """
        docstring
        """
        self._molecule = molecule
        self.uiTitleLineEdit.setText(self._molecule.GetName)
        self._methodWidgets[26].setValue(self._molecule.GetNetCharge)
        self.SetCahrgeMult()

    def MethodController(self, widgetNumber, selection):

        if widgetNumber == 1:
            if selection == 0:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 11) if i not in (5, 7, 9)]
                self._methodWidgets[5].setVisible(True)
            elif selection == 1:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 11) if i not in (4, 7, 9)]
                self._methodWidgets[4].setVisible(True)
                functional = self._semiempiricalFunctionals[self._methodWidgets[4].currentIndex()]
                self._keywordsLine[1] = f'{functional}'
                self._keywordsLine[2] = ''
            elif selection == 2:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 6)]
                [self._methodWidgets[i].setVisible(True) for i in range(6, 11) if i not in (7, 9)]
                self._keywordsLine[1] = 'hf/'
                self._keywordsLine[2] = self._basisSet[self._methodWidgets[8].currentIndex()]
            elif selection == 3:
                [self._methodWidgets[i].setVisible(False) for i in range(4, 6)]
                [self._methodWidgets[i].setVisible(True) for i in range(3, 11) if i not in (4, 5, 7, 9)]
                functional = self._dftFunctionals[self._methodWidgets[3].currentIndex()]
                self._keywordsLine[1] = f'{functional}/'
                self._keywordsLine[2] = self._basisSet[self._methodWidgets[8].currentIndex()]

        self.SetKeywords()
        self._parametersList = [self._methodWidgets[i].currentIndex() for i in range(0, 25) if i not in (6, 11, 13, 15, 16, 21, 23)]
        # print(self._parametersList)

    def JobTypeController(self, widgetNumber, selection):

        if widgetNumber == 0:
            if selection == 0:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(False)
            elif selection == 1:
                self.uiOptGroupBox.setVisible(True)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(False)
            elif selection == 2:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(True)
                self.uiIRCGroupBox.setVisible(False)
            elif selection == 3:
                self.uiOptGroupBox.setVisible(True)
                self.uiFreqGroupBox.setVisible(True)
                self.uiIRCGroupBox.setVisible(False)
            elif selection == 4:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(True)
            else:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(False)

    def SetTitle(self):

        self.uiTitleLabel.setText(self.uiTitleLineEdit.text())

    def SetKeywords(self):

        self.uiKeywordsLabel.setText(''.join(self._keywordsLine))

    def SetCahrgeMult(self):  # mejorar
        """
        docstring
        """
        self._chargeMultiplicityLine[0] = str(self._methodWidgets[26].value())

        combination = (self._molecule.GetValenceElectrons - self._methodWidgets[26].value() + self._molecule.GetNetCharge) % 2
        if combination == 0:
            self._methodWidgets[28].setVisible(True)
            self._methodWidgets[27].setVisible(False)
            self._methodWidgets[27].setValue(2)
            self._chargeMultiplicityLine[2] = str(self._methodWidgets[28].value())
        else:
            self._methodWidgets[28].setVisible(False)
            self._methodWidgets[27].setVisible(True)
            self._methodWidgets[28].setValue(1)
            self._chargeMultiplicityLine[2] = str(self._methodWidgets[27].value())
        
        self.uiChargeMultLabel.setText(''.join(self._chargeMultiplicityLine))        
# por aquí