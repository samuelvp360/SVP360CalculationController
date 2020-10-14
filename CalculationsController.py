#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
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
        self._keywordsLine = []
        for i in range(15):
            self._keywordsLine.append('')
        self._keywordsLine[0] = '# '
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
        self._shellType = ('Default Spin', 'Restricted', 'Unrestricted', 'Restricted-Open')
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
        self._fittingSet = (
            None, 'Default', 'AUTO', 'AUTO=N', 'AUTO=ALL', 'DEF2SV', 'DEF2TZV',
            'QZV', 'DGA1', 'DGA2', 'ASV', 'ATZ'
        )
        # self.chosenParameters = {
        #     'METHOD': None, 'FUNCTIONAL': None, 'BASIS': None,
        #     'BASIS2': None, 'SAVE FILE': False, 'STATUS': 'Pending'
        # }
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
            self.uiFittingSetLabel,  # 26
            self.uiFittingSetComboBox,  # 27
            self.uiFittingSet2Label,  # 28
            self.uiFittingSet2ComboBox,  # 29
            self.uiChargeSpinBox,  # 30
            self.uiEvenSpinSpinBox,  # 31
            self.uiOddSpinSpinBox,  # 32
            self.uiSparseCheckBox,  # 33
            self.uiAutoValueSpinBox  # 34
        ]

        # ------------------------------------POPULATE WIDGETS------------------------------------------------------
        self._methodWidgets[0].addItems(self._state)
        self._methodWidgets[1].addItems(self._methods)
        self._methodWidgets[1].setCurrentIndex(2)
        self._methodWidgets[2].addItems(self._shellType)
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
        self._methodWidgets[27].addItems(self._fittingSet)
        self._methodWidgets[29].addItems(self._fittingSet)
        [self._methodWidgets[i].setVisible(False) for i in range(3, 35) if i not in (6, 8, 10, 30, 31, 32)]

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

        # ------------------------------------SIGNALS---------------------------------------------------------------
        self._methodWidgets[1].currentIndexChanged.connect(lambda: self.SetKeywords(1, self._methodWidgets[1].currentIndex()))
        self._methodWidgets[2].currentIndexChanged.connect(lambda: self.SetKeywords(2, self._methodWidgets[2].currentIndex()))
        self._methodWidgets[3].currentIndexChanged.connect(lambda: self.SetKeywords(3, self._methodWidgets[3].currentIndex()))
        self._methodWidgets[4].currentIndexChanged.connect(lambda: self.SetKeywords(4, self._methodWidgets[4].currentIndex()))
        self._methodWidgets[5].currentIndexChanged.connect(lambda: self.SetKeywords(5, self._methodWidgets[5].currentIndex()))
        self._methodWidgets[7].currentIndexChanged.connect(lambda: self.SetKeywords(7, self._methodWidgets[7].currentIndex()))
        self._methodWidgets[8].currentIndexChanged.connect(lambda: self.SetKeywords(8, self._methodWidgets[8].currentIndex()))
        self._methodWidgets[9].currentIndexChanged.connect(lambda: self.SetKeywords(9, self._methodWidgets[9].currentIndex()))
        self._methodWidgets[10].currentIndexChanged.connect(lambda: self.SetKeywords(10, self._methodWidgets[10].currentIndex()))
        self._methodWidgets[12].currentIndexChanged.connect(lambda: self.SetKeywords(12, self._methodWidgets[12].currentIndex()))
        self._methodWidgets[14].currentIndexChanged.connect(lambda: self.SetKeywords(14, self._methodWidgets[14].currentIndex()))
        self._methodWidgets[27].currentIndexChanged.connect(lambda: self.SetKeywords(27, self._methodWidgets[27].currentIndex()))
        self._methodWidgets[33].stateChanged.connect(lambda: self.SetKeywords(33, self._methodWidgets[33].isChecked()))
        self._methodWidgets[34].valueChanged.connect(lambda: self.SetKeywords(34, self._methodWidgets[34].value()))
        self._jobsWidgets[0].currentIndexChanged.connect(lambda: self.JobTypeController(0, self._jobsWidgets[0].currentIndex()))
        self.uiTitleLineEdit.textChanged.connect(self.SetTitle)
        self._methodWidgets[30].valueChanged.connect(self.SetCahrgeMult)
        self._methodWidgets[31].valueChanged.connect(self.SetCahrgeMult)
        self._methodWidgets[32].valueChanged.connect(self.SetCahrgeMult)
        
        # ------------------------------------METHODS---------------------------------------------------------------

    def DetectMolecule(self, molecule):
        """
        docstring
        """
        self._molecule = molecule
        self.uiTitleLineEdit.setText(self._molecule.GetName)
        self._methodWidgets[30].setValue(self._molecule.GetNetCharge)
        self._methodWidgets[10].model().item(2).setEnabled(False)
        self.SetKeywords(1, self._methodWidgets[1].currentIndex())
        self.SetCahrgeMult()

    def SetKeywords(self, widgetNumber, selection):

        if widgetNumber == 1:
            if selection == 0:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 34) if i not in (5, 7, 30, 31, 32)]
                self._methodWidgets[5].setVisible(True)
                method = self._mechanics[self._methodWidgets[5].currentIndex()]
                self._keywordsLine[2] = f'{method}'
                for i in (3, 4, 9, 10, 11):
                    self._keywordsLine[i] = ''
            elif selection == 1:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 34) if i not in (4, 7, 30, 31, 32)]
                self._methodWidgets[4].setVisible(True)
                functional = self._semiempiricalFunctionals[self._methodWidgets[4].currentIndex()]
                self._keywordsLine[2] = f'{functional}'
                for i in (3, 4, 9, 10, 11):
                    self._keywordsLine[i] = ''
            elif selection == 2:
                [self._methodWidgets[i].setVisible(False) for i in range(3, 35) if i not in (6, 8, 10, 30, 31, 32)]
                [self._methodWidgets[i].setVisible(True) for i in (6, 8, 10)]
                self._methodWidgets[8].setCurrentIndex(0)
                self._methodWidgets[10].setCurrentIndex(0)
                self._keywordsLine[2] = 'HF/'
                self._keywordsLine[4] = self._basisSet[self._methodWidgets[8].currentIndex()]
                for i in range(9, 12):
                    self._keywordsLine[i] = ''
            elif selection == 3:
                [self._methodWidgets[i].setVisible(False) for i in range(4, 6)]
                [self._methodWidgets[i].setVisible(True) for i in (3, 6, 8, 33)]
                self._methodWidgets[33].setEnabled(False)
                self._methodWidgets[3].setCurrentIndex(2)
                functional = self._dftFunctionals[self._methodWidgets[3].currentIndex()]
                self._keywordsLine[2] = f'{functional}/'
                self._keywordsLine[4] = self._basisSet[self._methodWidgets[8].currentIndex()]
        elif widgetNumber == 2:
            if selection == 0:
                self._keywordsLine[1] = ''
            elif selection == 1:
                self._keywordsLine[1] = 'R'
            elif selection == 2:
                self._keywordsLine[1] = 'U'
            else:
                self._keywordsLine[1] = 'RO'
        elif widgetNumber == 3:  # faltan los últimos dos, faltan los del fitting set
            functional = self._dftFunctionals[selection]
            self._keywordsLine[2] = f'{functional}/'
            if selection in (0, 1, 6, 8, 9):
                [self._methodWidgets[i].setVisible(True) for i in (26, 27)]
                self._methodWidgets[33].setEnabled(True)
                self.SetKeywords(33, self._methodWidgets[33].isChecked())
            else:
                [self._methodWidgets[i].setVisible(False) for i in (26, 27, 34)]
                self._methodWidgets[33].setEnabled(False)
                self._methodWidgets[27].setCurrentIndex(0)
                for i in range(9, 12):
                    self._keywordsLine[i] = ''
        elif widgetNumber == 4:
            functional = self._semiempiricalFunctionals[selection]
            self._keywordsLine[2] = f'{functional}'
            self._keywordsLine[3] = ''
        elif widgetNumber == 5:
            method = self._mechanics[selection]
            self._keywordsLine[2] = f'{method}'
            self._keywordsLine[3] = ''
        elif widgetNumber == 7:
            if selection == 1:
                self._keywordsLine[3] = self._augmented[selection]
            else:
                self._keywordsLine[3] = ''
        elif widgetNumber == 8:
            self._keywordsLine[4] = self._basisSet[selection]

            if selection == 0:
                self._methodWidgets[10].setVisible(True)
                self._methodWidgets[9].setVisible(False)
                [self._methodWidgets[i].setCurrentIndex(0) for i in (10, 12, 14)]
                self._methodWidgets[10].model().item(2).setEnabled(False)
                for i in range(6, 12):
                    self._keywordsLine[i] = ''
            elif selection == 1:
                self._methodWidgets[9].setVisible(True)
                self._methodWidgets[10].setVisible(True)
                self._methodWidgets[9].setCurrentIndex(1)
                [self._methodWidgets[i].setCurrentIndex(0) for i in (10, 12, 14)]
                self._methodWidgets[9].model().item(2).setEnabled(False)
                self._methodWidgets[10].model().item(2).setEnabled(True)
                [self._methodWidgets[i].setVisible(False) for i in range(11, 16)]
                self._keywordsLine[4] = self._basisSet[selection].replace('G', '+G')
                for i in range(6, 9):
                    self._keywordsLine[i] = ''
            elif selection in (2, 3, 4):
                [self._methodWidgets[i].setVisible(True) for i in range(9, 16) if i != 10]
                self._methodWidgets[10].setVisible(False)
                self._methodWidgets[9].setCurrentIndex(1)
                self._methodWidgets[9].model().item(2).setEnabled(True)
                self._keywordsLine[4] = self._basisSet[selection].replace('G', '+G')
                for i in range(5, 9):
                    self._keywordsLine[i] = ''
            else:
                [self._methodWidgets[i].setVisible(False) for i in range(9, 16)]

            if selection in (5, 6, 7):
                self._methodWidgets[7].setVisible(True)
                [self._methodWidgets[i].setCurrentIndex(0) for i in (7, 12, 14)]
                for i in range(5, 9):
                    self._keywordsLine[i] = ''
            else:
                self._methodWidgets[7].setVisible(False)
                self._keywordsLine[3] = ''
        elif widgetNumber == 9:
            if selection == 1:
                self._keywordsLine[4] = self._basisSet[self._methodWidgets[8].currentIndex()].replace('G', '+G')
            elif selection == 2:
                self._keywordsLine[4] = self._basisSet[self._methodWidgets[8].currentIndex()].replace('G', '++G')
            else:
                self._keywordsLine[4] = self._basisSet[self._methodWidgets[8].currentIndex()]
        elif widgetNumber == 10:
            if selection in (1, 2):
                self._keywordsLine[5] = self._asterisk[selection]
            else:
                self._keywordsLine[5] = ''
        elif widgetNumber == 12:
            if selection != 0:
                self._keywordsLine[6] = '(' + self._firstPolarizedOrbitals[selection]
                self._keywordsLine[8] = ')'
            else:
                self._keywordsLine[6] = ''
                self._keywordsLine[8] = ''
        elif widgetNumber == 14:
            if selection != 0:
                self._keywordsLine[7] = ',' + self._secondPolarizedOrbitals[selection]
            else:
                self._keywordsLine[7] = ''
        elif widgetNumber == 27:
            if selection == 0:
                self._keywordsLine[9] = ''
            elif selection == 1:
                self._keywordsLine[9] = '/FIT'
            elif selection == 3:
                self._methodWidgets[34].setVisible(True)
                self._keywordsLine[9] = '/' + self._fittingSet[selection].replace('N', str(self._methodWidgets[34].value()))
            else:
                self._keywordsLine[9] = '/' + self._fittingSet[selection]
                self._methodWidgets[34].setVisible(False)
        elif widgetNumber == 33:
            if selection:
                self._keywordsLine[11] = ' SPARSE'
            else:
                self._keywordsLine[11] = ''
        elif widgetNumber == 34:
            self._keywordsLine[9] = re.sub(r'\d+', str(selection), self._keywordsLine[9])

        self.uiKeywordsLabel.setText(''.join(self._keywordsLine))

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

    def SetCahrgeMult(self):  # mejorar
        """
        docstring
        """
        self._chargeMultiplicityLine[0] = str(self._methodWidgets[30].value())

        self._currentElectrons = self._molecule.GetValenceElectrons - self._methodWidgets[30].value() + self._molecule.GetNetCharge
        combination = self._currentElectrons % 2
        if combination == 0:
            self._methodWidgets[32].setVisible(True)
            self._methodWidgets[31].setVisible(False)
            self._methodWidgets[31].setValue(2)
            self._chargeMultiplicityLine[2] = str(self._methodWidgets[32].value())
        else:
            self._methodWidgets[32].setVisible(False)
            self._methodWidgets[31].setVisible(True)
            self._methodWidgets[32].setValue(1)
            self._chargeMultiplicityLine[2] = str(self._methodWidgets[31].value())
        
        if combination == 0 and self._methodWidgets[32].value() == 1:
            self._methodWidgets[2].model().item(1).setEnabled(True)
        else:
            self._methodWidgets[2].setCurrentIndex(0)
            self._methodWidgets[2].model().item(1).setEnabled(False)
        
        self.uiChargeMultLabel.setText(''.join(self._chargeMultiplicityLine))
# por aquí