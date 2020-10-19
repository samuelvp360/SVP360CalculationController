#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from PyQt5 import QtWidgets as qtw
from PyQt5 import uic
# from Models import StatusModel
import multiprocessing
from psutil import virtual_memory


class CalculationsController(qtw.QWidget):

    def __init__(self):
        super(CalculationsController, self).__init__()
        uic.loadUi('Views/uiCalculationsWindow.ui', self)
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
        # self._addInputs = []

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
        self._state = ('Ground State', 'ZINDO', 'CIS') #, 'TD-SCF', 'TDA', 'EOM-CCSD')
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
        self._qeqCharges = (
            "Don't use", 'All atoms', 'Untyped atoms', 'Uncharged atoms'
        )
        self._includeExclude = ('Include Triples', 'Exclude Triples')
        self._includeExclude2 = ('Exclude Triples', 'Include Triples', 'Include MP4 Triples')
        self._states = ('Default', 'Singlet Only', 'Triplet Only', 'Singlets & Triplets')
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
            self.uiAutoValueSpinBox,  # 34
            self.uiQeqChargesLabel,  # 35
            self.uiQeqChargesComboBox,  # 36
            self.uiIncludeAllElectronsCheckBox,  # 37
            self.uiIncludeExcludeComboBox,  # 38
            self.uiReadAmpCheckBox,  # 39
            self.uiSaveAmpCheckBox,  # 40
            self.uiIncludeExclude2ComboBox,  # 41
            self.uiNumElectLabel,  # 42
            self.uiNumOrbLabel,  # 43
            self.uiNumElectSpinBox,  # 44
            self.uiNumOrbSpinBox,  # 45
            self.uiRFOCheckBox,  # 46
            self.uiStatesLabel,  # 47
            self.uiStatesComboBox,  # 48
            self.uiMoreStatesCheckBox,  # 49
            self.uiMoreStatesSpinBox,  # 50
            self.uiStateOfInterestCheckBox,  # 51
            self.uiStateOfInterestSpinBox,  # 52
        ]
        self.lastMethodWidget = len(self._methodWidgets)

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
        self._methodWidgets[36].addItems(self._qeqCharges)
        self._methodWidgets[38].addItems(self._includeExclude)
        self._methodWidgets[41].addItems(self._includeExclude2)
        self._methodWidgets[48].addItems(self._states)
        [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in (6, 8, 10, 30, 31, 32)]

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
        self._methodWidgets[0].currentIndexChanged.connect(lambda: self.SetKeywords(0, self._methodWidgets[0].currentIndex()))
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
        self._methodWidgets[36].currentIndexChanged.connect(lambda: self.SetKeywords(36, self._methodWidgets[36].currentIndex()))
        self._methodWidgets[37].stateChanged.connect(lambda: self.SetKeywords(37, self._methodWidgets[37].isChecked()))
        self._methodWidgets[38].currentIndexChanged.connect(lambda: self.SetKeywords(38, self._methodWidgets[38].currentIndex()))
        self._methodWidgets[39].stateChanged.connect(lambda: self.SetKeywords(39, self._methodWidgets[39].isChecked()))
        self._methodWidgets[40].stateChanged.connect(lambda: self.SetKeywords(40, self._methodWidgets[40].isChecked()))
        self._methodWidgets[41].currentIndexChanged.connect(lambda: self.SetKeywords(41, self._methodWidgets[41].currentIndex()))
        self._methodWidgets[44].valueChanged.connect(lambda: self.SetKeywords(44, self._methodWidgets[44].value()))
        self._methodWidgets[45].valueChanged.connect(lambda: self.SetKeywords(45, self._methodWidgets[45].value()))
        self._methodWidgets[46].stateChanged.connect(lambda: self.SetKeywords(46, self._methodWidgets[46].isChecked()))
        self._methodWidgets[48].currentIndexChanged.connect(lambda: self.SetKeywords(48, self._methodWidgets[48].currentIndex()))
        self._methodWidgets[49].stateChanged.connect(lambda: self.SetKeywords(49, self._methodWidgets[49].isChecked()))
        self._methodWidgets[50].valueChanged.connect(lambda: self.SetKeywords(50, self._methodWidgets[50].value()))
        self._methodWidgets[51].stateChanged.connect(lambda: self.SetKeywords(51, self._methodWidgets[51].isChecked()))
        self._methodWidgets[52].valueChanged.connect(lambda: self.SetKeywords(52, self._methodWidgets[52].value()))
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

        if widgetNumber == 0:
            if selection == 0:
                self._methodWidgets[1].setVisible(True)
                self.SetKeywords(1, self._methodWidgets[1].currentIndex())
            elif selection == 1:
                self._methodWidgets[2].model().item(3).setEnabled(False)
                self._methodWidgets[2].setCurrentIndex(0)
                for i in range(1, 12):
                    self._keywordsLine[i] = ''
                self._keywordsLine[2] = 'ZINDO'
                [self._methodWidgets[i].setVisible(False) for i in range(1, self.lastMethodWidget) if i not in (2, 30, 31, 32)]
                [self._methodWidgets[i].setVisible(True) for i in range(47, 53)]
                self.SetKeywords(48, self._methodWidgets[48].currentIndex())
            elif selection == 2:
                self._methodWidgets[2].model().item(3).setEnabled(False)
                self._methodWidgets[2].setCurrentIndex(0)
                for i in range(1, 12):
                    self._keywordsLine[i] = ''
                self._keywordsLine[2] = 'CIS'
                [self._methodWidgets[i].setVisible(False) for i in range(1, self.lastMethodWidget) if i not in (2, 30, 31, 32)]
                [self._methodWidgets[i].setVisible(True) for i in (6, 8)]
                [self._methodWidgets[i].setVisible(True) for i in range(47, 53)]
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (8, 48)]
            # elif selection == 3:
            #     self._methodWidgets[2].model().item(3).setEnabled(False)
            #     self._methodWidgets[2].setCurrentIndex(0)
            #     for i in range(1, 12):
            #         self._keywordsLine[i] = ''
            #     [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in (30, 31, 32)]
            #     [self._methodWidgets[i].setVisible(True) for i in (1, 2, 6, 8)]
            #     [self._methodWidgets[i].setVisible(True) for i in range(47, 53)]
            #     [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (8, 48)]
        elif widgetNumber == 1:
            self._methodWidgets[2].model().item(3).setEnabled(True)
            if selection == 0:
                for i in range(1, 12):
                    self._keywordsLine[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(2, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (5, 35, 36)]
                self.SetKeywords(5, self._methodWidgets[5].currentIndex())
            elif selection == 1:
                for i in range(2, 12):
                    self._keywordsLine[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 4, 33)]
                self._methodWidgets[33].setEnabled(True)
                # functional = self._semiempiricalFunctionals[self._methodWidgets[4].currentIndex()]
                # self._keywordsLine[2] = f'{functional}'
                self.SetKeywords(4, self._methodWidgets[4].currentIndex())
                self.SetKeywords(33, self._methodWidgets[33].isChecked())
            elif selection == 2:
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 10)]
                for i in range(2, 12):
                    self._keywordsLine[i] = ''
                self._keywordsLine[2] = 'HF/'
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
            elif selection == 3:
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 3, 6, 8, 33)]
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (3, 8)]
            elif selection == 4:
                for i in range(2, 12):
                    self._keywordsLine[i] = ''
                self._keywordsLine[2] = 'MP2'
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 10, 37)]
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(37, self._methodWidgets[37].isChecked())
            elif selection == 5:
                for i in range(2, 12):
                    self._keywordsLine[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 10, 37, 38)]
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (38, 8)]
                self.SetKeywords(37, self._methodWidgets[37].isChecked())
            elif selection == 6:
                for i in range(2, 12):
                    self._keywordsLine[i] = ''
                self._keywordsLine[2] = 'CCSD'
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 37, 39, 40, 41)]
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(39, 1)
            elif selection == 7:
                for i in range(2, 12):
                    self._keywordsLine[i] = ''
                self._keywordsLine[2] = 'BD'
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 39, 40)]
                self._methodWidgets[41].setCurrentIndex(0)
                self._methodWidgets[37].setChecked(False)
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(39, 1)
            elif selection == 8:
                for i in range(2, 12):
                    self._keywordsLine[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 42, 43, 44, 45, 46)]
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(46, 1)
        elif widgetNumber == 2:
            if selection == 0:
                self._keywordsLine[1] = ''
            elif selection == 1:
                self._keywordsLine[1] = 'R'
            elif selection == 2:
                self._keywordsLine[1] = 'U'
            else:
                self._keywordsLine[1] = 'RO'
        elif widgetNumber == 3:
            for i in range(2, 12):
                self._keywordsLine[i] = ''
            functional = self._dftFunctionals[selection]
            self._keywordsLine[2] = f'{functional}/'
            if selection in (0, 1, 6, 8, 9):
                [self._methodWidgets[i].setVisible(True) for i in (26, 27)]
                self._methodWidgets[33].setEnabled(True)
                self.SetKeywords(33, self._methodWidgets[33].isChecked())
                self.SetKeywords(27, self._methodWidgets[27].currentIndex())
            else:
                [self._methodWidgets[i].setVisible(False) for i in (26, 27, 34)]
                self._methodWidgets[33].setEnabled(False)
            self.SetKeywords(8, self._methodWidgets[8].currentIndex())
        elif widgetNumber == 4:
            functional = self._semiempiricalFunctionals[selection]
            self._keywordsLine[2] = f'{functional}'
            # self._keywordsLine[3] = ''
        elif widgetNumber == 5:
            method = self._mechanics[selection]
            self._keywordsLine[2] = f'{method}'
            self.SetKeywords(36, self._methodWidgets[36].currentIndex())
        elif widgetNumber == 7:
            if selection == 1:
                self._keywordsLine[3] = self._augmented[selection]
            else:
                self._keywordsLine[3] = ''
        elif widgetNumber == 8:
            for i in range(4, 9):
                self._keywordsLine[i] = ''
            self._keywordsLine[4] = self._basisSet[selection]

            if selection == 0:
                self._methodWidgets[7].setVisible(False)
                self._methodWidgets[10].setVisible(True)
                self._methodWidgets[9].setVisible(False)
                self._methodWidgets[10].setCurrentIndex(0)
                self._methodWidgets[10].model().item(2).setEnabled(False)
                self.SetKeywords(10, self._methodWidgets[10].currentIndex())
            elif selection == 1:
                self._methodWidgets[7].setVisible(False)
                [self._methodWidgets[i].setVisible(True) for i in (9, 10)]
                [self._methodWidgets[i].setVisible(False) for i in range(11, 16)]
                self._methodWidgets[9].setCurrentIndex(0)
                self._methodWidgets[9].model().item(2).setEnabled(False)
                self._methodWidgets[10].model().item(2).setEnabled(True)
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (9, 10)]
            elif selection in range(2, 5):
                self._methodWidgets[7].setVisible(False)
                [self._methodWidgets[i].setVisible(True) for i in range(9, 16) if i != 10]
                self._methodWidgets[10].setVisible(False)
                self._methodWidgets[9].model().item(2).setEnabled(True)
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (9, 12)]
            elif selection in range(5, 8):
                [self._methodWidgets[i].setVisible(False) for i in range(9, 16)]
                self._methodWidgets[7].setVisible(True)
                self.SetKeywords(7, self._methodWidgets[7].currentIndex())
            else:
                [self._methodWidgets[i].setVisible(False) for i in range(7, 16) if i != 8]
        elif widgetNumber == 9:
            basis = self._basisSet[self._methodWidgets[8].currentIndex()]
            if selection == 1:
                self._keywordsLine[4] = basis.replace('G', '+G')
            elif selection == 2:
                self._keywordsLine[4] = basis.replace('G', '++G')
            else:
                self._keywordsLine[4] = basis
        elif widgetNumber == 10:
            if selection in (1, 2):
                self._keywordsLine[5] = self._asterisk[selection]
            else:
                self._keywordsLine[5] = ''
        elif widgetNumber in (12, 14):
            first = self._methodWidgets[12].currentIndex()
            second = self._methodWidgets[14].currentIndex()
            if first != 0 or second != 0:
                self._keywordsLine[6] = '('
                self._keywordsLine[8] = ')'
                if first != 0 and second == 0:
                    self._keywordsLine[7] = self._firstPolarizedOrbitals[first]
                elif first == 0 and second != 0:
                    self._keywordsLine[7] = self._secondPolarizedOrbitals[second]
                elif first != 0 and second != 0:
                    self._keywordsLine[7] = self._firstPolarizedOrbitals[first] + ',' + self._secondPolarizedOrbitals[second]
            else:
                for i in range(6, 9):
                    self._keywordsLine[i] = ''
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
        elif widgetNumber == 36:
            if selection == 0:
                self._keywordsLine[3] = ''
            elif selection == 1:
                self._keywordsLine[3] = '=QEQ'
            elif selection == 2:
                self._keywordsLine[3] = '=UNTYPED'
            else:
                self._keywordsLine[3] = '=UNCHARGED'
        elif widgetNumber == 37:
            selectedMethod = self._methodWidgets[1].currentIndex()
            if selectedMethod == 4:
                if selection:
                    self._keywordsLine[3] = '=FULL/'
                else:
                    self._keywordsLine[3] = '/'
            elif selectedMethod == 5:
                if selection:
                    self._keywordsLine[3] = ',FULL)/'
                else:
                    self._keywordsLine[3] = ')/'
            elif selectedMethod == 6:
                self.SetKeywords(39, self._methodWidgets[39].isChecked())
        elif widgetNumber == 38:
            if selection == 0:
                self._keywordsLine[2] = 'MP4(SDTQ'
            else:
                self._keywordsLine[2] = 'MP4(SDQ'
        elif widgetNumber in range(39, 42):
            first = self._methodWidgets[39].isChecked()
            if first:
                firstVariable = 'READAMPLITUDES,'
            else:
                firstVariable = ''
            second = self._methodWidgets[40].isChecked()
            if second:
                secondVariable = 'SAVEAMPLITUDES,'
            else:
                secondVariable = ''
            third = self._methodWidgets[37].isChecked()
            if third:
                thirdVariable = 'FULL,'
            else:
                thirdVariable = ''
            fourth = self._methodWidgets[41].currentIndex()
            if fourth == 1:
                fourthVariable = 'T,'
            elif fourth == 2:
                fourthVariable = 'T,E4T,'
            else:
                fourthVariable = ''
            variables = (fourthVariable, thirdVariable, firstVariable, secondVariable)
            combination = (first, second, third, fourth)
            if sum(combination) > 0:
                if first + second == 1 and third + fourth == 0:
                    self._keywordsLine[3] = '=' + re.sub(r'\,$', '', ''.join(variables)) + '/'
                elif first == 0 and second == 0:
                    self._keywordsLine[3] = '(' + re.sub(r'\,$', '', ''.join(variables)) + ')/'
                else:
                    self._keywordsLine[3] = '=(' + re.sub(r'\,$', '', ''.join(variables)) + ')/'
            else:
                self._keywordsLine[3] = '/'
        elif widgetNumber in range(44, 47):
            self._keywordsLine[2] = f'CASSCF({self._methodWidgets[44].value()},{self._methodWidgets[45].value()}'
            if self._methodWidgets[46].isChecked():
                self._keywordsLine[3] = ',RFO)/'
            else:
                self._keywordsLine[3] = ')/'
        elif widgetNumber in range(48, 53):
            first = self._methodWidgets[48].currentIndex()
            if first == 0:
                firstVariable = ''
            elif first == 1:
                firstVariable = 'SINGLETS,'
            elif first == 2:
                firstVariable = 'TRIPLETS,'
            else:
                firstVariable = '50-50,'
            second = self._methodWidgets[49].isChecked()
            if second:
                self._methodWidgets[50].setEnabled(True)
                secondVariable = f'NSTATES={str(self._methodWidgets[50].value())},'
            else:
                self._methodWidgets[50].setEnabled(False)
                secondVariable = ''
            third = self._methodWidgets[51].isChecked()
            if third:
                self._methodWidgets[52].setEnabled(True)
                thirdVariable = f'ROOT={str(self._methodWidgets[52].value())},'
            else:
                self._methodWidgets[52].setEnabled(False)
                thirdVariable = ''
            variables = (firstVariable, secondVariable, thirdVariable)
            combination = (first, second, third)
            if sum(combination) > 0:
                if second or third:
                    self._keywordsLine[3] = '=(' + re.sub(r'\,$', '', ''.join(variables)) + ')/'
                else:
                    self._keywordsLine[3] = '=' + re.sub(r'\,$', '', ''.join(variables)) + '/'
            else:
                self._keywordsLine[3] = '/'
            if self._methodWidgets[0].currentIndex() == 1:
                self._keywordsLine[3] = re.sub(r'/$', '', self._keywordsLine[3])

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
