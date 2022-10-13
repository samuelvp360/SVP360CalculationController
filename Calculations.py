#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import subprocess
import os
import io
import shutil
import copy
from glob import glob
import multiprocessing
import numpy as np
import pandas as pd
from os import listdir
from PIL import Image, ImageQt
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from psutil import virtual_memory
# from vina import Vina
from Models import PandasModel, DisplayMFPModel, SimilarityModel
from Components import Similarity
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
# from openbabel import openbabel as ob
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
# from mpldatacursor import datacursor
from sklearn.metrics import r2_score
from scipy.cluster.hierarchy import dendrogram, linkage
# ob.obErrorLog.StopLogging()  # se puede cambiar por SetOutputLevel para logging
matplotlib.use('Qt5Agg')


class PlotCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, parent=None):
        self.fig = Figure(
            figsize=(12, 8), dpi=100,
            tight_layout=True
        )
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)


class Gaussian(qtw.QMainWindow):

    submitted = qtc.pyqtSignal(dict, str)

    def __init__(
        self, molecule, *, keywords=False, group=False,
        job_type=False, coordinates= False, project=False
    ):
        super().__init__()
        uic.loadUi('Views/uiCalculationsWindow.ui', self)
        self.uiOptGroupBox.setVisible(False)
        self.uiFreqGroupBox.setVisible(False)
        self.uiIRCGroupBox.setVisible(False)
    # ------------------------------------PROPERTIES-------------------------------------
        self.group = group
        self.job_type = job_type
        self.coordinates = coordinates
        self.project = project
        self.memory = virtual_memory().total // 1e9
        self.cpu = multiprocessing.cpu_count()
        self.keywords_line = []
        for i in range(20):
            self.keywords_line.append('')
        self.keywords_line[0] = '# '
        self._link0Line = ['', '', '', '']
        self._chargeMultiplicityLine = ['', ' ', '']
        self._lightBasis = ['', '', '', '', '', '']
        self.additional_input = ''
        self.populate_widgets()
    # ------------------------------------SIGNALS----------------------------------------
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
        self._methodWidgets[12].currentIndexChanged.connect(lambda: self.SetKeywords(12))
        self._methodWidgets[14].currentIndexChanged.connect(lambda: self.SetKeywords(14))
        self._methodWidgets[17].currentIndexChanged.connect(lambda: self.SetKeywords(17, self._methodWidgets[17].currentIndex()))
        self._methodWidgets[18].currentIndexChanged.connect(lambda: self.SetKeywords(18, self._methodWidgets[18].currentIndex()))
        self._methodWidgets[19].currentIndexChanged.connect(lambda: self.SetKeywords(19, self._methodWidgets[19].currentIndex()))
        self._methodWidgets[20].currentIndexChanged.connect(lambda: self.SetKeywords(20, self._methodWidgets[20].currentIndex()))
        self._methodWidgets[22].currentIndexChanged.connect(lambda: self.SetKeywords(22))
        self._methodWidgets[24].currentIndexChanged.connect(lambda: self.SetKeywords(24))
        self._methodWidgets[27].currentIndexChanged.connect(lambda: self.SetKeywords(27, self._methodWidgets[27].currentIndex()))
        self._methodWidgets[29].currentIndexChanged.connect(lambda: self.SetKeywords(29, self._methodWidgets[29].currentIndex()))
        self._methodWidgets[30].valueChanged.connect(self.SetCahrgeMult)
        self._methodWidgets[31].valueChanged.connect(self.SetCahrgeMult)
        self._methodWidgets[32].valueChanged.connect(self.SetCahrgeMult)
        self._methodWidgets[33].stateChanged.connect(lambda: self.SetKeywords(33, self._methodWidgets[33].isChecked()))
        self._methodWidgets[34].valueChanged.connect(lambda: self.SetKeywords(34, self._methodWidgets[34].value()))
        self._methodWidgets[36].currentIndexChanged.connect(lambda: self.SetKeywords(36, self._methodWidgets[36].currentIndex()))
        self._methodWidgets[37].stateChanged.connect(lambda: self.SetKeywords(37, self._methodWidgets[37].isChecked()))
        self._methodWidgets[38].currentIndexChanged.connect(lambda: self.SetKeywords(38, self._methodWidgets[38].currentIndex()))
        self._methodWidgets[39].stateChanged.connect(lambda: self.SetKeywords(39))
        self._methodWidgets[40].stateChanged.connect(lambda: self.SetKeywords(40))
        self._methodWidgets[41].currentIndexChanged.connect(lambda: self.SetKeywords(41))
        self._methodWidgets[44].valueChanged.connect(lambda: self.SetKeywords(44))
        self._methodWidgets[45].valueChanged.connect(lambda: self.SetKeywords(45))
        self._methodWidgets[46].stateChanged.connect(lambda: self.SetKeywords(46))
        self._methodWidgets[48].currentIndexChanged.connect(lambda: self.SetKeywords(48))
        self._methodWidgets[49].stateChanged.connect(lambda: self.SetKeywords(49))
        self._methodWidgets[50].valueChanged.connect(lambda: self.SetKeywords(50))
        self._methodWidgets[51].stateChanged.connect(lambda: self.SetKeywords(51))
        self._methodWidgets[52].valueChanged.connect(lambda: self.SetKeywords(52))
        self.uiAdditionalKeyLine.textChanged.connect(lambda: self.SetKeywords(60, self.uiAdditionalKeyLine.text()))
        self._jobsWidgets[0].currentIndexChanged.connect(lambda: self.SetJobType(0, self._jobsWidgets[0].currentIndex()))
        self._jobsWidgets[1].currentIndexChanged.connect(lambda: self.SetJobType(1, self._jobsWidgets[1].currentIndex()))
        self._jobsWidgets[2].currentIndexChanged.connect(lambda: self.SetJobType(2, self._jobsWidgets[2].currentIndex()))
        self._jobsWidgets[3].stateChanged.connect(lambda: self.SetJobType(3, self._jobsWidgets[3].isChecked()))
        self._jobsWidgets[4].stateChanged.connect(lambda: self.SetJobType(4, self._jobsWidgets[4].isChecked()))
        self._jobsWidgets[5].currentIndexChanged.connect(lambda: self.SetJobType(5, self._jobsWidgets[5].currentIndex()))
        self._jobsWidgets[6].currentIndexChanged.connect(lambda: self.SetJobType(6, self._jobsWidgets[6].currentIndex()))
        self._jobsWidgets[7].currentIndexChanged.connect(lambda: self.SetJobType(7, self._jobsWidgets[7].currentIndex()))
        self._jobsWidgets[8].stateChanged.connect(lambda: self.SetJobType(8, self._jobsWidgets[8].isChecked()))
        self._jobsWidgets[9].stateChanged.connect(lambda: self.SetJobType(9, self._jobsWidgets[9].isChecked()))
        self._jobsWidgets[10].stateChanged.connect(lambda: self.SetJobType(10, self._jobsWidgets[10].isChecked()))
        self._jobsWidgets[11].stateChanged.connect(lambda: self.SetJobType(11, self._jobsWidgets[11].isChecked()))
        self._jobsWidgets[12].stateChanged.connect(lambda: self.SetJobType(12, self._jobsWidgets[12].isChecked()))
        self._jobsWidgets[15].currentIndexChanged.connect(lambda: self.SetJobType(15, self._jobsWidgets[15].currentIndex()))
        self._jobsWidgets[16].currentIndexChanged.connect(lambda: self.SetJobType(16, self._jobsWidgets[16].currentIndex()))
        self._jobsWidgets[17].currentIndexChanged.connect(lambda: self.SetJobType(17, self._jobsWidgets[17].currentIndex()))
        self._jobsWidgets[18].stateChanged.connect(lambda: self.SetJobType(18, self._jobsWidgets[18].isChecked()))
        self._jobsWidgets[19].valueChanged.connect(lambda: self.SetJobType(19, self._jobsWidgets[19].value()))
        self._jobsWidgets[20].stateChanged.connect(lambda: self.SetJobType(20, self._jobsWidgets[20].isChecked()))
        self._jobsWidgets[21].valueChanged.connect(lambda: self.SetJobType(21, self._jobsWidgets[21].value()))
        self._jobsWidgets[22].stateChanged.connect(lambda: self.SetJobType(22, self._jobsWidgets[22].isChecked()))
        self.link_0_widgets[0].valueChanged.connect(lambda: self.SetLink0(0, self.link_0_widgets[0].value()))
        self.link_0_widgets[1].valueChanged.connect(lambda: self.SetLink0(1, self.link_0_widgets[1].value()))
        self.link_0_widgets[2].currentIndexChanged.connect(lambda: self.SetLink0(2, self.link_0_widgets[2].currentIndex()))
        self.link_0_widgets[3].textChanged.connect(lambda: self.SetLink0(3, self.link_0_widgets[3].text()))
        self.link_0_widgets[4].clicked.connect(lambda: self.SetLink0(4, 1))
        self.link_0_widgets[5].currentIndexChanged.connect(lambda: self.SetLink0(5, self.link_0_widgets[5].currentIndex()))
        self.link_0_widgets[6].textChanged.connect(lambda: self.SetLink0(6, self.link_0_widgets[6].text()))
        self.link_0_widgets[7].clicked.connect(lambda: self.SetLink0(7, 1))
        self._generalWidgets[0].stateChanged.connect(lambda: self.SetGeneral(0, self._generalWidgets[0].isChecked()))
        self._generalWidgets[1].stateChanged.connect(lambda: self.SetGeneral(1, self._generalWidgets[1].isChecked()))
        self._generalWidgets[2].stateChanged.connect(lambda: self.SetGeneral(2, self._generalWidgets[2].isChecked()))
        self._generalWidgets[7].stateChanged.connect(lambda: self.SetGeneral(7, self._generalWidgets[7].isChecked()))
        self._generalWidgets[9].stateChanged.connect(lambda: self.SetGeneral(9, self._generalWidgets[9].isChecked()))
        self._generalWidgets[11].stateChanged.connect(lambda: self.SetGeneral(11, self._generalWidgets[11].isChecked()))
        self._generalWidgets[12].valueChanged.connect(lambda: self.SetGeneral(12, self._generalWidgets[12].value()))
        self.uiTitleLineEdit.textChanged.connect(self.set_title)
        # self.uiAddInputCheckBox.stateChanged.connect(self.set_preview)
        # self.uiQueueCalcButton.clicked.connect(self.queu_calculation)

        self.prepare_molecule(molecule)
        if keywords:
            self.uiKeywordsLabel.setText(keywords)
            self.set_preview()
        else:
            self.set_preview()
            self.show()

    def prepare_molecule(self, molecule):
        self._molecule = molecule
        self._chk = self._molecule.get_name
        self._oldChk = self._molecule.get_name
        if not self.coordinates:
            coords = str(self._molecule.get_coordinates)
            self.coordinates = coords.replace('[', '').replace(']', '').replace('\n ', '\n')
        self.uiTitleLineEdit.setText(self._molecule.get_name)
        self.SetCahrgeMult()
        self.SetKeywords(1, self._methodWidgets[1].currentIndex())
        self._methodWidgets[30].setValue(self._molecule.get_formal_charge)
        self._methodWidgets[10].model().item(2).setEnabled(False)
        [self.SetLink0(i, self.link_0_widgets[i].value()) for i in range(2)]
        [self.SetLink0(i, self.link_0_widgets[i].currentIndex()) for i in (2, 5)]

    def populate_widgets(self):
        self._jobTypes = (
            'Energy', 'Optimization', 'Frequency', 'Opt+Freq', 'IRC', 'Scan', 'Stability', 'NMR'
        )
        self._optimizeToA = ('Minimum', 'TS (Berny)')
        self._forceConstants = (
            None, 'Calculate at First Point', 'Calculate at all Points', 'Read Internal', 'Read Cartesian'
        )
        self._computeRaman = (None, 'No', 'Yes')
        self._saveNormalModesOptions = (None, 'Yes', 'No')
        self._computeROA = ('No', 'Yes')
        self._readIncidentLightFreqs = (None, 'Yes', 'No')
        self._forceConstants2 = ('Calculate Once', 'Calculate Always', 'Read from .CHK')
        self._recorrect = (None, 'Never', 'Yes', 'Always', 'Test')
        self._followIRC = ('Both directions', 'Forward only', 'Reverse only')
        self._state = ('Ground State', 'ZINDO', 'CIS')  # , 'TD-SCF', 'TDA', 'EOM-CCSD')
        self._methods = (
            'Mechanics...', 'Semi-empirical...', 'Hartree-Fock', 'DFT...', 'MP2', 'MP4', 'CCSD', 'BD', 'CASSCF'
        )
        self._shellType = ('Default Spin', 'Restricted', 'Unrestricted', 'Restricted-Open')
        self._semiempiricalFunctionals = ('PM6', 'PDDG', 'AM1', 'PM3', 'PM3MM', 'INDO', 'CNDO')
        self._dftFunctionals = (
            'LSDA', 'BPV86', 'B3LYP', 'CAM-B3LYP', 'B3PW91', 'MPW1PW91', 'PBEPBE', 'HSEH1PBE',
            'HCTH', 'TPSSTPSS', 'WB97XD', 'APFD', 'MPWB1K', 'M06-2X'
        )
        self._mechanics = ('UFF', 'Dreiding', 'Amber')
        self._augmented = (None, 'aug-')
        self._basisSet = (
            'STO-3G', '3-21G', '6-31G', '6-31G\'', '6-311G', 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ',
            'LanL2DZ', 'LanL2MB', 'SDD', 'DGDZVP', 'DGDZVP2', 'DGTZVP', 'GEN', 'GENECP'
        )
        self._heavyBasisSet = ('LanL2DZ', 'LanL2MB', 'SDD', 'DGDZVP', 'DGDZVP2', 'DGTZVP')
        self._lightBasisSet = (
            'STO-3G', '3-21G', '6-31G', '6-31G\'', '6-311G', 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ'
        )
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
        self._chk = ("Don't save", 'Default name', 'Specify...')
        self._oldChk = ('No', 'Default name', 'Specify...')
        self._jobsWidgets = [
            self.uiJobTypeComboBox,  # 0
            self.uiOptimizeToAComboBox,  # 1
            self.uiForceConstComboBox,  # 2
            self.uiUseRFOStepsCheckBox,  # 3
            self.uiTightCheckBox,  # 4
            self.uiRamanComboBox,  # 5
            self.uiROAComboBox,  # 6
            self.uiNormalModeComboBox,  # 7
            self.uiVCDCheckBox,  # 8
            self.uiSkipCheckBox,  # 9
            self.uiAnharmonicCheckBox,  # 10
            self.uiAnharmonicModesCheckBox,  # 11
            self.uiProjectedFreqCheckBox,  # 12
            self.uiAnharmonicSpinBox,  # 13
            self.uiIncidentLightComboBox,  # 14
            self.uiFollowIRCComboBox,  # 15
            self.uiForceConst2ComboBox,  # 16
            self.uiRecorrectComboBox,  # 17
            self.uiMorePointsCheckBox,  # 18
            self.uiMorePointsSpinBox,  # 19
            self.uiRecalculateCheckBox,  # 20
            self.uiRecalculateSpinBox,  # 21
            self.uiIRCMaxCheckBox  # 22
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
            self.uiBasis3Label,  # 26
            self.uiBasis3ComboBox,  # 27
            self.uiFittingSetLabel,  # 28
            self.uiFittingSetComboBox,  # 29
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
        self.link_0_widgets = [
            self.uiMemorySpinBox,  # 0
            self.uiProcessorsSpinBox,  # 1
            self.uiChkComboBox,  # 2
            self.uiChkLine,  # 3
            self.uiSaveChkButton,  # 4
            self.uiOldChkComboBox,  # 5
            self.uiOldChkLine,  # 6
            self.uiSelectOldChkButton,  # 7
            self.uiLink0EditPlainText  # 8
        ]
        self._generalWidgets = [
            self.uiQuadraticallyCheckBox,  # 0
            self.uiIgnoreSymmetryCheckBox,  # 1
            self.uiWriteCartesiansCheckBox,  # 2
            self.uiModifiedRedundantCheckBox,  # 3
            self.uiCounterpoiseCheckBox,  # 4
            self.uiWriteConnectivityCheckBox,  # 5
            self.uiGaussianFragmentCheckBox,  # 6
            self.uiAdditionalPrintCheckBox,  # 7
            self.uiWritePDBCheckBox,  # 8
            self.uiPolarizabilitiesCheckBox,  # 9
            self.uiOpticalRotationsCheckBox,  # 10
            self.uiMaxDiskCheckBox,  # 11
            self.uiMaxDiskSpinBox  # 12
        ]

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
        self._methodWidgets[18].addItems(self._lightBasisSet)
        self._methodWidgets[19].addItems(self._diffuse)
        self._methodWidgets[20].addItems(self._asterisk)
        self._methodWidgets[22].addItems(self._firstPolarizedOrbitals)
        self._methodWidgets[24].addItems(self._secondPolarizedOrbitals)
        self._methodWidgets[27].addItems(self._heavyBasisSet)
        self._methodWidgets[29].addItems(self._fittingSet)
        self._methodWidgets[36].addItems(self._qeqCharges)
        self._methodWidgets[38].addItems(self._includeExclude)
        self._methodWidgets[41].addItems(self._includeExclude2)
        self._methodWidgets[48].addItems(self._states)
        [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in (6, 8, 10, 30, 31, 32)]

        self._jobsWidgets[0].addItems(self._jobTypes)
        self._jobsWidgets[1].addItems(self._optimizeToA)
        self._jobsWidgets[2].addItems(self._forceConstants)
        self._jobsWidgets[5].addItems(self._computeRaman)
        self._jobsWidgets[6].addItems(self._computeROA)
        self._jobsWidgets[7].addItems(self._saveNormalModesOptions)
        self._jobsWidgets[14].addItems(self._readIncidentLightFreqs)
        self._jobsWidgets[15].addItems(self._followIRC)
        self._jobsWidgets[16].addItems(self._forceConstants2)
        self._jobsWidgets[17].addItems(self._recorrect)

        self.link_0_widgets[0].setMaximum(self.memory)
        self.link_0_widgets[0].setValue(self.memory - 2)
        self.link_0_widgets[1].setMaximum(self.cpu)
        self.link_0_widgets[1].setValue(self.cpu - 1)
        self.link_0_widgets[2].addItems(self._chk)
        [self.link_0_widgets[i].setEnabled(False) for i in range(3, 8) if i != 5]
        self.link_0_widgets[5].addItems(self._oldChk)
        [self._generalWidgets[i].setEnabled(False) for i in (3, 4, 5, 6, 8, 10, 12)]

    def SetKeywords(self, widgetNumber=None, selection=None):

        if widgetNumber == 0:
            if selection == 0:
                self._methodWidgets[1].setVisible(True)
                self.SetKeywords(1, self._methodWidgets[1].currentIndex())
            elif selection == 1:
                self._methodWidgets[2].model().item(3).setEnabled(False)
                self._methodWidgets[2].setCurrentIndex(0)
                for i in range(1, 11):
                    self.keywords_line[i] = ''
                self.keywords_line[2] = 'ZINDO'
                [self._methodWidgets[i].setVisible(False) for i in range(1, self.lastMethodWidget) if i not in (2, 30, 31, 32)]
                [self._methodWidgets[i].setVisible(True) for i in range(47, 53)]
                self.SetKeywords(48, self._methodWidgets[48].currentIndex())
            elif selection == 2:
                self._methodWidgets[2].model().item(3).setEnabled(False)
                self._methodWidgets[2].setCurrentIndex(0)
                for i in range(1, 11):
                    self.keywords_line[i] = ''
                self.keywords_line[2] = 'CIS'
                [self._methodWidgets[i].setVisible(False) for i in range(1, self.lastMethodWidget) if i not in (2, 30, 31, 32)]
                [self._methodWidgets[i].setVisible(True) for i in (6, 8)]
                [self._methodWidgets[i].setVisible(True) for i in range(47, 53)]
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (8, 48)]
        elif widgetNumber == 1:
            self._methodWidgets[2].model().item(3).setEnabled(True)
            if selection == 0:
                for i in range(1, 11):
                    self.keywords_line[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(2, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (5, 35, 36)]
                self.SetKeywords(5, self._methodWidgets[5].currentIndex())
            elif selection == 1:
                for i in range(2, 11):
                    self.keywords_line[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 4, 33)]
                self._methodWidgets[33].setEnabled(True)
                self.SetKeywords(4, self._methodWidgets[4].currentIndex())
                self.SetKeywords(33, self._methodWidgets[33].isChecked())
            elif selection == 2:
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 10)]
                for i in range(2, 11):
                    self.keywords_line[i] = ''
                self.keywords_line[2] = 'HF/'
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
            elif selection == 3:
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 3, 6, 8, 33)]
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (3, 8)]
            elif selection == 4:
                for i in range(2, 11):
                    self.keywords_line[i] = ''
                self.keywords_line[2] = 'MP2'
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 10, 37)]
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(37, self._methodWidgets[37].isChecked())
            elif selection == 5:
                for i in range(2, 11):
                    self.keywords_line[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 10, 37, 38)]
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (38, 8)]
                self.SetKeywords(37, self._methodWidgets[37].isChecked())
            elif selection == 6:
                for i in range(2, 11):
                    self.keywords_line[i] = ''
                self.keywords_line[2] = 'CCSD'
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 37, 39, 40, 41)]
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(39, 1)
            elif selection == 7:
                for i in range(2, 11):
                    self.keywords_line[i] = ''
                self.keywords_line[2] = 'BD'
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 39, 40)]
                self._methodWidgets[41].setCurrentIndex(0)
                self._methodWidgets[37].setChecked(False)
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(39, 1)
            elif selection == 8:
                for i in range(2, 11):
                    self.keywords_line[i] = ''
                [self._methodWidgets[i].setVisible(False) for i in range(3, self.lastMethodWidget) if i not in range(30, 33)]
                [self._methodWidgets[i].setVisible(True) for i in (2, 6, 8, 42, 43, 44, 45, 46)]
                self.SetKeywords(8, self._methodWidgets[8].currentIndex())
                self.SetKeywords(46, 1)
        elif widgetNumber == 2:
            if selection == 0:
                self.keywords_line[1] = ''
            elif selection == 1:
                self.keywords_line[1] = 'R'
            elif selection == 2:
                self.keywords_line[1] = 'U'
            else:
                self.keywords_line[1] = 'RO'
        elif widgetNumber == 3:
            for i in range(2, 11):
                self.keywords_line[i] = ''
            functional = self._dftFunctionals[selection]
            self.keywords_line[2] = f'{functional}/'
            if selection in (0, 1, 6, 8, 9):
                [self._methodWidgets[i].setVisible(True) for i in (28, 29)]
                self._methodWidgets[33].setEnabled(True)
                self.SetKeywords(33, self._methodWidgets[33].isChecked())
                self.SetKeywords(29, self._methodWidgets[29].currentIndex())
            else:
                [self._methodWidgets[i].setVisible(False) for i in (28, 29, 34)]
                self._methodWidgets[33].setEnabled(False)
            self.SetKeywords(8, self._methodWidgets[8].currentIndex())
        elif widgetNumber == 4:
            functional = self._semiempiricalFunctionals[selection]
            self.keywords_line[2] = f'{functional}'
        elif widgetNumber == 5:
            method = self._mechanics[selection]
            self.keywords_line[2] = f'{method}'
            self.SetKeywords(36, self._methodWidgets[36].currentIndex())
        elif widgetNumber == 7:
            if selection == 1:
                self.keywords_line[3] = self._augmented[selection]
            else:
                self.keywords_line[3] = ''
        elif widgetNumber == 8:
            for i in range(3, 9):
                self.keywords_line[i] = ''
            self.keywords_line[4] = self._basisSet[selection]

            if selection == 0:
                [self._methodWidgets[i].setVisible(False) for i in range(7, 16) if i not in (8, 10)]
                [self._methodWidgets[i].setVisible(False) for i in range(16, 28)]
                self._methodWidgets[10].setVisible(True)
                self._methodWidgets[10].setCurrentIndex(0)
                self._methodWidgets[10].model().item(2).setEnabled(False)
                self.SetKeywords(10, self._methodWidgets[10].currentIndex())
            elif selection == 1:
                self._methodWidgets[7].setVisible(False)
                [self._methodWidgets[i].setVisible(True) for i in (9, 10)]
                [self._methodWidgets[i].setVisible(False) for i in range(11, 28)]
                self._methodWidgets[9].setCurrentIndex(0)
                self._methodWidgets[9].model().item(2).setEnabled(False)
                self._methodWidgets[10].model().item(2).setEnabled(True)
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (9, 10)]
            elif selection in range(2, 5):
                self._methodWidgets[7].setVisible(False)
                self._methodWidgets[10].setVisible(False)
                [self._methodWidgets[i].setVisible(False) for i in range(16, 28)]
                [self._methodWidgets[i].setVisible(True) for i in range(9, 16) if i != 10]
                self._methodWidgets[9].model().item(2).setEnabled(True)
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (9, 12)]
            elif selection in range(5, 8):
                [self._methodWidgets[i].setVisible(False) for i in range(9, 28)]
                self._methodWidgets[7].setVisible(True)
                self.SetKeywords(7, self._methodWidgets[7].currentIndex())
            elif selection in range(14, 16):
                [self._methodWidgets[i].setVisible(False) for i in range(9, 26) if i not in (16, 18, 20)]
                [self._methodWidgets[i].setVisible(True) for i in (16, 18, 20, 26, 27)]
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (18, 27)]
                self.SetAddInput()
            else:
                [self._methodWidgets[i].setVisible(False) for i in range(7, 28) if i != 8]
        elif widgetNumber == 9:
            basis = self._basisSet[self._methodWidgets[8].currentIndex()]
            if selection == 1:
                self.keywords_line[4] = basis.replace('G', '+G')
            elif selection == 2:
                self.keywords_line[4] = basis.replace('G', '++G')
            else:
                self.keywords_line[4] = basis
        elif widgetNumber == 10:
            if selection in (1, 2):
                self.keywords_line[5] = self._asterisk[selection]
            else:
                self.keywords_line[5] = ''
        elif widgetNumber in (12, 14):
            first = self._methodWidgets[12].currentIndex()
            second = self._methodWidgets[14].currentIndex()
            if first != 0 or second != 0:
                self.keywords_line[6] = '('
                self.keywords_line[8] = ')'
                if first != 0 and second == 0:
                    self.keywords_line[7] = self._firstPolarizedOrbitals[first]
                elif first == 0 and second != 0:
                    self.keywords_line[7] = self._secondPolarizedOrbitals[second]
                elif first != 0 and second != 0:
                    self.keywords_line[7] = self._firstPolarizedOrbitals[first] + ',' + self._secondPolarizedOrbitals[second]
            else:
                for i in range(6, 9):
                    self.keywords_line[i] = ''
        elif widgetNumber == 17:
            if selection == 1:
                self._lightBasis[0] = self._augmented[selection]
            else:
                self._lightBasis[0] = ''
            self.SetAddInput()
        elif widgetNumber == 18:
            for i in range(0, 6):
                self._lightBasis[i] = ''
            self._lightBasis[1] = self._lightBasisSet[selection]
            self.SetKeywords(27, self._methodWidgets[27].currentIndex())
            if selection == 0:
                [self._methodWidgets[i].setVisible(False) for i in range(17, 26) if i not in (18, 20)]
                self._methodWidgets[20].setVisible(True)
                self._methodWidgets[20].setCurrentIndex(0)
                self._methodWidgets[20].model().item(2).setEnabled(False)
                self.SetKeywords(20, self._methodWidgets[20].currentIndex())
            elif selection == 1:
                self._methodWidgets[17].setVisible(False)
                [self._methodWidgets[i].setVisible(True) for i in (19, 20)]
                [self._methodWidgets[i].setVisible(False) for i in range(21, 26)]
                self._methodWidgets[19].setCurrentIndex(0)
                self._methodWidgets[19].model().item(2).setEnabled(False)
                self._methodWidgets[20].model().item(2).setEnabled(True)
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (19, 20)]
            elif selection in range(2, 5):
                [self._methodWidgets[i].setVisible(False) for i in (17, 20)]
                [self._methodWidgets[i].setVisible(True) for i in range(21, 26) if i != 20]
                self._methodWidgets[19].model().item(2).setEnabled(True)
                [self.SetKeywords(i, self._methodWidgets[i].currentIndex()) for i in (19, 22)]
            elif selection in range(5, 8):
                [self._methodWidgets[i].setVisible(False) for i in range(19, 26)]
                self._methodWidgets[17].setVisible(True)
                self.SetKeywords(17, self._methodWidgets[17].currentIndex())
            self.SetAddInput()
        elif widgetNumber == 19:
            basis = self._lightBasisSet[self._methodWidgets[18].currentIndex()]
            if selection == 1:
                self._lightBasis[1] = basis.replace('G', '+G')
            elif selection == 2:
                self._lightBasis[1] = basis.replace('G', '++G')
            else:
                self._lightBasis[1] = basis
            self.SetAddInput()
        elif widgetNumber == 20:
            if selection in (1, 2):
                self._lightBasis[2] = self._asterisk[selection]
            else:
                self._lightBasis[2] = ''
            self.SetAddInput()
        elif widgetNumber in (22, 24):
            first = self._methodWidgets[22].currentIndex()
            second = self._methodWidgets[24].currentIndex()
            if first != 0 or second != 0:
                self._lightBasis[3] = '('
                self._lightBasis[5] = ')'
                if first != 0 and second == 0:
                    self._lightBasis[4] = self._firstPolarizedOrbitals[first]
                elif first == 0 and second != 0:
                    self._lightBasis[4] = self._secondPolarizedOrbitals[second]
                elif first != 0 and second != 0:
                    self._lightBasis[4] = self._firstPolarizedOrbitals[first] + ',' + self._secondPolarizedOrbitals[second]
            else:
                for i in range(3, 6):
                    self._lightBasis[i] = ''
            self.SetAddInput()
        elif widgetNumber == 27:
            self._heavyBasis = self._heavyBasisSet[selection]
            self.SetAddInput()
        elif widgetNumber == 29:
            if selection == 0:
                self.keywords_line[9] = ''
            elif selection == 1:
                self.keywords_line[9] = '/FIT'
            elif selection == 3:
                self._methodWidgets[34].setVisible(True)
                self.keywords_line[9] = '/' + self._fittingSet[selection].replace('N', str(self._methodWidgets[34].value()))
            else:
                self.keywords_line[9] = '/' + self._fittingSet[selection]
                self._methodWidgets[34].setVisible(False)
        elif widgetNumber == 33:
            if selection:
                self.keywords_line[10] = ' SPARSE'
            else:
                self.keywords_line[10] = ''
        elif widgetNumber == 34:
            self.keywords_line[9] = re.sub(r'\d+', str(selection), self.keywords_line[9])
        elif widgetNumber == 36:
            if selection == 0:
                self.keywords_line[3] = ''
            elif selection == 1:
                self.keywords_line[3] = '=QEQ'
            elif selection == 2:
                self.keywords_line[3] = '=UNTYPED'
            else:
                self.keywords_line[3] = '=UNCHARGED'
        elif widgetNumber == 37:
            selectedMethod = self._methodWidgets[1].currentIndex()
            if selectedMethod == 4:
                if selection:
                    self.keywords_line[3] = '=FULL/'
                else:
                    self.keywords_line[3] = '/'
            elif selectedMethod == 5:
                if selection:
                    self.keywords_line[3] = ',FULL)/'
                else:
                    self.keywords_line[3] = ')/'
            elif selectedMethod == 6:
                self.SetKeywords(39, self._methodWidgets[39].isChecked())
        elif widgetNumber == 38:
            if selection == 0:
                self.keywords_line[2] = 'MP4(SDTQ'
            else:
                self.keywords_line[2] = 'MP4(SDQ'
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
                    self.keywords_line[3] = '=' + re.sub(r'\,$', '', ''.join(variables)) + '/'
                elif first == 0 and second == 0:
                    self.keywords_line[3] = '(' + re.sub(r'\,$', '', ''.join(variables)) + ')/'
                else:
                    self.keywords_line[3] = '=(' + re.sub(r'\,$', '', ''.join(variables)) + ')/'
            else:
                self.keywords_line[3] = '/'
        elif widgetNumber in range(44, 47):
            self.keywords_line[2] = f'CASSCF({self._methodWidgets[44].value()},{self._methodWidgets[45].value()}'
            if self._methodWidgets[46].isChecked():
                self.keywords_line[3] = ',RFO)/'
            else:
                self.keywords_line[3] = ')/'
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
                    self.keywords_line[3] = '=(' + re.sub(r'\,$', '', ''.join(variables)) + ')/'
                else:
                    self.keywords_line[3] = '=' + re.sub(r'\,$', '', ''.join(variables)) + '/'
            else:
                self.keywords_line[3] = '/'
            if self._methodWidgets[0].currentIndex() == 1:
                self.keywords_line[3] = re.sub(r'/$', '', self.keywords_line[3])
        elif widgetNumber == 60:
            self.keywords_line[-1] = ' ' + selection.upper()

        self.uiKeywordsLabel.setText(''.join(self.keywords_line))
        self.set_preview()

    def SetJobType(self, widgetNumber, selection):

        if widgetNumber == 0:
            if selection == 0:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(False)
                self._methodWidgets[0].model().item(1).setEnabled(True)
                for i in range(-4, -1):
                    self.keywords_line[i] = ''
            elif selection == 1:
                self.uiOptGroupBox.setVisible(True)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(False)
                self._methodWidgets[0].setCurrentIndex(0)
                self._methodWidgets[0].model().item(1).setEnabled(False)
                self.keywords_line[-3] = ''
                self.keywords_line[-2] = ''
                self.SetJobType(1, 1)
            elif selection == 2:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(True)
                self.uiIRCGroupBox.setVisible(False)
                self._methodWidgets[0].model().item(1).setEnabled(True)
                self.keywords_line[-4] = ''
                self.SetJobType(5, 1)
            elif selection == 3:
                self.uiOptGroupBox.setVisible(True)
                self.uiFreqGroupBox.setVisible(True)
                self.uiIRCGroupBox.setVisible(False)
                self._methodWidgets[0].setCurrentIndex(0)
                self._methodWidgets[0].model().item(1).setEnabled(False)
                [self.SetJobType(i, 1) for i in (1, 5)]
            elif selection == 4:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(True)
                self._methodWidgets[0].setCurrentIndex(0)
                self._methodWidgets[0].model().item(1).setEnabled(False)
                for i in range(-4, -1):
                    self.keywords_line[i] = ''
                self.SetJobType(15, 1)
            else:
                self.uiOptGroupBox.setVisible(False)
                self.uiFreqGroupBox.setVisible(False)
                self.uiIRCGroupBox.setVisible(False)
        elif widgetNumber in range(1, 5):
            first = self._jobsWidgets[1].currentIndex()
            if first == 1:
                firstVariable = 'TS,'
            else:
                firstVariable = ''
            second = self._jobsWidgets[2].currentIndex()
            if second == 1:
                secondVariable = 'CALCFC,'
            elif second == 2:
                secondVariable = 'CALCALL,'
            elif second == 3:
                secondVariable = 'READFC,'
            elif second == 4:
                secondVariable = 'RCFC,'
            else:
                secondVariable = ''
            third = self._jobsWidgets[3].isChecked()
            if third:
                thirdVariable = 'RFO,'
            else:
                thirdVariable = ''
            fourth = self._jobsWidgets[4].isChecked()
            if fourth:
                fourthVariable = 'TIGHT,'
            else:
                fourthVariable = ''
            variables = (secondVariable, fourthVariable, thirdVariable, firstVariable)
            combination = (first, second, third, fourth)
            if sum(combination) > 0:
                if sum(combination) > 1:
                    self.keywords_line[-4] = ' OPT=(' + re.sub(r'\,$', '', ''.join(variables)) + ')'
                else:
                    self.keywords_line[-4] = ' OPT=' + re.sub(r'\,$', '', ''.join(variables))
            else:
                self.keywords_line[-4] = ' OPT'
        elif widgetNumber in range(5, 13):
            first = self._jobsWidgets[5].currentIndex()
            if first == 1:
                firstVariable = 'NORAMAN,'
            elif first == 2:
                if self._jobsWidgets[6].currentIndex() == 0:
                    firstVariable = 'RAMAN,'
                else:
                    firstVariable = ''
            else:
                firstVariable = ''
            second = self._jobsWidgets[6].currentIndex()
            if second == 1:
                secondVariable = 'ROA,'
            else:
                secondVariable = ''
            third = self._jobsWidgets[7].currentIndex()
            if third == 1:
                thirdVariable = 'SAVENORMALMODES,'
            elif third == 2:
                thirdVariable = 'NOSAVENORMALMODES,'
            else:
                thirdVariable = ''
            fourth = self._jobsWidgets[8].isChecked()
            if fourth:
                fourthVariable = 'VCD,'
            else:
                fourthVariable = ''
            fifth = self._jobsWidgets[9].isChecked()
            if fifth:
                fifthVariable = 'NODIAGFULL,'
            else:
                fifthVariable = ''
            sixth = self._jobsWidgets[10].isChecked()
            if sixth:
                sixthVariable = 'ANHARMONIC,'
                self._jobsWidgets[11].setEnabled(True)
            else:
                sixthVariable = ''
                self._jobsWidgets[11].setEnabled(False)
            seventh = self._jobsWidgets[11].isChecked()
            if seventh:
                seventhVariable = 'SELECTANHARMONICMODES,'
                self._jobsWidgets[13].setEnabled(True)
                self.SetJobType(13, self._jobsWidgets[13].value())
            else:
                seventhVariable = ''
                self._jobsWidgets[13].setEnabled(False)
            eighth = self._jobsWidgets[12].isChecked()
            if eighth:
                eighthVariable = 'PROJECTED,'
            else:
                eighthVariable = ''
            variables = (
                firstVariable, secondVariable, fourthVariable, fifthVariable,
                thirdVariable, sixthVariable, seventhVariable, eighthVariable
            )
            combination = (first, second, third, fourth, fifth, sixth, seventh, eighth)
            if sum(combination) > 0:
                if sum(combination) > 1:
                    self.keywords_line[-3] = ' FREQ=(' + re.sub(r'\,$', '', ''.join(variables)) + ')'
                else:
                    self.keywords_line[-3] = ' FREQ=' + re.sub(r'\,$', '', ''.join(variables))
            else:
                self.keywords_line[-3] = ' FREQ'
        elif widgetNumber == 14:
            if selection == 1:
                self.keywords_line[-2] = ' CPHF=RDFREQ'
            elif selection == 2:
                self.keywords_line[-2] = ' CPHF=NOREAD'
            else:
                self.keywords_line[-2] = ''
        elif widgetNumber in range(15, 23):
            first = self._jobsWidgets[15].currentIndex()
            if first == 1:
                firstVariable = 'FORWARD,'
            elif first == 2:
                firstVariable = 'REVERSE,'
            else:
                firstVariable = ''
            second = self._jobsWidgets[16].currentIndex()
            if second == 1:
                secondVariable = 'CALCALL,'
                self._jobsWidgets[20].setEnabled(False)
            elif second == 2:
                secondVariable = 'RCFC,'
                self._jobsWidgets[20].setEnabled(True)
                self._jobsWidgets[21].setEnabled(True)
            else:
                secondVariable = 'CALCFC,'
                self._jobsWidgets[20].setEnabled(True)
                self._jobsWidgets[21].setEnabled(True)
            third = self._jobsWidgets[17].currentIndex()
            if third == 1:
                thirdVariable = 'RECORRECT=NEVER,'
            elif third == 2:
                thirdVariable = 'RECORRECT=YES,'
            elif third == 3:
                thirdVariable = 'RECORRECT=ALWAYS,'
            elif third == 4:
                thirdVariable = 'RECORRECT=TEST,'
            else:
                thirdVariable = ''
            fourth = self._jobsWidgets[18].isChecked()
            if fourth:
                self._jobsWidgets[19].setEnabled(True)
                fourthVariable = f'MAXPOINTS={str(self._jobsWidgets[19].value())},'
            else:
                fourthVariable = ''
                self._jobsWidgets[19].setEnabled(False)
            fifth = self._jobsWidgets[20].isChecked()
            if fifth:
                self._jobsWidgets[21].setEnabled(True)
                self._jobsWidgets[16].model().item(1).setEnabled(False)
                fifthVariable = f'RECALC={str(self._jobsWidgets[21].value())},'
            else:
                self._jobsWidgets[21].setEnabled(False)
                self._jobsWidgets[16].model().item(1).setEnabled(True)
                fifthVariable = ''
            sixth = self._jobsWidgets[22].isChecked()
            if sixth:
                sixthVariable = 'HF/3-21G:HF/3-21G,'
                for i in range(1, 13):
                    self.keywords_line[i] = ''
                self.keywords_line[-4] = ' IRCMAX='
            else:
                sixthVariable = ''
                self.keywords_line[-4] = ' IRC='
                self.SetKeywords(0, self._methodWidgets[0].currentIndex())
            variables = (
                firstVariable, secondVariable, fourthVariable,
                fifthVariable, thirdVariable, sixthVariable
            )
            combination = (first, second, third, fourth, fifth, sixth)
            if sum(combination) >= 1:
                self.keywords_line[-3] = '(' + re.sub(r'\,$', '', ''.join(variables)) + ')'
            else:
                self.keywords_line[-3] = 'CALCFC'

        self.uiKeywordsLabel.setText(''.join(self.keywords_line))
        self.set_preview()

    def set_title(self, title):
        self.uiTitleLabel.setText(title)
        # self.set_preview()

    def SetCahrgeMult(self):
        """
        docstring
        """
        self._chargeMultiplicityLine[0] = str(self._methodWidgets[30].value())

        self._currentElectrons = self._molecule.get_valence_electrons - self._methodWidgets[30].value() + self._molecule.get_formal_charge
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
        self.set_preview()

    def SetLink0(self, widgetNumber, selection):
        """
        docstring
        """
        if widgetNumber == 0:
            if selection:
                self._link0Line[2] = f'%Mem={str(selection)}GB\n'
            else:
                self._link0Line[2] = ''
        elif widgetNumber == 1:
            if selection:
                self._link0Line[1] = f'%NProcShared={str(selection)}\n'
            else:
                self._link0Line[1] = ''
        elif widgetNumber == 2:
            if selection == 1:
                [self.link_0_widgets[i].setEnabled(False) for i in (3, 4)]
                self.link_0_widgets[3].setText(self._molecule.get_name)
            elif selection == 2:
                [self.link_0_widgets[i].setEnabled(True) for i in (3, 4)]
                self.link_0_widgets[3].clear()
                self.link_0_widgets[3].setPlaceholderText('Name_here')
            else:
                [self.link_0_widgets[i].setEnabled(False) for i in (3, 4)]
                self.link_0_widgets[3].clear()
            self.SetLink0(3, self.link_0_widgets[3].text())
        elif widgetNumber == 3:
            if selection == '':
                self._link0Line[3] = ''
            else:
                self._link0Line[3] = f'%CHK={selection}.chk\n'
        elif widgetNumber == 4:
            self._chk, _ = qtw.QFileDialog.getSaveFileName(
                self,
                'Save your chk file',
                options=qtw.QFileDialog.DontUseNativeDialog,
                filter='Chk files(*.chk)'
            )
            if self._chk:
                self.link_0_widgets[3].setText(f'{self._chk}.chk')
                self.SetLink0(3, self._chk)
        elif widgetNumber == 5:
            if selection == 1:
                [self.link_0_widgets[i].setEnabled(False) for i in (6, 7)]
                self.link_0_widgets[6].setText(f'{self._molecule.get_name}')
            elif selection == 2:
                [self.link_0_widgets[i].setEnabled(True) for i in (6, 7)]
                self.link_0_widgets[6].clear()
                self.link_0_widgets[6].setPlaceholderText('Name_here')
            else:
                [self.link_0_widgets[i].setEnabled(False) for i in (6, 7)]
                self.link_0_widgets[6].clear()
            self.SetLink0(6, self.link_0_widgets[6].text())
        elif widgetNumber == 6:
            if selection == '':
                self._link0Line[0] = ''
            else:
                self._link0Line[0] = f'%OLDCHK={selection}.chk\n'
        elif widgetNumber == 7:
            self._oldChk, _ = qtw.QFileDialog.getOpenFileName(
                self, 'Open an old chk file',
                options=qtw.QFileDialog.DontUseNativeDialog,
                filter='Chk files(*.chk)'
            )
            if self._oldChk:
                self.link_0_widgets[6].setText(f'{self._oldChk}.chk')
                self.SetLink0(6, self._oldChk)

        self.link_0_widgets[8].setPlainText(''.join(self._link0Line))
        self.set_preview()

    def SetGeneral(self, widgetNumber, selection):
        """
        docstring
        """
        if widgetNumber == 0:
            if selection:
                self.keywords_line[11] = ' SCF=QC'
            else:
                self.keywords_line[11] = ''
        if widgetNumber == 1:
            if selection:
                self.keywords_line[12] = ' NOSYMM'
            else:
                self.keywords_line[12] = ''
        if widgetNumber == 2:
            if selection:
                self.coordinates = self._molecule.get_coordinates
            else:
                self.coordinates = self._molecule.get_zmatrix
        if widgetNumber == 7:
            if selection:
                self.keywords_line[0] = '#P '
            else:
                self.keywords_line[0] = '# '
        if widgetNumber == 9:
            if selection:
                self.keywords_line[13] = ' POLAR'
            else:
                self.keywords_line[13] = ''
        if widgetNumber == 11:
            if selection:
                self._generalWidgets[12].setEnabled(True)
                self.SetGeneral(12, self._generalWidgets[12].value())
            else:
                self._generalWidgets[12].setEnabled(False)
                self.keywords_line[14] = ''
        if widgetNumber == 12:
            self.keywords_line[14] = f' MAXDISK={str(selection)}GB'

        self.uiKeywordsLabel.setText(''.join(self.keywords_line))
        self.set_preview()

    def SetAddInput(self):
        """
        docstring
        """
        heavy_atoms = ' '.join(self._molecule.get_heavy_atoms)
        light_atoms = ' '.join(self._molecule.get_light_atoms)
        self.uiAddInputPlainText.clear()
        self.uiAddInputPlainText.setPlainText(
            f"{heavy_atoms} 0\n{self._heavyBasis}\n****\n{light_atoms} 0\n{''.join(self._lightBasis)}\n****\n\n{heavy_atoms} 0\n{self._heavyBasis}"
        )

    def set_preview(self):
        """
        docstring
        """
        if self.uiAddInputCheckBox.isChecked():
            self.additional_input = '\n' + self.uiAddInputPlainText.toPlainText()
        else:
            self.additional_input = ''

        self.input = [
            self.link_0_widgets[8].toPlainText(),
            self.uiKeywordsLabel.text(),
            self.uiTitleLineEdit.text(),
            self.uiChargeMultLabel.text(),
            self.coordinates,
            self.additional_input
        ]
        self.uiPreviewPlainText.setPlainText(
            '{}{}\n\n {}\n\n{}\n{}\n{}\n'.format(*self.input)
        )

    def queue_calculation(self):
        """
        docstring
        """
        if not self.job_type:
            self.job_type = self._jobTypes[self._jobsWidgets[0].currentIndex()]
        name = f'molecules/{self._molecule.inchi_key}/{self.job_type}'
        existent_imputs = str(
            len(
                [i for i in listdir(f'molecules/{self._molecule.inchi_key}/') if re.sub(r'_\d+\.com', '', i) == self.job_type]
            )
        )
        input_file = f'{name}_{existent_imputs}.com'
        keywords = self.uiKeywordsLabel.text()
        charge_mult = self.uiChargeMultLabel.text().replace(' ', '/')
        # error por solucionar
        if not self.group:
            with open(input_file, 'w') as f:
                f.write(self.uiPreviewPlainText.toPlainText())
        output_file = input_file.replace('.com', '.log')
        calculation = {
            'type': self.job_type,
            'keywords': keywords,
            'charge_mult': charge_mult,
            'input_file': input_file,
            'output_file': output_file,
            'molecule': self._molecule.get_name,
            'molecule_id': self._molecule.inchi_key,
            'coordinates': self.coordinates,
            'status': 'Pending'
        }
        self.submitted.emit(calculation, self.project.name)
        self.close()


# class GaussianWorker():

    # def run(self, input_file, output_file):
        # # pasar al módulo Calculations, formando una nueva clase llamada
        # # GaussianWorker, que en la función run recibe el input y output y
        # # retorna si el cálculo fue exitoso o no.
        # # Primero debe hacer checkeo a ver si ya existe ese output y si ya el
        # # cálculo está terminado. Si es así, debe retornar True de una vez
        # '''
        # Parameters
        # ----------
        # input_file: the .com file as Gaussian input
        # output_file: the name of the .log file to be done

        # Returns
        # -------
        # True if the input file is successfully run in the shell
        # False if something went wrong with the gaussian executable
        # '''
        # env = os.environ.copy()
        # try:
            # gauss_exec = env['GAUSS_EXEDIR'].split('/')[-1]
        # except KeyError:
            # return False
            # print('no leyó el ejecutable de Gaussian')
        # cmd = f'{gauss_exec} < {input_file} > {output_file}'
        # subprocess.run(cmd, shell=True)
        # return True

    # def check(self, output_file):
        # if not os.path.exists(output_file):
            # return False
        # if 'Optimization' in output_file:
            # with open(output_file, 'r') as f:
                # file = f.read()
                # finished = re.findall('Optimization completed.', file)
                # if finished:
                    # return True
                # return False


class DockingPlotter(qtw.QWidget):

    # @logger.catch
    def __init__(self, projects, jobs):
        super().__init__()
        uic.loadUi('Views/uiDockingResults.ui', self)
        self.projects = projects
        self.jobs = jobs
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        self.selected_jobs = []
        self.selected_projects = []
        self.create_checks()
        self.create_df()
        self.show()

    def create_checks(self):
        self.checks = [
            qtw.QCheckBox(f'{project.name}') for project in self.projects
        ]
        for i, check in enumerate(self.checks):
            check.stateChanged.connect(self.select_project)
            self.uiProjectsLayout.addWidget(check, i)
        spacer = qtw.QSpacerItem(20, 40, qtw.QSizePolicy.Minimum, qtw.QSizePolicy.Expanding)
        self.uiProjectsLayout.addItem(spacer)

    def set_model(self):
        model = PandasModel(self.df)
        self.uiResultsTableView.setModel(model)
        self.uiResultsTableView.resizeColumnsToContents()
        self.uiResultsTableView.resizeRowsToContents()
        self.plot()

    def select_project(self):
        indexes = [p.isChecked() for p in self.checks]
        self.selected_projects = [
            p for check, p in zip(indexes, self.projects) if check
        ]
        self.create_df()

    def select_jobs(self):
        indexes = self.uiResultsTableView.selectedIndexes()
        if indexes:
            self.selected_jobs = [i.row() for i in indexes]
            self.plot()
        else:
            self.selected_jobs = []

    def create_df(self):
        dfs = []
        for project in self.selected_projects:
            if not project.job_ids:
                continue
            names = pd.Series(
                [mol.get_name for mol in project.molecules], name='Names'
            )
            energies = []
            names = []
            receptor = []
            box_size = []
            project_name = []
            for job in self.jobs:
                if job.id in project.job_ids:
                    if len(job.energies) == 0:
                        continue
                    energies.append(job.energies)
                    names.append(job.molecule)
                    receptor.append(job.receptor_name)
                    size_x = job.config.get('size_x')
                    size_y = job.config.get('size_y')
                    size_z = job.config.get('size_z')
                    box_size.append(f'{size_x}x{size_y}x{size_z}')
                    project_name.append(project.name)
            values = {
                'Molecule': names,
                'Receptor': receptor,
                'Box size (\u212B\u00B3)': box_size,
                'Binding energies': list(energies),
                'Project': project_name
            }
            dfs.append(pd.DataFrame(values))
        self.df = pd.concat(dfs) if dfs else pd.DataFrame([])
        self.set_model()

    def plot(self):
        self.canvas.ax.clear()
        if not self.selected_jobs or not self.df.shape[0]:
            self.uiSortedCheck.setChecked(False)
            self.uiCompareProjectButton.setEnabled(False)
            self.uiCompareReceptorButton.setEnabled(False)
            return
        self.sorted = self.uiSortedCheck.isChecked()
        self.uiSortedCheck.setEnabled(True)
        self.canvas.ax.set_title('Docking')
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.axhline(-7., linewidth=.5, color='red')
        to_plot_both = self.df.iloc[self.selected_jobs]
        receptors = pd.unique(to_plot_both['Receptor'])
        if len(receptors) > 1:
            self.uiCompareReceptorButton.setEnabled(True)
        else:
            self.uiCompareReceptorButton.setEnabled(False)
        projects = pd.unique(to_plot_both['Project'])
        if len(projects) > 1:
            self.uiCompareProjectButton.setEnabled(True)
        else:
            self.uiCompareProjectButton.setEnabled(False)
        total = to_plot_both.shape[0]
        if self.sorted:
            to_plot_both['Medians'] = to_plot_both['Binding energies'].apply(np.median)
            to_plot_both.sort_values(by='Medians', inplace=True)
        props={'color': 'blue', 'linewidth': 1.0}
        positions = np.arange(total)
        self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
        self.canvas.ax.boxplot(
            positions=positions,
            x=to_plot_both['Binding energies'],
            labels=to_plot_both['Molecule'],
            boxprops=props, medianprops=props,
            whiskerprops=props, capprops=props,
        )
        self.canvas.draw()

    def compare_project(self):
        self.canvas.ax.clear()
        self.canvas.ax.set_title('Docking')
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.axhline(-7., linewidth=.5, color='red')
        colors = ['blue', 'red', 'green', 'magenta', 'orange']
        to_plot_both = self.df.iloc[self.selected_jobs]
        projects = pd.unique(to_plot_both['Project'])
        if len(projects) != 2:
            return # here a critical warning
        total = to_plot_both.shape[0]
        both_df = []
        for i, project in enumerate(projects):
            mask = to_plot_both['Project'] == project
            to_plot = to_plot_both.loc[mask, ['Binding energies', 'Molecule']]
            both_df.append(to_plot)
        molecules_p1 = both_df[0]['Molecule'].values.tolist()
        molecules_p2 = both_df[1]['Molecule'].values.tolist()
        if both_df[0].shape[0] != both_df[1].shape[0]:
            return # avisar que no tienen el mismo número de compuestos
        elif molecules_p1 != molecules_p2:
            return # avisar que no están en el mismo orden
        for i, project in enumerate(projects):
            props = {'color': colors[i], 'linewidth': 1.0}
            total = to_plot.shape[0]
            positions = np.arange(total) * 2 - 0.4 * (-1) ** (i + 1)
            self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
            self.canvas.ax.plot([], c=colors[i], label=projects[i])
            self.canvas.ax.boxplot(
                positions=positions,
                x=both_df[i].loc[:, 'Binding energies'],
                labels=molecules_p1,
                boxprops=props, medianprops=props,
                whiskerprops=props, capprops=props,
            )
            self.canvas.ax.set_xticks(np.arange(total) * 2)
        self.canvas.ax.legend()
        self.canvas.draw()
        x1 = both_df[0]['Binding energies'].apply(np.median).values.tolist()
        x2 = both_df[1]['Binding energies'].apply(np.median).values.tolist()
        plt.figure(figsize=(10,10))
        plt.scatter(x1, x2, alpha=0.3)
        plt.title('Linear correlation between projects')
        plt.xlabel(projects[0])
        plt.ylabel(projects[1])
        fit = np.polyfit(x1, x2, 1)
        f = np.poly1d(fit)
        m, b = fit
        r2 = r2_score(x2, f(x1))
        x = np.linspace(min(x1), max(x1), num=100)
        plt.plot(
            x, f(x), c='orange', label='Trend line', linewidth=0.5
        )
        plt.annotate(
            f'r\u00B2 = {r2:.5f}\nm = {m:.5f}\nb = {b:.5f}', (min(x1), max(x2) - 0.5)
        )
        plt.legend(loc='lower right')
        plt.show()

    def compare_receptor(self):
        self.canvas.ax.clear()
        self.canvas.ax.set_title('Docking')
        self.canvas.ax.set_xlabel('Molecule')
        self.canvas.ax.set_ylabel('Binding Energy (kcal/mol)')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.axhline(-7., linewidth=.5, color='red')
        colors = ['blue', 'red', 'green', 'magenta', 'orange']
        to_plot_both = self.df.iloc[self.selected_jobs]
        receptors = pd.unique(to_plot_both['Receptor'])
        if len(receptors) != 2:
            return # here a critical warning
        total = to_plot_both.shape[0]
        both_df = []
        for i, receptor in enumerate(receptors):
            mask = to_plot_both['Receptor'] == receptor
            to_plot = to_plot_both.loc[mask, ['Binding energies', 'Molecule']]
            to_plot.reset_index(drop=True, inplace=True)
            both_df.append(to_plot)
        molecules_p1 = both_df[0]['Molecule'].values.tolist()
        molecules_p2 = both_df[1]['Molecule'].values.tolist()
        if both_df[0].shape[0] != both_df[1].shape[0]:
            return # avisar que no tienen el mismo número de compuestos
        elif molecules_p1 != molecules_p2:
            return # avisar que no están en el mismo orden
        for i, receptor in enumerate(receptors):
            props = {'color': colors[i], 'linewidth': 1.0}
            total = both_df[i].shape[0]
            positions = np.arange(total) * 2 - 0.4 * (-1) ** (i + 1)
            self.canvas.ax.set_xlim(left=-2, right=positions[-1] + 2)
            self.canvas.ax.plot([], c=colors[i], label=receptors[i])
            self.canvas.ax.boxplot(
                positions=positions,
                x=both_df[i].loc[:, 'Binding energies'],
                labels=molecules_p1,
                boxprops=props, medianprops=props,
                whiskerprops=props, capprops=props,
            )
            self.canvas.ax.set_xticks(np.arange(total) * 2)
        self.canvas.ax.legend()
        self.canvas.draw()
        # if working with all data is wanted
        # both_df = [df['Binding energies'].apply(np.sort) for df in both_df]
        # x1 = np.concatenate(both_df[0].values.tolist())
        # x2 = np.concatenate(both_df[1].values.tolist())
        # if working with the medians of each group of data is wanted
        x1 = both_df[0]['Binding energies'].apply(np.median).values.tolist()
        x2 = both_df[1]['Binding energies'].apply(np.median).values.tolist()
        plt.figure(figsize=(5,5))
        plt.scatter(x1, x2, alpha=0.3, label='Medians by molecule')
        plt.title('Linear correlation between receptors')
        plt.xlabel(receptors[0])
        plt.ylabel(receptors[1])
        fit = np.polyfit(x1, x2, 1)
        f = np.poly1d(fit)
        m, b = fit
        r2 = r2_score(x2, f(x1))
        x = np.linspace(min(x1), max(x1), num=100)
        plt.plot(
            x, f(x), c='orange', label='Trend line', linewidth=0.5
        )
        plt.annotate(
            f'r\u00B2 = {r2:.5f}\nm = {m:.5f}\nb = {b:.5f}', (min(x1), max(x2) - 0.5)
        )
        plt.legend(loc='lower right')
        plt.show()

    def export_to_excel(self):
        filename, _ = qtw.QFileDialog.getSaveFileName(
            self, 'Save results as excel data sheet',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Excel file (*.xlsx *.xls)'
        )
        if filename:
            if '.xlsx' not in filename or '.xls' not in filename:
                self.df.to_excel(filename + '.xlsx')
            else:
                self.df.to_excel(filename)

    def import_from_excel(self):
        filename, _ = qtw.QFileDialog.getOpenFileName(
            self, 'Save results as excel data sheet',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Excel file (*.xlsx *xls)'
        )
        if filename:
            new_df = pd.read_excel(filename)
            new_df['Binding energies'] = new_df['Binding energies'].apply(
                lambda x: np.array(
                    [float(i) for i in x.replace('[', '').replace(']', '').split()]
            ))
            new_df = new_df.iloc[:, 1:]
            self.df = pd.concat([self.df, new_df], ignore_index=True)
            self.set_model()


class RedockingPlotter(qtw.QWidget):

    def __init__(self, projects, jobs):
        super().__init__()
        uic.loadUi('Views/uiRedockingResults.ui', self)
        self.projects = projects
        self.jobs = jobs
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        self.selected_jobs = []
        self.selected_projects = []
        self.create_checks()
        self.create_df()
        self.show()

    def create_checks(self):
        self.checks = [
            qtw.QCheckBox(f'{project.name}') for project in self.projects
        ]
        for i, check in enumerate(self.checks):
            check.stateChanged.connect(self.select_project)
            self.uiProjectsLayout.addWidget(check, i)
        spacer = qtw.QSpacerItem(20, 40, qtw.QSizePolicy.Minimum, qtw.QSizePolicy.Expanding)
        self.uiProjectsLayout.addItem(spacer)

    def create_df(self):
        dfs = []
        for project in self.selected_projects:
            if not project.job_ids:
                continue
            names = pd.Series(
                [mol.get_name for mol in project.molecules], name='Names'
            )
            energies = []
            names = []
            receptor = []
            box_size = []
            rmsd = []
            project_name = []
            for job in self.jobs:
                if job.id in project.job_ids:
                    if len(job.energies) == 0:
                        continue
                    energies.append(job.energies)
                    names.append(job.molecule)
                    receptor.append(job.receptor_name)
                    size_x = job.config.get('size_x')
                    size_y = job.config.get('size_y')
                    size_z = job.config.get('size_z')
                    box_size.append(f'{size_x}x{size_y}x{size_z}')
                    project_name.append(project.name)
                    rmsd.append(job.rmsd)
            values = {
                'Molecule': names,
                'Receptor': receptor,
                'Box size (\u212B\u00B3)': box_size,
                'Binding energies': list(energies),
                'Project': project_name,
                'RMSD': rmsd
            }
            dfs.append(pd.DataFrame(values))
        self.df = pd.concat(dfs) if dfs else pd.DataFrame([])
        self.set_model()

    def set_model(self):
        model = PandasModel(self.df)
        self.uiResultsTableView.setModel(model)
        self.uiResultsTableView.resizeColumnsToContents()
        self.uiResultsTableView.resizeRowsToContents()
        self.plot()

    def plot(self):
        self.canvas.ax.clear()
        if not self.selected_jobs or not self.df.shape[0]:
            return
        self.ax2 = self.canvas.ax.twinx()
        self.ax2.clear()
        self.canvas.ax.set_ylabel('RMSD')
        self.canvas.ax.tick_params(axis='x', labelrotation=90)
        self.canvas.ax.set_ylim(bottom=0., top=2.5)
        props={'color': 'blue', 'linewidth': 1.0}
        to_plot = self.df.iloc[self.selected_jobs]
        self.canvas.ax.set_xlabel(to_plot['Molecule'])
        receptor_name = to_plot.loc[0, 'Receptor']
        self.canvas.ax.set_title(f'Docking on {receptor_name}')
        self.canvas.ax.boxplot(
            x=to_plot['RMSD'], positions=[-0.5],
            # labels=to_plot['Molecule'],
            boxprops=props, medianprops=props,
            whiskerprops=props, capprops=props,
        )
        self.ax2.set_ylabel('Binding Energy (kcal/mol)')
        props={'color': 'red', 'linewidth': 1.0}
        self.ax2.boxplot(
            x=to_plot['Binding energies'], positions=[0.5],
            labels=to_plot['Molecule'],
            boxprops=props, medianprops=props,
            whiskerprops=props, capprops=props
        )
        self.canvas.draw()

    def select_jobs(self):
        indexes = self.uiResultsTableView.selectedIndexes()
        if indexes:
            self.selected_jobs = [i.row() for i in indexes]
            self.plot()
        else:
            self.selected_jobs = []

    def select_project(self):
        indexes = [p.isChecked() for p in self.checks]
        self.selected_projects = [
            p for check, p in zip(indexes, self.projects) if check
        ]
        self.create_df()


class MyVina(qtw.QWidget):

    submitted = qtc.pyqtSignal(dict, str)

    def __init__(self, molecule, project=False, calc=False):
        super().__init__()
        uic.loadUi('Views/uiVina.ui', self)
        self.molecule = molecule
        self.project = project
        self.config = {}
        self.redocking = False
        self.nat_lig_path = False
        items = list(molecule.conf_dict.keys())
        if project and project.calculations:
            opt = [
                calc.get('keywords') for calc in project.calculations\
                   if calc.get('type') == 'Optimization'
            ]
            items.extend(opt)
            self.uiConformerCombo.addItems(items)
        self.set_conformer()
        self.set_auto_box_size()
        if calc:
            self.set_options(calc)
        else:
            self.show()

    def set_options(self, calc):
        calc = copy.deepcopy(calc)
        self.config = calc.get('config')
        self.times = calc.get('times')
        self.nat_lig_path = calc.get('nat_lig_path')
        self.redocking = calc.get('redocking')
        self.uiAutoBoxSize.setChecked(calc.get('auto_box_size'))
        self.conformer = calc.get('conformer')
        self.set_auto_box_size()

    def set_conformer(self):
        conf = self.uiConformerCombo.currentText()
        self.conformer = conf if conf else 'input'

    def set_auto_box_size(self):
        # according to a Rg to box size ratio of 0.35
        self.auto_box_size = self.uiAutoBoxSize.isChecked()
        if self.auto_box_size:
            box_size = 2.857 * self.molecule.descriptors.loc[0, 'Rg']
        else:
            box_size = 30.
        self.uiSizeXDouble.setValue(box_size)
        self.uiSizeYDouble.setValue(box_size)
        self.uiSizeZDouble.setValue(box_size)

    def open_receptor(self):
        rec_path, _ = qtw.QFileDialog.getOpenFileName(
            self, 'Selecciona el receptor',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Archivos pdb o pdbqt (*.pdb *.pdbqt)'
        )
        if rec_path:
            self.uiReceptorLine.setText(rec_path)
            self.uiSubmitButton.setEnabled(True)

    def receptor_prep(self, receptor_in, receptor_out):
        if not os.path.exists(receptor_out):
            # water suppression
            with open(receptor_in, 'r') as file:
                lines = file.readlines()
                new_lines = []
                for line in lines:
                    if line.startswith('HETATM'):
                        col = line.split()
                        if col[2] != 'O' and col[3] != 'HOH':
                            lines.append(line)
                    else:
                        new_lines.append(line)
            with open(receptor_in, 'w') as file:
                file.writelines(new_lines)
            # mol2 file generation in order to get the charges
            mol2_file = receptor_out.split(".")[0] + '.mol2'
            subprocess.run(
                f'obabel {receptor_in} -O {mol2_file} -h',
                shell=True,
                capture_output=True,
                text=True
            )
            # pdbqt file generation  of a rigid prot. with proper Hs at pH = 7.4
            subprocess.run(
                f'obabel {mol2_file} -O {receptor_out} -p 7.4 -xr',
                shell=True,
                capture_output=True,
                text=True
            )
            os.remove(mol2_file)

    def set_redocking(self):
        self.redocking = self.uiRedockingCheck.isChecked()

    def center_from_nat_ligand(self):
        nat_lig_path, _ = qtw.QFileDialog.getOpenFileName(
            self, 'Selecciona el ligando natural',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Archivos pdb o pdbqt (*.pdb *.pdbqt)'
        )
        if nat_lig_path:
            self.nat_lig_path = nat_lig_path
            self.uiRedockingCheck.setEnabled(True)
            with open(nat_lig_path, 'r') as file:
                lines = file.readlines()
                coords = np.empty((0, 3), dtype=float)
                for line in lines:
                    if line.startswith('HETATM'):
                        x, y, z = line.split()[-6:-3]
                        if not x.count('.'):
                            x, y, z = line.split()[-5:-2]
                    elif line.startswith('ATOM'):
                        x, y, z = line.split()[-7:-4]
                        if not x.count('.'):
                            x, y, z = line.split()[-6:-3]
                    else:
                        continue
                    coords = np.append(
                        coords, np.array([[float(x), float(y), float(z)]]), axis=0
                    )
                center = tuple(np.average(coords, axis=0))
                self.uiCenterXDouble.setValue(center[0])
                self.uiCenterYDouble.setValue(center[1])
                self.uiCenterZDouble.setValue(center[2])
                self.uiNatLigLine.setText(nat_lig_path)

    def queue_calculation(self):
        if not self.config:
            receptor = self.uiReceptorLine.text()
            *receptor_name, receptor_format = receptor.split('/')[-1].split('.')
            if isinstance(receptor_name, list):
                receptor_name = ''.join(receptor_name)
            if self.project:
                receptor_file = f'projects/{self.project.name}/{receptor_name}_receptor.pdbqt'
            else:
                receptor_file = f'molecules/{self.molecule.inchi_key}/{receptor_name}_receptor.pdbqt'
            if receptor_format == 'pdbqt':
                src = receptor
                dst = receptor_file
                shutil.copyfile(src, dst)
            else:
                self.receptor_prep(receptor, receptor_file)
            # center = (
                # self.uiCenterXDouble.value(),
                # self.uiCenterYDouble.value(),
                # self.uiCenterZDouble.value()
            # )
            self.times = self.uiTimesSpin.value()
            self.conformer = self.uiConformerCombo.currentText()
            self.config = {
                'receptor': receptor_file,
                'center_x': self.uiCenterXDouble.value(),
                'center_y': self.uiCenterYDouble.value(),
                'center_z': self.uiCenterZDouble.value(),
                'size_x': self.uiSizeXDouble.value(),
                'size_y': self.uiSizeYDouble.value(),
                'size_z': self.uiSizeZDouble.value(),
                'energy_range': self.uiDiffEnergySpin.value(),
                'num_modes': self.uiNumModesSpin.value(),
                'exhaustiveness': self.uiExhausSpin.value(),
            }
        # center = (
            # self.config['center_x'],
            # self.config['center_y'],
            # self.config['center_z']
        # )
        *receptor_name, receptor_format = self.config['receptor'].split('/')[-1].split('.')
        receptor_name = receptor_name[-1]
        ligand = f'molecules/{self.molecule.inchi_key}/{self.molecule.inchi_key}_{receptor_name}.pdbqt'
        previous = len(glob(ligand.split('.')[0] + '_*.pdbqt'))
        output_file = [
            f'{ligand.split(".")[0]}_{i}.pdbqt' for i in range(previous + 1, previous + self.times + 1)
        ]
        self.config.update(
            {
                'ligand': ligand,
                'size_x': self.uiSizeXDouble.value(),
                'size_y': self.uiSizeYDouble.value(),
                'size_z': self.uiSizeZDouble.value(),
            }
        )
        keywords = ''
        for k, v in self.config.items():
            if k in ('ligand', 'receptor'):
                v = v.split('/')[-1]
            keywords += f'{k}: {v}\n'
        calculation = {
            'type': 'Docking',
            'config': self.config,
            'receptor_in': self.config['receptor'],
            'receptor_out': self.config['receptor'],
            'receptor_name': receptor_name,
            'times': self.times,
            'charge_mult': f'{self.molecule.get_formal_charge}/1',
            'keywords': keywords,
            'input_file': f'{self.config.get("ligand").split(".")[0]}_conf.txt',
            'output_file': output_file,
            'molecule': self.molecule.get_name,
            'molecule_id': self.molecule.inchi_key,
            'auto_box_size': self.auto_box_size,
            'nat_lig_path': self.nat_lig_path,
            'redocking': self.redocking,
            'conformer': self.conformer if self.conformer else 'input',
            'status': 'Pending'
        }
        self.submitted.emit(calculation, self.project.name)
        self.close()


class FPExplorer(qtw.QWidget):

    def __init__(self, project):
        super().__init__()
        uic.loadUi('Views/uiDisplayFP.ui', self)
        self.project = project
        fp_types = self.project.fps[0].keys()
        self.uiFPTypeCombo.addItems(fp_types)
        self.uiFPTypeCombo.setCurrentText('Morgan2')

    def change_fp_type(self, fp_type):
        if fp_type:
            self.df = self.create_df(fp_type)
            self.set_model()

    def create_df(self, fp_type):
        fps = pd.DataFrame([fp.get(fp_type) for fp in self.project.fps])
        fps.set_axis([str(i) for i in range(2048)], axis=1, inplace=True)
        fps.set_axis(
            [m.get_name for m in self.project.molecules], axis=0, inplace=True
        )
        fps.replace(0, np.nan, inplace=True)
        fps.dropna(how='all', axis=1, inplace=True)
        fps.replace(np.nan, 0, inplace=True)
        fps = fps.applymap(int)
        return fps

    def set_model(self):
        model = DisplayMFPModel(self.df)
        self.uiFPTable.setModel(model)
        self.uiFPTable.resizeColumnsToContents()
        self.uiFPTable.resizeRowsToContents()

class DescriptorsExplorer(qtw.QWidget):

    def __init__(self, project):
        super().__init__()
        uic.loadUi('Views/uiDisplayDescriptors.ui', self)
        self.project = project
        self.set_model()

    def set_model(self):
        self.df = self.project.descriptors
        model = PandasModel(self.df)
        self.uiDescriptorsTable.setModel(model)
        self.uiDescriptorsTable.resizeColumnsToContents()
        self.uiDescriptorsTable.resizeRowsToContents()


class SimilarityExplorer(qtw.QWidget):

    def __init__(self, project):
        super().__init__()
        uic.loadUi('Views/uiSimilarity.ui', self)
        self.names = [m.get_name for m in project.molecules]
        self.molecules = project.molecules
        self.similarity = Similarity(project.fps, self.names)
        self.simil_metric = 'Tanimoto'
        self.fp_type = 'Morgan2'
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)
        simil_metrics = (
            'Tanimoto', 'Dice', 'Cosine', 'Sokal',
            'Russel', 'Kulczynski', 'Hamann'
        )
        self.uiSimilMetricCombo.addItems(simil_metrics)
        self.uiSimilMetricCombo.setCurrentText('Tanimoto')
        fp_types = project.fps[0].keys()
        self.uiFPTypeCombo.addItems(fp_types)
        self.uiFPTypeCombo.setCurrentText('Morgan2')
        self.uiShowSimilMapButton.clicked.connect(
            lambda: self.show_simil_map(self.fp_type)
        )

    def set_model(self):
        self.simil_df = self.similarity.similarity_matrix(
            self.fp_type, self.simil_metric
        )
        self.HCL = self.get_HCL()
        simil_model = SimilarityModel(self.simil_df)
        self.uiSimilarityTable.setModel(simil_model)
        self.uiSimilarityTable.resizeColumnsToContents()
        self.uiSimilarityTable.resizeRowsToContents()

    def set_simil_metric(self, simil_metric):
        self.simil_metric = simil_metric
        self.set_model()

    def set_fp_type(self, fp_type):
        self.fp_type = fp_type
        if fp_type in (
            'Morgan1', 'Morgan2', 'Morgan3',
            'Hasehd Topological Torsions', 'RDKit',
            'Hashed Atom Pairs (chiral)',
            'Hashed Atom Pairs (achiral)'
        ):
            self.uiShowSimilMapButton.setEnabled(True)
        else:
            self.uiShowSimilMapButton.setEnabled(False)
        self.set_model()

    def show_heatmap(self):
        self.heatmap = Heatmap(self.simil_df, self.fp_type, self.simil_metric)
        self.heatmap.plot()
        self.heatmap.show()

    def get_HCL(self):
        self.canvas.ax.clear()
        simil_matrix = self.simil_df.to_numpy()
        linked = linkage(simil_matrix, 'single')
        dend = dendrogram(
            linked, orientation='left', labels=self.names,
            distance_sort='descending', show_leaf_counts=True,
            ax=self.canvas.ax
        )
        self.canvas.ax.set_title(
            f'Similarity HCL cluster ({self.fp_type}; {self.simil_metric})'
        )
        self.canvas.ax.spines['left'].set_visible(False)
        self.canvas.ax.spines['top'].set_visible(False)
        self.canvas.ax.spines['right'].set_visible(False)
        self.canvas.draw()

    def show_simil_map(self, fp_type):
        indexes = self.uiSimilarityTable.selectedIndexes()
        if not indexes:
            return
        i = indexes[0].row()
        j = indexes[0].column()
        ref = self.molecules[i].mol
        ref_name = self.molecules[i].get_name
        ref = Chem.MolFromSmiles(Chem.MolToSmiles(ref))
        prob = self.molecules[j].mol
        prob_name = self.molecules[j].get_name
        prob = Chem.MolFromSmiles(Chem.MolToSmiles(prob))
        d = Draw.MolDraw2DCairo(550, 550)
        if 'Morgan' in fp_type:
            radius = int(fp_type[-1])
            simil_function = lambda m, i: SimilarityMaps.GetMorganFingerprint(
                m, i, radius=radius, fpType='bv'
            )
        elif fp_type == 'Hasehd Topological Torsions':
            simil_function = lambda m, i: SimilarityMaps.GetTTFingerprint(
                m, i, fpType='bv'
            )
        elif 'Hashed Atom Pairs' in fp_type:
            simil_function = lambda m, i: SimilarityMaps.GetAPFingerprint(
                m, i, fpType='bv'
            )
        elif fp_type == 'RDKit':
            simil_function = lambda m, i: SimilarityMaps.GetRDKFingerprint(
                m, i, fpType='bv'
            )
        d.DrawMolecule(ref)
        _, maxWeight = SimilarityMaps.GetSimilarityMapForFingerprint(
            ref, prob, simil_function, draw2d=d
        )
        d.FinishDrawing()
        data = d.GetDrawingText()
        bio = io.BytesIO(data)
        img = Image.open(bio)
        setattr(
            self, f'simil_map_{fp_type}_{i}_{j}',
            SimilMap(img, ref_name, prob_name, fp_type)
        )


class SimilMap(qtw.QWidget):

    def __init__(self, img, ref_name, prob_name, fp_type):
        super().__init__()
        uic.loadUi('Views/uiSimilMap.ui', self)
        self.img = img
        self.uiSimilMapLabel.setPixmap(qtg.QPixmap.fromImage(
            ImageQt.ImageQt(img)
        ))
        self.uiLabel.setText(
            f'Ref: {ref_name}\tProb: {prob_name}\nFP type: {fp_type}')
        self.show()

    def save_image(self):
        filename, extension = qtw.QFileDialog.getSaveFileName(
            self, 'Save image', 'Image.png', 'PNG (*.png)',
            options=qtw.QFileDialog.DontUseNativeDialog,
        )
        if filename:
            has_extension = '.' in filename
            if not has_extension:
                filename += '.png'
            self.img.save(filename)


class Heatmap(qtw.QWidget):

    def __init__(self, data, fp_type, simil_metric):
        super().__init__()
        uic.loadUi('Views/uiHeatmap.ui', self)
        self.data = data
        self.fp_type = fp_type
        self.simil_metric = simil_metric
        self.canvas = PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiPlotLayout.addWidget(self.canvas)

    def plot(self):
        self.canvas.ax.clear()
        simil_matrix = self.data.to_numpy()
        labels = self.data.columns.values.tolist()
        self.canvas.ax.set_title(
            f'Similarity Matrix ({self.fp_type}; {self.simil_metric})'
        )
        im = self.canvas.ax.imshow(simil_matrix, cmap='magma_r')
        # Create colorbar
        cbar = self.canvas.ax.figure.colorbar(im, ax=self.canvas.ax)
        cbar.ax.set_ylabel('', rotation=-90, va="bottom")
        # Show all ticks and label them with the respective list entries.
        self.canvas.ax.set_xticks(
            np.arange(simil_matrix.shape[1]), labels=labels,
            rotation=90
        )
        self.canvas.ax.set_yticks(np.arange(simil_matrix.shape[0]), labels=labels)
        # Let the horizontal axes labeling appear on top.
        self.canvas.ax.tick_params(
            top=True, bottom=False,
            labeltop=True, labelbottom=False
        )
        self.canvas.ax.tick_params(axis='x', labelsize=8)
        self.canvas.ax.tick_params(axis='y', labelsize=8)
        # Turn spines off and create white grid.
        self.canvas.ax.spines[:].set_visible(False)
        self.canvas.ax.set_xticks(np.arange(simil_matrix.shape[1] + 1) - .5, minor=True)
        self.canvas.ax.set_yticks(np.arange(simil_matrix.shape[0] + 1) - .5, minor=True)
        self.canvas.ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
        self.canvas.ax.tick_params(which="minor", bottom=False, left=False)
        self.canvas.draw()

