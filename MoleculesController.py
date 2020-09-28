#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Computational_calculations.Gaussian import Gaussian
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from Models import MoleculesModel, PropertiesModel, StatusModel
from Worker import WorkerThread
from rdkit import Chem
from PIL import ImageQt  # , Image
from Views import resources
from DB.MoleculesDB import MyZODB
import multiprocessing
from psutil import virtual_memory


# mainWindowBase, mainWindowForm = uic.loadUiType('Views/uiMainWindow.ui')
# calcWindowBase, calcWindowForm = uic.loadUiType('Views/uiCalculationsWindow.ui')


class MainWindow(qtw.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
        # mainWindowBase.__init__(self)
        # self.setupUi(self)
    # --------------------------------------------------PARAMETERS--------------------------------------------------
        self._filePathList = []
        self._molecules = []
        self._calcWidgets = []
        self.database = MyZODB()
        self._propertiesWidgets = [
            self.uiName,
            self.uiNameLabel,
            self.uiFormula,
            self.uiFormulaLabel,
            self.uiMolarMass,
            self.uiMolarMassLabel,
            self.uiInchikey,
            self.uiInchikeyLabel,
            self.uiSmiles,
            self.uiSmilesLabel,
            self.uiDisplayAvailableCalcButton,
            self.uiSubmitCalcButton,
            self.ui2DImage
            ]
    # --------------------------------------------------STYLE--------------------------------------------------
        self.uiOpen.setIcon(qtg.QIcon(':/icons/openIcon.png'))
        self.uiSave.setIcon(qtg.QIcon(':/icons/saveIcon.png'))
        self.uiDeleteMoleculeButton.setIcon(qtg.QIcon(':/icons/trashIcon.png'))
        self.uiNew.setIcon(qtg.QIcon(':/icons/newIcon.png'))
        self.uiResetButton.setIcon(qtg.QIcon(':/icons/resetIcon.png'))
        self.uiSaveToDB.setIcon(qtg.QIcon(':/icons/addDB.png'))
        self.uiStartCalcButton.setIcon(qtg.QIcon(':/icons/startIcon.png'))
        self.uiStopCalcButton.setIcon(qtg.QIcon(':/icons/cancelIcon.png'))
        self.uiLoadMoleculeButton.setIcon(qtg.QIcon(':/icons/loadIcon.png'))
        self.uiRemoveButton.setIcon(qtg.QIcon(':/icons/trashIcon.png'))
        self.statusBar().showMessage('Welcome to SVP360 Calculation Manager')
        for w in self._propertiesWidgets:
            w.setVisible(False)
        self.uiSaveToDB.setEnabled(False)
        self.uiSave.setEnabled(False)
        self.uiDeleteMoleculeButton.setEnabled(False)
        self.uiResetButton.setEnabled(False)
        self.uiRemoveButton.setEnabled(False)

    # --------------------------------------------------SIGNALS--------------------------------------------------
        self.uiLoadMoleculeButton.pressed.connect(self.LoadMolecules)
        self.uiSave.triggered.connect(self.SaveFiles)
        self.uiMoleculesList.clicked.connect(self.ShowProperties)
        self.uiDeleteMoleculeButton.pressed.connect(self.DeleteMolecule)
        self.uiResetButton.pressed.connect(self.ResetListOfMolecules)
        self.uiSubmitCalcButton.clicked.connect(self.StackedWidgetsManager)
        self.uiStartCalcButton.clicked.connect(self.RunCalculations)
        self.uiSaveToDB.triggered.connect(self.StoreMolecule)
        self.uiStopCalcButton.clicked.connect(lambda: self.uiStartCalcButton.setEnabled(True))
        self.uiRemoveButton.clicked.connect(self.RemoveFromQueue)

#     --------------------------------------------------FUNCTIONS------------------------------------------------
    def LoadMolecules(self):
        """
        Open File Dialog to get the files of the molecules

        Returns
        -------
        Feeds the self._filePathList
        """

        filePathList, _ = qtw.QFileDialog.getOpenFileNames(
            self,
            'Select your molecule(s)',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Mol files(*.mol2)'
            )

        key_names_loaded = []
        if len([i for i in filePathList if i in self._filePathList]) == 0:
            for index, key in enumerate(filePathList, start=1):
                exec(f'self.mol{index} = Gaussian("{key}")')
                exec(f"key_names_loaded.append(str(self.mol{index}.GetInchikey))")
                exec(f'self._filePathList.append("{key}")')

            for index, key in enumerate(key_names_loaded, start=1):
                widgetName = key+'CalcWidget'
                if self.database.SearchMolecule(key):
                    exec(f'fetched_{key} = self.database.FetchMolecule("{key}")')
                    exec(f'self._molecules.append(fetched_{key})')
                    exec(f'{widgetName} = CalculationsController(fetched_{key}.heavyAtomsPresent)')
                else:
                    exec(f'{key} = self.mol{index}')
                    exec(f'self._molecules.append({key})')
                    exec(f'{widgetName} = CalculationsController({key}.heavyAtomsPresent)')
                exec(f'del self.mol{index}')
                exec(f'self._calcWidgets.append({widgetName})')
                self.uiSaveToDB.setEnabled(True)
                self.uiSave.setEnabled(True)
                self.uiDeleteMoleculeButton.setEnabled(True)
                self.uiResetButton.setEnabled(True)
        else:
            qtw.QMessageBox.warning(self, 'Oops', 'This molecule has been already uploaded!')

        numMoleculesAdded = len(key_names_loaded)
        self.listModel = MoleculesModel(self._molecules)
        self.uiMoleculesList.setModel(self.listModel)
        if numMoleculesAdded == 1:
            self.statusBar().showMessage('1 Molecule loaded!')
        else:
            self.statusBar().showMessage(f'{numMoleculesAdded} Molecules loaded!')

    def StoreMolecule(self):

        self.GetSelectedMoleculeAndWidgets()
        if not self._selectedMolecule.stored:
            try:
                self.database.StoreMolecule(self._selectedMolecule)
                self.GetSelectedMoleculeAndWidgets(
                    self.database.FetchMolecule(self._selectedMolecule.GetInchikey)
                    )
                self.listModel.layoutChanged.emit()
                self.statusBar().showMessage(
                    f'{self._selectedMolecule.GetName} succesfully added to DB!'
                    )
                print(self.database.ShowAll())  # quitar
                self.uiSubmitCalcButton.setChecked(False)
                self.listModel = MoleculesModel(self._molecules)  # ojo
                self.uiMoleculesList.setModel(self.listModel)  # ojo
            except:
                print('algo salió mal')
        else:
            pass  # self.statusBar().showMessage(f'{self._selectedMolecule.GetName} is already in DB!')

    def SaveFiles(self):
        pass

        # nameSDFile, _ = qtw.QFileDialog.getSaveFileName(
        #     self,
        #     'Save the SDF file of your molecules',
        #     options=qtw.QFileDialog.DontUseNativeDialog,
        #     filter='SDF files(*.sdf)'
        #     )
        # name = nameSDFile.split('.')[0]
        # sdMolFile = Chem.SDWriter(f'{name}.sdf')
        # for index, mol in enumerate(self._filePathList):
        #     name = self._molecules[index].GetName
        #     exec(f'rdkitMol{index} = Chem.MolFromMol2File("{mol}")')
        #     exec(f'rdkitMol{index}.SetProp("_Name", "{name}")')
        #     exec(f'sdMolFile.write(rdkitMol{index})')

    def GetSelectedMoleculeAndWidgets(self, fetchedMolecule=None):

        indexes = self.uiMoleculesList.selectedIndexes()

        if indexes:
            index = indexes[0]
            if fetchedMolecule is not None:
                self._molecules[index.row()] = fetchedMolecule
                self._calcWidgets[index.row()] = CalculationsController(fetchedMolecule)
            else:
                self._selectedMolecule = self._molecules[index.row()]
                self._selectedCalcWidget = self._calcWidgets[index.row()]

    def ShowProperties(self):
        """
        Set the model and the mapper to feed the widgets with the properties
        of the molecule
        """
        for w in self._propertiesWidgets:
            w.setVisible(True)
        self.GetSelectedMoleculeAndWidgets()
        self.uiSubmitCalcButton.setChecked(False)
        self.propertiesModel = PropertiesModel(self._selectedMolecule)
        self._dataMapper = qtw.QDataWidgetMapper()
        self._dataMapper.setModel(self.propertiesModel)
        self._dataMapper.setOrientation(qtc.Qt.Vertical)
        self._dataMapper.addMapping(self.uiName, 0)
        self._dataMapper.addMapping(self.uiFormula, 1)
        self._dataMapper.addMapping(self.uiMolarMass, 2)
        self._dataMapper.addMapping(self.uiInchikey, 3)
        self._dataMapper.addMapping(self.uiSmiles, 4)
        self.ui2DImage.setPixmap(
            qtg.QPixmap.fromImage(ImageQt.ImageQt(self._selectedMolecule.Get2DImage))
            )
        self._dataMapper.toFirst()
        self.StackedWidgetsManager()

    def DeleteMolecule(self):
        """
        Remove the selected molecule object from the main list of them,
        along with the path to it and the widgets to handle its calculations.
        """
        indexes = self.uiMoleculesList.selectedIndexes()
        if indexes:
            thisIndex = indexes[0]
            del self._molecules[thisIndex.row()]
            del self._filePathList[thisIndex.row()]
            del self._calcWidgets[thisIndex.row()]
            newIndex = thisIndex.row() - 1
            self.listModel.layoutChanged.emit()
            if newIndex >= 0:
                self.uiMoleculesList.setCurrentIndex(self.listModel.index(newIndex))
            else:
                self.uiSaveToDB.setEnabled(False)
                self.uiSave.setEnabled(False)
                self.uiDeleteMoleculeButton.setEnabled(False)
                self.uiResetButton.setEnabled(False)
            for w in self._propertiesWidgets:
                w.setVisible(False)
            self.uiStackedWidget.setVisible(False)
            self.statusBar().showMessage('1 Molecule deleted!')
            self.QueueManager()
            self.statusModel.layoutChanged.emit()

    def ResetListOfMolecules(self):
        """
        Remove all items in the list of molecules
        """
        totalMolecules = len(self._molecules)
        del self._molecules[:]
        del self._filePathList[:]
        del self._calcWidgets[:]
        self.listModel.layoutChanged.emit()
        for w in self._propertiesWidgets:
            w.setVisible(False)
        self.uiStackedWidget.setVisible(False)
        self.uiSaveToDB.setEnabled(False)
        self.uiSave.setEnabled(False)
        self.uiDeleteMoleculeButton.setEnabled(False)
        self.uiResetButton.setEnabled(False)
        self.statusBar().showMessage(f'{totalMolecules} Molecule(s) deleted!')
        self.QueueManager()
        self.statusModel.layoutChanged.emit()

    def RemoveFromQueue(self):

        indexes = self.uiCalculationFlowTable.selectedIndexes()
        if indexes:
            thisIndex = indexes[0]
            del self.masterQueue[thisIndex.row()]
            self.statusModel.layoutChanged.emit()
        if len(self.masterQueue) == 0:
            self.uiRemoveButton.setEnabled(False)
            # falta quitarlo del queue de cada widget

    def StackedWidgetsManager(self):

        self.GetSelectedMoleculeAndWidgets()
        if self.uiSubmitCalcButton.isChecked():
            self.uiCalculationsLayout.addWidget(self.uiStackedWidget)
            self.uiStackedWidget.setVisible(True)
            for w in self._calcWidgets:
                self.uiStackedWidget.addWidget(w)
            self.uiStackedWidget.setCurrentWidget(self._selectedCalcWidget)
        else:
            self.uiStackedWidget.setParent(None)

    def QueueManager(self):

        self.masterQueue = []
        for index, w in enumerate(self._calcWidgets):
            self.masterQueue.extend([[self._molecules[index], i] for i in w.queue])
        self.statusModel = StatusModel(self.masterQueue)
        self.uiCalculationFlowTable.setModel(self.statusModel)
        if len(self.masterQueue) == 0:
            self.uiStartCalcButton.setEnabled(False)
            self.uiRemoveButton.setEnabled(False)
        else:
            self.uiRemoveButton.setEnabled(True)
            self.uiStartCalcButton.setEnabled(True)

    def RunCalculations(self):

        self.worker = WorkerThread(self.masterQueue)
        self.thread = qtc.QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.Calculate)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(lambda: self.uiStartCalcButton.setEnabled(True))
        self.worker.finished.connect(lambda: self.uiStopCalcButton.setEnabled(False))
        self.worker.okCalc.connect(self.ExecutionFlowManager)
        self.worker.startingCalc.connect(self.ExecutionFlowManager)
        self.uiStopCalcButton.clicked.connect(self.worker.Stop)
        self.worker.error.connect(self.ErrorManager)
        self.thread.start()
        self.uiStartCalcButton.setEnabled(False)
        self.uiStopCalcButton.setEnabled(True)

    @qtc.pyqtSlot(object, int, str)
    def ErrorManager(self, errorMol, index, status):  # aquí hay que mejorar

        self.thread.quit()
        self.statusBar().showMessage(
            f'Something get wrong with {errorMol.GetName} with index: {index}. Please check the log file'
            )
        self.masterQueue[index][1][6] = 'Failed'
        self.statusModel.layoutChanged.emit()
        if index != (len(self.masterQueue) - 1):
            self.RunCalculations()

    @qtc.pyqtSlot(object, int, str)
    def ExecutionFlowManager(self, procMol, index, status):
        if status == 'Running':
            self.masterQueue[index][1][6] = 'Running'
            self.statusModel.layoutChanged.emit()
            self.statusBar().showMessage(f'Currently processing {procMol.GetName}')
        elif status == 'Finished':
            self.masterQueue[index][1][6] = 'Finished'
            self.statusModel.layoutChanged.emit()
            self.statusBar().showMessage(f'Already processed {procMol.GetName}')

    def closeEvent(self, event):
        reply = qtw.QMessageBox.question(
            self, 'Window Close', 'Are you sure you want to close the window?',
            qtw.QMessageBox.Yes | qtw.QMessageBox.No, qtw.QMessageBox.No
            )

        if reply == qtw.QMessageBox.Yes:
            self.database.Close()
            event.accept()
            print('Window closed')
        else:
            event.ignore()


class CalculationsController(qtw.QWidget):

    def __init__(self, heavy):
        super(CalculationsController, self).__init__()
        uic.loadUi('Views/uiCalculationsWindow.ui', self)
        # self.setupUi(self)
        mem = virtual_memory()
        self.memory = mem.total//1000000 - 2000
        self.cpu = multiprocessing.cpu_count() - 1
        self._heavy = heavy
        self.queue = []
        self._methods = (None, 'DFT', 'Semiempirical')
        self._semiempiricalFunctionals = (None, 'AM1', 'PM3', 'PM6')
        self._dftFunctionals = (None, 'B3LYP', 'wB97XD', 'MPWB1K', 'M06-2X')
        self._basisSetsLightAtoms = (
            None, '6-31G(d)', '6-31G(d,p)', '6-31+G(d)', '6-31+G(d,p)',
            '6-311G(d)', '6-311G(d,p)', '6-311+G(d)', '6-311+G(d,p)'
            )
        self._basisSetsHeavyAtoms = (None, 'LANL2DZ', 'SDD')
        self.chosenParameters = {
            'METHOD': None, 'FUNCTIONAL': None, 'BASIS': None,
            'BASIS2': None, 'SAVE FILE': False, 'STATUS': 'Pending'
            }
        self._widgets = [
            self.uiGeoptMethodComboBox,
            self.uiGeoptDftFunctionalComboBox,
            self.uiGeoptSemiFunctionalComboBox,
            self.uiGeoptBasisComboBox,
            self.uiGeoptBasis2ComboBox,
            self.uiSaveOptFileCheckBox,
            self.uiGeoptQueueButton
            ]
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


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
