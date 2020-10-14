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
from CalculationsController import CalculationsController


class MainWindow(qtw.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
    # --------------------------------------------------PARAMETERS--------------------------------------------------
        self._filePathList = []
        self._molecules = []
        self.database = MyZODB()
        self.calculationSetupController = CalculationsController()
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
        self.uiCalculationsLayout.addWidget(self.calculationSetupController)
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
        self.calculationSetupController.setVisible(False)

    # --------------------------------------------------SIGNALS--------------------------------------------------
        self.uiLoadMoleculeButton.pressed.connect(self.LoadMolecules)
        self.uiSave.triggered.connect(self.SaveFiles)
        self.uiMoleculesList.clicked.connect(self.ShowProperties)
        self.uiDeleteMoleculeButton.pressed.connect(self.DeleteMolecule)
        self.uiResetButton.pressed.connect(self.ResetListOfMolecules)
        self.uiStartCalcButton.clicked.connect(self.RunCalculations)
        self.uiSaveToDB.triggered.connect(self.StoreMolecule)
        self.uiStopCalcButton.clicked.connect(lambda: self.uiStartCalcButton.setEnabled(True))
        self.uiRemoveButton.clicked.connect(self.RemoveFromQueue)
        self.uiSubmitCalcButton.pressed.connect(lambda: self.calculationSetupController.setVisible(True))
        self.uiSubmitCalcButton.pressed.connect(lambda: self.calculationSetupController.DetectMolecule(self._selectedMolecule))
        self.uiMoleculesList.clicked.connect(lambda: self.calculationSetupController.DetectMolecule(self._selectedMolecule))

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
                if self.database.SearchMolecule(key):
                    exec(f'fetched_{key} = self.database.FetchMolecule("{key}")')
                    exec(f'self._molecules.append(fetched_{key})')
                else:
                    exec(f'{key} = self.mol{index}')
                    exec(f'self._molecules.append({key})')
                exec(f'del self.mol{index}')
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

        nameSDFile, _ = qtw.QFileDialog.getSaveFileName(
            self,
            'Save the SDF file of your molecules',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='SDF files(*.sdf)'
        )
        name = nameSDFile.split('.')[0]
        sdMolFile = Chem.SDWriter(f'{name}.sdf')
        for index, mol in enumerate(self._filePathList):
            name = self._molecules[index].GetName
            exec(f'rdkitMol{index} = Chem.MolFromMol2File("{mol}")')
            exec(f'rdkitMol{index}.SetProp("_Name", "{name}")')
            exec(f'sdMolFile.write(rdkitMol{index})')

    def GetSelectedMolecule(self, fetchedMolecule=None):

        indexes = self.uiMoleculesList.selectedIndexes()

        if indexes:
            index = indexes[0]
            if fetchedMolecule is not None:
                self._molecules[index.row()] = fetchedMolecule
            else:
                self._selectedMolecule = self._molecules[index.row()]

    def ShowProperties(self):
        """
        Set the model and the mapper to feed the widgets with the properties
        of the molecule
        """
        for w in self._propertiesWidgets:
            w.setVisible(True)
        self.GetSelectedMolecule()
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
        self.listModel.layoutChanged.emit()
        for w in self._propertiesWidgets:
            w.setVisible(False)
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


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
