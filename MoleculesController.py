#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import transaction
from Molecules import Molecules
from CalculationsController import CalculationsController
from ResultsController import ResultsWidget
from IRController import SpectrumSelector
from Models import MoleculesModel, StatusModel, AvailableCalcModel
from Worker import WorkerThread
from Views import resources
from DB.MoleculesDB import MyZODB

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import uic

from rdkit import Chem
from PIL import ImageQt  # , Image


class MainWindow(qtw.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        uic.loadUi('Views/uiMainWindow.ui', self)
    # --------------------------------------------------PARAMETERS--------------------------------------------------
        self._filePathList = []
        self._molecules = []
        self.masterQueue = []
        self.database = MyZODB()
        self._propertiesWidgets = [
            self.uiName,
            self.uiDisplayAvailableCalcButton,
            self.uiSubmitCalcButton,
            self.ui2DImage,
            self.uiChangeNameButton
        ]
    # --------------------------------------------------STYLE-------------------------------------------------------
        # self.fontDB = qtg.QFontDatabase()
        # self.fontDB.addApplicationFont(":/fonts/CascadiaCode.ttf")
        # self.setFont(qtg.QFont("CascadiaCode", 11))
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
        self.uiSubmitCalcButton.setIcon(qtg.QIcon(':/icons/calculateIcon.png'))
        self.uiDisplayAvailableCalcButton.setIcon(qtg.QIcon(':/icons/displayCalcIcon.png'))
        self.statusBar().showMessage('Welcome to SVP360 Calculation Manager')
        self.uiSaveToDB.setEnabled(False)
        self.uiSave.setEnabled(False)
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
        self.uiSubmitCalcButton.clicked.connect(self.SubmitCalc)
        self._propertiesWidgets[1].clicked.connect(self.DisplayAvailableCalcs)
        self.uiShowResultsButton.clicked.connect(self.ShowResults)
        self.uiSpectraManagerButton.clicked.connect(self.DisplaySpectrum)
        self.uiUpdateDB.triggered.connect(self.StoreChanges)
        self._propertiesWidgets[4].clicked.connect(lambda: self._selectedMolecule.SetName(self._propertiesWidgets[0].text()))
    # ---------------------------------------------------METHODS--------------------------------------------------

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
                exec(f'self.mol{index} = Molecules("{key}")')
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

        self.GetSelectedMolecule()
        if not self._selectedMolecule.stored:
            try:
                self.database.StoreMolecule(self._selectedMolecule)
                self.GetSelectedMolecule(
                    self.database.FetchMolecule(self._selectedMolecule.GetInchikey)
                )
                self.listModel.layoutChanged.emit()
                self.statusBar().showMessage(
                    f'{self._selectedMolecule.GetName} succesfully added to DB!'
                )
                self.uiSubmitCalcButton.setChecked(False)
                self.listModel = MoleculesModel(self._molecules)
                self.uiMoleculesList.setModel(self.listModel)
            except:
                self.statusBar().showMessage('Something was wrong')
        else:
            self.statusBar().showMessage(f'{self._selectedMolecule.GetName} is already in DB!')

    def StoreChanges(self):

        self.GetSelectedMolecule()
        self._selectedMolecule.UpdateChanges()

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
        self.uiSpectraManagerButton.setEnabled(True)

    def ShowProperties(self):
        """
        Show the properties of the selected molecule in the status bar
        """
        self.GetSelectedMolecule()
        name = self._selectedMolecule.GetName
        weight = self._selectedMolecule.GetMolarMass
        formula = self._selectedMolecule.GetForm
        inchikey = self._selectedMolecule.GetInchikey
        smiles = self._selectedMolecule.GetSmiles
        [w.setEnabled(True) for w in self._propertiesWidgets if w != 1]
        if len(self._selectedMolecule.GetCalculations) > 0:
            self._propertiesWidgets[1].setEnabled(True)
        else:
            self._propertiesWidgets[1].setEnabled(False)
        self._propertiesWidgets[0].setText(name)
        self.statusBar().showMessage(f'{name}  {formula}  FW:{weight} uma  InchiKey:{inchikey}  Smiles:{smiles}')
        self.ui2DImage.setPixmap(
            qtg.QPixmap.fromImage(ImageQt.ImageQt(self._selectedMolecule.Get2DImage))
        )

    def SubmitCalc(self):
        """
        docstring
        """
        self.calculationSetupController = CalculationsController(self._selectedMolecule)
        self.calculationSetupController.submitted.connect(self.QueueManager)
        self.calculationSetupController.show()

    def DeleteMolecule(self):
        """
        Remove the selected molecule object from the main list of them,
        along with the path through.
        """
        indexes = self.uiMoleculesList.selectedIndexes()
        if indexes:
            thisIndex = indexes[0]
            newMasteQueue = []
            for idx, i in enumerate(self.masterQueue):
                if i[0].GetName == self._molecules[thisIndex.row()].GetName:
                    os.remove(self.masterQueue[idx][5])
                    del i
                else:
                    newMasteQueue.append(i)
            self.masterQueue = newMasteQueue
            del self._molecules[thisIndex.row()]
            del self._filePathList[thisIndex.row()]
            newIndex = thisIndex.row() - 1
            self.statusModel = StatusModel(self.masterQueue)
            self.uiCalculationFlowTable.setModel(self.statusModel)
            self.listModel.layoutChanged.emit()
            if len(self._molecules) > 0:
                self.uiMoleculesList.setCurrentIndex(self.listModel.index(newIndex))
                self.ShowProperties()
            else:
                self.ResetListOfMolecules()
                [w.setEnabled(False) for w in self._propertiesWidgets]

            self.statusBar().showMessage('1 Molecule deleted!')

    def ResetListOfMolecules(self):
        """
        Remove all items in the list of molecules
        """
        totalMolecules = len(self._molecules)
        del self._molecules[:]
        del self._filePathList[:]
        self.listModel.layoutChanged.emit()
        self.uiName.setText('')
        self.uiDisplayAvailableCalcButton.setEnabled(False)
        self.uiSubmitCalcButton.setEnabled(False)
        self.uiSpectraManagerButton.setEnabled(False)
        self.ui2DImage.clear()
        self.uiSaveToDB.setEnabled(False)
        self.uiSave.setEnabled(False)
        self.uiDeleteMoleculeButton.setEnabled(False)
        self.uiResetButton.setEnabled(False)
        self.statusBar().showMessage(f'{totalMolecules} Molecule(s) deleted!')
        self.masterQueue = []
        self.uiCalculationFlowTable.setModel(None)
        try:
            self.statusModel.layoutChanged.emit()
        except AttributeError:
            pass

    def RemoveFromQueue(self):

        indexes = self.uiCalculationFlowTable.selectedIndexes()
        if indexes:
            thisIndex = indexes[0]
            os.remove(self.masterQueue[thisIndex.row()][5])
            del self.masterQueue[thisIndex.row()]
            self.statusModel.layoutChanged.emit()
        if len(self.masterQueue) == 0:
            self.uiRemoveButton.setEnabled(False)
            self.uiStopCalcButton.setEnabled(False)

    @qtc.pyqtSlot(list)
    def QueueManager(self, calculation):

        self.masterQueue.append(calculation)
        self.statusModel = StatusModel(self.masterQueue)
        self.uiCalculationFlowTable.setModel(self.statusModel)
        self.uiCalculationFlowTable.resizeColumnsToContents()
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
        self._propertiesWidgets[1].setChecked(False)
        self.uiAvailableCalcsTableView.setModel(None)

    @qtc.pyqtSlot(object, int, str)
    def ErrorManager(self, errorMol, index, status):  # aquí hay que mejorar

        self.thread.quit()
        self.statusBar().showMessage(
            f'Something get wrong with {errorMol.GetName} with index: {index}. Please check the log file'
        )
        self.masterQueue[index][3] = 'Failed'
        self.statusModel.layoutChanged.emit()
        if index != (len(self.masterQueue) - 1):
            self.RunCalculations()

    @qtc.pyqtSlot(object, int, str)
    def ExecutionFlowManager(self, procMol, index, status):
        if status == 'Running':
            self.masterQueue[index][3] = 'Running'
            self.statusModel.layoutChanged.emit()
            self.statusBar().showMessage(f'Currently processing {procMol.GetName}')
        elif status == 'Finished':
            self.masterQueue[index][3] = 'Finished'
            self.statusModel.layoutChanged.emit()
            self.statusBar().showMessage(f'Already processed {procMol.GetName}')

    def DisplayAvailableCalcs(self):
        """
        docstring
        """
        if self._propertiesWidgets[1].isChecked():
            availableCalcModel = AvailableCalcModel(self._selectedMolecule)
            self.uiAvailableCalcsTableView.setModel(availableCalcModel)
            self.uiAvailableCalcsTableView.resizeColumnsToContents()
            self.uiShowResultsButton.setEnabled(True)
        else:
            self.uiAvailableCalcsTableView.setModel(None)
            self.uiShowResultsButton.setEnabled(False)

    def ShowResults(self):
        """
        docstring
        """
        indexes = self.uiAvailableCalcsTableView.selectedIndexes()
        if indexes:
            index = indexes[0]
            self.resultsWidget = ResultsWidget(self._selectedMolecule.GetCalculations[index.row()]['RESULTS'])
            self.resultsWidget.show()

    def DisplaySpectrum(self):
        """
        docstring
        """
        self.SpectrumSelector = SpectrumSelector(self._selectedMolecule)
        self.SpectrumSelector.show()

    def CloseAll(self):
        self.database.Close()
        try:
            self.SpectrumSelector.close()
        except AttributeError:
            pass

    def closeEvent(self, event):

        try:
            if self._selectedMolecule._p_changed:
                reply = qtw.QMessageBox.question(
                    self, 'Window Close',
                    'Some changes have not been stored yet, do you want to commit this changes?',
                    qtw.QMessageBox.Yes | qtw.QMessageBox.No | qtw.QMessageBox.Cancel,
                    qtw.QMessageBox.No
                )

                if reply == qtw.QMessageBox.Yes:
                    transaction.commit()
                    self._selectedMolecule._p_changed = False
                    transaction.commit()
                    self.CloseAll()
                    event.accept()
                    print('Window closed')
                elif reply == qtw.QMessageBox.No:
                    transaction.abort()
                    self.CloseAll()
                    event.accept()
                    print('Window closed')
                else:
                    event.ignore()
            else:
                self.CloseAll()
                event.accept()
        except AttributeError:
            self.CloseAll()
            event.accept()


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
