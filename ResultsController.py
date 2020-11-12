#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import uic
from Models import ResultsModel


class ResultsWidget(qtw.QWidget):
    """
    docstring
    """
    def __init__(self, calculation):
        super(ResultsWidget, self).__init__()
        uic.loadUi('Views/uiResultsWidget.ui', self)
        model = ResultsModel(calculation)
        self.uiResultsTableView.setModel(model)
