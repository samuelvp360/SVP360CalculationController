#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import uic
from Models import ResultsModel
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema


class ResultsWidget(qtw.QWidget):
    """
    docstring
    """
    def __init__(self, calculation):
        super(ResultsWidget, self).__init__()
        uic.loadUi('Views/uiResultsWidget.ui', self)
        model = ResultsModel(calculation)
        self.uiResultsTableView.setModel(model)

    def IRPlotter(self, csvData1, threshold=0.9975, orderMin=25, orderMax=5):
        """
        docstring
        """
        label = csvData1.replace('.csv', '')
        df = pd.read_csv(csvData1, names=('Wavelenght', 'Transmitance'))
        maxValue = df['Transmitance'].max()
        zeroCorrection = 1 - maxValue
        df['Transmitance'] = df['Transmitance'] + zeroCorrection

        df['min'] = df.loc[argrelextrema(df['Transmitance'], np.less_equal, order=orderMin)[0], 'Transmitance']
        df['max'] = df.loc[argrelextrema(df['Transmitance'], np.less_equal, order=orderMax)[0], 'Transmitance']

        mask1 = (df['max'] < threshold)
        df.loc[mask1, 'max'] = np.nan
        mask2 = (df['max'].isnull() is False)
        baselineX = df.loc[mask2, 'Wavelenght']
        baselineY = df.loc[mask2, 'Wavelenght']

        plt.scatter(df['Wavelenght'], df['min'], c='r')
        coef = np.polyfit(baselineX, baselineY, 4)
        baselineFit = np.poly1d(coef)
        plt.plot(df['Wavelenght'], baselineFit(df['Wavelenght']), c='k')
        plt.scatter(baselineX, baselineY, c='k')
        plt.plot(df['Wavelenght'], df['Transmitance'], label=label)
        plt.xlim(4000, 1000)
        plt.xlabel('Wavelenght')
        plt.ylabel('Transmitance')
        plt.title('IR')
        plt.legend()
        plt.show()
