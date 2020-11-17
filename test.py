#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from scipy.signal import argrelextrema

matplotlib.use('Qt5Agg')


class IRCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, width=5, height=4, dpi=200):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(IRCanvas, self).__init__(fig)
        

        # 1 - baselineFit(self.df['Wavelenght']) + self.df['Transmitance'] if baselineFit(self.df['Wavelenght']) > self.df['Transmitance'] else self.df['Transmitance']
        plt.scatter(self.df['Wavelenght'], self.df['min'], c='r')
        plt.scatter(baselineX, baselineY, c='k')
        plt.plot(self.df['Wavelenght'], self.df['Transmitance'], label=label)
        plt.xlim(4000, 1000)
        plt.xlabel('Wavelenght')
        plt.ylabel('Transmitance')
        plt.title('IR')
        plt.legend()
        plt.show()


class IRPlotter(qtw.QWidget):
    """
    docstring
    """
    def __init__(self, csvData1, threshold=0.9975, orderMin=20, orderMax=50):
        """
        docstring
        """
        label = csvData1.replace('.csv', '')
        self.df = pd.read_csv(csvData1, names=('Wavelenght', 'Transmitance'))
        # zero correction
        maxValue = self.df['Transmitance'].max()
        zeroCorrection = 1 - maxValue
        self.df['Transmitance'] = self.df['Transmitance'] + zeroCorrection
        # baseline detection
        self.df['max'] = self.df.loc[argrelextrema(self.df.Transmitance.values, np.greater_equal, order=orderMax)[0], 'Transmitance']
        mask1 = (self.df['max'] < threshold)
        self.df.loc[mask1, 'max'] = np.nan

        mask2 = (self.df['max'].isnull() == False)
        baselineX = self.df.loc[mask2, 'Wavelenght']
        baselineY = self.df.loc[mask2, 'Transmitance']

        coef = np.polyfit(baselineX, baselineY, 1)
        self.baselineFit = np.poly1d(coef)
        self.df['Transmitance'] = pd.Series([1 - self.baselineFit(i) + i for i in self.df['Transmitance']])

        self.df['min'] = self.df.loc[argrelextrema(self.df.Transmitance.values, np.less_equal, order=orderMin)[0], 'Transmitance']
        self.df['max'] = self.df.loc[argrelextrema(self.df.Transmitance.values, np.greater_equal, order=orderMax)[0], 'Transmitance']

        mask3 = (self.df['min'] > threshold)
        self.df.loc[mask3, 'min'] = np.nan

        canvas = IRCanvas(self, width=5, height=4, dpi=100)
        self.df.plot(ax=canvas.axes)
        self.setCentralWidget(canvas)
        self.show()
