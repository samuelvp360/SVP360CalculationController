#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PyQt5 import QtWidgets as qtw
from PyQt5 import uic
import pandas as pd
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from scipy.signal import argrelextrema

matplotlib.use('Qt5Agg')


class IRCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(11, 4), dpi=100, tight_layout=True, facecolor='#f8b4a5')
        self.axes = self.fig.add_subplot(111)
        super(IRCanvas, self).__init__(self.fig)


class IRPlotter(qtw.QWidget):
    """
    docstring
    """
    def __init__(self, csvData1):
        super(IRPlotter, self).__init__()
        uic.loadUi('Views/uiIRPlotter.ui', self)
        self.df = pd.read_csv(csvData1, names=('Wavelenght', 'Transmitance'))
        self.title = csvData1.replace('.csv', '')
        self.threshold = 0.9975
        self.canvas = IRCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiSpectraLayout.addWidget(self.canvas)
        self.Plot()
# ------------------------------------SIGNALS---------------------------------------------------------------
        self.uiZeroButton.clicked.connect(self.ZeroCorrection)
        self.uiDetectBaselineButton.clicked.connect(self.BaselineDetection)
        self.uiSubsBaselineButton.clicked.connect(self.BaselineSubs)
        self.uiThresholdSpinBox.valueChanged.connect(self.SetThreshold)
# ------------------------------------METHODS---------------------------------------------------------------

    def Plot(self, yList='Transmitance', showGrid=True):
        """
        docstring
        """
        self.canvas.axes.clear()
        self.canvas.axes.plot(self.df['Wavelenght'], self.df[yList])
        self.canvas.axes.set_title(self.title)
        self.canvas.axes.set_xlabel('Wavelenght [cm-1]')
        self.canvas.axes.set_ylabel('Transmitance')
        self.canvas.axes.set_xlim(4000, 1000)
        self.canvas.axes.grid(True)
        self.canvas.draw()

    def SetThreshold(self):
        """
        docstring
        """
        self.threshold = float(self.uiThresholdSpinBox.value())

    def ZeroCorrection(self):
        """
        docstring
        """
        maxValue = self.df['Transmitance'].max()
        zeroCorrection = 1 - maxValue
        self.df['Transmitance'] = self.df['Transmitance'] + zeroCorrection
        self.Plot()

    def BaselineDetection(self):
        """
        docstring
        """
        orderBase = 50
        self.df['baseline'] = self.df.loc[argrelextrema(self.df.Transmitance.values, np.greater_equal, order=orderBase)[0], 'Transmitance']
        mask1 = (self.df['baseline'] < self.threshold)
        self.df.loc[mask1, 'baseline'] = np.nan

        mask2 = (self.df['baseline'].isnull() == False)
        baselineX = self.df.loc[mask2, 'Wavelenght']
        baselineY = self.df.loc[mask2, 'Transmitance']
        coef = np.polyfit(baselineX, baselineY, 1)
        self.baselineFit = np.poly1d(coef)
        self.df['baseline'] = pd.Series([self.baselineFit(i) for i in self.df['Wavelenght']])
        self.Plot(yList=['Transmitance', 'baseline'])

    def BaselineSubs(self):
        """
        docstring
        """
        self.df['Transmitance'] = 1 - self.df['baseline'] + self.df['Transmitance']
        self.Plot()

    def PeakDetection(self):
        """
        docstring
        """
        orderPeaks = 20
        self.df['peaks'] = self.df.loc[argrelextrema(self.df.Transmitance.values, np.less_equal, order=orderPeaks)[0], 'Transmitance']
        mask3 = (self.df['peaks'] > self.threshold)
        self.df.loc[mask3, 'peaks'] = np.nan
        self.Plot()
