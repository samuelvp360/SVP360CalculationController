#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.ticker as ticker
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import uic
from io import StringIO
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from scipy.signal import argrelextrema, savgol_filter
from Models import IRDataModel
from Worker import IRWorkerThread
from Views import resources

matplotlib.use('Qt5Agg')


class SnappingCursor:
    """
    A cross hair cursor that snaps to the data point of a line, which is
    closest to the *x* position of the cursor.

    For simplicity, this assumes that *x* values of the data are sorted.
    """
    def __init__(self, ax, line):
        self.ax = ax
        self.horizontal_line = ax.axhline(color='k', lw=0.8, ls='--')
        self.vertical_line = ax.axvline(color='k', lw=0.8, ls='--')
        self.x, self.y = line.get_data()
        self._last_index = None
        # text location in axes coords
        self.text = ax.text(0.72, 0.9, '', transform=ax.transAxes)

    def set_cross_hair_visible(self, visible):
        need_redraw = self.horizontal_line.get_visible() != visible
        self.horizontal_line.set_visible(visible)
        self.vertical_line.set_visible(visible)
        self.text.set_visible(visible)
        return need_redraw

    def on_mouse_move(self, event):
        if not event.inaxes:
            self._last_index = None
            need_redraw = self.set_cross_hair_visible(False)
            if need_redraw:
                self.ax.figure.canvas.draw()
        else:
            self.set_cross_hair_visible(True)
            x, y = event.xdata, event.ydata
            index = len(self.x) - np.searchsorted(np.flip(self.x), x)
            # index = min(np.searchsorted(self.x, x), len(self.x) - 1)
            if index == self._last_index:
                return  # still on the same data point. Nothing to do.
            self._last_index = index
            x = self.x[index]
            y = self.y[index]
            # update the line positions
            self.horizontal_line.set_ydata(y)
            self.vertical_line.set_xdata(x)
            self.ax.figure.canvas.draw()


class IRCanvas(FigureCanvasQTAgg):
    """
    docstring
    """
    def __init__(self):
        self.fig = Figure(figsize=(6, 4), dpi=100)  #, facecolor='#2d2a2e')
        self.grid = GridSpec(17, 1, left=0.1, bottom=0.15, right=0.94, top=0.94, wspace=0.3, hspace=0.3)
        self.ax = self.fig.add_subplot(self.grid[0:, 0])
        self.ax2 = self.ax.twinx()
        super().__init__(self.fig)


class IRPlotter(qtw.QMainWindow):
    """
    docstring
    """
    data_sent = qtc.pyqtSignal(dict)

    def __init__(
        self,
        title='',
        solvent='',
        df=None,
        bands_df=None,
        predicted_bands=None
    ):
        super().__init__()
        uic.loadUi('Views/uiIRPlotter_new.ui', self)
        self.uiTransButton.setIcon(qtg.QIcon(':/icons/transIcon.png'))
        self.uiAbsorbanceButton.setIcon(qtg.QIcon(':/icons/absIcon.png'))
        self.uiPercentTransButton.setIcon(qtg.QIcon(':/icons/percentIcon.png'))
        self.uiDetectBaselineButton.setIcon(qtg.QIcon(':/icons/baselineIcon.png'))
        self.uiNormalizeButton.setIcon(qtg.QIcon(':/icons/normalizeIcon.png'))
        self.uiShowGridButton.setIcon(qtg.QIcon(':/icons/gridIcon.png'))
        self.canvas = IRCanvas()
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.uiToolbarLayout.addWidget(self.toolbar)
        self.uiSpectraLayout.addWidget(self.canvas)
        self.cid = self.canvas.fig.canvas.mpl_connect(
            'button_press_event', self.bands_detect_manual
        )
        self.title = title
        self.solvent = solvent
        self.predicted_bands = predicted_bands
        self.bands_fitted = False if not predicted_bands else True
        if bands_df is not None:
            self.uiBandFittingAutoButton.setEnabled(True)
            self.uiRemoveBandButton.setEnabled(True)

        self.base_threshold = self.uiBaseThresholdSpinBox.value()
        self.band_threshold = self.uiBandsThresholdSpinBox.value()
        self.order_base = self.uiOrderBaseSpinBox.value()
        self.order_bands = self.uiOrderBandsSpinBox.value()
        self.show_grid = True
        self.normalized = False
        self.snap_cursor = None
        self.baseline = False

        self.df = df
        self.copied_cells = None
        self.cwd = ''
        self.y_axis = '%Transmittance'
        self.bands_df = pd.DataFrame({
            'Wavenumber': [],
            'Transmittance': [],
            '%Transmittance': [],
            'Absorbance': [],
            'Transmittance(N)': [],
            '%Transmittance(N)': [],
            'Absorbance(N)': [],
            'HWHM': []
        }) if not bands_df else bands_df
# ------------------------------------SIGNALS---------------------------------------------------------------
        self.uiTransButton.clicked.connect(lambda: self.select_y_axis('Transmittance'))
        self.uiPercentTransButton.clicked.connect(lambda: self.select_y_axis('%Transmittance'))
        self.uiAbsorbanceButton.clicked.connect(lambda: self.select_y_axis('Absorbance'))

    def set_title(self, value):
        self.title = value
        self.statusBar().showMessage(
            f'Showing IR spectrum from {self.title} ({self.solvent})'
        )
        self.plot()

    def set_solvent(self, value):
        self.solvent = value
        self.statusBar().showMessage(
            f'Showing IR spectrum from {self.title} ({self.solvent})'
        )
        self.plot()

    def open_file(self):
        filename, _ = qtw.QFileDialog.getOpenFileName(
            self, 'Selecciona la molécula a cargar',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='CSV files(*.csv)'
        )
        if filename:
            *path, name = filename.split('/')
            path = '/'.join(path)
            name = ''.join(name.split('.')[:-1])
            self.cwd = (path)
            self.title = name
            self.bands_df = None
            with open(filename, 'r') as file:
                lines = file.readlines()
                new_lines = []
                for line in lines:
                    if line[0].isdigit():
                        new_lines.append(line)

                csv_str = ''.join(new_lines)

            try:
                self.df = pd.read_csv(
                    StringIO(csv_str), dtype=float,
                    names=('Wavenumber', 'Raw Data')
                )
            except ValueError:
                self.df = pd.read_csv(
                    StringIO(csv_str), delimiter=';',
                    decimal=',', dtype=float,
                    names=('Wavenumber', 'Raw Data')
                )
            self.set_y_axis()
            self.select_y_axis('%Transmittance')

    def plot(self):
        """
        docstring
        """
        self.canvas.ax.clear()
        self.canvas.ax2.clear()
        self.line = self.canvas.ax.plot(self.df['Wavenumber'], self.df[self.y_axis], linewidth=1.5, color='#272822')
        self.canvas.ax.set_title(f'{self.title} ({self.solvent})')
        self.canvas.ax.set_xlabel('Wavenumber [cm-1]')
        self.canvas.ax.set_ylabel(self.y_axis)  #, color='#f92672')
        self.canvas.ax.tick_params(axis='x')  #, colors='#66d9ef')
        self.canvas.ax.tick_params(axis='y')  #, colors='#66d9ef')
        self.canvas.ax2.tick_params(axis='y')  #, colors='#66d9ef')
        self.canvas.ax.set_xlim(4000, 500)
        self.canvas.ax.grid(
            True, linestyle=':', linewidth=0.5, color='#2d2a2e'
        ) if self.show_grid else self.canvas.ax.grid(False)
        if self.uiBandManualButton.isChecked():
            self.snap_cursor = SnappingCursor(self.canvas.ax, self.line[0])
            self.canvas.mpl_connect('motion_notify_event', self.snap_cursor.on_mouse_move)
        else:
            self.snap_cursor = None
        if self.baseline:
            self.canvas.ax.plot(
                self.df['Wavenumber'], self.baseline_df[self.y_axis],
                linewidth=1.0, color='#f92672'
            )
        if self.bands_fitted:
            self.canvas.ax2.set_ylabel('Residuals')  #, color='#f92672')
            for _, predicted in self.predicted_bands[self.y_axis].iteritems():
                self.canvas.ax.plot(self.df['Wavenumber'], predicted, linewidth=0.0)
                axis = self.y_axis.split('(')[0]
                if axis == 'Transmittance' or axis == '%Transmittance':
                    self.canvas.ax2.set_ylim(0.1, -0.005)
                else:
                    self.canvas.ax2.set_ylim(-0.005, 0.1)
                self.canvas.ax2.plot(
                    self.df['Wavenumber'], self.predicted_bands['Residuals'],
                    '-r', markersize=0.2, linewidth=0.1
                )
                self.canvas.ax.plot(
                    self.df['Wavenumber'], self.predicted_bands[self.y_axis]['Lorentzian'],
                    color='#ae81ff', linewidth=0.5, linestyle=':'
                )
                if self.y_axis == 'Transmittance' or self.y_axis == 'Transmittance(N)':
                    self.canvas.ax.fill_between(self.df['Wavenumber'], predicted, 1, alpha=0.2)
                elif self.y_axis == '%Transmittance' or self.y_axis == '%Transmittance(N)':
                    self.canvas.ax.fill_between(self.df['Wavenumber'], predicted, 100, alpha=0.2)
                else:
                    self.canvas.ax.fill_between(
                        self.df['Wavenumber'], predicted, 0, alpha=0.2
                    )
        if self.bands_df is not None and self.bands_df.shape[0]:
            for index, row in self.bands_df.iterrows():
                axis = self.y_axis.split('(')[0]
                if axis == 'Transmittance' or axis == '%Transmittance':
                    offset = -50
                else:
                    offset = 15
                self.canvas.ax.annotate(
                    str(round(row['Wavenumber'], 0)), xy=(row['Wavenumber'], row[self.y_axis]),
                    xytext=(0, offset), textcoords='offset points', ha='center', va='bottom',
                    rotation=90, bbox=dict(boxstyle='round', color='white', alpha=0.5, ec="white"),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                )
        self.canvas.draw()
        self.set_model()

    def set_model(self):
        if self.bands_df is not None:
            self.bands_model = IRDataModel(self.bands_df, self.y_axis)
            self.uiIRBandsTableView.setModel(self.bands_model)
            self.uiIRBandsTableView.resizeColumnsToContents()
            self.uiIRBandsTableView.resizeRowsToContents()

    def set_y_axis(self):
        """
        docstring
        """
        if 0.0 < self.df['Raw Data'].median() < 1.0:
            self.df['Transmittance'] = self.df['Raw Data']
            self.df['Absorbance'] = self.df['Transmittance'].apply(self.trans_to_abs)
            self.df['%Transmittance'] = self.df['Transmittance'].apply(self.trans_to_percent)
        elif self.df['Raw Data'].min() < 1.0 and 2 > self.df['Raw Data'].max() > 0:
            self.df['Absorbance'] = self.df['Raw Data']
            self.df['Transmittance'] = self.df['Absorbance'].apply(self.abs_to_trans)
            self.df['%Transmittance'] = self.df['Transmittance'].apply(self.trans_to_percent)
        else:
            self.df['%Transmittance'] = self.df['Raw Data']
            self.df['Transmittance'] = self.df['%Transmittance'].apply(self.trans_to_percent)
            self.df['Absorbance'] = self.df['%Transmittance'].apply(self.percent_to_abs)
        self.min_abs = self.df['Absorbance'].min()
        self.max_abs = self.df['Absorbance'].max()
        self.y_range_abs = self.max_abs - self.min_abs
        self.min_trans = self.df['Transmittance'].min()
        self.max_trans = self.df['Transmittance'].max()
        self.y_range_trans = self.max_trans - self.min_trans
        self.df['Absorbance(N)'] = self.df['Absorbance'].apply(self.norm_A)
        self.df['Transmittance(N)'] = self.df['Transmittance'].apply(self.norm_T)
        self.df['%Transmittance(N)'] = self.df['Transmittance(N)'] * 100

    def set_base_thres(self, value):
        self.base_threshold = value
        if self.df.get('baseline') is not None:
            self.baseline_detection()

    def set_bands_thres(self, value):
        self.band_threshold = value
        if self.df.get('Bands') is not None:
            del self.df['Bands']
            self.bands_detect_auto()

    def set_order_base(self, value):
        self.order_base = value
        if self.df.get('baseline') is not None:
            self.baseline_detection()

    def set_order_bands(self, value):
        self.order_bands = value
        if self.df.get('Bands') is not None:
            self.bands_detect_auto()

    def set_grid(self):
        """set_grid."""
        self.show_grid = self.uiShowGridButton.isChecked()
        self.plot()

    def set_normal(self):
        """set_normal.
        """
        self.normalized = self.uiNormalizeButton.isChecked()
        self.select_y_axis(self.y_axis) if self.normalized else self.select_y_axis(self.y_axis.replace('(N)', ''))

    def select_y_axis(self, unit):
        self.uiTransButton.setChecked(True) if unit == 'Transmittance' else self.uiTransButton.setChecked(False)
        self.uiPercentTransButton.setChecked(True) if unit == '%Transmittance' else self.uiPercentTransButton.setChecked(False)
        self.uiAbsorbanceButton.setChecked(True) if unit == 'Absorbance' else self.uiAbsorbanceButton.setChecked(False)
        self.y_axis = unit + '(N)' if self.normalized else unit
        self.plot()

    def baseline_detection(self):
        """ baseline_detection
        """
        self.baseline_df = pd.DataFrame({'Wavenumber': self.df['Wavenumber'].to_list()})
        self.baseline_df['Transmittance'] = self.df.loc[
            argrelextrema(
                self.df.Transmittance.values,
                np.greater_equal,
                order=self.order_base
            )[0], 'Transmittance'
        ]
        mask_1 = (self.baseline_df['Transmittance'] < self.base_threshold)
        self.baseline_df.loc[mask_1, 'Transmittance'] = np.nan
        mask_2 = (self.baseline_df['Transmittance'].isnull() == False)
        baseline_x = self.df.loc[mask_2, 'Wavenumber']
        baseline_y = self.df.loc[mask_2, 'Transmittance']
        coef = np.polyfit(baseline_x, baseline_y, 1)
        baseline_fit = np.poly1d(coef)
        self.baseline_df['Transmittance'] = pd.Series([baseline_fit(i) for i in self.df['Wavenumber']])
        self.baseline_df['Transmittance(N)'] = self.baseline_df['Transmittance'].apply(self.norm_T)
        self.baseline_df['%Transmittance'] = self.baseline_df['Transmittance'] * 100
        self.baseline_df['%Transmittance(N)'] = self.baseline_df['Transmittance(N)'].apply(self.trans_to_percent)
        self.baseline_df['Absorbance'] = self.baseline_df['Transmittance'].apply(self.trans_to_abs)
        self.baseline_df['Absorbance(N)'] = self.baseline_df['Absorbance'].apply(self.norm_A)
        self.baseline = True
        self.uiSubsBaselineButton.setEnabled(True)
        self.uiDetectBaselineButton.setEnabled(False)
        self.plot()

    def baseline_subs(self):
        self.df['Transmittance'] = 1 - self.baseline_df['Transmittance'] + self.df['Transmittance']
        self.df['%Transmittance'] = self.df['Transmittance'].apply(self.trans_to_percent)
        self.df['Absorbance'] = self.df['Transmittance'].apply(self.trans_to_abs)
        self.min_abs = self.df['Absorbance'].min()
        self.max_abs = self.df['Absorbance'].max()
        self.y_range_abs = self.max_abs - self.min_abs
        self.min_trans = self.df['Transmittance'].min()
        self.max_trans = self.df['Transmittance'].max()
        self.y_range_trans = self.max_trans - self.min_trans
        self.df['Transmittance(N)'] = self.df['Transmittance'].apply(self.norm_T)
        self.df['%Transmittance(N)'] = self.df['Transmittance(N)'].apply(self.trans_to_percent)
        self.df['Absorbance(N)'] = self.df['Absorbance'].apply(self.norm_A)
        self.uiSubsBaselineButton.setEnabled(False)
        self.uiBaseThresholdSpinBox.setEnabled(False)
        self.uiOrderBaseSpinBox.setEnabled(False)
        self.baseline = False
        if self.bands_df.shape[0] > 0:
            self.bands_detect_auto()
        else:
            self.plot()

    def trans_to_percent(self, y):
        return round(y * 100, 6)

    def trans_to_abs(self, y):
        return round(2.0 - math.log10(y * 100), 6)

    def PercentToTrans(self, y):
        return round(y / 100, 6)

    def percent_to_abs(self, y):
        return round(2.0 - math.log10(y), 6)

    def abs_to_trans(self, y):
        return round((10 ** (2 - y)) / 100, 6)

    def abs_to_percent(self, y):
        return round((10 ** (2 - y)), 6)

    def norm_T(self, y):
        return round((1.0 - (self.max_trans - y) / self.y_range_trans), 6)

    def norm_A(self, y):
        return round(((y - self.min_abs) / self.y_range_abs), 6)

    def restore(self):
        """
        docstring
        """
        self.uiDetectBaselineButton.setEnabled(True)
        self.uiBaseThresholdSpinBox.setEnabled(True)
        self.uiOrderBaseSpinBox.setEnabled(True)
        self.uiSubsBaselineButton.setEnabled(False)
        self.uiBandManualButton.setChecked(False)
        self.uiBandFittingAutoButton.setEnabled(False)
        self.uiRemoveBandButton.setEnabled(False)
        self.uiNormalizeButton.setChecked(False)
        self.bands_fitted = False
        self.normalized = False
        if self.df.get('Bands') is not None:
            del self.df['Bands']
        self.bands_df = pd.DataFrame({
            'Wavenumber': [],
            'Transmittance': [],
            '%Transmittance': [],
            'Absorbance': [],
            'Transmittance(N)': [],
            '%Transmittance(N)': [],
            'Absorbance(N)': [],
        })
        self.predicted_bands = {}
        self.set_y_axis()
        self.select_y_axis('%Transmittance(N)')

    def bands_detect_auto(self):
        """
        docstring
        """
        self.df['Bands'] = self.df.loc[
            argrelextrema(
                self.df['%Transmittance(N)'].values, np.less_equal,
                order=self.order_bands
            )[0], '%Transmittance(N)'
        ]
        mask3 = (self.df['Bands'] < self.band_threshold)
        self.df.loc[mask3, 'Bands'] = np.nan
        self.bands_df = self.df.loc[
            mask3, [
                'Wavenumber', 'Transmittance', '%Transmittance', 'Absorbance',
                'Transmittance(N)', '%Transmittance(N)', 'Absorbance(N)'
            ]
        ].copy()
        self.bands_df.reset_index(drop=True, inplace=True)
        self.bands_df['HWHM'] = pd.Series(
            map(
                self.HWHM, self.bands_df['Wavenumber'],
                self.bands_df['%Transmittance(N)']
            )
        )
        self.bands_df = self.bands_df.round(6)
        self.uiBandFittingAutoButton.setEnabled(True)
        self.uiRemoveBandButton.setEnabled(True)
        self.uiToolBox.setCurrentIndex(3)
        self.plot()

    def bands_detect_manual(self, event):
        """
        docstring
        """
        if self.uiBandManualButton.isChecked():
            ix, _ = event.xdata, event.ydata
            tolerance = 1  # this value can be moved to avoid errors if the x
            # sample data changes.
            mask = self.df['Wavenumber'].between(ix - tolerance, ix + tolerance)
            try:
                transmittance = self.df.loc[mask, 'Transmittance'].values[0]
                transmittance_N = self.df.loc[mask, 'Transmittance(N)'].values[0]
                percent = self.df.loc[mask, '%Transmittance'].values[0]
                percent_N = self.df.loc[mask, '%Transmittance(N)'].values[0]
                absorbance = self.df.loc[mask, 'Absorbance'].values[0]
                absorbance_N = self.df.loc[mask, 'Absorbance(N)'].values[0]
                new_band = pd.DataFrame({
                    'Wavenumber': [ix],
                    'Transmittance': [transmittance],
                    '%Transmittance': [percent],
                    'Absorbance': [absorbance],
                    'Transmittance(N)': [transmittance_N],
                    '%Transmittance(N)': [percent_N],
                    'Absorbance(N)': [absorbance_N],
                    'HWHM': [self.HWHM(ix, percent_N)]
                })
                if self.bands_df is None:
                    self.bands_df = pd.DataFrame({
                        'Wavenumber': [],
                        'Transmittance': [],
                        '%Transmittance': [],
                        'Absorbance': [],
                        'Transmittance(N)': [],
                        '%Transmittance(N)': [],
                        'Absorbance(N)': [],
                    })
                self.bands_df = self.bands_df.append(new_band, ignore_index=True)
                self.bands_df = self.bands_df.round(6)
            except IndexError:
                pass
            self.uiBandFittingAutoButton.setEnabled(True)
            self.uiRemoveBandButton.setEnabled(True)
            self.plot()

    def remove_band(self):
        indexes = self.uiIRBandsTableView.selectedIndexes()
        if indexes is not None:
            self.bands_df.drop([i.row() for i in indexes], inplace=True)
            self.bands_df.reset_index(drop=True, inplace=True)
            self.set_model()
            self.plot()

        if self.bands_fitted:
            reply = qtw.QMessageBox.question(
                self, 'Redo Fitting',
                'You have removed a band, do you want to redo the fitting?',
                qtw.QMessageBox.Yes | qtw.QMessageBox.No,
                qtw.QMessageBox.No
            )

            if reply == qtw.QMessageBox.Yes:
                self.bands_fitting()

    def HWHM(self, x, y):
        """HWHM.

        Parameters
        ----------
        x :
            wavenumber value
        y :
            normalized transmittance value
        """
        halfMax = round((100 + y) / 2, 6)

        mask = self.df['Wavenumber'].between(x, x + 10)
        mask_2 = self.df['Wavenumber'].between(x - 10, x)

        bandsToFitX = self.df.loc[mask, 'Wavenumber'].values
        bandsToFitY = self.df.loc[mask, '%Transmittance(N)'].values
        bandsToFitX2 = self.df.loc[mask_2, 'Wavenumber'].values
        bandsToFitY2 = self.df.loc[mask_2, '%Transmittance(N)'].values

        coef = np.polyfit(bandsToFitX, bandsToFitY, 1)
        bandSteepFit = np.poly1d(coef)
        coef2 = np.polyfit(bandsToFitX2, bandsToFitY2, 1)
        bandSteepFit2 = np.poly1d(coef2)

        xRange = np.linspace(x - 500, x + 500, 100000)
        fitted = pd.DataFrame({'x': xRange, 'y1': bandSteepFit(xRange), 'y2': bandSteepFit2(xRange)})

        tolerance = 1e-1

        mask3 = fitted['y1'].between(halfMax - tolerance, halfMax + tolerance)
        mask4 = fitted['y2'].between(halfMax - tolerance, halfMax + tolerance)

        found = fitted.loc[mask3, 'x'].values
        found2 = fitted.loc[mask4, 'x'].values
        if len(found) > 0 and len(found2) > 0:
            hwhm1 = abs(x - found[0])
            hwhm2 = abs(x - found2[0])
        else:
            tolerance = 1
            mask3 = fitted['y1'].between(
                halfMax - tolerance, halfMax + tolerance
            )
            mask4 = fitted['y2'].between(halfMax - tolerance, halfMax + tolerance)
            found = fitted.loc[mask3, 'x'].values
            found2 = fitted.loc[mask4, 'x'].values
            if len(found) > 0 and len(found2) > 0:
                hwhm1 = abs(x - found[0])
                hwhm2 = abs(x - found2[0])
            else:
                return 1.000

        mean = (hwhm1 + hwhm2) / 2
        return mean

    def bands_fitting(self):
        """
        docstring
        """
        self.uiBandFittingAutoButton.setEnabled(False)
        self.uiRemoveBandButton.setEnabled(False)
        self.worker = IRWorkerThread(self.bands_df, self.df)
        self.thread = qtc.QThread()
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.fit)
        self.worker.finished.connect(self.thread.quit)
        self.statusBar().showMessage(f'Now fitting {self.title} IR spectrum')
        self.worker.ok_fit.connect(self.ok_fit_manager)
        self.worker.error.connect(self.error_fit_manager)
        self.thread.start()

    @qtc.pyqtSlot(dict, object, object)
    def ok_fit_manager(self, predicted_bands, bands, df):
        self.predicted_bands = predicted_bands
        self.df = df
        self.bands_df = bands
        self.bands_df = self.bands_df.round(6)
        self.predicted_bands['Absorbance(N)'] = self.predicted_bands['Absorbance'].applymap(self.norm_A)
        self.predicted_bands['Transmittance'] = self.predicted_bands['Absorbance'].applymap(self.abs_to_trans)
        self.predicted_bands['Transmittance(N)'] = self.predicted_bands['Transmittance'].applymap(self.norm_T)
        self.predicted_bands['%Transmittance'] = self.predicted_bands['Transmittance'].applymap(self.trans_to_percent)
        self.predicted_bands['%Transmittance(N)'] = self.predicted_bands['Transmittance(N)'].applymap(lambda x: x * 100)
        self.predicted_bands['Absorbance']['Lorentzian'] = self.predicted_bands['Absorbance'].sum(axis=1)
        self.predicted_bands['Absorbance(N)']['Lorentzian'] = self.predicted_bands['Absorbance']['Lorentzian'].apply(self.norm_A)
        self.predicted_bands['Transmittance']['Lorentzian'] = self.predicted_bands['Absorbance']['Lorentzian'].apply(self.abs_to_trans)
        self.predicted_bands['Transmittance(N)']['Lorentzian'] = self.predicted_bands['Transmittance']['Lorentzian'].apply(self.norm_T)
        self.predicted_bands['%Transmittance']['Lorentzian'] = self.predicted_bands['Transmittance']['Lorentzian'] * 100
        self.predicted_bands['%Transmittance(N)']['Lorentzian'] = self.predicted_bands['Transmittance(N)']['Lorentzian'] * 100
        self.predicted_bands['Residuals'] = self.df['Absorbance'] - self.predicted_bands['Absorbance']['Lorentzian']
        self.statusBar().showMessage(f'Successful fitting of {self.title}')
        self.bands_fitted = True
        self.plot()

    @qtc.pyqtSlot()
    def error_fit_manager(self):
        self.statusBar().showMessage(f'Unsuccessful fitting for {self.title}')
        self.uiBandFittingAutoButton.setEnabled(True)
        self.uiRemoveBandButton.setEnabled(True)

    def save_data(self):
        filename, _ = qtw.QFileDialog.getSaveFileName(
            self, 'Save bands report in txt format',
            options=qtw.QFileDialog.DontUseNativeDialog,
            filter='Txt file (*.txt)',
            directory=self.cwd
        )
        if filename:
            new_filename = filename + '.txt' if not '.txt' in filename else filename
            general_info = f'Bands report for: {self.title} ({self.solvent})\n'
            with open(new_filename, 'w') as file:
                df_str = self.bands_df.to_string(header=True, index=True)
                file.write(general_info + df_str)

    # def send_data(self):
        # date = datetime.now()
        # bands_df = self.bands_df if self.bands_df.shape[0] else None
        # predicted = self.predicted_bands if self.predicted_bands.items() else None
        # self.data_sent.emit({
            # 'TYPE': 'FTIR',
            # 'DATE': str(date),
            # 'SOLVENT': self.solvent,
            # 'DATA FRAME': self.df,
            # 'BANDS': bands_df,
            # 'PREDICTED': predicted
        # })
        # reply = qtw.QMessageBox.question(
            # self, 'Window Close',
            # 'The spectrum has been stored, do you want to close the window?',
            # qtw.QMessageBox.Yes | qtw.QMessageBox.No,
            # qtw.QMessageBox.No
        # )

        # if reply == qtw.QMessageBox.Yes:
            # self.close()

    def keyPressEvent(self, event):
        ctrl_pressed = (event.modifiers() & qtc.Qt.ControlModifier)
        selected = self.uiIRBandsTableView.selectedIndexes()
        if event.key() == qtc.Qt.Key_C and ctrl_pressed:
            self.copied_cells = sorted(selected)
        elif event.key() == qtc.Qt.Key_V and ctrl_pressed and self.copied_cells:
            for cell in self.copied_cells:
                print(cell.data())
        elif event.key() == qtc.Qt.Key_M:
            if self.uiBandManualButton.isChecked():
                self.uiBandManualButton.setChecked(False)
                self.plot()
            else:
                self.uiBandManualButton.setChecked(True)
                self.plot()
        elif event.key() == qtc.Qt.Key_A:
            self.bands_detect_auto()
        elif event.key() == qtc.Qt.Key_Delete and selected:
            self.remove_band()
        if event.key() == qtc.Qt.Key_S and ctrl_pressed:
            self.save_data()
        if event.key() == qtc.Qt.Key_O and ctrl_pressed:
            self.open_file()


