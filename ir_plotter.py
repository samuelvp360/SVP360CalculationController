#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from PyQt5 import QtWidgets as qtw
from IRController import IRPlotter

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = IRPlotter()
    window.show()
    sys.exit(app.exec_())
