#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from PyQt5 import QtWidgets as qtw
from IRController import IRPlotter

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    window = IRPlotter(
        csv_path='/home/samuelvip/Documentos/Bayreuth/IR/JAN2IM.csv'
    )
    window.show()
    sys.exit(app.exec_())
