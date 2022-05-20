# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'station_gui.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.resize(1114, 764)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setMinimumSize(QtCore.QSize(0, 40))
        self.frame.setMaximumSize(QtCore.QSize(16777214, 40))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.frame)
        
        self.backwardButton = QtWidgets.QPushButton(self.frame)
        self.backwardButton.setText("<<<")
        self.horizontalLayout.addWidget(self.backwardButton)
        
        self.intervalButton = QtWidgets.QPushButton(self.frame)
        self.intervalButton.setText("Step")
        self.horizontalLayout.addWidget(self.intervalButton)
        
        self.forwardButton = QtWidgets.QPushButton(self.frame)
        self.forwardButton.setText(">>>")
        self.horizontalLayout.addWidget(self.forwardButton)
        
        self.gotoButton = QtWidgets.QPushButton(self.frame)
        self.gotoButton.setIconSize(QtCore.QSize(20, 20))
        self.gotoButton.setText("Go To")
        self.horizontalLayout.addWidget(self.gotoButton)

        self.filtButton = QtWidgets.QPushButton(self.frame)
        self.filtButton.setIconSize(QtCore.QSize(20, 20))
        self.filtButton.setText("FILT")
        self.horizontalLayout.addWidget(self.filtButton)

        self.specButton = QtWidgets.QPushButton(self.frame)
        self.specButton.setIconSize(QtCore.QSize(20, 20))
        self.specButton.setText("SPEC")
        self.horizontalLayout.addWidget(self.specButton)

        self.respButton = QtWidgets.QPushButton(self.frame)
        self.respButton.setIconSize(QtCore.QSize(20, 20))
        self.respButton.setText("RESP")
        self.horizontalLayout.addWidget(self.respButton)
        
        self.saveButton = QtWidgets.QPushButton(self.frame)
        self.saveButton.setText("USGS")
        self.horizontalLayout.addWidget(self.saveButton)
        
        self.psdButton = QtWidgets.QPushButton(self.frame)
        self.psdButton.setIconSize(QtCore.QSize(16, 16))
        self.psdButton.setText("PSD")
        self.horizontalLayout.addWidget(self.psdButton)
        
        self.exportButton = QtWidgets.QPushButton(self.frame)
        self.exportButton.setIconSize(QtCore.QSize(16, 16))
        self.exportButton.setText("Export")
        self.horizontalLayout.addWidget(self.exportButton)
        
        # more buttons
        self.toolButton = QtWidgets.QToolButton(self.frame)
        self.horizontalLayout.addWidget(self.toolButton)
        spacerItem = QtWidgets.QSpacerItem(452, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.gridLayout.addWidget(self.frame, 0, 0, 1, 2)
        self.widget = QtWidgets.QWidget(self.centralwidget)
        self.widget.setMinimumSize(QtCore.QSize(200, 16777215))
        self.widget.setMaximumSize(QtCore.QSize(200, 16777215))
        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget)
        self.tabWidget = QtWidgets.QTabWidget(self.widget)
        self.infotab = QtWidgets.QWidget()
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.infotab)
        self.plainTextEdit = QtWidgets.QTextBrowser(self.infotab)
        self.verticalLayout_2.addWidget(self.plainTextEdit)
        self.tabWidget.addTab(self.infotab, "")
        self.historytab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.historytab, "")
        self.verticalLayout.addWidget(self.tabWidget)
        self.groupBox = QtWidgets.QGroupBox(self.widget)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.groupBox)
        self.listWidget = QtWidgets.QListWidget(self.groupBox)
        self.horizontalLayout_2.addWidget(self.listWidget)
        self.verticalLayout.addWidget(self.groupBox)
        self.gridLayout.addWidget(self.widget, 1, 0, 1, 1)
        # self.widget_2 = QtWidgets.QWidget(self.centralwidget)
        # self.widget_2.setObjectName("widget_2")
        # self.gridLayout.addWidget(self.widget_2, 1, 1, 1, 1)
        
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1114, 22))
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuDatabase = QtWidgets.QMenu(self.menubar)
        self.menuSettings = QtWidgets.QMenu(self.menubar)
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        MainWindow.setStatusBar(self.statusbar)
        self.actionForward = QtWidgets.QAction(MainWindow)
        self.actionBackward = QtWidgets.QAction(MainWindow)
        self.actionGo_to = QtWidgets.QAction(MainWindow)
        self.actionSet_interval = QtWidgets.QAction(MainWindow)
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionOpen = QtWidgets.QAction(MainWindow)
        self.actionSave_event = QtWidgets.QAction(MainWindow)
        self.actionShow_info = QtWidgets.QAction(MainWindow)
        self.actionConfig = QtWidgets.QAction(MainWindow)
        self.actionPreferences = QtWidgets.QAction(MainWindow)
        self.actionRemove_Response = QtWidgets.QAction(MainWindow)
        self.actionRemove_Response.setCheckable(True)

        self.actionMatrix_Return = QtWidgets.QAction(MainWindow)
        self.actionMatrix_Return.setCheckable(True)

        self.action3component = QtWidgets.QAction(MainWindow)
        self.action3component.setCheckable(True)
        self.actionPolargram = QtWidgets.QAction(MainWindow)
        self.actionPolargram.setCheckable(True)
        self.actionSpectrogram = QtWidgets.QAction(MainWindow)
        self.actionSpectrogram.setCheckable(True)
        
        self.menuFile.addAction(self.actionForward)
        self.menuFile.addAction(self.actionBackward)
        self.menuFile.addAction(self.actionGo_to)
        self.menuFile.addAction(self.actionSet_interval)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionQuit)
        self.menuDatabase.addAction(self.actionOpen)
        self.menuDatabase.addAction(self.actionShow_info)
        self.menuDatabase.addAction(self.actionSave_event)
        self.menuSettings.addAction(self.action3component)
        self.menuSettings.addAction(self.actionRemove_Response)
        self.menuSettings.addAction(self.actionSpectrogram)
        self.menuSettings.addAction(self.actionPolargram)
        self.menuSettings.addAction(self.actionMatrix_Return)
        self.menuSettings.addAction(self.actionPreferences)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuDatabase.menuAction())
        self.menubar.addAction(self.menuSettings.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.gotoButton.setText(_translate("MainWindow", "Go to..."))
        self.saveButton.setText(_translate("MainWindow", "USGS"))
        self.toolButton.setText(_translate("MainWindow", "..."))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.infotab), _translate("MainWindow", "Info"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.historytab), _translate("MainWindow", "Historial"))
        self.groupBox.setTitle(_translate("MainWindow", "Events"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuDatabase.setTitle(_translate("MainWindow", "Database"))
        self.menuSettings.setTitle(_translate("MainWindow", "Settings"))
        self.actionForward.setText(_translate("MainWindow", "Forward"))
        self.actionBackward.setText(_translate("MainWindow", "Backward"))
        self.actionGo_to.setText(_translate("MainWindow", "Go to..."))
        self.actionSet_interval.setText(_translate("MainWindow", "Set interval"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionSave_event.setText(_translate("MainWindow", "Save event"))
        self.actionShow_info.setText(_translate("MainWindow", "Show info"))
        self.actionConfig.setText(_translate("MainWindow", "Config"))
        self.actionPreferences.setText(_translate("MainWindow", "Settings..."))
        self.actionRemove_Response.setText(
            _translate("MainWindow", "Remove Response"))
        self.action3component.setText(
            _translate("MainWindow", "Z-N-E component"))
        self.actionPolargram.setText(
            _translate("MainWindow", "Polargram"))
        self.actionMatrix_Return.setText(
            _translate("MainWindow", "Polargram Extend"))
        self.actionSpectrogram.setText(
            _translate("MainWindow", "Spectrogram"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())