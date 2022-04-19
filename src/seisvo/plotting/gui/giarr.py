#!/usr/bin/python3
# coding=utf-8


import numpy as np
from datetime import timedelta
import pyqtgraph
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvas
from obspy import UTCDateTime
from seisvo import Array, Trace2, Stream2, DB
from seisvo.signal.spectrum import power_density_spectrum
from seisvo.gui import Navigation
from seisvo.gui import notify
from seisvo.utils.plotting import get_colors


plt.rc('axes', labelsize=10)
plt.rc('axes', labelpad=4.0)
plt.rc('axes', titlepad=6.0)
plt.rc('axes', titlesize=10)
plt.rc('xtick', labelsize=10)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('lines', linewidth=0.5)
plt.rc('lines', linewidth=0.5)


icons_path = '/home/ifm/.local/lib/python3.8/site-packages/seisvo/gui/icons'


def gplot(array=None, air_file=None, **kwargs):
    qApp = QtWidgets.QApplication([])

    if array:
        app = iArrayGui(array=array, **kwargs)

    elif air_file:
        app = iArrayGui(air_file=air_file, **kwargs)

    else:
        app = None

    if app:
        app.show()
        qApp.exec_()


def removeID(db, parent=None):
    id_list = db.get_id_list()

    if not id_list:
        notify('seisvo', 'No events in database', status='error')
    else:
        id, ok = QtWidgets.QInputDialog.getInt(parent, "Select ID", "Enter ID to remove",
                                               max(id_list), min(id_list), max(id_list), 1)
        if ok:
            try:
                db.remove_id(id)
                notify('seisvo', 'ID %s removed!' % id, status='info')
            except:
                notify('seisvo', 'Error removing ID', status='error')


class Ui_array(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1200, 700)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tool_widget = QtWidgets.QWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tool_widget.sizePolicy().hasHeightForWidth())
        self.tool_widget.setSizePolicy(sizePolicy)
        self.tool_widget.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.tool_widget.setObjectName("tool_widget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.tool_widget)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.frame = QtWidgets.QFrame(self.tool_widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame.sizePolicy().hasHeightForWidth())
        self.frame.setSizePolicy(sizePolicy)
        self.frame.setMinimumSize(QtCore.QSize(900, 74))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.time_label = QtWidgets.QLabel(self.frame)
        self.time_label.setGeometry(QtCore.QRect(180, 10, 60, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.time_label.sizePolicy().hasHeightForWidth())
        self.time_label.setSizePolicy(sizePolicy)
        self.time_label.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.time_label.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.time_label.setFrameShadow(QtWidgets.QFrame.Plain)
        self.time_label.setTextFormat(QtCore.Qt.AutoText)
        self.time_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.time_label.setObjectName("time_label")
        self.timeBox = QtWidgets.QSpinBox(self.frame)
        self.timeBox.setGeometry(QtCore.QRect(250, 10, 51, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.timeBox.sizePolicy().hasHeightForWidth())
        self.timeBox.setSizePolicy(sizePolicy)
        self.timeBox.setMinimum(0)
        self.timeBox.setObjectName("timeBox")
        self.azimuth_label = QtWidgets.QLabel(self.frame)
        self.azimuth_label.setGeometry(QtCore.QRect(500, 10, 56, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.azimuth_label.sizePolicy().hasHeightForWidth())
        self.azimuth_label.setSizePolicy(sizePolicy)
        self.azimuth_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.azimuth_label.setObjectName("azimuth_label")
        self.azimuthBox = QtWidgets.QSpinBox(self.frame)
        self.azimuthBox.setGeometry(QtCore.QRect(560, 10, 54, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.azimuthBox.sizePolicy().hasHeightForWidth())
        self.azimuthBox.setSizePolicy(sizePolicy)
        self.azimuthBox.setMaximum(360)
        self.azimuthBox.setObjectName("azimuthBox")
        self.azm_cheqBox = QtWidgets.QCheckBox(self.frame)
        self.azm_cheqBox.setGeometry(QtCore.QRect(350, 10, 110, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.azm_cheqBox.sizePolicy().hasHeightForWidth())
        self.azm_cheqBox.setSizePolicy(sizePolicy)
        self.azm_cheqBox.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.azm_cheqBox.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.azm_cheqBox.setChecked(True)
        self.azm_cheqBox.setTristate(False)
        self.azm_cheqBox.setObjectName("azm_cheqBox")
        self.plotButton = QtWidgets.QPushButton(self.frame)
        self.plotButton.setGeometry(QtCore.QRect(800, 10, 71, 54))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotButton.sizePolicy().hasHeightForWidth())
        self.plotButton.setSizePolicy(sizePolicy)
        self.plotButton.setMinimumSize(QtCore.QSize(71, 54))
        self.plotButton.setToolTip("")
        self.plotButton.setText("")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("%s/plot.svg" % icons_path), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.plotButton.setIcon(icon)
        self.plotButton.setIconSize(QtCore.QSize(24, 24))
        self.plotButton.setAutoDefault(False)
        self.plotButton.setDefault(False)
        self.plotButton.setFlat(False)
        self.plotButton.setObjectName("plotButton")
        self.label = QtWidgets.QLabel(self.frame)
        self.label.setGeometry(QtCore.QRect(640, 20, 36, 31))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.frame)
        self.label_2.setGeometry(QtCore.QRect(680, 40, 28, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName("label_2")
        self.freqminSpin = QtWidgets.QDoubleSpinBox(self.frame)
        self.freqminSpin.setGeometry(QtCore.QRect(720, 40, 54, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.freqminSpin.sizePolicy().hasHeightForWidth())
        self.freqminSpin.setSizePolicy(sizePolicy)
        self.freqminSpin.setDecimals(1)
        self.freqminSpin.setMaximum(50.0)
        self.freqminSpin.setSingleStep(0.1)
        self.freqminSpin.setProperty("value", 0.5)
        self.freqminSpin.setObjectName("freqminSpin")
        self.label_3 = QtWidgets.QLabel(self.frame)
        self.label_3.setGeometry(QtCore.QRect(680, 10, 32, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setObjectName("label_3")
        self.freqmaxSpin = QtWidgets.QDoubleSpinBox(self.frame)
        self.freqmaxSpin.setGeometry(QtCore.QRect(720, 10, 54, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.freqmaxSpin.sizePolicy().hasHeightForWidth())
        self.freqmaxSpin.setSizePolicy(sizePolicy)
        self.freqmaxSpin.setDecimals(1)
        self.freqmaxSpin.setMaximum(50.0)
        self.freqmaxSpin.setSingleStep(0.1)
        self.freqmaxSpin.setProperty("value", 5.0)
        self.freqmaxSpin.setObjectName("freqmaxSpin")
        self.label_4 = QtWidgets.QLabel(self.frame)
        self.label_4.setGeometry(QtCore.QRect(180, 40, 63, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        self.label_4.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_4.setObjectName("label_4")
        self.timepadBox = QtWidgets.QSpinBox(self.frame)
        self.timepadBox.setGeometry(QtCore.QRect(250, 40, 61, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.timepadBox.sizePolicy().hasHeightForWidth())
        self.timepadBox.setSizePolicy(sizePolicy)
        self.timepadBox.setMinimum(5)
        self.timepadBox.setMaximum(3600)
        self.timepadBox.setProperty("value", 30)
        self.timepadBox.setObjectName("timepadBox")
        self.retrocess_button = QtWidgets.QPushButton(self.frame)
        self.retrocess_button.setGeometry(QtCore.QRect(20, 10, 61, 51))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.retrocess_button.sizePolicy().hasHeightForWidth())
        self.retrocess_button.setSizePolicy(sizePolicy)
        self.retrocess_button.setMinimumSize(QtCore.QSize(40, 26))
        self.retrocess_button.setText("")
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("%s/left_arr.svg" % icons_path), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.retrocess_button.setIcon(icon1)
        self.retrocess_button.setObjectName("retrocess_button")
        self.advance_button = QtWidgets.QPushButton(self.frame)
        self.advance_button.setGeometry(QtCore.QRect(100, 10, 61, 51))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.advance_button.sizePolicy().hasHeightForWidth())
        self.advance_button.setSizePolicy(sizePolicy)
        self.advance_button.setMinimumSize(QtCore.QSize(61, 26))
        self.advance_button.setText("")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("%s/right_arr.svg" % icons_path), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.advance_button.setIcon(icon2)
        self.advance_button.setObjectName("advance_button")
        self.smb_label = QtWidgets.QLabel(self.frame)
        self.smb_label.setGeometry(QtCore.QRect(340, 40, 61, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.smb_label.sizePolicy().hasHeightForWidth())
        self.smb_label.setSizePolicy(sizePolicy)
        self.smb_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.smb_label.setObjectName("smb_label")
        self.smb_spinBox = QtWidgets.QDoubleSpinBox(self.frame)
        self.smb_spinBox.setGeometry(QtCore.QRect(410, 40, 54, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.smb_spinBox.sizePolicy().hasHeightForWidth())
        self.smb_spinBox.setSizePolicy(sizePolicy)
        self.smb_spinBox.setDecimals(1)
        self.smb_spinBox.setMaximum(1.0)
        self.smb_spinBox.setSingleStep(0.1)
        self.smb_spinBox.setProperty("value", 0.5)
        self.smb_spinBox.setObjectName("smb_spinBox")
        self.correlation_label = QtWidgets.QLabel(self.frame)
        self.correlation_label.setGeometry(QtCore.QRect(475, 40, 81, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.correlation_label.sizePolicy().hasHeightForWidth())
        self.correlation_label.setSizePolicy(sizePolicy)
        self.correlation_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.correlation_label.setObjectName("correlation_label")
        self.correlationSpin = QtWidgets.QDoubleSpinBox(self.frame)
        self.correlationSpin.setGeometry(QtCore.QRect(560, 40, 54, 21))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.correlationSpin.sizePolicy().hasHeightForWidth())
        self.correlationSpin.setSizePolicy(sizePolicy)
        self.correlationSpin.setDecimals(2)
        self.correlationSpin.setMaximum(1.00)
        self.correlationSpin.setSingleStep(0.1)
        self.correlationSpin.setProperty("value", 0.0)
        self.correlationSpin.setObjectName("correlationSpin")
        self.correlationSpin.setEnabled(False)
        self.gridLayout_2.addWidget(self.frame, 0, 0, 1, 1)
        self.verticalLayout.addWidget(self.tool_widget)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 934, 25))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuDB = QtWidgets.QMenu(self.menubar)
        self.menuDB.setObjectName("menuDB")
        self.menuView = QtWidgets.QMenu(self.menubar)
        self.menuView.setObjectName("menuView")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionInfo = QtWidgets.QAction(MainWindow)
        self.actionInfo.setObjectName("actionInfo")
        self.actionClose = QtWidgets.QAction(MainWindow)
        self.actionClose.setObjectName("actionClose")
        self.actionPrint_image = QtWidgets.QAction(MainWindow)
        self.actionPrint_image.setObjectName("actionPrint_image")
        self.actionOpen_DB = QtWidgets.QAction(MainWindow)
        self.actionOpen_DB.setObjectName("actionOpen_DB")
        self.actionRemoveDB = QtWidgets.QAction(MainWindow)
        self.actionRemoveDB.setObjectName("actionOpen_DB")
        self.actionModel = QtWidgets.QAction(MainWindow)
        self.actionModel.setObjectName("actionModel_2")
        self.menuFile.addAction(self.actionOpen_DB)
        self.menuFile.addAction(self.actionPrint_image)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionClose)
        self.menuDB.addAction(self.actionRemoveDB)
        self.menuView.addAction(self.actionInfo)
        self.menuView.addAction(self.actionModel)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuDB.menuAction())
        self.menubar.addAction(self.menuView.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.time_label.setText(_translate("MainWindow", "Time Bin:"))
        self.azimuth_label.setText(_translate("MainWindow", "Azimuth:"))
        self.azm_cheqBox.setText(_translate("MainWindow", "max. azimuth"))
        self.label.setText(_translate("MainWindow", "Filter:"))
        self.label_2.setText(_translate("MainWindow", "min."))
        self.label_3.setText(_translate("MainWindow", "max."))
        self.label_4.setText(_translate("MainWindow", "Time pad:"))
        self.retrocess_button.setShortcut(_translate("MainWindow", "Left"))
        self.advance_button.setShortcut(_translate("MainWindow", "Right"))
        self.smb_label.setText(_translate("MainWindow", "SMB min:"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuView.setTitle(_translate("MainWindow", "View"))
        self.actionInfo.setText(_translate("MainWindow", "Info..."))
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.actionPrint_image.setText(_translate("MainWindow", "Print image..."))
        self.actionOpen_DB.setText(_translate("MainWindow", "Open DB"))
        self.menuDB.setTitle(_translate("MainWindow", "Database"))
        self.actionRemoveDB.setText(_translate("MainWindow", "Remove ID"))
        self.actionModel.setText(_translate("MainWindow", "Model..."))
        self.correlation_label.setText(_translate("MainWindow", "Correlation:"))


class Ui_wvfm(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1000, 700)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 934, 25))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionClose = QtWidgets.QAction(MainWindow)
        self.actionClose.setObjectName("actionClose")
        self.actionPrint_image = QtWidgets.QAction(MainWindow)
        self.actionPrint_image.setObjectName("actionPrint_image")
        self.menuFile.addAction(self.actionPrint_image)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionClose)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "iWvfm Gui"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.actionPrint_image.setText(_translate("MainWindow", "Print image..."))


class iArrayGui(QtWidgets.QMainWindow, Ui_array):
    def __init__(self, array=None, air_file=None, **kwargs):
        """
        kwargs info
            -- general --
                starttime: datetime object
                length_interval: int in min. default 15
                pad_time: int in sec. default 30
                maxSMB: float 0--1. default 0.5
            -- iarr --
                time_width:
                overlap:
                fq_band:
        """

        if not array and not air_file:
            raise ValueError (' array or air_file must be defined')
        
        if array and air_file:
            print( 'warn: iArrayGui initialized with array')

        super(iArrayGui, self).__init__()
        self.setupUi(self)
        self.setWindowTitle('iArray Gui')

        self.air_file = air_file
        self.starttime = kwargs.get('starttime', None)
        self.length_interval = kwargs.get('interval', 15)
        self.maxSMB = kwargs.get('maxSMB', 0.5)
        self.pad_time = kwargs.get('pad_time', 30)
        self.wvfm_main = QtWidgets.QMainWindow()
        self.wvfm_app = None
        self.fq_band = kwargs.get('fq_band',(0.5,5))

        if array:
            self.iarr = array
            self.air_file = False

            time_width = kwargs.get('time_width', None)
            overlap = kwargs.get('overlap', None)

            if not time_width or not overlap:
                raise ValueError(' time_width and overlap must be stated!')

            self.time_width = time_width
            self.overlap = overlap
            self.prefilt = kwargs.get('fq_band', (0.5, 5))

            self.azimuths = range(self.iarr.model.get('src_dgr')[0],
                                  self.iarr.model.get('src_dgr')[1]+1)

            if not self.starttime:
                raise ValueError('starttime must be stated!')

            self._check_database()
            self.initFrames()
            self.compute_ccorr()

        if self.air_file:
            if not self.starttime:
                self.starttime = self.air_file.stats.starttime

            endtime = self.starttime + timedelta(minutes=self.length_interval)
            self.out = self.air_file.get_dict(starttime=self.starttime, endtime=endtime)
            self.time_bin_max = len(self.out['time'])
            self.time_width = self.air_file.stats.time_width
            self.overlap = self.air_file.stats.overlap
            self.prefilt = self.air_file.stats.fq_band

            # loading array
            net_code = self.air_file.stats.code.split('.')[0]
            sta_code = self.air_file.stats.code.split('.')[1]
            self.iarr = Array(net_code, sta_code=sta_code, model=self.out['model'])
            self.azimuths = range(self.air_file.stats.src_dgr[0],
                                  self.air_file.stats.src_dgr[1]+1)

            self._check_database()
            self.initFrames()
            self.set_default()
            self.plot()


    def _check_database(self):
        if self.iarr.db_infrasound is None:
            net_code = self.iarr.info.code
            self.iarr.db_infrasound = DB.create_database(net_code, arr=True)
            notify('database', 'created >> infrasound.sql', status='info')

        self.db = self.iarr.db_infrasound


    def initFrames(self):
        self.actionClose.triggered.connect(self.close)
        self.actionClose.setShortcut('Esc')

        self.actionModel.triggered.connect(self.show_model)
        self.actionModel.setShortcut('Alt+M')

        self.actionInfo.triggered.connect(self.show_info)
        self.actionInfo.setShortcut('Alt+I')

        self.actionPrint_image.triggered.connect(self.save_fig)
        self.actionPrint_image.setShortcut('Alt+P')

        self.actionOpen_DB.triggered.connect(self.openDB)
        self.actionOpen_DB.setShortcut('F1')

        self.actionRemoveDB.triggered.connect(self.removeID)
        self.actionRemoveDB.setToolTip('Remove ID from database')
        self.actionRemoveDB.setShortcut('DEL')

        if self.db is None:
            self.actionOpen_DB.setEnabled(False)
            self.actionRemoveDB.setEnabled(False)

        self.retrocess_button.clicked.connect(self.__retrocess)
        self.retrocess_button.setToolTip('Retrocede')
        self.advance_button.clicked.connect(self.__advance)
        self.advance_button.setToolTip('Advance')
        self.plotButton.clicked.connect(self.__plot_waveform)
        self.plotButton.setToolTip('Plot waveform')

        self.timeBox.valueChanged.connect(self.set_maxazm)
        self.azm_cheqBox.stateChanged.connect(self.set_maxazm)
        self.smb_spinBox.valueChanged.connect(self.update_plot)

        self.freqminSpin.setValue(self.fq_band[0])
        self.freqmaxSpin.setValue(self.fq_band[1])

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setParent(self.centralwidget)
        self.verticalLayout.addWidget(self.canvas)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.canvas.setSizePolicy(sizePolicy)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.callbacks.connect('button_press_event', self.on_click)

        self.starttime_info_bar = QtWidgets.QLabel()
        self.statusbar.addWidget(self.starttime_info_bar)


    def compute_ccorr(self):
        start_time = self.starttime
        with pyqtgraph.BusyCursor():
            self.out = self.iarr.ripepe_ccorr(
                self.time_width, self.overlap, start_time, interval=self.length_interval, fq_band=self.prefilt)
            self.time_bin_max = len(self.out['time'])

        self.set_default()
        self.plot()


    def azimuthFormatter(self, x, pos):
        try:
            val = self.azimuths[int(x)]
            return str(val)
        except:
            pass


    def plot(self):
        import matplotlib.ticker as mtick
        import matplotlib.colors as mcolor

        ans = self.iarr.get_max_values(self.out)
        time_bin = self.out['time']

        self.figure.clear()
        self.figure.suptitle('%s   --   %s' % (time_bin[0], time_bin[-1]))

        ax1 = self.figure.add_axes([0.1, 0.75, 0.8, 0.15])
        ax2 = self.figure.add_axes([0.1, 0.55, 0.8, 0.15])
        ax3 = self.figure.add_axes([0.1, 0.35, 0.8, 0.15])
        ax4 = self.figure.add_axes([0.1, 0.15, 0.8, 0.15])
        cbar_ax = self.figure.add_axes([0.91, 0.15, 0.025, 0.15])

        self.axes = [ax1, ax2, ax3]
        self.red_lines = []

        # pressure plot
        p_max_r = ans[2][np.where(ans[0] >= self.maxSMB)]
        t_az = np.where(ans[0] >= self.maxSMB)
        ax1.scatter(t_az, p_max_r, color='darkred', label=r'$P_{max}$')
        ax1.plot(range(len(time_bin)), ans[3], color='k', label=r'$P_{avg}$')
        ax1.legend(bbox_to_anchor=(1.1, 1.05))
        ax1.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        ax1.set_ylabel(r'P [Pa]')
        ax1.set_xlim(0, len(time_bin))
        ax1.xaxis.set_major_formatter(mtick.NullFormatter())
        ax1.grid(True)

        # corr plot
        ax2.plot(range(len(time_bin)), ans[0], color='darkgreen')
        ax2.axhline(y=self.maxSMB, color='red', linestyle='--')
        ax2.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        ax2.set_ylabel(r'$r_{max}$')
        ax2.set_xlim(0, len(time_bin))
        ax2.xaxis.set_major_formatter(mtick.NullFormatter())
        ax2.grid(True)

        # azimuth plot
        az_r = ans[1][np.where(ans[0]>=self.maxSMB)]
        ax3.scatter(t_az, az_r)
        ax3.set_ylabel(r'$azm_{r_{max}}$ [º]')
        ax3.set_xlim(0, len(time_bin))
        ax3.set_ylim(0, 360)
        ax3.yaxis.set_major_locator(mtick.FixedLocator([0,90,180,270,360]))
        ax3.xaxis.set_major_formatter(mtick.NullFormatter())
        ax3.grid(True)

        # crosscorrelation plot
        mcorr = self.out['mcorr']
        mcorr = np.flipud(mcorr.T)
        azmbin = range(mcorr.shape[0])
        halfbin_time = (time_bin[1] - time_bin[0]) / 2.0
        halfbin_azmbin = (azmbin[1] - azmbin[0]) / 2.0
        extent = (mdates.date2num(time_bin[0] - halfbin_time),
                  mdates.date2num(time_bin[-1] + halfbin_time),
                  azmbin[0] - halfbin_azmbin,
                  azmbin[-1] + halfbin_azmbin)
        cmap = plt.get_cmap('Spectral_r')
        norm = mcolor.Normalize(0, 1)
        im = ax4.imshow(mcorr, cmap=cmap, norm=norm, interpolation='gaussian', extent=extent, aspect='auto')
        ax4.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        ax4.yaxis.set_major_formatter(mtick.FuncFormatter(self.azimuthFormatter))
        ax4.axis('tight')
        ax4.xaxis_date()
        ax4.set_xlabel('Time')
        ax4.set_ylabel('azm [º]')
        ax4.set_xlim(time_bin[0],time_bin[-1])
        self.figure.colorbar(im, cax=cbar_ax, ticks=[0, 0.5, 1], orientation='vertical',
                             format='%.1f')
        
        self.nav = Navigation(self.axes)
        self.canvas.draw()


    def update_plot(self):
        self.maxSMB=float(self.smb_spinBox.value())
        self.plot()


    def __advance(self):
        self.starttime += timedelta(minutes=self.length_interval)

        if not self.air_file:
            self.compute_ccorr()
        else:
            endtime = self.starttime + timedelta(minutes=self.length_interval)
            self.out = self.air_file.get_dict(starttime=self.starttime, endtime=endtime)
            self.time_bin_max = len(self.out['time'])
            self.plot()


    def __retrocess(self):
        self.starttime -= timedelta(minutes=self.length_interval)
        if not self.air_file:
            self.compute_ccorr()
        else:
            endtime = self.starttime + timedelta(minutes=self.length_interval)
            self.out = self.air_file.get_dict(starttime=self.starttime, endtime=endtime)
            self.time_bin_max = len(self.out['time'])
            self.plot()


    def get_values(self):
        time_pad = int(self.timepadBox.value())
        freq_min = float(self.freqminSpin.value())
        freq_max = float(self.freqmaxSpin.value())
        azm = int(self.azimuthBox.value())
        time_bin = int(self.timeBox.value())

        return [(time_bin, azm), time_pad, (freq_min, freq_max)]


    def set_default(self):
        self.azimuthBox.setMaximum(self.azimuths[-1])
        self.azimuthBox.setMinimum(self.azimuths[0])
        self.timeBox.setMaximum(len(self.out['time']))
        self.set_maxazm()


    def set_maxazm(self, time_bin=None):
        if time_bin:
            self.timeBox.setValue(time_bin)

        if self.azm_cheqBox.isChecked():
            self.azimuthBox.setEnabled(False)
            with pyqtgraph.BusyCursor():
                time_pad = int(self.timeBox.value())
                ans = self.iarr.get_max_values(self.out, time_pad)
                self.azimuthBox.setValue(ans[1][0])

        else:
            self.azimuthBox.setEnabled(True)

        azm = int(self.azimuthBox.value())
        time_bin = int(self.timeBox.value())
        corr_value = self.out['mcorr'][time_bin, self.azimuths.index(azm)]
        self.correlationSpin.setValue(corr_value)

        dict_out = self.iarr.get_starttime(self.out, time_bin, azm=azm)
        self.starttime_info_bar.setText('  Time Bin: %s   |   Starttime : %s' % (time_bin,
            dict_out[list(dict_out.keys())[0]][0].strftime('%Y-%m-%d %H:%M:%S')))


    def on_click(self, event):
        for ax in self.axes:
            if event.inaxes == ax:

                # print(mdates.num2date(event.xdata))
                if self.red_lines:
                    for var in self.red_lines:
                        var.remove()

                xlim = ax.get_xlim()
                ylim = ax.get_ylim()
                hvar = ax.axhline(y = event.ydata, color='k')
                vvar = ax.axvline(x = event.xdata, color='k')
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                self.red_lines = [hvar, vvar]

                self.set_maxazm(time_bin = int(event.xdata))
                self.canvas.draw()

                if event.dblclick:
                    self.__plot_waveform()


    def __plot_waveform(self):
        self._values = self.get_values()
        print('tick --> ', self._values)

        if self.wvfm_app:
            resize = True
            width = self.wvfm_app.main.width()
            height = self.wvfm_app.main.height()
        else:
            resize = False

        self.wvfm_app = iWvfmGui(self)

        if resize:
            self.wvfm_main.resize(width, height)

        self.wvfm_main.show()


    def show_model(self):
        print('model --> ', self.out['model'])


    def show_info(self):
        print('info --> ', self.out['info'])


    def save_fig(self):
        name_def = self.out['time'][0].strftime('%Y.%m.%d_%H.%M.%S_') + str(self.length_interval)
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', '%s.png' % name_def, filter="PNG Files (*.png)")
        if name[0]:
            self.figure.savefig(name[0])


    def openDB(self):
        self.db.open_db()


    def removeID(self):
        removeID(self.db, parent=self)


class iWvfmGui(QtWidgets.QMainWindow, Ui_wvfm):
    def __init__(self, main):
        super(iWvfmGui, self).__init__(main.wvfm_main)
        self.main = main
        
        self.time_bin = main._values[0][0]
        self.azm = main._values[0][1]
        self.pad_time = main._values[1]
        self.prefilt = main._values[2]

        self._initazm()

        # get station list
        self.station_list = [main.iarr.station[x] for x in range(len(main.iarr.station))]
        self.station_spec = self.station_list[0]
        self.db = main.iarr.db_infrasound

        # crate frame
        self.init_dict()
        self.setupUi(main.wvfm_main)
        self.setWindowTitle('iWvfm Gui')
        self.show_text()
        self.initFrames()

        self.plot()


    def initFrames(self):
        self.actionClose.triggered.connect(self.close)
        self.actionPrint_image.triggered.connect(self.save_fig)

        self.figure = Figure(figsize=(17,9), constrained_layout=True)
        self.canvas = FigureCanvas(self.figure)
        #self.canvas.setParent(self.centralwidget)
        self.verticalLayout.addWidget(self.canvas)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.canvas.setSizePolicy(sizePolicy)
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.callbacks.connect('button_release_event', self.on_click)
        self.canvas.callbacks.connect('key_press_event', self.on_key)


    def _initazm(self):
        self.azm_bin = Array.get_azimuth(self.main.out, self.azm)
        self.p_max = self.main.out['p_max'][self.time_bin, self.azm_bin]
        self.p_avg = self.main.out['p_avg'][self.time_bin, self.azm_bin]
        self.corr_coef = self.main.out['mcorr'][self.time_bin, self.azm_bin]
        self.fq_dominant = None
        self.v_bars = self.main.iarr.get_starttime(
            self.main.out, self.time_bin, azm=self.azm)
        self.stream = self.main.iarr.get_waveform(
            self.main.out, self.time_bin, azm=self.azm, time_pad=self.pad_time, fq_band=self.prefilt)


    def init_dict(self):
        self.dict_out = {
            'duration':None,
            'starttime':None,
            'freq':None,
            'label':None,
            'corr': self.corr_coef,
            'azm': self.azm,
            'p_avg':None,
            'p_max':None,
            'fq_dominant':None,
            'wvfm':{}
        }


    def show_text(self):
        string = ' Pmax:  {:.3f}    Corr:  {:.2f} '.format(self.p_max, self.corr_coef)
        string += ' Pavg:  {:.3f}    Azm:  {:d} '.format(self.p_avg, self.azm)

        if not self.dict_out['freq']:
            self.dict_out['freq'] = self.fq_dominant

        if not self.dict_out['duration'] and not self.dict_out['freq']:
            string += ' Dur.: ~   Freq: ~'

        elif not self.dict_out['duration'] and self.dict_out['freq']:
            string += ' Dur.: ~    Freq: {:.1f}'.format(self.dict_out['freq'])

        elif self.dict_out['duration'] and not self.dict_out['freq']:
            string += ' Dur.: {:.1f}    Freq: ~'.format(
                self.dict_out['duration'])
        else:
            string += ' Dur.: {:.1f}    Freq: {:.1f}'.format(self.dict_out['duration'], self.dict_out['freq'])

        self.statusbar.showMessage(string)


    def on_key(self, event):
        if event.key in ('t', 'x', 'm', 'i', 'e', '1', '2', '3', '4') and (event.inaxes in self.tr_axes or event.inaxes == self.shift_axes):
            if event.key == 't':
                self.dict_out['label'] = 'TR'
            
            elif event.key == 'x':
                self.dict_out['label'] = 'XX'
            
            elif event.key == 'm':
                self.dict_out['label'] = 'MB'
            
            elif event.key == 'e':
                self.dict_out['label'] = 'EX'
            
            elif event.key == 'i':
                self.dict_out['label'] ='IMB'
            
            else:
                self.dict_out['label'] = 'L'+event.key

            self.fill_dict()
            self.save_db()

        
        if event.key == 'delete' and (event.inaxes in self.tr_axes or event.inaxes == self.shift_axes) and self._drawevents:
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            id_to_remove = []
            for id, key in self._drawevents.items():
                evnt = key['event']
                if evnt.starttime < time < evnt.starttime + timedelta(seconds=evnt.duration):
                    self.db.remove_id(id)
                    id_to_remove.append(id)
                    [txt.remove() for txt in self._drawevents[id]['txt']]
                    [span.remove() for span in self._drawevents[id]['span']]
            
            self.canvas.draw()
            for id in id_to_remove:
                del self._drawevents[id]
                    
                
        if event.key == 'a':
            azm_bin, ok = QtWidgets.QInputDialog.getInt(
                self, "Get Azimuth", "Azimuth:", self.azm, 0, 359, 1)
            if ok:
                self.azm = azm_bin
                self._initazm()
                self.init_dict()
                self.show_text()
                self.plot()
        

        if event.key == 'backspace':
            station_id_list = [sta.info.id for sta in self.station_list]
            item, ok = QtWidgets.QInputDialog.getItem(
                self, "Select Specgram", "list of stations", [sta.info.id for sta in self.station_list], 0, False)
                    
            if ok and item:
                station_spec = self.station_list[station_id_list.index(item)]
                if station_spec != self.station_spec:
                    self.station_spec = station_spec
                    self.plot_specgram()
                    self.canvas.draw()
        

        if event.key == 'S' and (event.inaxes in self.tr_axes or event.inaxes == self.shift_axes):
            text, ok = QtWidgets.QInputDialog.getText(self, "Save event", "Label:", QtWidgets.QLineEdit.Normal, "")
            if ok and text != '':
                label = text.upper()
                self.dict_out['label'] = label
                self.fill_dict()
                self.save_db()


    def save_db(self):
        event = self.db.add_event(self.dict_out)
        notify('seisvo', 'ID %s saved!' % event.id, status='info')
        print('ID saved --> ', event.id, '-- time :', event.starttime.strftime('%Y.%m.%d %H:%M:%S'))
        # self.draw_events(self.time[0], self.time[-1])
        starttime = self.dict_out['starttime']
        endtime = starttime + timedelta(seconds=self.dict_out['duration'])
        self.draw_events(self.time[0], self.time[-1])
        self.draw_events(starttime, endtime, shift=True)
        self.canvas.draw()


    def save_fig(self):
        name_def = self.out['time'][0].strftime('%Y.%m.%d_%H.%M.%S_')
        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', '%s.png' % name_def, filter="PNG Files (*.png)")
        if name[0]:
            self.figure.savefig(name[0])


    def set_axes(self):
        self.figure.clear()
        nrows = len(self.stream) + 3
        ncols = 2
        grid = {'hspace':0.3, 'width_ratios':[1, 0.01], 'hspace': 0.02, 'wspace': 0.01, 'height_ratios':[2]+[2]*len(self.stream)+[0.5,5]}
        return self.figure.subplots(nrows, ncols, gridspec_kw=grid)


    def plot_specgram(self):
        self.spec_ax.cla()
        sta_code = self.station_spec.info.code
        sta_loc = self.station_spec.info.loc
        trace = Trace2(self.stream.select(station=sta_code, location=sta_loc)[0])
        trace.specgram(window_length=5, fq_band=self.prefilt, axes=self.spec_ax, axis_bar=self.spec_bar_ax)
        self.spec_ax.xaxis.set_major_formatter(mtick.NullFormatter())
        self.spec_ax.axes.get_xaxis().set_visible(False)
        time0 = trace.stats.starttime.datetime.strftime('%Y/%m/%d %H:%M:%S')
        self.spec_ax.set_title('Start time: %s' % time0)
        self.spec_ax.set_ylabel(self.station_spec.info.id)


    def plot(self):
        self.nav_trace = None
        self.nav_shift = None

        self.dt = self.stream[0].stats.delta
        axes = self.set_axes()
        self.spec_ax = axes[0,0]
        self.spec_bar_ax = axes[0, 1]

        # self.colors = ['navy', 'black', 'darkorange', 'blueviolet']
        self.colors = get_colors('zesty')

        self.tr_axes = []
        self.time = None
        k = 1
        for sta in self.main.iarr.station:
            sta_code = sta.info.code
            sta_loc = sta.info.loc
            trace = Trace2(self.stream.select(station=sta_code, location=sta_loc)[0])

            # plot specgram
            if sta == self.station_spec:
                self.plot_specgram()

            copied_trace = trace.copy()
            copied_trace.resample(20.0) # resample the trace to 20 sps

            if not self.time:
                self.time = copied_trace.get_time()
                starttime = self.time[0]
                endtime = self.time[-1]

            axes[k,0].plot(self.time, copied_trace.data, color=self.colors[k-1])
            st1 = self.time[np.argmin(np.abs(np.array(self.time)-self.v_bars[sta.info.id][0]))]
            et1 = self.time[np.argmin(np.abs(np.array(self.time)-self.v_bars[sta.info.id][1]))]
            axes[k,0].axvspan(st1, et1, alpha=0.1, color='red')
            axes[k,0].grid(True)
            axes[k,0].set_ylabel(sta.info.id)
            axes[k,0].set_xlim(starttime, endtime)
            axes[k, 1].axes.get_xaxis().set_visible(False)
            axes[k, 1].axes.get_yaxis().set_visible(False)
            axes[k, 1].set_frame_on(False)
            self.tr_axes += [axes[k,0]]
            k += 1

        for ax in self.tr_axes[:-1]:
            ax.xaxis.set_major_formatter(mtick.NullFormatter())
        
        self.nav_trace = Navigation(self.tr_axes, imshow_axes=self.spec_ax)

        # following in white
        axes[k, 1].axes.get_xaxis().set_visible(False)
        axes[k, 1].axes.get_yaxis().set_visible(False)
        axes[k, 1].set_frame_on(False)
        axes[k, 0].axes.get_xaxis().set_visible(False)
        axes[k, 0].axes.get_yaxis().set_visible(False)
        axes[k, 0].set_frame_on(False)

        axes[k+1, 1].axes.get_xaxis().set_visible(False)
        axes[k+1, 1].axes.get_yaxis().set_visible(False)
        axes[k+1, 1].set_frame_on(False)
        axes[k+1, 0].axes.get_xaxis().set_visible(False)
        axes[k+1, 0].axes.get_yaxis().set_visible(False)
        axes[k+1, 0].set_frame_on(False)

        self.shift_axes = axes[k+1, 0]

        self._drawevents = {}
        self.draw_events(starttime, endtime)
        self.canvas.draw()


    def clear_shift_plot(self):
        self.shift_axes.cla()
        self.shift_axes.axes.get_xaxis().set_visible(False)
        self.shift_axes.axes.get_yaxis().set_visible(False)
        self.shift_axes.set_frame_on(False)


    def on_click(self, event):
        if event.inaxes == self.spec_ax:
            self.dict_out['freq'] = float(event.ydata)
            self.show_text()
        
        if event.inaxes in self.tr_axes:
            if self.nav_trace.ticks['right'][0] and self.nav_trace.ticks['left'][0]:
                times = (self.nav_trace.ticks['right'][0], self.nav_trace.ticks['left'][0])
                starttime = mdates.num2date(min(times)).replace(tzinfo=None)
                endtime = mdates.num2date(max(times)).replace(tzinfo=None)
                self.plot_shift(starttime, endtime)
                self.dict_out['duration'] = (endtime-starttime).total_seconds()
                self.dict_out['starttime'] = starttime
                self.show_text()
            
        if event.inaxes == self.shift_axes:
            if self.nav_shift.ticks['right'][0] and self.nav_shift.ticks['left'][0]:
                times = (self.nav_shift.ticks['right'][0], self.nav_shift.ticks['left'][0])
                starttime = mdates.num2date(min(times)).replace(tzinfo=None)
                endtime = mdates.num2date(max(times)).replace(tzinfo=None)
                self.dict_out['duration'] = (endtime-starttime).total_seconds()
                self.dict_out['starttime'] = starttime
                self.show_text()


    def plot_shift(self, starttime, endtime):
        self.shift_axes.axes.get_xaxis().set_visible(True)
        self.shift_axes.axes.get_yaxis().set_visible(True)
        self.shift_axes.set_frame_on(True)
        self.shift_axes.cla()

        self.st_shift = Stream2()
        for sta in self.main.iarr.station:
            st = self.stream.select(station=sta.info.code, location=sta.info.loc)
            sec = float((int(self.main.iarr.dt_times[sta.info.id][self.azm_bin]) - 1) * self.dt)
            start_time = UTCDateTime(starttime + timedelta(seconds=sec))
            end_time = UTCDateTime(endtime + timedelta(seconds=sec))
            stream = Stream2(st.slice(start_time, end_time))
            self.st_shift += stream

            if self.main.iarr.central_sta == sta.info.code:
                self.shift_plot_time = stream[0].get_time()

        min_npts = min([tr.stats.npts for tr in stream])
        for n, trace in enumerate(self.st_shift):
            self.shift_axes.plot(self.shift_plot_time[:min_npts], trace.data[:min_npts], color=self.colors[n])

        self.shift_axes.set_xlabel('Seconds\n%s' % starttime.strftime('%H:%M:%S'))
        self.shift_axes.set_xlim(self.shift_plot_time[0], self.shift_plot_time[-1])
        self.shift_axes.xaxis.set_major_formatter(mdates.DateFormatter('%S.%f'))
        self.shift_axes.xaxis.set_major_locator(mtick.MaxNLocator(nbins=10, min_n_ticks=5))
        self.shift_axes.xaxis.set_minor_locator(mtick.AutoMinorLocator(4))
        self.nav_shift = Navigation(self.shift_axes)

        saved_txt = 't(TR) x(XX) m(MB) e(EX) i(IMB)'
        self.shift_axes.annotate(saved_txt, xy=(0.01,0.01), xycoords='figure fraction', fontsize=8, fontfamily='monospace')


        self.draw_events(self.shift_plot_time[0], self.shift_plot_time[-1], shift=True)
        self.canvas.draw()


    def compute_fq(self, data, sps):
        psd, fq = power_density_spectrum(data, sps, fq_band=self.prefilt)
        return fq[np.argmax(psd)]


    def fill_dict(self):
        starttime = self.dict_out['starttime']
        endtime = starttime + timedelta(seconds=self.dict_out['duration'])
        p_avg = []
        p_max = []
        fq_dominant = []
        for trace in self.st_shift:
            trace_id = '%s.%s.%s' % (trace.stats.network, trace.stats.station, trace.stats.location)
            self.dict_out['wvfm'][trace_id] = {}

            sec = float((int(self.main.iarr.dt_times[trace_id][self.azm_bin]) - 1) * self.dt)
            self.dict_out['wvfm'][trace_id]['delta'] = sec

            data = trace.get_data(starttime=starttime, endtime=endtime, detrend=True)

            data_abs = np.abs(data)
            pmax = data_abs.max()
            self.dict_out['wvfm'][trace_id]['p_max'] = pmax
            p_max += [pmax]

            pavg = data_abs.mean()
            p_avg += [pavg]
            self.dict_out['wvfm'][trace_id]['p_avg'] = pavg

            fqd = self.compute_fq(data, trace.stats.sampling_rate)
            fq_dominant += [fqd]
            self.dict_out['wvfm'][trace_id]['fq_dominant'] = fqd

        self.dict_out['fq_dominant'] = np.mean(fq_dominant)
        self.dict_out['p_avg'] = np.mean(p_avg)
        self.dict_out['p_max'] = np.mean(p_max)


    def draw_events(self, starttime, endtime, shift=False):
        events_list = self.db.get_event_list(starttime=starttime, endtime=endtime)

        for event in events_list:
            if event.id not in list(self._drawevents.keys()):
                exist_in_main = False
                self._drawevents[event.id] = {
                    'span': [],
                    'txt':[],
                    'event':event
                }
            else:
                exist_in_main = True
            st = event.starttime
            et = event.starttime + timedelta(seconds=event.duration)

            cond_in = st >= starttime and et <= endtime # start and end in
            cond_st_in = endtime > st > starttime and et > endtime # start in, end out
            cond_et_in = starttime > st and starttime < et < endtime # start out, end in
            cond_full_in = st < starttime and endtime < et # start and end out

            if cond_in or cond_st_in or cond_et_in or cond_full_in:

                if cond_in:
                    st1 = self.time[np.argmin(np.abs(np.array(self.time)-st))]
                    et1 = self.time[np.argmin(np.abs(np.array(self.time)-et))]

                if cond_st_in:
                    st1 = self.time[np.argmin(np.abs(np.array(self.time)-st))]
                    et1 = endtime

                if cond_et_in:
                    st1 = starttime
                    et1 = self.time[np.argmin(np.abs(np.array(self.time)-et))]

                if cond_full_in:
                    st1 = starttime
                    et1 = endtime

                if shift:
                    sp = self.shift_axes.axvspan(st1, et1, color='teal', alpha=0.2)
                    text_pos = mdates.date2num(st1 + timedelta(seconds=0.3))
                    tx = self.shift_axes.annotate(
                        'ID:%s(%s)' % (event.id, event.label), 
                        (text_pos, self.shift_axes.get_ylim()[1]),
                        color='teal', fontsize=9)
                    self._drawevents[event.id]['span'].append(sp)
                    self._drawevents[event.id]['txt'].append(tx)
                else:
                    if not exist_in_main:
                        for n, ax in enumerate(self.tr_axes):
                            sp = ax.axvspan(st1, et1, color='teal', alpha=0.2)
                            self._drawevents[event.id]['span'].append(sp)
                            if n == 0:
                                text_pos = mdates.date2num(st1 + timedelta(seconds=0.3))
                                tx = ax.annotate(
                                    'ID:%s(%s)' % (event.id, event.label), 
                                    (text_pos, ax.get_ylim()[1]),
                                    color='teal', fontsize=9)
                                self._drawevents[event.id]['txt'].append(tx)
