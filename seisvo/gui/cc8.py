#!/usr/bin/env python3
# coding=utf-8


import pyqtgraph
import datetime as dt

# matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.text import Annotation
from .utils import DateDialog, Navigate, notify
from ..database import CCE

# utils

def save_fig(fig, filename, parent=None):
    files_types = 'Image (*.png);; Document (*.pdf)'
    fileName, _ = QtWidgets.QFileDialog.getSaveFileName(parent, 'Save File', filename, files_types)
    if fileName != '':
        fig.savefig(fileName, bbox_inches='tight')
        print(f" [info] {fileName} created")


class CC8ConfigDialog(QtWidgets.QDialog):
    def __init__(self, ValsDict, minmax_times, parent=None):
        super(CC8ConfigDialog, self).__init__(parent)
        self.setWindowTitle("CC8 GUI - Config")

        # fixed size
        self.resize(300, 400)
        self.setMinimumSize(QtCore.QSize(300, 400))
        self.setMaximumSize(QtCore.QSize(300, 400))

        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self)
        self.verticalLayout   = QtWidgets.QVBoxLayout()
        self.gridLayout       = QtWidgets.QGridLayout()
        self.verticalLayout.addLayout(self.gridLayout)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.setLabels()
        self.setSpaces()
        self.minmax_times = minmax_times
        self.setEmbedings(ValsDict)

        self.buttonBox = QtWidgets.QDialogButtonBox(self)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
        self.verticalLayout_2.addWidget(self.buttonBox)

        self.buttonBox.accepted.connect(self.accept) # type: ignore
        self.buttonBox.rejected.connect(self.reject) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(self)
    

    def setLabels(self):
        self.label_1 = QtWidgets.QLabel(self)
        self.label_1.setText("maac_th")
        self.gridLayout.addWidget(self.label_1, 5, 1, 1, 1)
        
        self.label_2 = QtWidgets.QLabel(self)
        self.label_2.setText("error_th")
        self.gridLayout.addWidget(self.label_2, 6, 1, 1, 1)
        
        self.label_3 = QtWidgets.QLabel(self)
        self.label_3.setText("rms_th")
        self.gridLayout.addWidget(self.label_3, 7, 1, 1, 1)
        
        self.label_4 = QtWidgets.QLabel(self)
        self.label_4.setText("slow_max")
        self.gridLayout.addWidget(self.label_4, 10, 1, 1, 1)
        
        self.label_5 = QtWidgets.QLabel(self)
        self.label_5.setText("slow_int")
        self.gridLayout.addWidget(self.label_5, 11, 1, 1, 1)
        
        self.label_6 = QtWidgets.QLabel(self)
        self.label_6.setText("fq_idx")
        self.gridLayout.addWidget(self.label_6, 1, 1, 1, 1)
        
        self.label_7 = QtWidgets.QLabel(self)
        self.label_7.setText("For recomputed slowness map")
        self.label_7.setWordWrap(True)
        self.gridLayout.addWidget(self.label_7, 9, 0, 1, 3)

        self.label_8 = QtWidgets.QLabel(self)
        self.label_8.setText("time")
        self.gridLayout.addWidget(self.label_8, 2, 1, 1, 1)
        
        self.label_9 = QtWidgets.QLabel(self)
        self.label_9.setText("interval")
        self.gridLayout.addWidget(self.label_9, 3, 1, 1, 1)


    def setSpaces(self):

        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        self.gridLayout.addItem(spacerItem, 8, 1, 1, 1)
        
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem1, 12, 1, 1, 1)

        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem2, 0, 1, 1, 1)
        
        spacerItem3 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        self.gridLayout.addItem(spacerItem3, 4, 1, 1, 1)

        spacerItem4 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem4, 7, 0, 1, 1)
        
        spacerItem5 = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem5, 7, 3, 1, 1)

        spacerItem6 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem6)


    def setEmbedings(self, vdict):
        self.qtimebox = QtWidgets.QDateTimeEdit(self)
        self.qtimebox.setCalendarPopup(True)
        self.qtimebox.setDateTime(QtCore.QDateTime(vdict["time"]))
        self.qtimebox.setMinimumDateTime(QtCore.QDateTime(self.minmax_times[0]))
        self.qtimebox.setMaximumDateTime(QtCore.QDateTime(self.minmax_times[1]))
        self.gridLayout.addWidget(self.qtimebox, 2, 2, 1, 2)

        self.interval_dsb = QtWidgets.QDoubleSpinBox(self)
        self.interval_dsb.setDecimals(0)
        self.interval_dsb.setMinimum(10.0)
        self.interval_dsb.setMaximum(120.0)
        self.interval_dsb.setProperty("value", vdict["interval"])
        self.gridLayout.addWidget(self.interval_dsb, 3, 2, 1, 1)

        self.fq_cbox = QtWidgets.QComboBox(self)
        self.fq_cbox.addItems(vdict["fq_idx_list"])
        self.fq_cbox.setCurrentIndex(vdict["fq_idx_list"].index(vdict["fq_idx"]))
        self.gridLayout.addWidget(self.fq_cbox, 1, 2, 1, 1)

        self.maac_dsb = QtWidgets.QDoubleSpinBox(self)
        self.maac_dsb.setDecimals(2)
        self.maac_dsb.setMaximum(1.0)
        self.maac_dsb.setSingleStep(0.05)
        self.maac_dsb.setProperty("value", vdict["maac"])
        self.gridLayout.addWidget(self.maac_dsb, 5, 2, 1, 1)

        self.error_dsb = QtWidgets.QDoubleSpinBox(self)
        self.error_dsb.setMaximum(2.0)
        self.error_dsb.setSingleStep(0.05)
        self.error_dsb.setProperty("value", vdict["error"])
        self.gridLayout.addWidget(self.error_dsb, 6, 2, 1, 1)

        self.rms_dsb = QtWidgets.QDoubleSpinBox(self)
        self.rms_dsb.setDecimals(0)
        self.rms_dsb.setProperty("value", vdict["rms"])
        self.gridLayout.addWidget(self.rms_dsb, 7, 2, 1, 1)

        self.slomax_dsb = QtWidgets.QDoubleSpinBox(self)
        self.slomax_dsb.setDecimals(1)
        self.slomax_dsb.setMinimum(0.5)
        self.slomax_dsb.setMaximum(5.0)
        self.slomax_dsb.setSingleStep(0.5)
        self.slomax_dsb.setProperty("value", vdict["slomax"])
        self.gridLayout.addWidget(self.slomax_dsb, 10, 2, 1, 1)
        
        self.sloint_dsb = QtWidgets.QDoubleSpinBox(self)
        self.sloint_dsb.setMinimum(0.01)
        self.sloint_dsb.setMaximum(1.0)
        self.sloint_dsb.setSingleStep(0.01)
        self.sloint_dsb.setProperty("value", vdict["sloint"])
        self.gridLayout.addWidget(self.sloint_dsb, 11, 2, 1, 1)


    @staticmethod
    def getVals(ValsDict, parent=None):
        dialog = CC8ConfigDialog(ValsDict, parent)
        result = dialog.exec_()
        dvals = {
          "time":dialog.qtimebox.dateTime().toPyDateTime(),
          "interval":dialog.interval_dsb.value(),
          "fq_idx":dialog.fq_cbox.currentText(),
          "maac":dialog.maac_dsb.value(),
          "rms":dialog.rms_dsb.value(),
          "error":dialog.error_dsb.value(),
          "slomax":dialog.slomax_dsb.value(),
          "sloint":dialog.sloint_dsb.value(),
        }
        return (dvals, result==QtWidgets.QDialog.Accepted)

# SLOWMAP Widget
 
class CC8nidxWidget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS y controla los botones
    """

    def __init__(self, cc8out, nidx_list, parent=None):
        QtWidgets.QWidget.__init__(self)
        self.cc8out    = cc8out
        self.cc8canvas = parent  # main canvas
        self.nidx_list = nidx_list
        self.max_nidx  = len(self.nidx_list)
        self.nidx_idx  = None
        self.layout    = QtWidgets.QVBoxLayout()
        self.canvas    = CC8nidxCanvas(self)
        self.layout.addWidget(self.canvas)

        self.setLayout(self.layout)
        self.setCursor(QtCore.Qt.CrossCursor)
    

    def plot(self, nidx=None):
        if not nidx:
            nidx = self.nidx_list[self.nidx_idx][0]
        else:
            self.nidx_idx, = np.where(self.nidx_list == nidx)

        self.canvas.plot(int(nidx))
        self.canvas.setFocus()
    

    def update(self, cc8out, nidx_list):
        self.cc8out    = cc8out
        self.nidx_list = nidx_list
        self.max_nidx  = len(self.nidx_list)
        self.nidx_idx  = None
        self.close()


    def on_key(self, key):
        if key == QtCore.Qt.Key_Right:
            if self.nidx_idx + 1 ==  self.max_nidx:
                notify("CC8", "The bound of the file was reached!")
            else:
                self.nidx_idx += 1
                self.plot()
            
        if key == QtCore.Qt.Key_Left:
            if self.nidx_idx - 1 < 0:
                notify("CC8", "The bound of the file was reached!")
            else:
                self.nidx_idx -= 1
                self.plot()


    def save_fig(self, filename):
        save_fig(self.canvas.fig, filename, parent=self)


class CC8nidxCanvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.filter = True # by default, active the filter
        self.fig = Figure(figsize=(9,12))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)


    def on_key(self, event):
        if event.key == "h":
            print("   key   --   info   ")
            print(" --------  ----------")
            print("   ->    --  advance nidx")
            print("   <-    --  retroces nidx")
            print("   f1    --  save fig")
            print("    p    --  show particle motion of window")
            print("    f    --  remove/add filter")
            print("    a    --  save nidx into database")
            print("    l    --  add label nidx in database")
            print("   del   --  delete nidx from database")


        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)


        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)
        

        if event.key == 'f1':
            filename = f"Slowmap_1"
            self.parent.save_fig(filename)


        if event.key == 'p':
            # plot particle motion
            ans, _ = self.parent.cc8out.cc8stats.get_array()
            if ans:
                sta01   = ans.get_sta("01")
                half_w  = dt.timedelta(seconds=float(self.parent.cc8out.cc8stats.window/2))
                start   = self.parent.cc8canvas.hover_time - half_w
                end     = start + half_w + half_w
                fig     = sta01.particle_motion(start, end, baz=self.parent.cc8canvas.baz0)


        if event.key == 'f':
            if self.filter:
                self.filter = False
            else:
                self.filter = True
            self.plot(self.nidx)
        

        if event.key in ('a' 'l') and self.parent.cc8canvas.parent.is_db:
            
            if event.key == 'l':
                if self.row:
                    #relabel
                    rowid = self.row.id
                    label , ok = QtWidgets.QInputDialog.getText(self,\
                            f'Evento {rowid}', 'Define Etiqueta:', text=f"{self.row.label}")
                    if ok:
                        self.parent.cc8canvas.parent.db.relabel_event(rowid, label)
                        notify("CC8", f" Nueva etiqueta en ID {rowid}")
                        self.axes["ax2"].set_title(f" ID [{rowid}] LABEL [{label}]")
                        self.draw()
                
                else:
                    label , ok = QtWidgets.QInputDialog.getText(self,\
                            f'Nuevo Evento', 'Define Etiqueta:', text="")
                    if ok:
                        rowid    = self.parent.cc8canvas.add_event(nidx=self.nidx, label=label)
                        self.row = self.parent.cc8canvas.get_row(self.nidx)
                        self.axes["ax2"].set_title(f" ID [{rowid}] LABEL [{label}]")
                        self.draw()

            if event.key == 'a':
                if self.row:
                    # nidx is saved in database!!
                    return
                else:
                    rowid    = self.parent.cc8canvas.add_event(nidx=self.nidx)
                    self.row = self.parent.cc8canvas.get_row(self.nidx)
                    self.axes["ax2"].set_title(f" ID [{rowid}] LABEL [None]")
                    self.draw()
        

        if event.key == "delete" and self.row:
            self.parent.cc8canvas.parent.rm_nidx(self.nidx, self.row.id)
            self.axes["ax2"].set_title(f"")
            self.draw()
        

    def print_ticks(self):
        print("\n :: Ticks info :: ")
        for position, time in self.ticks.items():
            if time:
                print(f"  {position} >> {time} ")
        
        if self.ticks['right'] and self.ticks['left']:
            dist = abs((self.ticks['left'] - self.ticks['right']))
            print(f" distance |R-L|  >>  {dist:.2f} s  [{1/dist:.2f} Hz]")
    

    def on_click(self, event):
        if event.inaxes:
            if event.inaxes in self.nav.ax:
                
                # t = mdates.num2date(float(event.xdata))
                # t = t.replace(tzinfo=None)

                if event.button == 1:
                    self.ticks['left'] = event.xdata
                
                if event.button == 3:
                    self.ticks['right'] = event.xdata
                
                self.print_ticks()


    def plot(self, nidx):
        # CHECK IF NIDX is in DB and take ROW
        self.row = self.parent.cc8canvas.get_row(nidx)

        with pyqtgraph.BusyCursor():
            self.fig.clf()
            self.ticks = dict(right=None, left=None)
            gs = gridspec.GridSpec(4, 5, figure=self.fig, wspace=0.1, height_ratios=[1,0.1,0.75,0.75], width_ratios=[1,0.02,0.2,1,0.02])

            ax1  = self.fig.add_subplot(gs[0, 0])
            ax1b = self.fig.add_subplot(gs[0, 1])
            ax2  = self.fig.add_subplot(gs[0, 3])
            # ax2b = self.fig.add_subplot(gs[0, 4])
            ax3  = self.fig.add_subplot(gs[2, :])
            ax4  = self.fig.add_subplot(gs[3, :])
            self.axes = {
                "ax1":ax1,
                "ax2":ax2,
                "ax3":ax3,
                "ax4":ax4
            }

            if self.parent.cc8out.cc8stats.slowmap:
                self.parent.cc8out.plot_smap(nidx, axis=ax1, bar_axis=ax1b, fig=self.fig)

            try:
                # self.parent.cc8canvas.ccout.plot_detailed_smap(nidx, fq_idx=self.parent.cc8canvas.parent.fq_idx, show_title=False, slomax=self.parent.cc8canvas.parent.r_slomax, sloint=self.parent.cc8canvas.parent.r_sloint, axis=ax2, bar_axis=ax2b, fig=self.fig)
                # ax2.xaxis.set_major_formatter(mtick.NullFormatter())
                # ax2.yaxis.set_major_formatter(mtick.NullFormatter())
                # ax2.set_xlabel("")
                # ax2.set_ylabel("")

                if self.filter:
                    fqband = [1., 5.]
                else:
                    fqband = [] 
                
                self.parent.cc8out.plot_beamform(ntime=nidx, fq_band=fqband, off_sec=2, axes=[ax3, ax4], fig=self.fig, show_title=False)
                self.nav = Navigate([ax3, ax4], self, color='red', linewidth=0.5, alpha=0.5)
            except Exception as e:
                print(" [warn] .... ")
                print(e)

            if self.row:
                ax2.set_title(f" ID [{self.row.id}] LABEL [{self.row.label}]")

            self.draw()
            self.parent.show()
        
        if self.parent.cc8canvas:
            self.parent.cc8canvas.show_green(nidx)

        self.nidx = nidx

# CC8 Widget

class CC8Widget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS y controla los botones
    """

    def __init__(self, cc8, starttime, interval, fq_idx, db, **kwargs):

        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()
        self.cc8       = cc8
        self.fq_idx    = fq_idx
        
        # load database
        if isinstance(db, str):
            self.db    = CCE(db)
            self.is_db = True
        else:
            self.db    = None
            self.is_db = False
        
        self.olap_pct  = kwargs.get("olap", 0.1)
        self.maac_th   = kwargs.get("maac_th", 0.6)
        self.max_err   = kwargs.get("max_err", 0.7)
        self.rms_th    = kwargs.get("rms_th", 0.0)
        rms_lim        = kwargs.get("rms_lim",[])
        self.r_slomax  = 0.5
        self.r_sloint  = 0.01
        self.starttime = starttime
        self.interval  = dt.timedelta(minutes=interval)
        self.olap      = dt.timedelta(minutes=interval*self.olap_pct)

        # before plot, compute the bounds of the RMS
        if not rms_lim:
            rms_data = self.cc8.get(attr="rms", fq_idx=self.fq_idx).get_data("rms")
            rms_lim = [np.nanmin(rms_data), np.nanmax(rms_data)]

        self.rms_min = rms_lim[0]
        self.rms_max = rms_lim[1]

        # add canvas
        self.canvas = CC8Canvas(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        self.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.setFocus()

    
    def on_key(self, key):
        if key == QtCore.Qt.Key_Right:
            # print(self.starttime + self.interval, self.cc8.time_[-1])

            if self.starttime + self.interval > self.cc8.time_[-1]:
                notify("CC8", "The bound of the file was reached!")
                endtime = self.cc8.time_[-1]
                new_starttime = endtime - self.interval - self.olap
            else:
                endtime = None
                new_starttime = self.starttime + self.interval - self.olap
            
            self.plot(new_starttime, endtime=endtime)

        if key == QtCore.Qt.Key_Left:
            if self.starttime - self.interval < self.cc8.time_[0]:
                notify("CC8", "The bound of the file was reached!")
                new_starttime = self.cc8.time_[0]
            else:
                new_starttime = self.starttime - self.interval + self.olap
            self.plot(new_starttime)


    def plot(self, starttime, endtime=None):
        self.starttime = starttime
        self.canvas.plot(endtime=endtime) 


    def update(self, new_vals, plot):
        self.r_slomax  = new_vals["slomax"]
        self.r_sloint  = new_vals["sloint"]
        self.maac_th   = new_vals["maac"]
        self.max_err   = new_vals["error"]
        self.rms_th    = new_vals["rms"]
        self.interval  = dt.timedelta(minutes=new_vals["interval"])
        self.olap      = dt.timedelta(minutes=new_vals["interval"]*self.olap_pct)

        if self.fq_idx != new_vals["fq_idx"]:
            self.fq_idx = new_vals["fq_idx"]
            rms_data = self.cc8.get(attr="rms", fq_idx=self.fq_idx).get_data("rms")
            rms_lim = [np.nanmin(rms_data), np.nanmax(rms_data)]
            self.rms_min = rms_lim[0]
            self.rms_max = rms_lim[1]

        if plot:
            self.plot(new_vals["time"])


    def save_fig(self, filename):
        save_fig(self.canvas.fig, filename, parent=self)


    def save_nidx(self, nidx, edict, label=None):
        if self.is_db:
            rowid = self.db._add_row(edict)
            self.canvas.colorlist[nidx] = "red"
            self.canvas.nidx_db[rowid] = nidx

            if label:
                self.db.relabel_event(rowid, label)

            notify("CC8", f"new ID [{rowid}] in database!")
            return rowid


    def rm_nidx(self, nidx, eid):
        self.db.remove_event(eid)
        self.canvas.colorlist[nidx] = "blue"
        del self.canvas.nidx_db[eid]
        notify("CC8", f" ID [{eid}] removed from database!")


class CC8Canvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.WidgetNidex = None
        self.fig = Figure(figsize=(12,9))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
        self.callbacks.connect("motion_notify_event", self.on_hover)
        self.plot()


    def plot(self, endtime=None):
        self.fig.clf()
        self.ticks = dict(right=None, left=None)
        self.baz0  = None
        self.slow0 = None

        if not endtime:
            self.endtime = self.parent.starttime + self.parent.interval
        else:
            self.endtime = endtime

        with pyqtgraph.BusyCursor():
            # get out object
            self.ccout = self.parent.cc8.get(starttime=self.parent.starttime, endtime=self.endtime, slowmap=True, fq_idx=self.parent.fq_idx)
            
            # set nidx widget
            self.nidx_list = self.ccout.get_nidx(return_full=False, max_err=self.parent.max_err, maac_th=self.parent.maac_th, rms_th=self.parent.rms_th)
            
            if not self.WidgetNidex:
                # open a widget
                self.WidgetNidex = CC8nidxWidget(self.ccout, self.nidx_list, parent=self)
            else:
                self.WidgetNidex.update(self.ccout, self.nidx_list)

            self.time = np.array(self.ccout._dout["dtime"])
            self.fig_dict = self.ccout.plot(max_err=self.parent.max_err, maac_th=self.parent.maac_th, rms_th=self.parent.rms_th, fig=self.fig, x_time=True, rms_lim=[self.parent.rms_min, self.parent.rms_max], return_fdict=True)

            # add horizontal bar in maac_th and rms
            self.fig_dict["maac"]["axis"].axhline(self.parent.maac_th, color="r", ls="--", alpha=0.7, zorder=1)
            self.fig_dict["rms"]["axis"].axhline(self.parent.rms_th, color="r", ls="--", alpha=0.7, zorder=1)
            
            # add navigation
            self.axes_ = [ax["axis"] for _, ax in self.fig_dict.items()]
            self.nav = Navigate(self.axes_, self, color='red', active_cursor=False, linewidth=0.5, alpha=0.5)

            # search in database
            self.colorlist = ["blue"]*len(self.ccout._dout["dtime"])
            self.set_nidxdb()

            for _, fdict in self.fig_dict.items():
                fdict["sc"].set_facecolor(self.colorlist)

            self.draw()
    

    def set_nidxdb(self):
        if self.parent.is_db:
            saved_points = self.parent.db.get_episode_list(time_interval=(self.parent.starttime, self.endtime))
        else:
            saved_points = []
        
        self.nidx_db = {}

        if saved_points:
            for sp in saved_points:
                row = self.parent.db[sp]
                if row.time in self.ccout._dout["dtime"]:
                    sp_nidx = list(np.where(self.time == row.time)[0])[0]
                    self.nidx_db[sp] = sp_nidx
                    self.colorlist[sp_nidx] = "red"
        
        return
    

    def get_row(self, nidx):
        if self.parent.is_db:
            self.set_nidxdb()
            row = None

            for rowid, nidx_db in self.nidx_db.items():
                if nidx_db == nidx:
                    row = self.parent.db[rowid]
            
            return row

        return None


    def add_event(self, nidx, label=None):
        if not self.parent.is_db:
            notify("CC8", f" no database loaded! ")
            return
        
        # add event
        time = self.time[nidx]
        maac = self.fig_dict["maac"]["sc"].get_offsets()[nidx][1]
        slow = self.fig_dict["slow"]["sc"].get_offsets()[nidx][1]
        baz = self.fig_dict["baz"]["sc"].get_offsets()[nidx][1]
        rms  = self.fig_dict["rms"]["sc"].get_offsets()[nidx][1]
        
        event_dict = {
            "network":self.parent.cc8.stats.id.split(".")[0],
            "station":self.parent.cc8.stats.id.split(".")[1],
            "time":time,
            "slow":slow,
            "baz":baz,
            "maac":maac,
            "rms":rms,
            "cc8_file":self.parent.cc8.stats.file,
            "fqidx":self.parent.fq_idx,
        }

        return self.parent.save_nidx(nidx, event_dict, label=label)


    def get_axis_attr(self, axis):
        if axis in self.axes_:
            for attr, fdict in self.fig_dict.items():
                if fdict["axis"] == axis:
                    return attr
        return None


    def print_ticks(self):
        print("\n :: Ticks info :: ")
        for position, time in self.ticks.items():
            if time:
                print(f"  {position} >> {time} ")
        
        if self.ticks['right'] and self.ticks['left']:
            dist = abs((self.ticks['left'] - self.ticks['right']).total_seconds())
            print(f" distance |R-L|  >>  {dist}  [sec]")


    def on_click(self, event):
        if event.inaxes:
            if event.inaxes in self.axes_:
                
                t = mdates.num2date(float(event.xdata))
                t = t.replace(tzinfo=None)

                if event.button == 1:
                    self.ticks['left'] = t
                
                if event.button == 3:
                    self.ticks['right'] = t
                
                self.print_ticks()
    

    def show_green(self, nidx):
        colorlist       = self.colorlist.copy()
        colorlist[nidx] = "green"
        
        for _, fdict in self.fig_dict.items():
            fdict["sc"].set_facecolor(colorlist)
        
        self.draw()


    def on_hover(self, event):
        attr = self.get_axis_attr(event.inaxes)

        if attr:
            sc = self.fig_dict[attr]["sc"]
            cont, ind = sc.contains(event)
            if cont:
                x, _  = sc.get_offsets()[ind["ind"][0]]
                x     = mdates.num2date(float(x)).replace(tzinfo=None)
                nidx  = np.argmin(np.abs(x - self.time))
                self.hover_time = x
                self.hover_nidx = nidx
                self.baz0  = self.fig_dict["baz"]["sc"].get_offsets()[nidx][1]
                self.slow0 = self.fig_dict["slow"]["sc"].get_offsets()[nidx][1]
                self.show_green(nidx)


    def on_key(self, event):
        if event.key == "h":
            print("   key   --   info   ")
            print(" --------  ----------")
            print("   ->    --  advance interval")
            print("   <-    --  retroces interval")
            print("   f1    --  save fig")
            print("    h    --  show key info")
            print("    s    --  show slowmap of selected point")
            print("    t    --  show SLOW and BAZ of selected point")
            print("    c    --  open config widget")
            print("    p    --  show proble slowmap")
            print("    w    --  show beam form between ticks for\n\
                         SLOW and BAZ of selected point")
            print("   tab   --  goto saved ID")


        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)


        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)


        if event.key == "s" and self.hover_nidx:
            self.WidgetNidex.plot(nidx=self.hover_nidx)
        

        if event.key == "f1":
            filename = "Detecciones_1"
            self.parent.save_fig(filename)


        if event.key == "t" and self.hover_nidx:
            print(f"\n >> {self.hover_time.strftime('%Y %b %d %H:%M:%S.%f')} :: BAZ {self.baz0:5.1f} SLOW {self.slow0:4.2f}")


        if event.key == "c":
            cvals = {
                "time":self.parent.starttime,
                "interval":int(self.parent.interval.total_seconds()/60),
                "fq_idx_list":self.parent.cc8.stats.fqidx,
                "fq_idx":self.parent.fq_idx,
                "maac":self.parent.maac_th,
                "rms":self.parent.rms_th,
                "error":self.parent.max_err,
                "slomax":self.parent.r_slomax,
                "sloint":self.parent.r_sloint,
            }
            minmax_times = (self.parent.cc8.time_[0],self.parent.cc8.time_[-1])
            newvals, ans = CC8ConfigDialog.getVals(cvals, minmax_times)
            
            if ans:
                c1 = self.parent.max_err != newvals["error"] 
                c2 = self.parent.rms_th != newvals["rms"]
                c3 = self.parent.maac_th != newvals["maac"]
                c4 = self.parent.fq_idx != newvals["fq_idx"]
                c5 = self.parent.starttime != newvals["time"]
                c6 = int(self.parent.interval.total_seconds()/60) != newvals["interval"]
                
            
                if any([c1,c2,c3,c4,c5,c6]):
                    replot = True

                else:
                    replot = False

                self.parent.update(newvals, replot)


        if event.key == "w" and self.baz0 and self.slow0: 
            # particle motio
            if self.ticks['right'] and self.ticks['left']:
                starttime = min([self.ticks['right'],self.ticks['left']])
                endtime   = max([self.ticks['right'],self.ticks['left']])
                with pyqtgraph.BusyCursor():
                    fig = self.ccout.get_beamform(starttime, endtime, self.slow0, self.baz0)
                fig.close()
        

        if event.key == "p":
            # probability slowmap
            with pyqtgraph.BusyCursor():
                fig = self.ccout.prob_slowmap(fq_idx=self.parent.fq_idx, max_err=self.parent.max_err, maac_th=self.parent.maac_th)
            plt.close(fig)


        if event.key == "tab":
            id_list = self.parent.db.get_episode_list()
            print(id_list)
            # combo box get ID and change starttime/endtime
