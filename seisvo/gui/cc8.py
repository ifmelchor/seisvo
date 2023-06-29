#!/usr/bin/env python3
# coding=utf-8


import pyqtgraph
import datetime as dt

# matplotlib
import numpy as np
import matplotlib.dates as mdates
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.text import Annotation
from .utils import Navigate, notify

 
class CC8nidxWidget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS y controla los botones
    """

    def __init__(self, cc8out, nidx):
        QtWidgets.QWidget.__init__(self)
        self.cc8out = cc8out
        self.nidx   = nidx
        self.layout = QtWidgets.QVBoxLayout()
        # add canvas
        self.canvas = CC8nidxCanvas(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        self.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.setFocus()


class CC8nidxCanvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.fig = Figure(figsize=(9,9))
        FigureCanvas.__init__(self, self.fig)
        self.plot()


    def plot(self):
        with pyqtgraph.BusyCursor():
            self.fig.clf()
            axes = self.fig.subplots(2,2, gridspec_kw={"height_ratios":[1,0.75], "width_ratios":[1,0.02]})
            axes[1,1].axis("off")
            self.parent.cc8out.plot_smap(self.parent.nidx, axis=axes[0,0], bar_axis=axes[0,1], fig=self.fig)
            self.parent.cc8out.plot_wvfm(self.parent.nidx, off_sec=5, axis=axes[1,0], fig=self.fig, show_title=False)
            self.nav = Navigate([axes[1,0]], self, color='red', linewidth=0.5, alpha=0.5)
            self.draw()


class CC8Widget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS y controla los botones
    """

    def __init__(self, cc8, starttime, interval, fq_idx, slow_idx,\
        olap=0.1, maac_th=0.6, baz_int=[]):

        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()
        
        self.cc8       = cc8
        self.fq_idx    = fq_idx
        self.slow_idx  = slow_idx
        self.maac_th   = maac_th
        self.baz_int   = baz_int
        self.starttime = starttime
        self.interval  = dt.timedelta(minutes=interval)
        self.olap      = dt.timedelta(minutes=interval*olap)

        # add canvas
        self.canvas = CC8Canvas(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        self.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.setFocus()

        # init nidex widget

    

    def on_key(self, key):
        if key == QtCore.Qt.Key_Right:
            if self.starttime + self.interval > self.cc8.stats.endtime:
                notify("CC8", "The bound of the file was reached!")
                endtime = self.cc8.stats.endtime
                self.starttime = endtime - self.interval
            else:
                endtime = None
                self.starttime = self.starttime + self.interval
            self.canvas.plot(endtime=endtime)

        if key == QtCore.Qt.Key_Left:
            if self.starttime - self.interval + self.olap < self.cc8.starttime:
                notify("CC8", "The bound of the file was reached!")
                self.starttime = self.cc8.starttime
            else:
                self.starttime = self.starttime - self.interval + self.olap
            self.canvas.plot()


class CC8Canvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.fig = Figure(figsize=(12,9), subplotpars=SubplotParams(left=0.08,\
            right=0.92, wspace=0.1, top=0.95, bottom=0.05))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
        self.callbacks.connect("motion_notify_event", self.on_hover)
        self.plot()

        # init nidex widget
        self.WidgetNidex = None


    def plot(self, endtime=None):
        self.fig.clf()
        self.ticks = dict(right=None, left=None)

        if not endtime:
            self.endtime = self.parent.starttime + self.parent.interval
        else:
            self.endtime = endtime
        
        with pyqtgraph.BusyCursor():

            self.ccout = self.parent.cc8.get(starttime=self.parent.starttime,\
                endtime=self.endtime, slowmap=True, fq_idx=self.parent.fq_idx,\
                slow_idx=self.parent.slow_idx)
            
            self.fig_dict = self.ccout.plot(maac_th=self.parent.maac_th, fig=self.fig, return_fig_dict=True)
            
            # add navigation
            self.axes_ = [ax["axis"] for _, ax in self.fig_dict.items()]
            self.nav = Navigate(self.axes_, self, color='red',\
            linewidth=0.5, alpha=0.5)

            self.draw()


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
            dist = abs(self.ticks['left'] - self.ticks['right'])*60
            print(f" distance |R-L|  >>  {dist}  [sec]")


    def on_click(self, event):
        if event.inaxes:
            if event.inaxes in self.axes_:
                
                # t = mdates.num2date(float(event.xdata))
                # t = t.replace(tzinfo=None)

                if event.button == 1:
                    self.ticks['left'] = event.xdata
                
                if event.button == 3:
                    self.ticks['right'] = event.xdata
                
                self.print_ticks()
    

    def nidx_from_time(self, time_in_minutes):
        dtime    = self.ccout._dout["dtime"]
        duration = (dtime[-1]-dtime[0]).total_seconds()/60
        npts     = len(dtime)
        time     = np.linspace(0, duration, npts)
        return np.argmin(np.abs(time_in_minutes-time))


    def on_hover(self, event):
        attr = self.get_axis_attr(event.inaxes)

        if attr:
            sc = self.fig_dict[attr]["sc"]
            cont, ind = sc.contains(event)
            if cont:
                colorlist = ["blue"]*len(self.ccout._dout["dtime"])
                x, _ = sc.get_offsets()[ind["ind"][0]]
                nidx = self.nidx_from_time(x)
                colorlist[nidx] = "green"
                
                for _, fdict in self.fig_dict.items():
                    fdict["sc"].set_facecolor(colorlist)

                self.hover_nidx = nidx
        
                self.draw()


    def on_key(self, event):
        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)
        

        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)

        if event.key == "s" and self.hover_nidx:
            self.WidgetNidex = CC8nidxWidget(self.ccout, self.hover_nidx)
            self.WidgetNidex.show()



