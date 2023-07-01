#!/usr/bin/env python3
# coding=utf-8


import pyqtgraph
import datetime as dt

# matplotlib
import numpy as np
import matplotlib.dates as mdates
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.text import Annotation
from .utils import Navigate, notify

 
class CC8nidxWidget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS y controla los botones
    """

    def __init__(self, cc8out, nidx_list, parent=None):
        QtWidgets.QWidget.__init__(self)
        self.cc8out    = cc8out
        self.parent    = parent  # main canvas
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

        self.canvas.plot(nidx)
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


class CC8nidxCanvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.fig = Figure(figsize=(9,9))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
    
    def on_key(self, event):
        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)

        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)


    def print_ticks(self):
        print("\n :: Ticks info :: ")
        for position, time in self.ticks.items():
            if time:
                print(f"  {position} >> {time} ")
        
        if self.ticks['right'] and self.ticks['left']:
            dist = abs((self.ticks['left'] - self.ticks['right']))
            print(f" distance |R-L|  >>  {dist}  [sec]")
    

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
        with pyqtgraph.BusyCursor():
            self.fig.clf()
            self.ticks = dict(right=None, left=None)
            axes = self.fig.subplots(4,2, gridspec_kw={"height_ratios":[1,0.1,0.75,0.75], "width_ratios":[1,0.02]})
            axes[1,0].axis("off")
            axes[1,1].axis("off")
            axes[2,1].axis("off")
            axes[3,1].axis("off")
            self.parent.cc8out.plot_smap(nidx, axis=axes[0,0], bar_axis=axes[0,1], fig=self.fig)
            self.parent.cc8out.plot_wvfm(nidx, off_sec=5, axes=[axes[2,0], axes[3,0]], fig=self.fig, show_title=False)
            self.nav = Navigate([axes[2,0], axes[3,0]], self, color='red', linewidth=0.5, alpha=0.5)
            self.draw()
            self.parent.show()
        
        if self.parent.parent:
            self.parent.parent.show_green(nidx)


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

    
    def on_key(self, key):
        if key == QtCore.Qt.Key_Right:
            if self.starttime + self.interval > self.cc8.time_[-1]:
                notify("CC8", "The bound of the file was reached!")
                endtime = self.cc8.time_[-1]
                self.starttime = endtime - self.interval
            else:
                endtime = None
                self.starttime = self.starttime + self.interval
            
            self.canvas.plot(endtime=endtime)

        if key == QtCore.Qt.Key_Left:
            if self.starttime - self.interval + self.olap < self.cc8.time_[0]:
                notify("CC8", "The bound of the file was reached!")
                self.starttime = self.cc8.time_[0]
            else:
                self.starttime = self.starttime - self.interval + self.olap
            self.canvas.plot()


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

        if not endtime:
            self.endtime = self.parent.starttime + self.parent.interval
        else:
            self.endtime = endtime
        
        with pyqtgraph.BusyCursor():

            self.ccout = self.parent.cc8.get(starttime=self.parent.starttime,\
                endtime=self.endtime, slowmap=True, fq_idx=self.parent.fq_idx,\
                slow_idx=self.parent.slow_idx)
            
            # set nidx widget
            self.nidx_list = self.ccout.get_nidx(maac_th=self.parent.maac_th, baz_int=self.parent.baz_int)
            
            if not self.WidgetNidex:
                # open a widget
                self.WidgetNidex = CC8nidxWidget(self.ccout, self.nidx_list, parent=self)
            else:
                self.WidgetNidex.update(self.ccout, self.nidx_list)

            self.time = np.array(self.ccout._dout["dtime"])
            self.fig_dict = self.ccout.plot(maac_th=self.parent.maac_th, baz_int=self.parent.baz_int, fig=self.fig, x_time=True, return_fdict=True)
            
            # add navigation
            self.axes_ = [ax["axis"] for _, ax in self.fig_dict.items()]
            self.nav = Navigate(self.axes_, self, color='red', active_cursor=False, linewidth=0.5, alpha=0.5)

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
        colorlist = ["blue"]*len(self.ccout._dout["dtime"])
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
                x, _ = sc.get_offsets()[ind["ind"][0]]
                x = mdates.num2date(float(x))
                x = x.replace(tzinfo=None)
                nidx  = np.argmin(np.abs(x - self.time))
                self.hover_nidx = nidx
                self.show_green(nidx)


    def on_key(self, event):
        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)

        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)

        if event.key == "s" and self.hover_nidx:
            self.WidgetNidex.plot(nidx=self.hover_nidx)
            


