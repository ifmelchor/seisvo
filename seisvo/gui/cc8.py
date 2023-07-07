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
        fig = self.canvas.fig
        files_types = 'Image (*.png)'
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', filename, files_types)
        if fileName == '':
            return
        else:
            fig.savefig(fileName+".pdf", bbox_inches='tight', dpi = 300)


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
        
        if event.key == 'f1':
            filename = f"Slowmap_1"
            self.parent.save_fig(filename)


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

            try:
                self.parent.cc8out.plot_wvfm(nidx, off_sec=5, axes=[axes[2,0], axes[3,0]], fig=self.fig, show_title=False)
                self.nav = Navigate([axes[2,0], axes[3,0]], self, color='red', linewidth=0.5, alpha=0.5)
            except:
                print(" No data found to plot waveform.")

            self.draw()
            self.parent.show()
        
        if self.parent.parent:
            self.parent.parent.show_green(nidx)


class CC8Widget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS y controla los botones
    """

    def __init__(self, cc8, starttime, interval, fq_idx, olap=0.1, maac_th=0.6, max_err=0.7, rms_lim=[], baz_int=[]):

        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()
        
        self.cc8       = cc8
        self.fq_idx    = fq_idx
        self.maac_th   = maac_th
        self.max_err   = max_err
        self.baz_int   = baz_int
        self.starttime = starttime
        self.interval  = dt.timedelta(minutes=interval)
        self.olap      = dt.timedelta(minutes=interval*olap)

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
            if self.starttime + self.interval > self.cc8.time_[-1]:
                notify("CC8", "The bound of the file was reached!")
                endtime = self.cc8.time_[-1]
                self.starttime = endtime - self.interval - self.olap
            else:
                endtime = None
                self.starttime = self.starttime + self.interval - self.olap

            self.canvas.plot(endtime=endtime)

        if key == QtCore.Qt.Key_Left:
            if self.starttime - self.interval < self.cc8.time_[0]:
                notify("CC8", "The bound of the file was reached!")
                self.starttime = self.cc8.time_[0]
            else:
                self.starttime = self.starttime - self.interval + self.olap
            
            self.canvas.plot()


    def update_maac(self, new_maac):
        self.maac_th   = new_maac
        self.canvas.plot()
    

    def update_error(self, new_error):
        self.max_err   = new_error
        self.canvas.plot()


    def save_fig(self, filename):
        fig = self.canvas.fig
        files_types = 'Image (*.png)'
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', filename, files_types)
        if fileName == '':
            return
        else:
            fig.savefig(fileName, bbox_inches='tight', dpi = 300)


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

            self.ccout = self.parent.cc8.get(starttime=self.parent.starttime, endtime=self.endtime, slowmap=True, fq_idx=self.parent.fq_idx)
            
            # set nidx widget
            self.nidx_list = self.ccout.get_nidx(return_full=False, max_err=self.parent.max_err, maac_th=self.parent.maac_th, baz_int=self.parent.baz_int)
            
            if not self.WidgetNidex:
                # open a widget
                self.WidgetNidex = CC8nidxWidget(self.ccout, self.nidx_list, parent=self)
            else:
                self.WidgetNidex.update(self.ccout, self.nidx_list)

            self.time = np.array(self.ccout._dout["dtime"])
            self.fig_dict = self.ccout.plot(max_err=self.parent.max_err, maac_th=self.parent.maac_th, baz_int=self.parent.baz_int, fig=self.fig, x_time=True, rms_lim=[self.parent.rms_min, self.parent.rms_max], return_fdict=True)

            # add horizontal bar in maac_th
            self.fig_dict["maac"]["axis"].axhline(self.parent.maac_th, color="r", ls="--", alpha=0.7, zorder=5)
            
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
                self.hover_time = x
                self.baz0  = self.fig_dict["bazm"]["sc"].get_offsets()[self.hover_nidx][1]
                self.slow0 = self.fig_dict["slow"]["sc"].get_offsets()[self.hover_nidx][1]
                self.show_green(nidx)


    def on_key(self, event):
        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)


        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)


        if event.key == "s" and self.hover_nidx:
            self.WidgetNidex.plot(nidx=self.hover_nidx)
        

        if event.key == "f1":
            filename = "Detecciones_1.png"
            self.parent.save_fig(filename)


        if event.key == "t" and self.hover_nidx:
            print(f"\n >> {self.hover_time.strftime('%Y %b %d %H:%M:%S.%f')} :: BAZ {self.baz0:5.1f} SLOW {self.slow0:4.2f}")


        if event.key == "m":
            new_maac, ok = QtWidgets.QInputDialog.getDouble(self, "MAAC threshold","Value:", self.parent.maac_th, 0.0, 1.0, 2, QtCore.Qt.WindowStaysOnTopHint, 0.05)
            if ok:
                self.parent.update_maac(new_maac)
        

        if event.key == "e":
            new_error, ok = QtWidgets.QInputDialog.getDouble(self, "Error max. threshold","Value:", self.parent.max_err, 0.0, 10.0, 1, QtCore.Qt.WindowStaysOnTopHint, 0.1)
            if ok:
                self.parent.update_error(new_error)


        if event.key == "p" and self.baz0 and self.slow0:
            if self.ticks['right'] and self.ticks['left']:
                starttime = min([self.ticks['right'],self.ticks['left']])
                endtime   = max([self.ticks['right'],self.ticks['left']])
                fig = self.ccout.get_beamform(starttime, endtime, self.slow0, self.baz0)
                fig.savefig(f"{self.parent.cc8.stats.id}-{starttime.strftime('%y%m%d-%H%M%S')}-BAZ{self.baz0:.0f}.png")
        

        if event.key == "z":
            fig = self.ccout.prob_slowmap(fq_idx=self.parent.fq_idx, max_err=self.parent.max_err, maac_th=self.parent.maac_th, baz_int=self.parent.baz_int)

