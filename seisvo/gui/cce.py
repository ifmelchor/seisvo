#!/usr/bin/env python3
# coding=utf-8


import operator
import pyqtgraph
from matplotlib.gridspec import GridSpec
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from seisvo.gui.utils import Navigate


class _CCEiter:
    def __init__(self, cce, eid_list, eid=None):
        self.eid_list = eid_list
        self.cce = cce

        if not eid:
            self.n = 0
        else:
            self.n = self.eid_list.index(eid)


    def take(self):
        return self.cce[self.eid_list[self.n]]


    def get_eidlist(self):
        return [str(x) for x in self.eid_list]


    def next(self):
        if self.n == len(self.eid_list)-1:
            self.n = 0
        else:
            self.n += 1

        return self.take()


    def prev(self):
        if self.n == 0:
            self.n = len(self.eid_list)-1
        else:
            self.n -= 1

        return self.take()


class CCEWidget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS y controla los botones
    """

    def __init__(self, cce, eid_list, path_to_cc8file="./"):
        QtWidgets.QWidget.__init__(self)
        
        self.cce  = cce
        self.ptc  = path_to_cc8file
        self.iter = _CCEiter(cce, eid_list)
        self.eid  = self.iter.take()
        
        self.layout = QtWidgets.QVBoxLayout()
        self.canvas = CCECanvas(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        self.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.setFocus()


    def select_eid(self):
        eid , ok = QtWidgets.QInputDialog.getItem(self.canvas, 'Cambiar de EID', 'Selecciona:', self.iter.get_eidlist(), self.iter.n, editable=False)
        eid = int(eid)
        if ok and eid != self.eid:
            self.update(self.iter.eid_list, eid)


    def update(self, eid_list, init_eid):
        self.iter = _CCEiter(self.cce, eid_list, eid=init_eid)
        self.eid  = self.iter.take()
        self.canvas.plot()


    def relabel_eid(self, label):
        self.cce.relabel(self.eid.id, label.upper())
        self.eid  = self.iter.take()


    def onKey(self, key):
        if key == QtCore.Qt.Key_Right:
            self.eid = self.iter.next()
            self.canvas.plot()
            
        if key == QtCore.Qt.Key_Left:
            self.eid = self.iter.prev()
            self.canvas.plot()

        self.canvas.setFocus()


class CCECanvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.fig    = Figure(figsize=(12,9))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
        self.taper   = False
        self.off_sec = 5 # sec   
        self.plot()
        self.setFocus()


    def plot(self):
    
        self.parent.setWindowTitle(f" ID [{self.parent.eid.id}] LABEL [{self.parent.eid.label}]")
        
        with pyqtgraph.BusyCursor():

            self.fig.clf()
            self.ticks = dict(right=None, left=None)
            
            gs  = GridSpec(3, 4, wspace=0.5, width_ratios=[1, 2, 0.1, 1], height_ratios=[2, 1, 1])
            ax1 = self.fig.add_subplot(gs[0,1])
            ax2 = self.fig.add_subplot(gs[0,2])
            ax3 = self.fig.add_subplot(gs[1,:])
            ax4 = self.fig.add_subplot(gs[2,:])

            self.axes = {
                "ax1":ax1,
                "ax2":ax2,
                "ax3":ax3,
                "ax4":ax4
            }

            self.parent.eid.slowmap(path_to_cc8file=self.parent.ptc, axis=ax1, fig=self.fig, bar_axis=ax2)

            self.parent.eid.beamform(off_sec=self.off_sec, path_to_cc8file=self.parent.ptc, taper=self.taper, plot=True, fig=self.fig, axes=[ax3, ax4])
            
            self.nav = Navigate([ax3, ax4], self, color='red', linewidth=0.5, alpha=0.5)

            self.draw()


    def print_ticks(self):
        try:
            print("\n :: Ticks info :: ")
            for position, time in self.ticks.items():
                if time:
                    print(f"  {position} >> {time} ")
            
            if self.ticks['right'] and self.ticks['left']:
                dist = abs((self.ticks['left'] - self.ticks['right']))
                print(f" distance |R-L|  >>  {dist:.2f} sec [{1/dist:.2f} Hz]")
        except:
            self.ticks = dict(right=None, left=None)


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


    def on_key(self, event):

        if event.key == 'right':
            self.parent.onKey(QtCore.Qt.Key_Right)


        if event.key == 'left':
            self.parent.onKey(QtCore.Qt.Key_Left)


        if event.key == "t":
            self.taper = operator.not_(self.taper)
            self.plot()


        if event.key == "s":
            ans, ok = QtWidgets.QInputDialog.getDouble(self, 'Change offset time', 'Enter new window in seconds:', self.off_sec, 0, 600, 2)
            if ok:
                self.off_sec = ans
                self.plot()


        if event.key == "l":
            label , ok = QtWidgets.QInputDialog.getText(self,\
                            f'Evento {self.parent.eid.id}', 'Define Etiqueta:', text=f"{self.parent.eid.label}")
            if ok:
                self.parent.relabel_eid(label.upper())
                self.parent.setWindowTitle(f" ID [{self.parent.eid.id}] LABEL [{self.parent.eid.label}]")


        if event.key == "tab":
            self.parent.select_eid()
            

