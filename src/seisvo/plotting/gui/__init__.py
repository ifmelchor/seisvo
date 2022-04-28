#!/usr/bin/env python
# coding=utf-8

import os
import numpy as np
import matplotlib.colors as mcolor

from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib.widgets import AxesWidget

from seisvo.utils.plotting import plot_multiple_psd
from seisvo.gui.frames import psdfull_gui

from pynotifier import Notification
import pyqtgraph
import pyautogui

path = os.path.dirname(os.path.realpath(__file__))
icons_path = os.path.join(path, 'icons')


def notify(title, description, status='info', duration=2):
    """
    status: error, info, warn
    """
    n = Notification(
        title=title, 
        description=description, 
        icon_path='%s/%s.png' %(icons_path, status),
        duration=duration
        )
    n.send()


class Navigation(object):
    def __init__(self, ax, imshow_axes=None, base_scale=2, parent=None, **kwargs):
        
        self.parent = parent
        self.ax = ax
        self.imshow_axes = imshow_axes
        self.base_scale = base_scale

        if not (isinstance(ax, np.ndarray) or isinstance(ax, list)):
            self.ax = [ax]
        
        if not (isinstance(imshow_axes, np.ndarray) or isinstance(imshow_axes, list)) and imshow_axes:
            self.imshow_axes = [imshow_axes]
        
        self.fig = self.ax[0].get_figure()
        self.cursors = [Cursor(ax, useblit=True, color='red', linewidth=0.5, alpha=0.5) for ax in self.ax]
        self.max_xlim = self.ax[0].get_xlim()

        self.press = None
        self.cur_xlim = None
        self.x0 = None
        self.x1 = None
        self.xpress = None
        self.new_xlim = None

        self.ticks = dict(left=[None, []], right=[None, []]) # tick data and list of axvline
        self.fig.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.fig.canvas.setFocus()
        self.fig.canvas.mpl_connect('scroll_event', self.onZoom)
        self.fig.canvas.mpl_connect('button_press_event', self.onClkPress)
        self.fig.canvas.mpl_connect('button_release_event', self.onClkRelease)
        self.fig.canvas.mpl_connect('motion_notify_event', self.onMotion)


    def reset(self):
        for ax in self.ax:
            ax.set_xlim(self.max_xlim)
        
        if self.imshow_axes:
            for imax in self.imshow_axes:
                imax.set_xlim(self.max_xlim)


    def reset_ticks(self):
        if self.ticks['left'][0]:
            [line.remove() for line in self.ticks['left'][1]]
            self.ticks['left'][0] = None
        
        if self.ticks['right'][0]:
            [line.remove() for line in self.ticks['right'][1]]
            self.ticks['right'][0] = None

    def draw(self):
        with pyqtgraph.BusyCursor():
            self.fig.canvas.draw()


    def set_xlim(self, new_lim):

        if new_lim[0] < self.max_xlim[0]:
            new_lim[0] = self.max_xlim[0]

        if new_lim[1] > self.max_xlim[1]:
            new_lim[1] = self.max_xlim[1]

        for ax in self.ax:
            ax.set_xlim(new_lim)
        
        if self.imshow_axes:
            for imax in self.imshow_axes:
                imax.set_xlim(new_lim)


    def onZoom(self, event):
        if event.inaxes not in self.ax: 
            return
        else:
            cur_xlim = event.inaxes.get_xlim()
            xdata = event.xdata # get event x location

            if event.button == 'down':
                # deal with zoom in
                scale_factor = 1 / self.base_scale
            
            elif event.button == 'up':
                # deal with zoom out
                scale_factor = self.base_scale
            
            else:
                return
        
            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            new_lim = [xdata - new_width * (1-relx), xdata + new_width * (relx)]

            self.set_xlim(new_lim)
            self.draw()


    def onClkPress(self, event):
        if event.inaxes not in self.ax:
            return

        if event.dblclick:
            self.reset()
            return

        if event.button == 1:
            if self.ticks['left'][0]:
                self.old_tick = self.ticks['left'][0]
                [line.remove() for line in self.ticks['left'][1]]
            else:
                self.old_tick = None

            self.ticks['left'][0] = event.xdata
            self.ticks['left'][1] = []
            for ax in self.ax:
                self.ticks['left'][1].append(ax.axvline(event.xdata, color='k', alpha=0.7, ls='--', lw=1))

            self.cur_xlim = event.inaxes.get_xlim()
            self.new_xlim = ()

            # this is for release
            self.press = self.x0, event.xdata
            self.x0, self.xpress = self.press

        if event.button == 3:
            #right click
            if self.ticks['right'][0]:
                [line.remove() for line in self.ticks['right'][1]]

            self.ticks['right'][0] = event.xdata
            self.ticks['right'][1] = []
            for ax in self.ax:
                self.ticks['right'][1].append(ax.axvline(event.xdata, color='g', alpha=0.7, ls='--', lw=1))


    def onClkRelease(self, event):
        if event.inaxes not in self.ax: 
            return
        
        if self.new_xlim != ():
            # motion was True, remove new_tick
            self.ticks['left'][0] = None
            [line.remove() for line in self.ticks['left'][1]]

            if self.old_tick:
                self.ticks['left'][0] = self.old_tick
                self.ticks['left'][1] = []
                for ax in self.ax:
                    self.ticks['left'][1].append(ax.axvline(self.old_tick, color='k', alpha=0.7, ls='--', lw=1))

        self.press = None
        self.draw()


    def onMotion(self, event):
        if event.inaxes not in self.ax: 
            return

        if self.press is None:
            return

        dx = event.xdata - self.xpress

        # if dx == 0:
            # remove left click
            # self.ticks['left'][0] = None
            # self.ticks['left'][1].remove()
            
        self.new_xlim = self.cur_xlim - dx
        self.set_xlim(self.new_xlim)


class Cursor(AxesWidget):
    def __init__(self, ax, horizOn=True, vertOn=True, useblit=False,
                 **lineprops):
        
        AxesWidget.__init__(self, ax)

        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('draw_event', self.clear)

        self.visible = True
        self.horizOn = horizOn
        self.vertOn = vertOn
        self.useblit = useblit and self.canvas.supports_blit

        if self.useblit:
            lineprops['animated'] = True
        self.lineh = ax.axhline(ax.get_ybound()[0], visible=False, **lineprops)
        self.linev = ax.axvline(ax.get_xbound()[0], visible=False, **lineprops)

        self.background = None
        self.needclear = False

    def clear(self, event):
        """Internal event handler to clear the cursor."""
        if self.ignore(event):
            return
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        
        self.linev.set_visible(False)
        self.lineh.set_visible(False)

    def onmove(self, event):
        """Internal event handler to draw the cursor when the mouse moves."""
        if self.ignore(event):
            return
        if not self.canvas.widgetlock.available(self):
            return
        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)

            if self.needclear:
                self.canvas.draw()
                self.needclear = False
            return
        
        self.needclear = True
        
        if not self.visible:
            return
        
        self.linev.set_xdata((event.xdata, event.xdata))
        self.linev.set_visible(self.visible and self.vertOn)

        self.lineh.set_ydata((event.ydata, event.ydata))
        self.lineh.set_visible(self.visible and self.horizOn)

        self._update()

    def _update(self):
        if self.useblit:
            
            if self.background is not None:
                self.canvas.restore_region(self.background)
            
            self.ax.draw_artist(self.linev)
            self.ax.draw_artist(self.lineh)
            self.canvas.blit(self.ax.bbox)

        else:
            self.canvas.draw_idle()
        return False


class LDEfilter(QtWidgets.QDialog):
    def __init__(self, labels, default=[], **kwargs):
        QtWidgets.QDialog.__init__(self)
        self.ui = ldegui_filter.Ui_Dialog()
        self.ui.setupUi(self)
        self.items = []
        self.set_labels(labels)
        self.check_labels(default)
        self.ui.pushButton.clicked.connect(self.set_check)
        self.ui.pushButton_2.clicked.connect(self.set_uncheck)
        self.show()

    
    def set_labels(self, labels):
        j = 0
        for i in range(int(len(labels)/2)):
            self.add_items((labels[j], labels[j+1]))
            j += 2

        if len(labels) % 2 == 1:
            self.add_items(labels[-1])


    def add_items(self, item):
        root = self.ui.treeWidget.invisibleRootItem()
        child_count = root.childCount()

        if isinstance(item, tuple):
            # add item[0] and item[1] in two-colunm
            item_0 = QtWidgets.QTreeWidgetItem(self.ui.treeWidget)
            item_0.setCheckState(0, QtCore.Qt.Unchecked)
            item_0.setCheckState(1, QtCore.Qt.Unchecked)
            self.items += [(item_0, item)]
            __sortingEnabled = self.ui.treeWidget.isSortingEnabled()
            self.ui.treeWidget.setSortingEnabled(False)
            self.ui.treeWidget.topLevelItem(child_count).setText(0, item[0])
            self.ui.treeWidget.topLevelItem(child_count).setText(1, item[1])
            self.ui.treeWidget.setSortingEnabled(__sortingEnabled)

        else:
            # add item in first column
            item_0 = QtWidgets.QTreeWidgetItem(self.ui.treeWidget)
            item_0.setCheckState(0, QtCore.Qt.Unchecked)
            self.items += [(item_0, item)]
            __sortingEnabled = self.ui.treeWidget.isSortingEnabled()
            self.ui.treeWidget.setSortingEnabled(False)
            self.ui.treeWidget.topLevelItem(child_count).setText(0, item)
            self.ui.treeWidget.setSortingEnabled(__sortingEnabled)


    def check_labels(self, items):
        if items:
            for it in self.items:
                if isinstance(it[1], tuple):
                    if it[1][0] in items:
                        it[0].setCheckState(0, QtCore.Qt.Checked)

                    if it[1][1] in items:
                        it[0].setCheckState(1, QtCore.Qt.Checked)
                else:
                    if it[1] in items:
                        it[0].setCheckState(0, QtCore.Qt.Checked)

        else:
            for it in self.items:
                if isinstance(it[1], tuple):
                    it[0].setCheckState(0, QtCore.Qt.Checked)
                    it[0].setCheckState(1, QtCore.Qt.Checked)
                else:
                    it[0].setCheckState(0, QtCore.Qt.Checked)


    def set_check(self):
        for it in self.items:
            if isinstance(it[1], tuple):
                it[0].setCheckState(0, QtCore.Qt.Checked)
                it[0].setCheckState(1, QtCore.Qt.Checked)
            else:
                it[0].setCheckState(0, QtCore.Qt.Checked)


    def set_uncheck(self):
        for it in self.items:
            if isinstance(it[1], tuple):
                it[0].setCheckState(0, QtCore.Qt.Unchecked)
                it[0].setCheckState(1, QtCore.Qt.Unchecked)
            else:
                it[0].setCheckState(0, QtCore.Qt.Unchecked)


    def get_labels(self):
        checked_labels = []
        for it in self.items:
            if isinstance(it[1], tuple):
                if it[0].checkState(0) == QtCore.Qt.Checked:
                    checked_labels += [it[1][0]]
                
                if it[0].checkState(1) == QtCore.Qt.Checked:
                    checked_labels += [it[1][1]]
            else:
                if it[0].checkState(0) == QtCore.Qt.Checked:
                    checked_labels += [it[1]]

        return checked_labels


class PSD_GUI(QtWidgets.QMainWindow):
    def __init__(self, freq, psd_array, pd_array=None, **kwargs):
        QtWidgets.QMainWindow.__init__(self)
        self.ui = psdfull_gui.Ui_MainWindow()
        self.ui.setupUi(self)
        self.canvas = PSDCanvas(freq, psd_array, pd_array, parent=self, **kwargs)
        self.ui.verticalLayout.addWidget(self.canvas)
        self.ui.pushButton.clicked.connect(self.canvas.logPSD)
        self.ui.pushButton2.clicked.connect(self.canvas.logFreq)

        self.setWindowTitle(kwargs.get('title', 'PSD--PD'))


class PSDCanvas(FigureCanvas):
    def __init__(self, freq, psd_array, pd_array=None, parent=None, **kwargs):

        figsize = kwargs.get("figsize", (8, 4))
        dpi = kwargs.get("dpi", 100)
        plot_prob = kwargs.get('plot_prob', True)

        self.freq = freq
        self.psd_array = psd_array
        self.pd_array = pd_array
        self.log_psd = kwargs.get('log_psd', False)
        self.log_fq = kwargs.get('log_fq', False)
        
        self.plot_kwargs = dict(figsize = figsize,
                                dpi = dpi,
                                plot_prob = plot_prob,
                                log_psd = kwargs.get('log_psd', False),
                                log_fq = kwargs.get('log_fq', False),
                                colors = kwargs.get('colors', []),
                                )

        self.fig = Figure(figsize=figsize, dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, 
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        # self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()
        self.mpl_connect('button_release_event', self.print_tickinfo)
        self.__plot__()

    def __plot__(self):
        self.fig.clf()
        with pyqtgraph.BusyCursor():
            _, self.axes = plot_multiple_psd(self.freq, self.psd_array, self.pd_array, fig=self.fig, 
                plot=False, **self.plot_kwargs)

        self.nav_n = Navigation(self.axes, parent=self.fig)
        self.draw()

    def logPSD(self):
        if self.log_psd:
            self.log_psd = False
        else:
            self.log_psd = True

        self.plot_kwargs['log_psd'] = self.log_psd
        self.__plot__()

    def logFreq(self):
        if self.log_fq:
            self.log_fq = False
        else:
            self.log_fq = True

        self.plot_kwargs['log_fq'] = self.log_fq
        self.__plot__()


    def print_tickinfo(self, event):
        if isinstance(self.axes, np.ndarray):
            if event.inaxes == self.axes[0]:
                print('--------Tick-Info---------')
                try:
                    bfq = self.nav_n.ticks['left'][0]
                    print('   Black tick: %2.2f Hz' % bfq)
                except:
                    bfq = None

                try:
                    gfq = self.nav_n.ticks['right'][0]
                    print('   Green tick: %2.2f Hz' % gfq)
                except:
                    gfq = None

                if bfq and gfq:
                    print('   Delta: %2.2f Hz' % np.abs(bfq-gfq))
                print('-------------------------')
        else:
            if event.inaxes == self.axes:
                print('--------Tick-Info---------')
                try:
                    bfq = self.nav_n.ticks['left'][0]
                    print('   Black tick: %2.2f Hz' % bfq)
                except:
                    bfq = None

                try:
                    gfq = self.nav_n.ticks['right'][0]
                    print('   Green tick: %2.2f Hz' % gfq)
                except:
                    gfq = None

                if bfq and gfq:
                    print('   Delta: %2.2f Hz' % np.abs(bfq-gfq))
                print('-------------------------')
