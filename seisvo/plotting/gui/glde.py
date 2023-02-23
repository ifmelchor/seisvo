#!/usr/bin/env python
# coding=utf-8

import sys
import pyqtgraph
import numpy as np
import datetime as dt
import collections
import calendar
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
import matplotlib.colors as mcolor
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.qt_compat import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from seisvo import Network
from seisvo.signal.pdf import get_max_prob
from seisvo.gui import Navigation, PSD_GUI, LDEfilter
from seisvo.gui.frames import ldewidget, ldegui_filter, event_polar
from seisvo.gui import notify, get_norm
from functools import partial

plt.rc('axes', labelsize=9)
plt.rc('axes', labelpad=4.0)
plt.rc('axes', titlepad=6.0)
plt.rc('axes', titlesize=10)
plt.rc('xtick', labelsize=9)
plt.rc('ytick', labelsize=9)
plt.rc('lines', linewidth=0.5)
plt.rc('lines', color='black')
plt.rc('lines', markersize=4)
plt.rc('legend', markerscale=1.5)

# TOF classification
TYPES_EVENTS = ('SPT', 'MPT', 'BBT',
                'BT', 'ST', 'LPS', 'UK')
MEAN_COLOR = 'navy'
GAS_COLOR = 'lightskyblue'
ASH_COLOR = 'dimgrey'
GEYSER_COLOR = 'b'
SO2_COLOR = 'k'
T_COLOR = 'darkorange'
IC_COLOR = 'firebrick'
EXP_COLOR = 'blueviolet'
S02_label = r'S0$_2$'+'\nanomaly'
TA_label = 'Thermal\nanomaly'
IC_label = 'Incandescence'

COLOR_DICT = {'Z':'k', 'E':'navy', 'N':'firebrick'}


class LDEWindow(QtWidgets.QMainWindow):
    def __init__(self, lde, **kwargs):
        QtWidgets.QMainWindow.__init__(self)
        self.ui = ldewidget.Ui_MainWindow()
        self.ui.setupUi(self)

        self.canvas = LDEevent(self, lde, **kwargs)

        self.ui.verticalLayout.addWidget(self.canvas)
        self.ui.action_Open.triggered.connect(lde.open)
        
        self.ui.actionPlot_Scatter.triggered.connect(lde.scatter)

        self.ui.actionRelabel.triggered.connect(self.canvas._relabel)
        self.ui.actionRelabel.setShortcut('R')
        
        self.ui.actionRemove.triggered.connect(self.canvas._remove)
        self.ui.actionRemove.setShortcut('del')

        self.ui.actionPlotPolar.triggered.connect(self.canvas._plot_polar)
        self.ui.actionPlotPolar.setShortcut('ctrl+p')

        self.ui.actionGo_to_ID.triggered.connect(self.canvas._goto)
        self.ui.actionGo_to_ID.setShortcut('tab')
        
        self.ui.actionNext.triggered.connect(self.canvas._next)
        self.ui.actionNext.setShortcut('shift+right')
        self.ui.actionPrevious.triggered.connect(self.canvas._previous)
        self.ui.actionPrevious.setShortcut('shift+left')

        self.ui.actionODP.triggered.connect(partial(self.canvas._change_modeview, "odp"))
        self.ui.actionODP.setShortcut('ctrl+o')
        self.ui.actionStep_Range.triggered.connect(partial(self.canvas._change_modeview, "steprange"))
        self.ui.actionStep_Range.setShortcut('ctrl+s')
        self.ui.actionFull_Range.triggered.connect(partial(self.canvas._change_modeview, "fullrange"))
        self.ui.actionFull_Range.setShortcut('ctrl+f')

        self.ui.actionFilter.triggered.connect(self.canvas._filter)
        self.ui.actionFilter.setShortcut('f')

        self.ui.actionPlot_3D.setCheckable(True)
        self.ui.actionPlot_3D.setChecked(int(kwargs.get('three_component', False)))
        self.ui.actionPlot_3D.triggered.connect(self.canvas._three_component)

        self.ui.actionOpen_LTE.triggered.connect(self.canvas.openLTE)
        self.ui.actionOpen_LTE.setShortcut('F1')

        self.lte_parent = kwargs.get("lte_parent", None)
        if self.lte_parent:
            self.ui.actionOpen_LTE.setEnabled(False)

        self.ui.actionExit.triggered.connect(self.close)
        self.ui.actionExit.setShortcut('esc')


    def closeEvent(self, event):
        if self.lte_parent:
            self.lte_parent.setEnabled(True)
            self.lte_parent.canvas.plot()


class LDEevent(FigureCanvas):
    def __init__(self, parent, lde, **kwargs):
        self.fig = Figure(figsize=(20,9), subplotpars=SubplotParams(left=0.05, right=0.95, wspace=0.1))
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        self.lde = lde
        self.parent = parent
        self.plot_kwargs = dict(vmin=kwargs.get("vmin", -20),
                                vmax=kwargs.get("vmax", 60),
                                interval=kwargs.get("interval", 60),
                                base_scale=kwargs.get("base_scale", 2),
                                three_component=kwargs.get("three_component", False),
                                fig=self.fig)

        event_id = kwargs.get("event_id", None)
        self.psd_frame = None

        self.event_labels = kwargs.get("event_labels", [])
        self.event_time_interval = kwargs.get("event_time", ())
        self._set_event(id=event_id)

        self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self.callbacks.connect('key_press_event', self.on_key)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('button_release_event', self.show_tickinfo)
        self.callbacks.connect('scroll_event', self.on_zoom)

        self.spec_info = {}
        self.plot(default=True)


    def _set_statusbar(self):
        text = '\tID: %s | Duration: %.1f hr' % (self.evnt.id, self.evnt.duration)
        self.parent.ui.statusbar.showMessage(text, 0)


    def _set_modeview(self):
        if self.mode == 'fullrange':
            self.parent.ui.actionFull_Range.setChecked(True)
            self.parent.ui.actionFull_Range.setEnabled(False)
            self.parent.ui.actionStep_Range.setEnabled(True)
            self.parent.ui.actionStep_Range.setChecked(False)
            self.parent.ui.actionODP.setEnabled(True)
            self.parent.ui.actionODP.setChecked(False)

        if self.mode == 'steprange':
            self.parent.ui.actionFull_Range.setChecked(False)
            self.parent.ui.actionFull_Range.setEnabled(True)
            self.parent.ui.actionStep_Range.setEnabled(False)
            self.parent.ui.actionStep_Range.setChecked(True)
            self.parent.ui.actionODP.setEnabled(True)
            self.parent.ui.actionODP.setChecked(False)

        if self.mode == 'odp':
            self.parent.ui.actionFull_Range.setChecked(False)
            self.parent.ui.actionFull_Range.setEnabled(True)
            self.parent.ui.actionStep_Range.setEnabled(True)
            self.parent.ui.actionStep_Range.setChecked(False)
            self.parent.ui.actionODP.setEnabled(False)
            self.parent.ui.actionODP.setChecked(True)


    def _set_nroevent(self):
        text = '%s/%s' % (self.event_index+1, len(self.event_list))
        self.psd_ax[0].annotate(text, xy=(0.95, 0.95),
            xycoords='figure fraction', bbox=dict(boxstyle="round", fc=(1.0, 0.7, 0.7), ec="none"))


    def _is_event(self):
        event_list = self.lde.get_events(label=self.event_labels,
        time_interval=self.event_time_interval)

        if event_list:
            return True
        else:
            return False


    def _change_modeview(self, mode):
        if mode == 'steprange':
            step, ok = QtWidgets.QInputDialog.getInt(self, "Set step time", "Time in min: ", value=20)
            if ok:
                self.step_time = step
            else:
                return

        self.start = None
        self.mode = mode
        self._set_modeview()
        self.plot()


    def _set_event(self, id=None):
        if id and self.lde.is_event(id):
            evnt = self.lde.get_event(id)
        else:
            evnt = None

        self.event_list = self.lde.get_events(label=self.event_labels,
        time_interval=self.event_time_interval)

        self.nro_events = len(self.event_list)
        self.start = None

        if evnt:
            for e in self.event_list:
                if e.id == evnt.id:
                    self.event_index = self.event_list.index(e)
                    return
                
            notify('seisvo', '(**) Event ID not found', status='warn')

        self.event_index = 0

        if self.event_list[self.event_index].lte_file_sup:
            self.parent.ui.actionPlotPolar.setEnabled(True)
        else:
            self.parent.ui.actionPlotPolar.setEnabled(False)


    def _plot_polar(self):
        lde_window = LDEpolarWindow(self.evnt, polar_min=0.9)
        lde_window.show()


    def _three_component(self):
        self.plot_kwargs['three_component'] = self.parent.ui.actionPlot_3D.isChecked()
        self.plot()


    def plot(self, default=False, clean=False):
        if clean:
            self.fig.clf()
            self.draw()
            return

        self.axv = []
        self.nav = []
        self.evnt = self.event_list[self.event_index]

        if default:
            if self.evnt.duration > 12:
                self.mode = 'odp'
            else:
                self.mode = 'fullrange'

        else:
            if self.evnt.duration > 20 and self.mode == 'fullrange':
                self.mode = 'odp'
                notify('seisvo', 'mode ODP activated', status='info')

        self.plot_kwargs['mode'] = self.mode

        if self.mode == 'steprange':
            if not self.start:
                self.start = self.evnt.starttime

            self.plot_kwargs['start'] = self.start
            self.plot_kwargs['step'] = self.step_time
            axes, self.values, start = plot_event(self.evnt, **self.plot_kwargs)
            self.start = start
        
        else:
            self.plot_kwargs['start'] = None
            self.plot_kwargs['step'] = None
            axes, self.values, _ = plot_event(self.evnt, **self.plot_kwargs)

        (self.tr_ax, self.psd_ax, spec_ax) = axes
        self.fq_range = (self.values[0][1].min(), self.values[0][1].max())

        if self.mode in ('fullrange', 'steprange'):
            nav1 = Navigation(self.tr_ax, parent=self, specgram_axis=spec_ax)
            self.nav += [nav1]
        
        self._set_modeview()
        self._set_statusbar()
        self._set_nroevent()
        
        with pyqtgraph.BusyCursor():
            self.draw()


    def _filter(self):
        label_dict = self.lde.get_labels()
        label_list = list(label_dict.keys())
        label_list.sort()
        filt = LDEfilter(label_list, default=self.event_labels)
        ans = filt.exec_()

        if ans == QtWidgets.QDialog.Accepted:
            self.event_labels = filt.get_labels()
            self._set_event()
            if self.mode == 'steprange':
                self.plot(default=True)
            else:
                self.plot()


    def _relabel(self):
        # relabel event
        new_label, ok = QtWidgets.QInputDialog.getText(self, "Set new label", "Label :", text=self.evnt.label)

        if ok and new_label.upper() != self.evnt.label:
            self.lde.relabel(self.evnt, new_label.upper())
            notify('seisvo', 'ID %s relabeled!' % self.evnt.id, status='info')

            if self.event_labels and new_label.upper() not in self.event_labels:
                self._set_event()

            else:
                self._set_event(id=self.evnt.id)

            if self._is_event():
                if self.mode == 'steprange':
                    self.plot(default=True)
                else:
                    self.plot()
            else:
                self.plot(clean=True)


    def _remove(self):
        box = QtWidgets.QMessageBox()
        box.setIcon(QtWidgets.QMessageBox.Question)
        box.setWindowTitle('Delete Event %s' % self.evnt.id)
        box.setText('Confirm you really want to delete the event')
        box.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        box.setDefaultButton(QtWidgets.QMessageBox.No)
        buttonYes = box.button(QtWidgets.QMessageBox.Yes)
        buttonYes.setText('Delete')
        buttonNo = box.button(QtWidgets.QMessageBox.No)
        buttonNo.setText('Cancel')
        box.exec_()

        if box.clickedButton() == buttonYes:
            self.lde.remove(event=self.evnt)
            self._set_event()

            if self._is_event():
                if self.mode == 'steprange':
                    self.plot(default=True)
                else:
                    self.plot()
            else:
                self.plot(clean=True)


    def _goto(self):
        new_id, ok = QtWidgets.QInputDialog.getInt(self, "Set Event ID", "ID :", value=self.evnt.id)
        list_ids = [e.id for e in self.event_list]
        list_ids.remove(self.evnt.id)

        if ok and new_id in list_ids:
            try:
                self._set_event(id=new_id)
                if self.mode == 'steprange':
                    self.plot(default=True)
                else:
                    self.plot()
            except:
                notify('seisvo', 'ID not found', status='warn')


    def _next(self):
        self.event_index += 1
        if self.event_index > self.nro_events - 1:
            self.event_index = 0
        
        self.start = None
        if self.mode == 'steprange':
            self.plot(default=True)
        else:
            self.plot()


    def _previous(self):
        self.event_index -= 1
        if self.event_index < 0:
            self.event_index = self.nro_events - 1
        
        self.start = None
        if self.mode == 'steprange':
            self.plot(default=True)
        else:
            self.plot()


    def _nextstep(self):
        self.start += dt.timedelta(minutes=self.step_time)

        if self.start > self.evnt.starttime + dt.timedelta(hours=self.evnt.duration):
            notify('seisvo', 'Start time', status='warn')
            self.start -= dt.timedelta(minutes=self.step_time)
            return

        self.plot()


    def _previousstep(self):
        self.start -= dt.timedelta(minutes=self.step_time)
        
        if self.start < self.evnt.starttime:
            notify('seisvo', 'End time', status='warn')
            self.start += dt.timedelta(minutes=self.step_time)
            return

        self.plot()


    def on_key(self, event):
        if self.mode == 'steprange':
            if event.key == 'right':
                self._nextstep()

            if event.key == 'left':
                self._previousstep()

        if self.mode in ('fullrange', 'steprange'):
            if event.key == 'p':
                try:
                    times = [mdates.num2date(self.nav[0].ticks['left'][0]), 
                             mdates.num2date(self.nav[0].ticks['right'][0])
                             ]
                except:
                    notify('seisvo', 'No clicks found', status='warn')
                    return

                if self.psd_frame:
                    self.psd_frame.close()
                    self.psd_frame = None

                start_tick = min(times)
                end_tick = max(times)

                trace = self.evnt.get_stream(start_tick, end_tick)[0]
                delta = int((end_tick - start_tick).total_seconds()/60)

                if delta >= 5:
                    avg_step = 1
                else:
                    avg_step = None

                psd, freq = trace.psd(fq_band=list(self.evnt.get_lte().stats.fq_band),
                                       avg_step=avg_step, return_fig=False, plot=False)
                psd = np.array([psd])
                self.psd_frame = PSD_GUI(freq, psd, None, plot_prob=False, title='PSD')
                self.psd_frame.show()


    def on_click(self, event):
        if event.inaxes in self.psd_ax:
            (psd, freq0), (degree, freq1) = self.values
            spec = psd[np.argmin(np.abs(np.array(freq0)-event.xdata))]
            pd = degree[np.argmin(np.abs(np.array(freq1)-event.xdata))]

            text = 'Freq.: %.2f\n' % event.xdata
            text +='  PSD: %.2f\n' % spec
            text +='   PD: %.2f'   % pd
            
            if self.spec_info:
                self.spec_info['x_vline'].remove()
                self.spec_info['psd_vline'].remove()
                self.spec_info['pd_vline'].remove()
                self.spec_info['text'].remove()
                self.spec_info = {}

            self.spec_info['x_vline'] = self.psd_ax[0].axvline(event.xdata, color='navy')
            self.spec_info['psd_vline'] = self.psd_ax[0].axhline(spec, color='k')
            self.spec_info['pd_vline'] = self.psd_ax[1].axhline(pd, color='red')
            self.spec_info['text'] = self.psd_ax[0].annotate(text, xy=(0.75, 0.85), 
                xycoords='axes fraction', bbox=dict(boxstyle="round4,pad=.5", fc="0.8"))

            self.draw()


    def show_tickinfo(self, event):
        if event.inaxes == self.tr_ax and self.mode in ('fullrange', 'steprange'):
            tickinfo = ' %s: \n' % self.nav[0].sta_info
            
            try:
                ltime = mdates.num2date(self.nav[0].ticks['left'][1].get_xdata()[0])
                tickinfo += '  B: %s\n' % ltime.strftime("%Y-%m-%d %H:%M:%S")
            except:
                ltime = None
                tickinfo += '  B: None\n'

            try:
                rtime = mdates.num2date(self.nav[0].ticks['right'][1].get_xdata()[0])
                tickinfo += '  G: %s\n' % rtime.strftime("%Y-%m-%d %H:%M:%S")
            except:
                rtime = None
                tickinfo += '  G: None\n'

            if rtime and ltime:
                delta = (max([ltime, rtime]) - min([ltime, rtime])).total_seconds()
                tickinfo += '  Delta: %2.2f sec.\n' % delta
            else:
                delta = None
                tickinfo += '  Delta: None\n'

            print(tickinfo)


    def on_zoom(self, event):
        if event.inaxes in self.psd_ax:

            cur_xlim = self.psd_ax[0].get_xlim()
            xdata = event.xdata

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
            new_min = xdata - new_width * (1-relx)
            new_max = xdata + new_width * (relx)

            if new_min < self.fq_range[0]:
                new_min = self.fq_range[0]

            if new_max > self.fq_range[1]:
                new_max = self.fq_range[1]

            self.psd_ax[0].set_xlim(new_min, new_max)
            self.draw()


    def openLTE(self):
        from seisvo.gui.glte import LTEWindow
        startday = self.evnt.starttime.replace(hour=0, minute=0)
        interval = 15
        lte = self.evnt.get_lte()
        lte_window = LTEWindow(lte, starttime=startday, interval=interval, lde_parent=self.parent)
        lte_window.show()


def plot_event(lde_event, return_fig=False, **kwargs):
    """
    This code crate the plot frame of the event
    """

    mode = kwargs.get('mode', 'fullrange')
    if mode not in ('fullrange', 'steprange', 'odp'):
        mode = 'fullrange'
        print('Warninig: mode changed to fullrange')

    starttime = lde_event.starttime
    endtime = starttime + dt.timedelta(hours=lde_event.duration)
    
    if mode == 'steprange':
        time_delta = kwargs.get('step', 20)
        start = kwargs.get('start', None)
        if not start:
            start = starttime
        end = start + dt.timedelta(minutes=time_delta)
        if end > endtime:
            end = endtime
        starttime = start
        endtime = end

    # get stream
    st = lde_event.get_stream(starttime, endtime)

    # get psd
    psd, freq0 = lde_event.get_specgram(starttime, endtime)
    fq_band = (freq0[0], freq0[-1])
    psd = 10*np.log10(psd)
    psd_pdf = np.array(get_max_prob(psd))
    psd_pdf = psd_pdf/np.nanmax(psd_pdf)

    # get degree
    degree, freq1 = lde_event.get_polargram(starttime, endtime)
    degree = np.array(get_max_prob(degree))

    # create figure
    fig = kwargs.get('fig', None)
    if not fig:
        fig = plt.figure(figsize=kwargs.get('figsize', ()))
        return_axes = False
    else:
        fig.clf()
        return_axes = True

    title = "Event: %s ID: %s\n Time: %s" % (lde_event.label, lde_event.id, lde_event.starttime.strftime('%d %b %Y'))
    fig.suptitle(title)
    
    if mode == 'odp':
        interval = kwargs.get('interval', 30) # pameter for odp plot
        ncol = 2
        nrow = 1
        grid_opt = None

    else:
        ncol = 4
        nrow = 2
        grid_opt = {'width_ratios':[40,1,3,30], 'height_ratios':[3,2]}
        vmin = kwargs.get('vmin', None)
        vmax = kwargs.get('vmax', None)

    axes = fig.subplots(nrow, ncol, gridspec_kw=grid_opt)

    z_trace = st.get_component('Z')
    if mode in ('fullrange', 'steprange'):
        three_component = kwargs.get('three_component', False)

        gs = axes[1,2].get_gridspec()
        for ax in axes[:, -1]:
            ax.remove()

        psd_ax = fig.add_subplot(gs[:, -1])
        tr_ax = axes[0, 0]
        spec_ax = axes[1,0]

        z_waveform = z_trace.data
        time = z_trace.get_time()

        if three_component:
            z_waveform = z_waveform/np.nanmax(z_waveform)
            trace_list = [(z_waveform, 'Z')]
            yticks, yticklabels = [], []

            try:
                e_trace = st.get_component('E')
                e_waveform = e_trace.data
                e_waveform = e_waveform/np.nanmax(e_waveform)
                trace_list += [(e_waveform, 'E')]
            except:
                pass

            try:
                n_trace = st.get_component('N')
                n_waveform = n_trace.data
                n_waveform = n_waveform/np.nanmax(n_waveform)
                trace_list += [(n_waveform, 'N')]
            except:
                pass

            for n, trace in enumerate(trace_list):
                tr_ax.plot(time, trace[0] - n*7./4, color=COLOR_DICT[trace[1]])
                yticks += [- n*7./4]
                yticklabels += [trace[1]]

            tr_ax.set_yticks(yticks)
            tr_ax.set_yticklabels(yticklabels)
            tr_ax.set_ylim(yticks[-1]-1, +1)

        else:
            tr_ax.plot(time, z_waveform, 'k')
        
        tr_ax.set_xlim(time[0], time[-1])
        axes[0, 1].axes.get_xaxis().set_visible(False)
        axes[0, 1].axes.get_yaxis().set_visible(False)
        axes[0, 1].set_frame_on(False)
        axes[0, 2].axes.get_xaxis().set_visible(False)
        axes[0, 2].axes.get_yaxis().set_visible(False)
        axes[0, 2].set_frame_on(False)
        axes[1, 2].axes.get_xaxis().set_visible(False)
        axes[1, 2].axes.get_yaxis().set_visible(False)
        axes[1, 2].set_frame_on(False)

        im, _ = z_trace.specgram(axes=spec_ax, fq_band=fq_band, per_lap=0.25, v_min=vmin, v_max=vmax)
        cbar = fig.colorbar(im, cax=axes[1,1], label='PSD\n'+r'dB[cnts$^2$/Hz]', orientation='vertical')
        cbar.locator = mtick.MaxNLocator(nbins=4)
        cbar.update_ticks()

        tr_ax.set_xticks(axes[1,0].get_xticks())
        tr_ax.xaxis.set_major_formatter(mtick.NullFormatter())

    else:
        tr_ax = axes[0]
        spec_ax = None
        psd_ax = axes[1]
        z_trace.odp2(tr_ax, interval=interval)

    psd_ax.plot(freq0, psd_pdf, 'k', label='psd norm')
    psd_ax.grid(axis='both', which='major', color='k', linestyle='--', alpha=0.3)
    psd_ax.set_ylabel('Normalized PSD [dB]')
    psd_ax.set_xlabel('Freq. [Hz]')
    psd_ax.set_ylim(0,1)
    # psd_ax.set_yscale('log')

    pd_ax = psd_ax.twinx()
    pd_ax.plot(freq1, degree, 'r', alpha=0.5)
    pd_ax.set_ylim(0,1)
    pd_ax.tick_params(axis='y', colors='red')
    pd_ax.yaxis.label.set_color('red')
    pd_ax.set_ylabel('PD')

    psd_ax.set_xlim(fq_band)

    if return_fig:
        return fig, [[tr_ax, (psd_ax, pd_ax), spec_ax], [(psd, freq0), (degree, freq1)], starttime]

    if return_axes:
        return [[tr_ax, (psd_ax, pd_ax), spec_ax], [(psd, freq0), (degree, freq1)], starttime]

    else:
        plt.show()


def plot_event_gui(lde, **kwargs):
    app = QtWidgets.QApplication(sys.argv)
    lde_window = LDEWindow(lde, **kwargs)
    lde_window.show()
    sys.exit(app.exec_())


def __get_color__(event_label):

    if event_label == 'LPS':
        return 'sienna'
    
    if event_label == 'BT':
        return 'orangered'

    elif event_label == 'MPT':
        return 'teal'

    elif event_label == 'SPT':
        return 'navy'

    if event_label == 'BBT':
        return 'olive'

    elif event_label == 'ST':
        return 'purple'

    return 'firebrick'


def __get_marker__(event_label):
    if event_label[0] in ('C', 'S'):
        if event_label[-1] == 'B':
            return 'o'

        elif event_label[-1] == 'P':
            return 'v'

        elif event_label[-1] == 'M':
            return '*'

    else:
        return '^'


def _get_list(events, attr, r=False, gaussian=False):
    list_values = []
    for evnt in events:
        dout = evnt.get_values(gaussian=gaussian, r=r)
        list_values += [dout[attr]]

    return list_values


def _get_axis_label(attr):
    if attr == 'energy':
        return r'$\overline{e}$'+'\n[dB]'

    elif attr == 'pentropy':
        return r'$\overline{h}$'

    elif attr == 'fq_dominant':
        return r'$\overline{f_d}$'+'\n[Hz]'

    elif attr == 'fq_centroid':
        return r'$\overline{f_c}$'+'\n[Hz]'

    elif attr == 'degree_max':
        return r'$\overline{p_d}$'

    elif attr == 'degree_wavg':
        return r'$\overline{p_{avg}}$'

    elif attr == 'fq_polar':
        return r'$\overline{f_p}$'+'\n[Hz]'

    elif attr == 'r':
        return r'$r_{eh}$'

    else:
        return None


def _get_points(points, rect_min):
    freq = []
    rect = []
    azm = []
    dip = []
    spectrum = []
    for (f, r, s, d, a) in points:
        if r >= rect_min:
            freq += [f]
            rect += [r]
            azm += [a]
            dip += [d]
            spectrum += [s]
        # else:
        #     freq += [f]
        #     rect += [r]
        #     azm += [a]
        #     dip += [d]
        #     pentropy += [s]

    return np.array(freq), np.array(rect), np.array(azm), np.array(dip), np.array(spectrum)


def plot_activity(lde, activity=None, eruptive_activity=[], attrs=['e', 'pe', 'fd', 'fc', 'pm', 'fd'], 
    show=True, **kwargs):
    figsize = kwargs.get('figsize', (18,9))
    time_delta = kwargs.get('time_delta', 20)
    time_interval = kwargs.get('time_interval', ())
    e_mean = kwargs.get('e_mean', None)
    pe_mean = kwargs.get('pe_mean', None)


    def _get_axes(axes, activity, attr):
        i = attrs.index(attr)
        if activity:
            i += 3
        return axes[i]


    def _get_values(attr):
        if attr == 'e':
            return (5, 0.13)

        if attr == 'pe':
            return (0.05, 0.115)

        if attr in ('fd', 'fc', 'fp'):
            return (0.5, 0.1)

        if attr == 'pm':
            return (0.1, 0.115)

        return (0,0)


    fig = plt.figure(figsize=figsize)
    ncol = 1

    if activity:
        nrows = len(attrs) + 2
        grid_opt = {'height_ratios':[1,6]+[10]*5}
    
    else:
        nrows = len(attrs)
        grid_opt = None

    axes = fig.subplots(nrows, ncol, gridspec_kw=grid_opt)

    if activity:
        ax_light = axes[0]
        ax_activity = axes[1]
        ax_seis = axes[2]

    evnts = lde.get_events(time_interval=time_interval)
    
    if not time_interval:
        starttime = min([e.starttime for e in evnts]).date()
        endtime = max([e.starttime for e in evnts]).date()

    else:
        (starttime, endtime) = time_interval

    starttime = starttime.replace(day=1)
    endtime = endtime.replace(day=calendar.monthrange(endtime.year, endtime.month)[1])
    time_delta_in_sec = 60*time_delta
    npts = int((endtime-starttime).total_seconds()/time_delta_in_sec)
    time = [dt.datetime.fromordinal(starttime.toordinal())\
    + dt.timedelta(minutes=time_delta*x) for x in range(npts+1)]

    itertime = starttime
    months_years = []
    while itertime <= endtime:
        txt = (itertime.year, itertime.month)
        if txt not in months_years:
            months_years += [txt]
        itertime += dt.timedelta(days=1)

    # set major ticks
    major_text = [str(calendar.month_abbr[x[1]])+str(str(x[0])[-2:]) for x in months_years]
    list_months = [dt.datetime(x[0],x[1],13) for x in months_years]
    list_months.sort()
    major_locator = mdates.MonthLocator()
    minor_locator = mdates.DayLocator(bymonthday=None, interval=7)
    minor_formatt = mdates.DateFormatter('%d')

    # init handle list
    handles = []

    for n, attr in enumerate(attrs):
        axis = _get_axes(axes, activity, attr)
        vlim, prc_val = _get_values(attr)

        if attr == 'e':
            mean = e_mean

        elif attr == 'pe':
            mean = pe_mean

        else:
            mean = None

        if n == len(attrs)-1:
            minor = True

        else:
            minor= False

        values = _get_list(evnts, attr)
        v_min, v_max, v_mean = np.nanmin(values), np.nanmax(values), np.nanmean(values)
        v_lims = (v_min-vlim, v_max+vlim)
        v_diff = v_max-v_min

        for e, val  in zip(evnts, values):
            x_min = time[np.argmin(np.abs(np.array(time)-e.starttime))]
            x_width = dt.timedelta(hours=e.duration)
            x_ranges = [(x_min, x_width)]
            y_ranges = (val-(prc_val/2)*v_diff, prc_val*v_diff)
            color = __get_color__(e.label)
            h = axis.broken_barh(x_ranges, y_ranges, facecolors=color)
            if not handles or e.label not in [x1[1] for x1 in handles]:
                handles += [(h, e.label)]

        if mean:
            axis.axhline(y=mean, color=MEAN_COLOR, lw=0.5, ls='-.', alpha=0.5)

        axis.set_ylim(v_lims)
        axis.set_xlim(time[0], time[-1])
        axis.set_ylabel(_get_axis_label(attr))
        axis.xaxis.set_major_locator(major_locator)
        axis.xaxis.set_minor_locator(minor_locator)
        axis.xaxis.set_minor_formatter(mtick.NullFormatter())
        axis.xaxis.set_major_formatter(mtick.NullFormatter())

        if attr == 'e':
            fmtt = '%i'
        else:
            fmtt = '%.1f'

        axis.yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        axis.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
        axis.yaxis.set_major_formatter(mtick.FormatStrFormatter(fmtt))

        axis.grid(axis='x', which='major', color='k', linestyle='-', alpha=0.5)
        axis.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.3)

        if minor:
            axis.xaxis.set_minor_formatter(minor_formatt)
            for txt, tick in zip(major_text, list_months):
                axis.text(tick, v_min-3, s=txt, fontsize=10)

    if eruptive_activity:
        for n, ax in enumerate(axes):
            for erup in eruptive_activity:
                ax.axvspan(erup[0], erup[1], alpha=0.3, color=erup[2], zorder=0)

    if activity:
        if time_interval:
            info_activity = activity.get_activity_dict(starttime=time_interval[0],
                endtime=time_interval[1])
        else:
            info_activity = activity.get_activity_dict()

        if info_activity['explosion']:
            for (st, et) in info_activity['explosion']:
                st = dt.datetime(st.year, st.month, st.day, 12, 00, 00)
                st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                
                if not et:
                    h = ax_activity.scatter(st_bin, 1.5, s=9**2, marker='*', color=EXP_COLOR, zorder=5)
                
                else:
                    day_count = (et - st).days + 1
                    for day in (st + dt.timedelta(n) for n in range(day_count)):
                        st = dt.datetime(day.year, day.month, day.day, 12, 00, 00)
                        st_bin = time[np.argmin(np.abs(np.array(self.time)-st))]
                        h = ax_activity.scatter(st_bin, 1.5, s=9**2, marker='*', color=EXP_COLOR, zorder=5)
                
                if not handles or 'Explosion' not in [x1[1] for x1 in handles]:
                    handles += [(h, 'Explosion')]

        if info_activity['geyser']:
            for (st, et) in info_activity['geyser']:
                st = dt.datetime(st.year, st.month, st.day, 12, 00, 00)
                st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                if not et:
                    h = ax_activity.scatter(st_bin, 1.5, s=9**2, marker='*', color=GEYSER_COLOR, zorder=5)
                else:
                    day_count = (et - st).days + 1
                    for day in (st + dt.timedelta(n) for n in range(day_count)):
                        st = dt.datetime(day.year, day.month, day.day, 12, 00, 00)
                        st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                        h = ax_activity.scatter(st_bin, 1.5, s=9**2, marker='+', color=GEYSER_COLOR, zorder=5)
                
                if not handles or 'Geyser' not in [x1[1] for x1 in handles]:
                    handles += [(h, 'Geyser')]

        if info_activity['so2_anomaly']:
            for (st, et) in info_activity['so2_anomaly']:
                st = dt.datetime(st.year, st.month, st.day, 12, 00, 00)
                st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                
                if not et:
                    h = ax_activity.scatter(st_bin, 1, marker='o', color=SO2_COLOR, zorder=4)
                
                else:
                    day_count = (et - st).days + 1
                    for day in (st + dt.timedelta(n) for n in range(day_count)):
                        st = dt.datetime(day.year, day.month, day.day, 12, 00, 00)
                        st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                        h = ax_activity.scatter(st_bin, 1, marker='o', color=SO2_COLOR, zorder=4)
                
                if not handles or S02_label not in [x1[1] for x1 in handles]:
                    handles += [(h, S02_label)]

        if info_activity['temp_anomaly']:
            for (st, et) in info_activity['temp_anomaly']:
                st = dt.datetime(st.year, st.month, st.day, 12, 00, 00)
                st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                
                if not et:
                    h = ax_activity.scatter(st_bin, 0.5, s=5**2, marker='^', color=T_COLOR, zorder=3)
                
                else:
                    day_count = (et - st).days + 1
                    for day in (st + dt.timedelta(n) for n in range(day_count)):
                        st = dt.datetime(day.year, day.month, day.day, 12, 00, 00)
                        st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                        h = ax_activity.scatter(st_bin, 0.5, s=5**2, marker='^', color=T_COLOR, zorder=3)
                
                if not handles or TA_label not in [x1[1] for x1 in handles]:
                    handles += [(h, TA_label)]

        if info_activity['incandescence']:
            for (st, et) in info_activity['incandescence']:
                st = dt.datetime(st.year, st.month, st.day, 12, 00, 00)
                st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                
                if not et:
                    h = ax_activity.scatter(st_bin, 0.75, s=5**2, marker='s', color=IC_COLOR, zorder=3)
                
                else:
                    day_count = (et - st).days + 1
                    for day in (st + dt.timedelta(n) for n in range(day_count)):
                        st = dt.datetime(day.year, day.month, day.day, 12, 00, 00)
                        st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                        h = ax_activity.scatter(st_bin, 0.75, s=5**2, marker='s', color=IC_COLOR, zorder=3)
                
                if not handles or IC_label not in [x1[1] for x1 in handles]:
                    handles += [(h, IC_label)]

        if info_activity['ash']:
            for (st, et) in info_activity['ash']:
                if et:
                    et = et + dt.timedelta(days=1) - dt.timedelta(minutes=time_delta)
                
                else:
                    delta = dt.timedelta(days=1) - dt.timedelta(minutes=time_delta)
                    et = st + delta
                st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                et_bin = time[np.argmin(np.abs(np.array(time)-et))]
                h = ax_activity.axvspan(st_bin, et_bin, alpha=0.7, color=ASH_COLOR, zorder=2)
                
                if not handles or 'Ash' not in [x1[1] for x1 in handles]:
                    handles += [(h, 'Ash')]

        if info_activity['gas']:
            for (st, et) in info_activity['gas']:
                if et:
                    et = et + dt.timedelta(days=1) - dt.timedelta(minutes=time_delta)
                
                else:
                    delta = dt.timedelta(days=1) - dt.timedelta(minutes=time_delta)
                    et = st + delta
                
                st_bin = time[np.argmin(np.abs(np.array(time)-st))]
                et_bin = time[np.argmin(np.abs(np.array(time)-et))]
                h = ax_activity.axvspan(st_bin, et_bin, alpha=0.2, color=GAS_COLOR, zorder=1)
                
                if not handles or 'Gas' not in [x1[1] for x1 in handles]:
                    handles += [(h, 'Gas')]

        ax_activity.set_xlim(time[0], time[-1])
        ax_activity.xaxis.set_major_locator(major_locator)
        ax_activity.xaxis.set_minor_locator(minor_locator)
        ax_activity.xaxis.set_minor_formatter(mtick.NullFormatter())
        ax_activity.xaxis.set_major_formatter(mtick.NullFormatter())
        ax_activity.yaxis.set_major_locator(mtick.NullLocator())
        ax_activity.grid(axis='x', which='major', color='k', linestyle='-', alpha=0.5)
        ax_activity.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.3)
        ax_activity.set_ylim(0.2, 1.8)

        for (st, et, c) in info_activity['color_light']:
            if et:
                et = et + dt.timedelta(days=1) - dt.timedelta(minutes=time_delta)
            else:
                delta = dt.timedelta(days=1) - dt.timedelta(minutes=time_delta)
                et = st + delta
            st_bin = time[np.argmin(np.abs(np.array(time)-st))]
            et_bin = time[np.argmin(np.abs(np.array(time)-et))]
            width = et_bin - st_bin

            if c == 'o':
                c = 'orange'

            if c == 'y':
                c = 'yellow'

            ax_light.broken_barh([(st_bin, width)], (0,1), alpha=1, color=c)

        ax_light.set_xlim(time[0], time[-1])
        ax_light.xaxis.set_major_locator(major_locator)
        ax_light.xaxis.set_minor_locator(minor_locator)
        ax_light.xaxis.set_minor_formatter(mtick.NullFormatter())
        ax_light.xaxis.set_major_formatter(mtick.NullFormatter())
        ax_light.yaxis.set_major_locator(mtick.NullLocator())
        ax_light.grid(axis='x', which='major', color='k', linestyle='-', alpha=0.5)
        sax_light.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.3)
        ax_light.set_ylim(0, 1)

        max_val = 0
        min_val = np.inf
        for (st, et, nro) in info_activity['seismicity']:

            if et:
                et = et + dt.timedelta(days=1) - dt.timedelta(minutes=time_delta)
                delta = (et-st).days
                nro = nro/delta

            else:
                delta = dt.timedelta(days=1)
                et = st + delta

            st_bin = time[np.argmin(np.abs(np.array(time)-st))]
            et_bin = time[np.argmin(np.abs(np.array(time)-et))]
            width = et_bin - st_bin
            ax_seis.broken_barh([(st_bin, width)], (0,nro), facecolor='royalblue', edgecolor='None', alpha=0.7)
            
            if nro > max_val:
                max_val = nro
            if nro < min_val:
                min_val = nro

        ax_seis.set_xlim(time[0], time[-1])
        # ax_seis.set_yscale('log')
        # ax_seis.set_ylim(1, max_val)
        ax_seis.set_ylabel('events/day')
        #ax_seis.yaxis.set_major_locator(mtick.LogLocator(base=10.0, subs=None))
        ax_seis.xaxis.set_major_locator(major_locator)
        ax_seis.xaxis.set_minor_locator(minor_locator)
        ax_seis.xaxis.set_major_formatter(mtick.NullFormatter())
        ax_seis.xaxis.set_minor_formatter(mtick.NullFormatter())
        ax_seis.grid(axis='x', which='major', color='k', linestyle='-', alpha=0.5)
        ax_seis.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.3)

        if info_activity['so2_doas']:
            ax_so2doas = ax_seis.twinx()
            activity['so2_doas'].sort(key=lambda x: x[0])
            time = []
            so2 = []
            for (st, et, nro) in activity['so2_doas']:
                time += [time[np.argmin(np.abs(np.array(time)-st))],
                time[np.argmin(np.abs(np.array(time)-et))]]
                so2 += [nro, nro]

            ax_so2doas.plot(time, so2, 'k-')
            ax_so2doas.set_ylabel(r'S0$_2$'+'\nTon/day')
            ax_so2doas.ticklabel_format(axis='y',useMathText=True)
            ax_so2doas.xaxis.set_major_formatter(mtick.NullFormatter())
            ax_so2doas.xaxis.set_minor_formatter(mtick.NullFormatter())


    # plot the legend
    fig.legend(handles=[x1[0] for x1 in handles], labels=[x1[1] for x1 in handles],
        bbox_to_anchor=(0.70, 0.95), ncol=len(handles))
    fig.align_ylabels()

    if show:
        plt.show()

    else:
        return fig


def plot_scatter(lde, time_interval=(), y_axis=['pentropy', 'degree_max'], x_axis=['energy', 'fq_centroid'], show=True, **kwargs):
    figsize = kwargs.get('figsize', (12,12))
    fontsize = kwargs.get('fontsize', 10)
    events = kwargs.get('events', True)
    gaussian = kwargs.get('gaussian', False)
    umbral = kwargs.get('umbral', None)
    limits = kwargs.get('limits', None)
    umbral_color = kwargs.get('umbral_color', 'red')

    if events:
        events = lde.get_events(time_interval=time_interval)
        detections = overdetections = None

    else:
        detections = kwargs.get('detections', None)
        overdetections = kwargs.get('overdetections', None)

        if not detections and not overdetections:
            print('No events to plot')
            exit()

    if 'r' in y_axis or x_axis:
        r = True
    else:
        r = False
    
    # check y_axis and x_axis are ok!
    fig = plt.figure(figsize=figsize)
    nrows, ncol = len(x_axis), len(y_axis)
    axes = fig.subplots(nrows, ncol)

    if nrows + ncol == 2:
        axes = [axes]

    else:
        axes = axes.reshape(nrows*ncol,)

    pe_mean = kwargs.get('pe_mean', None)

    def get_d(attr):
        if attr == 'energy': 
            return 0.25
        if attr == 'pentropy': 
            return 0.0005
        if attr in ('fq_centroid', 'fq_dominant', 'fq_polar'): 
            return 0.025
        if attr in ('degree_max', 'degree_wavg', 'r'):
            return 0.005

    # combination of plots
    plots = []
    handles = []

    if events:
        labels = [e.label for e in events]
        ids = [e.id for e in events]
    else:
        labels = []
        ids =[]

    n = 0
    for y in y_axis:
        for x in x_axis:
            ax = axes[n]

            if events:
                valx = _get_list(events, x, r=r, gaussian=gaussian)
                valy = _get_list(events, y, r=r, gaussian=gaussian)
            else:
                valx = []
                valy = []
            
            if overdetections:
                for ov_id, ov_dict in overdetections.items():
                    valx += [ov_dict[x]]
                    valy += [ov_dict[y]]
                    labels += ['OV']
                    ids += [ov_id[:2]+r'$_{%i}$' % (int(ov_id[-1])+1)]

            if detections:
                for dt_id, dt_dict in detections.items():
                    valx += [dt_dict[x]]
                    valy += [dt_dict[y]]
                    labels += ['DT']
                    ids += [dt_id[:2]+r'$_{%i}$'% (int(dt_id[-1])+1)]

            dy = get_d(y)
            dx = get_d(x)

            for lbl, e_id, vx, vy  in zip(labels, ids, valx, valy):
                h = ax.scatter(vx, vy, s=40, marker=__get_marker__(lbl), c=__get_color__(lbl))
                ax.text(vx+dx, vy+dy, s=e_id, fontsize=fontsize, c=__get_color__(lbl))

                if lbl == 'OV':
                    lbl = 'Overdetection'

                if lbl == 'DT':
                    lbl = 'Detection'

                ax.grid(axis='both', which='major', color='k', linestyle='--', alpha=0.3)
                if not handles or lbl not in [x1[1] for x1 in handles]:
                    handles += [(h, lbl)]

            y_label = _get_axis_label(y)
            ax.set_ylabel(y_label)
            x_label = _get_axis_label(x)
            ax.set_xlabel(x_label)
            ax.xaxis.set_major_locator(mtick.LinearLocator(numticks=5))
            ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
            ax.yaxis.set_major_locator(mtick.LinearLocator(numticks=5))
            ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))

            if umbral:
                if y in umbral.keys():
                    h = ax.axhline(y=umbral[y], color=umbral_color, lw=0.5, ls='--', alpha=0.5)
                    if not handles or 'Threshold' not in [x1[1] for x1 in handles]:
                        handles += [(h, 'Threshold')]

                if x in umbral.keys():
                    h = ax.axvline(x=umbral[x], color=umbral_color, lw=0.5, ls='--', alpha=0.5)
                    if not handles or 'Threshold' not in [x1[1] for x1 in handles]:
                        handles += [(h, 'Threshold')]

            if limits:
                if y in limits.keys():
                    ax.set_ylim(limits[y])

                if x in limits.keys():
                    ax.set_xlim(limits[x])

            n += 1



    # e_types = tuple(collections.Counter([e.label for e in events]))
    ncol = len(handles)
    # if 'threshold' in [x1[1] for x1 in handles]:
    #     ncol -= 1

    fig.legend(handles=[x1[0] for x1 in handles], labels=[x1[1] for x1 in handles],
        loc=1, ncol=1)

    if show:
        plt.show()

    return fig


class LDEpolarWindow(QtWidgets.QMainWindow):
    def __init__(self, event, **kwargs):
        QtWidgets.QMainWindow.__init__(self)
        self.ui = event_polar.Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle("Event %s (%s) - Polar attribute" % (event.id, event.label))
        self.canvas = LDEpolar(event, parent=self, **kwargs)
        self.ui.verticalLayout.addWidget(self.canvas)

        self.ui.doubleSpinBox_3.setMaximum(self.ui.doubleSpinBox_4.value())
        self.ui.doubleSpinBox_4.setMinimum(self.ui.doubleSpinBox_3.value())

        # change the criteria for define rectiliniarity values
        self.ui.doubleSpinBox_2.valueChanged.connect(self.canvas._replot_rectiliniarity)
        
        self.ui.doubleSpinBox.valueChanged.connect(self.canvas._replot_polar_attributes)
        self.ui.doubleSpinBox_3.valueChanged.connect(self.canvas._replot_polar_attributes)
        self.ui.doubleSpinBox_4.valueChanged.connect(self.canvas._replot_polar_attributes)


class LDEpolar(FigureCanvas):
    def __init__(self, event, parent=None, **kwargs):
        self.fig = Figure(subplotpars=SubplotParams(left=0.1, right=0.9, top=0.95, bottom=0.05, hspace=0.5, wspace=0.1))
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(parent)

        self.parent = parent
        self.event = event
        self.lte = event.get_lte_sup()
        self.start = self.lte.stats.starttime
        self.end = self.lte.stats.endtime

        self.polar_min = kwargs.get('polar_min', self.lte.stats.threshold)
        self.parent.ui.doubleSpinBox_2.setValue(self.polar_min)
        self.parent.ui.doubleSpinBox_2.setMinimum(self.lte.stats.threshold)
        self.parent.ui.doubleSpinBox_2.setMaximum(1)
        
        self.freq_min, self.freq_max = self.lte.stats.fq_band
        self.parent.ui.doubleSpinBox_3.setValue(self.lte.stats.fq_band[0])
        self.parent.ui.doubleSpinBox_3.setMinimum(self.lte.stats.fq_band[0])
        self.parent.ui.doubleSpinBox_4.setValue(self.lte.stats.fq_band[1])
        self.parent.ui.doubleSpinBox_4.setMaximum(self.lte.stats.fq_band[1])
        
        self.rect_min = kwargs.get('rect_min', 0.75)
        self.parent.ui.doubleSpinBox.setValue(self.rect_min)

        self.rect_cmap = kwargs.get('rect_cmap', 'inferno')
        self.spec_min = kwargs.get('spec_min', None)
        self.spec_max = kwargs.get('spec_max', None)

        self.setParent(parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self.rect_hline = None
        self.clicks = dict(right=dict(tick=[], time=None), left=dict(tick=[], time=None))
        self.plot()
        
        self.parent.ui.pushButton.clicked.connect(self.reset_clicks)
        self.parent.ui.pushButton_2.clicked.connect(self.save_figure)
        self.callbacks.connect('button_press_event', self.on_click)


    def _set_axes(self):
        gs0 = self.fig.add_gridspec(3, 3, height_ratios=[1, 1, 2], width_ratios=[1, 1, 0.015])
        self.polargram_ax = self.fig.add_subplot(gs0[0,:2])
        self.polargram_bar_ax = self.fig.add_subplot(gs0[0,-1])
        self.rect_ax = self.fig.add_subplot(gs0[1,:2]) # first plot: frequency--rectiliniarity--psd
        self.rect_bar_ax = self.fig.add_subplot(gs0[1:,-1]) # first plot: colorb bar
        self.azm_ax = self.fig.add_subplot(gs0[2,0], projection='polar')
        self.dip_ax = self.fig.add_subplot(gs0[2,1], projection='polar')


    def _set_values(self):
        self.rect_min = self.parent.ui.doubleSpinBox.value()
        self.polar_min = self.parent.ui.doubleSpinBox_2.value()
        self.freq_min = self.parent.ui.doubleSpinBox_3.value()
        self.freq_max = self.parent.ui.doubleSpinBox_4.value()
        self.parent.ui.doubleSpinBox_3.setMaximum(self.freq_max)
        self.parent.ui.doubleSpinBox_4.setMinimum(self.freq_min)


    def _set_points(self):
        self.points = []
        specgram, freq = self.lte.get_attr('specgram', self.start, self.end)
        color = 10*np.log10(specgram)

        degree, freq = self.lte.get_attr('degree', self.start, self.end)
        rect, freq = self.lte.get_attr('rect', self.start, self.end)
        dip, freq = self.lte.get_attr('dip', self.start, self.end)
        azm, freq = self.lte.get_attr('azm', self.start, self.end)

        freq_k = [] # frequency
        rect_k = [] # rectiliniarity
        dip_k = [] # dip
        azm_k = [] # azm
        color_k = [] # energy

        for fq in range(rect.shape[1]):
            for tm in range(rect.shape[0]):
                if np.isfinite(rect[tm, fq]):
                    freq_k = freq[fq]
                    d = degree[tm, fq]
                    if self.freq_min <= freq_k <= self.freq_max and d >= self.polar_min:
                        rect_k = rect[tm, fq]
                        dip_k = dip[tm, fq]
                        azm_k = azm[tm, fq]
                        if len(color.shape) > 1:
                            color_k = color[tm, fq]
                        else:
                            color_k = color[tm]
                        self.points += [(freq_k, rect_k, color_k, dip_k, azm_k)]


    def _set_polargram(self):
        self.polargram_ax.cla()
        self.polargram_bar_ax.cla()
        
        data = self.lte.get_attr('degree')
        zdata = data[0].T
        ydata = data[1]
        xdata = self.time

        zdata = np.flipud(zdata)
        cmap = plt.get_cmap('inferno')
        cmap.set_bad(color='azure', alpha = 0.5)
        norm, _ =  get_norm(cmap, zdata, v_min=0, v_max=1)

        halfbin_time = (mdates.date2num(xdata[1]) - mdates.date2num(xdata[0])) / 2.0
        halfbin_freq = (ydata[1] - ydata[0]) / 2.0
        extent = (mdates.date2num(xdata[0]) + halfbin_time,
                  mdates.date2num(xdata[-1]) - halfbin_time,
                  ydata[0] + halfbin_freq,
                  ydata[-1] - halfbin_freq)

        interpolation = 'gaussian'
        im = self.polargram_ax.imshow(zdata, cmap=cmap, norm=norm, interpolation=interpolation, extent=extent, aspect='auto')
        self.polargram_ax.xaxis_date()
        self.polargram_ax.axis('tight')
        self.polargram_ax.set_xlim(self.time[0], self.time[-1])
        self.polargram_ax.set_ylim(self.lte.stats.fq_band)
        self.polargram_ax.set_xlabel('Time (UTC)')
        self.polargram_ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M \n %d %b'))

        self.freq_majorticks = mtick.MultipleLocator(1)
        self.freq_minorticks = mtick.AutoMinorLocator(2)

        self.polargram_ax.yaxis.set_major_locator(self.freq_majorticks)
        self.polargram_ax.yaxis.set_minor_locator(self.freq_minorticks)

        # add colorbar
        cbar = self.fig.colorbar(im, cax=self.polargram_bar_ax, orientation='vertical')
        cbar.locator = mtick.MaxNLocator(nbins=4)
        cbar.set_label('Polarization\nDegree')
        cbar.update_ticks()


    def _set_rectiliniarity(self):
        self.rect_ax.cla()
        self.rect_bar_ax.cla()

        self._set_points()
        x_freq, y_rect, _, _, spectrum = _get_points(self.points, 0)
        # status bar add ('num points (total) >>>>', len(points))

        if not self.spec_min:
            self.spec_min = spectrum.min()
        
        if not self.spec_max:
            self.spec_max = spectrum.max()

        cmap = plt.get_cmap(self.rect_cmap)
        levels = mtick.MaxNLocator(nbins=500).tick_values(self.spec_min, self.spec_max)
        self.norm = mcolor.BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        im = self.rect_ax.scatter(x_freq, y_rect, c=spectrum, norm=self.norm, cmap=self.rect_cmap, alpha=0.75)
        self.rect_ax.set_xlabel('Frequency [Hz]')
        self.rect_ax.set_ylabel('Rectiliniarity')
        self.rect_ax.set_xlim(self.lte.stats.fq_band)
        self.rect_ax.xaxis.set_major_locator(self.freq_majorticks)
        self.rect_ax.xaxis.set_minor_locator(self.freq_minorticks)
        self.rect_ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        self.rect_ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))

        cbar = self.fig.colorbar(im, cax=self.rect_bar_ax, orientation='vertical')
        cbar.locator = mtick.MaxNLocator(nbins=4)
        cbar.set_label('PSD')
        cbar.update_ticks()


    def _set_polar_attributes(self):
        self.azm_ax.cla()
        self.dip_ax.cla()

        self._set_hline()

        x_freq, y_rect, y_azm, y_dip, spectrum = _get_points(self.points, self.rect_min)
        self.parent.ui.statusbar.showMessage('Total points: %i' % len(x_freq), 0)

        im_a = self.azm_ax.scatter(y_azm*np.pi/180, x_freq, c=spectrum, s=(10*y_rect)**2, norm=self.norm, cmap=self.rect_cmap, alpha=0.75, linewidths=0.3, edgecolors='k')
        freq_min, freq_max = self.lte.stats.fq_band
        self.azm_ax.set_rlim(self.lte.stats.fq_band)
        self.azm_ax.set_rticks(np.arange(math.ceil(freq_min),int(freq_max),1))
        self.azm_ax.set_xticklabels(['0', '45', '90', '135', '180', '-135', '-90', '-45'])
        self.azm_ax.set_rorigin(-0.5)
        self.azm_ax.set_rlabel_position(40)
        self.azm_ax.set_theta_direction(+1)
        self.azm_ax.set_theta_zero_location('N', offset=0)
        self.azm_ax.tick_params(labelleft=True, labelright=False,
                        labeltop=False, labelbottom=True)
        trans, _ , _ = self.azm_ax.get_xaxis_text1_transform(-10)
        self.azm_ax.set_ylabel('Azimuth', fontsize=12)

        im_d = self.dip_ax.scatter((y_dip)*np.pi/180, x_freq, c=spectrum, s=(10*y_rect)**2, norm=self.norm, cmap=self.rect_cmap, alpha=0.75, linewidths=0.3, edgecolors='k')
        self.dip_ax.set_rorigin(-0.5)
        self.dip_ax.set_rlim(self.lte.stats.fq_band)
        self.dip_ax.set_rticks(np.arange(math.ceil(freq_min),int(freq_max),1))
        self.dip_ax.set_rlabel_position(90)
        self.dip_ax.set_ylabel('Dip', fontsize=12)
        self.dip_ax.set_theta_direction(-1)
        self.dip_ax.set_thetamin(0)
        self.dip_ax.set_thetamax(90)
        self.dip_ax.set_theta_zero_location('E', offset=0)
        self.dip_ax.tick_params(labelleft=True, labelright=False,
                        labeltop=False, labelbottom=True)
        trans, _ , _ = self.dip_ax.get_xaxis_text1_transform(-10)


    def plot(self):
        self._set_axes()

        self.time = self.lte.get_time()

        self._set_polargram()

        self._set_rectiliniarity()

        self._set_polar_attributes()


    def save_figure(self):
        file_format = self.parent.ui.comboBox.currentText()
        file_name = './' + self.event.lte_file_sup
        file_name = file_name.replace('.lte', '.'+file_format)

        if file_format == 'svg':
            transparent = True
        else:
            transparent = False

        self.fig.savefig(file_name, transparent=transparent)


    def _set_hline(self):
        if self.rect_hline:
            self.rect_hline.remove()
        self.rect_hline = self.rect_ax.axhline(y=self.rect_min, color='gray', linestyle='--')


    def _replot_polar_attributes(self):
        self._set_values()
        self._set_points()
        self._set_polar_attributes()
        self.draw()


    def _replot_rectiliniarity(self):
        self._set_values()
        self._set_rectiliniarity()
        self._set_polar_attributes()
        self.draw()


    def reset_clicks(self):
        if self.clicks['left']['tick']:
            for ac_ver_bar in self.clicks['left']['tick']:
                ac_ver_bar.remove()
            self.clicks['left']['tick'] = []
            self.clicks['left']['time'] = None

        if self.clicks['right']['tick']:
            for ac_ver_bar in self.clicks['right']['tick']:
                ac_ver_bar.remove()
            self.clicks['right']['tick'] = []
            self.clicks['right']['time'] = None

        self.start = self.lte.stats.starttime
        self.end = self.lte.stats.endtime
        self._replot_rectiliniarity()
        self._set_polar_attributes()
        self.draw()


    def on_click(self, event):
        if event.inaxes == self.polargram_ax:
            t = mdates.num2date(float(event.xdata))
            t = t.replace(tzinfo=None)
            t = self.time[np.argmin(np.abs(np.array(self.time)-t))]

            if event.button == 1:
                if self.clicks['left']['tick']:
                    for ac_ver_bar in self.clicks['left']['tick']:
                        ac_ver_bar.remove()
                    self.clicks['left']['tick'] = []

                self.clicks['left']['time'] = t
                tax = self.polargram_ax.axvline(t, color='r')
                self.clicks['left']['tick'].append(tax)

            if event.button == 3:
                if self.clicks['right']['tick']:
                    for ac_ver_bar in self.clicks['right']['tick']:
                        ac_ver_bar.remove()
                    self.clicks['right']['tick'] = []

                self.clicks['right']['time'] = t
                tax = self.polargram_ax.axvline(t, color='darkgreen')
                self.clicks['right']['tick'].append(tax)

            if event.dblclick:
                self.reset_clicks()


            if self.clicks['left']['time'] and self.clicks['right']['time']:
                times = (self.clicks['left']['time'], self.clicks['right']['time'])
                self.start = min(times)
                self.end = max(times)
                self._replot_rectiliniarity()
                self._set_polar_attributes()
            else:
                self.draw()


def plot_event_polar(event, **kwargs):
    app = QtWidgets.QApplication(sys.argv)
    lde_window = LDEpolarWindow(event, **kwargs)
    lde_window.show()
    sys.exit(app.exec_())