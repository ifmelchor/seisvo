#!/usr/bin/env python
# coding=utf-8

import os, sys
import numpy as np
import datetime as dt
import pyqtgraph

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.qt_compat import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvas

from seisvo import LDE, Network
from seisvo.gui import notify, get_norm, PSD_GUI
from seisvo.gui.gstation import StationWindow
from seisvo.gui.frames import ltegui2
from seisvo.gui.glde import LDEWindow, __get_color__


def sort_list(llist):
    sort_list = [None] * 10
    for item in llist:
        if item == 'energy':
            sort_list[0] = item
        if item == 'pentropy':
            sort_list[1] = item
        if item == 'fq_dominant':
            sort_list[2] = item
        if item == 'fq_centroid':
            sort_list[3] = item
        if item == 'specgram':
            sort_list[4] = item
        if item == 'degree_max':
            sort_list[5] = item
        if item == 'degree_wavg':
            sort_list[6] = item
        if item == 'fq_polar':
            sort_list[7] = item
        if item == 'degree':
            sort_list[8] = item
        if item == 'rect':
            sort_list[9] = item
    return list(filter(None.__ne__, sort_list))


def get_ticks(time_step, time, tint, day_off=1, yday=False):
    one_day = int((60/time_step)*24)
    one_day_bins = np.arange(0, len(time), one_day)
    major_ticks = one_day_bins[day_off::tint]
    minor_ticks = list(set(one_day_bins).difference(set(major_ticks)))
    minor_ticks.sort()
    if yday:
        major_ticks_pos = [str(time[x].timetuple().tm_yday) for x in major_ticks]
    else:
        major_ticks_pos = [str(time[x].day) for x in major_ticks]
    return major_ticks, minor_ticks, major_ticks_pos


class LTEWindow(QtWidgets.QMainWindow):
    def __init__(self, lte, starttime, endtime, interval, list_attr, lde_parent=None, **kwargs):
        QtWidgets.QMainWindow.__init__(self)
        self.ui = ltegui2.Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle('LTE: '+lte.stats.id)

        self.canvas = LTECanvas(lte, starttime, endtime, interval, list_attr, parent=self, **kwargs)
        self.ui.horizontalLayout.addWidget(self.canvas)

        self.ui.actionLTE_info.triggered.connect(self.canvas.show_info)
        self.ui.actionSave_Figure.triggered.connect(self.canvas.save_fig)

        self.ui.actionReset_Ticks.triggered.connect(self.canvas.reset_ticks)
        self.ui.buttonReset_Ticks.clicked.connect(self.canvas.reset_ticks)
        
        self.ui.actionOpen_Database.triggered.connect(self.canvas.lde.open)
        self.ui.actionOpen_Database.setShortcut('F1')

        self.ui.actionGo_to.triggered.connect(self.canvas.set_starttime)
        self.ui.buttonGo_to.clicked.connect(self.canvas.set_starttime)
        self.ui.actionGo_to.setShortcut('tab')

        self.ui.actionSet_interval.triggered.connect(self.canvas.set_interval)
        self.ui.buttonSet_interval.clicked.connect(self.canvas.set_interval)
        self.ui.actionSet_interval.setShortcut('shift+tab')

        self.ui.actionWaveform.triggered.connect(self.canvas.plot_seismogram)
        self.ui.buttonWaveform.clicked.connect(self.canvas.plot_seismogram)
        self.ui.actionWaveform.setShortcut('p')

        self.ui.actionPSD.triggered.connect(self.canvas.plot_psd)
        self.ui.buttonPSD.clicked.connect(self.canvas.plot_psd)
        self.ui.actionPSD.setShortcut('shift+p')

        self.ui.actionSave_Event.triggered.connect(self.canvas.write_event)
        self.ui.buttonSave_Event.clicked.connect(self.canvas.write_event)
        self.ui.actionSave_Event.setShortcut('W')

        self.ui.actionPlotEvent.triggered.connect(self.canvas.plot_event)
        self.ui.actionPlotEvent.setShortcut('ctrl+E')

        self.ui.actionRemove.triggered.connect(self.canvas.remove_event)
        self.ui.actionRemove.setShortcut('ctrl+R')

        self.ui.actionScatter.triggered.connect(lambda: self.canvas.plot_lde(scatter=True))
        self.ui.actionPlot.triggered.connect(self.canvas.plot_lde)
        
        self.ui.actionWaveform.setEnabled(False)
        self.ui.buttonWaveform.setEnabled(False)

        self.ui.actionPSD.setEnabled(False)
        self.ui.buttonPSD.setEnabled(False)

        self.ui.actionSave_Event.setEnabled(False)
        self.ui.buttonSave_Event.setEnabled(False)

        self.ui.buttonForward.clicked.connect(self.canvas.move_forward)
        self.ui.buttonForward.setShortcut('right')
        self.ui.buttonBackwards.clicked.connect(self.canvas.move_backwards)
        self.ui.buttonBackwards.setShortcut('left')

        self.ui.eventWidget.installEventFilter(self)
        self.ui.eventWidget.itemClicked.connect(self.canvas.event_selected)

        # if LDE gui is the parent, hide it and when close LTE, turn on.
        self.lde_parent = lde_parent
        if self.lde_parent:
            self.lde_parent.setEnabled(False)

    def eventFilter(self, source, event):
        if (event.type() == QtCore.QEvent.ContextMenu and
            source is self.ui.eventWidget):
            item = source.itemAt(event.pos())
            id_event = int(item.text().split(' ')[3])
            
            menu = QtWidgets.QMenu()
            
            infoButton = QtWidgets.QAction('See Info', self)
            infoButton.triggered.connect(lambda: self.canvas.show_eventinfo(i=id_event))
            menu.addAction(infoButton)

            plotButton = QtWidgets.QAction('Plot', self)
            plotButton.setStatusTip('Plot event')
            plotButton.triggered.connect(lambda: self.canvas.plot_event(i=id_event))
            menu.addAction(plotButton)

            psdButton = QtWidgets.QAction('PSD', self)
            psdButton.setStatusTip('Plot PSD')
            psdButton.triggered.connect(lambda: self.canvas.plot_event_psd(i=id_event))
            menu.addAction(psdButton)

            openButton = QtWidgets.QAction('Open', self)
            openButton.setStatusTip('Open in Station GUI')
            openButton.triggered.connect(lambda: self.canvas.plot_event(i=id_event, sta=True))
            menu.addAction(openButton)        

            relabelButton = QtWidgets.QAction('Relabel', self)
            relabelButton.setStatusTip('Relabel event')
            relabelButton.triggered.connect(lambda: self.canvas.relabel_event(i=id_event))
            menu.addAction(relabelButton)

            removeButton = QtWidgets.QAction('Remove', self)
            removeButton.setStatusTip('Remove event')
            removeButton.triggered.connect(lambda: self.canvas.remove_event(i=id_event))
            menu.addAction(removeButton)
            
            menu.exec_(event.globalPos())
                
            return True
        else:
            return False

    def closeEvent(self, event):
        if self.lde_parent:
            self.lde_parent.setEnabled(True)
            with pyqtgraph.BusyCursor():
                if self.lde_parent.canvas.lde.is_event(self.lde_parent.canvas.evnt.id):
                    self.lde_parent.canvas._set_event(self.lde_parent.canvas.evnt.id)

                if self.lde_parent.canvas.mode == 'steprange':
                    self.lde_parent.canvas.plot(default=True)
                
                else:
                    self.lde_parent.canvas.plot()


class LTECanvas(FigureCanvas):
    def __init__(self, lte, starttime, endtime, interval, list_attr, parent=None, **kwargs):
        self.lte = lte
        self.lde = LDE(self.lte.stats.id)
        
        self.delta = dt.timedelta(days=interval)
        self.interval = interval

        self.full_start = self.lte.stats.starttime
        self.starttime = starttime
        
        self.full_end = self.lte.stats.endtime
        self.endtime = endtime

        self.list_attr = list_attr

        self.parent = parent
        self.kwargs = kwargs
        self.psd_frame = None
        self.sta_frame = None
        self.lde_frame = None

        figsize = kwargs.get("figsize", (20,9))
        self.fig = Figure(figsize=figsize, dpi=100,
            subplotpars=SubplotParams(left=0.08, right=0.92, wspace=0.1, top=0.95, bottom=0.05))
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self.cliks = dict(right=dict(tick=[], time=None), left=dict(tick=[], time=None))
        self.callbacks.connect('button_press_event', self.on_click)

        self.plot()


    def set_statusbar(self, interval):
        text = '\t %s  --- %s  | Interval: %i' % (self.starttime.strftime('%d %B %Y'), self.endtime.strftime('%d %B %Y'), interval)
        self.parent.ui.statusbar.showMessage(text, 0)


    def show_info(self):
        lte_info = " LTE File: %s\n ID: %s\n Start Time: %s\n End Time: %s\n Time Step: %s\n" % (
            os.path.basename(self.lte.lte_file), self.lte.stats.id, self.lte.stats.starttime.strftime("%Y-%m-%d"),
            self.lte.stats.endtime.strftime("%Y-%m-%d"), self.lte.stats.time_step)

        lte_info += " Sampling Rate: %s\n Remove Response: %s\n Freq. Band: %s\n Avg. Step: %s\n" % (
            self.lte.stats.sampling_rate, self.lte.stats.remove_response, self.lte.stats.fq_band, self.lte.stats.avg_step)

        lte_info += " Polargram: %s\n Matrix Return: %s\n PE (d,t): (%s,%s)\n" % (self.lte.stats.polargram,
            self.lte.stats.matrix_return, self.lte.stats.p_order, self.lte.stats.tau)

        self.parent.ui.textTickInfo.setText(lte_info)


    def get_attr_values(self, starttime, endtime, r=False):
        """
        Get the mean values of the period data specified
        If end is not specified, take the value at time of start.
        """

        if r:
            attr_list = self.list_attr + ['r']
        else:
            attr_list = self.list_attr

        attr_mean_dict = {}
        for attr in attr_list:
            if not self.lte.is_matrix(attr):
                if endtime:
                    attr_mean_dict[attr] = self.lte.get_values(attr, starttime, endtime)[2]
                else:
                    attr_mean_dict[attr] = self.lte.get_attr(attr, starttime, None)

        return attr_mean_dict


    def plot(self):
        self.fig.clf()

        self.delta = dt.timedelta(days=self.interval)

        if self.starttime < self.full_start:
            self.starttime = self.full_start
            notify('seisvo', 'End of the LTE file', status='warn')

        self.endtime = self.starttime + self.delta

        interval = self.interval
        if self.endtime > self.full_end:
            self.endtime = self.full_end
            interval = (self.endtime - self.starttime).days
            notify('seisvo', 'End of the LTE file', status='warn')

        self.time = self.lte.get_time(self.starttime, self.endtime)

        _, self.axes = get_fig(self.lte, self.starttime, self.endtime, interval, self.list_attr, 
            fig=self.fig, settitle=False, **self.kwargs)

        self.events = self.lde.get_events(time_interval=(self.starttime, self.endtime))
        
        self.parent.ui.eventWidget.clear()
        
        self.__eventlist = {}
        if self.events:
            for event in self.events:
                self.__eventlist[event.id] = []

            for i, attr in enumerate(self.list_attr):
                self.add_events(i)
            
            self.add_info_events()

        with pyqtgraph.BusyCursor():
            self.draw()

        self.set_statusbar(interval)


    def add_info_events(self):
        for e in self.events:
            info_evnt = '  ID: %s | Label: %s' % (e.id, e.label)
            self.parent.ui.eventWidget.addItem(info_evnt)


    def add_events(self, i):
        for event in self.events:
            st = event.starttime
            et = st + dt.timedelta(hours=event.duration)
            color = __get_color__(event.label)
            code = event.id

            cond1 = st >= self.starttime and et <= self.endtime
            cond2 = self.endtime > st > self.starttime and et > self.endtime
            cond3 = self.starttime > st and self.starttime < et < self.endtime
            cond4 = st < self.starttime and self.endtime < et

            if cond1:
                width = dt.timedelta(hours=event.duration)

            if cond2:
                width = self.endtime - st

            if cond3:
                st = self.starttime
                width = et - st

            if cond4:
                st = self.starttime
                width = self.endtime - st

            if i == 0:
                mid_bin = st + dt.timedelta(seconds=width.total_seconds()/2)
                x_fraction_num = mdates.date2num(self.starttime) - mdates.date2num(mid_bin)
                x_fraction_den = mdates.date2num(self.starttime) - mdates.date2num(self.endtime)
                x_fraction = x_fraction_num/x_fraction_den
                txt = self.axes[i, 0].annotate(code, xy=(x_fraction, 1.05), xycoords='axes fraction', c=color, fontsize=6)
                self.__eventlist[code] += [txt]

            span = self.axes[i, 0].axvspan(st, st + width, alpha=0.2, color=color, label=code)
            self.__eventlist[code] += [span]


    def event_selected(self):
        # search for golds colors
        for key in self.__eventlist.keys():
            if self.__eventlist[key][0].get_color() == 'gold':
                for item in self.__eventlist[key]:
                    event = self.lde.get_event(key)
                    item.set_color(__get_color__(event.label))

        item = self.parent.ui.eventWidget.selectedItems()[0]
        event_id = int(item.text().split(' ')[3])

        for item in self.__eventlist[event_id]:
            item.set_color('gold')

        with pyqtgraph.BusyCursor():
            self.draw()


    def remove_event(self, i=None):
        if not i:
            i, ok = QtWidgets.QInputDialog.getInt(self, "Remove Event","Event ID:")
        else:
            msgBox = QtWidgets.QMessageBox()
            msgBox.setIcon(QtWidgets.QMessageBox.Information)
            msgBox.setText("Are you sure you want to remove the event (ID) %s?" % i)
            msgBox.setWindowTitle("Remove Event")
            msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel)
            returnValue = msgBox.exec()
            if returnValue == QtWidgets.QMessageBox.Ok:
                ok = True
            else:
                return

        if ok and self.lde.is_event(i):
            self.lde.remove(id=i)
            self.plot()


    def relabel_event(self, i):
        evnt = self.lde.get_event(id=i)
        current_label = evnt.label
        i, ok = QtWidgets.QInputDialog.getText(self, "Relabel Event","Event Label: ", text=current_label)

        if ok:
            new_label = str(i).upper()
            self.lde.relabel(evnt, new_label)
            with pyqtgraph.BusyCursor():
                self.plot()


    def plot_event(self, i=None, sta=False):
        if not i:
            i, ok = QtWidgets.QInputDialog.getInt(self, "Plot Event","Event ID:")
        
        else:
            ok = True

        if ok and self.lde.is_event(i):
            if sta:
                event = self.lde.get_event(id=i)
                sta = self.lde.station

                with pyqtgraph.BusyCursor():
                    if self.sta_frame:
                        self.sta_frame.close()
                        self.sta_frame = None
                    self.sta_frame = StationWindow(sta, event.starttime, self.lte.stats.id.split('.')[3], 
                        30, specgram=True, specgram_lim=(-20, 60), 
                        fq_band=self.lte.stats.fq_band, one_day=False)
                    self.sta_frame.show()

            else:
                with pyqtgraph.BusyCursor():
                    if self.lde_frame:
                        self.lde_frame.close()
                        self.lde_frame = None
                    self.lde_frame = LDEWindow(self.lde, event_id=i, lte_parent=self.parent, vmin=-20, vmax=60)
                    self.lde_frame.show()
                self.parent.setEnabled(False)
        else:
            notify('seisvo', 'Event ID not found', status='warn')


    def plot_event_psd(self, i):
        if self.psd_frame:
            self.psd_frame.close()
            self.psd_frame = None

        event = self.lde.get_event(id=i)
        start = event.starttime
        end = start + dt.timedelta(hours=event.duration)

        data = self.lte.get_attr('specgram', start, end)
        psdlist = data[0]

        data2 = self.lte.get_attr('degree', start, end)
        pdlist = data2[0]

        with pyqtgraph.BusyCursor():
            self.psd_frame = PSD_GUI(data[1], psdlist, pd_array=pdlist, 
                plot_prob=True, title='PSD/PD -- %s' % self.lte.get_station_id())
            self.psd_frame.show()


    def reset_ticks(self):
        if self.cliks['left']['tick']:
            for ac_ver_bar in self.cliks['left']['tick']:
                ac_ver_bar.remove()
            self.cliks['left']['tick'] = []
            self.cliks['left']['time'] = None

        if self.cliks['right']['tick']:
            for ac_ver_bar in self.cliks['right']['tick']:
                ac_ver_bar.remove()
            self.cliks['right']['tick'] = []
            self.cliks['right']['time'] = None

        self.parent.ui.actionWaveform.setEnabled(False)
        self.parent.ui.buttonWaveform.setEnabled(False)
        self.parent.ui.actionPSD.setEnabled(False)
        self.parent.ui.buttonPSD.setEnabled(False)
        self.parent.ui.actionSave_Event.setEnabled(False)
        self.parent.ui.buttonSave_Event.setEnabled(False)

        self.draw()
        self.show_infoticks()


    def get_infoticks(self, r=False):
        if self.cliks['left']['time']:
            red_tick = self.cliks['left']['time']
        else:
            red_tick = None

        if self.cliks['right']['time']:
            green_tick = self.cliks['right']['time']
        else:
            green_tick = None

        if not green_tick and not red_tick:
            return red_tick, green_tick, np.nan, [np.nan]*5

        if red_tick and green_tick:
            times = [red_tick, green_tick]
            hr_dur = (max(times)-min(times)).total_seconds()/3600.0
            dict_ans = self.get_attr_values(min(times), max(times), r=r)

        else:
            if red_tick:
                time = red_tick
            else:
                time = green_tick

            hr_dur = np.nan
            dict_ans = self.get_attr_values(time, None)

        return red_tick, green_tick, hr_dur, dict_ans


    def show_infoticks(self):
        red_time, green_time, hr_dur, dict_ans = self.get_infoticks(r=True)

        if red_time:
            red_time = red_time.strftime("%Y-%m-%d %H:%M")

        if green_time:
            green_time = green_time.strftime("%Y-%m-%d %H:%M")

        text_tick = " Red click : %s\n" % red_time
        text_tick += " Green click : %s\n" % green_time
        text_tick += " Time delta  : %2.1f hr\n" % hr_dur

        for key in dict_ans.keys():
            text_tick += " %s : %2.2f\n" % (key, dict_ans[key])

        self.parent.ui.textTickInfo.setText(text_tick)


    def show_eventinfo(self, i):
        evnt = self.lde.get_event(i)
        average_values = evnt.get_values(r=True)

        info_text = evnt.__str__()
        info_text += '\n -- Mean values ----- \n' 
        for key in average_values:
            info_text += '    %s : %2.2f\n' % (key, average_values[key])

        self.parent.ui.textTickInfo.setText(info_text)


    def move_forward(self):
        self.starttime += self.delta
        self.plot()


    def move_backwards(self):
        self.starttime -= self.delta
        self.plot()


    def set_interval(self):
        i, ok = QtWidgets.QInputDialog.getInt(self, "Set interval","Interval:", self.interval, 5, 150, 5)
        if ok and int(i) != self.interval:
            self.interval = i
            self.plot()


    def set_starttime(self):
        current_day = self.starttime.strftime("%Y%m%d")
        time, ok = QtWidgets.QInputDialog.getText(self, "Set Starttime", "Time (YYYYMMDD):", text=current_day)
        if ok and str(time) != current_day:
            try:
                time = dt.datetime.strptime(str(time), "%Y%m%d")
            except:
                time = dt.datetime.now()

            if self.full_start <= time < self.full_end:
                self.starttime = time
                self.plot()
            
            else:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("Time out of the bounds.")
                msg.setWindowTitle("Error")
                msg.exec_()


    def plot_seismogram(self):

        if self.sta_frame:
            self.sta_frame.close()
            self.sta_frame = None

        net = Network(self.lte.stats.id.split('.')[0])
        sta = net.get_sta(self.lte.stats.id.split('.')[1], self.lte.stats.id.split('.')[2])

        times = [self.cliks['right']['time'], self.cliks['left']['time']]
        starttime = min(times)
        endtime = max(times)
        delta = endtime - starttime
        delta_minutes = int((endtime - starttime).total_seconds()/60)

        if delta.days > 1:
            one_day = True
            specgram = False
            interval = 60 #min

        else:
            one_day = False
            specgram = True
            interval = None #min

        with pyqtgraph.BusyCursor():
            self.sta_frame = StationWindow(sta, starttime, self.lte.stats.id.split('.')[3], delta_minutes, specgram=specgram,
                specgram_lim=(-20, 60), fq_band=self.lte.stats.fq_band, one_day=one_day,
                interval=interval)
            self.sta_frame.show()


    def plot_psd(self):
        if self.psd_frame:
            self.psd_frame.close()
            self.psd_frame = None

        times = [self.cliks['right']['time'], self.cliks['left']['time']]
        starttime = min(times)
        endtime = max(times)
        data = self.lte.get_attr('specgram', starttime, endtime)
        psdlist = data[0]
        data2 = self.lte.get_attr('degree', starttime, endtime)
        pdlist = data2[0]

        with pyqtgraph.BusyCursor():
            self.psd_frame = PSD_GUI(data[1], psdlist, pd_array=pdlist, 
                plot_prob=True, title='PSD/PD -- %s' % self.lte.get_station_id())
            self.psd_frame.show()


    def write_event(self):
        red_time, green_time, hr_dur, _ = self.get_infoticks()
        if red_time and green_time:
            lbl, ok = QtWidgets.QInputDialog.getText(self, "Long Duration Event type", "Label:")
            if ok:
                time = min([red_time, green_time])
                dict_out = {}
                dict_out['starttime'] = time
                dict_out['duration'] = hr_dur
                dict_out['label'] = lbl.upper()
                dict_out['lte_file'] = os.path.join(self.lte.stats.id.split('.')[0], os.path.basename(self.lte.lte_file))
                self.lde.save_event(dict_out)
                self.plot()


    def save_fig(self):
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save F:xile", "./", "Images (*.png *.svg *.pdf)")
        if fileName:
            self.fig.savefig(fileName)


    def on_click(self, event):
        try:
            t = mdates.num2date(float(event.xdata))
            t = t.replace(tzinfo=None)
            t = self.time[np.argmin(np.abs(np.array(self.time)-t))]
        except:
            return

        if event.button == 1:
            if self.cliks['left']['tick']:
                for ac_ver_bar in self.cliks['left']['tick']:
                    ac_ver_bar.remove()
                self.cliks['left']['tick'] = []


        if event.button == 3:
            if self.cliks['right']['tick']:
                for ac_ver_bar in self.cliks['right']['tick']:
                    ac_ver_bar.remove()
                self.cliks['right']['tick'] = []


        for iax in range(self.axes.shape[0]):
            if event.button == 1:
                self.cliks['left']['time'] = t
                tax = self.axes[iax,0].axvline(t, color='r')
                self.cliks['left']['tick'].append(tax)

            if event.button == 3 and not event.dblclick:
                self.cliks['right']['time'] = t
                tax = self.axes[iax,0].axvline(t, color='darkgreen')
                self.cliks['right']['tick'].append(tax)

            if event.button == 3 and event.dblclick:
                self.cliks['right']['time'] = None


        if self.cliks['right']['time'] and self.cliks['left']['time']:
            self.parent.ui.actionWaveform.setEnabled(True)
            self.parent.ui.buttonWaveform.setEnabled(True)

            self.parent.ui.actionPSD.setEnabled(True)
            self.parent.ui.buttonPSD.setEnabled(True)

            self.parent.ui.actionSave_Event.setEnabled(True)
            self.parent.ui.buttonSave_Event.setEnabled(True)

        else:
            self.parent.ui.actionWaveform.setEnabled(False)
            self.parent.ui.buttonWaveform.setEnabled(False)

            self.parent.ui.actionPSD.setEnabled(False)
            self.parent.ui.buttonPSD.setEnabled(False)

            self.parent.ui.actionSave_Event.setEnabled(False)
            self.parent.ui.buttonSave_Event.setEnabled(False)


        self.draw()
        self.show_infoticks()


    def plot_lde(self, scatter=False):
        if scatter:
            fig = plt.figure(figsize=(12,8))
            self.lde.scatter(fig=fig)
        else:
            if not self.lde.activity:
                f_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select the file containing activity', "./", "Text (*.lde *.txt)")
                if os.path.isfile(str(f_path)):
                    try:
                        self.lde.load_activity(f_path)
                    except:
                        self.lde.activity = None
                        msg = QtWidgets.QMessageBox()
                        msg.setIcon(QtWidgets.QMessageBox.Critical)
                        msg.setText("Error reading file")
                        msg.setWindowTitle("Error")
                        msg.exec_()
            self.lde.plot()


def plot_gui(lte, starttime, endtime, interval, list_attr, **kwargs):
    app = QtWidgets.QApplication(sys.argv)
    lte_window = LTEWindow(lte, starttime, endtime, interval, list_attr, **kwargs)
    lte_window.show()
    sys.exit(app.exec_())





