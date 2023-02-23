#!/usr/bin/python3
# seisvo
'''

 This module plot station object

'''

from functools import partial
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas

import numpy as np
import datetime as dt
import pyqtgraph 
import os

from seisvo import SDE
from seisvo.utils import get_catalog
from seisvo.plotting import get_colors
from seisvo.plotting.gui import Navigation, PSD_GUI
from seisvo.plotting.gui.frames import MainGUI


SHORTCUT_COLORS = get_colors('tolb')
SHORTCUT_LABELS = {
    0:'REG',
    1:'VT',
    2:'LP',
    3:'VLP',
    4:'TR',
    5:'EXP',
    6:'UK'
}
PHASE_COLORS = {
    'p':get_colors('okabe')[3],
    's':get_colors('okabe')[0],
    'f':get_colors('okabe')[6]
    }



class StationWindow(QtWidgets.QMainWindow, MainGUI):
    def __init__(self, station, starttime, channel, delta, sde, **kwargs):
        super(StationWindow, self).__init__()
        self.setupUi(self)
        self.canvas = StationCanvas(station, starttime, channel, delta, sde, parent=self, **kwargs)
        
        self.plotPSD = self.canvas.plotPSDeid
        self.plotEID = self.canvas.plotEID
        self.removeEID = self.canvas.removeEID
        self.relabelEID = self.canvas.relabelEID

        self.gridLayout.addWidget(self.canvas, 1, 1, 1, 1)
        self.setWindowTitle(station.stats.id + ' [{}]'.format(self.canvas.sde.sql_path))
        self.actionOpen.triggered.connect(self.canvas.sde.open)

        # add actions of the windows menu
        # self.ui.forwardButton.clicked.connect(self.canvas.forward_step)
        # self.ui.backwardButton.clicked.connect(self.canvas.backward_step)
        # self.ui.gotoButton.clicked.connect(self.canvas.set_starttime)
        # self.ui.intervalButton.clicked.connect(self.canvas.set_interval)
        self.saveButton.clicked.connect(self.canvas.show_catalog)
        # self.ui.saveButton.setShortcut('S')
        # self.ui.psdButton.clicked.connect(self.canvas.plot_psd)
        # self.ui.exportButton.clicked.connect(self.canvas.export_to_mlf)

        # self.actionRemove_Response.setChecked(int(kwargs.get('remove_response', False)))
        # self.actionRemove_Response.triggered.connect(self.canvas.set_removeresponse)
        
        # self.actionSpectrogram.setChecked(int(kwargs.get('specgram', False)))
        # self.actionSpectrogram.triggered.connect(self.canvas.set_spectrogram)
        
        # self.actionPolargram.setChecked(int(kwargs.get('polargram', False)))
        # self.actionPolargram.triggered.connect(self.canvas.set_polargram)
        
        # self.action3component.setChecked(int(kwargs.get('three_component', False)))
        # self.action3component.triggered.connect(self.canvas.set_3component)
        
        # self.actionPreferences.triggered.connect(self.canvas.set_settings)
        # self.toolButton.clicked.connect(self.canvas.set_settings)
        self.listWidget.installEventFilter(self)
        # self.listWidget.itemClicked.connect(self.canvas.event_selected)


    def eventFilter(self, source, event):
        if (event.type() == QtCore.QEvent.ContextMenu and
            source is self.listWidget):
            item = source.itemAt(event.pos())
            id_event = int(item.text().split('[')[0].strip())
            
            menu = QtWidgets.QMenu()

            plotPSD_Button = QtWidgets.QAction('PSD', self)
            displayPSD = partial(self.plotPSD, id_event)
            plotPSD_Button.triggered.connect(displayPSD)
            menu.addAction(plotPSD_Button)

            plotButton = QtWidgets.QAction('Plot', self)
            plotButton.setStatusTip('Plot')
            displayEID = partial(self.plotEID, id_event)
            plotButton.triggered.connect(displayEID)
            menu.addAction(plotButton)

            relabelButton = QtWidgets.QAction('Relabel', self)
            relabelButton.setStatusTip('Relabel')
            relabelEID = partial(self.relabelEID, [id_event])
            relabelButton.triggered.connect(relabelEID)
            menu.addAction(relabelButton)

            removeButton = QtWidgets.QAction('Remove', self)
            removeButton.setStatusTip('Remove')
            removePSD = partial(self.removeEID, [id_event])
            removeButton.triggered.connect(removePSD)
            menu.addAction(removeButton)
            
            menu.exec_(event.globalPos())
            return True

        else:
            return False


    def closeEvent(self, event):
        if self.canvas.psd_frame:
            self.canvas.psd_frame.close()
        
        if self.canvas.EID_frame:
            self.canvas.EID_frame.close()


class StationCanvas(FigureCanvas):
    def __init__(self, station, starttime, channel, delta, sde=None, parent=None, **kwargs):

        # load objects
        self.starttime = starttime
        self.station = station
        self.channel = channel
        self.delta = delta
        self.endtime = starttime + dt.timedelta(minutes=delta)
        self.parent = parent
        self.olap = kwargs.get("delta_olap", 0.2)
        self.add_stations = kwargs.get("add_stations", [])
        
        self.plot_kwargs = dict(
            specgram=kwargs.get("specgram", 0),
            interval=kwargs.get("interval", 60),
            fq_band=kwargs.get("fq_band", [0.5, 10]),
            remove_response=kwargs.get("remove_response", True),
            color_name=kwargs.get("color_name", 'zesty')
        )
        
        self.psd_frame = None
        self.EID_frame = None
        self.info_dict = {
            'left':None,
            'right':None,
            'diff':None,
            'freq':None,
            'channel':None
        }
        
        # change starttime
        self.parent.gotoButton.clicked.connect(self.set_starttime)
        self.parent.gotoButton.setShortcut('g')

        # change delta
        self.parent.intervalButton.clicked.connect(self.set_interval)
        self.parent.intervalButton.setShortcut('d')

        # navigate
        self.parent.forwardButton.clicked.connect(self.forward_step)
        self.parent.forwardButton.setShortcut('right')
        self.parent.backwardButton.clicked.connect(self.backward_step)
        self.parent.backwardButton.setShortcut('left')

        # filter, spectrogam, and response
        self.parent.specButton.clicked.connect(self.set_spectrogram)
        self.parent.specButton.setShortcut('backspace')

        self.parent.respButton.clicked.connect(self.remove_response)
        self.parent.respButton.setShortcut('r')
        
        self.parent.filtButton.clicked.connect(self.set_fqband)
        self.parent.filtButton.setShortcut('a')

        self.parent.psdButton.clicked.connect(self.plot_psd)

        if not sde:
            sde = './%s' % station.stats.id
        
        if isinstance(sde, SDE):
            self.sde = sde
        else:
            self.sde = SDE(sde)
        
        # load canvas
        figsize = kwargs.get("figsize", (17,5))
        dpi = kwargs.get("dpi", 100)
        self.fig = Figure(figsize=figsize, dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, 
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self.mpl_connect('button_release_event', self.on_click)
        self.mpl_connect('key_press_event', self.on_key)
        self.__plot()


    def __plot(self):
        self.fig.clf()
        self.parent.plainTextEdit.setText('')

        with pyqtgraph.BusyCursor():
            ans = get_fig(
                self.station, 
                self.starttime,
                self.channel,
                self.endtime,
                return_axes=True,
                fig=self.fig,
                **self.plot_kwargs
                )
        
        self.chan_axes, self.specgram_axes, self.time, self.channel_list = ans

        if self.specgram_axes:
            self.nav = Navigation(self.chan_axes, imshow_axes=self.specgram_axes[0], parent=self)
        else:
            self.nav = Navigation(self.chan_axes, parent=self)

        self.__eventlist = {}
        self.__show_events()

        self.showCatalog = False
        with pyqtgraph.BusyCursor():
            self.draw()


    def show_tickinfo(self):
        text = f" {self.info_dict['channel']}:\n"
        text += f" L: {self.info_dict['left']}\n"
        text += f" R: {self.info_dict['right']}\n"

        if self.info_dict['diff']:
            text += f" R-L: {self.info_dict['diff']:.2f} sec\n"
        
        if self.info_dict['freq']:
            text += f" Freq.: {self.info_dict['freq']:.2f}\n"
        self.parent.plainTextEdit.setText(text)


    def __update_phase(self, phase, item):
        info_dict = {}
        time = item['phases'][phase]['time']
        
        if time:
            if phase == 'p':
                info_dict['time_P'] = (time - item['event'].starttime).total_seconds()

            elif phase == 's':
                info_dict['time_S'] = (time - item['event'].starttime).total_seconds()
            
            else:
                if item['phases']['p']['time']:
                    info_dict['event_duration'] = (time - item['phases']['p']['time']).total_seconds()

        else:
            if phase == 'p':
                info_dict['time_P'] = None 
                info_dict['event_duration'] = None
            
            elif phase == 's':
                info_dict['time_S'] = None
            
            else:
                info_dict['event_duration'] = None
        
        row = item['event'].get_row(self.station.stats.id)
        self.sde.update_row(row.id, info_dict)
        self.__show_events(draw=True)


    def on_key(self, event):
        if event.inaxes in self.chan_axes and event.key in ('p', 'P', 's', 'S', 'f', 'F') and self.__eventlist:
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            list_eid = []
            for id, item in self.__eventlist.items():
                evnt = item['event']
                if evnt.starttime < time < evnt.endtime:
                    list_eid.append(id)
            
            if list_eid:
                if len(list_eid) == 1: # item and id are defined
                    if event.key in ('p', 's', 'f'): # add/change phase
                        if item['phases'][event.key]['time']: # change phase
                            item['phases'][event.key]['time'] = time
                            [tick.remove() for tick in item['phases'][event.key]['tick']]
                            item['phases'][event.key]['txt'].remove()
                        
                        else: # add phase
                            item['phases'][event.key]['time'] = time
                        
                        self.__update_phase(event.key, item)
                        
                    
                    if event.key in ('P', 'S', 'F'):
                        if item['phases'][event.key.lower()]['time']: # remove phase
                            item['phases'][event.key.lower()]['time'] = None
                            [tick.remove() for tick in item['phases'][event.key.lower()]['tick']]
                            item['phases'][event.key.lower()]['txt'].remove()
                        
                        if event.key == 'P':
                            item['phases']['f']['time'] = None
                            [tick.remove() for tick in item['phases']['f']['tick']]
                            item['phases']['f']['txt'].remove()
                        
                        self.__update_phase(event.key.lower(), item)
                
                else:
                    print(' Two IDs found.')
            
            return True
        
        # relabel event
        if event.inaxes in self.chan_axes and event.key == 'L' and self.__eventlist:
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            id_to_relabel = []
            for id, item in self.__eventlist.items():
                evnt = item['event']
                if evnt.starttime < time < evnt.endtime:
                    id_to_relabel.append(id)

            if id_to_relabel:
                self.relabelEID(id_to_relabel)

        # remove event
        if event.inaxes in self.chan_axes and event.key == 'delete' and self.__eventlist:
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            id_to_remove = []
            for id, item in self.__eventlist.items():
                evnt = item['event']
                if evnt.starttime < time < evnt.endtime:
                    id_to_remove += [id]
                    [s.remove() for s in item['span']]
                    item['txt'].remove()
            
            if id_to_remove:
                self.removeEID(id_to_remove)


        # add event
        if self.info_dict['left'] and self.info_dict['right']:
            if event.key in map(str, SHORTCUT_LABELS.keys()) or event.key == 'l':
                event_to_save = {}
                
                if event.key == 'l':
                    lbl, ok = QtWidgets.QInputDialog.getText(self, "Save event", "Label:")
                    if ok:
                        event_to_save['label'] = lbl.upper()
                    else:
                        return
                
                else:
                    event_to_save['label'] = SHORTCUT_LABELS[int(event.key)]
                
                if self.info_dict['channel'][-1] == 'P':
                    event_to_save['event_type'] = 'P'
                else:
                    event_to_save['event_type'] = 'S'
                
                event_to_save['starttime'] = min(self.info_dict['right'], self.info_dict['left'])
                event_to_save['duration'] = self.info_dict['diff']
                event_to_save['network'] = self.station.stats.id.split('.')[0]
                event_to_save['station'] = self.station.stats.id.split('.')[1]
                event_to_save['location'] = self.station.stats.id.split('.')[2]
                event_to_save['event_id'] = self.sde.last_eid() + 1

                self.sde.add_row(event_to_save)
                self.__show_events(draw=True)


    def on_click(self, event):
        if event.inaxes in self.chan_axes:
            for ax in self.chan_axes:
                if ax == event.inaxes:
                    channel = ax.yaxis.get_label().get_text()
                    self.info_dict['channel'] = channel

                    if self.nav.ticks['right'][0]:
                        r_time = mdates.num2date(self.nav.ticks['right'][0]).replace(tzinfo=None)
                    else:
                        r_time = None
                    
                    self.info_dict['right'] = r_time

                    if self.nav.ticks['left'][0]:
                        l_time = mdates.num2date(self.nav.ticks['left'][0]).replace(tzinfo=None)
                    else:
                        l_time = None
                    
                    self.info_dict['left'] = l_time
                    
                    if r_time and l_time:
                        tuple_time = (r_time, l_time)
                        diff = (max(tuple_time) - min(tuple_time)).total_seconds()
                    else:
                        diff = None
                    
                    self.info_dict['diff'] = diff
                    self.show_tickinfo()

        if event.inaxes == self.specgram_axes[0]:
            channel = event.inaxes.yaxis.get_label().get_text()
            freq = event.ydata
            self.info_dict['freq'] = freq
            self.show_tickinfo()


    def __show_events(self, draw=False):
        self.parent.listWidget.clear()
        
        # remove all span
        if self.__eventlist:
            for _, item in self.__eventlist.items():
                [s.remove() for s in item['span']]
                item['txt'].remove()
            self.__eventlist = {}

        # create __eventlist
        eid_list = self.sde.get_eid_list(time_interval=(self.starttime, self.endtime))
        if eid_list:
            for eid in eid_list:
                event = self.sde[eid]
                row = event.get_row(self.station.stats.id)
                
                if row.time_P:
                    Ptime = {
                        'time':event.starttime + dt.timedelta(seconds=row.time_P),
                        'tick':[],
                        'txt':None
                    }
                else:
                    Ptime = {
                        'time':None,
                        'tick':[],
                        'txt':None
                    }
                
                if row.time_S:
                    Stime = {
                        'time':event.starttime + dt.timedelta(seconds=row.time_S),
                        'tick':None,
                        'txt':None
                    }
                else:
                    Stime = {
                        'time':None,
                        'tick':[],
                        'txt':None
                    }
                
                if row.event_duration and row.time_P:
                    Ftime = {
                        'time':event.starttime + dt.timedelta(seconds=row.time_P + row.event_duration),
                        'tick':[],
                        'txt':None
                    }
                else:
                    Ftime = {
                        'time':None,
                        'tick':[],
                        'txt':None
                    }
                
                self.__eventlist[eid] = {
                    'span': [],
                    'txt': None,
                    'event': event,
                    'phases': {
                        'p': Ptime,
                        's': Stime,
                        'f': Ftime
                    }
                }

        # print events
        for eid, item in self.__eventlist.items():
            st = item['event'].starttime
            et = item['event'].endtime

            cond_in = st >= self.starttime and et <= self.endtime # start and end in
            cond_st_in = self.endtime > st > self.starttime and et > self.endtime # start in, end out
            cond_et_in = self.starttime > st and self.starttime < et < self.endtime # start out, end in
            # cond_full_in = st < self.starttime and endtime < et # start and end out

            if cond_in:
                st1 = st
                et1 = et

            elif cond_st_in:
                st1 = st
                et1 = self.endtime

            elif cond_et_in:
                st1 = self.starttime
                et1 = et

            else:# cond_full_in:
                st1 = self.starttime
                et1 = self.endtime

            for a, ax in enumerate(self.chan_axes):
                if item['event'].label in [l for _, l in SHORTCUT_LABELS.items()]:
                    n = [l for _, l in SHORTCUT_LABELS.items()].index(item['event'].label)
                    c = SHORTCUT_COLORS[n]
                else:
                    c = 'teal'
                self.__eventlist[eid]['span'] += [ax.axvspan(st1, et1, color=c, alpha=0.5)]

                if a == 0:
                    text_pos = mdates.date2num(st1 + dt.timedelta(seconds=0.3))
                    tx = ax.annotate(
                        'EID:%s(%s)' % (eid, item['event'].label), 
                        (text_pos, 1),
                        color='k',
                        fontsize=9)
                    self.__eventlist[eid]['txt'] = tx
            
            info_evnt = '  %s [%s]' % (eid, item['event'].label)
            self.parent.listWidget.addItem(info_evnt)
        
        # print event phases
        off = dt.timedelta(seconds=0.25)
        for eid, item in self.__eventlist.items():
            for n, ax in enumerate(self.chan_axes):
                if item['phases']['p']['time']:
                    tick = ax.axvline(item['phases']['p']['time'], color=PHASE_COLORS['p'], lw=1.1)
                    item['phases']['p']['tick'].append(tick)
                    if n == 0:
                        txt = ax.annotate('P', xy=(item['phases']['p']['time'] + off, 0.8), color=PHASE_COLORS['p'], fontsize=9)
                        item['phases']['p']['txt'] = txt
                
                if item['phases']['s']['time']:
                    tick = ax.axvline(item['phases']['s']['time'], color=PHASE_COLORS['s'], lw=1.1)
                    item['phases']['s']['tick'].append(tick)
                    if n == 0:
                        txt = ax.annotate('S', xy=(item['phases']['s']['time'] + off, 0.8), color=PHASE_COLORS['s'], fontsize=9)
                        item['phases']['s']['txt'] = txt
                
                if item['phases']['f']['time']:
                    tick = ax.axvline(item['phases']['f']['time'], color=PHASE_COLORS['f'], lw=1.1)
                    item['phases']['f']['tick'].append(tick)
                    if n == 0:
                        txt = ax.annotate('F', xy=(item['phases']['f']['time'] + off, 0.8), color=PHASE_COLORS['f'], fontsize=9)
                        item['phases']['f']['txt'] = txt

        if draw:
            self.draw()


    def plotPSDeid(self, eid):
        event = self.sde[eid]
        st = event.starttime
        et = event.endtime
        self.plot_psd(times=[st, et])
        

    def plotEID(self, eid):
        self.EID_frame = self.sde.plot_gui(
            eid=eid, 
            remove_response=self.plot_kwargs['remove_response'],
            fq_band=self.plot_kwargs['fq_band'],
            stations=self.add_stations,
            app=True)
        

    def removeEID(self, id_to_remove):
        for id in id_to_remove:
                self.sde.remove_event(id)
                del self.__eventlist[id]
            
        self.parent.listWidget.clear()
        for id, item in self.__eventlist.items():
            info_evnt = '  %s [%s]' % (id, item['event'].label)
            self.parent.listWidget.addItem(info_evnt)
        self.draw()
    

    def relabelEID(self, id_to_relabel):
        for id in id_to_relabel:
            evnt = self.sde[id]
            new_label, ok = QtWidgets.QInputDialog.getText(self, "Relabel Event","Event Label: ", text=evnt.label)
            if ok and new_label != evnt.label:
                self.sde.relabel_event(id, new_label.upper())
        
        self.__show_events(draw=True)


    def export_to_mlf(self):
        path = './dataset'
        mseed_path = '%s/data' % path
        zmlffile = '%s/BTs_Z.mlf' % path

        if not os.path.isdir(path):
            os.mkdir(path)

        if not os.path.isdir(mseed_path):
            os.mkdir(mseed_path)

        if not os.path.isfile(zmlffile):
            with open(zmlffile, 'w') as f:
                f.write('#!MLF!#\n')
                f.write('# Copahue Volcano\n')
                f.write('# Ivan Melchor - ifmelchor@unrn.edu.ar\n')
                f.write('.\n')

        lbl, ok = QtWidgets.QInputDialog.getText(self, "Event type: S(eismic)", "Label:")
        if ok:
            times = [self.nav[0].ticks['left'][0], self.nav[0].ticks['right'][0]]
            start_tick = mdates.num2date(min(times))
            start_tick = start_tick.replace(tzinfo=None)
            start_tick = self.time[np.argmin(np.abs(np.array(self.time)-start_tick))]
            end_tick = mdates.num2date(max(times))
            end_tick = end_tick.replace(tzinfo=None)
            end_tick = self.time[np.argmin(np.abs(np.array(self.time)-end_tick))]
            stream = self.station.get_stream(start_tick, end_tick, remove_response=True)
            vtrace = stream.get_component('Z')
            etrace = stream.get_component('E')
            mseed_efile = './data/%s%sE.mseed' % (self.code_base, start_tick.strftime('%y%m%d_%H%M%S'))
            factor = 10000


            with pyqtgraph.BusyCursor():

                with open(zmlffile, 'a') as f:
                    mseed_vfile_to_write = './data/%s%sZ.mseed' % (self.code_base, start_tick.strftime('%y%m%d_%H%M%S'))
                    mseed_vfile = '%s/%s%sZ.mseed' % (mseed_path, self.code_base, start_tick.strftime('%y%m%d_%H%M%S'))
                    vtrace.write(mseed_vfile, format='MSEED')
                    f.write('.\n"'+mseed_vfile_to_write+'"\n')

                    sub_evnts = self.sde.get_events(time_interval=(start_tick, end_tick))
                    for n, (ei, ej) in enumerate(zip(sub_evnts[:-1], sub_evnts[1:])):
                        start_i = int(abs((start_tick - ei.starttime).total_seconds()))*factor
                        if n == 0 and start_i > 0:
                            f.write("{0:7s} {1:7s} NS\n".format('0', str(start_i)))
                        
                        end_event = ei.starttime + dt.timedelta(minutes=ei.duration)
                        end_i = int(abs((start_tick - end_event).total_seconds()))*factor
                        f.write("{0:7s} {1:7s} {2:s}\n".format(str(start_i), str(end_i), lbl.upper()))

                        start_j = int(abs((start_tick - ej.starttime).total_seconds()))*factor
                        if start_j > end_i:
                            f.write("{0:7s} {1:7s} NS\n".format(str(end_i), str(start_j)))

                        if n == len(sub_evnts)-2:
                            end_event = ej.starttime + dt.timedelta(minutes=ej.duration)
                            end_j = int(abs((start_tick - end_event).total_seconds()))*factor
                            f.write("{0:7s} {1:7s} {2:s}\n".format(str(start_j), str(end_j), lbl.upper()))

                            endtick = int(abs((start_tick - end_tick).total_seconds()))*factor
                            if endtick > end_j:
                                f.write("{0:7s} {1:7s} NS\n".format(str(end_j), str(endtick)))


    def forward_step(self):
        step = dt.timedelta(minutes=self.delta)
        olap = dt.timedelta(minutes=self.delta*self.olap)
        self.starttime = self.starttime - olap + step
        self.endtime = self.starttime + step
        self.__plot()


    def backward_step(self):
        step = dt.timedelta(minutes=self.delta)
        olap = dt.timedelta(minutes=self.delta*self.olap)
        self.starttime = self.starttime + olap - step
        self.endtime = self.starttime + step
        self.__plot()


    def set_interval(self):
        i, ok = QtWidgets.QInputDialog.getInt(self, "Set interval","Interval:", self.delta, 5, 150)
        if ok and int(i) != self.delta:
            self.delta = i
            step = dt.timedelta(minutes=self.delta)
            self.endtime = self.starttime + step
            self.__plot()


    def set_starttime(self):
        current_day = self.starttime.strftime("%Y%m%d%H%M")
        time, ok = QtWidgets.QInputDialog.getText(self, "Set Starttime", "Time (YYYYMMDDHHMM):", text=current_day)
        if ok and str(time) != current_day:
            try:
                time = dt.datetime.strptime(str(time), "%Y%m%d%H%M")
                self.starttime = time
                self.endtime = self.starttime + dt.timedelta(minutes=self.delta)
                self.__plot()

            except:
                return False


    def plot_psd(self, times=[]):
        if self.psd_frame:
            self.psd_frame.close()
            self.psd_frame = None

        tick_times = (self.nav.ticks['right'][0], self.nav.ticks['left'][0])
        if not times and not all(tick_times):
            return
        
        if times:
            start = min(times)
            end = max(times)
        else:
            r_time = mdates.num2date(self.nav.ticks['right'][0]).replace(tzinfo=None)
            l_time = mdates.num2date(self.nav.ticks['left'][0]).replace(tzinfo=None)
            start = min([r_time, l_time])
            end = max([r_time, l_time])
        
        if (end - start).total_seconds()/60 >= 5:
                avg_step = 1
        else:
            avg_step = None

        with pyqtgraph.BusyCursor():
            stream = self.station.get_stream(start, end, remove_response=self.plot_kwargs['remove_response'])

            l_trace = []
            l_color = []
            for c, chan in enumerate(self.channel_list):
                trace = stream.get(channel=chan)
                
                psd, freq = trace.psd(
                    fq_band=self.plot_kwargs['fq_band'], 
                    avg_step=avg_step, 
                    olap=0.25, 
                    return_fig=False, 
                    plot=False
                    )
                l_trace.append(psd)
                l_color.append(get_colors(self.plot_kwargs['color_name'])[c])
                # psd = np.array([psd])
            
            if len(l_trace) == 3:
                pd, _ = stream.polardegree(fq_band=self.plot_kwargs['fq_band'], avg_step=avg_step, 
                    olap=0.25, plot=False)
            else:
                pd = None

            # add matrix return and rect/dip/azm plots
            self.psd_frame = PSD_GUI(freq, np.vstack(l_trace), pd, plot_prob=False, title='PSD/PD', colors=l_color)
            self.psd_frame.show()


    def remove_response(self):
        if self.plot_kwargs['remove_response']:
            self.plot_kwargs['remove_response'] = False
        else:
            self.plot_kwargs['remove_response'] = True
        self.__plot()


    def set_spectrogram(self):
        if self.specgram_axes:
            chn = self.specgram_axes[0].yaxis.get_label().get_text()
            items = self.channel_list
            index = items.index(chn)
            item, ok = QtWidgets.QInputDialog.getItem(self, 
                "select channel", 
                "list of channels", 
                items, 
                index, 
                False
            )
            if ok:
                self.plot_kwargs['specgram'] = item
                self.__plot()
    

    def set_fqband(self):
        if self.plot_kwargs['fq_band']:
            current_filter = ';'.join(map(str, self.plot_kwargs['fq_band']))
        else:
            current_filter = ''
        
        filt_text, ok = QtWidgets.QInputDialog.getText(self, "Set filter", "fq_min;fq_max:", text=current_filter)
        if ok and str(filt_text) != current_filter:
            fq_min = float(filt_text.split(';')[0])
            fq_max = float(filt_text.split(';')[1])
            self.plot_kwargs['fq_band'] = (fq_min, fq_max)
            self.__plot()


    def show_catalog(self):
        print('\n ============== CATALOG ==============')
        print(f" n |  {'Time':^8}  {'M':^6}  {'Distance'}  {'Depth'}")
        with pyqtgraph.BusyCursor():
            if self.showCatalog:
                [arrival.remove() for arrival in self.USGSarrivals['ticks']]
                [arrival.remove() for arrival in self.USGSarrivals['text']]
                self.showCatalog = False
                self.draw()
                
            else:
                self.USGSarrivals = {'ticks':[], 'text':[]}
                off = dt.timedelta(minutes=10)
                cat = get_catalog(self.starttime-off, self.endtime, receiver=(self.station.stats.lat, self.station.stats.lon))
                if cat:
                    for n, c in enumerate(cat):
                        if c[-1]:
                            if self.starttime <= c[-1] <= self.endtime:
                                print(f" {n:>2} | {c[0].strftime('%H:%M:%S')}  {c[4]}  {c[2]:.1f} km {c[3]:.0f} km")

                                if c[2] > 500:
                                    col = 'k'
                                else:
                                    col = 'r'

                                for ax in self.chan_axes:
                                    arrival = ax.axvline(c[-1], color=col)
                                    self.USGSarrivals['ticks'].append(arrival)
                                text_pos = mdates.date2num(c[-1] + dt.timedelta(seconds=0.3))
                                self.USGSarrivals['text'].append(self.chan_axes[0].annotate(str(n), (text_pos, 1), color='k', fontsize=9))
                    
                    if self.USGSarrivals['ticks']:
                        self.showCatalog = True
                        self.draw()


def get_fig(station, starttime, channel, endtime, return_axes=False, **kwargs):
    """[summary]

    Args:
        station: Station object
        starttime: datetime
        channel (string): the channel to plot
        delta (int): step time
        return_fig (bool, optional): Return fig and axes. Defaults to False.
    """

    # define some variables
    fq_band = kwargs.get('fq_band', (0.5, 10))
    remove_response = kwargs.get('remove_response', True)
    sample_rate = kwargs.get('sample_rate', None)
    color_name = kwargs.get('color_name', 'zesty')
    specgram = kwargs.get('specgram', 0)

    colors = get_colors(color_name)
    
    if isinstance(specgram, str):
        if specgram not in station.stats.chan:
            specgram = None
    
    elif isinstance(specgram, int):
        if specgram < len(station.stats.chan)-1:
            specgram = station.stats.chan[specgram]
        else:
            specgram = None
    
    else:
        specgram = None

    # define the stream
    if channel != 'all':
        if not channel in station.stats.chan:
            raise ValueError('channel %s not in %s' % (channel, station.stats.chan))
        else:
            stream = station.get_stream(starttime, endtime, sample_rate=sample_rate, chan=channel, remove_response=remove_response)
            channel_list = [channel]
            if not stream:
                raise ValueError(' no data for channel %s' %channel)
    else:
        stream = station.get_stream(starttime, endtime, sample_rate=sample_rate, remove_response=remove_response)
        channel_list = [tr.stats.channel for tr in stream]

    fig = kwargs.get('fig', None)
    grid = {'hspace':0.1, 'left':0.08, 'right':0.92, 'wspace':0.05, 'top':0.95, 'bottom':0.05}
    
    if specgram:
        nrows = 1 + len(channel_list)
        ncols = 2
        grid['width_ratios'] = [1, 0.01]
    else:
        nrows = len(channel_list)
        ncols = 1
    
    if not fig:
        dpi = kwargs.get('dpi', 100)
        figsize = kwargs.get('figsize', (20,9))
        fig, axes = plt.subplots(nrows, ncols, gridspec_kw=grid, figsize=figsize, dpi=dpi)
    else:
        axes = fig.subplots(nrows, ncols, gridspec_kw=grid)
    
    if isinstance(axes, np.ndarray):
        axes = axes.reshape(nrows, ncols)
    else:
        axes = np.array([axes]).reshape(nrows, ncols)

    # if specgram, plot first
    if specgram:
        spec_ax = axes[0,0]
        spec_bar_ax = axes[0,1]
        spec_stream = station.get_stream(starttime, endtime, chan=specgram, remove_response=False)
        spec_stream[0].specgram(
            axes=spec_ax, 
            fq_band=fq_band, 
            per_lap=0.75,
            xlabel=' ',
            axis_bar=spec_bar_ax, 
            axis_bar_label='PSD\n'+r'[dB cnts$^2$/Hz]')
        chn = spec_stream[0].stats.channel
        spec_ax.set_ylabel(chn, fontsize=12, color=colors[channel_list.index(chn)])
        spec_ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        spec_ax.xaxis.set_major_formatter(mtick.NullFormatter())
        spec_return = (spec_ax, spec_bar_ax)
    else:
        spec_return = ()

    # plot channels
    chan_axes = []
    time = None
    for i, channel in enumerate(channel_list):
        if specgram:
            i += 1
            axes[i, 1].axes.get_xaxis().set_visible(False)
            axes[i, 1].axes.get_yaxis().set_visible(False)
            axes[i, 1].set_frame_on(False)
        
        ax = axes[i,0]
        
        if channel == channel_list[0]:
            if remove_response:
                max_amplitude_text = r'Max. amplitude in $\mu m$'
            else:
                max_amplitude_text = r'Max. amplitude in cnts'
            ax.annotate(
                max_amplitude_text, 
                xy=(0,0.9), 
                xycoords='axes fraction', 
                color='k'
                )
        
        trace = stream.get(channel=channel)
        if time == None:
            time = trace.get_time()

        data = trace.get_data(detrend=True, fq_band=fq_band)
        max_ampl = np.nanmax(np.abs(data))
        norm_data = data/max_ampl
        ax.plot(time, norm_data, color=colors[channel_list.index(channel)], lw=1.0)
        ax.set_ylabel(channel, fontsize=12, color=colors[channel_list.index(channel)])
        ax.set_ylim(-1, 1)

        if remove_response:
            max_ampl = max_ampl*10**6 # micrometers/sec

        ax.annotate(
            f'{max_ampl:.2f}',
            xy=(-0.01,0.75),
            xycoords='axes fraction',
            color=colors[channel_list.index(channel)],
            bbox=dict(boxstyle="round", fc="w", alpha=1)
        )

        if channel == channel_list[-1]:
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        else:
            ax.xaxis.set_major_formatter(mtick.NullFormatter())

        ax.yaxis.set_major_locator(mtick.NullLocator())
        ax.grid(axis='x', which='major', color='k', ls='--', alpha=0.4)
        ax.grid(axis='x', which='minor', color='k', ls='--', alpha=0.2)
        chan_axes.append(ax)

    if specgram:
        xaxis = [spec_ax] + chan_axes
    else:
        xaxis = chan_axes
    
    [ax.set_xlim(time[0], time[-1]) for ax in xaxis]
    [ax.xaxis.set_major_locator(mtick.MaxNLocator(nbins=6)) for ax in xaxis]
    [ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(4)) for ax in xaxis]

    if (time[-1] - time[0]).total_seconds() > 86400:
        xaxis[0].set_title('%s--%s', time[0].strftime('%d %b %Y'), time[-1].strftime('%d %b %Y'))
    else:
        xaxis[0].set_title(time[0].strftime('%d %b %Y'))
    
    fig.align_ylabels()

    if return_axes:
        return chan_axes, spec_return, time, channel_list

    else:
        return fig, axes


def plot_station_gui(station, starttime, sde, channel='all', delta=30, app=False, **kwargs):
    if not app:
        app = QtWidgets.QApplication([])
        _exec = True
    else:
        _exec = False
    
    sta_frame = StationWindow(station, starttime, channel, delta, sde, **kwargs)
    sta_frame.show()
    
    if _exec:
        app.exec_()
    
    else:
        return sta_frame



