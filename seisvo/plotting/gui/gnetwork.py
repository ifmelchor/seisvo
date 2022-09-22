

import numpy as np
import datetime as dt
import pyqtgraph
from functools import partial

from seisvo import SDE
from seisvo.utils import get_times_bounds
from seisvo.core.obspyext import Stream2
from seisvo.plotting import get_colors
from seisvo.plotting.gui import Navigation, PSD_GUI
from seisvo.plotting.gui.frames import MainGUI


import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
from matplotlib.gridspec import GridSpec
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas

SHORTCUT_COLORS = get_colors('tolb')

DEFAULT_LABEL_COLOR = "teal"
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

def get_label_color(label):
    if label in [l for _, l in SHORTCUT_LABELS.items()]:
        n = [l for _, l in SHORTCUT_LABELS.items()].index(label)
        return SHORTCUT_COLORS[n]
    else:
        return DEFAULT_LABEL_COLOR


class NetworkCanvas(FigureCanvas):
    def __init__(self, network, sta_list, starttime, delta, component, sde_file, parent=None, **kwargs):

        # load objects
        self.network = network
        self.sta_list = sta_list
        self.starttime = starttime
        self.component = component
        self.delta = delta
        self.parent = parent # mainwindow of the GUI

        self.olap = kwargs.get("delta_olap", 0.2)
        self.plot_kwargs = dict(
            fq_band = kwargs.get("fq_band", ()),
            color_name = kwargs.get("color_name", "tolm"),
            remove_response = kwargs.get("remove_response", False),
            sample_rate = kwargs.get("sample_rate", None),
            specgram = kwargs.get("specgram")
        )

        self.sde = SDE(sde_file)
        
        self.loadGui()
        self.setPlot()


    def loadGui(self):

        self.PSDframe = None
        self.EIDframe = None
        
        self.info_dict = {
            'left':None,
            'right':None,
            'diff':None,
            'freq':None,
            'trace':None
        }
        
        # navigate
        self.parent.gotoButton.clicked.connect(self.setStarttime)
        self.parent.gotoButton.setShortcut('g')

        self.parent.intervalButton.clicked.connect(self.setInterval)
        self.parent.intervalButton.setShortcut('d')

        self.parent.forwardButton.clicked.connect(self.forwardStep)
        self.parent.forwardButton.setShortcut('right')

        self.parent.backwardButton.clicked.connect(self.backwardStep)
        self.parent.backwardButton.setShortcut('left')

        # filter, spectrogam, and response
        self.parent.specButton.clicked.connect(self.setResponse)
        self.parent.specButton.setShortcut('backspace')

        self.parent.respButton.clicked.connect(self.setResponse)
        self.parent.respButton.setShortcut('r')
        
        self.parent.filtButton.clicked.connect(self.setFQband)
        self.parent.filtButton.setShortcut('a')

        # filter, spectrogam, and response

        self.parent.psdButton.clicked.connect(self.plotPSD)

        self.fig = Figure(figsize=(20,9))
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, 
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self.mpl_connect('button_release_event', self.onClick)
        self.mpl_connect('key_press_event', self.onKey)


    def setPlot(self):
        self.fig.clf()
        self.parent.plainTextEdit.setText('')

        self.endtime = self.starttime + dt.timedelta(minutes=self.delta)
        with pyqtgraph.BusyCursor():
            (self.trace_axes, self.specgram_axes, self.time, self.trace_list) = network_plot(self.network, self.sta_list, self.starttime, self.component, self.endtime, full_return=True, fig=self.fig, **self.plot_kwargs)
        
        if self.specgram_axes:
            self.nav = Navigation(self.trace_axes, imshow_axes=self.specgram_axes[0], parent=self)
        else:
            self.nav = Navigation(self.trace_axes, parent=self)

        self._sta_ids = '.'.join([self.network.stats.code, tr_id] for tr_id in self.trace_list)
        self.eventlist = {}
        
        self.drawEvents()
        self.showCatalog = False
        
        with pyqtgraph.BusyCursor():
            self.draw()
    

    def showTickInfo(self):
        text = f" {self.info_dict['trace']}:\n"
        text += f" L: {self.info_dict['left']}\n"
        text += f" R: {self.info_dict['right']}\n"

        if self.info_dict['diff']:
            text += f" R-L: {self.info_dict['diff']:.2f} sec\n"
        
        if self.info_dict['freq']:
            text += f" Freq.: {self.info_dict['freq']:.2f}\n"
        
        self.parent.plainTextEdit.setText(text)

    # --------------
    #   DATABASE
    # --------------

    def getEIDforTime(self, time, sta_id):
        eid_list = []
        for eid, item in self.eventlist.items():
            event = item['event']
            if sta_id not in event.stations:
                if event.starttime <= time <= event.endtime:
                    eid_list += [eid]
        
        if len(eid_list) == 1:
            return eid_list[0]
        
        elif len(eid_list) > 1:
            return eid_list
        
        else:
            return None
        

    def drawEvents(self, draw=False):
        self.parent.listWidget.clear()
        
        # clean previous eventlist
        if self.eventlist:
            for _, item in self.eventlist.items():
                [s.remove() for s in item['span']]
                item['txt'].remove()
            self.eventlist = {}
            
        # get the EID list
        eid_list = self.sde.get_eid_list(time_interval=(self.starttime, self.endtime))
        if eid_list:
            for eid in eid_list:
                event = self.sde[eid]
                # create empty dict for eid
                self.eventlist[eid] = {'span': [], 'txt': None, 'event': event} 
                # define the times
                st1, et1 = get_times_bounds(self.starttime, self.endtime, event.starttime, event.endtime)
                # define the color
                c = get_label_color(event.label)

                n = 0
                for sta_id in event.stations:
                    if sta_id in self._sta_ids:
                        ax = self.trace_axes[self._sta_ids.index(sta_id)]
                        self.eventlist[eid]['span'] += [ax.axvspan(st1, et1, color=c, alpha=0.5)]

                        if n == 0:
                            text_pos = mdates.date2num(st1 + dt.timedelta(seconds=0.3))
                            tx = ax.annotate(
                                'EID:%s(%s)' % (eid, event.label), 
                                (text_pos, 1),
                                color='k',
                                fontsize=9)
                            self.eventlist[eid]['txt'] = tx
                        
                        n += 1
                
                info_evnt = '  %s [%s]' % (eid, event.label)
                self.parent.listWidget.addItem(info_evnt)
        
        if draw:
            self.draw()


    def onKey(self, event):
        if event.inaxes in self.trace_axes and event.key in ('p', 'P', 's', 'S', 'f', 'F') and self.eventlist:
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            list_eid = []
            for id, item in self.eventlist.items():
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
                        
                        self.updatePhase(event.key, item)
                        
                    
                    if event.key in ('P', 'S', 'F'):
                        if item['phases'][event.key.lower()]['time']: # remove phase
                            item['phases'][event.key.lower()]['time'] = None
                            [tick.remove() for tick in item['phases'][event.key.lower()]['tick']]
                            item['phases'][event.key.lower()]['txt'].remove()
                        
                        if event.key == 'P':
                            item['phases']['f']['time'] = None
                            [tick.remove() for tick in item['phases']['f']['tick']]
                            item['phases']['f']['txt'].remove()
                        
                        self.updatePhase(event.key.lower(), item)
                
                else:
                    print(' Two IDs found.')
            
            return True
        
        # relabel event
        if event.inaxes in self.trace_axes and event.key == 'L' and self.eventlist:
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            id_to_relabel = []
            for id, item in self.eventlist.items():
                evnt = item['event']
                if evnt.starttime < time < evnt.endtime:
                    id_to_relabel.append(id)

            if id_to_relabel:
                self.relabelEID(id_to_relabel)

        # remove event
        if event.inaxes in self.trace_axes and event.key == 'delete' and self.eventlist:
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            id_to_remove = []
            for id, item in self.eventlist.items():
                evnt = item['event']
                if evnt.starttime < time < evnt.endtime:
                    id_to_remove += [id]
                    [s.remove() for s in item['span']]
                    item['txt'].remove()
            
            if id_to_remove:
                self.removeEID(id_to_remove)
        
        # add station to eid
        if event.inaxes in self.trace_axes and event.key == '+' and self.eventlist:
            trace = ax.yaxis.get_label().get_text()
            sta_id = '.'.join([self.network.stats.code, trace])
            time = mdates.num2date(event.xdata).replace(tzinfo=None)
            eid = self.getEIDforTime(time, sta_id)
            # if isinstance(eid, list): choose eid
            if eid:
                self.addStaEID(eid, sta_id)

        # remove station to eid
        if event.inaxes in self.trace_axes and event.key == '-' and self.eventlist:
            trace = ax.yaxis.get_label().get_text()
            sta_id = '.'.join([self.network.stats.code, trace])
            time = mdates.num2date(event.xdata).replace(tzinfo=None)


        # add event
        if self.info_dict['left'] and self.info_dict['right']:
            if event.key in map(str, SHORTCUT_LABELS.keys()) or event.key == 'l':
                if event.key == 'l':
                    lbl, ok = QtWidgets.QInputDialog.getText(self, "Save event", "Label:")
                    if ok:
                        label = lbl.upper()
                    else:
                        return
                else:
                    label = SHORTCUT_LABELS[int(event.key)]
                self.addEID(label)


    def onClick(self, event):
        if event.inaxes in self.trace_axes:
            for ax in self.trace_axes:
                if ax == event.inaxes:
                    trace = ax.yaxis.get_label().get_text()
                    self.info_dict['trace'] = trace

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
            freq = event.ydata
            self.info_dict['freq'] = freq
            self.show_tickinfo()


    def addEID(self, label):
        event_to_save = {}
        event_to_save['event_type'] = 'S' # only for seismic
        event_to_save['label'] = label
        event_to_save['starttime'] = min(self.info_dict['right'], self.info_dict['left'])
        event_to_save['duration'] = self.info_dict['diff']
        event_to_save['network'] = self.network.stats.code
        event_to_save['station'] = self.info_dict['trace'].split('.')[0]
        event_to_save['location'] = self.info_dict['trace'].split('.')[1]
        event_to_save['event_id'] = self.sde.last_eid() + 1
        self.sde.add_row(event_to_save)
        self.drawEvents(draw=True)


    def addStaEID(self, eid, sta_id):
        network  = sta_id.split('.')[0]
        station  = sta_id.split('.')[1]
        location = sta_id.split('.')[2]
        self.sde.append_station(eid, network, station, location, "S")
        self.drawEvents(draw=True)

    
    def removeStaEID(self, eid, sta_id):


    def removeEID(self, id_to_remove):
        for id in id_to_remove:
            self.sde.remove_event(id)
            del self.eventlist[id]
            
        self.parent.listWidget.clear()
        for id, item in self.eventlist.items():
            info_evnt = '  %s [%s]' % (id, item['event'].label)
            self.parent.listWidget.addItem(info_evnt)

        self.draw()
    

    def relabelEID(self, sta_id, id_to_relabel):
        for id in id_to_relabel:
            evnt = self.sde[id]
            new_label, ok = QtWidgets.QInputDialog.getText(self, "Relabel Event", "Event Label: ", text=evnt.label)
            
            if ok and new_label != evnt.label:
                self.sde.relabel_event(id, new_label.upper())
        
        self.drawEvents(draw=True)
    
    # -------------
    #   NAVIGATE
    # -------------

    def forwardStep(self):
        step = dt.timedelta(minutes=self.delta)
        olap = dt.timedelta(minutes=self.delta*self.olap)
        self.starttime = self.starttime - olap + step
        self.endtime = self.starttime + step
        self.setPlot()


    def backwardStep(self):
        step = dt.timedelta(minutes=self.delta)
        olap = dt.timedelta(minutes=self.delta*self.olap)
        self.starttime = self.starttime + olap - step
        self.endtime = self.starttime + step
        self.setPlot()


    def setInterval(self):
        i, ok = QtWidgets.QInputDialog.getInt(self, "Set interval","Interval:", self.delta, 5, 150)
        if ok and int(i) != self.delta:
            self.delta = i
            step = dt.timedelta(minutes=self.delta)
            self.endtime = self.starttime + step
            self.setPlot()


    def setStarttime(self):
        current_day = self.starttime.strftime("%Y%m%d%H%M")
        time, ok = QtWidgets.QInputDialog.getText(self, "Set Starttime", "Time (YYYYMMDDHHMM):", text=current_day)
        if ok and str(time) != current_day:
            try:
                time = dt.datetime.strptime(str(time), "%Y%m%d%H%M")
                self.starttime = time
                self.endtime = self.starttime + dt.timedelta(minutes=self.delta)
                self.setPlot()

            except:
                return False

    # --------------
    # CHANGE INPUTS
    # --------------

    def setResponse(self):
        if self.plot_kwargs['remove_response']:
            self.plot_kwargs['remove_response'] = False
        else:
            self.plot_kwargs['remove_response'] = True
        self.setPlot()


    def setSpecgram(self):
        if self.specgram_axes:
            chn = self.specgram_axes[0].yaxis.get_label().get_text()
            index = self.trace_list.index(chn)
            item, ok = QtWidgets.QInputDialog.getItem(self, 
                "select channel", 
                "list of channels", 
                self.trace_list,
                index, 
                False
            )
            if ok:
                self.plot_kwargs['specgram'] = item
                self.setPlot()
    

    def setFQband(self):
        if self.plot_kwargs['fq_band']:
            current_filter = ';'.join(map(str, self.plot_kwargs['fq_band']))
        else:
            current_filter = ''
        
        filt_text, ok = QtWidgets.QInputDialog.getText(self, "Set filter", "fq_min;fq_max:", text=current_filter)
        if ok and str(filt_text) != current_filter:
            fq_min = float(filt_text.split(';')[0])
            fq_max = float(filt_text.split(';')[1])
            self.plot_kwargs['fq_band'] = (fq_min, fq_max)
            self.setPlot()

    # --------------
    #     UTILS
    # --------------

    def plotPSD(self, times=[]):
        if self.PSDframe:
            self.PSDframe.close()
            self.PSDframe = None

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
        
        cmap = get_colors(self.plot_kwargs["color_name"])

        with pyqtgraph.BusyCursor():
            
            dtrace = {"PSDs":[], "colors":[], "labels":[]}

            for n, tr_id in enumerate(self.trace_list):
                sta_code = tr_id.split('.')[0]
                sta_loc = tr_id.split('.')[1]
                sta = self.network.get_sta(sta_code, loc=sta_loc)
                channel = sta.get_chan(self.component)
                trace = sta.get_stream(start, end, channel=channel, remove_response=self.plot_kwargs['remove_response'])[0]
                psd, freq = trace.psd(fq_band=self.plot_kwargs['fq_band'], avg_step=avg_step, olap=0.25, return_fig=False, plot=False)

                if n == 0:
                    dtrace["freq"] = freq

                dtrace["PSDs"] += [psd]
                dtrace["colors"] += [cmap[n]]
                dtrace["labels"] += [tr_id]

            # add matrix return and rect/dip/azm plots
            self.PSDframe = PSD_GUI(dtrace["freq"], psds=np.vstack(dtrace["PSDs"]), title='PSD', colors=dtrace["colors"], labels=dtrace["labels"])
            self.PSDframe.show()

    
    def plotEID(self, eid):

        station_list = [".".join([self.network.code, i]) for i in self.trace_list]

        self.EIDframe = self.sde.plot_gui(
            eid=eid, 
            remove_response=self.plot_kwargs['remove_response'],
            fq_band=self.plot_kwargs['fq_band'],
            stations=station_list,
            app=True)


class NetworkWindow(QtWidgets.QMainWindow, MainGUI):
    def __init__(self, network, sta_list, starttime, delta, component, sde_file, **kwargs):
        super(NetworkWindow, self).__init__()
        self.setupUi(self)
        self.canvas = NetworkCanvas(network, sta_list, starttime, delta, component, sde_file, parent=self, **kwargs)
        
        # self.plotPSD = self.canvas.plotPSDeid
        self.plotEID = self.canvas.plotEID
        self.removeEID = self.canvas.removeEID
        self.relabelEID = self.canvas.relabelEID

        self.gridLayout.addWidget(self.canvas, 1, 1, 1, 1)
        self.setWindowTitle(network.stats.code + ' [{}]'.format(self.canvas.sde.sql_path))
        self.actionOpen.triggered.connect(self.canvas.sde.open)

        # add actions of the windows menu
        # self.saveButton.clicked.connect(self.canvas.show_catalog)
        # self.listWidget.itemClicked.connect(self.canvas.event_selected)
        self.listWidget.installEventFilter(self)


    def eventFilter(self, source, event):
        if (event.type() == QtCore.QEvent.ContextMenu and
            source is self.listWidget):
            item = source.itemAt(event.pos())
            id_event = int(item.text().split('[')[0].strip())
            
            menu = QtWidgets.QMenu()

            # plotPSD_Button = QtWidgets.QAction('PSD', self)
            # displayPSD = partial(self.plotPSD, id_event)
            # plotPSD_Button.triggered.connect(displayPSD)
            # menu.addAction(plotPSD_Button)

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
        if self.canvas.PSDframe:
            self.canvas.PSDframe.close()
        
        if self.canvas.EIDframe:
            self.canvas.EIDframe.close()


def get_text(remove_response, x):
    if remove_response:
        txt = f'Max.: {x*10**6:.2f}'
        return txt + r' $\mu m$'
    else:
        txt = f'Max.: {x:.2f}'
        return txt + ' cnts'


def network_plot(network, station_list, starttime, component, endtime, full_return=False, **kwargs):
    """[summary]

    Args:
        station: List of station objects
        starttime: datetime
        channel: the channel/s to plot
        endtime: datetime
    """

    # define some variables
    fq_band = kwargs.get('fq_band', ())
    remove_response = kwargs.get('remove_response', False)
    sample_rate = kwargs.get('sample_rate', None)
    color_name = kwargs.get('color_name', "tolm")
    specgram = kwargs.get('specgram', station_list[0])
    colors = get_colors(color_name)

    # define the stream
    stream = Stream2()
    trace_ids = [] # for color index and plot
    for sta_str in station_list:
        sta_code = sta_str.split('.')[0]
        staloc_code = sta_str.split('.')[1]
        sta = network.get_sta(sta_code, loc=staloc_code)
        channel = sta.get_chan(component)
        try:
            stream += sta.get_stream(starttime, endtime, sample_rate=sample_rate, chan=channel, remove_response=remove_response)
            trace_ids += [sta_str]
        except:
            print(f" warn: error loading station {sta.stats.id}")

    # define the figure frame
    fig = kwargs.get('fig', None)
    if not fig:
        fig = plt.figure(constrained_layout=True)
    
    true_trace_list = list(set(['.'.join(s.id.split('.')[1:3]) for s in stream]))
    nrows = len(true_trace_list)
    ncols = 1
    
    if specgram in true_trace_list:
        nrows += 1
        ncols += 1
        gs = GridSpec(nrows, ncols, figure=fig, hspace=0.1, left=0.08, right=0.92, wspace=0.05, top=0.95, bottom=0.05, width_ratios=[1, 0.01])
    else:
        gs = GridSpec(nrows, ncols, figure=fig, hspace=0.1, left=0.08, right=0.92, wspace=0.05, top=0.95, bottom=0.05)
        
    # if specgram, plot first
    i = 0
    spec_return = ()
    specgram_trace = stream.get_component(component, station=specgram.split('.')[0], loc=specgram.split('.')[1])
    if specgram_trace:
        
        if remove_response:
            label = r"[dB cnts$^2$/Hz]"
        else:
            label = r"[dB m$^2$s$^{-2}$/Hz]"
        
        spec_ax = fig.add_subplot(gs[i, 0])
        spec_bar_ax = fig.add_subplot(gs[i, 1])
        specgram_trace.plot_specgram(axes=spec_ax, fq_band=fq_band, per_lap=0.75, xlabel="", axis_bar=spec_bar_ax, axis_bar_label=label)
        spec_ax.set_ylabel(specgram, fontsize=12, color=colors[trace_ids.index(specgram)])
        spec_ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        spec_ax.xaxis.set_major_formatter(mtick.NullFormatter())
        spec_return = (spec_ax, spec_bar_ax)
        
        i += 1

        # move the specgram trace to the first position of trace_ids
        spec_pos = trace_ids.index(specgram)
        trace_ids.insert(0, trace_ids.pop(spec_pos))

    # plot traces
    axlist = []
    for n, tr_id in enumerate(trace_ids):
        color = colors[n]
        ax =  fig.add_subplot(gs[i, 0])
        sta_code = tr_id.split('.')[0]
        sta_loc = tr_id.split('.')[1]
        trace = stream.get_component(component, station=sta_code, loc=sta_loc)
        data = trace.get_data(detrend=True, fq_band=fq_band)
        max_ampl = np.nanmax(np.abs(data))
        norm_data = data/max_ampl

        if n == 0:
            time = trace.get_time()
        
        ax.plot(time, norm_data, color=color, lw=1.0)
        ax.set_ylabel(tr_id, fontsize=12, color=color)
        ax.set_ylim(-1.05, 1.05)
        
        ax.annotate(
                get_text(remove_response, max_ampl),
                xy=(-0.01,0.75),
                xycoords='axes fraction',
                color=color,
                bbox=dict(boxstyle="round", fc="w", alpha=1)
            )

        if n == len(tr_id)-1:
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        else:
            ax.xaxis.set_major_formatter(mtick.NullFormatter())
        
        ax.yaxis.set_major_locator(mtick.NullLocator())
        ax.grid(axis='x', which='major', color='k', ls='--', alpha=0.4)
        ax.grid(axis='x', which='minor', color='k', ls='--', alpha=0.2)
            
        axlist.append(ax)
        i += 1
    
    if spec_return:
        xaxis = [spec_ax] + axlist
    else:
        xaxis = axlist
    
    [ax.set_xlim(time[0], time[-1]) for ax in xaxis]
    [ax.xaxis.set_major_locator(mtick.MaxNLocator(nbins=6)) for ax in xaxis]
    [ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(4)) for ax in xaxis]

    if (time[-1] - time[0]).total_seconds() > 86400:
        xaxis[0].set_title('%s--%s', time[0].strftime('%d %b %Y'), time[-1].strftime('%d %b %Y'))
    else:
        xaxis[0].set_title(time[0].strftime('%d %b %Y'))
    
    fig.align_ylabels()

    if full_return:
        return (axlist, spec_return, time, trace_ids)

    else:
        return fig
        

def init_network_gui(network, sta_list, starttime, delta, component, sde_file, init_app=True, **kwargs):
    
    if init_app:
        app = QtWidgets.QApplication([])
    
    NETframe = NetworkWindow(network, sta_list, starttime, delta, component, sde_file, **kwargs)
    NETframe.show()
    
    if init_app:
        app.exec_()
    else:
        return NETframe
        

