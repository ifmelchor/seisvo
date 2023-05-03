#!/usr/bin/env python3
# coding=utf-8

import pyqtgraph
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
from seisvo.plotting import get_colors
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from .utils import Navigate, Picker, getYesNo


phase_colors = {
    'P':get_colors('okabe')[3],
    'S':get_colors('okabe')[0],
    'F':get_colors('okabe')[6]
}

comp_colors = {
    'Z':get_colors('zesty')[0],
    'E':get_colors('zesty')[1],
    'N':get_colors('zesty')[3]
}


def _plot_row(row, fig=None, off_sec=0, fq_band=(0.5,15)):
    stream = row.get_stream(off_seconds=off_sec)

    if not stream:
        return None
    
    if not fig:
        fig = plt.figure(figsize=(12,9))
        return_fig = True
    else:
        return_fig = False

    default_sort = "ZNE"
    components = ''.join([tr.stats.channel[-1] for tr in stream])
    diff = len(default_sort) - len(components)
    axes = fig.subplots(len(stream), 1)

    if len(stream) == 1:
        axes = [axes]

    time = stream[0].get_time()
    for axn, trace in enumerate(stream):
        comp = trace.stats.channel[-1]
        n = default_sort.index(comp)-diff
        if n < 0:
            n = 0
        y = trace.get_data(detrend=True, norm=True, fq_band=fq_band)
        axes[n].plot(time, y, color=comp_colors[comp])
        axes[n].set_xlim(time[0], time[-1])
        axes[n].set_ylim(-1.1, 1.1)
        axes[n].yaxis.set_major_formatter(mtick.NullFormatter())
        axes[n].set_ylabel(trace.stats.channel)
        axes[n].grid(axis='x',which='major',ls='--',color='k',alpha=0.2)
        # axes[n].annotate(txt, xy=(0,1.1), xycoords='axes fraction', color='k')

    axes[0].set_title(row.get_station_id())
    # build phase dictionary
    phases = {
        "P":{
            "time":None,
            "weight":None,
            "onset":None,
            "artist":[],
            "artist_text":None
        },
        "S":{
            "time":None,
            "weight":None,
            "artist":[],
            "artist_text":None
        },
        "F":{
            "time":None,
            "artist":[],
            "artist_text":None
        }
    }

    if row.time_P:
        phases['P']["time"] = row.time_P
        if row.weight_P:
            phases["P"]['weight'] = row.weight_P
        if row.onset_P:
            phases["P"]['onset'] = row.onset_P
        
    if row.time_S:
        phases['S']["time"] = row.time_S
        if row.weight_S:
            phases["S"]['weight'] = row.weight_S
        
    if row.time_F:
        phases['F']["time"] = row.time_F
    
    # draw phases
    for wave, phase in phases.items():
        if phase["time"]:
            for ax in axes:
                phase['artist'].append(ax.axvline(phase["time"], color=phase_colors[wave], lw=1.1))
            txt = Picker.phase_text(wave, phase)
            phase['artist_text'] = axes[0].annotate(txt, xy=(phase["time"], 1),color=phase_colors[wave])


    if return_fig:
        return fig

    return axes, phases



# def _plot_multirow(row_list, fig, off_sec, component=None, fq_band=()):



class _Canvas(FigureCanvas):
    def __init__(self, event, station_id):

        # add Fig frame and linked with canvas
        self.fig = Figure(figsize=(12,9),
            subplotpars=SubplotParams(left=0.08, right=0.92, wspace=0.1, top=0.95, bottom=0.05))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
        
        self.event = event
        self.max_row = len(event.rows_)
        
        try:
            for n, row in enumerate(self.event):
                if row.get_station_id() == station_id:
                    self.load_row(n)
        except:
            print(" warn :: station {station_id} not found")
            self.load_row(0)

        self.plot()


    def load_row(self, r):
        self.row_index  = r
        self.station_id = self.event[r].get_station_id()
        self.row        = self.event[r]


    def plot(self):
        self.fig.clf()
        self.ticks = dict(right=None, left=None)
        
        with pyqtgraph.BusyCursor():
            self.axes_, self.phase_ = _plot_row(self.row, fig=self.fig)
        
        # load picker
        self.picker_ = Picker(self.axes_, self.phase_, self.event.sde, self.row.id, self, phase_colors=phase_colors)

        # load navigation
        self.nav_ = Navigate(self.axes_, self, color='red', linewidth=0.5, alpha=0.5)

        self.draw()


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
                self.nav_.reset_ticks()
                
                t = mdates.num2date(float(event.xdata))
                t = t.replace(tzinfo=None)

                if event.button == 1:
                    self.ticks['left'] = t
                
                if event.button == 3:
                    self.ticks['right'] = t
                
                self.print_ticks()
    

    def on_key(self, event):
        if event.key =='up':
            r = self.row_index + 1

            if r == self.max_row:
                self.load_row(0)
            else:
                self.load_row(r)
            
            self.plot()
        

        if event.key =='down':
            r = self.row_index - 1

            if r < 0:
                self.load_row(self.max_row-1)
            else:
                self.load_row(r)
            
            self.plot()
        

        if event.key =='-':
            ans = getYesNo("Clear phases", "Estás seguro de que quieres eliminas las fases?")
            if ans.exec():
                self.picker_.clear()
                self.picker_.save()
                self.draw()


        if event.key == 'backspace':
            print("seleccionar station_id")
            sid , ok = QtWidgets.QInputDialog.getItem(None, 'Cambiar de StationID', 'Selecciona:', self.event.stations, self.station_id, False)
            if ok and sid != self.station_id:
                for n, row in enumerate(self.event):
                    if row.get_station_id() == sid:
                        self.load_row(n)
                self.plot()
        


class EventWidget(QtWidgets.QWidget):
    def __init__(self, sde, event_id, station_id):
        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()
        self.canvas = None

        self.sde = sde
        self.event_list = sde.get_event_list()
        self.max_index = len(self.event_list)

        if event_id not in self.event_list:
            self.load_event(self.event_list[0], station_id)
        else:
            self.load_event(event_id, station_id)
        

    def load_event(self, eid, station_id):
        self.event_idx = self.event_list.index(eid)
        self.event = self.sde[eid]

        if station_id in self.event.stations:
            self.station_id = station_id
        else:
            self.station_id = self.event[0].get_station_id()
        
        self.setWindowTitle(f"Evento ID :: {eid}")
        self.load_canvas()


    def load_canvas(self):
        if self.canvas: # clear pre-existing canvas
            self.layout.removeWidget(self.canvas)
            self.canvas.deleteLater()
        
        self.print_event()    
        self.canvas = _Canvas(self.event, self.station_id)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)


    def print_event(self):
        print(self.event)
    

    def print_event_phases(self):
        print(" --  PHASE EVENT info  -- ")
        print("     EID     |       TIME      |  NPHASE/NSTATION  ")
        for event in self.sde:
            time = event.starttime
            nsta = len(event)
            nphase = event.get_phases(nro_phases=True)
            txt = f"  {event.id:^13.0f}  {time:^17.0f}   {nphase}/{nsta} "
            print(txt)


    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_Right:
            idx = self.event_idx + 1
            if idx == self.max_index:
                idx = 0
            self.load_event(self.event_list[idx], self.station_id)


        if event.key() == QtCore.Qt.Key_Left:
            idx = self.event_idx - 1
            if idx < 0:
                idx = self.max_index-1
            self.load_event(self.event_list[idx], self.station_id)
        

        if event.key() == QtCore.Qt.Key_Tab:
            eid , ok = QtWidgets.QInputDialog.getItem(None, 'Cambiar de evento', 'Selecciona ID:', self.event_list, self.event.id, False)
            if ok and int(eid) != self.event.id:
                self.load_event(eid, self.station_id)


        if event.key() == QtCore.Qt.Key_Delete:
            print("DELETE event")

        
        if event.key() == QtCore.Qt.Key_Insert:
            print("CLONE event")
           

        if event.key() == QtCore.Qt.Key_R:
            text , ok = QtWidgets.QInputDialog.getText(None, f'Nueva etiqueta Evento {self.event.id}', 'Nueva Etiqueta:', text=self.event.label)
            if ok:
                self.event.relabel(text.upper())
        

        if event.key() == QtCore.Qt.Key_Plus:
            # Print event phases
            self.print_event_phases()


        if event.key() == QtCore.Qt.Key_Enter:
            # Print event info
            self.print_event()
        
         
