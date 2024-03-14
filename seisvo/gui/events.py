#!/usr/bin/env python3
# coding=utf-8

import pyqtgraph
import os
import datetime as dt
import numpy as np
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
from obspy.signal.invsim import estimate_magnitude
from seisvo.plotting import get_colors
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from .utils import Navigate, Picker


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

default_fqband = {
    "f1" : (),
    "f2" : (0.5,15),
    "f3" : (1,10),
    "f4" : (2,5),
    "f5" : (1,3),
}


def _plot_row(row, fig=None, off_sec=0, fq_band=default_fqband["f1"], focus=None):
    stream = row.get_stream(off_seconds=off_sec)

    if not stream:
        return None
    
    if not fig:
        fig = plt.figure(figsize=(12,9))
        return_fig = True
    else:
        return_fig = False

    default_sort = "ZNE"
    components = ''.join([tr.stats.channel[-1] for tr in stream if tr.stats.channel[-1] in ("Z", "N", "E")])
    diff = len(default_sort) - len(components)
    axes = fig.subplots(len(components), 1)

    if len(stream) == 1:
        axes = [axes]

    axes_dict = {}
    time = stream[0].get_time()
    for axn, trace in enumerate(stream):
        comp = trace.stats.channel[-1]
        if comp in ("Z", "N", "E"):
            n = default_sort.index(comp)-diff
            if n < 0:
                n = 0

            if focus:
                y = trace.get_data(detrend=True, norm=False, fq_band=fq_band)
                yfocus = trace.get_data(starttime=focus[0], endtime=focus[1], detrend=True, fq_band=fq_band)
                max_value = np.abs(yfocus).max()
                y /= max_value

            else:
                y = trace.get_data(detrend=True, norm=False, fq_band=fq_band)
                max_value = np.abs(y).max()
                y /= max_value

            axes_dict[trace.stats.channel] = axes[n]
            axes[n].plot(time, y, color=comp_colors[comp])
            axes[n].annotate(f"max count {max_value:.1f}", xy=(0,0.9), xycoords='axes fraction', color="k")
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

    return axes, phases, axes_dict


class _rowIter(object):
    def __init__(self, event, station_id=None):
        self.event  = event
        self.nrows  =len(event.rows_)

        self.n = 0 # by default take first

        if station_id:
            for n, row in enumerate(event):
                if row.get_station_id() == station_id:
                    self.n = n
                    break


    def eid(self):
        return self.event.id


    def rowid(self):
        return self.event[self.n].id


    def stations(self):
        return self.event.stations


    def row(self):
        return self.event[self.n]


    def station_id(self):
        return self.event[self.n].get_station_id()


    def to_staid(self, station_id):
        for n, row in enumerate(self.event):
            if row.get_station_id() == station_id:
                self.n = n


    def next(self):
        if self.n + 1 == self.nrows:
            self.n = 0
        else:
            self.n += 1


    def prev(self):
        if self.n - 1 < 0:
            self.n = self.nrows - 1
        else:
            self.n -= 1


    def update(self, event):
        self.event = event


class _eventIter(object):
    def __init__(self, sde, eid):
        self.sde  = sde
        self.elst = sde.get_event_list()
        self.neid = len(self.elst)
        
        # by default take the first event
        self.n    = 0
        if eid and eid in self.elst:
            self.n = self.elst.index(eid)


    def get_elst_str(self):
        return list(map(str, self.elst))


    def update(self):
        self.elst = self.sde.get_event_list()
        self.neid = len(self.elst)
        if self.n >= self.neid:
            self.n = self.neid - 1


    def event(self):
        return self.sde[self.elst[self.n]]

    
    def to_eid(self, eid):
        self.update()
        if eid in self.elst:
            self.n = self.elst.index(eid)


    def next(self):
        if self.n + 1 == self.neid:
            self.n = 0
        else:
            self.n += 1


    def prev(self):
        if self.n - 1 < 0:
            self.n = self.neid - 1
        else:
            self.n -= 1


    def __str__(self):
        return self.event().__str__()
    

class EventCanvas(FigureCanvas):
    def __init__(self, event, station_id, fb="f1", parent=None):

        # add Fig frame and linked with canvas
        self.parent = parent
        self.fig = Figure(figsize=(12,9), subplotpars=SubplotParams(left=0.08, right=0.92, wspace=0.1, top=0.95, bottom=0.05))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)

        # intit focus dict
        self.focus = {"active":False, "ticks":[]}
        self.fb = fb

        # init row iterator
        self.rowIter = _rowIter(event, station_id=station_id)
        self.plot()
        

    def reload_event(self):
        event = self.parent.sde[self.rowIter.event.id]
        self.rowIter.update(event)


    def plot(self):
        self.fig.clf()
        self.ticks = dict(chan="", right=None, left=None)
        self.reload_event()

        row = self.rowIter.row()
        
        with pyqtgraph.BusyCursor():

            self.axes_, phase_, self.axes_dict = _plot_row(row, fig=self.fig, fq_band=default_fqband[self.fb], focus=self.focus["ticks"])

            # load picker
            self.picker_ = Picker(self.axes_, phase_, self.parent.sde, row.id, self, phase_colors=phase_colors)

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


    def get_channel_from_axes(self, axes):
        for chan, ax in self.axes_dict.items():
            if ax == axes:
                return chan
        return None


    def on_click(self, event):
        if event.inaxes:
            if event.inaxes in self.axes_:
                self.ticks['chan'] = self.get_channel_from_axes(event.inaxes)
                t = mdates.num2date(float(event.xdata))
                t = t.replace(tzinfo=None)

                if event.button == 1:
                    self.ticks['left'] = t
                
                if event.button == 3:
                    self.ticks['right'] = t
                
                self.print_ticks()
    

    def change_row(self, arg):
        if arg == "up":
            self.rowIter.next()
            
        if arg == "down":    
            self.rowIter.prev()

        if arg in self.rowIter.stations():
            self.rowIter.to_staid(arg)


    def on_key(self, event):

        if event.key in list(default_fqband.keys()):
            if event.key != self.fb:
                print(f"\n  [info]  >>>  new freq. band  set to  {default_fqband[self.fb]}")
                self.fb = event.key
                self.plot()
            return


        if event.key =='escape':
            print("\n  <<< [Info]  EXIT of the Insert/Picker mode ")
            event = self.rowIter.event
            self.parent.setWindowTitle(f"Event {event.id} ({event.label})")
            self.parent.setCursor(QtCore.Qt.ArrowCursor)
            self.parent.setFocus()
            return
        

        if event.key in ("up", "down"):
            event = self.rowIter.event
            self.change_row(event.key)
            self.plot()
            return

        
        if event.key == "o":
            if not self.focus["active"]:

                if self.ticks['right'] and self.ticks['left']:
                    rtick = self.ticks["right"]
                    ltick = self.ticks["left"]

                    self.focus = {
                    "active":True,
                    "ticks":[min([rtick, ltick]), max([rtick, ltick])]
                    }
                    print("\n  <<< [Info]  entry FOCUS mode between ticks")
                    self.plot()

            else:
                self.focus = {
                "active":False,
                "ticks":[]
                }
                print("\n  <<< [Info]  exit FOCUS mode ")
                self.plot()
            
            return


        if event.key =='delete':
            change = False
            for wave in ["P", "S", "F"]:
                if self.picker_.phase[wave]["artist"]:
                    change = True
                    self.picker_.clear(wave)
            
            if change:
                self.picker_.save()
                self.draw()
            return
        

        if event.key == "tab":
            curr     = self.rowIter.station_id()
            sid , ok = QtWidgets.QInputDialog.getItem(self, 'Cambiar de StationID', 'Selecciona:', self.rowIter.stations(), 
                self.rowIter.stations().index(curr), editable=False)
            
            if ok and sid != curr:
                self.change_row(sid)
                self.plot()
            return


        if event.key == "+":
            self.rowIter.event.info()
            return


        if event.key == "a":
            if self.ticks['right'] and self.ticks['left']:
                p2p_amp, p2p_delta  = self.get_p2p(self.ticks['chan'], self.ticks['left'], self.ticks['right'])
                # print(f" Peak to peak [info]\n    Amplitude: {p2p_amp:.2f} [counts]\n    Time  :{p2p_delta:.2f} [sec]")
                self.parent.eventIter.sde._update_row(self.rowIter.rowid(), {"value_1":p2p_amp, "value_2":p2p_delta, "string_1":self.ticks['chan']})
                self.rowIter.event.info()
            return
        

        if event.key == "h":
            print("\n ------ INFO -----")
            print(" <-- / --> :: change station")
            print("    tab    :: select station")
            print("  f1---f5  :: change bandpass filter")
            print("     a     :: save peak2peak")
            print("    supr   :: clean all picked phases")
            print(" p/1/2/3/4 :: pick P phase")
            print("     P     :: clean P phase")
            print(" s/6/7/8/9 :: pick P phase")
            print("     S     :: clean P phase")
            print("    f/F    :: pick/clean F phase")
            print("    esc    :: exit PICKER mode")
            print("")
            return


    def get_p2p(self, channel, time1, time2):
        # compute peak to peak amplitude and period

        sta = self.rowIter.row().get_station()
        st  = sta.get_stream(
            starttime=min([time1, time2]), 
            endtime=max([time1, time2]), 
            channel=channel
            )
        data = st[0].data
        fs   = st[0].stats.sampling_rate

        # amplitude
        peak_max = abs(data.max())
        peak_min = abs(data.min())
        p2p_amp  = peak_min + peak_max

        # timespan
        t0 = np.argmin(data)
        t1 = np.argmax(data)
        p2p_time = 2*abs(t1-t0)/fs

        return p2p_amp, p2p_time


class EventWidget(QtWidgets.QWidget):
    def __init__(self, sde, event_id, station_id, parent=False):
        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()
        self.canvas = None
        self.parent = parent
        self.sde    = sde
        self.setFocus() 
        # by default the focus is on the Widget, 
        # to change the Focus of the Canvas, press key ``I''
        self.eventIter = _eventIter(sde, event_id)
        self.load_canvas(station_id)
        

    def load_canvas(self, station_id=None):
        # load event
        event = self.eventIter.event()
        
        # clear pre-existing canvas
        if self.canvas:
            # take station_id
            station_id = self.canvas.rowIter.station_id()
            self.layout.removeWidget(self.canvas)
            self.canvas.deleteLater()

        else:
            if station_id in event.stations:
                station_id = station_id
            else:
                station_id = event[0].get_station_id()
            
        # print event
        print(self.eventIter) 
        
        # init canvas
        self.setWindowTitle(f"Event {event.id} ({event.label})")
        self.canvas = EventCanvas(event, station_id, parent=self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
    

    def update(self, station_id=None, eid=None):
        self.eventIter.update()

        if eid:
            self.eventIter.to_eid(eid)
        
        if not station_id:
            station_id = self.canvas.rowIter.station_id()

        self.load_canvas(station_id=station_id)


    def change_event(self, event_key):
        if event_key == QtCore.Qt.Key_Right:
            self.eventIter.next()

        if event_key == QtCore.Qt.Key_Left:
            self.eventIter.prev()


    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key_H:
            print("\n ------ INFO -----")
            print("  <--/-->  :: change event")
            print("    ^/v    :: change station")
            print("    tab    :: select event")
            print("    supr   :: delete event")
            print("     t     :: move StationGUI to current event time")
            print("     i     :: entry PICKER mode")
            print("     r     :: relabel event")
            print("     w     :: write to hypo71 format")
            # print("     m     :: compute local magnitude")
            print("     +     :: print full info")
            print("")
            return


        if event.key() in [QtCore.Qt.Key_Right, QtCore.Qt.Key_Left]:
            self.change_event(event.key())
            self.load_canvas()
            return


        if event.key() in [QtCore.Qt.Key_Up, QtCore.Qt.Key_Down]:
            action = ["up", "down"][[QtCore.Qt.Key_Up, QtCore.Qt.Key_Down].index(event.key())]
            self.canvas.change_row(action)
            self.canvas.plot()
            return


        if event.key() == QtCore.Qt.Key_Tab:
            eid , ok = QtWidgets.QInputDialog.getItem(
                None, 'Cambiar de evento', 'Selecciona ID:',
                self.eventIter.get_elst_str(), self.eventIter.n, False)
            if ok:
                self.eventIter.to_eid(int(eid))
                self.load_canvas()
            return
        

        if event.key() == QtCore.Qt.Key_I:
            event = self.eventIter.event()
            print("\n  >>> [Info]  ENTRY on Insert/Picker mode ")
            self.setWindowTitle(f"[PICKER MODE] Event {event.id} ({event.label})")
            # self.setWindowTitle(f"[P]  Evento ID [{self.event.label}] {self.event.id}")
            self.setCursor(QtCore.Qt.CrossCursor)
            self.canvas.setFocus()
            return

           
        if event.key() == QtCore.Qt.Key_R:
            event = self.eventIter.event()
            text , ok = QtWidgets.QInputDialog.getText(
                None, f'Nueva etiqueta Evento', 'Nueva Etiqueta:', 
                text=event.label)
            if ok:
                event.relabel(text.upper())
            return
        

        if event.key() == QtCore.Qt.Key_Plus:
            self.eventIter.event().info()
            return
        

        if event.key() == QtCore.Qt.Key_W:
            event = self.eventIter.event()
            # create file output
            out  = f"location/{event.id}/phase.in"

            if not os.path.isdir(f"./location/{event.id}"):
                os.makedirs("./phase")

            fout = open(out, 'w')
            try:
                event.to_hypo71(phsfile=fout, show=False)
                fout.close()
                print(f" OK :: file {out} created!")
            
            except Exception as e:
                fout.close()
                print(e) 
                print(" FAIL :: error found")
            return


        if self.parent:
            if event.key() == QtCore.Qt.Key_Delete:
                ans = self.eventIter.event().id
                with pyqtgraph.BusyCursor():
                    self.parent._delete_event(ans)
                return

            if event.key() == QtCore.Qt.Key_R:
                event = self.eventIter.event()
                label, ok = QtWidgets.QInputDialog.getText(self, "Relabel event","Label: ", text=event.label)
                label = label.upper()
                if ok and label != event.label:
                    with pyqtgraph.BusyCursor():
                        self.parent._relabel_event(event.id, label)
                return
            
            if event.key() == QtCore.Qt.Key_T:
                event = self.eventIter.event()
                interval = self.parent.parent.interval.total_seconds()
                new_time = event.starttime - dt.timedelta(seconds=interval/2)
                with pyqtgraph.BusyCursor():
                    self.parent.parent.plot(new_time)
                return
            
            # if event.key() == QtCore.Qt.Key_M:
                # print("compute local magnitude")
                # estimate_magnitude(resp, p2p, T, distancia_Hipo)
                # return


    def closeEvent(self, event):
        if self.parent:
            self.parent.eGUI = None
        event.accept()


