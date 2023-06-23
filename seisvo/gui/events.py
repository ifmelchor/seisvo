#!/usr/bin/env python3
# coding=utf-8

import pyqtgraph
import numpy as np
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

        if focus:
            y = trace.get_data(detrend=True, norm=False, fq_band=fq_band)
            yfocus = trace.get_data(starttime=focus[0], endtime=focus[1], detrend=True, fq_band=fq_band)
            max_value = np.abs(yfocus).max()
            y /= max_value

        else:
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



class EventCanvas(FigureCanvas):
    def __init__(self, event, station_id, fb="f1", parent=None):

        # add Fig frame and linked with canvas
        self.parent = parent
        self.fig = Figure(figsize=(12,9), subplotpars=SubplotParams(left=0.08, right=0.92, wspace=0.1, top=0.95, bottom=0.05))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
        self.event = event
        self.max_row = len(event.rows_)
        self.fb = fb

        self.focus = {
        "active":False,
        "ticks":[]
        }
        
        try:
            for n, row in enumerate(self.event):
                if row.get_station_id() == station_id:
                    self.load_row(n)
        except:
            print(" warn :: station {station_id} not found")
            self.load_row(0)


    def load_row(self, r):
        self.row_index  = r
        self.station_id = self.event[r].get_station_id()
        self.row        = self.event[r]
        self.plot()


    def reload_event(self):
        self.event = self.parent.sde[self.event.id]


    def plot(self):
        self.fig.clf()
        self.ticks = dict(right=None, left=None)
        self.reload_event()
        
        with pyqtgraph.BusyCursor():

            self.axes_, phase_ = _plot_row(self.row, fig=self.fig, fq_band=default_fqband[self.fb], focus=self.focus["ticks"])

            # load picker
            self.picker_ = Picker(self.axes_, phase_, self.event.sde, self.row.id, self, phase_colors=phase_colors)

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
                # self.nav_.reset_ticks()
                
                t = mdates.num2date(float(event.xdata))
                t = t.replace(tzinfo=None)

                if event.button == 1:
                    self.ticks['left'] = t
                
                if event.button == 3:
                    self.ticks['right'] = t
                
                self.print_ticks()
    

    def change_row(self, arg):

        if arg == "up":
            r = self.row_index + 1
            if r == self.max_row:
                self.load_row(0)
            else:
                self.load_row(r)
        

        if arg == "down":    
            r = self.row_index - 1
            if r < 0:
                self.load_row(self.max_row-1)
            else:
                self.load_row(r)


        if arg in self.event.stations:
            for n, row in enumerate(self.event):
                if row.get_station_id() == arg:
                    self.load_row(n)


    def on_key(self, event):
        # print(event.key)

        if event.key in list(default_fqband.keys()):
            if event.key != self.fb:
                print(f"\n  [info]  >>>  new freq. band  set to  {default_fqband[self.fb]}")
                self.fb = event.key
                self.plot()


        if event.key =='escape':
            print("\n  <<< [Info]  EXIT of the Insert/Picker mode ")
            self.parent.setWindowTitle(f"Evento ID :: {self.event.id}")
            self.parent.setCursor(QtCore.Qt.ArrowCursor)
            self.parent.setFocus()
        

        if event.key in ("up", "down"):
            self.change_row(event.key)

        
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

        #     key = [QtCore.Qt.Key_Right, QtCore.Qt.Key_Left][["right", "left"].index(event.key)]
        #     self.parent.change_event(key)
        

        if event.key =='delete':
            change = False
            for wave in ["P", "S", "F"]:
                if self.picker_.phase[wave]["artist"]:
                    change = True
                    self.picker_.clear(wave)
            
            if change:
                self.picker_.save()
                self.draw()
        

        if event.key == "tab":
            sid , ok = QtWidgets.QInputDialog.getItem(None, 'Cambiar de StationID', 'Selecciona:', self.event.stations, self.event.stations.index(self.station_id), editable=False)
            if ok and sid != self.station_id:
                self.change_row(sid)

        if event.key == "+":
            self.parent.print_event_phases()
        

        if event.key == " ":
            print("\n ------ INFO -----")
            print(" <-- / --> :: cambiar estación")
            print("    tab    :: elige estación")
            print("  f1---f5  :: cambia filtro")
            print("    supr   :: elimina picks")
            print(" p/1/2/3/4 :: pica fase P")
            print("     P     :: elimina fase P")
            print(" s/6/7/8/9 :: pica fase S")
            print("     S     :: elimina fase S")
            print("    f/F    :: pica/elimina fase F")
            print("    esc    :: salir del modo PICKER")
            print("")


class EventWidget(QtWidgets.QWidget):
    def __init__(self, sde, event_id, station_id):
        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()
        self.canvas = None
        self.setFocus() 
        # by default the focus is on the Widget, to change the Focus of the Canvas, press key I(nsert)

        self.sde = sde
        self.event_list = sde.get_event_list()
        self.max_index = len(self.event_list)

        if event_id not in self.event_list:
            self.load_event(self.event_list[0], station_id)
        else:
            self.load_event(event_id, station_id)
        

    def load_event(self, eid, station_id):
        # self.idx = eid
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
        self.canvas = EventCanvas(self.event, self.station_id, parent=self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)


    def print_event(self):
        print(self.event)


    def change_event(self, event_key):
        if event_key == QtCore.Qt.Key_Right:
            idx = self.event_idx + 1
            if idx == self.max_index:
                idx = 0
            self.load_event(self.event_list[idx], self.station_id)

        if event_key == QtCore.Qt.Key_Left:
            idx = self.event_idx - 1
            if idx < 0:
                idx = self.max_index-1
            self.load_event(self.event_list[idx], self.station_id)

    

    def print_event_phases(self):
        print(f" --  PHASE EVENT [{self.event.id}] info  -- ")
        print("   STA   |   P   |   S   |   P-S   ")
        for row in self.sde[self.event.id]:
            text = f"{row.station:^9}"

            if row.time_P:
                p_t = (row.time_P - self.event.starttime).total_seconds()
                text += f" {p_t:^7.2f}"
            else:
                p_t = None
                text += f" {'-':^7}"

            if row.time_S:
                s_t = (row.time_S - self.event.starttime).total_seconds()
                text += f" {s_t:^7.2f}"
            else:
                s_t = None
                text += f" {'-':^7}"

            if p_t and s_t:
                slp = s_t - p_t
                text += f" {slp:^7.2f}"
            else:
                slp = None
                text += f" {'-':^7}"

            print(text)


    def keyPressEvent(self, event):
        # print(event.key())

        if event.key() in [QtCore.Qt.Key_Right, QtCore.Qt.Key_Left]:
            self.change_event(event.key())


        if event.key() in [QtCore.Qt.Key_Up, QtCore.Qt.Key_Down]:
            action = ["up", "down"][[QtCore.Qt.Key_Up, QtCore.Qt.Key_Down].index(event.key())]
            self.canvas.change_row(action)
        

        if event.key() == QtCore.Qt.Key_Tab:
            eid , ok = QtWidgets.QInputDialog.getItem(None, 'Cambiar de evento', 'Selecciona ID:', list(map(str, self.event_list)), self.event_idx, False)
            if ok and int(eid) != self.event.id:
                self.load_event(int(eid), self.station_id)


        if event.key() == QtCore.Qt.Key_Delete:
            print("DELETE event")
        

        if event.key() == QtCore.Qt.Key_I:
            print("\n  >>> [Info]  ENTRY on Insert/Picker mode ")
            self.setWindowTitle(f"[P]  Evento ID :: {self.event.id}")
            self.setCursor(QtCore.Qt.CrossCursor)
            self.canvas.setFocus()

        
        if event.key() == QtCore.Qt.Key_W:
            print("Create new Event")
            # if self.canvas.ticks['right'] and self.canvas.ticks['left']:
            #     rtick = self.canvas.ticks["right"]
            #     ltick = self.canvas.ticks["left"]
                
            #     end = max([rtick, ltick])
            #     start = min([rtick, ltick])

            #     # ask for label?

            #     # save into database
            #     self.sde.add_event()

                # reload events

                # notify 

           

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

        
        if event.key() == QtCore.Qt.Key_Space:
            print("\n ------ INFO -----")
            print(" <-- / --> :: cambiar evento")
            print("    tab    :: elige evento")
            print("    supr   :: elimina evento")
            print("     i     :: entra en modo PICKER")
            print("    ins    :: clonar evento")
            print("     r     :: re-etiqueta evento")
            print("     +     :: info fases")
            print("   enter   :: info evento")
            print("")
        
        