#!/usr/bin/env python3
# coding=utf-8

from seisvo.plotting import get_colors
from seisvo.plotting.gui import Navigation
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas

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


def _phasetext(phase_dict, wave):

    if wave == "P":



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
        axes[n].yaxis.set_major_formatter(mtick.NullFormatter)
        axes[n].set_ylabel(trace.stats.channel)
        axes[n].grid(axis='x',which='major',ls='--',color='k',alpha=0.2)
        # axes[n].annotate(txt, xy=(0,1.1), xycoords='axes fraction', color='k')

    # draw phases
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
            "artist":None,
            "artist_text":None
        },
        "F":{
            "time":None,
            "artist":None,
            "artist_text":None
        }
    }

    if row.time_P:
        phases['P']["time"] = row.get_Pdate()
        if row.weight_P:
            phases["P"]['weight'] = row.weight_P
        if row.onset_P:
            phases["P"]['onset'] = row.onset_P
        for ax in axes:
            phases["P"]['artist'].append(ax.axvline(phases['P']["time"], color=phase_colors["P"], lw=1.1))
        txt = _phasetext(phases["P"], wave="P")
        phases["P"]['artist_text'] = axes[0].annotate(txt)

    if row.time_S:
        phases['S']["time"] = row.get_Sdate()
        if row.weight_S:
            phases["S"]['weight'] = row.weight_S
        for ax in axes:
            phases["S"]['artist'].append(ax.axvline(phases['S']["time"], color=phase_colors["S"], lw=1.1))
        txt = _phasetext(phases["S"], wave="S")
        phases["S"]['artist_text'] = axes[0].annotate(txt)
    
    if row.time_F:
        phases['F']["time"] = row.time_F
        for ax in axes:
            phases["F"]['artist'].append(ax.axvline(phases['F']["time"], color=phase_colors["F"], lw=1.1))
        phases["F"]['artist_text'] = axes[0].annotate("F")

    if return_fig:
        return fig

    return axes, phases


class _Canvas(FigureCanvas):
    def __init__(self, event, station_id=None, parent=None):
        
        self.parent = parent
        self.event = event
        
        self.max_row = len(event.rows_)
        if not station_id:
            self.load_row(0)
        else:
            for n, row in enumerate(self.event):
                if row.get_station_id() == station_id:
                    self.load_row(n)
        
        self.load_frame()
        self.plot()


    def load_frame(self):
        self.fig = Figure(figsize=(12,9),
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
        self.callbacks.connect('key_press_event', self.on_key)


    def load_row(self, r):
        self.row_index  = r
        self.station_id = self.event[r].get_station_id()
        self.row        = self.event[r]


    def plot(self):
        self.fig.clf()
        self.ticks = dict(right=[], left=[])
        with pyqtgraph.BusyCursor():
            self.axes, self.phases = _plot_row(self.row, fig=self.fig)
            self.nav  = [Navigation(ax, parent=self) for ax in self.axes]

        self.draw()
    

    




    # def on_click(self):
    # def on_key(self):
        

