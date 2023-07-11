#!/usr/bin/env python3
# coding=utf-8


import pyqtgraph
import datetime as dt

# matplotlib
import matplotlib.dates as mdates
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.text import Annotation

# seisvo
from ..utils import in_interval
from ..plotting.lte import LTESTAplot
from .utils import Navigate, notify
from ..database import LDE


class DrawEpisodes:
    def __init__(self, lde, axes, interval):
        self.lde = lde
        self.axes = axes
        self.delta = interval.total_seconds()/60
        self.edict = {}


    def draw(self, starttime, endtime):
        eid_list = self.lde.get_episode_list(time_interval=(starttime, endtime))

        if eid_list:
            self.edict = {}

            for eid in eid_list:
                e = self.lde[eid] 
                self.edict[eid] = []

                dur = (e.starttime+dt.timedelta(minutes=e.duration*0.35) - starttime).total_seconds()/60
                x_pos = dur/self.delta
                for n, ax in enumerate(self.axes):
                    a1  = ax.axvline(e.starttime,  alpha=0.5, color="g", ls='dashed')
                    a2  = ax.axvline(e.endtime,  alpha=0.5, color="g", ls='dashed')
                    a3  = ax.axvspan(e.starttime,  e.endtime, alpha=0.15, color="g")
                    self.edict[eid] += [a1, a2, a3]
                    
                    if n == 0:
                        a4  = ax.annotate(f"{e.label.upper()} [{e.id}]", (x_pos,1.05), xycoords="axes fraction")
                        self.edict[eid].append(a4)
                        
    
    def draw_eid(self, eid, axis_starttime):
        if self.lde.is_eid(eid):
            e = self.lde[eid]
            self.edict[eid] = []

            dur = (e.starttime+dt.timedelta(minutes=e.duration*0.35) - axis_starttime).total_seconds()/60
            x_pos = dur/self.delta
            for n, ax in enumerate(self.axes):
                a1  = ax.axvline(e.starttime,  alpha=0.5, color="g", ls='dashed')
                a2  = ax.axvline(e.endtime,  alpha=0.5, color="g", ls='dashed')
                a3  = ax.axvspan(e.starttime,  e.endtime, alpha=0.15, color="g")
                self.edict[eid] += [a1, a2, a3]
                
                if n == 0:
                    a4  = ax.annotate(f"{e.label.upper()} [{e.id}]", (x_pos,1.05), xycoords="axes fraction")
                    self.edict[eid].append(a4)


    def clear(self, eid=None):
        if self.edict:
            if eid:
                for id, artist_list in self.edict.items():
                    if id == eid:
                        [a.remove() for a in artist_list]
                del self.edict[eid]                    
            
            else:
                for id, artist_list in self.edict.items():
                    [a.remove() for a in artist_list]
                self.edict = {}


    def detect(self, tick_date):
        eid_list = []
        if self.edict:
            for eid in list(self.edict.keys()):
                e = self.lde[eid]
                if in_interval(tick_date, tick_date, (e.starttime, e.endtime)):
                    eid_list.append(eid)
        return eid_list


    def relabel(self, eid):
        e = self.lde[eid]
        for a in self.edict[eid]:
            if isinstance(a, Annotation):
                a.set_text(f"{e.label.upper()} [{e.id}]")


class LTEWidget(QtWidgets.QWidget):
    """
    Genera el espacio físico donde se va a alojar el CANVAS (LTECANVAS) y controla los botones
    """
    def __init__(self, lte, lde, starttime, interval, attr_list, chan_list, olap=0.1, **kwargs):
        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()
        
        if self.lde:
            self.lde = LDE(lde)
        else:
            self.lde
        self.lte = lte 
        self.starttime = starttime
        self.interval  = dt.timedelta(hours=interval)
        self.olap      = dt.timedelta(hours=interval*olap)
        self.attr_list = attr_list
        self.chan_list = chan_list
        self.fig_kwargs = kwargs

        # add canvas
        self.canvas = LTECanvas(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        self.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.setFocus()
    

    def on_key(self, key):
        if key == QtCore.Qt.Key_Right:
            if self.starttime + self.interval > self.lte.endtime:
                notify("LTE", "The bound of the file was reached!")
                endtime = self.lte.endtime
                self.starttime = endtime - self.interval
            else:
                endtime = None
                self.starttime = self.starttime + self.interval
            self.canvas.plot(endtime=endtime)

        if key == QtCore.Qt.Key_Left:
            if self.starttime - self.interval + self.olap < self.lte.starttime:
                notify("LTE", "The bound of the file was reached!")
                self.starttime = self.lte.starttime
            else:
                self.starttime = self.starttime - self.interval + self.olap
            self.canvas.plot()


    def relabel_event(self, eid, net_label):
        if self.lde:
            self.lde.relabel_event(eid, net_label)
            notify("LTE", f"Episode {eid} relabed!")


    def save_episode(self, episode_dict):
        if self.lde:
            eid = self.lde.add_episode(episode_dict)
            notify("LTE", f"Episode {eid} save in database")
            return eid


    def remove_episode(self, eid):
        if self.lde:
            self.lde.remove_event(eid)
            notify("LTE", f"Episode {eid} removed!")


class LTECanvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.fig = Figure(figsize=(12,9), subplotpars=SubplotParams(left=0.08,\
            right=0.92, wspace=0.1, top=0.95, bottom=0.05))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
        self.plot()


    def plot(self, endtime=None):
        self.fig.clf()
        self.ticks = dict(right=None, left=None)

        if not endtime:
            self.endtime = self.parent.starttime + self.parent.interval
        else:
            self.endtime = endtime
        
        with pyqtgraph.BusyCursor():

            out = self.parent.lte.get(self.parent.starttime, self.endtime,\
                attr=self.parent.attr_list, chan=self.parent.chan_list,\
                db_scale=True, azimuth_ambiguity=True)
            
            _, self.axes = LTESTAplot(out, self.parent.chan_list, self.parent.attr_list,\
                fig=self.fig, return_stats=False, plot=False, **self.parent.fig_kwargs)
            
            # add navigation
            self.axes_   = [ax[0] for _, ax in self.axes.items() if len(ax) == 1]
            imaxes_list = [ax[0] for _, ax in self.axes.items() if len(ax) == 2]
            self.nav = Navigate(self.axes_, self, im_axis=imaxes_list, color='red',\
            linewidth=0.5, alpha=0.5)

            # get episodes draw object
            if self.parent.lde:
                self.draw_ = DrawEpisodes(self.parent.lde, self.axes_, self.parent.interval)
                self.draw_.draw(self.parent.starttime, self.endtime)
            else:
                self.draw_ = None

            self.draw()

    
    def print_ticks(self):
        print("\n :: Ticks info :: ")
        for position, time in self.ticks.items():
            if time:
                print(f"  {position} >> {time} ")
        
        if self.ticks['right'] and self.ticks['left']:
            dist = abs((self.ticks['left'] - self.ticks['right']).total_seconds()/60)
            print(f" distance |R-L|  >>  {dist}  [min]")


    def on_click(self, event):
        if event.inaxes:
            if event.inaxes in self.axes_:
                
                t = mdates.num2date(float(event.xdata))
                t = t.replace(tzinfo=None)

                if event.button == 1:
                    self.ticks['left'] = t
                
                if event.button == 3:
                    self.ticks['right'] = t
                
                self.print_ticks()
        

    def on_key(self, event):
        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)
        

        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)
        

        if event.key == "s":
            if self.parent.lde:
                if self.ticks['right'] and self.ticks['left']:
                    starttime = min(self.ticks['right'], self.ticks['left'])
                    duration  = abs((self.ticks['left'] - self.ticks['right']).total_seconds()/60)
                    (net_code, sta_code, loc) = self.parent.lte.stats.id.split(".")
                    label , ok = QtWidgets.QInputDialog.getText(None,\
                        'Nuevo Episodio', 'Define Etiqueta:', text="TR")
                    if ok:
                        edict = {
                            "network":net_code,
                            "station":sta_code,
                            "location":loc,
                            "starttime":starttime,
                            "duration":duration,
                            "label":label.upper(),
                            "lte_file":self.parent.lte.stats.file
                        }
                        new_eid = self.parent.save_episode(edict)
                        self.draw_.draw_eid(new_eid, self.parent.starttime)
                        self.draw()


        if event.key in ("delete", "l"):
            if self.parent.lde:
                if self.draw_.edict:
                    tick_date = mdates.num2date(event.xdata).replace(tzinfo=None)
                    eid_list = self.draw_.detect(tick_date)
                    
                    if len(eid_list) == 1:
                        if event.key == "delete":
                            self.parent.remove_episode(eid_list[0])
                            self.draw_.clear(eid=eid_list[0])
                        else:
                            e = self.parent.lde[eid_list[0]]
                            label , ok = QtWidgets.QInputDialog.getText(None,\
                        f'Episodio {e.id}', 'Define Etiqueta:', text=f"{e.label}")
                            if ok:
                                self.parent.relabel_event(eid_list[0], label)
                                self.draw_.relabel(eid_list[0])
                                self.draw()
                    else:
                        # sid , ok = QtWidgets.QInputDialog.getItem(None, 'Cambiar de StationID',\
                        #     'Selecciona:', self.event.stations, self.event.stations.index(self.station_id),\
                        #     editable=False)
                        print(eid_list)
                    
                    self.draw()

