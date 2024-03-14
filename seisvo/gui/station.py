#!/usr/bin/env python
# coding=utf-8


from obspy import UTCDateTime
from obspy.imaging.waveform import WaveformPlotting
import datetime as dt
import os
import numpy as np
import matplotlib.gridspec as plg
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
import pyqtgraph
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from .utils import Navigate, DateDialog, SimplePlotDialog
from .events import EventWidget
from ..database import SDE
from ..utils import HILoc

EVENTS_KEYS = {
    "1":"VT",
    "2":"LP",
    "3":"RG",
    "4":"TR",
    "5":"VLP",
}


def _plot_stream(fig, station, starttime, endtime, spec_v, wlen, **kw_stream):
    """[summary]

    Args:
        station: Station object
        starttime: datetime
        channel (string): the channel to plot
        delta (int): step time
        return_fig (bool, optional): Return fig and axes. Defaults to False.
    """

    # get stream and trace
    stream = station.get_stream(starttime, endtime, **kw_stream).select2(component="Z")

    # clean figure and add axes
    fig.clf()
    gs = plg.GridSpec(2, 2, wspace=0.05, width_ratios=[1,0.025])
    ax_trac = fig.add_subplot(gs[0,0])
    ax_spec = fig.add_subplot(gs[1,0])
    ax_spec_cbar = fig.add_subplot(gs[1,1])

    # check if there are data
    trace  = stream[0]

    # plot trace
    data = trace.plot_trace(ax_trac)
    
    # plot specgram
    fq_band = kw_stream.get("fq_band", ())
    trace.plot_specgram(ax_spec, window_length=wlen, fq_band=fq_band, date_list=trace.get_time(),
        axis_bar=ax_spec_cbar, axis_bar_label='PSD', y_label="Freq. [Hz]",
        x_label="UTC Time", v_max=spec_v[1], v_min=spec_v[0])

    # define ticks
    ax_trac.xaxis.set_minor_locator(mtick.AutoMinorLocator(4))
    ax_trac.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    fig.suptitle(starttime.strftime("%d %B %Y"))

    ax_spec.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
    ax_spec.xaxis.set_minor_locator(mtick.AutoMinorLocator(4))
    ax_spec.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

    fig.align_ylabels()

    return (ax_trac, ax_spec)


def _day_plot(fig, station, starttime, endtime, events, fq_band=(), **kwargs):

    # clean fig
    fig.clf()

    # get stream
    stream = station.get_stream(starttime, endtime, fq_band=fq_band, sample_rate=10).select2(component="Z")

    # configure plot
    kwargs["stream"] = stream[0]
    kwargs["fig"]    = fig
    kwargs["show"]   = False
    kwargs["type"]   = "dayplot"
    kwargs["tick_format"] = "%d/%m/%Y \n %H:%M:%S"

    wp  = WaveformPlotting(**kwargs)

    # draw
    wp.plot_waveform()

    for event in events:
        wp._plot_event(event)

    return fig.axes[0]


def _net_ztrace(network, starttime, endtime, fq_band=()):
    
    fig  = Figure(figsize=(9,4))
    axis = fig.add_subplot(1, 1, 1)

    # get stream
    stream = network.get_stream(starttime, endtime, component="Z", fq_band=fq_band, norm=True)
    arrayz  = stream.to_array(norm=True)
    ncom    = arrayz.shape[0]
    arrayz += 2*np.array(list(range(ncom))).reshape(ncom,1)

    # get time
    start = min([stream[n].stats.starttime for n in range(ncom)]).datetime
    end   = max([stream[n].stats.endtime   for n in range(ncom)]).datetime
    # delta = dt.timedelta(seconds=1/stream[0].stats.sampling_rate)
    time  = np.linspace(mdates.date2num(start), mdates.date2num(end), arrayz.shape[1])
    time  = mdates.num2date(time)

    labels = [f"{stream[n].stats.station}.{stream[n].stats.channel}" for n in range(ncom)]

    axis.plot(time, arrayz.T, color="k", lw=0.8)
    axis.set_xlim(time[0],time[-1])
    axis.grid(axis='x', which='major', color='k', ls='--', alpha=0.4)
    axis.grid(axis='x', which='minor', color='k', ls='--', alpha=0.2)
    axis.xaxis.set_minor_locator(mtick.AutoMinorLocator(4))
    axis.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    axis.set_title(start.strftime("%d %B %Y"))
    axis.set_yticks(2*np.array(list(range(ncom))))
    axis.set_yticklabels(labels)

    return fig


def _zne_polar(stream, win, rl_th=0.7):

    with pyqtgraph.BusyCursor():
        ans = stream.time_polar(win)

    if ans:
        timePolar, rl, baz, inc = ans

        fig  = Figure(figsize=(9,9))
        gs   = plg.GridSpec(4, 1, wspace=0.5, height_ratios=[1,0.75,0.75,0.75])
        lax1 = fig.add_subplot(gs[0,0])
        lax2 = fig.add_subplot(gs[1,0], sharex=lax1) # Rl
        lax3 = fig.add_subplot(gs[2,0], sharex=lax1) # AZ
        lax4 = fig.add_subplot(gs[3,0], sharex=lax1) # INC
        axes = [lax1, lax2, lax3, lax4]

        # plot ZNE
        start = stream[0].stats.starttime.datetime
        end   = stream[0].stats.endtime.datetime
        duration  = (end - start).total_seconds()
        timeWvfm  = np.linspace(0, duration, stream[0].stats.npts)
        
        Z = stream.select(component="Z")[0].data
        Zmax = np.abs(Z).max()
        N = stream.select(component="N")[0].data
        Nmax = np.abs(N).max()
        E = stream.select(component="E")[0].data
        Emax = np.abs(E).max()

        lax1.plot(timeWvfm, Z/Zmax, color="k", lw=0.8)
        lax1.plot(timeWvfm, N/Nmax - 2, color="r", lw=0.8)
        lax1.plot(timeWvfm, E/Emax - 4, color="g", lw=0.8)
        lax1.set_ylim(-6,2)
        lax1.set_yticks([-4, -2, 0])
        lax1.set_yticklabels(["E", "N", "Z"])
        lax1.xaxis.set_major_formatter(mtick.NullFormatter())

        rrl1 = np.where(rl >= rl_th)
        rrl0 = np.where(rl < rl_th)
        lax2.scatter(timePolar[rrl1], rl[rrl1], ec="k", c="w")
        lax2.scatter(timePolar[rrl0], rl[rrl0], ec="k", c="w", alpha=0.5)
        lax2.set_ylim(0,1)
        lax2.axhline(rl_th, color="r", ls="--", lw=0.5, alpha=0.5)
        lax2.set_ylabel(r"R$_L$")
        lax2.xaxis.set_major_formatter(mtick.NullFormatter())

        # baz
        lax3.scatter(timePolar[rrl1], baz[rrl1], ec="k", c="w")
        lax3.scatter(timePolar[rrl0], baz[rrl0], ec="k", c="w", alpha=0.5)
        lax3.set_ylim(0,180)
        lax3.set_ylabel("PAZ")
        lax3.set_yticks([45,90,135])
        lax3.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        lax3.xaxis.set_major_formatter(mtick.NullFormatter())

        # inc
        lax4.scatter(timePolar[rrl1], inc[rrl1], ec="k", c="w")
        lax4.scatter(timePolar[rrl0], inc[rrl0], ec="k", c="w", alpha=0.5)
        lax4.set_ylim(-5,95)
        lax4.set_yticks([0,45,90])
        lax4.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        lax4.set_ylabel("PIN")
        lax4.set_xlabel("Time [s]")

        for ax in axes:
            ax.set_xlim(0, duration)
            ax.grid(color="k", ls='--', alpha=0.5)

        return fig

    else:
        return None
    # return fig


class EventDrawer(object):
    def __init__(self, parent):
        self.eventList = []
        self.eventDict = {}
        self.parent    = parent


    def export(self, time_interval):
        eventList = self.parent.parent.db.get_event_list(time_interval=time_interval)
        expList   = []

        for eid in eventList:
            event = self.parent.parent.db[eid]
            if time_interval[0] <= event.starttime <= time_interval[1]:
                event_dict = {
                    "time": UTCDateTime(event.starttime),
                    "text": f"{event.label}"}
                expList.append(event_dict)

        return expList


    def add_event(self, eid):
        self.eventList.append(eid)
        event = self.parent.parent.db[eid]

        t0  = event.starttime
        t1  = event.endtime

        span = self.parent.axes[0].axvspan(t0, t1, color="k", alpha=0.2, zorder=1)
        tp   = t0 + dt.timedelta(seconds=event.duration/2)
        yl   = self.parent.axes[0].get_ylim()
        txt  = self.parent.axes[0].annotate('%s(%s)' % (eid, event.label), (tp, yl[1]-0.1*abs(yl[1]-yl[0])), color='r', weight="bold", fontfamily="monospace", stretch='condensed')

        self.eventDict[eid] = {
            "event": event,
            "span": span,
            "text": txt
        }

        self.parent.draw()


    def update(self):
        self.eventList = self.parent.parent.db.get_event_list(time_interval=(self.parent.parent.starttime, self.parent.endtime))
        for eid in self.eventList:
            if self.eventDict.get(eid, None):
                self.eventDict[eid]["event"] = self.parent.parent.db[eid]


    def draw(self):
        if self.eventList:
            starttime = self.parent.parent.starttime
            endtime   = self.parent.endtime
            axis      = self.parent.axes[0]

            # event database
            self.eventDict = {}
            for eid in self.eventList:
                event = self.parent.parent.db[eid]

                st = event.starttime
                et = event.endtime
                
                # start and end in
                cond_in = st >= starttime and et <= endtime

                # start in, end out
                cond_st_in = endtime > st > starttime and et > endtime
                
                # start out, end in
                cond_et_in = starttime > st and starttime < et < endtime

                if cond_in:
                    t0 = st
                    t1 = et

                elif cond_st_in:
                    t0 = st
                    t1 = endtime

                elif cond_et_in:
                    t0 = starttime
                    t1 = et

                # cond_full_in
                else:
                    t0 = starttime
                    t1 = endtime

                span = axis.axvspan(t0, t1, color="k", alpha=0.2, zorder=1)

                tp = t0 + dt.timedelta(seconds=event.duration/2)
                yl = axis.get_ylim()
                txt  = axis.annotate('%s(%s)' % (eid, event.label), (tp, yl[1]-0.1*abs(yl[1]-yl[0])), color='r', weight="bold", fontfamily="monospace", stretch='condensed')

                self.eventDict[eid] = {
                    "event": event,
                    "span": span,
                    "text": txt
                }

        self.parent.draw()


    def get_eid(self, time):
        if self.eventDict:
            for eid, dd in self.eventDict.items():
                event = dd["event"]
                if event.starttime<=time<=event.endtime:
                    return eid

        return None


    def remove_event(self, eid):
        self.eventDict[eid]["span"].remove()
        self.eventDict[eid]["text"].remove()
        del self.eventDict[eid]
        self.eventList.remove(eid)
        self.parent.draw()


    def relabel_event(self, eid):
        event = self.parent.parent.db[eid]
        self.eventDict[eid][event] = event
        self.eventDict[eid]["text"].set_text('%s(%s)' % (eid, event.label))
        self.parent.draw()


class StationMainWidget(QtWidgets.QWidget):
    def __init__(self, network, station, starttime, interval, olap, sde_db, hyp, spec_v, spec_wlen=None, **sta_kwargs):

        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()

        self.network    = network
        self.station    = station
        self.sta_kwargs = sta_kwargs
        self.firsttime  = station.stats.starttime
        self.lasttime   = station.stats.endtime

        # load database
        self.db    = None
        self.is_db = False
        self.hyp   = None

        if isinstance(sde_db, str):
            self.db    = SDE(sde_db)
            self.is_db = True
            if isinstance(hyp, HILoc):
                self.hyp = hyp

        self.olap_pct  = olap
        self.starttime = starttime
        self.spec_v    = spec_v

        if not spec_wlen:
            spec_wlen = 50

        self.spec_wlen = spec_wlen
        self.interval  = dt.timedelta(minutes=interval)
        self.olap      = dt.timedelta(minutes=interval*self.olap_pct)

        # add canvas
        self.canvas = StationMainCanvas(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        self.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.setFocus()


    def all_stations(self, starttime, endtime):
        st  = self.network.get_stream(starttime, endtime, component=None)
        return st.get_stations_id()


    def on_key(self, key):
        if key == QtCore.Qt.Key_Right:
            if self.starttime + self.interval > self.lasttime:
                notify("SEISVO", "The bound of the file was reached!")
                endtime = self.lasttime
                new_starttime = endtime - self.interval - self.olap
            else:
                endtime = None
                new_starttime = self.starttime + self.interval - self.olap
            
            self.plot(new_starttime, endtime=endtime)

        if key == QtCore.Qt.Key_Left:
            if self.starttime - self.interval < self.firsttime:
                notify("SEISVO", "The bound of the file was reached!")
                new_starttime = self.firsttime
            else:
                new_starttime = self.starttime - self.interval + self.olap
            
            self.plot(new_starttime)

        day_interval = dt.timedelta(hours=24)
        
        if key == QtCore.Qt.Key_Up:
            if self.starttime + day_interval > self.lasttime:
                notify("SEISVO", "The bound of the file was reached!")
                endtime = self.lasttime
                new_starttime = endtime - day_interval - self.olap
            else:
                endtime = None
                new_starttime = self.starttime + day_interval - self.olap

            self.plot(new_starttime, endtime=endtime)

        if key == QtCore.Qt.Key_Down:
            if self.starttime - day_interval < self.firsttime:
                notify("SEISVO", "The bound of the file was reached!")
                new_starttime = self.firsttime
            else:
                new_starttime = self.starttime - day_interval + self.olap
            
            self.plot(new_starttime)


    def plot(self, starttime, endtime=None):
        self.starttime = starttime
        self.canvas.plot(endtime=endtime)

        if self.canvas.DayPlotWidget:
            self.canvas.DayPlotWidget.update(self)


    def closeEvent(self, event):
        close = QtWidgets.QMessageBox.question(self, "QUIT", "Are you sure want to exit?", QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        
        if close == QtWidgets.QMessageBox.Yes:
            try:
                self.canvas.DayPlotWidget.close()
            except:
                pass

            try:
                self.canvas.eGUI.close()
            except:
                pass
        
            event.accept()
        else:
            event.ignore()


class StationMainCanvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.DayPlotWidget = None
        
        # init event drawer object
        if self.parent.is_db:
            self.eDraw = EventDrawer(self)
            self.eGUI  = None

        self.fig = Figure(figsize=(12,9))
        FigureCanvas.__init__(self, self.fig)

        self.callbacks.connect('key_press_event', self.on_key)
        self.callbacks.connect('button_press_event', self.on_click)
        self.plot()


    def plot(self, endtime=None):

        self.ticks = dict(right=None, left=None)

        if not endtime:
            self.endtime = self.parent.starttime + self.parent.interval
        else:
            self.endtime = endtime

        with pyqtgraph.BusyCursor():

            # put events of database into a dictionary
            if self.parent.is_db:
                self.eDraw.update()
            
            # plot data
            self.axes = _plot_stream(self.fig, self.parent.station, self.parent.starttime, self.endtime, self.parent.spec_v, self.parent.spec_wlen, **self.parent.sta_kwargs)

            # update widget
            if self.DayPlotWidget:
                self.DayPlotWidget.update(self.parent)

        self.nav = Navigate([self.axes[0]], self, im_axis=[self.axes[1]], color='red', active_cursor=True, linewidth=0.5, alpha=0.5)
        
        # plot events
        if self.parent.is_db:
            self.eDraw.draw()


    def get_time_ticks(self):
        tL = self.nav.ticks['left'][0]
        if tL:
            tLeft = mdates.num2date(tL).replace(tzinfo=None)
        else:
            tLeft = None

        tR = self.nav.ticks['right'][0]
        if tR:
            tRight = mdates.num2date(tR).replace(tzinfo=None)
        else:
            tRight = None

        return (tLeft, tRight)


    def print_info(self):
        print("\n -------------------")
        print("      INFO      ")
        print(" -------------------")

        if self.parent.is_db:
            print("\n :: SDE :: ")
            print(self.parent.db)

        tLeft, tRight = self.get_time_ticks()

        print("\n :: TICKS :: ")
        if tLeft:
            print(f"  L >> {tLeft} ")

        if tRight:
            print(f"  R >> {tRight} ")
        
        if tLeft and tRight:
            dist   = abs((tLeft - tRight).total_seconds())
            print(f" distance |R-L|  >>  {dist}  [sec]\n")

            print("\n :: STATIONS :: ")
            start = min((tLeft, tRight))
            end   = max((tLeft, tRight))
            print(self.parent.all_stations(start, end))

        else:
            print("\n :: STATIONS :: ")
            print(self.parent.all_stations(self.parent.starttime, self.endtime))

        print("\n ------------------- \n")


    def on_click(self, event):
        if event.inaxes == self.axes[1]:
            print(f"\n Freq :: {event.ydata:.2f} Hz")


    def goto_time(self, parent):
        minmax = (self.parent.firsttime, self.parent.lasttime)
        new_time, ok = DateDialog.getDateTime(self.parent.starttime, minmax, parent)
        if ok:
            self.parent.plot(new_time)


    def _delete_event(self, eid):
        self.parent.db.remove_event(eid)
        self.eDraw.remove_event(eid)
        
        if self.DayPlotWidget:
            self.DayPlotWidget.update(self.parent)
        
        if self.eGUI:
            if len(self.parent.db) > 1:
                self.eGUI.update()
            else:
                self.eGUI.close()
                self.eGUI = None


    def _relabel_event(self, eid, label):
        self.parent.db.relabel_event(eid, label)
        self.eDraw.relabel_event(eid)
        
        if self.DayPlotWidget:
            self.DayPlotWidget.update(self.parent)
        
        if self.eGUI:
            self.eGUI.update()


    def on_key(self, event):
        if event.key == "h":
            print("   key   --   info   ")
            print(" --------  ----------")
            print(" NAVIGATE ")
            print("   ->    --  advance interval")
            print("   <-    --  retroces interval")
            print("  ^(up)  --  advance day")
            print(" v(down) --  retroces day")
            print("   tab   --  select day")
            print("backspace--  select station")
            print("\n TOOLS ")
            print("    i    --  print info")
            print("    d    --  interact with DayPlot GUI")
            print("    z    --  plot network Z component")
            print("    p    --  plot polar GUI")
            print("    f    --  change bandpass filter")
            print("    v    --  change specgram params")
            print("\n EVENTS ")
            print("    s    --  save event into database. Ask label")
            print("    1    --  save VT into database")
            print("    2    --  save LP into database")
            print("    3    --  save RG into database")
            print("    4    --  save TR into database")
            print("    5    --  save VLP into database")
            print("    r    --  relabel saved event")
            print("    +    --  print phase info")
            print("    e    --  init EventWidget GUI")
            print("   l/L   --  locate event / clean location")
            print("   del   --  remove saved event")
            return


        if event.key == 'i':
            self.print_info()
            return


        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)
            return


        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)
            return


        if event.key == 'up':
            self.parent.on_key(QtCore.Qt.Key_Up)
            return


        if event.key == 'down':
            self.parent.on_key(QtCore.Qt.Key_Down)
            return


        if event.key == "d":
            if not self.DayPlotWidget:
                self.DayPlotWidget = DayPlotWidget(self.parent)
                self.DayPlotWidget.show()
            else:
                self.DayPlotWidget.update(self.parent)
            return


        if event.key == "z":
            tLeft, tRight = self.get_time_ticks()
            fq_band       = self.parent.sta_kwargs.get("fq_band", ())
            if tLeft and tRight:
                start = min((tLeft, tRight))
                end   = max((tLeft, tRight))
                fig = _net_ztrace(self.parent.network, start, end, fq_band=fq_band)
                figdialog = SimplePlotDialog(fig)
                figdialog.exec_()
            return


        if event.key == "p":
            tLeft, tRight = self.get_time_ticks()
            fq_band       = self.parent.sta_kwargs.get("fq_band", ())
            if tLeft and tRight:
                self.print_info()
                start  = min((tLeft, tRight))
                end    = max((tLeft, tRight))
                stream = self.parent.station.get_stream(start, end, fq_band=fq_band)
                max_win = int(0.8*stream[0].stats.npts)
                win, ok = QtWidgets.QInputDialog.getInt(self, "Set window length", "Interval:", 50, 5, max_win)
                if ok:
                    try:
                        fig    = _zne_polar(stream, win, rl_th=0.7)
                        figdialog = SimplePlotDialog(fig)
                        figdialog.exec_()
                    except:
                        print(" warn :: error by computing polarization")
            return


        if event.key == "f":
            fq_band = self.parent.sta_kwargs.get("fq_band", ())

            if fq_band:
                curr_freq = f"{fq_band[0]};{fq_band[1]}"
            else:
                curr_freq = ";"
            
            new_freq, ok = QtWidgets.QInputDialog.getText(self, "Change frequency","Input new frequency range: ", text=curr_freq)

            if ok:
                try:
                    fq0, fq1 = new_freq.split(";")
                    fq0 = float(fq0)
                    fq1 = float(fq1)
                    if fq0 < fq1:
                        self.parent.sta_kwargs["fq_band"] = [fq0, fq1]
                        self.plot()
                except:
                    print("\n warn :: Frequency range could not be changed!")
            return


        if event.key == "v":
            vmin, vmax = self.parent.spec_v
            wlen       = self.parent.spec_wlen

            curr_vspec = f"{vmin};{vmax}/{wlen}"
            new_vspec, ok = QtWidgets.QInputDialog.getText(self, "Change specgram limits","Input new range (vmin;vmax/vlen): ", text=curr_vspec)

            if ok:
                try:
                    v, w =  new_vspec.split("/")
                    v0, v1 = v.split(";")
                    v0 = float(v0)
                    v1 = float(v1)
                    wl = float(w)
                    if v0 < v1:
                        self.parent.spec_v = [v0, v1]
                        self.parent.spec_wlen = wl
                        self.plot()
                except:
                    print("\n warn :: Specgram parameters could not be changed!")
            return


        if event.key == "tab":
            self.goto_time(self)
            return

        # EVENT DATABASE #
        if event.key in ('s', '1', '2', '3', '4', '5', "delete", "r", "l", "L", "e", "+") and not self.parent.is_db:
            print(" warn :: No SDE database found!")

        else:
            if event.key in ('s', '1', '2', '3', '4', '5'):
                tLeft, tRight = self.get_time_ticks()
                if tLeft and tRight:
                    start = min((tLeft, tRight))
                    end   = max((tLeft, tRight))
                    duration = (end-start).total_seconds()
                    
                    if event.key == "s":
                        # take label
                        label, ok = QtWidgets.QInputDialog.getText(self, "Input label","Label: ", text="")
                        if ok and label != "":
                            label = label.upper()
                        else:
                            return
                    else:
                        label = EVENTS_KEYS[event.key]

                    # search all station id availables
                    station_id_list = self.parent.all_stations(start, end)

                    # save into database
                    eid = self.parent.db.add_event(station_id_list, label, start, duration)

                    # reset
                    self.nav.reset_ticks()
                    self.eDraw.add_event(eid)

                    if self.DayPlotWidget:
                        self.DayPlotWidget.update(self.parent)
                    
                    if self.eGUI:
                        self.eGUI.update()
                return


            if event.inaxes == self.axes[0]:
                xtime = mdates.num2date(event.xdata).replace(tzinfo=None)
                ans = self.eDraw.get_eid(xtime)
                if ans:
                    if event.key == "delete":
                        self._delete_event(ans)
                        return

                    if event.key == "r":
                        event = self.eDraw.eventDict[ans]["event"]
                        label, ok = QtWidgets.QInputDialog.getText(self, "Relabel event","Label: ", text=event.label)
                        if ok and label != event.label:
                            self._relabel_event(ans, label.upper())
                        return
                            
                    if event.key == "e":
                        ans        = self.eDraw.get_eid(xtime)
                        station_id = self.parent.station.stats.id
                        if self.eGUI:
                            self.eGUI.update(station_id, ans)
                        else:
                            self.eGUI = EventWidget(self.parent.db, ans, station_id, parent=self)
                            self.eGUI.show()
                        return

                    if event.key == "+":
                        event = self.eDraw.eventDict[ans]["event"]
                        event.info()
                        return

                    if event.key in ("l", "L") and self.parent.hyp:
                        # locate event in nro phases is 4
                        evnt = self.eDraw.eventDict[ans]["event"]

                        if event.key == "L":
                            return
                            # clean location
                            # remove string_2
                            # remove folder and files in ./loc

                        else:
                            ok = evnt.locate(self.parent.hyp, './loc')
                            if ok:
                                evnt.loc.show() # show KML file
                                self.eDraw.update()



class DayPlotWidget(QtWidgets.QWidget):
    def __init__(self, parent):

        QtWidgets.QWidget.__init__(self)
        self.layout = QtWidgets.QVBoxLayout()

        self.parent  = parent
        self.fq_band = parent.sta_kwargs.get("fq_band", ())
        self.starttime = parent.starttime
        self.interval = dt.timedelta(hours=12)

        # add canvas
        self.canvas = DayPlotCanvas(self)
        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)
        # self.setCursor(QtCore.Qt.CrossCursor)
        self.canvas.setFocus()


    def on_key(self, key):
        if key == QtCore.Qt.Key_Right:
            if self.starttime + self.interval > self.parent.lasttime:
                notify("SEISVO", "The bound of the file was reached!")
                endtime = self.parent.lasttime
                new_starttime = endtime - self.interval - self.parent.olap
            else:
                endtime = None
                new_starttime = self.starttime + self.interval - self.parent.olap
            
            self.plot(new_starttime, endtime=endtime)

        if key == QtCore.Qt.Key_Left:
            if self.starttime - self.interval < self.parent.firsttime:
                notify("SEISVO", "The bound of the file was reached!")
                new_starttime = self.parent.firsttime
            else:
                new_starttime = self.starttime - self.interval + self.parent.olap
            
            self.plot(new_starttime)


    def plot(self, starttime, endtime=None):
        self.starttime = starttime
        self.canvas.plot(endtime=endtime)


    def update(self, parent):
        self.parent = parent
        self.plot(parent.starttime)


    def closeEvent(self, event):
        self.parent.canvas.DayPlotWidget = None
        event.accept()


class DayPlotCanvas(FigureCanvas):
    def __init__(self, parent):
        self.parent = parent
        self.fig = Figure(figsize=(12,9))
        FigureCanvas.__init__(self, self.fig)
        self.callbacks.connect('button_press_event', self.on_click)
        self.callbacks.connect('key_press_event', self.on_key)
        self.plot()


    def plot(self, endtime=None):
        if not endtime:
            self.endtime = self.parent.starttime + self.parent.interval
        else:
            self.endtime = endtime

        # get events
        events = self.parent.parent.canvas.eDraw.export((self.parent.starttime, self.endtime))

        with pyqtgraph.BusyCursor():
            self.axis = _day_plot(self.fig, self.parent.parent.station, self.parent.starttime, self.endtime, events, fq_band=self.parent.fq_band, interval=60, right_vertical_labels=False, vertical_scaling_range=15e3, one_tick_per_line=True, color=['k', 'r'], show_y_UTC_label=False)

        self.draw()


    def on_click(self, event):
        if event.inaxes:
            if event.inaxes == self.axis:
                tlist = self.axis.yaxis.get_ticklabels()
                pos   = np.array([t.get_position()[1] for t in tlist])
                time  = tlist[np.argmin(np.abs(pos-event.ydata))].get_text()
                starttime = dt.datetime.strptime(time, "%d/%m/%Y \n %H:%M:%S")
                self.parent.parent.plot(starttime)


    def on_key(self, event):
        if event.key == 'right':
            self.parent.on_key(QtCore.Qt.Key_Right)

        if event.key == 'left':
            self.parent.on_key(QtCore.Qt.Key_Left)

        if event.key == "tab":
            self.parent.parent.canvas.goto_time(self)

