#!/usr/bin/env python
# coding=utf-8

import os, sys
import numpy as np
import datetime as dt
import pyqtgraph

import matplotlib.dates as mdates
from matplotlib.figure import Figure
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas

from seisvo.plotting import plot_gram
from seisvo.plotting.gui import notify
from seisvo.plotting.gui.frames import lte_base, lte_pdf
from seisvo.lte.plotting import ltaoutsta_plot, get_axes_dict

# from seisvo.plotting.base.lte import plotLTE, defaultLTE_color, default_labels

default_LTE_color = get_colors('zesty')[1]
default_color_dictlabels = {
    'NB':'k',
    'BB':'k',
    'HT':'k',
    'MT':'k'
}

CLICK_COLOR = ('r', 'darkgreen')


def plot_gui(lte, list_chan, list_attr, starttime, endtime, interval, lde, lde_parent, **kwargs):
    app = QtWidgets.QApplication(sys.argv)
    lte_window = LTEWindow(chan, lte, starttime, endtime, interval, list_attr, lde, lde_parent, **kwargs)
    lte_window.show()
    sys.exit(app.exec_())


class LTEWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        QtWidgets.QMainWindow.__init__(self)
        self.ui = lte_base.Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle('LTE: '+ lte.stats.id)

        self.lte        = args[0]
        self.list_chan  = args[1]
        self.list_attr  = args[2]
        self.starttime  = args[3]
        self.endtime    = args[4]
        self.interval   = args[5]
        self.lde        = args[6]
        self.lde_parent = args[7]
        self.ckwargs    = kwargs

        self.lteout = self.lte.get(attr=self.list_attr, chan=self.list_chan, starttime=self.starttime, endtime=self.endtime)
        self.set_canvas()


    def set_canvas(self):
        self.canvas = LTECanvas(self)
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
    def __init__(self, parent):
        self.parent    = parent
        self.pdf_frame = None
        self.events_   = {} 
        self.__figinit__()


    def __figinit__(self):
        self.fig = Figure(figsize=self.parent.ckwargs.get("figsize", (20,9)), dpi=100)
        self.axes = get_axes_dict(self.parent.lteout, self.parent.attr_list, self.parent.chan_list, self.fig)
        
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()
        
        self.cliks = dict(right=dict(tick=[], time=None), left=dict(tick=[], time=None))
        self.callbacks.connect('button_press_event', self.on_click)

        self.update_interval() # init interval
        self.plot()


    def set_statusbar(self, interval):
        text = f'\t {self.parent.starttime.strftime('%d %B %Y')}  --- {self.parent.endtime.strftime('%d %B %Y')}  | Interval: {interval:.0f}'
        self.parent.ui.statusbar.showMessage(text, 0)


    def show_info(self):
        lte_info = self.parent.lte.__str__()
        self.parent.ui.textTickInfo.setText(lte_info)


    def update_interval(self):
        self.delta = dt.timedelta(days=self.parent.interval)
        self.set_statusbar(self.interval)


    def plot(self):
        self.fig.clf()

        if self.parent.starttime < self.parent.lte.stats.starttime:
            self.parent.starttime = self.parent.lte.stats.starttime
            notify('LTE', 'starttime set to LTE.starttime', status='info')

        self.parent.endtime = self.parent.starttime + self.delta

        if self.parent.endtime > self.parent.lte.stats.endtime:
            self.parent.endtime = self.parent.lte.stats.endtime
            notify('LTE', 'endtime set to LTE.endtime', status='info')
        
        interval = (self.parent.endtime - self.parent.starttime).days

        self.parent.ckwargs['settitle'] = False
        self.col_dict = self.parent.ckwargs.get('col_dict', default_color_dictlabels)

        self.stats = ltaoutsta_plot(
                        self.parent.lteout, 
                        self.parent.chan_list, 
                        self.parent.attr_list, 
                        fig=self.fig, 
                        axes=self.axes, 
                        plot=False,
                        return_stats=True,
                        datetime=self.parent.ckwargs.get("datetime"), 
                        **kwargs)

        # show episodes in LDE
        ids = self.parent.lde.get_episodes_id(time_interval=(self.parent.starttime, self.parent.endtime))

        if ids:
            self.show_events(ids)
        
        with pyqtgraph.BusyCursor():
            self.reset_ticks()


#----------------------------------------------
#-----------  EVENTS --------------
#----------------------------------------------

    def show_events(self, ids, reshow=False):        
        if not reshow:
            # self.plte.show_events(self.parent.lde, col_dict=self.col_dict)
            for eid in ids:
                self.events_[eid] = {'ticks':[], 'event':self.parent.lde[eid]}
                for key in self.axes:
                    if len(key) == 1: # scalar parameters
                        ax  = self.axes[key][0]
                        col = self.col_dict.get(episode.label, defaultLTE_color)
                        v1  = ax.axvline(episode.starttime,  alpha=0.5, color=col, ls='dashed')
                        v2  = ax.axvline(episode.endtime,  alpha=0.5, color=col, ls='dashed')
                        v3  = ax.avspan(episode.starttime,  episode.endtime, alpha=0.15, color=col)
                        self.events_[eid]['ticks'] += [v1,v2,v3]
                        
                        if i == 0:
                            ymax = self.get_ymax(attr)
                            txt = f'ID:{episode.id}[{episode.label}]'
                            mid_bin = episode.starttime + dt.timedelta(hours=episode.duration/2)
                            t = ax.annotate(txt, (mid_bin, ymax), color=col, xycoords='data')
                            self.events_[eid]['text'] = t

            


        self.parent.ui.eventWidget.clear()
        if self.plte.events_:
            for _, deid in self.plte.events_.items():
                e = deid['event']
                info_evnt = '  ID: %s | Label: %s' % (e.id, e.label)
                self.parent.ui.eventWidget.addItem(info_evnt)


    def deselect_event(self):
        for _, deid in self.plte.events_.items():
            if deid['ticks'][0].get_color() == 'gold':
                deid['text'].set_color(self.col_dict.get(deid['event'].label, defaultLTE_color))
                for artist in deid['ticks']:
                    artist.set_color(self.col_dict.get(deid['event'].label, defaultLTE_color))
                return


    def select_event(self):
        self.deselect_event()

        item = self.parent.ui.eventWidget.selectedItems()[0]
        event_id = int(item.text().split(' ')[3])

        for item in self.plte.events_[event_id]['ticks']:
            item.set_color('gold')
    
        self.plte.events_[event_id]['text'].set_color('gold')
        self.draw()


    def remove_event(self, i=None):
        if not i:
            i, ok = QtWidgets.QInputDialog.getInt(self, "Remove Event", "Event ID:")
        
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

        if ok and self.parent.lde.is_event(i):
            self.parent.lde.remove_row(i)
            with pyqtgraph.BusyCursor():
                self.plte.clear_events(i)
                self.show_events(reshow=True)
                self.draw()


    def relabel_event(self, i):
        if self.parent.lde.is_event(i):
            evnt = self.parent.lde.get(i)
            current_label = evnt.label
            new_label, ok = QtWidgets.QInputDialog.getText(self, "Relabel Event", "Event Label: ", text=current_label)

            if ok:
                new_label = str(new_label).upper()
                self.parent.lde.relabel_row(i, new_label)
                
                with pyqtgraph.BusyCursor():
                    self.plte.clear_events(i)
                    self.show_events(reshow=True)
                    self.draw()


    def write_event(self):
        lticktime, rticktime, hr_dur, _ = self.get_infoticks()
        ticktimes = [lticktime, rticktime]

        if all(ticktimes):
            label, ok = QtWidgets.QInputDialog.getText(self, "Long Duration Event type", "Label:")
            if ok:
                net_code = self.parent.lte.stats.id.split('.')[0]
                sta_code = self.parent.lte.stats.id.split('.')[1]
                loc_code = self.parent.lte.stats.id.split('.')[2]
                self.parent.lde.__add_episode__(
                    net_code,
                    sta_code,
                    loc_code,
                    label.upper(),
                    min(ticktimes),
                    hr_dur,
                    self.parent.lte.lte_file
                )

                notify('New event', 'ID saved')

                with pyqtgraph.BusyCursor():
                    self.plte.clear_events()
                    self.show_events()
                    self.reset_ticks()

#----------------------------------------------
#-----------  TICKS   --------------
#----------------------------------------------

    def clear_click(self, button):
        # left click
        if button == 1:
            if self.cliks['left']['tick']:
                for ac_ver_bar in self.cliks['left']['tick']:
                    ac_ver_bar.remove()
                self.cliks['left']['tick'] = []

        # right click
        if button == 3:
            if self.cliks['right']['tick']:
                for ac_ver_bar in self.cliks['right']['tick']:
                    ac_ver_bar.remove()
                self.cliks['right']['tick'] = []


    def on_click(self, event):
        if event.inaxes in self.axes:
            t = mdates.num2date(float(event.xdata))
            t = t.replace(tzinfo=None)
            t = self.plte.xtime[np.argmin(np.abs(np.array(self.plte.xtime)-t))]
        else:
            return

        # double click clean clicks!
        if event.dblclick:
            self.clear_click(1)
            self.clear_click(3)
        
        else:
            if event.button in (1, 3):
                self.clear_click(event.button)
                for iax in self.axes:
                    if event.button == 1:
                        self.cliks['left']['time'] = t
                        tax = iax.axvline(t, color=CLICK_COLOR[0])
                        self.cliks['left']['tick'].append(tax)

                    if event.button == 3:
                        self.cliks['right']['time'] = t
                        tax = iax.axvline(t, color=CLICK_COLOR[1])
                        self.cliks['right']['tick'].append(tax)
            else:
                return

        self.draw()
        self.update_buttons()
        self.show_infoticks()


    def reset_ticks(self):
        self.clear_click(1)
        self.clear_click(3)
        self.update_buttons()
        self.show_infoticks()
        self.draw()


    def get_infoticks(self):
        lticktime = self.cliks['left'].get('time', None)
        rticktime = self.cliks['right'].get('time', None)
        
        ticktimes = [lticktime, rticktime]
        if any(ticktimes):
            if all(ticktimes):
                hr_dur = (max(ticktimes)-min(ticktimes)).total_seconds()/3600.0
                dict_ans = self.parent.lte.get_stats(self.parent.list_attr, chan=self.chan, starttime=min(ticktimes), endtime=max(ticktimes))
            else:
                ticktime = list(filter(None, ticktimes))[0]
                hr_dur = np.nan
                dict_ans = self.parent.lte.get_stats(self.parent.list_attr, chan=self.chan, starttime=ticktime, endtime=None)
        else:
            hr_dur = np.nan
            dict_ans = None
        
        return lticktime, rticktime, hr_dur, dict_ans


    def show_infoticks(self, eid=None):
        if eid:
            episode = self.parent.lde.get(eid)
            lticktime, rticktime = episode.starttime, episode.endtime
            hr_dur = episode.duration
            dict_ans = self.parent.lte.get_stats(self.parent.list_attr, chan=self.chan, starttime=episode.starttime, endtime=episode.endtime)
        
        else:
            lticktime, rticktime, hr_dur, dict_ans = self.get_infoticks()

            if lticktime:
                lticktime = lticktime.strftime("%Y-%m-%d %H:%M")

            if rticktime:
                rticktime = rticktime.strftime("%Y-%m-%d %H:%M")

        text_tick = " Left click : %s\n" % lticktime
        text_tick += " Right click : %s\n" % rticktime
        text_tick += " Time delta [hr] : %2.1f hr\n" % hr_dur

        for attr, values in dict_ans.items():
            text_tick += " %s ::\n" % attr
            for txt, val in zip(['min', 'max', 'mean', 'mode'], values):
                text_tick += "   %s : %2.2f\n" % (txt, val)

        self.parent.ui.textTickInfo.setText(text_tick)

#----------------------------------------------
#-----------  GUI & NAVIGATE   --------------
#----------------------------------------------

    def update_buttons(self):
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


    def move_forward(self):
        self.starttime += self.delta
        self.plot()


    def move_backwards(self):
        self.starttime -= self.delta
        self.plot()


    def set_interval(self):
        max_int = (self.parent.lte.stats.endtime - self.parent.lte.stats.starttime).days
        new_interval, ok = QtWidgets.QInputDialog.getInt(self, "Set interval", "Interval [days]:", self.interval, 5, int(max_int), 5)
        if ok and int(new_interval) != self.interval:
            self.interval = int(new_interval)
            self.update_interval()
            self.plot()


    def set_starttime(self):
        current_day = self.starttime.strftime("%Y%m%d")
        time, ok = QtWidgets.QInputDialog.getText(self, "Set Starttime", "Time (YYYYMMDD):", text=current_day)
        if ok and str(time) != current_day:
            try:
                time = dt.datetime.strptime(str(time), "%Y%m%d")
            except:
                time = dt.datetime.now()

            if self.parent.lte.stats.starttime <= time < self.parent.lte.stats.endtime:
                self.starttime = time
                self.plot()
            
            else:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("Time out of the bounds.")
                msg.setWindowTitle("Error")
                msg.exec_()


    def show_pdf(self, eid=None):
        if eid:
            episode = self.parent.lde.get(eid)
            starttime = episode.starttime
            endtime = episode.endtime
        
        else:
            ticktimes = [self.cliks['right']['time'], self.cliks['left']['time']]
            starttime = min(ticktimes)
            endtime = max(ticktimes)
            
        with pyqtgraph.BusyCursor():
            if self.pdf_frame:
                self.pdf_frame.update(starttime, endtime)
                self.pdf_frame.show()
            else:
                self.pdf_frame = PDFWindow(self.parent, starttime, endtime)
                self.pdf_frame.show()


    def save_fig(self):
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save F:xile", "./", "Images (*.png *.svg *.pdf)")
        if fileName:
            self.fig.savefig(fileName)


class PDFWindow(QtWidgets.QMainWindow):
    def __init__(self, main, starttime, endtime):
        QtWidgets.QMainWindow.__init__(self)
        self.ui = lte_pdf.Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.main = main
        self.starttime = starttime
        self.endtime = endtime

        attrs = self.main.list_attr
        self.scalar_attrs = self.main.lte.is_attr(attrs, only_scalars=True)
        self.vector_attrs = self.main.lte.is_attr(attrs, only_vectors=True)

        if self.scalar_attrs:
            self.scalar_canvas = PDFCanvas(self, 1, init_attr=self.scalar_attrs[0])
            self.ui.horizontalLayout.addWidget(self.scalar_canvas)
        else:
            self.scalar_canvas = None

        if self.vector_attrs:
            self.vector_canvas = PDFCanvas(self, 2, init_attr=self.vector_attrs[0])
            self.ui.horizontalLayout_2.addWidget(self.vector_canvas)
        else:
            self.vector_canvas = None

        self.set_window_title()
        self.add_buttons()
    

    def set_window_title(self):
        self.setWindowTitle('PDF %s -- %s' % (
            self.starttime.strftime('%d %b %Y'),
            self.endtime.strftime('%d %b %Y'))
        )


    def add_buttons(self):
        self.button_ = {}
        
        # scalar attributes
        for sat in self.scalar_attrs:
            pushButton = QtWidgets.QPushButton(self.ui.frame)
            self.ui.horizontalLayout.addWidget(pushButton)
            pushButton.setText(sat)
            pushButton.clicked.connect(lambda: self.scalar_canvas.plot(sat))
            self.button_[sat] = pushButton

        # vector attributes
        for vat in self.vector_attrs:
            pushButton = QtWidgets.QPushButton(self.ui.frame_2)
            self.ui.horizontalLayout.addWidget(pushButton)
            pushButton.setText(vat)
            pushButton.clicked.connect(lambda: self.vector_canvas.plot(vat))
            self.button_[vat] = pushButton


    def update(self, starttime, endtime):
        self.starttime = starttime
        self.endtime = endtime
        self.set_window_title()

        if self.scalar_canvas:
            self.scalar_canvas.plot(self.scalar_canvas.last_attr)
        
        if self.vector_canvas:
            self.vector_canvas.plot(self.vector_canvas.last_attr)


class PDFCanvas(FigureCanvas):
    def __init__(self, parent, ncol, init_attr=None):
        self.parent = parent
        self.lte = parent.main.parent.lte

        self.fig = Figure(figsize=(9,3.5), dpi=100)
    
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,
            QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.setParent(self.parent)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setFocus()

        self.set_frame(ncol)
        self.plot(init_attr)
    

    def set_frame(self, ncol):
        grid = {'hspace':0.3, 'left':0.08, 'right':0.92, 'wspace':0.1, 'top':0.95, 'bottom':0.05}
        
        if ncol == 2:
            grid['width_ratios'] = [1, 0.01]
            grid['wspace'] = 0.01
            self.vect_attr = True
        else:
            self.vect_attr = False
    
        self.axes = self.fig.subplots(1, ncol, gridspec_kw=grid)
    

    def plot(self, attr):
        self.last_attr = attr

        if self.vect_attr:
            ax = self.axes[0, 0]
            ax.cla()
            cax = self.axes[0, 1]
            cax.cla()

            x_space, y_space, pdf = self.lte.pdf(attr, chan=self.parent.main.chan, starttime=self.parent.starttime, endtime=self.parent.endtime)
            
            pltkwargs = {
                'y_label': r'$f$ [Hz]', 
                'x_label': default_labels.get(attr, attr), 
                'v_max':np.percentile(pdf, 95), 
                'v_min':np.percentile(pdf, 5), 
                'axis_bar':cax, 
                'bar_label':'PDF'
                }

            plot_gram(y_space, pdf, x_space, ax, **pltkwargs)

        else:
            self.axes.cla()
            x_space, pdf = self.lte.pdf(attr, chan=self.parent.main.chan, starttime=self.parent.starttime, endtime=self.parent.endtime)
            self.axes.plot(x_space, pdf, color='k')
            self.axes.set_xlabel(default_labels.get(attr, attr))
            self.axes.set_ylabel('PDF')




        







    
