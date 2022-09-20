#!/usr/bin/env python3
# coding=utf-8

import pyqtgraph
import sys
from seisvo.plotting import get_colors
from seisvo.plotting.gui import Navigation
from matplotlib.figure import Figure, SubplotParams
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas


phase_colors={'p':get_colors('okabe')[3],
              's':get_colors('okabe')[0],
              'f':get_colors('okabe')[6]}

class MplCanvas(FigureCanvas):
    def __init__(self, sde, eid_list, eid, parent=None, **kwargs):
        figsize = (12,9)
        self.parent = parent
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
        self.callbacks.connect('key_press_event', self.on_key)

        self.sde = sde
        self.eid_list = eid_list
        self.index = self.eid_list.index(eid)

        self.remove_response = kwargs.get('remove_response', True)
        self.fq_band = kwargs.get('fq_band', (1,10))
        self.stations = kwargs.get('stations', None)
        self.off_time = kwargs.get('off_time', 0)

        self.load_event()
        self.__plot__()


    def load_event(self):
        self.eid = self.eid_list[self.index]
        self.event = self.sde[self.eid]


    def __plot__(self):
        self.fig.clf()
        self.ticks = dict(right=[], left=[])
        with pyqtgraph.BusyCursor():
            self.axes1, self.axes2 = self.event.plot(
                stations=self.stations,
                fq_band=self.fq_band, 
                remove_response=self.remove_response, 
                fig=self.fig, 
                return_axes=True,
                off_time=self.off_time
                )
            self.nav1 = [Navigation(ax, parent=self) for ax in self.axes1]
            self.nav2 = [Navigation(ax, parent=self) for ax in self.axes2]
            self.__initphases__()
            self.__drawphases__()


    def __initphases__(self):
        self.event = self.sde[self.eid] # reload in case that new stations were added
        self.phases = {}
        self._phases = {}
        
        for staid in self.event.stations:
            self.phases[staid] = {
                'time_P':None,
                'time_S':None,
                'time_F':None,
                'onset_P':None,
                'mov_P':None,
                'weight_P':None,
                'weight_S':None,
                'max_period':None
                }

            self._phases[staid] = {
                'p' : [['','P','',''], None, None],
                's' : [['S',''], None, None],
                'f' : ['F', None, None]
            }
            
            row = self.event.get_row(staid)
            
            if row.time_P:
                self.phases[staid]['time_P'] = row.time_P

            if row.onset_P:
                self.phases[staid]['onset_P'] = row.onset_P

            if row.first_motion_P:
                self.phases[staid]['mov_P'] = row.first_motion_P

            if isinstance(row.weight_P, int):
                self.phases[staid]['weight_P'] = row.weight_P

            if row.time_S:
                self.phases[staid]['time_S'] = row.time_S

            if isinstance(row.weight_S, int):
                self.phases[staid]['weight_S'] = row.weight_S

            if row.event_duration and row.time_P:
                self.phases[staid]['time_F'] = row.time_P + row.event_duration

            self.__updatephasetxt__(station=staid)


    def __drawphases__(self, station=None, phase=None):
        for ax in self.axes1:
            staid = ax.yaxis.get_label().get_text()
            
            if not station or station==staid:
                pass
            else:
                continue
            
            if (phase == 'p' or not phase) and self.phases[staid]['time_P']: # draw P
                if self._phases[staid]['p'][1]: # if exist, remove old P
                    self._phases[staid]['p'][1].remove()
                    self._phases[staid]['p'][2].remove()
                    self._phases[staid]['p'] = [['','P','',''], None, None] 

                self.__updatephasetxt__(staid, 'p')
                txt = ''.join(self._phases[staid]['p'][0])
                self._phases[staid]['p'][1] = ax.axvline(self.phases[staid]['time_P'], color=phase_colors['p'], lw=1.1)
                self._phases[staid]['p'][2] = ax.annotate(txt, xy=(self.phases[staid]['time_P'], 6.8), color=phase_colors['p'], fontsize=12)
            
            if (phase == 's' or not phase) and self.phases[staid]['time_S']: # draw S
                if self._phases[staid]['s'][1]: # if exist, remove old S
                    self._phases[staid]['s'][1].remove()
                    self._phases[staid]['s'][2].remove()
                    self._phases[staid]['s'] = [['S',''], None, None]

                self.__updatephasetxt__(staid, 's')
                txt = ''.join(self._phases[staid]['s'][0])
                self._phases[staid]['s'][1] = ax.axvline(self.phases[staid]['time_S'], color=phase_colors['s'], lw=1.1)
                self._phases[staid]['s'][2] = ax.annotate(txt, xy=(self.phases[staid]['time_S'], 6.8), color=phase_colors['s'], fontsize=12)
            
            if (phase == 'f' or not phase) and self.phases[staid]['time_F']: # draw F
                if self._phases[staid]['f'][1]: # if exist, remove old F
                    self._phases[staid]['f'][1].remove()
                    self._phases[staid]['f'][2].remove()
                    self._phases[staid]['f'] = ['F', None, None]

                txt = self._phases[staid]['f'][0]
                self._phases[staid]['f'][1] = ax.axvline(self.phases[staid]['time_F'], color=phase_colors['f'], lw=1.1)
                self._phases[staid]['f'][2] = ax.annotate(txt, xy=(self.phases[staid]['time_F'], 6.8), color=phase_colors['f'], fontsize=12)

        self.draw()


    def __clearphase__(self, staid, phase):
        if phase == 'p':
            self.phases[staid]['time_P'] = None
            self.phases[staid]['weight_P'] = None
            self.phases[staid]['onset_P'] = None
            self.phases[staid]['mov_P'] = None

            if self._phases[staid]['p'][1]:
                self._phases[staid]['p'][1].remove()
                self._phases[staid]['p'][2].remove()
                self._phases[staid]['p'] = [['','P','',''], None, None]
        
        if phase == 's':
            self.phases[staid]['time_S'] = None
            self.phases[staid]['weight_S'] = None

            if self._phases[staid]['s'][1]:
                self._phases[staid]['s'][1].remove()
                self._phases[staid]['s'][2].remove()
                self._phases[staid]['s'] = [['S',''], None, None]
        
        if phase == 'f':
            self.phases[staid]['time_F'] = None

            if self._phases[staid]['f'][1]:
                self._phases[staid]['f'][1].remove()
                self._phases[staid]['f'][2].remove()
                self._phases[staid]['f'] = ['F', None, None]


    def __updatephasetxt__(self, station=None, phase=None):
        for staid in self.event.stations:
            if not station or station==staid:
                pass
            else:
                continue
            
            if phase=='p' or not phase:
                if self.phases[staid]['onset_P']:
                    self._phases[staid]['p'][0][0] = self.phases[staid]['onset_P']
                else:
                    self._phases[staid]['p'][0][0] = ''
                
                if isinstance(self.phases[staid]['weight_P'], int):
                    self._phases[staid]['p'][0][2] = str(self.phases[staid]['weight_P'])
                
                if self.phases[staid]['mov_P']:
                    self._phases[staid]['p'][0][3] = self.phases[staid]['mov_P']
                else:
                    self._phases[staid]['p'][0][3] = ''
            
            if phase=='s' or not phase:
                if isinstance(self.phases[staid]['weight_S'], int):
                    self._phases[staid]['s'][0][1] = str(self.phases[staid]['weight_S'])


    def __todatabase__(self, station):
        info_dict = {
            'time_P': self.phases[station]['time_P'],
            'time_S': self.phases[station]['time_S'],
            'onset_P': self.phases[station]['onset_P'],
            'first_motion_P': self.phases[station]['mov_P'],
            'weight_P': self.phases[station]['weight_P'],
            'weight_S': self.phases[station]['weight_S']
            }
        if self.phases[station]['time_F'] and self.phases[station]['time_P']:
            info_dict['event_duration'] = self.phases[station]['time_F']-self.phases[station]['time_P']
        else:
            info_dict['event_duration'] = None

        row = self.event.get_row(station)
        self.sde.update_row(row.id, info_dict)


    def on_click(self, event):
        if event.inaxes:
            for ax in self.axes1:
                if ax == event.inaxes:
                    station = ax.yaxis.get_label().get_text()
                    y = int(event.ydata)
                    if y == 0:
                        component = 'N'
                    elif y == 2:
                        component = 'E'
                    elif y == 4:
                        component = 'Z'
                    else:
                        component = '?'
                    
                    self.clean_ticks_on_nav(ax)
                    if event.button == 1:
                        self.ticks['left'] = [station, component, event.xdata]
                    
                    if event.button == 3:
                        self.ticks['right'] = [station, component, event.xdata]
                    
                    self.print_click_info()
                
            if event.inaxes in self.axes2:
                print(f' Click: {event.xdata:.2f} Hz')


    def clean_ticks_on_nav(self, ax):
        for n, axi in enumerate(self.axes1):
            if axi != ax:
                try:
                    self.nav1[n].reset_ticks()
                except:
                    pass
    

    def print_click_info(self):
        if self.ticks['right'] and self.ticks['left']: # if left and right click: print diff
            R_station = self.ticks['right'][0]
            L_station = self.ticks['left'][0]

            if R_station == L_station:
                diff = False
            else:
                diff =True
            
            if self.ticks['right'][2] < self.ticks['left'][2]:
                Ltick, Rtick = self.ticks['right'][2], self.ticks['left'][2]
            else:
                Rtick, Ltick = self.ticks['right'][2], self.ticks['left'][2]
            
            if diff:
                print(f' Click info :: {L_station} Left: {Ltick:.2f} sec | {R_station} Right: {Rtick:.2f} sec  || R-L: {Rtick-Ltick:.2f} sec')
            else:
                print(f' Click info {L_station} :: Left: {Ltick:.2f} sec | Right: {Rtick:.2f} sec  || R-L: {Rtick-Ltick:.2f} sec')
        
        else:
            if self.ticks['right']:
                print(f" Click info :: {self.ticks['right'][0]} Right: {self.ticks['right'][2]:.2f} sec")
            else:
                print(f" Click info :: {self.ticks['left'][0]} Left: {self.ticks['left'][2]:.2f} sec")


    def on_key(self, event):
        if event.key in ('left', 'backspace'):
            self.index -= 1
            if len(self.eid_list) < 0:
                self.index = len(self.eid_list) - 1
            
            self.load_event()
            self.__plot__()


        if event.key in ('right', 'enter'):
            self.index += 1
            if len(self.eid_list) == self.index:
                self.index = 0
            
            self.load_event()
            self.__plot__()


        if event.key == 'g': # go to EID
            eid , ok = QtWidgets.QInputDialog.getInt(None, 'Ir a EID', 'EID:', value=self.id_event, min=min(self.eid_list), max=max(self.eid_list))
            if ok and eid in self.eid_list:
                self.index = self.eid_list.index(self.eid)
                self.load_event()
                self.__plot__()

        
        if event.key == 'r': #relabel
            text , ok = QtWidgets.QInputDialog.getText(None, 'Re-etiqueta ID', 'Nueva Etiqueta:', text=self.event.label)
            if ok:
                self.sde.relabel_event(self.event.id, text.upper())
                self.load_event()
                self.__plot__()


        if event.key in ('p', '0', '1', '2', '3', '4', 'P'): # phase picking P
            for ax in self.axes1:
                if ax == event.inaxes:
                    station = ax.yaxis.get_label().get_text()
                    self.__clearphase__(station, 'p')
                    self.draw()
                    
                    if event.key != 'P':
                        self.phases[station]['time_P'] = event.xdata
                        
                        if event.key in ('0', '1', '2', '3', '4'):
                            self.phases[station]['weight_P'] = ['0', '1', '2', '3', '4'].index(event.key)
                        
                        else:
                            self.phases[station]['weight_P'] = ''
                        
                        self.__updatephasetxt__(station=station, phase='p')
                        self.__drawphases__(station=station, phase='p')
                    
                    self.__todatabase__(station)
        

        if event.key in ('i', 'e', 'c', 'd'):
            for ax in self.axes1:
                if ax == event.inaxes:
                    station = ax.yaxis.get_label().get_text()
                    
                    if event.key == self.phases[station]['onset_P'] or event.key == self.phases[station]['mov_P']:
                        
                        if event.key in ('c', 'd'): # remove first mov
                            self.phases[station]['mov_P'] = None
                        
                        else: # remove onset
                            self.phases[station]['onset_P'] = None
                    
                    else:
                        if event.key in ('c', 'd'): # update first mov
                            self.phases[station]['mov_P'] = event.key
                        
                        else: # update onset
                            self.phases[station]['onset_P'] = event.key
                    
                    self.__drawphases__(station=station, phase='p')
                    self.__todatabase__(station)


        if event.key in ('s', '=', '!', '"', '·', '$', 'S'): # phase picking S
            for ax in self.axes1:
                if ax == event.inaxes:
                    station = ax.yaxis.get_label().get_text()
                    self.__clearphase__(station, 's')
                    self.draw()
                    
                    if event.key != 'S':
                        self.phases[station]['time_S'] = event.xdata
                        
                        if event.key in ('=', '!', '"', '·', '$'):
                            self.phases[station]['weight_S'] = ['=', '!', '"', '·', '$'].index(event.key)
                        
                        else:
                            self.phases[station]['weight_S'] = ''
                        
                        self.__updatephasetxt__(station=station, phase='s')
                        self.__drawphases__(station=station, phase='s')
                    
                    self.__todatabase__(station)
        

        if event.key in ('f',  'F'): # phase picking F
            for ax in self.axes1:
                if ax == event.inaxes:
                    station = ax.yaxis.get_label().get_text()
                    self.__clearphase__(station, 'f')
                    self.draw()
                    
                    if event.key != 'F':
                        self.phases[station]['time_F'] = event.xdata
                        self.__drawphases__(station=station, phase='f')
                    
                    self.__todatabase__(station)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, sde, eid_list, eid=None, **kwargs):
        QtWidgets.QMainWindow.__init__(self)
        self.setWindowTitle('%s' %sde.sql_path)
        self.centralwidget = QtWidgets.QWidget(self)
        self.verticalLayout = QtWidgets.QHBoxLayout(self.centralwidget)

        if not eid:
            # show the first ID
            eid = eid_list[0]

        # self.msg = BoxWindow()
        self.canvas = MplCanvas(sde, eid_list, eid, **kwargs)
        self.verticalLayout.addWidget(self.canvas)
        self.setCentralWidget(self.centralwidget)

    # def closeEvent(self, event):
    #     self.msg.close()


def plot_SDE_event(sde, eid_list, eid=None, app=False, **kwargs):

    if not app:
        app = QtWidgets.QApplication(sys.argv)
        _exec = True
    else:
        _exec = False
    
    w = MainWindow(sde, eid_list, eid, **kwargs)
    w.show()
    
    if _exec:
        app.exec_()
    
    else:
        return w
