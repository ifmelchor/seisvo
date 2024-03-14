#!/usr/bin/env python3
# coding=utf-8

from pynotifier import Notification, NotificationClient
from pynotifier.backends import platform

from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas, NavigationToolbar2QT
from matplotlib.widgets import AxesWidget
import matplotlib.dates as mdates

P_PICKER_KEYS = ['P', 'p', '1', '2', '3', '4']
S_PICKER_KEYS = ['S', 's', '6', '7', '8', '9']


def notify(title, message, duration=2):
    c = NotificationClient()
    c.register_backend(platform.Backend())
    n = Notification(title=title, message=message, duration=duration)
    c.notify_all(n)


class SimplePlotWidget(QtWidgets.QWidget):
    def __init__(self, fig, parent=None):
        super(SimplePlotWidget, self).__init__(parent)
        self.canvas = FigureCanvas(fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)


class SimplePlotDialog(QtWidgets.QDialog):
  def __init__(self, fig, parent=None):
    super(SimplePlotDialog, self).__init__(parent)
    self.setWindowTitle("Plot Figure")
    self.plot_widget = SimplePlotWidget(fig, parent)
    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(self.plot_widget)
    self.setLayout(layout)
    self.plot_widget.canvas.draw()


class DateDialog(QtWidgets.QDialog):
    def __init__(self, current_time, minmax_times, parent=None):
        super(DateDialog, self).__init__(parent)

        layout = QtWidgets.QVBoxLayout(self)

        # nice widget for editing the date
        self.datetime = QtWidgets.QDateTimeEdit(self)
        self.datetime.setCalendarPopup(True)
        self.datetime.setDateTime(QtCore.QDateTime(current_time))
        self.datetime.setMinimumDateTime(QtCore.QDateTime(minmax_times[0]))
        self.datetime.setMaximumDateTime(QtCore.QDateTime(minmax_times[1]))
        layout.addWidget(self.datetime)

        # OK and Cancel buttons
        self.buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal, self)
        layout.addWidget(self.buttons)

        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

    # get current date and time from the dialog
    def dateTime(self):
        return self.datetime.dateTime()

    # static method to create the dialog and return (date, time, accepted)
    @staticmethod
    def getDateTime(current_time, minmax_times, parent=None):
        dialog = DateDialog(current_time, minmax_times, parent)
        result = dialog.exec_()
        date = dialog.dateTime()
        return (date.toPyDateTime(), result==QtWidgets.QDialog.Accepted)


class Cursor(AxesWidget):
    def __init__(self, ax, horizOn=True, vertOn=True, useblit=False, **lineprops):
        
        AxesWidget.__init__(self, ax)

        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('draw_event', self.clear)

        self.visible = True
        self.horizOn = horizOn
        self.vertOn = vertOn
        self.useblit = useblit and self.canvas.supports_blit

        if self.useblit:
            lineprops['animated'] = True

        self.lineh = ax.axhline(ax.get_ybound()[0], visible=False, **lineprops)
        self.linev = ax.axvline(ax.get_xbound()[0], visible=False, **lineprops)

        self.background = None
        self.needclear = False


    def clear(self, event):
        """Internal event handler to clear the cursor."""
        if self.ignore(event):
            return
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        
        self.linev.set_visible(False)
        self.lineh.set_visible(False)


    def onmove(self, event):
        """Internal event handler to draw the cursor when the mouse moves."""
        
        if self.ignore(event):
            return
        
        if not self.canvas.widgetlock.available(self):
            return
        
        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)

            if self.needclear:
                self.canvas.draw()
                self.needclear = False
            return
        
        self.needclear = True
        
        if not self.visible:
            return
        
        self.linev.set_xdata((event.xdata, event.xdata))
        self.linev.set_visible(self.visible and self.vertOn)

        self.lineh.set_ydata((event.ydata, event.ydata))
        self.lineh.set_visible(self.visible and self.horizOn)

        self._update()


    def _update(self):
        if self.useblit:    
            if self.background is not None:
                self.canvas.restore_region(self.background)
            
            self.ax.draw_artist(self.linev)
            self.ax.draw_artist(self.lineh)
            self.canvas.blit(self.ax.bbox)

        else:
            self.canvas.draw_idle()
        
        return False


class Navigate(object):
    def __init__(self, axes, canvas, im_axis=None, base_scale=2, active_cursor=True, **lineprops):
        # lineprops --> color='red', linewidth=0.5, alpha=0.5
        """
        axes is a list of Axes
        """

        self.ax  = axes
        self.canvas = canvas
        self.base_scale = base_scale
        self.imshow_axes = im_axis

        if active_cursor:
            self.cursors = []
            for ax in axes:
                cursor = Cursor(ax, useblit=True, **lineprops)
                self.cursors.append(cursor)
            
        self.clicked = None
        self.xtick = None

        self.cur_xlim = None
        self.new_xlim = None
        self.max_xlim = axes[0].get_xlim()

        self.ticks = dict(left=[None, []], right=[None, []]) # tick data and list of axvline
        self.canvas.mpl_connect('scroll_event', self.onZoom)
        self.canvas.mpl_connect('motion_notify_event', self.onMotion)
        self.canvas.mpl_connect('button_press_event', self.onClkPress)
        self.canvas.mpl_connect('button_release_event', self.onClkRelease)


    def reset(self):
        for ax in self.ax:
            ax.set_xlim(self.max_xlim)
        
        if self.imshow_axes:
            for imax in self.imshow_axes:
                imax.set_xlim(self.max_xlim)
        
        self.canvas.draw()


    def reset_ticks(self):
        if self.ticks['left'][0]:
            [line.remove() for line in self.ticks['left'][1]]
            self.ticks['left'][0] = None
        
        if self.ticks['right'][0]:
            [line.remove() for line in self.ticks['right'][1]]
            self.ticks['right'][0] = None


    def set_xlim(self, new_lim):
        if new_lim[0] < self.max_xlim[0]:
            new_lim[0] = self.max_xlim[0]

        if new_lim[1] > self.max_xlim[1]:
            new_lim[1] = self.max_xlim[1]

        for ax in self.ax:
            ax.set_xlim(new_lim)
        
        if self.imshow_axes:
            for imax in self.imshow_axes:
                imax.set_xlim(new_lim)


    def onZoom(self, event):
        if event.inaxes not in self.ax: 
            return
        
        cur_xlim = event.inaxes.get_xlim()
        xdata = event.xdata # get event x location

        if event.button == 'down':
            # deal with zoom in
            scale_factor = 1 / self.base_scale
        
        elif event.button == 'up':
            # deal with zoom out
            scale_factor = self.base_scale
        
        else:
            return
    
        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
        new_lim = [xdata - new_width * (1-relx), xdata + new_width * (relx)]

        self.set_xlim(new_lim)
        self.canvas.draw()


    def onClkPress(self, event):
        if event.inaxes not in self.ax:
            return

        if event.dblclick:
            self.reset()
            return
        
        # save click info
        self.clicked = event.button
        self.xtick   = event.xdata

        # get current xlim
        self.cur_xlim = event.inaxes.get_xlim()
        self.new_xlim = ()


    def onClkRelease(self, event):
        if event.inaxes not in self.ax: 
            return

        if event.xdata == self.xtick:
            # draw left/right ticks
            
            if self.clicked == 1:
                # remove old tick
                if self.ticks['left'][0]:
                    [line.remove() for line in self.ticks['left'][1]]
                
                self.ticks['left'][0] = self.xtick
                self.ticks['left'][1] = []

                tick_artist = []
                for ax in self.ax:
                    tick_artist.append(ax.axvline(event.xdata, color='r', alpha=1, ls='-', lw=1))

                self.ticks['left'][1] = tick_artist
            
            if self.clicked == 3:
                if self.ticks['right'][0]:
                    [line.remove() for line in self.ticks['right'][1]]
                
                self.ticks['right'][0] = self.xtick

                tick_artist = []
                for ax in self.ax:
                    tick_artist.append(ax.axvline(event.xdata, color='g', alpha=1, ls='-', lw=1))

                self.ticks['right'][1] = tick_artist
        
        else:
            # draw from onMotion
            self.set_xlim(self.new_xlim)
        
        self.canvas.draw()


    def onMotion(self, event):
        if event.inaxes not in self.ax: 
            return

        if self.clicked == 1:
            # compute new xlim
            dx = event.xdata - self.xtick
            self.new_xlim = self.cur_xlim - dx
            # draw in onClkRelease


class Picker(object):
    def __init__(self, axes, phase, sde, row_id, canvas, **kwargs):
        self.db  = sde
        self.row_id = row_id
        self.axes = axes
        self.phase  = phase
        self.canvas = canvas
        self.phase_colors = kwargs.get("phase_colors",{"P":"r", "S":"g", "F":"b"})
        self.canvas.mpl_connect('key_press_event', self.on_key)

    def clear(self, wave):
        assert wave in ["P", "S", "F"]

        if self.phase[wave]["artist"]:
            for artist in self.phase[wave]["artist"]:
                artist.remove()
            self.phase[wave]["artist_text"].remove()

        if wave == "P":
            self.phase[wave] = {
                "time":None,
                "weight":None,
                "onset":None,
                "artist":[],
                "artist_text":None
            }

        if wave == "S":
            self.phase[wave] = {
                "time":None,
                "weight":None,
                "artist":[],
                "artist_text":None
            }

        if wave == "F":
            self.phase[wave] = {
                "time":None,
                "artist":[],
                "artist_text":None
            }
        

    def draw(self):
        "draw the phases without artist"
        for wave, phase in self.phase.items():
            if phase["time"] and not phase["artist"]:
                for ax in self.axes:
                    phase['artist'].append(ax.axvline(phase["time"], color=self.phase_colors[wave], lw=1.1))
        
                txt = self.phase_text(wave, phase)
                phase['artist_text'] = self.axes[0].annotate(txt, xy=(phase["time"], 1), color=self.phase_colors[wave])

        self.canvas.draw()


    def update_text(self):
        "do not draw again, just change the artist text of each phase"

        self.phase["P"]["artist_text"].remove()
        txt = self.phase_text("P", self.phase["P"])
        self.phase["P"]['artist_text'] = self.axes[0].annotate(txt, xy=(self.phase["P"]["time"],1), color=self.phase_colors["P"])
        self.canvas.draw()


    def save(self):
        idict = {
            'time_P':self.phase["P"]['time'],
            'weight_P':self.phase["P"]['weight'],
            'onset_P':self.phase["P"]['onset'],
            'time_S':self.phase["S"]['time'],
            'weight_S':self.phase["S"]['weight'],
            'time_F':self.phase["F"]['time']
        }
        self.db._update_row(self.row_id, idict)


    def on_key(self, event):
        if event.inaxes:
            if event.inaxes in self.axes:
                t = mdates.num2date(float(event.xdata))
                t = t.replace(tzinfo=None)
                

                if event.key in P_PICKER_KEYS:
                    if event.key == "P": # remove P
                        self.clear("P")
                    else: # draw P
                        if event.key == 'p':
                            w = 0
                        else:
                            w = event.key
                        self.clear("P")
                        self.phase["P"]["time"] = t
                        self.phase["P"]["weight"] = w
                    self.draw()
                    self.save()
                

                if event.key in S_PICKER_KEYS:
                    if event.key == "S": # remove P
                        self.clear("S")
                    else: # draw P
                        if event.key == 's':
                            w = 0
                        else:
                            w = P_PICKER_KEYS[2:][S_PICKER_KEYS[2:].index(event.key)]
                        self.clear("S")
                        t = mdates.num2date(float(event.xdata))
                        t = t.replace(tzinfo=None)
                        self.phase["S"]["time"] = t
                        self.phase["S"]["weight"] = w
                    self.draw()
                    self.save()
                

                if event.key in ['f',  'F']:
                    if event.key == "F": # remove P
                        self.clear("F")
                    else: # draw P
                        self.clear("F")
                        self.phase["F"]["time"] = t
                    self.draw()
                    self.save()


                if event.key in ['c', 'C', 'd', 'D']:
                    if self.phase["P"]["time"]:
                        if event.key.lower() == self.phase["P"]['onset']:
                            self.phase["P"]['onset'] = None
                        else:
                            self.phase["P"]['onset'] = event.key
                        self.update_text()
                        self.save()


    @staticmethod
    def phase_text(wave, phase_dict):
        assert wave in ["P", "S", "F"]
        # print(phase_dict)
        a = phase_dict.get("onset", None)
        if not a:
            a = ""
        else:
            a = str(a)
        b = phase_dict.get("weight", None)
        if not b:
            b = ""
        else:
            b = str(b)
        txt = "".join([a,wave,b])
        return txt


