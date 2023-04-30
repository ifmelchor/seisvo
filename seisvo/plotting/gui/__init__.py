#!/usr/bin/env python
# coding=utf-8

import os
from pynotifier import Notification
from matplotlib.widgets import AxesWidget
import matplotlib.dates as mdates

path = os.path.dirname(os.path.realpath(__file__))
icons_path = os.path.join(path, 'icons')

P_PICKER_KEYS = ['P', 'p', '1', '2', '3', '4']
S_PICKER_KEYS = ['S', 's', '!', '"', '·', '$']


def notify(title, description, status='info', duration=2):
    """
    status: error, info, warn
    """
    n = Notification(
        title=title, 
        description=description, 
        icon_path='%s/%s.png' %(icons_path, status),
        duration=duration
        )
    n.send()


class _Cursor(AxesWidget):
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



class _Navigate(object):
    def __init__(self, axes, canvas, im_axis=None, base_scale=2, **lineprops):
        # lineprops --> color='red', linewidth=0.5, alpha=0.5
        """
        axes is a list of Axes
        """

        self.ax  = axes
        self.canvas = canvas
        self.base_scale = base_scale
        self.imshow_axes = im_axis
        self.cursors = []
        for ax in axes:
            cursor = _Cursor(ax, useblit=True, **lineprops)
            self.cursors.append(cursor)
        
        self.x0 = None
        self.x1 = None
        self.press    = None
        self.xpress   = None
        self.cur_xlim = None
        self.new_xlim = None
        self.max_xlim = ax[0].get_xlim()

        self.ticks = dict(left=[None, []], right=[None, []]) # tick data and list of axvline
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()
        self.canvas.mpl_connect('scroll_event', self.onZoom)
        self.canvas.mpl_connect('button_press_event', self.onClkPress)
        self.canvas.mpl_connect('button_release_event', self.onClkRelease)
        self.canvas.mpl_connect('motion_notify_event', self.onMotion)


    def reset(self):
        for ax in self.ax:
            ax.set_xlim(self.max_xlim)
        
        if self.imshow_axes:
            for imax in self.imshow_axes:
                imax.set_xlim(self.max_xlim)


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
        else:
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

        if event.button == 1:
            if self.ticks['left'][0]:
                self.old_tick = self.ticks['left'][0]
                [line.remove() for line in self.ticks['left'][1]]
            else:
                self.old_tick = None

            self.ticks['left'][0] = event.xdata
            self.ticks['left'][1] = []
            for ax in self.ax:
                self.ticks['left'][1].append(ax.axvline(event.xdata, color='k', alpha=0.7, ls='--', lw=1))

            self.cur_xlim = event.inaxes.get_xlim()
            self.new_xlim = ()

            # this is for release
            self.press = self.x0, event.xdata
            self.x0, self.xpress = self.press

        if event.button == 3:
            #right click
            if self.ticks['right'][0]:
                [line.remove() for line in self.ticks['right'][1]]

            self.ticks['right'][0] = event.xdata
            self.ticks['right'][1] = []
            for ax in self.ax:
                self.ticks['right'][1].append(ax.axvline(event.xdata, color='g', alpha=0.7, ls='--', lw=1))


    def onClkRelease(self, event):
        if event.inaxes not in self.ax: 
            return
        
        if self.new_xlim != ():
            # motion was True, remove new_tick
            self.ticks['left'][0] = None
            [line.remove() for line in self.ticks['left'][1]]

            if self.old_tick:
                self.ticks['left'][0] = self.old_tick
                self.ticks['left'][1] = []
                for ax in self.ax:
                    self.ticks['left'][1].append(ax.axvline(self.old_tick, color='k', alpha=0.7, ls='--', lw=1))

        self.press = None
        self.canvas.draw()


    def onMotion(self, event):
        if event.inaxes not in self.ax: 
            return

        if self.press is None:
            return

        dx = event.xdata - self.xpress

        # if dx == 0:
            # remove left click
            # self.ticks['left'][0] = None
            # self.ticks['left'][1].remove()
            
        self.new_xlim = self.cur_xlim - dx
        self.set_xlim(self.new_xlim)



class _Picker(object):
    def __init__(self, axes, phase, sde, row_id, canvas, **kwargs):
        self.db  = sde
        self.row_id = row_id
        self.axes = axes
        self.phase  = phase
        self.canvas = canvas
        self.phase_colors = kwargs.get("phase_colors",{"P":"r", "S":"g", "F":"b"})
        canvas.callbacks.connect('key_press_event', self.on_key)
    

    def clear(self, wave):
        assert wave in ["P", "S", "F"]

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
        
                txt = self.phase_text(wave, self.phase[wave])
                phase['artist_text'] = self.axes[0].annotate(txt, xy=(), color=self.phase_colors[wave])

        self.canvas.draw()


    def update_text(self):
        "do not draw again, just change the artist text of each phase"

        self.phase["P"]["artist_text"].remove()
        txt = self.phase_text("P", self.phase["P"])
        phase['artist_text'] = self.axes[0].annotate(txt)
        self.canvas.draw()


    def save(self):
        idict = {
            'time_P':self.phase["P"]['time'],
            'weight_P':self.phase["P"]['weight'],
            'onset_P':self.phase["P"]['onset'],
            'time_S':self.phase["S"]['time_S'],
            'weight_S':self.phase["S"]['weight'],
            'time_F':self.phase["F"]['time']
        }
        self.sde.update_row(self.row_id, idict)


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
        a = phase_dict.get("onset", None)
        b = phase_dict.get("weigth", None)
        txt = "".join([a,wave,b])
        return txt

