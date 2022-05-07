#!/usr/bin/env python3
# coding=utf-8

import numpy as np
import datetime as dt
from seisvo.file.lte import VECTORAL_PARAMS, VECTORAL_OPT_PARAMS
from seisvo.plotting import plot_gram, get_colors
import matplotlib.ticker as mtick
import matplotlib.dates as mdates

defaultLTE_color = get_colors('zesty')[1]


default_labels = {
    'energy':r'$e$'+'[dB]',
    'specgram':'PSD [dB]',
    'degree':'PD',
    'rect':'R',
    'azm':r'$\Theta_H$ [$\degree$]',
    'dip':r'$\Theta_V$ [$\degree$]',
    'pentropy': 'PE',
    'fq_dominant': r'$f_d$'+'[Hz]',
    'fq_centroid': r'$f_c$'+'[Hz]',
}


class plotLTE(object):
    def __init__(self, fig, lte, start, end, interval, list_attr, **kwargs):
        self.lte = lte
        self.fig = fig
        self.plotkwargs = kwargs
        
        # check times
        if start < lte.stats.starttime:
            raise ValueError('warn: start < lte.stats.starttime')
        
        if end > lte.stats.endtime:
            raise ValueError('warn: end > lte.stats.endtime')

        self.starttime = start
        self.endtime = end
        self.interval = interval

        # check that all attr are available of lte!
        self.list_attr = lte.check_list_attr(list_attr, return_list=True)
        if not self.list_attr:
            raise ValueError('not attr available')

        self.axes_ = {}

        # build the frame
        self.set_frame()
        self.plot()


    def is_twocolumns(self):
        a = [attr for attr in self.list_attr if attr in VECTORAL_PARAMS+VECTORAL_OPT_PARAMS]
        return any(a)
    

    def set_frame(self):
        grid = {'hspace':0.3, 'left':0.08, 'right':0.92, 'wspace':0.1, 'top':0.95, 'bottom':0.05}
        
        if self.is_twocolumns():
            ncols = 2
            grid['width_ratios'] = [1, 0.01]
            grid['wspace'] = 0.01
        else:
            ncols = 1
        
        nrows = len(self.list_attr)
        self.axes = self.fig.subplots(nrows, ncols, gridspec_kw=grid)

        if isinstance(self.axes, np.ndarray):
            self.axes = self.axes.reshape(nrows, ncols)
        else:
            self.axes = np.array([self.axes]).reshape(nrows, ncols)


    def set_xformat(self):

        if self.interval <= 1:
            major_locator = mdates.HourLocator(interval=1)
            major_formatt = mdates.DateFormatter('%d %b\n%H:%M')
            minor_locator = mdates.MinuteLocator(byminute=[15, 30, 45])
            minor_formatt = mtick.NullFormatter()

        elif self.interval <= 3:
            major_locator = mdates.HourLocator(interval=2)
            major_formatt = mdates.DateFormatter('%d\n%H:%M')
            minor_locator = mdates.MinuteLocator(byminute=[15, 30, 45])
            minor_formatt = mtick.NullFormatter()

        elif self.interval <= 10:
            major_locator = mdates.DayLocator(interval=1)
            major_formatt = mdates.DateFormatter('%d %H')
            minor_locator = mdates.HourLocator(byhour=[6, 12, 18, 24])
            minor_formatt = mtick.NullFormatter()

        elif 45 >= self.interval > 10 :
            major_locator = mdates.DayLocator(interval=7)
            major_formatt = mdates.DateFormatter('%d')
            minor_locator = mdates.DayLocator(interval=1)
            minor_formatt = mtick.NullFormatter()

        else:
            major_locator = mdates.WeekdayLocator(interval=2)
            major_formatt = mdates.DateFormatter('%d-%m')
            minor_locator = mdates.DayLocator(interval=7)
            minor_formatt = mtick.NullFormatter()
        
        return (major_locator, major_formatt), (minor_locator, minor_formatt)


    def get_vlim(self, attr):
        if attr == 'specgram':
            v_min = self.plotkwargs.get('specgram_vmin', None)
            v_max = self.plotkwargs.get('specgram_vmax', None)

        if attr == 'degree':
            v_min = 0
            v_max = 1
        
        if attr == 'rect':
            v_min = 0
            v_max = 1

        if attr == 'azm':
            v_min = 0
            v_max = 180
        
        if attr == 'dip':
            v_min = 0
            v_max = 90
        
        return (v_min, v_max)


    def get_ymax(self, attr):
        _, (_, y_max) = self.get_ylim(attr)
        if not y_max:
            data = self.lte.get_attr(attr, self.starttime, self.endtime)
            y = data
            if attr == 'energy':
                y[np.where(y == 0)] = np.nan
                y = 10*np.log10(y)
            y_max = y[np.isfinite(y)].max()
        return y_max 


    def get_ylim(self, attr):
        if attr == 'energy':
            vmin = self.plotkwargs.get('db_min', None)
            vmax = self.plotkwargs.get('db_max', None)

        if attr == 'pentropy':
            vmin = self.plotkwargs.get('pe_min', 0)
            vmax = self.plotkwargs.get('pe_max', 1)

        if attr == 'fq_dominant':
            vmin = self.plotkwargs.get('fd_min', None)
            vmax = self.plotkwargs.get('fd_max', None)

        if attr == 'fq_centroid':
            vmin = self.plotkwargs.get('fc_min', None)
            vmax = self.plotkwargs.get('fc_max', None)
        
        return (vmin, vmax)


    def __show_data__(self, i, attr, show_tickslabels, xformat, **kwargs):
        data = self.lte.get_attr(attr, self.starttime, self.endtime)
        ax = self.axes[i, 0]
        self.axes_[attr].append(ax)
        ax.cla()

        if i == 0 and self.plotkwargs.get('settitle', True):
            ax.set_title('%s %i' % (self.starttime.strftime('%B'), self.starttime.year))

        if self.lte.is_matrix(attr):
            cax = self.axes[i, 1]
            self.axes_[attr].append(cax)
            cax.cla()

            z = data[0].T
            y = data[1]

            if attr == 'specgram':
                z[np.where(z == 0)] = np.nan
                z = 10*np.log10(z)
            
            kwargs['axis_bar'] = cax
            plot_gram(y, z, self.xtime, ax, kwargs)
            ax.set_ylim(self.lte.stats.fq_band)
                    
        else:
            y = data
            if attr == 'energy':
                y[np.where(y == 0)] = np.nan
                y = 10*np.log10(y)

            ax.plot(self.xtime, y, 'k')
            ax.set_ylabel(kwargs.get('y_label'))
            ylim = kwargs.get('ylim', None)

            if ylim:
                ax.set_ylim(ylim)
            else:
                y_min = y[np.isfinite(y)].min()
                y_max = y[np.isfinite(y)].max()
                ax.set_ylim(y_min, y_max)
            
            if self.is_twocolumns():
                self.axes_[i, 1].axes.get_xaxis().set_visible(False)
                self.axes_[i, 1].axes.get_yaxis().set_visible(False)
                self.axes_[i, 1].set_frame_on(False)
            
            # show mode values
            if self.plotkwargs.get('show_mode', True):
                ans = self.lte.get_stats(attr)
                ax.axhline(y=ans[3], color='r', lw=0.5, ls='-.', alpha=0.7, zorder=7)
        
        major_locator, major_formatt = xformat[0]
        minor_locator, minor_formatt = xformat[1]

        if attr == 'energy':
            fmtt = '%i'
        else:
            fmtt = '%.1f'

        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_minor_locator(minor_locator)
        ax.xaxis.set_minor_formatter(mtick.NullFormatter())
        ax.xaxis.set_major_formatter(mtick.NullFormatter())
        ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter(fmtt))

        if show_tickslabels:
            ax.xaxis.set_minor_formatter(minor_formatt)
            ax.xaxis.set_major_formatter(major_formatt)


    def show_events(self, lde, col_dict={}):
        self.events_ = {}
        ids = lde.get_episodes_id(time_interval=(self.starttime, self.endtime))

        if ids:
            for eid in ids:
                episode = lde[eid]
                self.events_[eid] = {'ticks':[], 'event':episode}
                for i, attr in enumerate(self.list_attr):
                    if not self.lte.is_matrix(attr):
                        ax = self.axes_[attr][0]
                        col = col_dict.get(episode.label, defaultLTE_color)
                        v1 = ax.axvline(episode.starttime,  alpha=0.5, color=col, ls='dashed')
                        v2 = ax.axvline(episode.endtime,  alpha=0.5, color=col, ls='dashed')
                        v3 = ax.avspan(episode.starttime,  episode.endtime, alpha=0.15, color=col)
                        self.events_[eid]['ticks'] += [v1,v2,v3]
                        if i == 0:
                            ymax = self.get_ymax(attr)
                            txt = f'ID:{episode.id}[{episode.label}]'
                            mid_bin = episode.starttime + dt.timedelta(hours=episode.duration/2)
                            t = ax.annotate(txt, (mid_bin, ymax), color=col, xycoords='data')
                            self.events_[eid]['text'] = t


    def plot(self):
        num_i = len(self.list_attr)-1
        self.xtime = self.lte.get_time(self.starttime, self.endtime)
        xformat = self.set_xformat()

        for i, attr in enumerate(self.list_attr):
            self.axes_[attr] = []

            if i == num_i:
                show_tickslabels = True
            else:
                show_tickslabels = False

            label = default_labels.get(attr, attr)

            if self.lte.is_matrix(attr):
                ylabel = r'$f$ [Hz]'
                (v_min, v_max) = self.get_vlim(attr)
                self.__show_data__(i, attr, show_tickslabels, xformat, v_max=v_max, v_min=v_min, y_label=ylabel, bar_label=label)
            
            else:
                (vmin, vmax) = self.get_ylim(attr)
                self.__show_data__(i, attr, show_tickslabels, xformat, ylim=(vmin, vmax), y_label=ylabel)
    

    def clear_events(self, id=None):
        def clear_id(eid):
            idict = self.events_.get(eid)
            if idict:
                for artist in idict['ticks']:
                    artist.remove()
                idict['text'].remove()
                del self.events_[eid]
        
        if id:
            clear_id(id)
        
        else:
            for id in self.events_.keys():
                clear_id(id)

    
    def set_text(self, id, txt):
        idict = self.events_.get(id)
        idict['text'].set_text(txt)

