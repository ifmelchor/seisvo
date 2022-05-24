#!/usr/bin/env python3
# coding=utf-8

import numpy as np
import datetime as dt
from seisvo.plotting import plot_gram, get_colors

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates

defaultLTE_color = get_colors('zesty')[1]


default_labels = {
    'energy':r'$e$'+'[dB]',
    'specgram':'PSD [dB]',
    'degree':'PD',
    'rect':'R',
    'azimuth':r'$\Theta_H$ [$\degree$]',
    'elevation':r'$\Theta_V$ [$\degree$]',
    'pentropy': 'PE',
    'fq_dominant': r'$f_d$'+'[Hz]',
    'fq_centroid': r'$f_c$'+'[Hz]',
    'dsar': 'DSAR'
}


class plotLTE(object):
    def __init__(self, chan, fig, lte, start, end, interval, list_attr, **kwargs):
        self.lte = lte
        self.fig = fig
        self.chan = chan
        self.plotkwargs = kwargs
        
        self.starttime = start
        self.endtime = end
        self.interval = interval
        
        self.list_attr = list_attr
        
        self.axes_ = {}
        self.events_ = {}

        # build the frame
        self.set_frame()
        self.plot()


    def any_vector(self):
        ans_list = self.lte.is_attr(self.list_attr, only_vectors=True)
        if ans_list:
            return True
        else:
            return False
    

    def set_frame(self):
        grid = {'hspace':0.15, 'left':0.08, 'right':0.92, 'wspace':0.1, 'top':0.95, 'bottom':0.05}
        
        if self.any_vector():
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
        
        if isinstance(self.xtime[0], dt.datetime):
            if self.interval <= 1:
                major_locator = mdates.HourLocator(interval=1)
                major_formatt = mdates.DateFormatter('%d %b\n%H:%M')
                minor_locator = mdates.MinuteLocator(byminute=[15, 30, 45])
                minor_formatt = mtick.NullFormatter()

            elif self.interval <= 10:
                major_locator = mdates.DayLocator(interval=1)
                major_formatt = mdates.DateFormatter('%d %b %H:%M')
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
        
        else:
            major_locator = mtick.LinearLocator(10)
            major_formatt = mtick.FormatStrFormatter('%i')
            minor_locator = mtick.AutoMinorLocator(2)
            minor_formatt = None
        
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

        if attr == 'azimuth':
            v_min = 0
            v_max = 180
        
        if attr == 'elevation':
            v_min = 0
            v_max = 90
        
        return (v_min, v_max)


    def get_ymax(self, attr):
        _, (_, y_max) = self.get_ylim(attr)
        if not y_max:
            data = self.ddata[attr]
            y = data
            if attr == 'energy':
                y[np.where(y == 0)] = np.nan
                y = 10*np.log10(y)
            y_max = y[np.isfinite(y)].max()
        return y_max 


    def get_ylim(self, attr):
        # vmin, vmax = None, None

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
        data = self.ddata[attr]
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
            plot_gram(y, z, self.xtime, ax, **kwargs)
            ax.set_ylim(self.lte.stats.fq_band)
                    
        else:
            y = data
            if attr == 'energy':
                y[np.where(y == 0)] = np.nan
                y = 10*np.log10(y)

            ax.plot(self.xtime, y, 'k')
            ax.set_xlim(self.xtime[0], self.xtime[-1])
            ax.set_ylabel(kwargs.get('y_label'))
            ylim = kwargs.get('ylim', None)

            if ylim:
                ax.set_ylim(ylim)
            else:
                y_min = y[np.isfinite(y)].min()
                y_max = y[np.isfinite(y)].max()
                ax.set_ylim(y_min, y_max)
            
            if self.any_vector():
                self.axes[i, 1].axes.get_xaxis().set_visible(False)
                self.axes[i, 1].axes.get_yaxis().set_visible(False)
                self.axes[i, 1].set_frame_on(False)
            
            # show mode values
            if self.plotkwargs.get('show_mode', False):
                ans = self.lte.get_stats(attr, chan=self.chan)
                ax.axhline(y=ans[attr][3], color='r', lw=0.5, ls='-.', alpha=0.7, zorder=7)
        
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
            ax.xaxis.set_major_formatter(major_formatt)
            if minor_formatt:
                ax.xaxis.set_minor_formatter(minor_formatt)


    def show_events(self, lde, col_dict={}):
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
        self.xtime, self.ddata = self.lte.get(self.list_attr, chan=self.chan, starttime=self.starttime, endtime=self.endtime)
        
        if self.plotkwargs.get('xaxis_in_hours', True):
            total_hr = (self.endtime - self.starttime).total_seconds()/3600
            self.xtime = np.linspace(0, total_hr, self.lte.stats.nro_time_bins)
        
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
                self.__show_data__(i, attr, show_tickslabels, xformat, ylim=(vmin, vmax), y_label=label)
    

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


# class plotMultiLTE(object)


def plotPeaksSpecPDF(peak, show=True, **kwargs):
    grid = {'left':0.15, 'right':0.95, 'top':0.90, 'bottom':0.2}
    fig, ax = plt.subplots(1,1, figsize=kwargs.get('figsize',(4,3)), gridspec_kw=grid)

    norm = max([info['fq_prob'] for _, info in peak.peaks_.items()])
    
    ax.plot(peak.fq_space_, peak.fq_pdf_/norm, color='k') # plot normalized PDF

    if peak.threshold:
        ax.axhline(y=peak.threshold, color='r', ls='--')
    

    for p_i, info in peak.peaks_.items():
        ax.scatter(info['fq'], info['fq_prob']/norm, color='r', marker='o', ec='k')
        half_width = info['width']
        y1_index = np.argmin(np.abs(peak.fq_space_-(info['fq']-half_width)))
        y2_index = np.argmin(np.abs(peak.fq_space_-(info['fq']+half_width)))
        x_fill = peak.fq_space_[y1_index:y2_index].reshape(-1,)
        y_fill = peak.fq_pdf_[y1_index:y2_index]/norm
        ax.fill_between(x_fill, y_fill, color='k', alpha=0.1)
        # ax.annotate(f'#{p_i}', xy=(info['fq']+0.05, info['fq_prob']+0.01), bbox=dict(boxstyle="round", fc="w", ec="k", lw=0.8), fontsize=8, family='monospace')
    
    ax.set_xlim(peak.fq_space_.min(), peak.fq_space_.max())
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('Normalized PDF($\mathcal{P}$)', fontsize=kwargs.get('fs', 10))
    ax.set_xlabel(r'$f$ [Hz]', fontsize=kwargs.get('fs', 10))

    ax.set_xticks([1,2,3,4,5])
    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(mtick.FixedLocator([0.25,0.5,0.75,1.0]))

    ax.grid(which='major', axis='both', ls='--', alpha=0.5)
    ax.grid(which='minor', axis='both', ls=':', alpha=0.25)

    ax.tick_params(axis='y', labelsize=kwargs.get('fs', 10))
    ax.tick_params(axis='x', labelsize=kwargs.get('fs', 10))

    if show:
        plt.show()

    return fig


def plotPeaksPDF(peak, n, show=True, **kwargs):
    nPeak = peak.peaks_[n]

    PDFs = {}
    if nPeak['sp']:
        sp_space = np.linspace(-170, -100, 1000).reshape(-1,1)
        sp_pdf = np.exp(nPeak['sp']['kde'].score_samples(sp_space))
        PDFs['sp'] = (sp_space, sp_pdf)

    if nPeak['rect']:
        r_space = np.linspace(0, 1, 1000).reshape(-1,1)
        r_pdf = np.exp(nPeak['rect']['kde'].score_samples(r_space))
        PDFs['rect'] = (r_space, r_pdf)
    
    if nPeak['thH']:
        azm_space = np.linspace(0, 180, 1000).reshape(-1,1)
        azm_pdf = np.exp(nPeak['thH']['kde'].score_samples(azm_space))
        PDFs['thH'] = (azm_space, azm_pdf)
    
    if nPeak['thV']:
        dip_space = np.linspace(0, 90, 1000).reshape(-1,1)
        dip_pdf = np.exp(nPeak['thV']['kde'].score_samples(dip_space))
        PDFs['thV'] = (dip_space, dip_pdf)
    
    n_axis = len(PDFs)

    if n_axis == 0:
        print(' nothing to plot')
        return
    
    grid = {'hspace':0.1, 'wspace':0.1, 'left':0.15, 'right':0.95, 'top':0.90, 'bottom':0.2}
    fig, axis = plt.subplots(1, n_axis, figsize=(n_axis*2,2.5), gridspec_kw=grid)

    fs = kwargs.get('fs', 10)
    spec_th = kwargs.get('spec_th', 0.1)

    if not isinstance(axis, np.ndarray):
        axis = [axis]

    for ax, att in zip(axis, PDFs.keys()):
        space, pdf = PDFs[att]
        pdf /= pdf.max()
        ax.plot(space, pdf, color='k')

        if att == 'rect':
            rect_th = peak.peak_thresholds['rect_th']
            ax.fill_between(space[space>rect_th], pdf[np.where(space>rect_th)[0]], alpha=0.1, color='k')
            ax.set_xlim(0,1)
            ax.set_xlabel(r'R', fontsize=fs)
            ax.set_xticks([0.3, 0.5, 0.7])
        
        elif att == 'thH':
            ax.set_xlim(0,180)
            ax.set_xlabel(r'$\Theta_H$ [$\degree$]', fontsize=fs)
            ax.set_xticks([45, 90, 135])

        elif att == 'thV':
            ax.set_xlim(0,90)
            ax.set_xlabel(r'$\Theta_V$ [$\degree$]', fontsize=fs)
            ax.set_xticks([22, 45, 70])
        
        else:
            ax.set_xlim(-170,-110)
            ax.set_xlabel(r'PSD [dB]', fontsize=fs)
            ax.axhline(spec_th, color='k', ls='--')
        
        ax.set_ylim(0,1.1)
            
        ax.tick_params(axis='x', labelsize=fs)
        ax.set_yticks([0.25, 0.5, 0.75])
        ax.grid(which='major', axis='both', ls='--', alpha=0.5)
        ax.grid(which='minor', axis='both', ls=':', alpha=0.25)
        
        if ax == axis[0]:
            ax.tick_params(axis='y', labelsize=fs)
            ax.set_ylabel(r'Normalized PDF', fontsize=fs)
        else:
            ax.yaxis.set_major_formatter(mtick.NullFormatter())
    
    if show:
        plt.show()

    return fig


# revise
# peak should be a list of peaks asoaicted to a different channel, so you can plot all dominant frequencies found in the three component channel.
class plotDPeakTEVO(object):
    def _init__(self, peak, attr_list, fq_range_dict, fig, marker_dict={}, color_dict={}, **kwargs):
        
        self.fig = fig
        self.peak = peak
        self.attr_list = peak.is_attr(attr_list, only_vectors=True)
        self.top_attr = kwargs.get('top_scalar', 'energy')
        self.top_chan = kwargs.get('top_scalar_chan', peak.stats.channel[0])

        self.fq_range_dict = fq_range_dict
        self.marker_dict = marker_dict
        self.color_dict = color_dict

        self.set_frame()
        self.plot()


    def set_frame(self):

        gridspec_dict = {
            'left':0.1,
            'right':0.95,
            'height_ratios':[0.75,1,1],
            'hspace':0.1,
            'wspace':0.1
            }

        nrows = len(self.attr_list) + 1
        self.axes = self.fig.subplots(nrows, 1, gridspec_kw=gridspec_dict)


    def plot(self):

        time, dout = self.peak.get(self.top_attr, chan=self.top_chan, starttime=self.peak.starttime, endtime=self.peak.endtime)

        scalar_data = dout[self.top_attr]
        self.axes[0].plot(time, scalar_data, color='k')
        self.axes[0].set_ylabel(default_labels.get(self.top_attr, self.top_attr))
        prc5 = np.percentile(scalar_data, 5)
        prc95 = np.percentile(scalar_data, 95)
        self.axes[0].set_ylim(prc5, prc95)

        show_ylabel = False

        for group, fq_range in self.fq_range_dict.items():
            dpeaks_dict = self.peak.get_dominant_peaks(fq_range)
            color = self.color_dict.get(group, 'k')
            marker = self.marker_dict.get(group, 'o')

            for t, data_dict in dpeaks_dict.items():
                time_t = time[t]

                for n, attr in enumerate(self.attr_list):
                    data = data_dict[attr]
                    ax = self.axes[n+1]
                    ax.scatter([time_t]*len(data), data, color=color, ec='k', marker=marker)

                    if not show_ylabel:
                        ax.set_ylabel(default_labels.get(attr, attr))

                        if attr == 'azimuth':
                            ax.set_ylim(0,180)
                        
                        if attr in ('degree', 'rect'):
                            ax.set_ylim(0,1)
                        
                        if attr == 'elevation':
                            ax.set_ylim(0,90)
                
                show_ylabel = True

        
        [ax.set_xlim(time[0], time[-1]) for ax in self.axes]
        [ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2)) for ax in self.axes]
        [ax.xaxis.set_major_formatter(mtick.NullFormatter()) for ax in self.axes[:-1]]
        [ax.grid(which='both', axis='both', color='k', alpha=0.35, ls='--', zorder=1) for ax in self.axes]

        # add legend
        # for group, fq_range in self.fq_range_dict.items():


                    







