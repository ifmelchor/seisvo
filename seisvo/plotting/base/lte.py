#!/usr/bin/env python3
# coding=utf-8

import numpy as np
import datetime as dt
from seisvo.plotting import plot_gram, get_colors
from seisvo.file.lte import POLAR_PARAMS

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates

defaultLTE_color = get_colors('zesty')[1]


default_labels = {
    'energy':'Energy\n' + r' [dB//($m^2 s^{-2}$)]',
    'specgram':'PSD\n' + r'[dB//($m^2 s^{-2}$/Hz)]',
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

        # take data and time
        ans = lte.get(list_attr, chan=chan, starttime=start, endtime=end, **kwargs)
        self.xtime, self.ddata = ans

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
        grid = {'hspace':0.15, 'left':0.08, 'right':0.92, 'wspace':0.05, 'top':0.95, 'bottom':0.05}
        
        if self.any_vector():
            ncols = 2
            grid['width_ratios'] = [1, 0.01]
        else:
            ncols = 1
        
        self.nrows = 0
        for attr in self.list_attr:
            if attr in POLAR_PARAMS:
                self.nrows += 1
            else:
                self.nrows += 1*len(self.chan)

        self.axes = self.fig.subplots(self.nrows, ncols, gridspec_kw=grid)


        if isinstance(self.axes, np.ndarray):
            self.axes = self.axes.reshape(self.nrows, ncols)
        else:
            self.axes = np.array([self.axes]).reshape(self.nrows, ncols)


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
            cmap = 'Spectral_r'

        if attr == 'degree':
            v_min = 0
            v_max = 1
            cmap = 'Spectral_r'
        
        if attr == 'rect':
            v_min = 0
            v_max = 1
            cmap = 'Spectral_r'

        if attr == 'azimuth':
            v_min = 0
            v_max = 180
            cmap = 'twilight'
        
        if attr == 'elevation':
            v_min = 0
            v_max = 90
            cmap = 'Spectral_r'
        
        return cmap, (v_min, v_max)


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


    def __show_data__(self, chan, attr, show_tickslabels, xformat, **kwargs):
        
        if len(self.chan) > 1:
            if attr in POLAR_PARAMS:
                data = self.ddata[attr]
                attr_key = attr
            else:
                data = self.ddata[chan][attr]
                attr_key = '/'.join([chan,attr])
        else:
            data = self.ddata[attr]
            attr_key = attr

        if attr_key in list(self.axes_.keys()):
            return
        
        else:
            ax = self.axes[self.i, 0]
            self.axes_[attr_key] = ax
            ax.cla()

            if self.i == 0 and self.plotkwargs.get('settitle', True):
                ax.set_title('%s %i' % (self.starttime.strftime('%B'), self.starttime.year))

            if self.lte.is_matrix(attr):
                cax = self.axes[self.i, 1]
                self.axes_[attr_key+'/cax'] = cax
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
                    self.axes[self.i, 1].axes.get_xaxis().set_visible(False)
                    self.axes[self.i, 1].axes.get_yaxis().set_visible(False)
                    self.axes[self.i, 1].set_frame_on(False)
                
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
                ax.set_xlabel('Time [h]')
                if minor_formatt:
                    ax.xaxis.set_minor_formatter(minor_formatt)
            
            # change axis
            self.i += 1


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
        
        xformat = self.set_xformat()
        self.i = 0

        for chan in self.chan:
            for attr in self.list_attr:
                if self.i == self.nrows-1:
                    show_tickslabels = True
                else:
                    show_tickslabels = False

                label = default_labels.get(attr, attr)

                if self.lte.is_matrix(attr):
                    ylabel = r'$f$ [Hz]'
                    cmap, (v_min, v_max) = self.get_vlim(attr)
                    self.__show_data__(chan, attr, show_tickslabels, xformat, v_max=v_max, v_min=v_min, cmap=cmap, y_label=ylabel, bar_label=label)
                
                else:
                    (vmin, vmax) = self.get_ylim(attr)
                    self.__show_data__(chan, attr, show_tickslabels, xformat, ylim=(vmin, vmax), y_label=label)
                

    

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
    fig, ax = plt.subplots(1,1, figsize=kwargs.get('figsize',(9,3)), gridspec_kw=grid)

    norm = max([info['fq_prob'] for _, info in peak.peaks_.items()])
    
    ax.plot(peak.fq_space_, peak.fq_pdf_/norm, color='k') # plot normalized PDF

    if peak.threshold:
        ax.axhline(y=peak.threshold, color='r', ls='--')
    
    for _, info in peak.peaks_.items():
        ax.scatter(info['fq'], info['fq_prob']/norm, color='r', marker='o', ec='k')
        half_width = info['width']
        y1_index = np.argmin(np.abs(peak.fq_space_-(info['fq']-half_width)))
        y2_index = np.argmin(np.abs(peak.fq_space_-(info['fq']+half_width)))
        x_fill = peak.fq_space_[y1_index:y2_index].reshape(-1,)
        y_fill = peak.fq_pdf_[y1_index:y2_index]/norm
        ax.fill_between(x_fill, y_fill, color='k', alpha=0.1)
        # ax.annotate(f'#{p_i}', xy=(info['fq']+0.05, info['fq_prob']+0.01), bbox=dict(boxstyle="round", fc="w", ec="k", lw=0.8), fontsize=8, family='monospace')
    
    fs = kwargs.get('fs', 10)
    ax.set_xlim(peak.fq_space_.min(), peak.fq_space_.max())
    ax.set_ylim(0, 1.1)
    ax.set_ylabel('$P(\mathcal{P}_{wp}$)', fontsize=fs)
    ax.set_xlabel(r'$f$ [Hz]', fontsize=fs)

    ax.set_xticks([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(mtick.FixedLocator([0, 0.5, 1.0]))

    ax.grid(which='major', axis='both', ls='--', alpha=0.5)
    ax.grid(which='minor', axis='both', ls=':', alpha=0.25)

    ax.tick_params(axis='y', labelsize=fs)
    ax.tick_params(axis='x', labelsize=fs)

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
    
    grid = {'hspace':0.1, 'wspace':0.3, 'left':0.15, 'right':0.95, 'top':0.90, 'bottom':0.2}
    fig, axis = plt.subplots(1, n_axis, figsize=(n_axis*2,2.5), gridspec_kw=grid)

    fs = kwargs.get('fs', 10)

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
            ax.set_xticks([0, 0.5, 1])
            ax.set_ylabel(r'$P(\mathcal{R}_%s$)' %n, fontsize=fs)
            ax.set_title(r'$N_R$ = %s' % nPeak['rect']['n'])
        
        elif att == 'thH':
            ax.set_xlim(0,180)
            ax.set_xlabel(r'$\Theta_H$ [$\degree$]', fontsize=fs)
            ax.set_xticks([0, 90, 180])
            ax.set_ylabel(r'$P(\mathcal{H}_%s$)' %n, fontsize=fs)
            ax.set_title(r'$N_H$ = %s' % nPeak['thH']['n'])

        elif att == 'thV':
            ax.set_xlim(0,90)
            ax.set_ylabel(r'$P(\mathcal{V}_%s$)' %n, fontsize=fs)
            ax.set_xlabel(r'$\Theta_V$ [$\degree$]', fontsize=fs)
            ax.set_xticks([0, 45, 90])
            ax.set_title(r'$N_V$ = %s' % nPeak['thV']['n'])
        
        else:
            ax.set_xlim(-170,-110)
            ax.set_xlabel('PSD\n' +r'[dB//(m$^2$s$^{-2}$/Hz)]', fontsize=fs)
            ax.set_ylabel(r'$P(\mathcal{S}_%s$)' %n, fontsize=fs)
            ax.set_title(r'$N_S$ = %s' % nPeak['sp']['n'])
        
        ax.set_ylim(0,1.1)
        ax.tick_params(axis='x', labelsize=fs)
        ax.set_yticks([0, 0.5, 1])
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.grid(which='major', axis='both', ls='--', alpha=0.5)
        ax.grid(which='minor', axis='both', ls=':', alpha=0.25)
        
        if ax == axis[0]:
            ax.tick_params(axis='y', labelsize=fs)
        else:
            ax.yaxis.set_major_formatter(mtick.NullFormatter())
    
    if show:
        plt.show()

    return fig


class plotDPeakTEVO(object):
    def __init__(self, lte, starttime, endtime, chan_list, attr_list, fq_range_dict, fig, **kwargs):

        if not fig:
            self.fig = plt.figure(figsize=kwargs.get('figsize', (8,9)))
        else:
            self.fig = fig
        
        self.lte = lte
        self.attr_list = lte.is_attr(attr_list, only_vectors=True) # move to main call
        self.top_attr = kwargs.get('top_scalar', 'energy')
        self.top_chan = kwargs.get('top_scalar_chan', chan_list[0])

        self.chan_list = chan_list
        
        self.fqs = fq_range_dict
        
        self.starttime = starttime
        self.endtime = endtime

        self.set_frame()
        self.plot(**kwargs)

        if kwargs.get('title', None):
            self.axes[0].set_title(kwargs.get('title'))

        self.fig.align_ylabels()


    def set_frame(self):
        gridspec_dict = {
            'left':0.1,
            'right':0.95,
            'height_ratios':[0.75]+[1]*len(self.attr_list),
            'hspace':0.1,
            'wspace':0.1
            }
        nrows = len(self.attr_list) + 1
        self.axes = self.fig.subplots(nrows, 1, gridspec_kw=gridspec_dict)
        self.axes[-1].set_xlabel('Time [hr]')


    def plot(self, **kwargs):

        time, dout = self.lte.get(self.top_attr, chan=self.top_chan, starttime=self.starttime, endtime=self.endtime, **kwargs)

        if self.top_attr == 'energy':
            scalar_data = 10*np.log10(dout[self.top_attr])
        else:
            scalar_data = dout[self.top_attr]

        self.axes[0].plot(time, scalar_data, color='k')

        # define the limits
        y_lim = kwargs.get('y_lim', None)
        y_ticks = kwargs.get('y_ticks', None)
        y_label = kwargs.get('y_label', True)

        if y_label:
            self.axes[0].set_ylabel(default_labels.get(self.top_attr, self.top_attr))
        
        if y_lim:
            self.axes[0].set_ylim(y_lim)
        else:
            prc5 = np.percentile(scalar_data[np.isfinite(scalar_data)], 1)
            prc95 = np.percentile(scalar_data[np.isfinite(scalar_data)], 99)
            self.axes[0].set_ylim(prc5, prc95)
        
        if y_ticks:
            self.axes[0].set_yticks(y_ticks)

        self.axes[0].yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        self.axes[0].xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        
        show_ylabel = True
        # add legend
        legend_ = {'h':[], 't':[]}
        for chan in self.chan_list:
                peak = self.lte.get_Peaks(chan, starttime=self.starttime, endtime=self.endtime, peak_thresholds=kwargs.get('peak_thresholds', {}))
                
                for fq, fq_dict in self.fqs.items():
                    fq_range = [fq-fq_dict['width'], fq+fq_dict['width']]
                    fq_peaksdict = peak.get_dominant_peaks(fq_range)
                    fq_color = fq_dict.get('color', 'k')
                    fq_marker = fq_dict.get('marker', 'o')
                    fq_label = fq_dict.get('label', str(round(fq, 1)))

                    for t, data_dict in fq_peaksdict.items():
                        time_t = time[t]
                
                        for n, attr in enumerate(self.attr_list):
                            data = data_dict[attr]
                            ax = self.axes[n+1]

                            if data:
                                h = ax.scatter([time_t]*len(data), data, color=fq_color, ec='k', marker=fq_marker, alpha=0.5)

                                if fq_label not in legend_['t']:
                                    legend_['h'] += [h]
                                    legend_['t'] += [fq_label]

                            if show_ylabel:

                                if y_label:
                                    ax.set_ylabel(default_labels.get(attr, attr))

                                if attr == 'azimuth':
                                    ax.set_ylim(0, 180)
                                    ax.invert_yaxis()
                                    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
                                    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
                            
                                if attr in ('degree', 'rect'):
                                    ax.set_ylim(0, 1)
                                    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
                                    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
                            
                                if attr == 'elevation':
                                    ax.set_ylim(0, 90)
                                    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
                                    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
                    
                        show_ylabel = False

        [ax.set_xlim(time[0], time[-1]) for ax in self.axes]
        [ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2)) for ax in self.axes]
        [ax.xaxis.set_major_formatter(mtick.NullFormatter()) for ax in self.axes[:-1]]
        [ax.grid(which='both', axis='both', color='k', alpha=0.35, ls='--', zorder=1) for ax in self.axes]

        self.fig.legend(legend_['h'], legend_['t'], title=r'$f$ [Hz]', ncol=1)




                    







