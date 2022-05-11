#!/usr/bin/env python
# coding=utf-8

import math, os
import numpy as np
import datetime as dt
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
import matplotlib.colors as mcolor
import matplotlib._color_data as mcd

plt.rc('axes', labelsize=10)
plt.rc('axes', labelpad=4.0)
plt.rc('axes', titlepad=6.0)
plt.rc('axes', titlesize=10)
plt.rc('xtick', labelsize=10)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('lines', linewidth=0.5)
plt.rc('lines', linewidth=0.5)

def truncate_cmap(cmap, minval=0.0, maxval=1.0, n=100):
    '''
    https://stackoverflow.com/a/18926541
    '''
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    new_cmap = mcolor.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    
    return new_cmap


def count_consec(listrand):
    count=1
    consec_list=[]

    first = False
    for i in range(len(listrand[:-1])):
        if first is False:
            first=listrand[i]
        
        if listrand[i]+1 == listrand[i+1]:
            count+=1
        
        else:
            consec_list.append((first, count))
            first = False
            count=1

    # Account for the last iteration
    if first is False:
        first=listrand[-1]
    
    consec_list.append((first, count))
    return consec_list


def pplot_control(title, date_list, nro_traces, sample_rate, npts, max_value):

    fig, axs = plt.subplots(4,1, constrained_layout=True, figsize=(16,9), sharex=True)
    fig.suptitle(title, fontsize=12)

    mk = '+'
    axs[0].scatter(date_list, nro_traces, c='r', marker=mk)
    axs[0].set_ylabel('nro_traces')
    y_min = min(filter(lambda k: k is not None, nro_traces))
    y_max = max(filter(lambda k: k is not None, nro_traces))
    axs[0].set_ylim(y_min - 1, y_max + 1)
    axs[0].axes.get_xaxis().set_visible(False)

    axs[1].scatter(date_list, sample_rate, c='darkgreen', marker=mk)
    axs[1].set_ylabel('sample_rate')
    y_min = min(filter(lambda k: k is not None, sample_rate))
    y_max = max(filter(lambda k: k is not None, sample_rate))
    axs[1].set_ylim(y_min - 1, y_max + 1)
    axs[1].axes.get_xaxis().set_visible(False)
    
    axs[2].scatter(date_list, npts, c='grey', marker=mk)
    axs[2].set_ylabel('npts_total')
    y_min = min(filter(lambda k: k is not None, npts))
    y_max = max(filter(lambda k: k is not None, npts))
    delta = (y_max*10/100)
    axs[2].set_ylim(y_min - delta, y_max + delta)
    axs[2].axes.get_xaxis().set_visible(False)

    axs[3].scatter(date_list, max_value, c='k', marker=mk)
    axs[3].set_ylabel('max_count')
    y_min = min(filter(lambda k: k is not None, max_value))
    y_max = max(filter(lambda k: k is not None, max_value))
    axs[3].set_ylim(y_min - 1, y_max + 1)
    axs[3].set_xlim([date_list[0], date_list[-1]])
    axs[3].xaxis.set_major_formatter(mdates.DateFormatter('%Y.%j'))

    plt.show()


def plot_check(*args):

    title = args[0]
    sta_list = args[1]
    availability = args[2]
    day_list = args[3]
    
    fig = plt.figure(figsize=(10,7))

    ax = fig.add_subplot(111)
    n = len(sta_list)
    ticks = range(1,2*n,2)
    pos = range(0,availability.shape[1])
    width = 0.5
    major_bin = 5
    minor_bin = 1

    for l, k in zip(ticks, pos):
        klist = availability[:,k]

        green = list(np.where(klist==1)[0])
        if green:
            green = count_consec(green)
            ax.broken_barh(green, (l-(width/2), width), facecolors=('tab:green'))
    
        red = list(np.where(klist==0)[0])
        if red:
            red = count_consec(red)
            ax.broken_barh(red, (l-(width/2), width), facecolors=('tab:red'))

    ax.set_ylim(0,2*n)
    ax.set_yticks(ticks)

    # remove net code of the ticklabels
    yticks = []
    for y in sta_list:
        yticks.append('.'.join(y.split('.')[1:]))
    ax.set_yticklabels(yticks)

    ax.set_xlim(0, len(day_list))
    ax.xaxis.set_major_locator(mtick.MultipleLocator(major_bin))
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(minor_bin))

    ax.grid(axis='x', which='major', color='k', linestyle='-', alpha=0.7)
    ax.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.5)

    xticklabels = []
    for d in day_list[::major_bin]:
        xticklabels.append('%s.%s' % (d.year, d.timetuple().tm_yday))

    ax.set_xticks(range(0,len(day_list),major_bin))
    ax.set_xticklabels(xticklabels)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.set_title(title)
    fig.tight_layout()
    plt.show()


def plot_tr_detection(dict_out, tremors=None, figsize=(17,5), plot=True, fill_white=[]):
    fig, (ax1, ax3, ax2) = plt.subplots(3,1, figsize=figsize, gridspec_kw=dict(hspace=0.2))

    major_locator = mdates.MonthLocator()
    minor_locator = mdates.WeekdayLocator()
    major_formatt = mdates.DateFormatter('%b')
    minor_formatt = mtick.NullFormatter()

    time = np.array(dict_out['time'])
    r_eh = np.array(dict_out['r_eh'])
    pavg = np.array(dict_out['pavg'])
    pmax = np.array(dict_out['pmax'])

    if fill_white:
        for item in fill_white:
            bin1 = np.argmin(np.abs(time-item[0]))
            bin2 = np.argmin(np.abs(time-item[1]))
            r_eh[bin1:bin2] = np.nan
            pavg[bin1:bin2] = np.nan
            pmax[bin1:bin2] = np.nan

    ax1.plot(time, r_eh, color='k')
    h1 = ax1.axhline(y=dict_out['r_eh_thr'], color='red', lw=0.5, ls='--', alpha=0.5)
    ax1.set_ylabel(r'$r_{eh}$')
    ax1.set_ylim(-1,1)
    ax1.xaxis.set_major_locator(major_locator)
    ax1.xaxis.set_major_formatter(minor_formatt)
    ax1.xaxis.set_minor_locator(minor_locator)
    # ax1.yaxis.set_major_locator(mtick.LinearLocator(numticks=2))
    # ax1.yaxis.set_major_formatter(mtick.StrMethodFormatter("{x:.2f}"))
    ax1.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax1.grid(axis='x', which='major', color='k', linestyle='--', alpha=0.3)
    ax1.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.1)
    ax1.set_xlim(time[0], time[-1])

    ax2.plot(time, pavg, color='k')
    ax2.axhline(y=dict_out['pavg_thr'], color='red', lw=0.5, ls='--', alpha=0.5)
    ax2.set_ylabel(r'$\overline{p_{avg}}$')
    ax2.set_ylim(0,1)
    ax2.xaxis.set_major_locator(major_locator)
    ax2.xaxis.set_major_formatter(major_formatt)
    ax2.xaxis.set_minor_locator(minor_locator)
    # ax2.yaxis.set_major_locator(mtick.LinearLocator(numticks=2))
    # ax2.yaxis.set_major_formatter(mtick.StrMethodFormatter("{x:.2f}"))
    ax2.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax2.grid(axis='x', which='major', color='k', linestyle='--', alpha=0.3)
    ax2.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.1)
    ax2.set_xlim(time[0], time[-1])
    ax2.set_xlabel('2012')


    ax3.plot(time, pmax, color='k')
    ax3.axhline(y=dict_out['pmax_thr'], color='red', lw=0.5, ls='--', alpha=0.5)
    ax3.set_ylabel(r'$\overline{p_{d}}$')
    ax3.set_ylim(0,1)
    ax3.xaxis.set_major_locator(major_locator)
    ax3.xaxis.set_major_formatter(minor_formatt)
    ax3.xaxis.set_minor_locator(minor_locator)
    # ax3.yaxis.set_major_locator(mtick.LinearLocator(numticks=2))
    # ax3.yaxis.set_major_formatter(mtick.StrMethodFormatter("{x:.2f}"))
    ax3.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax3.grid(axis='x', which='major', color='k', linestyle='--', alpha=0.3)
    ax3.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.1)
    ax3.set_xlim(time[0], time[-1])

    handles = [(h1, 'Threshold')]

    if tremors:
        for tr in tremors:
            for ax in [ax1, ax2, ax3]:
                h2 = ax.axvspan(tr[0], tr[0]+tr[1], alpha=0.2, color='green')
                if not 'Tremor' in [h[1] for h in handles]:
                    handles += [(h2, 'Tremor')]


    fig.legend(handles=[x1[0] for x1 in handles], labels=[x1[1] for x1 in handles],
        loc='upper right', ncol=len(handles))

    if plot:
        plt.show()
    
    return fig


def plot_multiple_psd(freq, psd_array, pd_array=None, **kwargs):
    fig = kwargs.get('fig', None)
    figsize = kwargs.get('figsize', (8, 4))
    dpi = kwargs.get('dpi', 100)
    plot = kwargs.get('plot', True)
    plot_prob = kwargs.get('plot_prob', True)
    log_psd = kwargs.get('log_psd', False)
    log_fq = kwargs.get('log_fq', False)
    colors = kwargs.get('colors', [])

    if isinstance(pd_array, np.ndarray):
        rows = 2
    else:
        rows = 1
    
    grid = {'hspace':0.1, 'left':0.12, 'right':0.95, 'wspace':0.1, 'top':0.90, 'bottom':0.1}

    if fig:
        axes = fig.subplots(rows, 1, gridspec_kw=grid)

    else:
        fig, axes = plt.subplots(rows, 1, gridspec_kw=grid, figsize=figsize, dpi=dpi)

    if isinstance(pd_array, np.ndarray):
        psd_axis, pd_axis = axes
    else:
        psd_axis = axes

    # psd_array = 10*np.log(psd_array)

    if plot_prob:
        psd_prob = get_max_prob(psd_array)
        psd_axis.plot(freq, psd_prob, 'r', ls='--', lw=2, zorder=10)

    if colors:
        for i in range(psd_array.shape[0]):
            psd_axis.plot(freq, psd_array[i,:], color=colors[i])
    else:
        psd_axis.plot(freq, psd_array.T, color='k')
    psd_axis.set_xlim(min(freq), max(freq))

    if plot_prob:
        psd_axis.set_ylim(-10, psd_prob.max()+20)

    psd_axis.set_ylabel(r'PSD [$m^2 \cdot s^{-2}/Hz$] dB')
    psd_axis.set_title('Number of curves: %i' % psd_array.shape[0])
    psd_axis.yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
    psd_axis.xaxis.set_minor_locator(mtick.AutoMinorLocator(3))
    psd_axis.grid(axis='x', which='major', color='k', linestyle='--', alpha=0.3)
    psd_axis.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.1)

    # if log_fq:
    #     psd_axis.set_xscale('log')

    if isinstance(pd_array, np.ndarray):
        psd_axis.xaxis.set_major_formatter(mtick.NullFormatter())

        if plot_prob:
            pd_prob = get_max_prob(pd_array)
            pd_axis.plot(freq, pd_prob, 'r', ls='--', lw=2, zorder=10)
        
        pd_axis.plot(freq, pd_array.T, 'k')
        pd_axis.set_xlim(min(freq), max(freq))
        pd_axis.set_ylabel('PD')
        pd_axis.yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        pd_axis.xaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        pd_axis.grid(axis='x', which='major', color='k', linestyle='--', alpha=0.3)
        pd_axis.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.1)

        if log_fq:
            pd_axis.set_xscale('log')

    if plot:
        plt.show()

    return fig, axes


def plot_pdf(pdf, y_bins, x_bins, axis=None, plot=True, **kwargs):
    fig = None

    if not axis:    
        figsize = kwargs.get('figsize', (8, 4))
        grid = {'hspace':0.3, 'left':0.12, 'right':0.90, 'wspace':0.1, 'top':0.95, 'bottom':0.15, 'width_ratios':[1, 0.02]}
        fig, axes = plt.subplots(1, 2, gridspec_kw=grid, figsize=figsize, dpi=100)
        axis = axes[0]
        kwargs['axis_bar'] = axes[1]

    kwargs['y_label'] = 'PSD\n'+r'dB[cnts$^2$/Hz]'
    kwargs['x_label'] = 'Freq. [Hz]'
    kwargs['bar_label'] = 'PDF'

    plot_gram(y_bins, pdf, x_bins, axis, **kwargs)

    if kwargs.get('show_models', False):
        from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

        _, nlnm = get_nlnm() # NLNM model
        p, nhnm = get_nhnm() # NHNM model

        axis.plot(1/p, nlnm, color='k', lw=1.2, label='NLNM', zorder=10)
        axis.plot(1/p, nhnm, color='k', lw=1.2, label='NHNM', zorder=10)

    if plot:
        plt.show()

    return fig


def plot_gram(y, array, x, axis, return_bar=False, **kwargs):
    """
    This code was designer to plot both spectrogram and polargram
    """

    # check input
    if isinstance(array, np.ndarray):
        if len(array.shape) == 2:
            pass
        else:
            raise ValueError('array must be a ndarray of M,N')
    else:
        raise ValueError('array must be a ndarray of M,N')

    if isinstance(x[0], dt.datetime):
        is_time = True
    else:
        is_time = False

    if is_time:
        x1 = mdates.date2num(x[0])
        x2 = mdates.date2num(x[1])
        xFI = mdates.date2num(x[-1])
    else:
        x1 = x[0]
        x2 = x[1]
        xFI = x[-1]

    halfbin_x = (x2 - x1) / 2.0
    halfbin_y = (y[1] - y[0]) / 2.0
    extent = (
        x1 + halfbin_x,
        xFI - halfbin_x,
        y[0] + halfbin_y,
        y[-1] - halfbin_y
    )
    rel_norm = kwargs.get('rel_norm', False)
    v_max = kwargs.get('v_max', None)
    v_min = kwargs.get('v_min', None)

    if rel_norm:
        v_max = 1
        v_min = 0
        array_norm = np.ones(shape=array.shape)
        for k in range(array.shape[1]):
            array_norm[:, k] = array[:, k] / array[:, k][np.isfinite(array[:, k])].max()
        array = array_norm

    else:
        if not v_max:
            v_max = np.nanmax(array[np.isfinite(array)])
    
        if not  v_min:
            v_min = np.nanmin(array[np.isfinite(array)])

    
    interpolation = kwargs.get('interpolation', 'gaussian')
    cmap = kwargs.get('cmap', 'Spectral_r')
    norm1 = kwargs.get('norm', None)

    if not norm1:
        norm1 = mcolor.Normalize(v_min, v_max)

    im = axis.imshow(np.flipud(array), cmap=cmap, norm=norm1, interpolation=interpolation, extent=extent, aspect='auto')
    axis.axis('tight')
    axis.grid(False)
    
    if is_time:
        axis.xaxis_date()

    fq_logscale = kwargs.get('fq_logscale', False)
    if fq_logscale:
        axis.set_yscale('log')

    label_size = kwargs.get('label_size', 10)
    axis_y_label = kwargs.get('y_label', None)
    axis_x_label = kwargs.get('x_label', None)
    axis.set_xlabel(axis_x_label, size=label_size)
    axis.set_ylabel(axis_y_label, size=label_size)

    axis_bar = kwargs.get('axis_bar', None)
    if axis_bar:
        axis_bar_label = kwargs.get('bar_label', None)
        fig = axis.get_figure()
        cbar = fig.colorbar(im, cax=axis_bar, orientation='vertical')
        cbar.locator = mtick.MaxNLocator(nbins=4)
        cbar.update_ticks()
        cbar.set_label(axis_bar_label, size=label_size)

        if return_bar:
            return cbar

    return im, (v_min, v_max)


def polargram_plot(polar_degree, freq, trace=None, fq_band=(), opt_param=None, savefig=False, c_map='inferno'):
    """
    Create the figure frame
    """

    def get_norm(cmap, matrix, v_min=None, v_max=None):

        if not v_min:
            v_min = np.nanmin(matrix[np.isfinite(matrix)])
        if not v_max:
            v_max = np.nanmax(matrix[np.isfinite(matrix)])
        
        levels = mtick.MaxNLocator(nbins=500).tick_values(v_min, v_max)
        norm = mcolor.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        return norm, (v_min, v_max)

    def set_ax(ax, fq_band):
        ax.set_ylabel("Freq. (Hz)")
        lims2 = ax.axis('tight')
        ax.set_xlim(0, lims2[1])
        if fq_band != ():
            ax.set_ylim(fq_band)
        ax.grid(False)
        ax.axes.get_xaxis().set_visible(False)

    polar_dgr = np.flipud(polar_degree)

    n = 1
    
    if opt_param:
        try:
            azm = np.flipud(opt_param['azm'])
            plot_azmith = True
            n += 1
        except:
            plot_azmith = False

        try:
            dip = np.flipud(opt_param['dip'])
            plot_dip = True
            n += 1
        except:
            plot_dip = False

        try:
            plot_rect = True
            rect = np.flipud(opt_param['rect'])
            n += 1
        except:
            plot_rect = False

    if trace:
        n += 2

    cmap = plt.get_cmap(c_map)
    xbins = range(polar_dgr.shape[1])
    ybins = freq
    halfbin_time = (xbins[1] - xbins[0]) / 2.0
    halfbin_freq = (ybins[1] - ybins[0]) / 2.0
    extent = (xbins[0] - halfbin_time,
              xbins[-1] + halfbin_time,
              ybins[0] - halfbin_freq,
              ybins[-1] + halfbin_freq)

    fig, axs = plt.subplots(n,2, 
        gridspec_kw={"width_ratios":[1, 0.05],
                     "wspace":0.05},
        constrained_layout=False,
        figsize=(16,9))
        
    i = iter(range(n))
    
    # trace/specgram
    if trace:
        tr_i = next(i)
        ax_trace = axs[tr_i,0]
        time = trace.get_time()
        ax_trace.plot(time, trace.data)
        starttime = trace.stats.starttime.datetime.strftime("%Y-%m-%d %H:%M:%S")
        endtime = trace.stats.endtime.datetime.strftime("%Y-%m-%d %H:%M:%S")
        ax_trace.set_title('%s -- %s' % (starttime, endtime))
        ax_trace.set_ylabel(trace.id)
        ax_trace.set_xlim(time[0], time[-1])
        axs[tr_i,1].axes.get_xaxis().set_visible(False)
        axs[tr_i,1].axes.get_yaxis().set_visible(False)
        axs[tr_i,1].axis('off')

        # specgram
        i_spec = next(i)
        ax_spec = axs[i_spec,0]
        ans = trace.specgram(axes=ax_spec, fq_band=fq_band)
        ax_spec.set_ylabel("Freq. (Hz)")
        ax_spec.axes.get_xaxis().set_visible(False)
        ax_spec_color = axs[i_spec,1]
        fig.colorbar(ans[0], cax=ax_spec_color, orientation='vertical')
        
    # polargram
    i_poldgr = next(i)
    ax_poldgr = axs[i_poldgr,0]
    norm, (v_min, v_max) = get_norm(cmap, polar_dgr, v_min=0, v_max=1)
    cax_poldgr = ax_poldgr.imshow(polar_dgr, cmap=cmap, norm=norm, interpolation='gaussian', extent=extent)
    set_ax(ax_poldgr, fq_band)
    ax_poldgr_color = axs[i_poldgr,1]
    cb_poldgr = fig.colorbar(cax_poldgr, cax=ax_poldgr_color, ticks=[v_min, (v_max-v_min)/2, v_max], orientation='vertical')
    cb_poldgr.set_label('Polarization degree')

    if opt_param:
        # azimuth
        if plot_azmith:
            i_azm = next(i)
            ax_azm = axs[i_azm,0]
            norm, (v_min, v_max) = get_norm(cmap, azm, v_min=-180, v_max=180)
            cax_azm = ax_azm.imshow(azm, cmap=cmap, norm=norm, interpolation='gaussian', extent=extent)
            set_ax(ax_azm, fq_band)
            ax_azm_color = axs[i_azm,1]
            cb_azm = fig.colorbar(cax_azm, cax=ax_azm_color, ticks=[v_min, (v_max-v_min)/2, v_max], orientation='vertical')
            cb_azm.set_label('azimuth')

        # dip
        if plot_dip:
            i_dip = next(i)
            ax_dip = axs[i_dip,0]
            norm, (v_min, v_max) = get_norm(cmap, dip, v_min=0, v_max=90)
            cax_dip = ax_dip.imshow(dip, cmap=cmap, norm=norm, interpolation='gaussian', extent=extent)
            set_ax(ax_dip, fq_band)
            ax_dip_color = axs[i_dip,1]
            cb_dip = fig.colorbar(cax_dip, cax=ax_dip_color, ticks=[v_min, (v_max-v_min)/2, v_max], orientation='vertical')
            cb_dip.set_label('dip')

        # rect
        if plot_rect:
            i_rect = next(i)
            ax_rect = axs[i_rect,0]
            norm, (v_min, v_max) = get_norm(cmap, rect, v_min=0, v_max=1)
            cax_rect = ax_rect.imshow(rect, cmap=cmap, norm=norm, interpolation='gaussian', extent=extent)
            set_ax(ax_rect, fq_band)
            ax_rect_color = axs[i_rect,1]
            cb_rect = fig.colorbar(cax_rect, cax=ax_rect_color, ticks=[v_min, (v_max-v_min)/2, v_max], orientation='vertical')
            cb_rect.set_label('rectiliniarity')

    cmap.set_bad(color='azure', alpha = 0.5)

    if savefig:
        fig_name = 'polar_%s_%s.png' % (time[0].strftime('%m-%d-%Y'), time[-1].strftime('%m-%d-%Y'))
        fig.savefig(os.path.join(savefig, fig_name), dpi=100, format='png')
        fig.clf()
        plt.close()
    else:
        plt.show()


def plot_activity(lde, afile, starttime, endtime, day_interval_plot=90, **kwargs):

    from seisvo.gui.glde import __get_color__

    pos = dict(ash=(0.5, 0.08), 
           so2=(0.4, 0.08), 
           temp_anomaly=(0.3, 0.08), 
           incandescence=(0.3, 0.08), 
           explosions=(0.2,0.1))

    colors = dict(ash=mcd.XKCD_COLORS["xkcd:black"].upper(),
                  so2=mcd.XKCD_COLORS["xkcd:green"].upper(),
                  temp_anomaly=mcd.XKCD_COLORS["xkcd:pink"].upper(),
                  incandescence=mcd.XKCD_COLORS["xkcd:red"].upper(),
                  explosions=mcd.XKCD_COLORS["xkcd:goldenrod"].upper(),)

    marker_edge = dict(ash=None,
                       so2=None,
                       temp_anomaly=None,
                       incandescence=None,
                       explosions='black')

    marker = dict(ash='s',
                  so2='s',
                  temp_anomaly='s',
                  incandescence='s',
                  explosions='*')

    marker_size = dict(ash=50,
                  so2=50,
                  temp_anomaly=50,
                  incandescence=50,
                  explosions=100)

    title = dict(ash='Ash',
                 so2=r'S0$_2$'+'anomaly',
                 temp_anomaly='Thermal anomaly',
                 incandescence='Incandescence',
                 explosions='Explosion')


    nro_days = (endtime-starttime).total_seconds()/(3600*24)
    nro_plots = math.ceil(nro_days/day_interval_plot)

    grid_opt = dict(left=0.1, right=0.9, top=0.9, bottom=0.1)
    figsize = kwargs.get('figsize', (11, nro_plots*5/3))
    fig, axes = plt.subplots(nro_plots, 1, figsize=figsize, gridspec_kw=grid_opt)

    start = starttime
    handles1 = []
    handles2 = []
    alpha_event = kwargs.get('alpha_event', 0.8)
    alpha_activity = kwargs.get('alpha_activity', 1)

    for n, ax in enumerate(axes):
        end = start + dt.timedelta(days=day_interval_plot) - \
                dt.timedelta(minutes=20)

        if end > endtime:
            end = endtime
        
        for evnt in lde.get_events(time_interval=(start, end)):
            if evnt.label in ('BBT', 'MPT', 'SPT'):
                event_name = evnt.label.lower()
            else:
                event_name = evnt.label
            he = ax.broken_barh([(evnt.starttime,  dt.timedelta(hours=evnt.duration))],
                                (-1.2, 0.8), alpha=alpha_event, ec=None, facecolor=__get_color__(evnt.label))
            mid_time = evnt.starttime + dt.timedelta(hours=evnt.duration/2)
            #ax.annotate(evnt.id, (mid_time, 1.1), color=__get_color__(evnt.label), fontsize=9)
            
            if not handles1 or event_name not in [h[1] for h in handles1]:
                handles1 += [(he, event_name)]
        
        act = afile.get_activity_dict(start, end)
        for key, times in act.items():
            if times:
                for titem in times:
                    if not titem[1]:
                        delta = dt.timedelta(hours=24)
                    
                    else:
                        delta = titem[1]-titem[0]

                    if key == 'explosions':
                        nroday = int(delta.total_seconds()/3600/24)
                        if nroday > 1:
                            for t in [titem[0] + dt.timedelta(hours=x*24) for x in range(nroday+1)]:
                                ha = ax.scatter(t, pos[key][0], s=marker_size[key], marker=marker[key],
                                    alpha=alpha_activity, ec=marker_edge[key], facecolor=colors[key])
                        else:
                            ha = ax.scatter(titem[0], pos[key][0], s=marker_size[key], marker=marker[key],
                                alpha=alpha_activity, ec=marker_edge[key], facecolor=colors[key])

                    else:
                        ha = ax.broken_barh([(titem[0], delta)], pos[key],
                            alpha=alpha_activity, ec=marker_edge[key], facecolor=colors[key])
                    
                    if not handles2 or title[key] not in [h[1] for h in handles2]:
                        handles2 += [(ha, title[key])]
            
        ax.set_xlim(start, end)
        ax.set_ylim(-1.5, 0.65)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
        ax.xaxis.set_minor_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b.%d'))
        ax.yaxis.set_major_locator(mtick.NullLocator())
        ax.grid(axis='x', which='major', color='k', linestyle='--', alpha=0.3)
        ax.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.1)

        # add shadow
        ax.broken_barh([(start, end-start)], (0.15,0.65),
                       alpha=0.2, ec=None, facecolor='gray', zorder=0)

        start += dt.timedelta(days=day_interval_plot)

    fig.legend(handles=[x1[0] for x1 in handles1], labels=[x1[1] for x1 in handles1],
               loc=1, ncol=len(handles1), title="Tremor episodes")
    fig.legend(handles=[x1[0] for x1 in handles2], labels=[x1[1] for x1 in handles2],
               loc=2, ncol=len(handles2), edgecolor='black', facecolor='gray', framealpha=0.2, title="Activity")

    return fig


def get_colors(*args):
    color_dict = {}
    color_dict['zesty'] = ['#F5793A', '#A95AA1', '#85C0F9', '#0F2080']
    color_dict['corp'] = ['#BDB8AD', '#EBE7E0', '#C6D4E1', '#44749D']
    color_dict['elegant'] = ['#ABC3C9', '#E0DCD3', '#CCBE9F', '#382119']
    color_dict['retro'] = ['#601A4A', '#EE442F', '#63ACBE', '#F9F4EC']
    color_dict['okabe'] = ['#009E73', '#56B4E9', '#E69F00', '#000000', '#CC79A7', '#D55E00', '#007282', '#F0E442']
    color_dict['tolb'] = ['#228833', '#EE6677', '#BBBBBB', '#AA3377', '#66CCEE', '#CCBB44', '#4477AA']
    color_dict['toll'] = ['#77AADD', '#AAAA00', '#BBCC33', '#DDDDDD', '#44BB99', '#99DDFF', '#FFAABB', '#EEDD88', '#EE8866']
    color_dict['tolm'] = ['#117733', '#44AA99', '#88CCEE', '#DDDDDD', '#AA4499', '#882255', '#CC6677', '#999933', '#DDCC77', '#332288']
    
    if args[0] in list(color_dict.keys()):
        return color_dict[args[0]]
    
    elif args[0] == 'plot':
        fig = plt.figure(figsize=(13, 8))
        ax1 = fig.add_subplot(241, projection='polar')
        ax2 = fig.add_subplot(242, projection='polar')
        ax3 = fig.add_subplot(243, projection='polar')
        ax4 = fig.add_subplot(244, projection='polar')
        ax5 = fig.add_subplot(245, projection='polar')
        ax6 = fig.add_subplot(246, projection='polar')
        ax7 = fig.add_subplot(247, projection='polar')
        ax8 = fig.add_subplot(248, projection='polar')
        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

        for ax, name in zip(axes, list(color_dict.keys())):
            w = 360/len(color_dict[name])
            pos = np.arange(0, 360, w)*(np.pi/180)
            ax.set_xlabel(name, size=14)
            ax.bar(pos, 1, width=w*(np.pi/180), color=color_dict[name])

            for i in range(len(pos)):
                ax.text(pos[i], 1, i, color='black', fontsize='large')
        
            ax.grid(False)
            ax.spines['polar'].set_visible(False)
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15, wspace=0.5, hspace=0.5)
        plt.show()
    else:
        print('warning: %s is not a colormap' % args[0])
    
    return color_dict


def plot_sde_event(event, **kwargs):
    """
    This function plots the event of SDE database.

    Parameters
    ----------
    event : Event (object) of SDE database
    """

    fq_band = kwargs.get('fq_band', (0.5,10))
    avg_step = kwargs.get('avg_step', None)
    olap = kwargs.get('olap', 0.25)
    add_sta = kwargs.get('stations', None)
    fig = kwargs.get('fig', None)
    return_axes = kwargs.get('return_axes', False)
    remove_response = kwargs.get('remove_response', True)
    sample_rate = kwargs.get('sample_rate', None)
    off_time = kwargs.get('off_time', 0)

    if sample_rate:
        del kwargs['sample_rate']

    del kwargs['remove_response']

    def __check__(sta_id):
        try:
            net_code = sta_id.split('.')[0]
            sta_code = sta_id.split('.')[1]
            loc_code = sta_id.split('.')[2]
            return True
        
        except IndexError:
            return False


    def __getstream__(station_list):
        # create row if station_id not exist but stream does
        
        starttime = event.starttime - dt.timedelta(minutes=5)
        endtime = event.endtime + dt.timedelta(minutes=5)

        stream = None
        for staid in station_list:
            row = event.get_row(staid)

            if row:
                sta = row.get_station()
            
            else:
                new_event = {}
                new_event['event_id'] = event.id
                new_event['event_type'] = 'S'
                new_event['network'] = staid.split('.')[0]
                new_event['station'] = staid.split('.')[1]
                new_event['location'] = staid.split('.')[2]
                n_row = event.sde.append_event(event.id, new_event)
                row = event.sde.get_id(n_row)
                sta = row.get_station()
            
            if not stream:
                stream = sta.get_stream(starttime, endtime, remove_response=remove_response, sample_rate=sample_rate, **kwargs)
            else:
                stream += sta.get_stream(starttime, endtime, remove_response=remove_response, sample_rate=sample_rate, **kwargs)

        return stream

    station_list = event.stations
    if add_sta:
        if isinstance(add_sta, str):
            if add_sta not in station_list and __check__(add_sta):
                station_list.append(add_sta)
        
        elif isinstance(add_sta, list) or isinstance(add_sta, tuple): 
            for s in add_sta:
                if s not in station_list and __check__(s):
                    station_list.append(s)
        else:
            raise ValueError(' add_stations should be list, sring or tuple')
    
    nro_rows = len(station_list)
    stream = __getstream__(station_list)
    
    if not fig:
        fig = plt.figure(figsize=(12,9))

    gs = GridSpec(nro_rows, 2, left=0.05, right=0.95, top=0.92, bottom=0.05, wspace=0.08, hspace=0.2, width_ratios=[1,0.25])
    
    axes = []
    s_axes = []
    for i in range(nro_rows):
        axes.append(fig.add_subplot(gs[i, 0]))
        s_axes.append(fig.add_subplot(gs[i, 1]))

    start = event.starttime - dt.timedelta(seconds=off_time)
    end = event.endtime + dt.timedelta(seconds=off_time)
    comp_c = dict(Z=get_colors('zesty')[0],
                  E=get_colors('zesty')[1],
                  N=get_colors('zesty')[3]
                  )
    
    for i, stid in enumerate(station_list):
        sta = stid.split('.')[1]
        loc = stid.split('.')[2]
        z_comp = stream.get_component('Z', station=sta, loc=loc)
        e_comp = stream.get_component('E', station=sta, loc=loc)
        n_comp = stream.get_component('N', station=sta, loc=loc)

        psdZ, freq = z_comp.psd(start, end,
            fq_band=fq_band,
            avg_step=avg_step,
            olap=olap,
            return_fig=False,
            plot=False
            )
        psdZ = np.array([psdZ])

        psdE, freq = e_comp.psd(start, end,
            fq_band=fq_band,
            avg_step=avg_step,
            olap=olap,
            return_fig=False,
            plot=False
            )
        psdE = np.array([psdE])

        psdN, freq = n_comp.psd(start, end,
            fq_band=fq_band,
            avg_step=avg_step,
            olap=olap,
            return_fig=False,
            plot=False
            )
        psdN = np.array([psdN])

        dataZ = z_comp.get_data(start, end, detrend=True, fq_band=fq_band)
        dataE = e_comp.get_data(start, end, detrend=True, fq_band=fq_band)
        dataN = n_comp.get_data(start, end, detrend=True, fq_band=fq_band)
        time = np.arange(len(dataZ))*z_comp.stats.delta

        maxZ = np.nanmax(dataZ)
        maxE = np.nanmax(dataE)
        maxN = np.nanmax(dataN)

        dataZ /= maxZ
        axes[i].plot(time, dataZ + 5, lw=1.2, color=comp_c['Z'])
        s_axes[i].plot(freq, psdZ.T/psdZ.T.max(), lw=1.2, color=comp_c['Z'], label='%s_Z' % stid)        
        dataE /= maxE
        axes[i].plot(time, dataE + 2.5, lw=1.2, color=comp_c['E'])
        s_axes[i].plot(freq, psdE.T/psdE.T.max(), lw=1.2, color=comp_c['E'], label='%s_E' % stid)    
        dataN /= maxN
        axes[i].plot(time, dataN, lw=1.2, color=comp_c['N'])
        s_axes[i].plot(freq, psdN.T/psdN.T.max(), lw=1.2, color=comp_c['N'], label='%s_N' % stid)   

        axes[i].set_xlim(0, time[-1])
        axes[i].set_ylim(-1.2, 6.8)
        axes[i].set_yticks([5, 2.5, 0])
        axes[i].set_yticklabels(['Z', 'E', 'N'])
        axes[i].xaxis.set_minor_locator(mtick.AutoMinorLocator(2))

        if not i == len(station_list)-1:
            axes[i].xaxis.set_major_formatter(mtick.NullFormatter())
            s_axes[i].xaxis.set_major_formatter(mtick.NullFormatter())
        else:
            fig.suptitle(('%s [EID: %i] %s' % (event.label, event.id, event.starttime.strftime('%Y %b %d %H:%M'))), fontsize=12)
            axes[i].set_xlabel('Time [sec]')
            s_axes[i].set_xlabel('Freq. [Hz]')
            
        axes[i].grid(axis='x', which='major', ls='--', color='k', alpha=0.2)
        axes[i].set_ylabel(stid)
        
        if remove_response:
            axes[i].annotate(f'{maxZ*10**6:.2f}', xy=(0,0.75), color=comp_c['Z'], 
            xycoords='axes fraction', bbox=dict(boxstyle="Round, pad=0.2", fc='w', ec='k'))
            axes[i].annotate(f'{maxE*10**6:.2f}', xy=(0,0.45), color=comp_c['E'], 
            xycoords='axes fraction', bbox=dict(boxstyle="Round, pad=0.2", fc='w', ec='k'))
            axes[i].annotate(f'{maxN*10**6:.2f}', xy=(0,0.12), color=comp_c['N'], 
            xycoords='axes fraction', bbox=dict(boxstyle="Round, pad=0.2", fc='w', ec='k'))
        else:
            axes[i].annotate(f'{maxZ:.1f}', xy=(0,0.75), color=comp_c['Z'], 
            xycoords='axes fraction', bbox=dict(boxstyle="Round, pad=0.2", fc='w', ec='k'))
            axes[i].annotate(f'{maxE:.1f}', xy=(0,0.3), color=comp_c['E'], 
            xycoords='axes fraction', bbox=dict(boxstyle="Round, pad=0.2", fc='w', ec='k'))
            axes[i].annotate(f'{maxN:.1f}', xy=(0,0.12), color=comp_c['N'], 
            xycoords='axes fraction', bbox=dict(boxstyle="Round, pad=0.2", fc='w', ec='k'))
        
        s_axes[i].xaxis.set_major_locator(mtick.MaxNLocator(5))
        s_axes[i].yaxis.set_major_formatter(mtick.NullFormatter())
        s_axes[i].grid(axis='x', which='major', ls='--', color='k', alpha=0.2)
        s_axes[i].set_ylabel('PSD')
    
    if remove_response:
        txt = r'Max. amplitude in $\mu$m/s'
    else:
        txt = f'Max. amplitude in cnts'
    
    axes[0].annotate(txt, xy=(0,1.1), xycoords='axes fraction', color='k')

    if return_axes:
        return axes, s_axes
    
    else:
        return fig


def plot_lte_peaks_evo(peaks, fq_range=None, out='time_evo', format='pdf', spj=5, rj=3, pj=1):
    fig = plt.figure(figsize=(8, 9))
    gs1 = GridSpec(
        6,
        2,
        left=0.1,
        right=0.9,
        top=0.90,
        bottom=0.1,
        wspace=0.05,
        hspace=0.1,
        height_ratios=[0.75]+5*[1],
        width_ratios=[1, 0.02]
        )
    
    # energy
    erg = peaks.get_attr('energy')
    erg[np.where(erg == 0)] = np.nan
    erg = 10*np.log10(erg)

    # x-axis in hours
    bins = erg.shape[0]
    time_step_lte = peaks.stats.time_step
    time_in_hours = np.linspace(0, bins*time_step_lte/60, bins)

    # fq range
    if not fq_range:
        fq_range = peaks.stats.time_step
    
    # plot energy
    erg_ax = fig.add_subplot(gs1[0, 0])
    erg_ax.plot(time_in_hours, erg, 'k')
    erg_ax.set_ylabel('Energy\n' + r'[m$^2\cdot$s$^{-2}$ dB]')
    try:
        erg_ax.set_ylim(
            int(np.percentile(erg, 1)),
            int(np.percentile(erg, 99))
            )
    except:
        pass
    erg_ax.set_xlim(time_in_hours[0], time_in_hours[-1])
    erg_ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(3))
    erg_ax.xaxis.set_major_formatter(mtick.NullFormatter())
    erg_ax.grid(axis='x', which='major', color='k',
                linestyle='-', alpha=0.2)
    erg_ax.grid(axis='x', which='minor', color='k',
                linestyle='--', alpha=0.2)
    
    # get max/min spectral values
    vmin_list = []
    vmax_list = []
    for vlis in peaks.all_peaks:
        if vlis:
            vmin_list += [min(vlis['sp'])]
            vmax_list += [max(vlis['sp'])]
    vmin = np.floor(np.percentile(vmin_list, 10))
    vmax = np.ceil(np.percentile(vmin_list, 90))
    sp_norm = mcolor.Normalize(vmax=vmax, vmin=vmin)
    sp_cmap = 'Spectral_r'
    sp_ball = (abs(np.floor(min(vmin_list)))+10, 1.1)
    sp_bar = cm.ScalarMappable(norm=sp_norm, cmap=sp_cmap)

    sp_ax = fig.add_subplot(gs1[1, 0])
    sp_cax = fig.add_subplot(gs1[1, 1])

    pd_ax = fig.add_subplot(gs1[2, 0])
    pd_cax = fig.add_subplot(gs1[2, 1])

    r_ax = fig.add_subplot(gs1[3, 0])
    r_cax = fig.add_subplot(gs1[3, 1])

    azm_ax = fig.add_subplot(gs1[4, 0])
    azm_cax = fig.add_subplot(gs1[4, 1])

    dip_ax = fig.add_subplot(gs1[5, 0])
    dip_cax = fig.add_subplot(gs1[5, 1])

    axes = [sp_ax, pd_ax, r_ax, azm_ax, dip_ax]

    # colormap for PD
    pd_norm = mcolor.Normalize(vmax=1, vmin=0.8)
    pd_cmap = truncate_cmap('Greys', 0.2, 1)
    pd_bar = cm.ScalarMappable(norm=pd_norm, cmap=pd_cmap)

    # colormap for Rectilinearity
    zesty = get_colors('zesty')
    colors = [zesty[2]]*3 + [zesty[0]]*4 + [zesty[1]]*3
    r_norm = mcolor.Normalize(vmax=1, vmin=0)
    r_cmap = mcolor.LinearSegmentedColormap.from_list("zesty", colors, N=10)
    r_bar = cm.ScalarMappable(norm=r_norm, cmap=r_cmap)

    # colormap for azimuth
    azm_and_dip_cmap = 'Greys'
    azm_norm = mcolor.Normalize(vmax=180, vmin=0)
    azm_bar = cm.ScalarMappable(norm=azm_norm, cmap=azm_and_dip_cmap)

    dip_norm = mcolor.Normalize(vmax=90, vmin=0)
    dip_bar = cm.ScalarMappable(norm=dip_norm, cmap=azm_and_dip_cmap)

    spj = 5
    for n, peak in enumerate(peaks.all_peaks[::spj]):
        t = time_in_hours[::spj][n]
        if peak:
            for k, fp in enumerate(peak['fq']):
                sp = peak['sp'][k]
                pd = peak['pd'][k]
                s = (sp+sp_ball[0])**sp_ball[1]
                
                sp_ax.scatter(t, fp, c=sp, marker='.', norm=sp_norm, s=s, edgecolor='k', lw=0.5, cmap=sp_cmap)
                pd_ax.scatter(t, fp, c=pd, marker='.', norm=pd_norm, s=s, edgecolor='k', lw=0.5, cmap=pd_cmap)

    rj = 3
    for n, peak in enumerate(peaks.all_peaks[::rj]):
        t = time_in_hours[::rj][n]
        if peak:
            for k, fp in enumerate(peak['fq']):
                sp = peak['sp'][k]
                s = (sp+sp_ball[0])**sp_ball[1]
                r = peak['r'][k]
                r_ax.scatter(t, fp, c=r, marker='.', norm=r_norm, s=s, edgecolor='k', lw=0.5, cmap=r_cmap)

    pj = 2
    for n, peak in enumerate(peaks.all_peaks[::pj]):
        t = time_in_hours[::pj][n]
        if peak:
            for k, fp in enumerate(peak['fq']):
                sp = peak['sp'][k]
                s = (sp+sp_ball[0])**sp_ball[1]

                azm = peak['azm'][k]
                azm_ax.scatter(t, fp, c=azm, marker='.', norm=azm_norm, s=s, edgecolor='k', lw=0.5, cmap=azm_and_dip_cmap)

                dip = peak['dip'][k]
                dip_ax.scatter(t, fp, c=dip, marker='.', norm=dip_norm, s=s, edgecolor='k', lw=0.5, cmap=azm_and_dip_cmap)

    # add colorbars
    sp_label = 'PSD\n'+r'[10$\cdot\log_{10} \left(m^2s^{-2}/Hz\right)$ dB]'
    spbar = fig.colorbar(sp_bar, cax=sp_cax, orientation='vertical', label=sp_label)
    spbar.ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=3))

    pdbar = fig.colorbar(pd_bar, cax=pd_cax, orientation='vertical', label='PD')
    pdbar.ax.yaxis.set_major_locator(mtick.FixedLocator([0, 0.8, 1]))

    rbar = fig.colorbar(r_bar, cax=r_cax, orientation='vertical', label='R')
    rbar.ax.yaxis.set_major_locator(mtick.FixedLocator([0.3, 0.7]))

    azm_label = r'$\Theta_H$ [º]'
    azmbar = fig.colorbar(azm_bar, cax=azm_cax, orientation='vertical', label=azm_label)
    azmbar.ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))

    dip_label = r'$\Theta_V$ [º]'
    dipbar = fig.colorbar(dip_bar, cax=dip_cax, orientation='vertical', label=dip_label)
    dipbar.ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))

    # add labels and limits
    for n, ax in enumerate(axes):
        ax.set_ylim(fq_range)
        ax.set_ylabel('Freq. [Hz]')
        
        ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=3))
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))

        ax.set_xlim(time_in_hours[0], time_in_hours[-1])
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        
        ax.grid(axis='x', which='major', color='k', linestyle='-', alpha=0.2)
        ax.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.2)
        
        if n == len(axes)-1:
            ax.set_xlabel('Time [Hours]')
        else:
            ax.xaxis.set_major_formatter(mtick.NullFormatter())

    file_name = out+'.'+format
    if format == 'svg':
        fig.savefig(file_name, transparent=True)
    else:
        fig.savefig(file_name)


def plot_lte_peak_evo(peaks, fq, fq_off=0.1, pd_throld=0.8, r_throld=0.75, out=None, out_dir='./', format='pdf', show=False):
    fig = plt.figure(figsize=(8, 9))
    gs1 = GridSpec(
        4,
        1,
        left=0.1,
        right=0.9,
        top=0.90,
        bottom=0.1,
        wspace=0.05,
        hspace=0.1,
        # width_ratios=[1, 0.02],
        height_ratios=[0.75]+3*[1]
        )
    
    # energy
    erg = peaks.get_attr('energy')
    erg[np.where(erg == 0)] = np.nan
    erg = 10*np.log10(erg)

    # x-axis in hours
    bins = erg.shape[0]
    time_step_lte = peaks.stats.time_step
    time_in_hours = np.linspace(0, bins*time_step_lte/60, bins)
    
    # plot energy
    erg_ax = fig.add_subplot(gs1[0, 0])
    erg_ax.plot(time_in_hours, erg, 'k')
    erg_ax.set_ylabel('Energy\n' + r'[m$^2\cdot$s$^{-2}$ dB]')

    try:
        erg_ax.set_ylim(
            int(np.percentile(erg, 1)),
            int(np.percentile(erg, 99))
            )
    except:
        pass
    
    sp_ax = fig.add_subplot(gs1[1, 0])
    # r_ax = fig.add_subplot(gs1[2, 0])
    azm_ax = fig.add_subplot(gs1[2, 0])
    dip_ax = fig.add_subplot(gs1[3, 0])
    axes = [erg_ax, sp_ax, azm_ax, dip_ax]

    def add_points(ax, y, time_bin):
        if len(y) > 1:
            # ax.errorbar(time_bin, np.mean(y), ls='', yerr=[y.min(), y.max()], marker='s', color='k', alpha=0.7)
            ax.scatter(time_bin, np.mean(y), marker='o', color='k', alpha=0.7)
        else:
            ax.scatter(time_bin, y[0], marker='o', color='k', alpha=0.7)

    # get the dataset
    for n, peak in enumerate(peaks.all_peaks):
        t = time_in_hours[n]
        
        if peak:
            sp_t = np.array([peak['sp'][k] for k, fp in enumerate(peak['fq']) if fq-fq_off<fp<fq+fq_off])
            pd_t = np.array([peak['pd'][k] for k, fp in enumerate(peak['fq']) if fq-fq_off<fp<fq+fq_off])
            r_t = np.array([peak['r'][k] for k, fp in enumerate(peak['fq']) if fq-fq_off<fp<fq+fq_off])
            tH_t = np.array([peak['azm'][k] for k, fp in enumerate(peak['fq']) if fq-fq_off<fp<fq+fq_off])
            tV_t = np.array([peak['dip'][k] for k, fp in enumerate(peak['fq']) if fq-fq_off<fp<fq+fq_off])
            
            if sp_t.size > 0:
                add_points(sp_ax, sp_t, t)

                if any(pd_t[pd_t > pd_throld]):
                    # add_points(pd_ax, pd_t[pd_t > pd_throld], t)
                    r_pdth = r_t[np.where(pd_t > pd_throld)]
                    tH_pdth = tH_t[np.where(pd_t > pd_throld)]
                    tV_pdth = tV_t[np.where(pd_t > pd_throld)]
                    # add_points(r_ax, r_pdth, t)

                    if any(r_pdth[r_pdth > r_throld]):
                        tH_rth = tH_pdth[np.where(r_pdth > r_throld)]
                        tV_rth = tV_pdth[np.where(r_pdth > r_throld)]
                        add_points(azm_ax, tH_rth, t)
                        add_points(dip_ax, tV_rth, t)
    
    sp_label = 'PSD\n'+r'[10$\cdot\log_{10} \left(m^2s^{-2}/Hz\right)$ dB]'
    tH_label = r'$\Theta_H$ [º]'
    tV_label = r'$\Theta_V$ [º]'

    sp_ax.set_ylabel(sp_label)
    # r_ax.set_ylabel(r_label)
    azm_ax.set_ylabel(tH_label)
    dip_ax.set_ylabel(tV_label)
    [ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4)) for ax in axes[1:]]
    # r_ax.set_ylim(0, 1)
    azm_ax.set_ylim(0, 180)
    dip_ax.set_ylim(0, 90)
    
    [ax.set_xlim(time_in_hours[0], time_in_hours[-1]) for ax in axes]
    [ax.xaxis.set_major_formatter(mtick.NullFormatter()) for ax in axes[:-1]]
    [ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(3)) for ax in axes]
    [ax.grid(axis='x', which='major', color='k', linestyle='-', alpha=0.2) for ax in axes]
    [ax.grid(axis='x', which='minor', color='k', linestyle='--', alpha=0.2) for ax in axes]
    axes[0].set_title(f' Freq.: {fq:.2f} [Hz]')
    axes[-1].set_xlabel('Time [hr]')

    if show:
        plt.show()
    
    else:
        if not out:
            pk_l = [(n, peaks.peaks[n]['fq']) for n in range(1, peaks.nro+1)]
            idx = list(filter(lambda x: x[1] == fq, pk_l))[0][0]
            out = 'tevo_' + peaks.stats.id + '_' + str(idx)

        file_name = os.path.join(out_dir, out + '.' + format)
        
        if format == 'svg':
            fig.savefig(file_name, transparent=True)
        
        else:
            fig.savefig(file_name)


def plot_out_ccorr(out, maxSMB=0.5, save=False):

    from seisvo import Array

    def azimuthFmtt(x, pos):
            azimuths = range(out['model'].get('src_dgr')[0], out['model'].get('src_dgr')[1])
            try:
                return str(azimuths[int(x)])
            except:
                pass

    ans = Array.get_max_values(out)
    time_bin = out['time']
    time_bin_max = len(time_bin)

    fig = plt.figure(figsize=(8, 6))
    fig.suptitle('%s   --   %s' % (out['time'][0], out['time'][-1]))

    ax1 = fig.add_axes([0.1, 0.75, 0.8, 0.15])
    ax2 = fig.add_axes([0.1, 0.55, 0.8, 0.15])
    ax3 = fig.add_axes([0.1, 0.35, 0.8, 0.15])
    ax4 = fig.add_axes([0.1, 0.15, 0.8, 0.15])
    cbar_ax = fig.add_axes([0.91, 0.15, 0.025, 0.15])

    # pressure plot
    p_max_r = ans[2][np.where(ans[0] >= maxSMB)]
    t_az = time_bin[np.where(ans[0] >= maxSMB)]
    ax1.scatter(t_az, p_max_r, color='darkred', label=r'$P_{max}$')
    ax1.plot(time_bin, ans[3], color='k', label=r'$P_{avg}$')
    ax1.legend(bbox_to_anchor=(1.1, 1.05))
    ax1.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
    ax1.set_ylabel(r'P [Pa]')
    ax1.set_xlim(time_bin[0], time_bin[-1])
    ax1.grid(True)

    # corr plot
    ax2.plot(time_bin, ans[0], color='darkgreen')
    ax2.axhline(y=maxSMB, color='red', linestyle='--')
    ax2.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
    ax2.set_ylabel(r'$r_{max}$')
    ax2.set_xlim(time_bin[0], time_bin[-1])
    ax2.xaxis.set_major_formatter(mtick.NullFormatter())
    ax2.grid(True)

    # azimuth plot
    az_r = ans[1][np.where(ans[0]>=maxSMB)]
    ax3.scatter(t_az, az_r)
    ax3.set_ylabel(r'$azm_{r_{max}}$ [º]')
    ax3.set_xlim(time_bin[0], time_bin[-1])
    ax3.set_ylim(0, 359)
    ax3.yaxis.set_major_locator(mtick.FixedLocator([0,90,180,270]))
    ax3.xaxis.set_major_formatter(mtick.NullFormatter())
    ax3.grid(True)

    # crosscorrelation plot
    mcorr = out['mcorr']
    mcorr = np.flipud(mcorr.T)
    azmbin = range(mcorr.shape[0])
    halfbin_time = (time_bin[1] - time_bin[0]) / 2.0
    halfbin_azmbin = (azmbin[1] - azmbin[0]) / 2.0
    extent = (mdates.date2num(time_bin[0] - halfbin_time),
                mdates.date2num(time_bin[-1] + halfbin_time),
                azmbin[0] - halfbin_azmbin,
                azmbin[-1] + halfbin_azmbin)

    cmap = plt.get_cmap('Spectral_r')
    norm = mcolor.Normalize(0, 1)

    im = ax4.imshow(mcorr, cmap=cmap, norm=norm, interpolation='gaussian', extent=extent, aspect='auto')
    ax4.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
    ax4.yaxis.set_major_formatter(mtick.FuncFormatter(azimuthFmtt))

    fig.colorbar(im, cax=cbar_ax, ticks=[0, 0.5, 1], orientation='vertical',
                            format='%.1f')
    lims = ax4.axis('tight')
    ax4.xaxis_date()
    ax1.xaxis.set_major_formatter(mtick.NullFormatter())
    ax4.set_xlabel('Time')
    ax4.set_ylabel('azm [º]')

    ax4.annotate(f"window_length: {out['info']['time_width']}\noverlap: {out['info']['overlap']}\nmaxSMB: {maxSMB}", xy=(0.05,0.05), xycoords='figure fraction', fontsize=8, fontfamily='monospace')
    
    if save:
        date_fmt = out['time'][0].strftime('%Y%m%d')
        name = '%s.%s.png' % (out['info']['code'], date_fmt)
        fig.savefig(name, dpi=200)
        del fig
    
    else:
        plt.show()


def plot_peaks_spec_pdf(peak, plot=True, **kwargs):
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

    if plot:
        plt.show()

    return fig


def plot_peaks_pdf(peak, n, plot=True, **kwargs):
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
    
    if plot:
        plt.show()

    return fig