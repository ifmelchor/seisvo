#!/usr/bin/env python3
# coding=utf-8

import scipy
import numpy as np
import datetime as dt

from matplotlib import cm
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.ticker as mtick
import matplotlib.dates as mdates

from .utils import plot_gram, get_colors, get_time_format

default_LTE_color = get_colors('zesty')[1]

default_labels = {
    'energy':'Energy',
    'specgram':'PSD',
    'degree':'PD',
    'rect':'Rect.',
    'azimuth':'Azimuth',
    'elevation':'Elevation',
    'pentropy': 'PE',
    'dsar': 'DSAR'
}

def count_rows(lteout, chan_list, attr_list):
        nrow = 0
        for attr in attr_list:
            if attr in lteout.chan_attr_list:
                for chan in chan_list:
                    nrow += 1
            else:
                nrow += 1
        return nrow


def get_vector_bar(attr, **kwargs):
    if attr == 'specgram':
        v_min = kwargs.get('specgram_vmin', None)
        v_max = kwargs.get('specgram_vmax', None)
        cmap  = kwargs.get('specgram_cmap', 'Spectral_r')

    if attr == 'degree':
        v_min = 0
        v_max = 1
        cmap  = kwargs.get('degree_cmap', 'Spectral_r')
    
    if attr == 'rect':
        v_min = 0
        v_max = 1
        cmap  = kwargs.get('rect_cmap', 'Spectral_r')

    if attr == 'azimuth':
        v_min = 0
        v_max = 180
        cmap  = kwargs.get('azimuth_cmap', 'twilight')
    
    if attr == 'elev':
        v_min = 0
        v_max = 90
        cmap  = kwargs.get('elev_cmap', 'Spectral_r')
    
    return cmap, (v_min, v_max)


def get_axes_dict(lteout, attr_list, chan_list, fig):

    grid = {'hspace':0.15, 'left':0.08, 'right':0.92, 'wspace':0.05, 'top':0.95, 'bottom':0.05}
    
    if lteout.any_vector(attr_list):
        col = 2
        grid['width_ratios'] = [1, 0.01]
    else:
        col = 1
    
    row = count_rows(lteout, chan_list, attr_list)

    axes = fig.subplots(row, col, gridspec_kw=grid)
    
    if isinstance(axes, np.ndarray):
        axes = axes.reshape(row, col)
    else:
        axes = np.array([axes]).reshape(row, col)

    axes_dict = {}
    n = 0
    for attr in attr_list:
        if attr in lteout.chan_attr_list:
            for chan in chan_list:
                key = '/'.join([chan, attr])

                if lteout.any_vector([attr]):
                    axes_dict[key] = (axes[n,0], axes[n,1])
                
                else:
                    axes[n,1].axes.get_xaxis().set_visible(False)
                    axes[n,1].axes.get_yaxis().set_visible(False)
                    axes[n,1].set_frame_on(False)
                    axes_dict[key] = (axes[n,0],)
                
                n += 1
        else:
            key = attr
            if lteout.any_vector([key]):
                axes_dict[key] = (axes[n,0], axes[n,1])
            
            else:
                axes[n,1].axes.get_xaxis().set_visible(False)
                axes[n,1].axes.get_yaxis().set_visible(False)
                axes[n,1].set_frame_on(False)
                axes_dict[key] = (axes[n,0],)
            n += 1
    
    return axes_dict


def LTESTAplot(lteout, chan_list, attr_list, fig=None, return_stats=False, plot=False, **kwargs):

    if not fig:
        fig = plt.figure()
    
    axes = get_axes_dict(lteout, attr_list, chan_list, fig)
            
    # define time and x format
    hr_interval = (lteout.endtime - lteout.starttime).total_seconds()/3600
    time_format = get_time_format(True, hr_interval/24)
    time = lteout._dout["dtime"]
    freq = lteout._dout["freq"]
    
    # get stats of scalar values
    nrow = len(axes)
    n = 0
    stats = {}
    for attr in attr_list:
        if attr in lteout.chan_attr_list:
            for chan in chan_list:
                key = '/'.join([chan, attr])
                data = lteout._dout[chan][attr]
                ax = axes[key]
                
                if not lteout.any_vector([attr]): # scalar attribute
                    ax[0].plot(time, data, color='k', ls="--")
                    ax[0].scatter(time, data, color='blue', ec="k")
                    data_stats = lteout.get_stats(attr, chan)
                    
                    # set vlimits
                    vmin, vmax = (None, None)
                    if kwargs.get(attr+"_vmin", None): 
                        vmin = kwargs.get(attr+"_vmin")
                    
                    if kwargs.get(attr+"_vmax", None):
                        vmax = kwargs.get(attr+"_vmax")
                    
                    if not vmin or not vmax:
                        stats[key] = data_stats
                        if not vmin:
                            vmin = np.floor(data_stats[0])
                        else:
                            vmax = np.ceil(data_stats[1])
                    
                    ax[0].set_ylim(vmin, vmax)
                    ax[0].set_xlim(time[0], time[-1])
                    ax[0].set_ylabel(default_labels.get(attr, attr))
                
                else:
                    cmap, vlim = get_vector_bar(attr, **kwargs)
                    plot_gram(freq, data.T, time, ax[0], axis_bar=ax[1],\
                        v_max=vlim[1], v_min=vlim[0], cmap=cmap,\
                        bar_label=default_labels.get(attr, attr))

                    ax[0].set_ylim(freq[0], freq[-1])
                    ax[0].set_ylabel("Freq")
                
                ax[0].xaxis.set_major_locator(time_format[0][0])
                ax[0].xaxis.set_minor_locator(time_format[1][0])
                ax[0].xaxis.set_minor_formatter(mtick.NullFormatter())
                ax[0].xaxis.set_major_formatter(mtick.NullFormatter())
                ax[0].yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
                ax[0].yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
                ax[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
                
                n += 1
                if n == nrow:
                    # last row
                    ax[0].xaxis.set_major_formatter(time_format[0][1])
                    ax[0].set_xlabel('Time [h]')
                    
                    if time_format[1][1]:
                        ax[0].xaxis.set_minor_formatter(time_format[1][1])

        else:
            data = lteout._dout[attr]
            ax = axes[attr]
            cmap, vlim = get_vector_bar(attr, **kwargs)
            plot_gram(freq, data.T, time, ax[0], axis_bar=ax[1], v_max=vlim[1],\
                v_min=vlim[0], cmap=cmap, bar_label=default_labels.get(attr, attr))
            ax[0].set_ylim(freq[0], freq[-1])

            ax[0].xaxis.set_major_locator(time_format[0][0])
            ax[0].xaxis.set_minor_locator(time_format[1][0])
            ax[0].xaxis.set_minor_formatter(mtick.NullFormatter())
            ax[0].xaxis.set_major_formatter(mtick.NullFormatter())
            ax[0].yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
            ax[0].yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
            ax[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))

            n += 1
            if n == nrow:
                # last row
                ax[0].xaxis.set_major_formatter(time_format[0][1])
                ax[0].set_xlabel('Time [h]')
                
                if time_format[1][1]:
                    ax[0].xaxis.set_minor_formatter(time_format[1][1])
    
    if plot:
        plt.show()
    
    if return_stats:
        return stats
    else:
        return (fig, axes)



def plotPeaksSpecPDF(peak_dict, fq, fq_th=None, plot=True, **kwargs):

    grid = {'left':0.15, 'right':0.95, 'top':0.90, 'bottom':0.2}
    fig, ax = plt.subplots(1,1, figsize=kwargs.get('figsize',(9,3)), gridspec_kw=grid)

    norm = fq[1].max()
    ax.plot(fq[0], fq[1]/norm, color='k') # plot normalized PDF

    if fq_th:
        ax.axhline(y=fq_th, color='r', ls='--')
    
    for _, info in peak_dict.items():
        ax.scatter(info['fq'], info['fq_prob']/norm, color='r', marker='o', ec='k')
        half_width = info['width']
        y1_index = np.argmin(np.abs(fq[0]-(info['fq']-half_width)))
        y2_index = np.argmin(np.abs(fq[0]-(info['fq']+half_width)))
        x_fill = fq[0][y1_index:y2_index].reshape(-1,)
        y_fill = fq[1][y1_index:y2_index]/norm
        ax.fill_between(x_fill, y_fill, color='k', alpha=0.1)
        # ax.annotate(f'#{p_i}', xy=(info['fq']+0.05, info['fq_prob']+0.01), bbox=dict(boxstyle="round", fc="w", ec="k", lw=0.8), fontsize=8, family='monospace')
    
    ax.set_xlim(fq[0].min(), fq[0].max())
    ax.set_ylim(0, 1.1)
    ax.set_ylabel(r'$P(\mathcal{P}_{wp}$)')
    ax.set_xlabel(r'$f$ [Hz]')

    ax.set_xticks([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])
    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(mtick.FixedLocator([0, 0.5, 1.0]))

    ax.grid(which='major', axis='both', ls='--', alpha=0.5)
    ax.grid(which='minor', axis='both', ls=':', alpha=0.25)

    ax.tick_params(axis='y')
    ax.tick_params(axis='x')

    if plot:
        plt.show()

    return fig


def plotPeaksPDF(nPeak, plot=True):

    PDFs = {}
    if nPeak['sp']:
        sp_space = np.linspace(nPeak['sp']['min'], nPeak['sp']['max'], 1000).reshape(-1,1)
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

    if not isinstance(axis, np.ndarray):
        axis = [axis]

    for ax, att in zip(axis, PDFs.keys()):
        space, pdf = PDFs[att]
        pdf /= pdf.max()
        ax.plot(space, pdf, color='k')

        if att == 'rect':
            # rect_th = peak.peak_thresholds['rect_th']
            # ax.fill_between(space[space>rect_th], pdf[np.where(space>rect_th)[0]], alpha=0.1, color='k')
            ax.set_xlim(0,1)
            ax.set_xlabel(r'R')
            ax.set_xticks([0, 0.5, 1])
            ax.set_ylabel(r'$P(\mathcal{R}$)')
            ax.set_title(r'$N_R$ = %s' % nPeak['rect']['n'])
        
        elif att == 'thH':
            ax.set_xlim(0,180)
            ax.set_xlabel(r'$\Theta_H$ [$\degree$]')
            ax.set_xticks([0, 90, 180])
            ax.set_ylabel(r'$P(\mathcal{H}$)')
            ax.set_title(r'$N_H$ = %s' % nPeak['thH']['n'])

        elif att == 'thV':
            ax.set_xlim(0,90)
            ax.set_ylabel(r'$P(\mathcal{V}$)')
            ax.set_xlabel(r'$\Theta_V$ [$\degree$]')
            ax.set_xticks([0, 45, 90])
            ax.set_title(r'$N_V$ = %s' % nPeak['thV']['n'])
        
        else:
            ax.set_xlim(space[0],space[-1])
            ax.set_xlabel('PSD\n' +r'[dB//(m$^2$s$^{-2}$/Hz)]')
            ax.set_ylabel(r'$P(\mathcal{S}$)')
            ax.set_title(r'$N_S$ = %s' % nPeak['sp']['n'])
        
        ax.set_ylim(0,1.1)
        ax.tick_params(axis='x')
        ax.set_yticks([0, 0.5, 1])
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.grid(which='major', axis='both', ls='--', alpha=0.5)
        ax.grid(which='minor', axis='both', ls=':', alpha=0.25)
        
        if ax == axis[0]:
            ax.tick_params(axis='y')
        else:
            ax.yaxis.set_major_formatter(mtick.NullFormatter())
    
    if plot:
        plt.show()

    return fig


def plotDFreqTimeEvo(dfq, fhigh=None, flow=None, plot=True, **kwargs):

    fig = plt.figure(figsize=(13,8))
    gs = GridSpec(4, 2, figure=fig, hspace=0.1, left=0.08, right=0.92, wspace=0.05, top=0.95, bottom=0.05, width_ratios=[1, 0.02])

    sxx_ax = fig.add_subplot(gs[0, 0])
    rec_ax = fig.add_subplot(gs[1, 0])
    azi_ax = fig.add_subplot(gs[2, 0])
    ele_ax = fig.add_subplot(gs[3, 0])
    fq_cax = fig.add_subplot(gs[:, 1])

    if not fhigh:
        fhigh = np.ceil(np.nanmax(dfq["fq"]))

    if not flow:
        flow = np.floor(np.nanmin(dfq["fq"]))

    sxx_ax.set_title(kwargs.get("title", None))
    cmap = plt.get_cmap('Spectral_r')
    norm = mcolor.Normalize(flow, fhigh)

    sxx_ax.errorbar(dfq["time"], dfq["sxx"], yerr=np.vstack((dfq["sxx_r1"], dfq["sxx_r2"])), fmt="k.", capsize=5, zorder=1)
    sxx_ax.scatter(dfq["time"], dfq["sxx"], c=dfq["fq"], ec="k", norm=norm, cmap=cmap, zorder=2)
    sxx_min = np.floor(np.min(dfq["sxx"]-dfq["sxx_r1"]))
    sxx_max = np.ceil(np.max(dfq["sxx"]+dfq["sxx_r2"]))
    sxx_ax.set_ylim(sxx_min,sxx_max)
    sxx_ax.set_ylabel("Power [dB]")

    rec_ax.errorbar(dfq["time"], dfq["rect"], yerr=np.vstack((dfq["rect_r1"], dfq["rect_r2"])), fmt="k.", capsize=5, zorder=1)
    rec_ax.scatter(dfq["time"], dfq["rect"], c=dfq["fq"], ec="k", norm=norm, cmap=cmap, zorder=2)
    rec_ax.set_ylim(-0.1,1.1)
    rec_ax.set_ylabel("Rect.")

    azi_ax.errorbar(dfq["time"], dfq["azim"], yerr=np.vstack((dfq["azim_r1"], dfq["azim_r2"])), fmt="k.", capsize=5, zorder=1)
    azi_ax.scatter(dfq["time"], dfq["azim"], c=dfq["fq"], ec="k", norm=norm, cmap=cmap, zorder=2)
    azi_ax.set_ylim(-5,185)
    azi_ax.set_ylabel("Azimuth")

    ele_ax.errorbar(dfq["time"], dfq["elev"], yerr=np.vstack((dfq["elev_r1"], dfq["elev_r2"])), fmt="k.", capsize=5, zorder=1)
    ele_ax.scatter(dfq["time"], dfq["elev"], c=dfq["fq"], ec="k", norm=norm, cmap=cmap, zorder=2)
    ele_ax.set_ylim(-5,95)
    ele_ax.set_ylabel("Elevation")

    minor_locator = mdates.DayLocator(interval=1)
    major_locator = mdates.MonthLocator(interval=1)
    major_formatt = mdates.DateFormatter('%d %b %y')
    for ax in [sxx_ax, rec_ax, azi_ax, ele_ax]:
        ax.set_xlim(dfq["starttime"], dfq["endtime"])
        ax.xaxis.set_minor_locator(minor_locator)
        ax.xaxis.set_major_locator(major_locator)
        ax.grid(which="major", ls="--", alpha=0.3, color="k")
        ax.grid(which="minor", ls=":", alpha=0.2, color="k")
        if ax != ele_ax:
            ax.xaxis.set_major_formatter(mtick.NullFormatter())
        else:
            ax.xaxis.set_major_formatter(major_formatt)

    fq_im = cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(fq_im, cax=fq_cax, orientation='vertical', label="Hz")

    if plot:
        plt.show()
    
    return fig


def plotPeakTimeEvo(pkt, plot=True, **kwargs):

    title = kwargs.get("title", None)
    flow  = kwargs.get("flow", None)
    fhigh = kwargs.get("fhigh", None)
    t_min = kwargs.get("tmin", None)
    t_max = kwargs.get("tmax", None)
    ax_dict = kwargs.get("ax_dict", {})
    date_format = kwargs.get("date_format", mdates.DateFormatter('%d %b\n%H:%M'))
    
    if not ax_dict:
        fig = plt.figure(figsize=(13,8))
        gs = GridSpec(4, 3, figure=fig, hspace=0.1, left=0.08, right=0.92, wspace=0.05, top=0.95, bottom=0.05, width_ratios=[1, 0.2, 0.02])
        sxx_ax = fig.add_subplot(gs[0, 0])
        sxx_pdf_ax = fig.add_subplot(gs[0, 1])
        rec_ax = fig.add_subplot(gs[1, 0])
        rec_pdf_ax = fig.add_subplot(gs[1, 1])
        azi_ax = fig.add_subplot(gs[2, 0])
        azi_pdf_ax = fig.add_subplot(gs[2, 1])
        ele_ax = fig.add_subplot(gs[3, 0])
        ele_pdf_ax = fig.add_subplot(gs[3, 1])
        fq_cax = fig.add_subplot(gs[:, 2])
        ax_dict = {
            "sxx":sxx_ax,
            "sxx_pdf":sxx_pdf_ax,
            "rec":rec_ax,
            "rec_pdf":rec_pdf_ax,
            "azi":azi_ax,
            "azi_pdf":azi_pdf_ax,
            "ele":ele_ax,
            "ele_pdf":ele_pdf_ax,
            "fq_cax":fq_cax
        }
    else:
        plot = False
        fig = None

    # sxx norm/cmap
    fq_min = np.floor(np.min([pkt[ch]["stats"]["fq"][0] for ch in list(pkt.keys())]))
    fq_max = np.ceil(np.max([pkt[ch]["stats"]["fq"][1] for ch in list(pkt.keys())]))
    cmap = plt.get_cmap('Spectral_r')
    norm = mcolor.Normalize(fq_min, fq_max)

    chan = list(pkt.keys())

    ttt = np.array([])
    fqt = np.array([])
    sxt = np.array([])
    rct = np.array([])
    azt = np.array([])
    evt = np.array([])
    
    for ch in chan:
        nbin = pkt[ch]["time"].shape[0]
        for n in range(nbin):
            fqn = pkt[ch]["pks"]["fq"][n]
            sxn = pkt[ch]["pks"]["sxx"][n]
            rcn = pkt[ch]["pks"]["rect"][n]
            azn = pkt[ch]["pks"]["azim"][n]
            evn = pkt[ch]["pks"]["elev"][n]

            # filter in frequency
            if flow and fhigh:
                fql = np.where((fqn >= flow)&(fqn <= fhigh))
            elif flow and not fhigh:
                fql = np.where(fqn >= flow)
            elif fhigh and not flow:
                fql = np.where(fqn <= fhigh)
            else:
                fql = None

            if fql:
                fqn = fqn[fql]
                sxn = sxn[fql]
                rcn = rcn[fql]
                azn = azn[fql]
                evn = evn[fql]

            tnn = np.array([pkt[ch]["time"][n]]*len(fqn))
            ttt = np.hstack((ttt,tnn))
            fqt = np.hstack((fqt,fqn))
            sxt = np.hstack((sxt,sxn))
            rct = np.hstack((rct,rcn))
            azt = np.hstack((azt,azn))
            evt = np.hstack((evt,evn))

    ax_list = []
    
    if ax_dict["sxx"]:
        ax_dict["sxx"].scatter(ttt, sxt, c=fqt, ec="k", alpha=0.5, norm=norm, cmap=cmap)
        sxx_min = np.floor(np.min([pkt[ch]["stats"]["sxx"][0] for ch in list(pkt.keys())]))
        sxx_max = np.ceil(np.max([pkt[ch]["stats"]["sxx"][1] for ch in list(pkt.keys())]))
        ax_dict["sxx"].set_ylim(sxx_min, sxx_max)
        ax_dict["sxx"].set_ylabel("Power [dB]")
        ax_list.append(ax_dict["sxx"])
    
    if ax_dict["rec"]:
        ax_dict["rec"].scatter(ttt, rct, c=fqt, ec="k", alpha=0.5, norm=norm, cmap=cmap)
        ax_dict["rec"].set_ylim(-0.1, 1.1)
        ax_dict["rec"].set_ylabel("Rect.")
        ax_list.append(ax_dict["rec"])
    
    if ax_dict["azi"]:
        ax_dict["azi"].scatter(ttt, azt, c=fqt, ec="k", alpha=0.5, norm=norm, cmap=cmap)
        ax_dict["azi"].set_ylim(-5, 185)
        ax_dict["azi"].set_ylabel("Azimuth")
        ax_list.append(ax_dict["azi"])
    
    if ax_dict["ele"]:
        ax_dict["ele"].scatter(ttt, evt, c=fqt, ec="k", alpha=0.5, norm=norm, cmap=cmap)
        ax_dict["ele"].set_ylim(-2, 92)
        ax_dict["ele"].set_ylabel("Elevation")
        ax_list.append(ax_dict["ele"])

    ax_list[-1].xaxis.set_major_formatter(date_format)
    
    # time norm
    if not t_min:
        t_min = min([min(pkt[ch]["time"]) for ch in list(pkt.keys())])
    
    if not t_max:
        t_max = max([max(pkt[ch]["time"]) for ch in list(pkt.keys())])
    
    for ax in ax_list:
        ax.set_xlim(t_min, t_max)
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(4))
        ax.grid(which="major", ls="--", alpha=0.3, color="k")
        ax.grid(which="minor", ls=":", alpha=0.2, color="k")
        if ax != ax_list[-1]:
            ax.xaxis.set_major_formatter(mtick.NullFormatter())

    # add colorbar
    if ax_dict["fq_cax"]:
        fq_im = cm.ScalarMappable(norm=norm, cmap=cmap)
        plt.colorbar(fq_im, cax=ax_dict["fq_cax"], orientation='vertical', label="Hz")

    # compute probability
    ax_pdf_list = []
    if ax_dict["sxx_pdf"]:
        sxx_space = np.linspace(sxx_min, sxx_max, 500)
        sxx_kde = scipy.stats.gaussian_kde(sxt)
        sxx_pdf = sxx_kde(sxx_space)
        ax_dict["sxx_pdf"].plot(sxx_pdf, sxx_space, color='k')
        ax_dict["sxx_pdf"].set_ylim(sxx_min, sxx_max)
        ax_pdf_list.append(ax_dict["sxx_pdf"])
        # sxx_pdf_ax.set_title("PDF")
    
    if ax_dict["rec_pdf"]:
        ret_space = np.linspace(0, 1, 500)
        ret_kde = scipy.stats.gaussian_kde(rct[np.isfinite(rct)])
        ret_pdf = ret_kde(ret_space)
        ax_dict["rec_pdf"].plot(ret_pdf, ret_space, color='k')
        ax_dict["rec_pdf"].set_ylim(-0.1, 1.1)
        ax_pdf_list.append(ax_dict["rec_pdf"])
    
    if ax_dict["azi_pdf"]:
        azi_space = np.linspace(0, 180, 500)
        azi_kde = scipy.stats.gaussian_kde(azt[np.isfinite(azt)])
        azi_pdf = azi_kde(azi_space)
        ax_dict["azi_pdf"].plot(azi_pdf, azi_space, color='k')
        ax_dict["azi_pdf"].set_ylim(-5, 185)
        ax_pdf_list.append(ax_dict["azi_pdf"])
    
    if ax_dict["ele_pdf"]:
        ele_space = np.linspace(0, 90, 500)
        ele_kde = scipy.stats.gaussian_kde(evt[np.isfinite(evt)])
        ele_pdf = ele_kde(ele_space)
        ax_dict["ele_pdf"].plot(ele_pdf, ele_space, color='k')
        ax_dict["ele_pdf"].set_ylim(-2, 92)
        ax_pdf_list.append(ax_dict["ele_pdf"])

    for ax in ax_pdf_list:
        ax.xaxis.set_major_formatter(mtick.NullFormatter())
        ax.yaxis.set_major_formatter(mtick.NullFormatter())
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.grid(which="major", ls="--", alpha=0.3, color="k")
        ax.grid(which="minor", ls=":", alpha=0.2, color="k")

    if plot:
        plt.show()
    
    return fig


