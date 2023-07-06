
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
from .utils import get_colors

def location_map(arr, exclude_locs=[], show=True):

    fig, ax = plt.subplots(1,1, figsize=(4,4))

    for loc, utm in arr.utm.items():
        if loc not in exclude_locs:
            x, y = utm["easting"], utm["northing"]
            ax.scatter(x, y, marker="^", color="k")
            ax.annotate(loc, xy=(x, y))
    
    ax.set_ylabel("UTM northing")
    ax.set_xlabel("UTM easting")

    if show:
        plt.show()

    return fig


def traces_psd(psd_dict, freq, db_scale=True, vmin=None, vmax=None, show=False, title=None, colorname="tolm"):

    fig, ax = plt.subplots(1,1, figsize=(10,8))
    colorlist = get_colors(colorname)

    ymin = 999
    ymax = 0

    for n, (loc, psd) in enumerate(psd_dict.items()):

        if db_scale:
            psd = 10*np.log10(psd)

        if ymin > psd.min():
            ymin = psd.min()

        if ymax < psd.max():
            ymax = psd.max()

        if bool(n%2):
            ls = "--"
        else:
            ls = "-" 

        ax.plot(freq, psd, color=colorlist[n], label=loc, ls=ls)
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.grid(which="major", axis="x", color="k", ls="-", alpha=0.25)
        ax.grid(which="minor", axis="x", color="k", ls="--", alpha=0.15)
        ax.set_ylabel("PSD [dB]")
        ax.set_xlabel("Freq [Hz]")
        ax.set_title(title)

    if not vmin:
        vmin = ymin

    if not vmax:
        vmax = ymax

    ax.set_ylim(vmin,vmax)
    fig.legend(title="LOC", ncol=1)

    if show:
        plt.show()

    return fig
    

def simple_cc8_plot(dtime, datattr, bounds, slowpdf, bazmpdf, show=True, **kwargs):
    fig          = kwargs.get("fig", None)
    title        = kwargs.get("title", None)
    maac_rv      = kwargs.get("maac_rv", None)
    x_time       = kwargs.get("x_time", False)
    rms_max      = kwargs.get("rms_max", False)
    rms_min      = kwargs.get("rms_min", False)
    return_fdict = kwargs.get("return_fdict", False)

    grid = {'hspace':0.15, 'left':0.08, 'right':0.92, 'wspace':0.05, 'width_ratios':[1,0.1]}

    if fig:
        axes = fig.subplots(4, 2, gridspec_kw=grid)
        show = False
    else:
        fig, axes = plt.subplots(4, 2, figsize=(12,8), sharex='col', gridspec_kw=grid)

    axes[0,0].set_title(title)

    if not x_time:
        duration = (dtime[-1]-dtime[0]).total_seconds()/60
        npts     = len(dtime)
        time     = np.linspace(0, duration, npts)
    else:
        time     = dtime

    default_labels = {
    'rms':'RMS [dB]',
    'slow':'Slowness [s/km]',
    'bazm':'Back-azimuth' r'[$\degree$]',
    'maac':'MAAC',
    }

    fig_dict = {}

    for n, attr in enumerate(["maac", "rms", "slow", "bazm"]):
        fig_dict[attr] = {}
        fig_dict[attr]["axis"] = axes[n,0]
        colors   = ["blue"]*len(datattr[attr])
        fig_dict[attr]["sc"]   = axes[n,0].scatter(time, datattr[attr], facecolor=colors, edgecolor="k", alpha=0.7, zorder=3)

        axes[n,0].set_ylabel(default_labels[attr])

        if not x_time:
            axes[n,0].set_xlim(0, duration)
        else:
            axes[n,0].set_xlim(dtime[0],dtime[-1])
        
        if attr == "rms":
            axes[n,0].set_ylim(np.floor(rms_min), np.ceil(rms_max))

        if attr == "maac":
            axes[n,0].set_ylim(0, 1)
            if isinstance(maac_rv, np.ndarray):
                axes[n,0].scatter(time, maac_rv, facecolor="blue", edgecolor="k", alpha=0.2, zorder=1)

        if attr in ("slow", "bazm"):
            axes[n,0].errorbar(time, datattr[attr], yerr=bounds[attr].T, capsize=5, color="k", alpha=0.2, fmt="none", zorder=1)
            
            if attr == "slow":
                x, y = slowpdf
                axes[n,1].set_ylim(y[0], y[-1])
                axes[n,0].set_ylim(y[0], y[-1])
            else:
                x, y = bazmpdf
                axes[n,1].set_ylim(-5, 365)
                axes[n,0].set_ylim(-5, 365)
            
            normx = x/x.max()
            fig_dict[attr]["pdf"] = normx
            axes[n,1].plot(normx, y, color='k')
            # axes[n,1].grid(which="minor", axis="y", color="k", ls="-", alpha=0.35)
            axes[n,1].grid(which="major", axis="y", color="k", ls="--", alpha=0.20)
            axes[n,1].yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
            axes[n,1].yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
            axes[n,1].yaxis.set_major_formatter(mtick.NullFormatter())
            axes[n,1].xaxis.set_major_formatter(mtick.NullFormatter())
            axes[n,1].xaxis.set_major_locator(mtick.NullLocator())

        else:
            axes[n,1].axis("off")

        axes[n,0].grid(which="major", axis="x", color="k", ls="-", alpha=0.35)
        axes[n,0].grid(which="minor", axis="x", color="k", ls="--", alpha=0.20)
        axes[n,0].grid(which="major", axis="y", color="k", ls="--", alpha=0.20)

        axes[n,0].yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
        axes[n,0].yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        axes[n,0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))

        axes[n,0].xaxis.set_minor_locator(mtick.AutoMinorLocator(3))

        if attr == "bazm":
            axes[n,0].set_xlabel('Time [min]')
        else:
            axes[n,0].xaxis.set_major_formatter(mtick.NullFormatter())
        

    fig.align_labels()

    if show:
        plt.show()

    if return_fdict:
        return fig_dict
    
    else:
        return fig


def simple_slowmap(slomap, sloint, slomax, show=True, **kwargs):
    cmap  = kwargs.get("cmap", "Spectral_r")
    title = kwargs.get("title", None)
    fig   = kwargs.get("fig", None)
    axis  = kwargs.get("axis", None)
    vlim  = kwargs.get("vlim", [])
    bar_axis  = kwargs.get("bar_axis", None)
    bar_label = kwargs.get("bar_label", "CC")
    interpolation = kwargs.get("interpolation", "gaussian")

    if not axis or not fig:
        fig, axes = plt.subplots(1,2, figsize=(12,8), gridspec_kw={"width_ratios":[1,0.02]})
        axis = axes[0]
        bar_axis = axes[1]
    else:
        show = False
    
    halfbin = sloint / 2.0
    extent = (
        -slomax + halfbin,
         slomax - halfbin,
        -slomax + halfbin,
         slomax - halfbin
    )

    if not vlim:
        vmin, vmax = 0, 1
        ticks=[0,0.25,0.5,0.75,1]
    else:
        vmin, vmax = vlim
        ticks=None

    im = axis.imshow(np.flipud(slomap).T, cmap=cmap, interpolation=interpolation,\
        extent=extent, aspect='auto', vmin=vmin, vmax=vmax)

    axis.set_title(title)
    axis.set_xlabel("x [s/km]")
    axis.set_ylabel("y [s/km]")
    axis.grid(ls="--", color="k", alpha=0.3, zorder=3)

    fig.colorbar(im, cax=bar_axis, orientation='vertical', ticks=ticks, label=bar_label)

    if show:
        plt.show()

    return fig


def window_wvfm(wvfm_dict, time, startw, endw, show=True, **kwargs):

    title = kwargs.get("title", None)
    fig   = kwargs.get("fig", None)
    axes  = kwargs.get("axes", None)
    colorname = kwargs.get("colorname", "tolm")

    if not axes or not fig:
        fig, axes = plt.subplots(2,1, figsize=(9,4), sharex=True)
    else:
        show = False

    colorlist = get_colors(colorname)
    axes[0].set_title(title)

    avg_data = np.zeros(len(time))
    for n, (loc, data) in enumerate(wvfm_dict.items()):
        axes[0].plot(time, data, lw=1.2, color=colorlist[n], label=loc, alpha=0.7)
        avg_data += data/np.abs(data).max()

    # axes[0].xaxis.set_major_formatter(mtick.NullFormatter())
    
    # plot average wvfm
    axes[1].plot(time, avg_data, color="k")
    axes[1].set_xlabel("Time [sec]")
    
    for ax in axes:
        if startw and endw:
            ax.axvspan(startw, endw, color="k", alpha=0.05)
        ax.grid(which="major", axis="x", color="k", ls="-",  alpha=0.25)
        ax.grid(which="minor", axis="x", color="k", ls="--", alpha=0.15)
        ax.set_xlim(time[0], time[-1])
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.yaxis.set_major_formatter(mtick.NullFormatter())
            
    fig.legend(title="LOC", ncol=1, loc='lower right')

    if show:
        plt.show()

    return fig