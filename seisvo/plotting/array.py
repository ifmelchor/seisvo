
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from .utils import get_colors

def _slowbnds_error(x, th):
    xmax = np.argmax(x)
    x1 = np.argmin(np.abs(x[:xmax]-th))
    x2 = np.argmin(np.abs(x[xmax:]-th))
    return (x1, xmax+x2)


def _detections(count, time):
    count = np.where(count<1, np.nan, count)

    fig, ax = plt.subplots(1,1)
    ax.plot(time, count, color="k")
    ax.set_ylabel("nro_detection/hours")
    ax.set_xlim(time[0],time[-1])

    return fig


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


def traces_psd(psd_dict, freq, db_scale=False, vmin=None, vmax=None, show=False, title=None, colorname="tolm"):

    fig, ax = plt.subplots(1,1, figsize=(10,8))
    colorlist = get_colors(colorname)

    ymin = 999
    ymax = 0

    for n, (loc, psd) in enumerate(psd_dict.items()):

        if db_scale:
            psd = 10*np.log10(psd)

        if ymin > np.nanmin(psd):
            ymin = np.nanmin(psd)

        if ymax < np.nanmax(psd):
            ymax = np.nanmax(psd)

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
        ax.set_xscale("log")
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
    

def simple_cc8_plot(dtime, datattr, datapdf, show=True, **kwargs):
    fig          = kwargs.get("fig", None)
    title        = kwargs.get("title", None)
    maac_rv      = kwargs.get("maac_rv", None)
    rms_rv       = kwargs.get("rms_rv", None)
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
    'baz':'Back-azimuth' r'[$\degree$]',
    'maac':'MAAC',
    }

    fig_dict = {}
    for n, attr in enumerate(["maac", "rms", "slow", "baz"]):
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
            if isinstance(rms_rv, np.ndarray):
                axes[n,0].scatter(time, rms_rv, facecolor="w", edgecolor="k", alpha=0.2, zorder=1)
            if isinstance(rms_rv, list):
                for r in rms_rv:
                    axes[n,0].scatter(time, r, facecolor="w", edgecolor="k", alpha=0.2, zorder=1)

        if attr == "maac":
            axes[n,0].set_ylim(0, 1)
            if isinstance(maac_rv, np.ndarray):
                axes[n,0].scatter(time, maac_rv, facecolor="w", edgecolor="k", alpha=0.2, zorder=1)
            if isinstance(maac_rv, list):
                for m in maac_rv:
                    axes[n,0].scatter(time, m, facecolor="w", edgecolor="k", alpha=0.2, zorder=1)

        if attr in ("slow", "baz"):
            bound = datattr[attr+"bnd"]
            yerr = np.abs(bound[:,1]-bound[:,0])/2
            axes[n,0].errorbar(time, datattr[attr], yerr=yerr, capsize=5, color="k", alpha=0.2, fmt="none", zorder=1)
            x, y = datapdf[attr]
            
            if attr == "slow":
                axes[n,1].set_ylim(y[0], y[-1])
                axes[n,0].set_ylim(y[0], y[-1])
            else:
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

        if attr == "baz":
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


def simple_slowmap(slomap, sloint, slomax, slow0, ccerr=0.9, show=True, **kwargs):
    cmap  = kwargs.get("cmap", "Spectral_r")
    title = kwargs.get("title", None)
    fig   = kwargs.get("fig", None)
    axis  = kwargs.get("axis", None)
    vlim  = kwargs.get("vlim", [])
    bar_axis  = kwargs.get("bar_axis", None)
    bar_label = kwargs.get("bar_label", "C")
    figsize   = kwargs.get("figsize", (12,8))
    interpolation = kwargs.get("interpolation", "gaussian")

    if not axis or not fig:
        fig, axes = plt.subplots(1,2, figsize=figsize, gridspec_kw={"width_ratios":[1,0.02]})
        axis = axes[0]
        bar_axis = axes[1]
    else:
        show = False
    
    extent = (-abs(slow0[0])-slomax, 
               abs(slow0[0])+slomax, 
              -abs(slow0[0])-slomax, 
               abs(slow0[0])+slomax)
    
    slomap = np.flipud(slomap).T
    alphas = mcolor.Normalize(0, ccerr, clip=True)(np.abs(slomap))
    alphas = np.clip(alphas, ccerr, 1)

    nite = 1 + 2*int(slomax/sloint)
    slox = np.linspace(-slomax, slomax, nite)
    
    ticks = [0,0.25,0.5,0.75,1] # MAAC scale from 0 to 1
    
    maacth = slomap.max()*ccerr
    im = axis.imshow(slomap, cmap=cmap, interpolation=interpolation, extent=extent, aspect='auto', vmin=0, vmax=1, zorder=1)
    axis.contour(slox, -slox, slomap, levels=[maacth], colors="r")
    
    maxpos = np.where(slomap==slomap.max())
    slov0x = np.linspace(0, slox[maxpos[1]], 100)
    slov0y = -np.linspace(0, slox[maxpos[0]], 100)

    axis.plot(slov0x, slov0y, ls="--", lw=0.8, color="k", alpha=0.7, zorder=2)
    axis.scatter(slox[maxpos[1]],-slox[maxpos[0]], marker="o", color="r", ec="k", zorder=3)
    axis.scatter(0, 0, marker="o", color="k", ec="k", zorder=3)

    # plot error bars
    # y1, y2 = _slowbnds_error(slomap[maxpos[0],:].reshape(-1,), maacth)
    # x1, x2 = _slowbnds_error(slomap[:,maxpos[1]].reshape(-1,), maacth)
    # axis.errorbar(slox[maxpos[0]],slox[maxpos[1]], yerr=abs(slox[y2]-slox[y1]), xerr=abs(slox[x2]-slox[x1]), capsize=5, color="k", lw=0.8, fmt="none", zorder=2)
    # plot points of the bounds
    # axis.scatter(slox[maxpos[1]], -slox[y1], marker="o", color="k", ec="k", alpha=0.25, zorder=3)
    # axis.scatter(slox[maxpos[1]], -slox[y2], marker="o", color="k", ec="k", alpha=0.25, zorder=3)
    # axis.scatter(slox[x1], -slox[maxpos[0]], marker="o", color="k", ec="k", alpha=0.25, zorder=3)
    # axis.scatter(slox[x2], -slox[maxpos[0]], marker="o", color="k", ec="k", alpha=0.25, zorder=3)
    # plot error cross
    # slov1x = np.linspace(slox[x1],slox[x2],100)
    # slov1y = 100*[slox[maxpos[0]]]
    # axis.plot(slov1x, slov1y, ls="--", lw=0.8, color="k", alpha=0.25, zorder=2)
    # slov2x = 100*[slox[maxpos[1]]]
    # slov2y = np.linspace(slox[y1],slox[y2],100)
    # axis.plot(slov2x, slov2y, ls="--", lw=0.8, color="k", alpha=0.25, zorder=2)

    axis.set_title(title)
    axis.set_xlabel("x [s/km]")
    axis.set_ylabel("y [s/km]")
    axis.xaxis.set_minor_locator(mtick.MultipleLocator(sloint))
    axis.yaxis.set_minor_locator(mtick.MultipleLocator(sloint))
    axis.xaxis.set_major_locator(mtick.MultipleLocator(0.5))
    axis.yaxis.set_major_locator(mtick.MultipleLocator(0.5))
    axis.grid(which="major", ls="-", lw=0.8, color="k", alpha=0.3, zorder=3)
    axis.grid(which="minor", ls=":", lw=0.5, color="k", alpha=0.3, zorder=3)

    fig.colorbar(im, cax=bar_axis, orientation='vertical', ticks=ticks, label=bar_label)

    if show:
        plt.show()

    return fig


def beamform_wvfm(wvfm_dict, suma, time, show=True, **kwargs):

    title = kwargs.get("title", None)
    fig   = kwargs.get("fig", None)
    axes  = kwargs.get("axes", None)
    shadow_times = kwargs.get("shadow_times", ())
    colorname = kwargs.get("colorname", "tolm")

    if not axes or not fig:
        fig, axes = plt.subplots(2,1, figsize=(9,4), sharex=True)
    else:
        show = False

    colorlist = get_colors(colorname)
    axes[0].set_title(title)
    for n, (loc, data) in enumerate(wvfm_dict.items()):
        axes[0].plot(time, data, lw=1.2, color=colorlist[n], label=loc, alpha=0.7)
    
    # axes[0].xaxis.set_major_formatter(mtick.NullFormatter())
    
    # plot average wvfm
    axes[1].plot(time, suma, color="k")
    axes[1].set_xlabel("Time [sec]")
    
    for ax in axes:
        if shadow_times:
            ax.axvspan(shadow_times[0], shadow_times[1], color="k", alpha=0.2)
        ax.grid(which="major", axis="x", color="k", ls="-",  alpha=0.25)
        ax.grid(which="minor", axis="x", color="k", ls="--", alpha=0.15)
        ax.set_xlim(time[0], time[-1])
        ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(2))
        ax.yaxis.set_major_formatter(mtick.NullFormatter())
            
    fig.legend(title="LOC", ncol=1, loc='lower right')

    if show:
        plt.show()

    return fig


def plot_slowmap(array, starttime, window, offsec=3, taper=True, slowarg={}):

    # create frame
    fig  = plt.figure(figsize=(12,9))
    gs   = GridSpec(3, 4, wspace=0.5, width_ratios=[1, 2, 0.1, 1], height_ratios=[2, 1, 1])
    ax1  = fig.add_subplot(gs[0,1])
    ax2  = fig.add_subplot(gs[0,2])
    ax3  = fig.add_subplot(gs[1,:])
    ax4  = fig.add_subplot(gs[2,:])
    _, ans0 = array.slowmap(starttime, window, slowarg=slowarg, plot=True, axis=ax1, fig=fig, bar_axis=ax2)

    slow = ans0["slow"][0]
    baz  = ans0["baz"][0]

    time0 = starttime - dt.timedelta(seconds=offsec) 
    time1 = starttime + dt.timedelta(seconds=window+offsec)
    duration = (time1 - time0).total_seconds()
    startw = duration/2 - window/2
    endw   = startw + window
    array.beamform(time0, time1, slow, baz, shadow_times=(startw, endw), taper=taper, slowarg=slowarg, plot=True, fig=fig, axes=[ax3,ax4])

    plt.show()

    return fig
