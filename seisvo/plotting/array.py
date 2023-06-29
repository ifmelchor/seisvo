
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
import numpy as np
from .utils import get_colors


def location_map(arr, exclude_locs=[], show=True):

    fig, ax = plt.subplots(1,1)

    for loc, utm in arr.utm.items():
        if loc not in exclude_locs:
            x, y = utm["easting"], utm["northing"]
            xdisp = x+0.01*x
            ydisp = y+0.01*y
            ax.scatter(x, y, marker="^", color="k")
            ax.annotate(loc, xy=(xdisp, ydisp))
    
    ax.set_ylabel("UTM northing")
    ax.set_xlabel("UTM easting")

    if show:
        plt.show()

    return fig


def traces_psd(psd_dict, freq, db_scale=True,\
    vmin=None, vmax=None, show=True, title=None, colorname="tolm"):

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

    if show:
        plt.show()

    return fig
    

def simple_cc8_plot(dtime, datattr, slowpdf, bazmpdf, show=True):

    grid = {'hspace':0.15, 'left':0.08, 'right':0.92,\
        'wspace':0.05, 'top':0.95, 'width_ratios':[1,0.1],'bottom':0.05}
    
    fig, axes = plt.subplots(4, 2, figsize=(12,8), sharex='col', gridspec_kw=grid)

    duration = (dtime[-1]-dtime[0]).total_seconds()/60
    npts     = len(dtime)
    time     = np.linspace(0, duration, npts)
    axes[0,0].set_title(f"{dtime[0]} -- {dtime[-1]}")

    default_labels = {
    'rms':'RMS [dB]',
    'slow':'Slowness [s/km]',
    'bazm':'Back-azimuth' r'[$\degree$]',
    'maac':'MAAC',
    }

    for n, attr in enumerate(["rms", "maac", "slow", "bazm"]):

        axes[n,0].scatter(time, datattr[attr], color="blue", edgecolor="k", alpha=0.6)
        axes[n,0].set_ylabel(default_labels[attr])
        axes[n,0].set_xlim(0, duration)

        if attr in ("slow", "bazm"):
            if attr == "slow":
                x, y = slowpdf
            else:
                x, y = bazmpdf
            
            normx = x/x.max()
            axes[n,1].plot(normx, y, color='k')
            axes[n,1].grid(which="minor", axis="y", color="k", ls="-", alpha=0.35)
            axes[n,1].grid(which="major", axis="y", color="k", ls="--", alpha=0.20)

        else:
            axes[n,1].axis("off")

        axes[n,0].grid(which="major", axis="x", color="k", ls="-", alpha=0.35)
        axes[n,0].grid(which="minor", axis="x", color="k", ls="--", alpha=0.20)

        axes[n,0].yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
        axes[n,0].yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        axes[n,0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
    
        # last plot
        axes[n,0].xaxis.set_minor_locator(mtick.AutoMinorLocator(3))
        if attr == "bazm":
            axes[n,0].set_xlabel('Time [min]')
        else:
            axes[n,0].xaxis.set_major_formatter(mtick.NullFormatter())

    fig.align_labels()

    if show:
        plt.show()

    return fig


def simple_slowmap(slomap, sloint, slomax, show=True, **kwargs):

    cmap = kwargs.get("cmap", "Spectral_r")
    interpolation = kwargs.get("interpolation", "gaussian")
    title = kwargs.get("title", None)

    fig, axes = plt.subplots(1,2, figsize=(12,8), gridspec_kw={"width_ratios":[1,0.02]})
    
    halfbin = sloint / 2.0
    extent = (
        -slomax + halfbin,
         slomax - halfbin,
        -slomax + halfbin,
         slomax - halfbin
    )

    im = axes[0].imshow(np.flipud(slomap), cmap=cmap, interpolation=interpolation,\
        extent=extent, aspect='auto', vmin=0, vmax=1)

    axes[0].set_title(title)
    axes[0].set_xlabel("x [s/km]")
    axes[0].set_ylabel("y [s/km]")
    axes[0].grid(ls="--", color="k", alpha=0.3)

    fig.colorbar(im, cax=axes[1], orientation='vertical', ticks=[0,0.25,0.5,0.75,1], label="CC")

    if show:
        plt.show()

    return fig



