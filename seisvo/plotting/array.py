
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
    