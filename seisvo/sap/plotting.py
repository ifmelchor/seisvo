#!/usr/bin/env python3
# coding=utf-8

import numpy as np
import datetime as dt
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates
import matplotlib.colors as mcolor
import matplotlib.animation as mmation
from seisvo.plotting import plot_gram, get_colors
from seisvo.utils import get_time_format

default_LTE_color = get_colors('zesty')[1]


default_labels = {
    'rms':'RMS [dB]',
    'slow':'Slowness [km/s]',
    'bazm':'Back-azimuth' r'[$\degree$]',
    'maac':'MAAC',
}

default_colors = {
    "1":'k',
    "2":get_colors('zesty')[0],
    "3":get_colors('zesty')[2],
    "4":get_colors('zesty')[3]
}

default_linestyle = {
    "1":"-",
    "2":"--",
    "3":":"
}

def get_axes_dict(attr_list, fig):

    grid = {'hspace':0.15, 'left':0.08, 'right':0.92, 'wspace':0.05, 'top':0.95, 'bottom':0.05}
    col = 1
    row = len(attr_list)

    axes = fig.subplots(row, col, gridspec_kw=grid)
    
    if isinstance(axes, np.ndarray):
        axes = axes.reshape(row, col)
    
    else:
        axes = np.array([axes]).reshape(row, col)

    axes_dict = {}
    n = 0
    for attr in attr_list:
        axes_dict[attr] = axes[n,0]
        n += 1

    return axes_dict


def cc8_plot(cc8out, attr_list, fq_slo_idx, fig=None, axes=None, plot=False, return_stats=False, datetime=False, **kwargs):

    if not fig:
        fig = plt.figure(figsize=(12,8))
    
    if axes:
        clear_ax = True
    else:
        axes = get_axes_dict(attr_list, fig)
        clear_ax = False
            
    # define time and x format
    day_interval = (cc8out.endtime_ - cc8out.starttime_).total_seconds()/3600
    time_format = get_time_format(datetime, day_interval)

    if datetime:
        time = cc8out._dout["dtime"]
    else:
        time = cc8out._dout["time"]
    
    # get stats of scalar values
    fq_idx   = [fqsl.split('/')[0] for fqsl in fq_slo_idx]
    slow_idx = [fqsl.split('/')[1] for fqsl in fq_slo_idx]
    stats = cc8out.get_stats(attr_list=attr_list, fq_idx=fq_idx, slow_idx=slow_idx)
    
    for fqsl in fq_slo_idx:
        fq_idx   = fqsl.split('/')[0]
        color    = default_colors[fq_idx]
        slow_idx = fqsl.split('/')[1]
        linesty  = default_linestyle[slow_idx]

        for n, attr in enumerate(attr_list):
            ax = axes[attr]

            if clear_ax:
                ax.cla()

            key = "/".join([fqsl, attr])
            data = cc8out._dout[key]
            stdat = stats[key]

            if attr == "rms":
                data = 10*np.log10(data)
            
            ax.plot(time, data, color=color, ls=linesty, label=fqsl)
            data_min = np.floor(stdat[0])
            data_max = np.ceil(stdat[1])
            ax.set_ylim(data_min, data_max)
            ax.set_ylabel(default_labels[attr])
            ax.set_xlim(time[0], time[-1])

            ax.xaxis.set_major_locator(time_format[0][0])
            ax.xaxis.set_minor_locator(time_format[1][0])
            ax.xaxis.set_minor_formatter(mtick.NullFormatter())
            ax.xaxis.set_major_formatter(mtick.NullFormatter())
            ax.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4, min_n_ticks=3))
            ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(3))
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))

            if n == len(attr_list)-1:
                ax.xaxis.set_major_formatter(time_format[0][1])
                ax.set_xlabel('Time [h]')
                if time_format[1][1]:
                    ax.xaxis.set_minor_formatter(time_format[1][1])
            
    if plot:
        plt.show()
    
    if kwargs.get("save", False):
        filename = kwargs.get("filename", "./cc8_plot.png")
        fig.savefig(filename)
    
    if return_stats:
        return stats
    else:
        return (fig, axes)


def simple_slowness_plot(data, slomax, sloinc, title="", bar_label="CC", **kwargs):

    fig = plt.figure(figsize=(8,8))
    ax, barax = fig.subplots(1,2, gridspec_kw={"width_ratios":[1,0.02]})

    halfbin = sloinc / 2.0
    
    extent = (
        -slomax + halfbin,
         slomax - halfbin,
        -slomax + halfbin,
         slomax - halfbin
    )

    im = ax.imshow(np.flipud(data), cmap='Spectral_r', interpolation='gaussian', extent=extent, aspect='auto', vmin=0, vmax=1)
    ax.set_xlabel("px [s/km]")
    ax.set_ylabel("py [s/km]")
    ax.grid(ls="--", color="k", alpha=0.3)
    ax.set_title(title)
    fig.colorbar(im, cax=barax, orientation='vertical', ticks=[0,0.25,0.5,0.75,1], label=bar_label)

    if kwargs.get("plot", True):
        plt.show()
    
    fileout = kwargs.get("fileout", None)
    if fileout:
        fig.savefig(fileout)

    return fig 


def slowness_map_motion(cc8out, fq_slo_idx, fps, plot=True, starttime=None, endtime=None):

    fq_idx  = fq_slo_idx.split("/")[0]
    slo_idx = fq_slo_idx.split("/")[1]
    key = '/'.join([fq_slo_idx,"slowmap"])
    slowmap = cc8out._dout[key]
    dtimes  = cc8out._dout["dtime"]

    if starttime:
        n0 = np.argmin(np.abs(starttime - dtimes))
    else:
        n0 = 0
    
    if endtime:
        nf = np.argmin(np.abs(endtime - dtimes))+1
    else:
        nf = slowmap.shape[0]-1

    nf -= n0
    fig = plt.figure(figsize=(8,8))
    ax, barax = fig.subplots(1,2, gridspec_kw={"width_ratios":[1,0.02]})

    def set_title(i):
	    title = f"Fq band [{fq_idx}] : {cc8out.cc8.stats.fq_bands[int(fq_idx)-1]}  /  slow_idx : {slo_idx} \n  Time : {dtimes[i].strftime('%Y %b %d %H:%M:%S')}"
	    ax.set_title(title)
    
    slow_inc = cc8out.cc8.stats.slow_inc[int(slo_idx)-1]
    slow_max = cc8out.cc8.stats.slow_max[int(slo_idx)-1]
    halfbin = slow_inc / 2.0
    
    extent = (
        -slow_max + halfbin,
         slow_max - halfbin,
        -slow_max + halfbin,
         slow_max - halfbin
    )

    im = ax.imshow(np.flipud(slowmap[n0,:,:]), cmap='Spectral_r', interpolation='gaussian', extent=extent, aspect='auto', vmin=0, vmax=1)
    ax.set_xlabel("px [s/km]")
    ax.set_ylabel("py [s/km]")
    ax.grid(ls="--", color="k", alpha=0.3)
    set_title(n0)
    fig.colorbar(im, cax=barax, orientation='vertical', ticks=[0,0.25,0.5,0.75,1], label="CC")

    def afunc(i):
        im.set_array(np.flipud(slowmap[i,:,:]))
        set_title(i)
        return [im]
    
    anim = mmation.FuncAnimation(fig, afunc, frames=nf, interval=1000/fps)

    if plot:
        plt.show()

    return anim


def plot_slowbaz_tmap(data, title, imax=-1, fileout="slowbaz.png"):

    if imax <= 0:
        imax = len(data)
        
    fig = plt.figure(figsize=(8,8))
    ax, barax = fig.subplots(1, 2, gridspec_kw={"width_ratios":[1,0.02]})

    xlist = []
    ylist = []
    clist = []
    for i in range(imax):
        for x, y in data[i]:
            clist.append(i/imax)
            xlist.append(x)
            ylist.append(y)
    
    ax.scatter(np.array(xlist), np.array(ylist), s=15, c=np.array(clist), cmap="Spectral_r", ec="k", vmin=0, vmax=1)

    im = cm.ScalarMappable(norm=mcolor.Normalize(vmin=0,vmax=1), cmap="Spectral_r")
    fig.colorbar(im, cax=barax, orientation='vertical', ticks=[], label="time")
    ax.set_xlabel("slowness [s/km]")
    ax.set_ylabel("back azimuth")
    ax.set_ylim(0,360)
    ax.set_title(title)
    fig.savefig(fileout)


def plot_array_response(arf):
    slow_max = arf.slow_max_
    halfbin = arf.slow_inc_ / 2.0

    extent = (
        -slow_max + halfbin,
         slow_max - halfbin,
        -slow_max + halfbin,
         slow_max - halfbin
    )

    fig = plt.figure(figsize=(8,8))
    ax, barax = fig.subplots(1,2, gridspec_kw={"width_ratios":[1,0.02]})
    title = f"STA: {arf.sar_.sta_code} \n LOCS: {arf.sar_.locs} \n Fq band : {arf.fq_band_} Hz"
    ax.set_title(title)

    power = arf.power_
    power /= power.max()
    vmin = 0
    vmax = 1

    im = ax.imshow(np.flipud(power), cmap='Spectral_r', interpolation='gaussian', extent=extent, aspect='auto', vmin=vmin, vmax=vmax)

    fig.colorbar(im, cax=barax, orientation='vertical', label="Power")
    
    ax.set_xlabel("px [s/km]")
    ax.set_ylabel("py [s/km]")

    plt.show()

    return fig


def plot_array(sarray):

    fig = plt.figure(figsize=(8,8))
    ax  = fig.subplots(1,1)

    ax.scatter(sarray.yUTM, sarray.xUTM, color="k", marker="D")
    ax.set_title(sarray.sta_code)
    ax.set_ylabel("UTM [km]")
    ax.set_xlabel("UTM [km]")

    for n, loc in enumerate(sarray.locs):
        x = sarray.xUTM[n] + 10**(-6)*sarray.xUTM[n]
        y = sarray.yUTM[n] + 10**(-6)*sarray.yUTM[n]
        ax.annotate(loc, (y, x))#, bbox=dict(boxstyle="round,pad=0.3",))

    plt.show()

    return fig

    
        


            






