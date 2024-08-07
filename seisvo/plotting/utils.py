#!/usr/bin/env python
# coding=utf-8

import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.dates as mdates
import datetime as dt
import numpy as np


def get_time_format(datetime, day_interval):
    if datetime:
        if day_interval <= 1:
            major_locator = mdates.HourLocator(interval=1)
            major_formatt = mdates.DateFormatter('%d %b\n%H:%M')
            minor_locator = mdates.MinuteLocator(byminute=[15, 30, 45])
            minor_formatt = mtick.NullFormatter()

        elif day_interval <= 10:
            major_locator = mdates.DayLocator(interval=1)
            major_formatt = mdates.DateFormatter('%d %b\n%H:%M')
            minor_locator = mdates.HourLocator(byhour=[6, 12, 18, 24])
            minor_formatt = mtick.NullFormatter()

        elif 45 >= day_interval > 10 :
            major_locator = mdates.DayLocator(interval=3)
            major_formatt = mdates.DateFormatter('%d')
            minor_locator = mdates.DayLocator(interval=1)
            minor_formatt = mtick.NullFormatter()

        else:
            major_locator = mdates.WeekdayLocator(interval=2)
            major_formatt = mdates.DateFormatter('%d-%m\n%H:%M')
            minor_locator = mdates.DayLocator(interval=7)
            minor_formatt = mtick.NullFormatter()
    
    else:
        major_locator = mtick.LinearLocator(10)
        major_formatt = mtick.FormatStrFormatter('%i')
        minor_locator = mtick.AutoMinorLocator(2)
        minor_formatt = None
    
    return (major_locator, major_formatt), (minor_locator, minor_formatt)


def truncate_cmap(cmap, minval=0.0, maxval=1.0, n=100):
    '''
    from https://stackoverflow.com/a/18926541
    '''
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    new_cmap = mcolor.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    
    return new_cmap


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


def plot_gram(y, array, x, axis, **kwargs):
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

    dx2 = 0#(x2 - x1) / 2.0
    dy2 = 0#(y[1] - y[0]) / 2.0
    extent = (
        x1 + dx2,
        xFI - dx2,
        y[0] + dy2,
        y[-1] - dy2
    )

    norm = kwargs.get('norm', None)
    if not norm:
        v_max = kwargs.get('v_max', None)
        if not v_max:
            v_max = np.nanmax(array[np.isfinite(array)])

        v_min = kwargs.get('v_min', None)
        if not  v_min:
            v_min = np.nanmin(array[np.isfinite(array)])

        norm = mcolor.Normalize(v_min, v_max)

    interpolation = kwargs.get('interpolation', 'spline36')
    cmap = kwargs.get('cmap', 'Spectral_r')
    im = axis.imshow(np.flipud(array), cmap=cmap, norm=norm, interpolation=interpolation, extent=extent)#, aspect="auto")
    axis.axis('tight')
    axis.grid(which="major", color="k", ls="-", alpha=0.3, zorder=4)
    axis.grid(which="minor", color="k", ls="--", alpha=0.3, zorder=4)
    
    if is_time:
        axis.xaxis_date()

    fq_logscale = kwargs.get('fq_logscale', False)
    if fq_logscale:
        axis.set_xscale('log')

    y_label = kwargs.get('y_label', None)
    
    if y_label:
        axis.set_ylabel(y_label)
    x_label = kwargs.get('x_label', None)
    
    if x_label:
        axis.set_xlabel(x_label)

    axis_bar = kwargs.get('axis_bar', None)
    orientation = kwargs.get('orientation', 'vertical')
    if axis_bar:
        axis_bar_label = kwargs.get('bar_label', None)
        fig = axis.get_figure()
        cbar = fig.colorbar(im, cax=axis_bar, orientation=orientation)
        cbar.locator = mtick.MaxNLocator(nbins=4)
        cbar.update_ticks()
        cbar.set_label(axis_bar_label)

        return cbar