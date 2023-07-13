#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.dates as mdates


def particle_motion(zcomp, tcomp, rcomp, show=True, **fig_kwargs):

    baz     = fig_kwargs.get("baz", None)

    fig = plt.figure()
    ax  = fig.add_subplot(1,1,1, projection="3d")
    ax.plot(tcomp, rcomp, zcomp, color="k", label='particle motion')

    z_label = "Z"
    if baz:
        ax.set_title(f"BAZ = {baz:.1f}")
        t_label = "Transverse"
        r_label = "Radial"
    else:
        t_label = "NS"
        r_label = "EW"

    ax.set_ylabel(r_label)
    ax.xaxis.set_major_formatter(mtick.NullFormatter())
    ax.set_zlabel(z_label)
    ax.zaxis.set_major_formatter(mtick.NullFormatter())
    ax.set_xlabel(t_label)
    ax.yaxis.set_major_formatter(mtick.NullFormatter())
    
    ax.legend()

    if show:
        plt.show()

    return fig


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
