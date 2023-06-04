#!/usr/bin/env python
# coding=utf-8

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

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
