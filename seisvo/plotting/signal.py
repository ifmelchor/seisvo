

import scipy.signal as ss
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from seisvo.signal.utils import nearest_pow_2
from .utils import plot_gram


def spectrogram(data, sample_rate, axes, per_lap=0.75, window_length=None, fq_band=(), date_list=None, **kwargs):
    """
    Code to plot the spectrogram.
    The code returns AxesImage and the upper and lower value of the norm.
    kwargs: rel_norm (False), v_max (None), v_min (None), interpolation ('gaussian')
    fq_logscale (False), axis_bar, axis_bar_label
    """

    npts = len(data)

    if not window_length:
        window_length = int(sample_rate)//2

    # defining nfft (windows of the FFT transfrom)
    nfft = int(nearest_pow_2(window_length * sample_rate))
    if nfft > npts:
        nfft = int(nearest_pow_2(npts / 8.0))

    # pad with zeros to smooth the spectrogram
    mult = int(nearest_pow_2(8.0))
    mult = mult * nfft

    # overlap
    nlap = int(nfft * float(per_lap))

    # compute the spectrogram
    freq, time, sxx = ss.spectrogram(data, fs=sample_rate,
        nperseg=nfft, nfft=mult, noverlap=nlap, scaling='spectrum')

    if fq_band != ():
        fnptlo = np.argmin(np.abs(freq-fq_band[0]))
        fnpthi = np.argmin(np.abs(freq-fq_band[1]))
        freq = freq[fnptlo:fnpthi]
        sxx = sxx[fnptlo:fnpthi,:]

    if not date_list:
        date_list = time
        kwargs['is_time'] = False
    else:
        kwargs['is_time'] = True

    # normalize spectrogram
    db_norm = kwargs.get('db_norm', True)
    if db_norm:
        sxx = 10*np.log10(sxx)

    im, (v_min, v_max) = plot_gram(freq, sxx, date_list, axes, **kwargs)

    return im, (v_min, v_max)


def plotPDF(pdf, y_bins, x_bins, axis=None, plot=True, **kwargs):
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

    # masked = np.ma.masked_where(pdf<1e-06, pdf)
    plot_gram(y_bins, pdf, x_bins, axis, **kwargs)

    if kwargs.get('show_models', False):
        _, nlnm = get_nlnm() # NLNM model
        p, nhnm = get_nhnm() # NHNM model

        axis.plot(1/p, nlnm, color='k', lw=1.2, label='NLNM', zorder=10)
        axis.plot(1/p, nhnm, color='k', lw=1.2, label='NHNM', zorder=10)

    if plot:
        plt.show()

    return fig