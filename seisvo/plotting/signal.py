

import numpy as np
from scipy.signal import spectrogram as sss
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
    freq, time, sxx = sss(data, fs=sample_rate,
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

