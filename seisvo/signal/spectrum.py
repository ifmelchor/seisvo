#!/usr/bin/python3
# coding=utf-8

from __future__ import division

import numpy as np
from multitaper import MTSpec
from scipy import signal

from obspy.signal.invsim import cosine_taper
from seisvo.signal import freq_bins, _nearest_pow_2
from seisvo.plotting import plot_gram



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
    nfft = int(_nearest_pow_2(window_length * sample_rate))
    if nfft > npts:
        nfft = int(_nearest_pow_2(npts / 8.0))

    # pad with zeros to smooth the spectrogram
    mult = int(_nearest_pow_2(8.0))
    mult = mult * nfft

    # overlap
    nlap = int(nfft * float(per_lap))

    # compute the spectrogram
    freq, time, sxx = signal.spectrogram(data, fs=sample_rate,
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


def power_density_spectrum(data, sample_rate, fq_band=(), avg_step=None, olap=0, drm_params=False, **kwargs):
    """
    Compute the PSD using multitaper algortihm

    :param avg_step: step in minutes to average psd
    :type avg_step: float
    """

    time_bandwidth = kwargs.get("time_bandwidth", 3.5)
    taper = kwargs.get("taper", False)
    taper_p = kwargs.get("taper_p", 0.05)
    nfft = kwargs.get("nfft", 'uppest')
    npts = len(data)

    if avg_step:
        npts_step = int(avg_step*60*sample_rate)
        npts_olap = int(olap*npts_step)
        
        freq_n, (fnptlo, fnpthi), nfft_n = freq_bins(npts_step, sample_rate,
            fq_band=fq_band, nfft=nfft, get_freq=True)

        N = 0
        npts_start = 0
        psd_avg = np.zeros((len(freq_n),), dtype='float64')
        
        if drm_params:
            centroid = 0
            dominant = 0
            energy = 0

        while npts_start + npts_step <= npts:
            data_n = data[npts_start:npts_start+npts_step]
            data_n_mean = np.nanmean(data_n[np.isfinite(data_n)])
            data_n = data_n - data_n_mean

            if taper:
                tap = cosine_taper(npts_step, p=taper_p)
                data_n = data_n * tap

            MTSn = MTSpec(data_n, dt=1/sample_rate, nfft=nfft_n, nw=time_bandwidth).spec.reshape(-1,)
            
            psd_n = MTSn[fnptlo:fnpthi]
            psd_avg += psd_n
            
            if drm_params:
                centroid += get_centroid(freq_n, psd_n)
                dominant += get_dominant(freq_n, psd_n)
                energy += sum(psd_n)

            npts_start += npts_step - npts_olap
            N += 1
        
        to_return = [psd_avg/N, freq_n]

        if drm_params:
            to_return += [centroid/N, dominant/N, energy/N]
    
    else:
        freq, (fnptlo, fnpthi), nfft = freq_bins(npts, sample_rate,
        fq_band=fq_band, nfft=nfft, get_freq=True)

        if taper:
            tap = cosine_taper(npts, p=taper_p)
            data = data * tap

        MTSy = MTSpec(data, dt=1/sample_rate, nfft=nfft, nw=time_bandwidth).spec.reshape(-1,)
        psd = MTSy[fnptlo:fnpthi]

        to_return = [psd, freq]

        if drm_params:
            centroid = get_centroid(freq, psd)
            dominant = get_dominant(freq, psd)
            energy = sum(psd)
            to_return += [centroid, dominant, energy]

    return to_return


def cross_spectrum(xdata, ydata, sample_rate, avg_step=None, fq_band=(), **kwargs):
    """
    Compute the cross spectum between xdata and ydata using the multitaper algortihm

    :param avg_step: step in minutes to average psd
    :type avg_step: int
    """

    time_bandwidth = kwargs.get("time_bandwidth", 3.5)
    taper = kwargs.get("taper", False)
    taper_p = kwargs.get("taper_p", 0.05)
    nfft = kwargs.get("nfft", 'uppest')

    npts = len(xdata)
    delta = float(1/sample_rate)
    nro_tapers = int(2*time_bandwidth)-1

    if avg_step:
        npts_step = int(avg_step*60*sample_rate-1)
        if npts_step > npts:
            raise ValueError(' step is greater than data lenght!')

        N = 1
        npts_start = 0
        while npts_start + npts_step <= npts:
            xdata_n = xdata[npts_start:npts_start+npts_step]
            xdata_n_mean = np.nanmean(xdata_n[np.isfinite(xdata_n)])
            xdata_n = xdata_n - xdata_n_mean

            ydata_n = ydata[npts_start:npts_start+npts_step]
            ydata_n_mean = np.nanmean(ydata_n[np.isfinite(ydata_n)])
            ydata_n = ydata_n - ydata_n_mean

            if taper:
                tap = cosine_taper(npts_step, p=taper_p)
                xdata_n = xdata_n * tap
                ydata_n = ydata_n * tap

            freq_n, (fnptlo, fnpthi), nfft_n = freq_bins(npts_step, sample_rate,
                fq_band=fq_band, nfft=nfft, get_freq=True)

            MTSx = MTSpec(xdata_n, dt=delta, nfft=nfft_n, nw=time_bandwidth)
            x = MTSx.yk

            MTSy = MTSpec(ydata_n, dt=delta, nfft=nfft_n, nw=time_bandwidth)
            y = MTSy.yk

            psd_xy_n = np.zeros((len(freq_n),), dtype='complex128')
            for k in range(nro_tapers):
                psd_xy_n += x[fnptlo:fnpthi, k] * np.conj(y[fnptlo:fnpthi, k])
            psd_xy_n = psd_xy_n/nro_tapers

            if N == 1:
                psd_xy_avg = np.array(psd_xy_n, dtype='complex128')
            else:
                psd_xy_avg += psd_xy_n

            npts_start += npts_step
            N += 1

        return psd_xy_avg/N, freq_n

    else:
        xdata_mean = np.nanmean(xdata[np.isfinite(xdata)])
        xdata = xdata - xdata_mean

        ydata_mean = np.nanmean(ydata[np.isfinite(ydata)])
        ydata = ydata - ydata_mean

        if taper:
            tap = cosine_taper(npts, p=taper_p)
            xdata = xdata * tap
            ydata = ydata * tap

        freq, (fnptlo, fnpthi), nfft = freq_bins(npts, sample_rate,
            fq_band=fq_band, nfft=nfft, get_freq=True)

        MTSx = MTSpec(xdata, dt=delta, nfft=nfft_n, nw=time_bandwidth)
        x = MTSx.yk

        MTSy = MTSpec(ydata, dt=delta, nfft=nfft_n, nw=time_bandwidth)
        y = MTSy.yk

        psd_xy = np.zeros((len(freq),), dtype='complex128')
        for k in range(nro_tapers):
            psd_xy += x[fnptlo:fnpthi, k] * np.conj(y[fnptlo:fnpthi, k])
        
        # Following Prieto et al 2009, no division with K should be applied.
        # psd_xy = psd_xy/nro_tapers

        return psd_xy, freq


def cross_spectrum_matrix(idata, jdata, kdata, sample_rate, avg_step=None, fq_band=(), **kwargs):
    """
    This code compute the cross spectral matrix using multitaper algorithm
    """

    from itertools import product

    data = [idata, jdata, kdata]
    nfft = kwargs.get("nfft", 'uppest')

    if avg_step:
        npts = int(avg_step*60*sample_rate)
    else:
        npts = len(idata)

    freq, _, _ = freq_bins(npts, sample_rate, fq_band=fq_band, nfft=nfft, get_freq=True)
    spec_mat = np.zeros((3, 3, freq.shape[0]), dtype='complex128')

    index = set(product(set(range(3)), repeat=2))
    for i, j in list(index):
        ans = cross_spectrum(data[i], data[j], sample_rate, avg_step=avg_step, fq_band=fq_band, **kwargs)
        spec_mat[i,j] = ans[0]

    return spec_mat, freq


def get_centroid(freq, psd):
    """
       Get the centroid freq of the psd
    """

    if len(freq) != len(psd):
        raise ValueError(' freq and psd are of different length!')

    # remove nan values of psd and freq
    finite_args = np.isfinite(psd)
    freq = freq[finite_args]
    psd = psd[finite_args]

    numerator = 0
    denominator = 0

    for fq_i, psd_i in zip(freq, psd):
        numerator += fq_i * psd_i
        denominator += psd_i

    centroid = numerator/denominator
    return centroid


def get_dominant(freq, psd):
    """
       Get the dominant freq of the psd
    """

    if len(freq) != len(psd):
        raise ValueError(' freq and psd are of different length!')

    # look only at finite values
    finite_args = np.isfinite(psd)
    freq = freq[finite_args]
    psd = psd[finite_args]

    max_psd = np.argmax(psd)
    dominant = freq[max_psd]

    return dominant


def smooth_psd(psd, freq, fq_band=(0.1,10)):
    """
    This code smoothes the PSD in full octaves, given a fq_band, following Mcnamara and Boulond (2004).
    """

    def octave_intervals(fq_start, fq_end):
        fq = fq_start
        delta = 2**0.125
        full_octaves = [fq_start]
        while fq <= fq_end:
            fq = fq*2**(0.125)
            full_octaves += [fq]
        return full_octaves

    fc = []
    psd_avg = []
    psd_db = 10*np.log10(psd)
    fs = octave_intervals(fq_start=fq_band[0], fq_end=fq_band[1])

    for i in fs:
        fl = 2*i
        psd_avg += [np.mean(psd_db[np.where((freq>=i) & (freq<fl))])]
        fc += [np.sqrt(i*fl)]

    return psd_avg, fc
