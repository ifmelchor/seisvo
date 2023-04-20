#!/usr/bin/env python3
# coding=utf-8

import sys
import datetime as dt
import numpy as np
import scipy
import math as m
import multiprocessing as mp
from .proba import get_PDF

def _nearest_pow_2(x):
    """
    Find power of two nearest to x
    :type x: float
    :param x: Number
    :rtype: Int
    :return: Nearest power of 2 to x
    """
    a = int(m.pow(2, m.ceil(np.log2(x))))
    b = int(m.pow(2, m.floor(np.log2(x))))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b


def _uppest_pow_2(x):
    return int(2 ** (m.ceil(m.log(x, 2))))


def get_tseries_stats(x):
    v_min = x.min()
    v_max = x.max()
    v_mean = x.mean()

    x_range = np.linspace(v_min, v_max, 500)
    gkde = scipy.stats.gaussian_kde(x)
    kde = gkde(x_range)
    v_mode = x_range[np.argmax(kde)]

    return [v_min, v_max, v_mean, v_mode]


def get_pdf_data(x, bandwidth, db_scale=False, **kwargs):

    if len(x.shape) == 1:
        x = x.reshape(-1,1)
        
    masked_data = np.ma.masked_invalid(x)
    data = np.ma.compress_rowcols(masked_data, axis=0)

    if db_scale:
        data = 10*np.log10(data)
    
    xmin = kwargs.get("xmin", None)
    if not xmin:
        xmin = data.min()
    
    xmax = kwargs.get("xmax", None)
    if not xmax:
        xmax = data.max()

    xspace = kwargs.get("xspace", 1000)
    y_space = np.linspace(xmin, xmax, xspace).reshape(xspace, 1)
    pdf = get_PDF(data, y_space, bandwidth, **kwargs)

    return y_space.reshape(xspace,), pdf


def get_freq(npts, fs, fq_band=[], pad=1.0):
    assert pad >= 1, "pad must be >= 1.0"

    fftleng = round(pad*npts)
    
    if fftleng % 2 == 0:
        halffreq = int((fftleng)/2)+1
    else:
        halffreq = int(np.ceil((fftleng)/2))
    
    freq = fs*np.linspace(0,1,fftleng+1)[:halffreq+1]
    if list(fq_band):
        fnptlo = np.argmin(np.abs(freq-fq_band[0]))
        fnpthi = np.argmin(np.abs(freq-fq_band[1]))
        freq = freq[fnptlo:fnpthi+1]
    
    else:
        fnptlo = 0
        fnpthi = len(freq)
    
    return freq, (fnptlo, fnpthi)


def freq_bins(npts, sampling_rate, fq_band=[], nfft='uppest', get_freq=False):
    """
    Get the nro of freq bins for a nfft
    :type nspt: int
    :type sampling_rate: int
    :type fq_band: tuple of two floats
    :type nfft: str. 'uppest' or 'nearest'
    :return: Nro of bins
    """

    if not isinstance(nfft, int):
        if nfft == 'uppest':
            nfft = _uppest_pow_2(npts)

        elif nfft == 'nearest':
            nfft = _nearest_pow_2(npts)

        else:
            raise TypeError('nfft must be int, "uppest" or "nearest"')
    
    freq = np.fft.fftfreq(nfft, d=1/sampling_rate)[:int(nfft/2)]

    if list(fq_band):
        fnptlo = np.argmin(np.abs(freq-fq_band[0]))
        fnpthi = np.argmin(np.abs(freq-fq_band[1]))
        freq = freq[fnptlo:fnpthi]
    else:
        fnptlo = 0
        fnpthi = len(freq)

    if get_freq:
        return freq, (fnptlo, fnpthi), nfft
    
    else:
        return len(freq)


def time_bins(start_time, end_time, interval, olap):
    total_time = (end_time - start_time).total_seconds()/60
    n = ((total_time/interval) - olap) / (1 - olap)
    return int(n)


def degree_mean(x, ambiguity180=False):
    # convert to radian an computes sin(x)
    xrad = np.array(x)*np.pi/180
    
    if ambiguity180:
        xrad_mean = np.arctan(np.sin(xrad).sum()/np.abs(np.cos(xrad)).sum())
        xrad_std = np.sin(xrad).std()

    else:
        # the ambiguity is for 360 degree
        xrad_mean = np.arctan(np.sin(xrad).sum()/np.cos(xrad).sum())
        xrad_std = np.sin(xrad).std()

    return xrad_mean*180/np.pi, xrad_std*180/np.pi
    

class SSteps(object):
    """Class for fitting intervals in continuous data matching with window_length and overlap

    Parameters
    ----------
    start_time : datetime
    end_time   : datetime
    interval   : int [min]    interval in memory
    window     : float [sec]  time window of processing
    win_olap   : float [0-1)  by default is 0
    subwindow  : int [sec]    for moving average. by default is 0
    subw_olap  : float [0-1)  by default is 0
    
    """
    def __init__(self, start_time, end_time, interval, window, win_olap=0, subwindow=0, subw_olap=0, verbose=False, validate=False):
        
        if interval != -1:
            assert window/60 <= interval
        
        assert subwindow < window
        assert 0 <= win_olap < 1
        assert 0 <= subw_olap < 1
        assert end_time > start_time

        self.start_time = start_time
        self.end_time   = end_time

        # define all parameters in sec
        self.total_time = (end_time - start_time).total_seconds()

        if interval == -1:
            self.interval   = self.total_time
        
        else:
            self.interval   = interval*60 # in sec

        self.window     = window # in sec
        self.win_olap   = win_olap
        self.subwindow  = subwindow
        self.subw_olap  = subw_olap

        self.win_adv  = 1 - self.win_olap
        self.subw_adv = 1 - self.subw_olap

        self.fit()

        if validate:
            self.print()

    
    @staticmethod
    def nsteps(total_sec, window, olap):
        num = total_sec - olap*window
        denom = window * (1 - olap)
        steps = np.floor(num/denom)
        return steps
    

    def fit(self):
        # compute how many intervals you need
        self.nro_intervals = self.total_time/self.interval
        
        # compute steps in window
        self.total_nwin = self.nsteps(self.total_time, self.window, self.win_olap)
        self.int_nwin   = self.nsteps(self.interval, self.window, self.win_olap)
        
        # compute steps in subwindow
        if self.subwindow > 0:
            self.nsubwin = self.nsteps(self.window, self.subwindow, self.subw_olap)       
            self.total_subwin = self.nsubwin * self.int_nwin
        else:
            self.nsubwin = None
            self.total_subwin = None

        self.nwin_eff = self.int_nwin * self.nro_intervals
    
    
    def fit_interval(self, intmin, intmax, intrange):

        print( f"\n  starttime :: {self.start_time}  |  endtime :: {self.end_time} ")
        print( f" window [sec] :: {self.window}    |    overlap :: {self.win_olap}  ")
        print( f" total nwin :: {self.total_nwin} ") 
        print( "\n  int [min]  ::   nro int   ::  to compute  :: noComputed  (per int)")
        for interval in range(intmin, intmax+1, intrange):
            int_nwin = self.nsteps(interval*60, self.window, self.win_olap)
            nro_int  = self.total_time/(interval*60)
            to_compute   = int_nwin * nro_int
            not_computed = self.total_nwin-to_compute
            not_computed_per = (self.total_nwin-to_compute)/nro_int
            last_col = f"{not_computed:.1f} ({not_computed_per:.2f})"
            ans = f"{interval:^13.1f}::{nro_int:^13.1f}::{to_compute:^14.1f}::{last_col:^22s}"
            print(ans)
        print("")
    

    def print(self):
        print('')
        print(' SSteps INFO')
        print(' -------------')
        print(f'{"     Start time":^25s} :', self.start_time)
        print(f'{"       End time":^25s} :', self.end_time)
        print(f'{"   window [sec]":^25s} :', self.window)
        print(f'{"    window olap":^25s} :', self.win_olap)
        print(f'{"  Total windows":^25s} :', self.total_nwin)
        print(f'{"subwindow [sec]":^25s} :', self.subwindow)
        print(f'{" subwindow olap":^25s} :', self.subw_olap)

        print('\n    ------------ INTERVALS -------------------\n')
        
        for k, v in [
            (f'{"          interval [min]":^25s}', self.interval/60),
            (f'{"           Nro intervals":^25s}', self.nro_intervals),
            (f'{" Nro windows in interval":^25s}', self.int_nwin),
            (f'{"  Nro windows to compute":^25s}', self.nwin_eff),
            (f'{" Nro windows NO computed":^25s}', self.total_nwin-self.nwin_eff),
            (f'{"       |-> per intterval":^25s}', (self.total_nwin-self.nwin_eff)/self.nro_intervals),
            (f'{"Nro subwindows in window":^25s}', self.nsubwin),
            (f'{"        Total subwindows":^25s}', self.total_subwin)]:
            print(f'{k} :', v)
        
        print('')


