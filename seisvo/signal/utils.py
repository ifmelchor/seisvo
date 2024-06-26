#!/usr/bin/env python3
# coding=utf-8

import datetime as dt
import scipy.signal as ss
import numpy as np
import math

class SSteps(object):
    """Class for fitting intervals in continuous data matching with window_length and overlap

    Parameters
    ----------
    start_time : datetime
    end_time   : datetime
    interval   : int [min]    interval in memory
    window     : float [sec]  time window of processing
    win_olap   : float [0-1)  by default is 0
    subwindow  : float [sec]  for moving average. by default is 0, i.e., no moving average to apply
    subw_olap  : float [0-1)  by default is 0
    
    """
    def __init__(self, start_time, end_time, window, interval=None,\
        win_olap=0, subwindow=0, subw_olap=0, logfile=False):
        
        if interval:
            assert window <= interval*60
        
        assert subwindow < window
        assert 0 <= win_olap < 1
        assert 0 <= subw_olap < 1
        assert end_time > start_time

        self.start_time = start_time
        self.end_time   = end_time

        # define all parameters in sec
        self.total_time = (end_time - start_time).total_seconds()

        if interval:
            self.interval   = interval*60 # in sec
        else:
            self.interval   = self.total_time

        self.window     = window
        self.win_olap   = win_olap
        self.subwindow  = subwindow
        self.subw_olap  = subw_olap

        # make some calculation
        self.win_adv  = 1 - self.win_olap
        self.subw_adv = 1 - self.subw_olap

        # compute how many intervals you need
        nro_intervals = self.total_time/self.interval

        if nro_intervals.is_integer():
            self.nro_int  = int(np.ceil(nro_intervals)) - 1
        else:
            self.nro_int  = int(np.ceil(nro_intervals))
        
        # compute number of windows for total time
        total_nwin = self.nsteps(self.total_time, self.window, self.win_olap)
        
        # compute number the maximum windows per interval
        self.int_nwin   = int(np.ceil(total_nwin/nro_intervals))
        
        # toff seconds per interval
        self.int_toff = self.window*(1 + self.win_adv*(self.int_nwin-1)) - self.interval
        
        # compute steps in subwindow
        if self.subwindow > 0:
            self.nsubwin = self.nsteps(self.window, self.subwindow, self.subw_olap)       
            self.total_subwin = int(np.floor((self.nsubwin * self.int_nwin)*nro_intervals))
        else:
            self.nsubwin = None
            self.total_subwin = None
        
        if logfile:
            self.logfile = open("ssteps.log", "w") 

        else:
            self.logfile = None

        self.fit()
        self.print()


    @staticmethod
    def nsteps(total_sec, window, olap):
        num = total_sec - olap*window
        denom = window * (1 - olap)
        return int(np.floor(num/denom))
    

    def fit(self):
        """
        Function for testing the iteration procedure in LTE and CC8 files
        """

        toff  = dt.timedelta(seconds=self.int_toff)
        olap  = dt.timedelta(seconds=float(self.win_olap*self.window))
        total_nwin = 0

        if self.nro_int > 0:
            for nint in range(self.nro_int):
                if nint >= 1:
                    start_win = end_win - olap
                else:
                    start_win = self.start_time
                
                end_win   = start_win + dt.timedelta(seconds=self.interval) + toff

                if self.logfile:
                    self.logfile.write(f"\n\n  INTERVAL {nint:>3} [nwin = {total_nwin:>5}] ENDTIME = {end_win}\n")
                
                for i in range(self.int_nwin):
                    total_nwin += 1
                    ew  = start_win + dt.timedelta(seconds=float(self.window*self.win_adv*i))
                    ewd = ew + dt.timedelta(seconds=self.window)
                    
                    if self.logfile:
                        if i in (0,1):
                            self.logfile.write(f"\n{total_nwin:>7} :: {ew.strftime('%Y-%m-%d %H:%M:%S.%f')} -- {ewd.strftime('%Y-%m-%d %H:%M:%S.%f')}")
                        elif i == 2:
                            self.logfile.write("\n                             ...")
                        elif i in (self.int_nwin-2, self.int_nwin-1):
                            self.logfile.write(f"\n{total_nwin:>7} :: {ew.strftime('%Y-%m-%d %H:%M:%S.%f')} -- {ewd.strftime('%Y-%m-%d %H:%M:%S.%f')}")
                        else:
                            pass
            # last interval
            start_win  = end_win - olap
            last_delta = self.end_time - start_win
        else:
            nint = -1
            start_win  = self.start_time
            end_win    = self.end_time
            last_delta = end_win - start_win

        if last_delta.total_seconds() > 1:
            end_win = start_win + last_delta
            nint += 1
            if self.logfile:
                self.logfile.write(f"\n\n    LAST   {nint:>3} [nwin = {total_nwin:>5}] ENDTIME = {end_win}\n")
            
            i = 0
            while start_win + dt.timedelta(seconds=self.window) <= self.end_time:
                total_nwin += 1
                ewi = start_win + dt.timedelta(seconds=self.window)
                # print(f"{total_nwin:>7} :: {start} -- {start + dt.timedelta(seconds=self.window)}")
                if self.logfile:
                    if i in (0,1):
                        self.logfile.write(f"\n{total_nwin:>7} :: {start_win.strftime('%Y-%m-%d %H:%M:%S.%f')} -- {ewi.strftime('%Y-%m-%d %H:%M:%S.%f')}")
                    elif i == 2:
                        self.logfile.write("\n                             ...")
                    else:
                        pass
                
                if start_win + dt.timedelta(seconds=float(3*self.window)) >= self.end_time and self.logfile:
                    self.logfile.write(f"\n{total_nwin:>7} :: {start_win.strftime('%Y-%m-%d %H:%M:%S.%f')} -- {ewi.strftime('%Y-%m-%d %H:%M:%S.%f')}")
                
                start_win += dt.timedelta(seconds=float(self.window-(self.win_olap*self.window)))
                i += 1
        
        else:
            i = 0
        
        self.nro_intervals = nint+1
        self.total_nwin =  total_nwin
        self.last_nwin = i
        #true_end_time = start + dt.timedelta(seconds=self.window)

        if self.logfile:
            self.logfile.close()

        return True
    

    def print(self):
        print('\n --------------- SSTEPS INFO -------------------')
        print(f'{"     Start time":^25s} :', self.start_time)
        print(f'{"       End time":^25s} :', self.end_time)
        print(f'{"   window [sec]":^25s} :', self.window)
        print(f'{"    window olap":^25s} :', self.win_olap)
        print(f'{"  Total windows":^25s} :', self.total_nwin)
        print(f'{"subwindow [sec]":^25s} :', self.subwindow)
        print(f'{" subwindow olap":^25s} :', self.subw_olap)
        print('\n --------------- INTERVALS -------------------\n')
        for k, v in [
            (f'{"          interval [min]":^25s}', self.interval/60),
            (f'{"           Nro intervals":^25s}', self.nro_intervals),
            (f'{" Nro windows in interval":^25s}', self.int_nwin),
            (f'{" Nro windows in last int":^25s}', self.last_nwin),
            (f'{"Nro subwindows in window":^25s}', self.nsubwin),
            (f'{"        Total subwindows":^25s}', self.total_subwin)]:
            print(f'{k} :', v)
        
        print('')


    def to_dict(self):
        dout = {
            "int_extra_sec":self.int_toff,
            "nro_intervals":self.nro_intervals,
            "total_nwin":self.total_nwin,
            "nwin":self.int_nwin,
            "last_nwin":self.last_nwin,
            "wadv":self.win_adv,
            "nswin":self.nsubwin,
            "swadv":self.subw_adv
        }
        return dout


def get_freq(npts, fs, fq_band=[], pad=1.0):
    assert pad >= 1, "pad must be >= 1.0"

    fftleng = round(pad*npts)
    
    if fftleng % 2 == 0:
        halffreq = int((fftleng)/2)+1
    else:
        halffreq = int(np.ceil((fftleng)/2))
    
    freq = fs*np.linspace(0,1,fftleng+1)[:halffreq+1]
    if list(fq_band):
        fnptlo = int(np.argmin(np.abs(freq-fq_band[0])))
        fnpthi = int(np.argmin(np.abs(freq-fq_band[1])))
        freq = freq[fnptlo:fnpthi+1]
    
    else:
        fnptlo = 0
        fnpthi = len(freq)
    
    return freq, (fnptlo, fnpthi)


def nearest_pow_2(x):
    """
    Find power of two nearest to x
    :type x: float
    :param x: Number
    :rtype: Int
    :return: Nearest power of 2 to x
    """
    a = int(math.pow(2, math.ceil(np.log2(x))))
    b = int(math.pow(2, math.floor(np.log2(x))))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b


def uppest_pow_2(x):
    return int(2 ** (math.ceil(math.log(x, 2))))


def get_centroid(x, y):
    """
       Get the the centroid value of y = f(x).
       Typically x :: freq and y :: PSD
    """

    assert len(x) == len(y)

    # remove nan values of x and y
    finite_args = np.isfinite(y)
    x = x[finite_args]
    y = y[finite_args]

    numerator = 0
    denominator = 0

    for xi, yi in zip(x, y):
        numerator += xi * yi
        denominator += yi

    centroid = numerator/denominator
    return centroid


def get_dominant(x, y):
    """
       Get the the value of x of the maximum y = f(x).
       Typically x :: freq and y :: PSD
    """

    assert len(x) == len(y)

    # look only at finite values
    finite_args = np.isfinite(y)
    x = x[finite_args]
    y = y[finite_args]

    dominant = x[np.argmax(y)]

    return dominant


def octave_intervals(fq_start, fq_end, delta):
    
    fq = fq_start
    full_octaves = [fq_start]
    while fq <= fq_end:
        fq *= delta
        full_octaves += [fq]
    
    return full_octaves


def smooth_psd(psd, freq, fq_band=(0.1,10), delta=2**0.125):
    """
    This code smoothes the PSD in full octaves (delta=2**0.125), given a fq_band, following Mcnamara and Boulond (2004).
    """

    fc = []
    psd_avg = []
    psd_db = 10*np.log10(psd)
    fs = octave_intervals(fq_band[0], fq_band[1], delta)

    for i in fs:
        fl = 2*i
        psd_avg += [np.mean(psd_db[np.where((freq>=i) & (freq<fl))])]
        fc += [np.sqrt(i*fl)]

    return psd_avg, fc


def corr_coeff_index(nro):
    """
    Esta function me devuelve los indices de correlacion cruzada entre pares de sensores
    :return: lista de pares de tuplas
    """
    cci = []
    for i in range(nro):
        for j in range(i,nro):
            if i != j:
                cci.append((i,j))
    return cci