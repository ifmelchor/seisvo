#!/usr/bin/env python3
# coding=utf-8

import datetime as dt
import numpy as np
import math as m

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

    if fq_band:
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


class SSteps(object):
    def __init__(self, start_time, end_time, window_length, overlap, **kwargs):
        """Class for fitting intervals in continuous data matching with window_length and overlap

        Parameters
        ----------
        start_time : [datetime]
        end_time : [datetime]
        window_length : [int]
            [in sec]
        overlap : [int]
            [in sec]

        """

        if overlap >= window_length:
            raise ValueError ('Overlap must be smaller than window_length!')

        self.start_time = start_time
        self.end_time = end_time
        self.total_time = ((self.end_time - self.start_time).total_seconds()/60)
        default_interval = kwargs.get('interval', window_length) # in minutes!

        if self.total_time < default_interval:
            raise ValueError ('interval must be smaller than total time')

        self.overlap = overlap
        self.window_length = window_length
        self.sample_rate = kwargs.get('sample_rate', 100)
        
        # default interval
        self.make_interval(default_interval)
        self.fit(best_interval=False)
    

    def make_interval(self, interval):
        self.interval_ = dt.timedelta(minutes=interval)
        self.nro_intervals = self.total_time/interval


    def steps(self):
        '''Returns the number of steps of the moving window along the time interval. If appropiate, also returns the seconds that cannot be completely fitted in a window at the end of the interval.
        It is strongly recommended to use a combination of window length and overlap that return step[1]=0.\nParameters:
        window_length: an integer expressing length in seconds
        overlap: an integer expressing seconds
        sps:default value is 100 samples per second.
        '''
        
        interval_in_samples = self.interval_.total_seconds()*self.sample_rate
        distance_of_windows = (self.window_length - self.overlap)*self.sample_rate
        self.steps_ = int(interval_in_samples/distance_of_windows)
        self.rest_ = int(interval_in_samples%distance_of_windows)
    

    def last_intervalsetps(self):
        interval_in_samples = self.last_interval_.total_seconds()*self.sample_rate
        distance_of_windows = (self.window_length - self.overlap)*self.sample_rate
        self.laststeps_ = int(interval_in_samples/distance_of_windows)
        self.lastrest_ = int(interval_in_samples%distance_of_windows)


    def last_interval(self):
        '''Return the width of the last interval that will we applied to unpack the data. Parameters:
        time_series_starttime: datetime object
        time_series_endtime: datetime object
        window_length: an integer expressing length in seconds
        '''
        
        # Start time must be smaller than end time!

        new_start = self.start_time

        while new_start < self.end_time - self.interval_:
            new_start = new_start + self.interval_
        
        if new_start == self.end_time - self.interval_:
            self.last_interval_ = self.interval_
        
        else:
            self.last_interval_ = self.end_time - new_start
    

    def last_window(self):
        '''Return the width of the last window of processing, this is, the width of the last bin.\n Parameters:
        last_interval: width of the last interval (timedelta object)
        window_length: an integer expressing length in seconds
        overlap: an integer expressing seconds
        sps: default value is 100 samples per second.
        '''

        interval_in_samples = self.last_interval_.total_seconds()*self.sample_rate
        distance_of_windows = (self.window_length - self.overlap)*self.sample_rate

        rest = interval_in_samples % distance_of_windows
        if rest:
            self.last_window_ = self.window_length + rest/self.sample_rate
        
        else:
            self.last_window_ = 0


    def fit_interval(self, interval_range, verbose):
        for i in range(interval_range[0], interval_range[1]):
            self.make_interval(i)
            self.steps()
            self.last_interval()
            self.last_intervalsetps()
            self.last_window()
            if self.last_window_ == 0 and self.rest_ == 0 and self.lastrest_ == 0:
                if verbose:
                    print('\n  [fit] best interval: ', self.interval_)
                break


    def fit(self, best_interval=True, interval_range=(15,30), verbose=False):
        self.steps()
        self.last_interval()
        self.last_intervalsetps()
        self.last_window()

        if best_interval:
            self.fit_interval(interval_range, verbose)


    def total_steps(self):
        interval = self.interval_.total_seconds()/60
        last_interval = self.last_interval_.total_seconds()/60
        if last_interval != interval:
            return int(m.floor(self.nro_intervals)*self.steps_) + int(self.laststeps_)
        else:
            return int(m.floor(self.nro_intervals)*self.steps_)


    def print(self):
        print('')
        print('    -------------------------------')
        for k, v in [
            ('total_time [min]', self.total_time),
            ('interval [min]', self.interval_.total_seconds()/60),
            ('nro_intervals', self.nro_intervals),
            ('window_length [sec]', self.window_length),
            ('last_window_length [sec]', self.last_window_),
            ('overlap [sec]', self.overlap),
            ('steps_per_interval', self.steps_),
            ('steprest_per_interval', self.rest_),
            ('last_interval [min]', self.last_interval_.total_seconds()/60),
            ('last_interval_steps', self.laststeps_),
            ('last_interval_rest', self.lastrest_),
            ]:
            print(f'{k:>25s} :', v)
        print('    -------------------------------')
        print(f'{"total_steps":>25s} :', self.total_steps())
        print('')