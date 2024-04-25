#!/usr/bin/python3
# seisvo
'''
 This module extend the obspy modules
'''

import scipy
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import obspy.signal.invsim as osi
from obspy import Stream, Trace, read, UTCDateTime
from obspy.signal.util import _npts2nfft
from obspy.signal import filter as osf
from .signal import get_PSD, get_TimePolar
from .plotting.signal import spectrogram


def read2(*args, **kwargs):
    st = read(*args, **kwargs)
    return Stream2(st)


class Trace2(Trace):
    def __init__(self, trace):
        super().__init__(data=trace.data, header=trace.stats)
    

    def get_data(self, starttime=None, endtime=None, demean=True, detrend=False, fq_band=(), abs=False, sample_rate=None, norm=False):
        """
        This code returns a numpy array of the data
        """

        if starttime and endtime:
            st = UTCDateTime(starttime)
            et = UTCDateTime(endtime)
            self = self.slice(st, et)
        
        if sample_rate:
            self = self.resample(sample_rate)
        
        data = self.data

        if detrend:
            demean = True
            
        if demean:
            data = data - data.mean()

        if detrend:
            data = scipy.signal.detrend(data)

        if list(fq_band):
            tr_filt = self.filter2(fq_band)
            data = tr_filt.data
        
        if abs:
            data = np.abs(data)
        
        if norm:
            maxval = np.abs(data).max()
            data = data / maxval
        
        return data


    def get_time(self, starttime=None, endtime=None):
        if not starttime:
            starttime = self.stats.starttime.datetime

        if not endtime:
            endtime = self.stats.endtime.datetime

        delta = self.stats.delta
        npts = int((endtime-starttime).total_seconds()*self.stats.sampling_rate)+1
        timelist = [starttime + dt.timedelta(seconds=delta*x) for x in range(npts)]
        return np.array(timelist, dtype=dt.datetime)


    def plot_specgram(self, axis, date_list=None, window_length=None, fq_band=(), per_lap=0.75, show=False, **kw_plotgram):
        """
        This code plots the spectrogram of the trace at specific AxisPlot
        cax is the axis to show the colobar
        """

        data = self.get_data(detrend=True)
        spectrogram(data, self.stats.sampling_rate, axis, per_lap=per_lap,
            window_length=window_length, fq_band=fq_band, date_list=date_list, **kw_plotgram)
        
        if show:
            plt.show()

        return
    

    def plot_trace(self, axis, fq_band=(), show=False):
        data = self.get_data(fq_band=fq_band)
        time = self.get_time()

        axis.plot(time, data, color="k", lw=0.8)
        axis.set_xlim(time[0], time[-1])
        
        axis.annotate(f'{self.id}', xy=(0.05,0.9), xycoords='axes fraction', color="k", fontfamily="monospace", stretch='condensed')
        axis.grid(axis='both', which='major', color='k', ls='--', alpha=0.4)
        axis.grid(axis='both', which='minor', color='k', ls='--', alpha=0.2)

        if show:
            plt.show()

        return


    def filter2(self, fq_band, **kwargs):
        """
        Apply a filter
        """

        assert isinstance(fq_band, (list, tuple, np.ndarray))

        sample_rate = self.stats.sampling_rate
        new_trace   = self.copy()
        data        = self.data - self.data.mean()
        zerophase   = kwargs.get("zerophase", True)
        corners     = kwargs.get("corners", 2)

        taper = kwargs.get("taper", False)
        if taper:
            taper_p = kwargs.get("taper_p", 0.05)
            data = data * osi.cosine_taper(len(data), p=taper_p)

        if fq_band[0] and fq_band[1]:
            # data = osf.bandpass(data, freqmin=fq_band[0], freqmax=fq_band[1], df=sample_rate, corners=corners, zerophase=zerophase)
            b, a = scipy.signal.butter(corners, fq_band, fs=sample_rate, btype='band')
            data = scipy.signal.filtfilt(b, a, data)

        if fq_band[0] and not fq_band[1]:
            data = osf.highpass(data, freq=fq_band[0], df=sample_rate, zerophase=zerophase)

        if fq_band[1] and not fq_band[0]:
            data = osf.lowpass(data, freq=fq_band, df=sample_rate, zerophase=zerophase)

        new_trace.data = data

        return new_trace


    def remove_response2(self, response, **kwargs):
        """
        Remove instrument response from the config resp file.
        The function is based on Obspy.
        """

        taper          = kwargs.get('taper', False)
        taper_fraction = kwargs.get('taper_fraction', 0.05)
        water_level    = kwargs.get('water_level', 60)
        output         = kwargs.get('output', 'VEL')
        pre_filt       = kwargs.get('pre_filt', [0.01, 0.005, 40, 45]) #

        new_trace = self.copy()
        data = new_trace.data.astype(np.float64)
        npts = len(data)
                      
        # time domain pre-processing
        data -= data.mean()
        if taper:
            data = data * osi.cosine_taper(npts, taper_fraction, sactaper=True, halfcosine=False)
    
        # Transform data to Frequency domain
        nfft = _npts2nfft(npts)
        data = np.fft.rfft(data, n=nfft)
        freq_response, freqs = response.get_evalresp_response(self.stats.delta, nfft, output)

        if pre_filt:
            freq_domain_taper = osi.cosine_sac_taper(freqs, flimit=pre_filt)
            data *= freq_domain_taper
        
        if water_level is None:
            freq_response[0] = 0.0
            freq_response[1:] = 1.0 / freq_response[1:]
        else:
            osi.invert_spectrum(freq_response, water_level)

        data = data * freq_response
        data[-1] = abs(data[-1]) + 0.0j
        data = np.fft.irfft(data)
        new_trace.data = data[0:npts]

        return new_trace


    def remove_factor(self, factor, disp):
        """
        Divide the data with a factor (typically, the instrument sensitivity) from the config seisvo file.
        disp :: return displacement [m] instead of velocity [m/s]
        """

        tr = self.copy()

        tr.data = tr.data / factor
            
        if disp:
            tr.data -= tr.data.mean()
            tr_disp = scipy.integrate.cumtrapz(tr.data, dx=tr.stats.delta, initial=0)
            tr_disp = scipy.signal.resample(tr_disp, tr.stats.npts)
            tr.data = tr_disp
        
        return tr


    def psd(self, fq_band=()):
        data = self.get_data(demean=True)
        fs = self.stats.sampling_rate
        ans = get_PSD(data, fs, fq_band=fq_band, NW=3.5, pad=1.0, full_return=False)
        return ans

class Stream2(Stream):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def __iter__(self):
        tr = [Trace2(self.traces[item]) for item in range(len(self))]
        return tr.__iter__()


    def __getitem__(self, item):
        return Trace2(self.traces[item])


    def get_stations_id(self):
        station_id = []
        for tr in self:
            station_id.append('.'.join(tr.id.split('.')[:-1]))
        return list(set(station_id))


    def select2(self, **kwargs):
        st = self.select(**kwargs)
        if st:
            return Stream2(st)
        else:
            return None


    def filter2(self, fq_band, **kwargs):
        st = Stream2()
        for tr in self:
            st.append(Trace2(tr).filter2(fq_band, **kwargs))
        return st


    def remove_response2(self, resp, **kwargs):
        st = Stream2()
        for trace in self:
            if isinstance(resp, dict):
                resp_dict = resp[trace.stats.channel]
            else:
                resp_dict = resp
            tr = Trace2(trace).remove_response2(resp_dict, **kwargs)
            st.append(tr)
        return st
    

    def remove_factor(self, factor, disp=False):
        st = Stream2()
        for trace in self:
            if isinstance(factor, dict):
                fval = factor[trace.stats.channel]
            else:
                fval = factor
            tr = Trace2(trace).remove_factor(fval, disp)
            st.append(tr)
        return st
    

    def get_bounds(self):
        starttime = max([tr.stats.starttime.datetime for tr in self])
        endtime   = min([tr.stats.endtime.datetime for tr in self])
        return (starttime, endtime)


    def to_array(self, sort=None, return_info=False, **trace_kwargs):
        """
        Return a np.ndarray with the traces of the stream

        sort must be a string or none and is only valid for sort multicomponent stream ZNE
        """

        channel_list = [tr.stats.channel for tr in self]
        if len(channel_list) > 3:
            sort = None

        if sort:
            assert len(sort) <= len(self)
            # check that all components are available
            channel_list = [chan for chan in channel_list if chan[-1] in sort]
            sort = ''.join([chan[-1] for chan in channel_list])
        
        
        # build an empty array
        if not trace_kwargs.get("starttime", None):
            starttime = max([trace.stats.starttime.datetime for trace in self])
            trace_kwargs["starttime"] = starttime
            

        if not trace_kwargs.get("endtime", None):
            endtime = min([trace.stats.endtime.datetime for trace in self])
            trace_kwargs["endtime"] = endtime
            

        sample_rate = trace_kwargs.get("sample_rate", self[0].stats.sampling_rate)
        npts = int((endtime - starttime).total_seconds()*sample_rate + 1)
        mdata = np.empty((len(self), npts))

        info = []
        for n, trace in enumerate(self):
            if sort:
                n = sort.index(trace.stats.channel[-1])
                info.append(trace.stats.channel)
            else:
                info.append(trace.stats.station + trace.stats.channel[-1])

            try:
                trd = trace.get_data(**trace_kwargs)
                if len(trd) < npts:
                    print(f" [warn] trace {trace.stats.channel} resample {len(trd)} --> {npts}")
                    trd = scipy.signal.resample(trd, npts)
            except Exception as e:
                print(f"  error reading trace :: {e}")
                return None
            
            mdata[n,:] = trd[:npts]
        
        if return_info:
            return mdata, info
        else:
            return mdata


    def time_polar(self, window, olap=0.9, **trace_kwargs):

        data = self.to_array(sort="ZNE", **trace_kwargs)
        fs   = trace_kwargs.get("sample_rate", self[0].stats.sampling_rate)

        if data.shape[0] == 3:
            ans = get_TimePolar(data, fs, window, olap)
        else:
            ans = None
            
        return ans