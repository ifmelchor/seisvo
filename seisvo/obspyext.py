#!/usr/bin/python3
# seisvo
'''
 This module extend the obspy modules
'''

import scipy
import numpy as np
import datetime as dt
import obspy.signal.invsim as osi
from obspy import Stream, Trace, read, UTCDateTime
from obspy.signal.util import _npts2nfft
from obspy.signal import filter as osf
from .signal import get_PSD


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
        return [starttime + dt.timedelta(seconds=delta*x) for x in range(npts)]


    def plot_specgram(self, window_length=None, starttime=None, endtime=None, axes=None, fq_band=(), per_lap=0.75, returnfig=False, **kwargs):
        """
        This code plots the spectrogram

        :param data: array_like, data series
        :param window_length: int, the length of the window to compute the fft. By default samp_rate / 2
        :param rel_norm: True for relative normalization. By default, 'abs' normalization
        :param axes: axes object to plot the spectrogram
        :param cmap: 'jet' by default
        """

        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        from seisvo.signal.spectrum import spectrogram

        v_min = kwargs.get("v_min", None)
        v_max = kwargs.get("v_max", None)
        a_min = kwargs.get("a_min", None)
        a_max = kwargs.get("a_max", None)
        # cmap = kwargs.get("cmap", 'Spectral_r')
        # rel_norm = kwargs.get("rel_norm", False)
        # interpolation = kwargs.get("interpolation", 'gaussian')
        title = kwargs.get('title', None)
        major_ticks = kwargs.get('major_ticks', None)
        major_fmt = kwargs.get('major_fmt', None)
        minor_ticks = kwargs.get('minor_ticks', None)
        # minor_fmt = kwargs.get('minor_fmt', None)
        major_xtickslabels = kwargs.get('major_xtickslabels', None)
        xlabel = kwargs.get('xlabel', None)
        dateticks = kwargs.get("dateticks", 'byminute')

        if not starttime:
            starttime = self.stats.starttime.datetime

        if not axes:
            fig, axs = plt.subplots(2,2, gridspec_kw={"width_ratios":[1, 0.02],"wspace":0.02},
                constrained_layout=True, figsize=(8, 4))
            ax_trace = axs[0,0]
            ax_spect = axs[1,0]
            ax_spect_color = axs[1,1]
            axs[0,1].axes.get_xaxis().set_visible(False)
            axs[0,1].axes.get_yaxis().set_visible(False)
            axs[0,1].axis('off')

            data = self.get_data(starttime=starttime, endtime=endtime, detrend=True, fq_band=fq_band)
            time = self.get_time(starttime=starttime, endtime=endtime)

            if not starttime and not endtime:
                starttime = self.stats.starttime.datetime
                endtime = self.stats.endtime.datetime

            if not major_ticks:
                if dateticks == 'byminute':
                    major_ticks = mdates.MinuteLocator(byminute=range(0,60,20))
                    major_fmt = mdates.DateFormatter('%H:%M:%S')
                    minor_ticks = mdates.MinuteLocator(byminute=range(0,60,10))

                if dateticks == 'byhour':
                    major_ticks = mdates.HourLocator(byhour=range(0,24,2))
                    major_fmt = mdates.DateFormatter('%H:%M')
                    minor_ticks = mdates.HourLocator(byhour=range(0,24,1))

            ax_trace.plot(time, data, 'k')
            ax_trace.set_ylabel('cnts')
            ax_trace.set_xlim(time[0], time[-1])

            if a_min and a_max:
                ax_trace.set_ylim(a_min, a_max)

            if not title:
                title = '%s.%s.%s \n %s - %s' % (self.stats.station,
                                                 self.stats.location,
                                                 self.stats.channel,
                                                 starttime,
                                                 endtime)
            ax_trace.set_title(title)

            if minor_ticks:
                ax_trace.xaxis.set_minor_locator(minor_ticks)

            if major_ticks:
                ax_trace.xaxis.set_major_locator(major_ticks)
                ax_trace.xaxis.set_major_formatter(major_fmt)

            if major_xtickslabels:
                ax_trace.set_xticklabels(major_xtickslabels)

        else:
            ax_spect = axes

        data = self.get_data(starttime=starttime, endtime=endtime, detrend=True, taper=True)

        im, (v_min, v_max) = spectrogram(data, self.stats.sampling_rate, ax_spect, per_lap=per_lap,
            window_length=window_length, fq_band=fq_band, date_list=self.get_time(), **kwargs)

        if not xlabel:
            xlabel = 'UTC Time ['+starttime.strftime('%Y-%m-%d')+']'

        ax_spect.set_xlabel(xlabel)

        if minor_ticks:
            ax_spect.xaxis.set_minor_locator(minor_ticks)

        if major_ticks:
            ax_spect.xaxis.set_major_locator(major_ticks)
            ax_spect.xaxis.set_major_formatter(major_fmt)

        if major_xtickslabels:
            ax_spect.set_xticklabels(major_xtickslabels)

        #ax_trace.set_xticklabels([' ']*len(ax_trace.get_xticks()))

        if not axes:
            fig.colorbar(im, cax=ax_spect_color, label='PSD\n'+r'dB[cnts$^2$/Hz]',
                orientation='vertical')
            fig.align_labels()

            if returnfig:
                return fig
            
            else:
                plt.show()

        return im, (v_min, v_max)


    def filter2(self, fq_band, **kwargs):
        """
        Apply a filter
        """

        assert isinstance(fq_band, (list, tuple, np.ndarray))

        sample_rate = self.stats.sampling_rate
        new_trace = self.copy()
        data = self.data - self.data.mean()

        zerophase = kwargs.get("zerophase", True)
        corners = kwargs.get("corners", 2)

        taper = kwargs.get("taper", False)
        if taper:
            taper_p = kwargs.get("taper_p", 0.05)
            data = data * osi.cosine_taper(len(data), p=taper_p)

        if fq_band[0] and fq_band[1]:
            data = osf.bandpass(data, freqmin=fq_band[0], freqmax=fq_band[1], df=sample_rate, corners=corners, zerophase=zerophase)

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

        resp = response.load()
        if resp:                
            # time domain pre-processing
            data -= data.mean()
            if taper:
                data = data * osi.cosine_taper(npts, taper_fraction, 
                        sactaper=True, halfcosine=False)
        
            # Transform data to Frequency domain
            nfft = _npts2nfft(npts)
            data = np.fft.rfft(data, n=nfft)
            freq_response, freqs = resp.get_evalresp_response(self.stats.delta, nfft, output)

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

        else:
            if response.factor:
                new_trace.data = data*response.factor

            else:
                print(' [remove_response] No resp info loadded.')

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


class Stream2(Stream):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def __iter__(self):
        tr = [Trace2(self.traces[item]) for item in range(len(self))]
        return tr.__iter__()


    def __getitem__(self, item):
        return Trace2(self.traces[item])


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
