#!/usr/bin/python3
# seisvo
'''
 This module extend the obspy modules
'''

import numpy as np
from scipy import signal
import datetime as dt

from obspy import Stream, Trace, read, UTCDateTime
from obspy.signal import filter as osf

from seisvo.signal.spectrum import power_density_spectrum, cosine_taper

def read2(*args, **kwargs):
    st = read(*args, **kwargs)
    return Stream2(st)


class Trace2(Trace):
    def __init__(self, trace):
        super().__init__(data=trace.data, header=trace.stats)
    

    def get_data(self, starttime=None, endtime=None, detrend=True, fq_band=(), abs=False, sample_rate=None, **kwargs):
        """
        This code returns a numpy array of the data
        """

        if starttime:
            st = UTCDateTime(starttime)
        else:
            st = self.stats.starttime

        if endtime:
            et = UTCDateTime(endtime)
        else:
            et = self.stats.endtime

        tr = self.slice(st, et)

        if sample_rate:
            tr = self.resample(sample_rate)
        
        data = tr.data

        if detrend:
            data = signal.detrend(data)

        if list(fq_band):
            tr_filt = tr.filter2(fq_band, **kwargs)
            data = tr_filt.data
        
        if abs:
            data = np.abs(data)
        
        return data


    def get_time(self, starttime=None, endtime=None):
        if not starttime:
            starttime = self.stats.starttime.datetime

        if not endtime:
            endtime = self.stats.endtime.datetime

        dt = self.stats.delta
        npts = int((endtime-starttime).total_seconds()*self.stats.sampling_rate)+1
        return [starttime + dt.timedelta(seconds=dt*x) for x in range(npts)]


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


    def psd(self, starttime=None, endtime=None, fq_band=(), mov_avg_step=None, olap_step=0.0, drm_params=False, **kwargs):
        """Compute PSD using multitaper algorithm

        Parameters
        ----------
        starttime : _type_, optional
            _description_, by default None
        endtime : _type_, optional
            _description_, by default None
        fq_band : tuple, optional
            _description_, by default ()
        mov_avg_step : _type_, optional
            _description_, by default None
        olap_step : int, optional
            _description_, by default 0

        Returns
        -------
        _type_
            _description_
        """

        if not starttime:
            starttime = self.stats.starttime

        if not endtime:
            endtime = self.stats.endtime

        if starttime < self.stats.starttime:
            raise ValueError(' no data for this dates')
        
        if endtime > self.stats.endtime:
            raise ValueError(' no data for this dates')

    
        if isinstance(starttime, dt.datetime):
            starttime = UTCDateTime(starttime)
            
        if isinstance(endtime, dt.datetime):
            endtime = UTCDateTime(endtime)
            
        if mov_avg_step:
            interval_min = int(endtime - starttime)

            if interval_min < mov_avg_step:
                raise ValueError ('interval is less than mov_avg_step!')
        
            if olap_step > 1:
                raise ValueError ('overlap percent is greater than 1')

        data = self.get_data(starttime=starttime, endtime=endtime, detrend=True)
        ans = power_density_spectrum(data, self.stats.sampling_rate, fq_band=fq_band,
            avg_step=mov_avg_step, olap=olap_step, drm_params=drm_params, **kwargs)

        return ans


    def filter2(self, fq_band, **kwargs):
        """
        Apply a filter
        """

        sample_rate = self.stats.sampling_rate
        new_trace = self.copy()

        # detrend and taper
        data = self.data - self.data.mean()

        zerophase = kwargs.get("zerophase", True)
        taper = kwargs.get("taper", False)
        taper_p = kwargs.get("taper_p", 0.05)
        bandstop = kwargs.get("bandstop", False)

        if taper:
            data = data * cosine_taper(len(data), p=taper_p)

        if fq_band[0] and fq_band[1]:
            if bandstop:
                data = osf.bandstop(data, freqmin=fq_band[0], freqmax=fq_band[1], df=sample_rate, zerophase=zerophase)
            else:
                data = osf.bandpass(data, freqmin=fq_band[0], freqmax=fq_band[1], df=sample_rate, zerophase=zerophase)

        elif fq_band[0] and not fq_band[1]:
            data = osf.highpass(data, freq=fq_band[0], df=sample_rate, zerophase=zerophase)

        elif fq_band[1] and not fq_band[0]:
            data = osf.lowpass(data, freq=fq_band, df=sample_rate, zerophase=zerophase)

        else:
            raise TypeError(' fq_band must be a list')

        new_trace.data = data
        return new_trace


    def remove_response2(self, response, **kwargs):
        """
        Remove instrument response from the config resp file.
        The function is based on Obspy.
        """

        import obspy.signal.invsim as osi
        from obspy.signal.util import _npts2nfft

        taper = kwargs.get('taper', False)
        taper_fraction = kwargs.get('taper_fraction', 0.05)
        water_level = kwargs.get('water_level', 60)
        output = kwargs.get('output', 'VEL')
        pre_filt = kwargs.get('pre_filt', [0.01, 0.005, 40, 45]) #

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


class Stream2(Stream):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def __getitem__(self, item):
        return Trace2(self.traces[item])


    def get_component(self, component, station=None, loc=None):
        tr = self.select(station=station, location=loc, component=component)
        if tr:
            return Trace2(tr[0])
        else:
            raise ValueError('Component do not exist')


    def filter2(self, fq_band, **kwargs):
        st = Stream2()
        for tr in self:
            st.append(Trace2(tr).filter2(fq_band, **kwargs))
        return st


    def remove_response2(self, resp_dict, **kwargs):
        st = Stream2()
        for trace in self:
            tr = Trace2(trace).remove_response2(resp_dict, **kwargs)
            st.append(tr)
        return st