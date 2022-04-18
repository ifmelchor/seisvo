#!/usr/bin/python3
# seisvo
'''
 This module extend the obspy modules
'''

import os
import pickle
import numpy as np
from scipy import signal
from datetime import timedelta, datetime

# import matplotlib.pyplot as plt
# import matplotlib.ticker as mtick
# import matplotlib.dates as mdates

from obspy import Stream, Trace, read, UTCDateTime
from obspy.core.util import create_empty_data_chunk
from obspy.signal.invsim import cosine_taper
from seisvo.signal.spectrum import spectrogram


def read2(*args, **kwargs):
    st = read(*args, **kwargs)
    return Stream2(st)


class Trace2(Trace):
    def __init__(self, trace):
        super().__init__(data=trace.data, header=trace.stats)
        

    def get_data(self, starttime=None, endtime=None, detrend=True, fq_band=(), **kwargs):
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
        data = tr.data

        if detrend:
            data = signal.detrend(data)

        if fq_band != ():
            tr_filt = tr.filter2(fq_band, **kwargs)
            data = tr_filt.data

        return data


    def get_time(self, starttime=None, endtime=None):
        if not starttime:
            starttime = self.stats.starttime.datetime

        if not endtime:
            endtime = self.stats.endtime.datetime

        dt = self.stats.delta
        npts = int((endtime-starttime).total_seconds()*self.stats.sampling_rate)+1
        return [starttime + timedelta(seconds=dt*x) for x in range(npts)]


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


    def psd(self, starttime=None, endtime=None, fq_band=(), avg_step=None, axis=None,
        plot=True, log_scale='both', opt_return=False, return_fig=False, **kwargs):
        """
        This code get the power spectrum density.
        By default it return freq, psd arrays and prints the PSD plot.
        :param avg_step: step in minutes to average psd
        :type avg_step: int
        :param plot: boolean
        :param save: boolean, If True save the psd to png file as Station_Localization_Starttime.png,
        :param loglog: string. 'y' for yaxis, 'x' for xaxis and 'both' for both axis
        :param opt_return: boolean. If True, function also return the centroid and dominant of the PSD
        """

        from seisvo.signal.spectrum import power_density_spectrum

        if opt_return:
            from seisvo.signal.spectrum import get_centroid, get_dominant

        if starttime:
            start_time = UTCDateTime(starttime)
        else:
            start_time = self.stats.starttime
        if endtime:
            end_time = UTCDateTime(endtime)
        else:
            end_time = self.stats.endtime

        data = self.get_data(starttime=starttime, endtime=endtime, detrend=True, **kwargs)
        psd, freq = power_density_spectrum(data, self.stats.sampling_rate, fq_band=fq_band,
            avg_step=avg_step, **kwargs)

        if plot or return_fig or axis:

            if not axis:
                fig = plt.figure(figsize=(9,3))
                ax = fig.add_subplot(111)
            else:
                ax = axis
            
            if kwargs.get('db_scale', False):
                psd = 10*np.log10(psd)

            h = ax.plot(freq, psd, ls=kwargs.get('linestyle', '-'), color=kwargs.get('colorline', 'k'), label=kwargs.get('label', None))
            # ax.set_title('%s.%s.%s \n %s - %s' % (self.stats.station,
            #                                       self.stats.location,
            #                                       self.stats.channel,
            #                                       start_time.datetime,
            #                                       end_time.datetime))
            if kwargs.get('db_scale', False):
                ax.set_ylabel('PSD\n'+r'[m$^2$s$^{-4}$/Hz dB]')
            else:
                ax.set_ylabel('PSD\n'+r'[m$^2$s$^{-2}$/Hz]')

            ax.set_xlabel('Freq. [Hz]')
            ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(6))

            if log_scale in ('x', 'y', 'both'):
                if log_scale == 'x':
                    ax.set_xscale('log')

                if log_scale == 'y':
                    ax.set_yscale('log')

                if log_scale == 'both':
                    ax.set_yscale('log')
                    ax.set_xscale('log')

            if plot:
                plt.show()

            if axis:
                return h
            
            else:
                return fig

        list_to_return = [psd, freq]

        if opt_return:
            dominant = get_dominant(freq, psd)
            centroid = get_centroid(freq, psd)
            list_to_return += [(dominant, centroid)]

        return list_to_return


    def rms(self, starttime=None, endtime=None, fq_band=(), avg_step=None, **kwargs):
        """
        Get the RMS in the time and frequency domain
        :param step: in min
        :type step: int
        """

        from seisvo.signal.spectrum import power_density_spectrum

        if starttime:
            start_time = UTCDateTime(starttime)
        else:
            start_time = self.stats.starttime

        if endtime:
            end_time = UTCDateTime(endtime)
        else:
            end_time = self.stats.endtime

        data = self.get_data(starttime=starttime, endtime=endtime, detrend=True, fq_band=fq_band)
        psd, freq = power_density_spectrum(data, self.stats.sampling_rate, fq_band=fq_band,
            avg_step=avg_step, **kwargs)

        return sum(psd), sum(np.power(data,2))*self.stats.delta


    def filter2(self, fq_band, **kwargs):
        """
        Apply a filter
        """

        import obspy.signal.filter as osf
        from obspy.signal.invsim import cosine_taper

        sample_rate = self.stats.sampling_rate
        new_trace = self.copy()

        # detrend and taper
        data = self.data - self.data.mean()

        zerophase = kwargs.get("zerophase", True)
        taper = kwargs.get("taper", True)
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
    

    def get(self, *args, **kwargs):
        st = self.select( *args, **kwargs)
        if st:
            return Trace2(st[0])
        else:
            print(' warn: empty trace')
            return None


    def cross_spectrum_matrix(self, starttime=None, endtime=None, station=None, fq_band=(), avg_step=False, polar_return=False, matrix_return=False, threshold=0.75, **kwargs):
        """
        This functions return the cross_spectral denisty matrix using the multitaper algorithm
        :param station:
        :param freq_bands:
        :param time_width:
        :param overlap:
        :param time_bandwidth:
        :return:
        """

        from seisvo.signal.spectrum import cross_spectrum_matrix
        from seisvo.signal.polarization import polarization_content

        Z = self.get_component('Z', station=station)
        E = self.get_component('E', station=station)
        N = self.get_component('N', station=station)

        zdata = Z.get_data(starttime=starttime, endtime=endtime, detrend=True)
        edata = E.get_data(starttime=starttime, endtime=endtime, detrend=True)
        ndata = N.get_data(starttime=starttime, endtime=endtime, detrend=True)

        if len(zdata) == len(edata) == len(ndata):
            if Z.stats.sampling_rate == E.stats.sampling_rate == N.stats.sampling_rate:
                spec_mat = cross_spectrum_matrix(zdata, ndata, edata, Z.stats.sampling_rate, avg_step=avg_step, fq_band=fq_band, **kwargs)
                if polar_return:
                    polar_ans = polarization_content(spec_mat[0], matrix_return=matrix_return, threshold=threshold)

                    # get max degree
                    max_degree = polar_ans[0].max()

                    # get frequency at max degree
                    fq_polar = spec_mat[1][np.argmax(polar_ans[0])]

                    return polar_ans + [max_degree, fq_polar]
                else:
                    return spec_mat
            else:
                print('Sampling rate of all the traces must be the same')
                return
        else:
            print('Length of all the traces must be the same')
            return


    def filter2(self, fq_band, **kwargs):
        st = Stream2()
        for tr in self:
            tr.detrend()
            st.append(Trace2(tr).filter2(fq_band, **kwargs))
        return st


    def remove_response2(self, resp_dict, **kwargs):
        st = Stream2()
        for trace in self:
            tr = Trace2(trace).remove_response2(resp_dict, **kwargs)
            st.append(tr)
        return st


    def polardegree(self, starttime=None, endtime=None, station=None, fq_band=[], avg_step=None,
        olap=0.75, matrix_return=False, opt_return=False, **kwargs):
        """
        This function compute the polarization degree

        :param starttime: datetime (optional)
        :param endtime: datetime (optional)
        :param avg_step: int in minutes (optional)
        :param station: optional. string
        :return:
        """

        from seisvo.signal.polarization import polar_degree

        Z = self.get_component('Z', station=station)
        E = self.get_component('E', station=station)
        N = self.get_component('N', station=station)

        zdata = Z.get_data(starttime=starttime, endtime=endtime, detrend=True)
        edata = E.get_data(starttime=starttime, endtime=endtime, detrend=True)
        ndata = N.get_data(starttime=starttime, endtime=endtime, detrend=True)

        sample_rate = Z.stats.sampling_rate

        if not Z.stats.sampling_rate == E.stats.sampling_rate == N.stats.sampling_rate:
            raise ValueError('all traces must have same sampling rate')

        if not Z.stats.npts == E.stats.npts == N.stats.npts:
            raise ValueError('all traces must have same length')

        if matrix_return:
            freq, polar_dgr, azm, dip, rect = polar_degree(zdata, ndata, edata, sample_rate, per_lap=olap, avg_step=avg_step,
                fq_band=fq_band,  matrix_return=True, **kwargs)
        
        else:
            freq, polar_dgr = polar_degree(zdata, ndata, edata, sample_rate, per_lap=olap, avg_step=avg_step,
                fq_band=fq_band, **kwargs)

        for_return = [polar_dgr, freq]
        
        if opt_return:
            degree_dominant = polar_dgr.max()
            fq_dominant = freq[np.argmax(polar_dgr)]
            for_return += [(degree_dominant, fq_dominant)]

        if matrix_return:
            for_return += [azm, dip, rect]
                
        return for_return


    def polargram(self, window_length=None, starttime=None, endtime=None, station=None, fq_band=[],
        olap=0.75, matrix_return=False, axis=None, **kwargs):

        """
        This function plot the polarization degree

        :param starttime: datetime
        :param endtime: datetime
        :param time_width: seconds
        :param station: optional. string
        :param kwargs: cdm_overlap, cdm_width (sec)
        :return:
        """

        from seisvo.signal.polarization import polarization_content, polar_degree
        from seisvo.utils.plotting import plot_gram
        from seisvo.signal import freq_bins, _nearest_pow_2

        Z = self.get_component('Z', station=station)
        E = self.get_component('E', station=station)
        N = self.get_component('N', station=station)

        zdata = Z.get_data(starttime=starttime, endtime=endtime, detrend=True)
        edata = E.get_data(starttime=starttime, endtime=endtime, detrend=True)
        ndata = N.get_data(starttime=starttime, endtime=endtime, detrend=True)

        sample_rate = Z.stats.sampling_rate
        npts = Z.stats.npts

        if not Z.stats.sampling_rate == E.stats.sampling_rate == N.stats.sampling_rate:
            raise ValueError('all traces must have same sampling rate')

        if not Z.stats.npts == E.stats.npts == N.stats.npts:
            raise ValueError('all traces must have same length')

        if not window_length:
            window_length = int(sample_rate)//50

        # defining nfft (windows of the FFT transfrom)
        nfft = int(_nearest_pow_2(window_length * sample_rate))
        if nfft > npts:
            nfft = int(_nearest_pow_2(npts / 8.0))

        # pad with zeros to smooth the spectrogram
        mult = int(_nearest_pow_2(8.0))
        mult = mult * nfft

        # overlap
        nlap = int(mult * float(olap))
        step = mult - nlap
        
        shape = zdata.shape[:-1]+((zdata.shape[-1]-nlap)//step, mult)
        strides = zdata.strides[:-1]+(step*zdata.strides[-1], zdata.strides[-1])

        zresult = np.lib.stride_tricks.as_strided(zdata, shape=shape, strides=strides)
        nresult = np.lib.stride_tricks.as_strided(ndata, shape=shape, strides=strides)
        eresult = np.lib.stride_tricks.as_strided(edata, shape=shape, strides=strides)

        # freqshape
        freq, _, _ = freq_bins(mult, sample_rate, fq_band=fq_band, nfft='uppest', get_freq=True)
        polar_dgr = np.empty((freq.shape[0], shape[0]))
        
        if matrix_return:
            azm = np.empty((freq.shape[0], shape[0]))
            dip = np.empty((freq.shape[0], shape[0]))
            rect = np.empty((freq.shape[0], shape[0]))

            for n, (zi, ni, ei) in enumerate(zip(zresult, nresult, eresult)):
                freq, polar_dgr[:,n], azm[:,n], dip[:,n], rect[:,n] = polar_degree(zi, ni, ei, sample_rate, avg_step=None, fq_band=fq_band, matrix_return=True, **kwargs)
        
            for_return = [freq, polar_dgr, rect, azm, dip]
        
        else:
            for n, (zi, ni, ei) in enumerate(zip(zresult, nresult, eresult)):
                freq, polar_dgr[:,n]  = polar_degree(zi, ni, ei, sample_rate, avg_step=None, fq_band=fq_band, **kwargs)
            
            for_return = [freq, polar_dgr]
        
        return for_return

        