#!/usr/bin/python3
# coding=utf-8

import os
import utm
import scipy
import numpy as np
import datetime as dt
from obspy import UTCDateTime
from obspy.core.inventory.response import Response

from .utils import nCPU
from .obspyext import Stream2
from .lte.base import _new_LTE
from .signal import SSteps, get_freq, get_Polar, get_PSD, get_PDF
from .plotting import plotPDF


class Station(object):
    def __init__(self, StaFile):
        self.stats = StaFile


    def __str__(self):
        return self.stats.__str__()


    def is_component(self, component):
        """ 
        Check if the component exist
        """ 

        ans = False
        for true_chan in self.stats.channels:
            if true_chan[-1] == component:
                ans = True

        return ans


    def is_three_component(self):
        return all(list(map(self.is_component, ["Z", "N", "E"])))


    def get_chan(self, component):
        """ 
        Get all channels with the same component
        """

        if isinstance(component, str):
            component = [component]

        ans = []
        for true_chan in self.stats.channels:
            if true_chan[-1] in component:
                ans.append(true_chan)

        if len(ans) == 1:
            return ans[0]
        
        return ans


    def remove_response(self, stream, **kwargs):
        """
        Remove response of a stream2 using data loaded in RESP file of SEISVO.

        """

        resp = self.stats.get_response()

        if not isinstance(stream, Stream2):
            stream = Stream2(stream)

        if resp:
            stream_resp = stream.remove_response2(resp, **kwargs)

        else:
            print(" >>> warn: no response file info")
            stream_resp = stream
        
        return stream_resp
    

    def remove_factor(self, stream, disp=False):
        """
        Remove sensitivity of a stream2 using data loaded in RESP file of SEISVO.

        """

        factor = self.stats.get_factor()

        if not isinstance(stream, Stream2):
            stream = Stream2(stream)

        if factor:
            stream_resp = stream.remove_factor(factor, disp=disp)

        else:
            print(" >>> warn: no response file info")
            stream_resp = stream
        
        return stream_resp


    def get_stream(self, starttime, endtime, channel=None, **kwargs):
        """
        Get stream from station object
        :param starttime: datetime
        :param endtime: datetime
        :param chan: list or string e.j.:'SHZ' or 'HHN'. optional
        :param kwargs: remove_sensitivity, remove_response, prefilt, fill_values
        :return: stream2 object
        """

        assert starttime < endtime
        assert starttime >= self.stats.starttime
        assert endtime <= self.stats.endtime

        channel_list = self.stats.check_channel(channel)
        assert channel_list

        # most of MSEED files do not start in 00:00:00 and (finish in 23:59:59), 
        # this is why we need to define time delta to overcome this lack.
        offtime = kwargs.get('offtime', 10)
        time_delta = dt.timedelta(minutes=offtime)
        t0, tf = self.stats.get_offtimes(starttime, endtime, time_delta)
        day_diff = (tf.date() - t0.date()).days + 1
        date_list = [t0.date() + dt.timedelta(days=i) for i in range(day_diff)]

        stream, sample_rate_list = Stream2(), []
        for day in date_list:
            for channel in channel_list:
                st_day = self.stats.__read_file__(channel, date=day, stream=True, starttime=t0, endtime=tf)
                if st_day:
                    stream += st_day
                    sample_rate_list.append(st_day[0].stats.sampling_rate)
                else:
                    print(' warn [STA.ID: %s.%s]. No data for %s' % (self.stats.code, self.stats.location, day.strftime('%d %b %Y')))

        # if different traces with different sample rates are in date, discard stream
        if len(list(set(sample_rate_list))) > 1:
            print('error: stream with mixed sampling rates. Revise data.')
            return None

        # merge traces
        fill_value = kwargs.get('fill_value', None)
        method     = kwargs.get('method', 0)
        stream.merge(method=method, fill_value=fill_value)
        st = stream.slice(UTCDateTime(starttime), UTCDateTime(endtime))

        if st:
            if kwargs.get('remove_response', False):
                rrkwargs = kwargs.get('rrkwargs', {})
                st = self.remove_response(st, **rrkwargs)
                kwargs['remove_sensitivity'] = False
            
            if kwargs.get('remove_sensitivity', False):
                disp = kwargs.get('disp', False)
                st = self.remove_factor(st, disp=disp)
            
            sample_rate = kwargs.get('sample_rate', None)
            if sample_rate and sample_rate < sample_rate_list[0]:
                st = Stream2(st.resample(sample_rate))

                if st[0].stats.sampling_rate != sample_rate:
                    print("\n  [warn] stream data could not be resampled. Check your data.")
                    return None
        
            prefilt = kwargs.get('prefilt', [])

            if prefilt:
                st = st.filter2(fq_band=prefilt)
            
        return st


    def get_psd(self, starttime, endtime, channel, window=0, olap=0.75,\
        fq_band=(0.1,15), return_pdf=False, plot_pdf=False, **st_kwargs):
        """
        Compute the PSD between start and end time for a specific channel.
        if window [float, seconds] > 0, then a moving average is applied with overlap [olap] between [0,1)

        if return_pdf is True, return PDF of all window segments of the trace 
        """

        assert self.stats.check_channel(channel)
        stream = self.get_stream(starttime, endtime, channel=channel, **st_kwargs)
        
        assert stream.get_bounds() == (starttime, endtime)
        data = stream.to_array(detrend=True)
        fs = stream[0].stats.sampling_rate

        full_return, lwin = (False, None)
        if window > 0:
            lwin = window*fs
            if return_pdf:
                full_return = True
        
        ans = get_PSD(data[0,:], fs, lwin=lwin, olap=olap,\
            fq_band=fq_band, NW=3.5, pad=1.0, full_return=full_return)

        if return_pdf:
            array = 10*np.log(ans[0])
            start = np.floor(array.min())
            stop  = np.ceil(array.max())
            space = np.linspace(start, stop, num=1000).reshape(-1,1)
            pdf = get_PDF(array, space)
            ans = {
                "psd":ans[0],
                "pdf":pdf,
                "y":space.reshape(-1,),
                "freq":ans[1]
                }

            if plot_pdf:
                plotPDF(ans["pdf"], ans["y"], ans["freq"])

        return ans


    def get_polarization(self, starttime, endtime, window=0, olap=0.75, fq_band=(1.,5.),\
        full_analysis=True, return_pdf=False, azimuth_ambiguity=True, **st_kwargs):
        """
        Compute the polarization analysis between start and end time.
        if window [float, seconds] > 0, then a moving average is applied with overlap [olap] between [0,1)

        if return_pdf is True, return PDF of all window segments of the trace
        """

        assert self.is_three_component()
        channel_list = self.get_chan(["Z","N","E"])
        stream = self.get_stream(starttime, endtime, channel=channel_list, **st_kwargs)
        assert stream.get_bounds() == (starttime, endtime)
        
        data = stream.to_array(detrend=True, sort="ZNE")
        fs = stream[0].stats.sampling_rate

        if window > 0:
            lwin = window*fs
        else:
            lwin = None
            olap = None

        freq, polar = get_Polar(data, fs, lwin=lwin, olap=olap,\
            fq_band=fq_band, return_all=return_pdf, full_return=full_analysis)

        if full_analysis and azimuth_ambiguity:
            ts = polar["azimuth"]
            polar["azimuth"] = np.where(ts>180, ts-180, ts)

        if return_pdf:
            if full_analysis:
                for key, array in polar.items():
                    if key in ("degree", "rect"):
                        start, stop = (0, 1)
                    elif key == "azimuth":
                        if azimuth_ambiguity:
                            start, stop = (0, 180)
                        else:
                            start, stop = (0, 360)
                    elif key == "elev":
                        start, stop = (0, 90)
                    else:
                        start = np.floor(array.min())
                        stop  = np.ceil(array.max())
                    space = np.linspace(start, stop, num=1000).reshape(-1,1)
                    pdf = get_PDF(array, space)
                    polar[key]["pdf"] = {"pdf":pdf,"y":space.reshape(-1,)}
            
            else:
                start, stop = (0, 1)
                space = np.linspace(start, stop, num=1000).reshape(-1,1)
                pdf = get_PDF(array, space)
                polar = {"degree":polar, "pdf":pdf,"y":space.reshape(-1,)}
        
        return freq, polar


    def lte(self, starttime, endtime, window, subwindow, interval=None,\
        channel=None, window_olap=0, subwindow_olap=0, **kwargs):
        """ Compute (station) LTE file

        Parameters
        ----------
        starttime : datetime
        
        endtime : datetime
        
        window : float [sec]
            length of the time window (in min) to reduce
        
        subwindow : float [sec]
            length of the time window (in sec) for moving average over the window. 
            If ``subwindow=0`` no moving average is applied.
        
        interval : int [min]
            length of time window (in h) to store data in memory. By default is 1 hour.
        
        channel : _type_, optional
            _description_, by default None
        
        subwindow_olap : float
            overlap percent for moving average over the interval, by default 0.75

        Returns
        -------
        LTE
            LTE object
        """

        assert starttime < endtime
        assert starttime >= self.stats.starttime
        assert endtime <= self.stats.endtime

        # load parameters
        channel     = tuple(self.stats.check_channel(channel))
        sample_rate = kwargs.get('sample_rate', 40)
        njobs       = kwargs.get('njobs', 2)
        pad         = kwargs.get("pad",1.0)
        fq_band     = kwargs.get('fq_band', (0.5, 10))
        file_name   = kwargs.get("file_name", None)
        out_dir     = kwargs.get("out_dir", './')

        # defining base params
        ltebase = dict(
            id              = self.stats.id,
            type            = "station",
            channel         = channel,
            starttime       = starttime.strftime('%Y-%m-%d %H:%M:%S.%f'),
            endtime         = endtime.strftime('%Y-%m-%d %H:%M:%S.%f'),
            interval        = interval,
            window          = window,
            window_olap     = window_olap,
            subwindow       = subwindow,
            subwindow_olap  = subwindow_olap,
            sample_rate     = int(sample_rate),
            pad             = pad,
            fq_band         = fq_band,
            polar           = kwargs.get('polar', False),
            opt_params      = kwargs.get('opt_params', False),
            rm_sens         = kwargs.get('rm_sens', False),
            opt_twin        = kwargs.get('opt_twin', 150),
            opt_th          = kwargs.get('opt_th', 5.0),
            pe_tau          = int(kwargs.get("pe_tau", 2)),
            pe_order        = int(kwargs.get("pe_order", 7)),
            time_bandwidth  = kwargs.get('time_bandwidth', 3.5)
        )

        ss = SSteps(starttime, endtime, window, interval=interval, win_olap=window_olap, subwindow=subwindow, subw_olap=subwindow_olap)
        ssdict = ss.to_dict()

        ltebase["lwin"] = int(ss.window*sample_rate)
        ltebase["nwin"] = ssdict["nwin"]
        ltebase["last_nwin"] = ssdict["last_nwin"]
        ltebase["wadv"] = ssdict["wadv"]
        ltebase["nro_intervals"] = ssdict["nro_intervals"]
        ltebase["nro_time_bins"] = ssdict["total_nwin"]
        ltebase["int_extra_sec"] = ssdict["int_extra_sec"]

        if subwindow > 0:
            ltebase["lswin"] = int(ss.subwindow*sample_rate)
            ltebase["nswin"] = ssdict["nswin"]
            ltebase["swadv"] = ssdict["swadv"]
            nfs, _ = get_freq(ltebase["lswin"], sample_rate, fq_band=fq_band, pad=pad)
        else:
            ltebase["lswin"] = ltebase["nswin"] = ltebase["swadv"] = None
            nfs, _ = get_freq(ltebase["lwin"], sample_rate, fq_band=fq_band, pad=pad)

        ltebase["nro_freq_bins"] = len(nfs)

        if njobs >= nCPU or njobs == -1:
            njobs = nCPU - 2
        
        if njobs > ssdict["nro_intervals"]:
            njobs = ssdict["nro_intervals"]
            print(f"warn  ::  njobs set to {njobs}")
            
        # create hdf5 file and process data
        if not file_name:
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (self.stats.id, starttime.year, starttime.timetuple().tm_yday, endtime.year, endtime.timetuple().tm_yday, window)
        
        file_name_full = os.path.join(out_dir, file_name)
        if file_name_full.split('.')[-1] != 'lte':
            file_name_full += '.lte'

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)

        lte = _new_LTE(self, file_name_full, ltebase, njobs)
        
        return lte


    def plot_control(self, chan, startdate=None, enddate=None):
        """
        Plot information of your mseed files
        :param chan: channel e.j. 'SHZ'
        :param starttime: datetime object, optional
        :param endtime: datetime object. optinal
        """

        if chan not in self.stats.channels:
            raise TypeError('Channel not loaded')

        if not startdate:
            startdate = self.stats.starttime

        if not enddate:
            enddate = self.stats.endtime

        day_diff = (enddate - startdate).days
        date_list = [startdate + dt.timedelta(days=i) for i in range(day_diff+1)]

        npts = []
        sample_rate = []
        nro_traces = []
        filesize = []

        for i, item in enumerate(date_list):
            print ('  reading data... (%d%%)' % (100*i/len(date_list)), end='\r')
            st = self.stats.__read_file__(chan, date=item, stream=True, headonly=True)
            
            if st:
                npts += [sum([tr.stats.npts for tr in st])]
                sample_rate += [sum([tr.stats.sampling_rate for tr in st])/len(st)]
                nro_traces += [len(st)]
                filesize += [sum([tr.stats.filesize for tr in st])]
            
            else:
                npts += [None]
                sample_rate += [None]
                nro_traces += [None]
                filesize += [None]

        print('\t\t')

        title = '%s \n %s -- %s' % (self.stats.id,
            startdate.strftime('%Y.%m.%d'),
            enddate.strftime('%Y.%m.%d'))
        
        pplot_control(title, date_list, nro_traces, sample_rate, npts, filesize)


    # def plot(self, starttime, sde, channel='all', delta=30, app=False, **kwargs):
    #     """ Plot seismogram of the station in a simple GUI. 

    #     Args:
    #         channel (str or list): string or list of string with the channels to plot
    #         starttime (datetime): start time to plot
    #         delta (int, optional): in minutes. Defaults to 15.
    #         return_fig (bool, optional): Return fig and axes objects. Defaults to False.
    #     """

    #     from seisvo.plotting.gui.gstation import plot_station_gui

    #     window = plot_station_gui(self, starttime, sde, channel=channel, delta=delta, app=app, **kwargs)
        
    #     if app:
            
    #         return window




    
