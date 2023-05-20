#!/usr/bin/python3
# coding=utf-8

import os
import utm
import scipy
import pickle
import numpy as np
import datetime as dt
import multiprocessing as mp

from glob import glob
from obspy.core.inventory.response import Response

from seisvo import seisvo_paths
from .obspyext import UTCDateTime, read2, Stream2
from .lte import StationLTE
from .signal import SSteps, get_freq
from .signal.polarization import PolarAnalysis
from .plotting import pplot_control


class Station(object):
    def __init__(self, StaFile):
        self.stats = StaFile
        # self.resp_ = get_respfile(self.stats.net, self.stats.code, self.stats.loc)
        self.__set_dates__()


    def __str__(self):
        return self.stats.__str__()


    def __get_offtimes__(self, starttime, endtime, offtime):
        
        if starttime - offtime > self.starttime:
            start = starttime - offtime
        else:
            start = starttime

        if endtime + offtime < self.endtime:
            end = endtime + offtime
        else:
            end = endtime
        
        return (start, end)


    def __read_file__(self, chan, julian_date=(), date=None, stream=False, **readkwargs):
        """
        Read file if exists for a specific date of the station's channel

        :param chan: channel e.j. 'SHZ'
        :param julian_date: tuple of ints. year and julian day
        :param date: datetime object. optional
        :param stream: if True, return stream object
        :return: file_path or stream if exist, or None if not exist

        """

        if chan not in self.stats.channels:
            raise ValueError('Channel not loaded')

        if julian_date:
            year = julian_date[0]
            yday = julian_date[1]

        elif date:
            year = date.year
            yday = date.timetuple().tm_yday

        else:
            raise TypeError("A 'date' or 'julian_date' not specified")

        file_name = '%s.%s.%s.%s.D.%s.%03d' % (self.stats.network.code, self.stats.code, self.stats.location, chan, year, yday)
        file = os.path.join(self.stats.network.sds_path, str(year), self.stats.network.code, self.stats.code, '%s.D' % chan, file_name)
        
        if os.path.isfile(file):
            if stream:
                starttime = readkwargs.get("starttime", None)
                endtime   = readkwargs.get("endtime", None)
                headonly  = readkwargs.get("headonly", False)

                if isinstance(starttime, dt.datetime):
                    starttime = UTCDateTime(starttime)
                
                if isinstance(endtime, dt.datetime):
                    endtime = UTCDateTime(endtime)
                return read2(file, starttime=starttime, endtime=endtime, headonly=headonly)
            
            else:
                return file
        else:
            return False


    def __set_dates__(self):
        """
        Set startdate and enddate of the station object.
        """

        start = dt.datetime(2969,1,1)
        end   = dt.datetime(1969,1,1)

        sds_path = self.stats.network.sds_path
        net_code = self.stats.network.code
        sta_code = self.stats.code

        for chan in self.stats.channels:
            flist = glob(os.path.join(sds_path, '*', net_code, sta_code, f'{chan}.D', '*'))
            flist = list(filter(lambda x : len(x.split('.')[-1]) == 3, flist))
            datelist = [i[-8:] for i in flist]
            datelist.sort()

            if datelist:
                startdate = dt.datetime.strptime(datelist[0], '%Y.%j')
                sd_st = self.__read_file__(chan, date=startdate, stream=True, headonly=True)
                if sd_st:
                    starttime = sd_st[0].stats.starttime.datetime
                    if starttime < start:
                        start = starttime
                
                enddate = dt.datetime.strptime(datelist[-1], '%Y.%j')
                ed_st = self.__read_file__(chan, date=enddate, stream=True, headonly=True)
                if ed_st:
                    endtime   = ed_st[0].stats.endtime.datetime                    
                    if endtime > end:
                        end = endtime

        self.starttime = start
        self.endtime   = end


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

        if resp:
            stream_resp = stream.remove_factor(factor, disp=disp)

        else:
            print(" >>> warn: no response file info")
            stream_resp = stream
        
        return stream_resp


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

        ans = []
        for true_chan in self.stats.channels:
            if true_chan[-1] == component:
                ans.append(true_chan)

        return ans


    def check_channel(self, channel=None):
        """
        Check a channel list or string and return a list of available channels
        if channel is None, it returns stats.chan
        """

        if isinstance(channel, type(None)):
            return self.stats.channels

        if isinstance(channel, str):
            if channel not in self.stats.channels:
                print('channel %s not available' % channel)
                return
            else:
                true_chan = [channel]
            
        else:
            # channel is a list/tuple/array
            true_chan = [ch for ch in channel if ch in self.stats.channels]

        return true_chan


    def get_latlon(self, return_utm=False):
        """
        Get longitude coord. in 'degree' or 'utm'.
        """

        if lat not in self.stats.keys_ or lon not in self.stats.keys_:
            print(" lat/lon not defined in network JSON file")
            return
        
        lat = self.stats.lat
        lon = self.stats.lon

        if return_utm:
            return utm.from_latlon(lat, lon)
        
        else:
            return (lat, lon)


    def get_filelist(self, chan, startdate=None, enddate=None):
        """
        Return a list of files for a station's channel
        :param chan: channel e.j. 'SHZ'
        :param julian_date: tuple of ints. year and julian day
        :param date: datetime object. optinal
        :return: file_path or False
        """

        if chan not in self.stats.channels:
            raise TypeError('Channel not loaded')

        if not startdate:
            startdate = self.starttime

        if not enddate:
            enddate = self.endtime

        day_diff = (enddate - startdate).days
        date_list = [startdate + dt.timedelta(days=i) for i in range(day_diff+1)]
        file_list = [self.__read_file__(chan, date=i) for i in date_list]

        return file_list


    def get_stream(self, starttime, endtime, channel=None, **kwargs):
        """
        Get stream from station object
        :param starttime: datetime
        :param endtime: datetime
        :param chan: list or string e.j.:'SHZ' or 'HHN'. optional
        :param kwargs: remove_sensitivity, remove_response, prefilt, fill_values
        :return: stream2 object
        """

        verbose = kwargs.get('verbose', False)

        assert starttime < endtime, "starttime is greater than endtime"
        assert starttime >= self.starttime, "starttime not valid. No data available"
        assert endtime <= self.endtime, "endtime not valid. No data available"

        # most of MSEED files do not start in 00:00:00 and (finish in 23:59:59), 
        # this is why we need to define time delta to overcome this lack.
        offtime = kwargs.get('offtime', 10)
        time_delta = dt.timedelta(minutes=offtime)
        t0, tf = self.__get_offtimes__(starttime, endtime, time_delta)
        day_diff = (tf.date() - t0.date()).days + 1
        date_list = [t0.date() + dt.timedelta(days=i) for i in range(day_diff)]

        stream = Stream2()
        sample_rate_list = []
        for day in date_list:
            for ch in self.check_channel(channel):
                st_day = self.__read_file__(ch, date=day, stream=True, starttime=t0, endtime=tf)
                
                if st_day:
                    stream += st_day
                    sample_rate_list.append(st_day[0].stats.sampling_rate)
                
                else:
                    if verbose:
                        print(' warn [STA.ID: %s.%s]. No data for %s' % (self.stats.code, self.stats.location, day.strftime('%d %b %Y')))

        # if different traces with different sample rates are in date, discard stream
        if len(list(set(sample_rate_list))) > 1:
            if verbose:
                print('error: stream with mixed sampling rates. Revise your data.')
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

            prefilt = kwargs.get('prefilt', [])

            if prefilt:
                st = st.filter2(fq_band=prefilt)
            
            return st
        
        else:
            return None


    def get_mdata(self, start, end, channel_list, sort=None, verbose=True, **kwargs):
        """
        Return a numpy array with seismic data
        """

        if sort:
            assert len(sort) <= len(channel_list)
            channel_list = [chan for chan in channel_list if chan[-1] in sort]
        else:
            sort = ''.join([chan[-1] for chan in channel_list])

        # get stream
        st = self.get_stream(start, end, channel=channel_list, **kwargs)

        if st:
            # check if nro of channels and traces match up
            if len(channel_list) != len(st):
                print(" error: number of channels do not match with channel_list")
                return None

            # check if npts of all traces in stream match up
            npts_list = [s.stats.npts for s in st]
            if len(set(npts_list)) == 1:
                npts = npts_list[0]
            else:
                print(" error: npts of the channels do not match!")
                return None

            total_sec = (end-start).total_seconds()
            sample_rate = st[0].stats.sampling_rate
            true_npts = int(sample_rate*total_sec)
            
            # get indexes
            if true_npts+1 == npts:
                nin = 0
                nfi = true_npts

                # create the matrix
                mat = np.empty((len(st),true_npts))

                for tr in st:
                    if len(st) > 1:
                        n = sort.index(tr.stats.channel[-1]) 
                    else:
                        n = 0
                    mat[n,nin:nfi] = tr.data[:-1]
            
                return mat
            else:
                print(" error: no full data")
        
        return None


    def polar_analysis(self, starttime, endtime, window=0, olap=0.75, full_analysis=True, **kwargs):
        """ Compute LTE file

        Parameters
        ----------
        starttime : datetime
        
        endtime : datetime
        
        window : int [sec]
            length of the time window (in sec) for moving average over the window. 
            If ``subwindow=0`` no moving average is applied.
        
        olap : float
            overlap percent for moving average over the interval, by default 0.75

        Returns
        -------
        dict object
        """

        assert starttime < endtime
        assert starttime >= self.starttime
        assert endtime <= self.endtime
        assert window >= 0

        sample_rate = kwargs.get('sample_rate', 40)
        fq_band     = kwargs.get('fq_band', (0.5, 10))
        NW          = float(kwargs.get("NW", 3.5))
        pad         = float(kwargs.get("pad", 1.0))

        if window > 0:
            ss   = SSteps(starttime, endtime, -1, window, win_olap=olap, validate=True)
            nwin = int(ss.int_nwin)
            lwin = int(window*sample_rate)
            nadv = float(ss.win_adv)
        else:
            nwin = None
            lwin = None
            nadv = None

        data = self.get_mdata(starttime, endtime, self.stats.channels, sort="ZNE", **kwargs)

        from juliacall import Main as jl
        jl.seval("using LTE")
        polar_ans = jl.polar_run(jl.Array(data), jl.Tuple(fq_band), int(sample_rate), NW, pad, nwin, lwin, nadv, full_analysis)

        return polar_ans


    def lte(self, starttime, endtime, window, subwindow, interval=1, channel=None, subwindow_olap=0.75, **kwargs):
        """ Compute (station) LTE file

        Parameters
        ----------
        starttime : datetime
        
        endtime : datetime
        
        window : int [min]
            length of the time window (in min) to reduce
        
        subwindow : int [sec]
            length of the time window (in sec) for moving average over the window. 
            If ``subwindow=0`` no moving average is applied.
        
        interval : int [h]
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
        assert starttime >= self.starttime
        assert endtime <= self.endtime
        assert window/60 < interval

        # load parameters
        channel     = self.check_channel(channel)
        sample_rate = kwargs.get('sample_rate', 40)
        njobs       = kwargs.get('njobs', -1)
        pad         = kwargs.get("pad",1.0)
        fq_band     = kwargs.get('fq_band', (0.5, 10))
        validate    = kwargs.get("validate", True) # if True, ssteps will ask for confirmation

        # defining base params
        ltebase = dict(
            id              = self.stats.id,
            type            = "station",
            channel         = channel,
            starttime       = starttime.strftime('%Y-%m-%d %H:%M:%S'),
            endtime         = endtime.strftime('%Y-%m-%d %H:%M:%S'),
            interval        = interval,
            window          = window,
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

        ss = SSteps(starttime, endtime, interval*60, window*60, subwindow=subwindow, subw_olap=subwindow_olap, validate=validate)

        if njobs >= mp.cpu_count() or njobs == -1:
            njobs = mp.cpu_count() - 2
        
        if njobs > int(ss.nro_intervals):
            njobs = int(ss.nro_intervals)
            print(f"warn  ::  njobs set to {njobs}")
            
        ltebase["nwin"] = int(ss.int_nwin)
        ltebase["lwin"] = int(window*60*sample_rate)
        ltebase["nro_intervals"] = int(ss.nro_intervals)
        ltebase["nro_time_bins"] = int(ss.int_nwin*ss.nro_intervals)

        if subwindow > 0:
            # if subwindow, a moving average is applied
            lwin = int(subwindow*sample_rate)
            ltebase["lswin"] = lwin
            ltebase["nswin"] = int(ss.nsubwin)
            ltebase["nadv"]  = float(ss.subw_adv)
            ltebase["nro_freq_bins"] = len(get_freq(lwin, sample_rate, fq_band=fq_band, pad=pad)[0])
        else:
            ltebase["lswin"] = None
            ltebase["nswin"] = None
            ltebase["nadv"]  = None
            ltebase["nro_freq_bins"] = len(get_freq(ltebase["lwin"], sample_rate, fq_band=fq_band, pad=pad)[0])
        

        # create hdf5 file and process data
        file_name = kwargs.get("file_name", None)
        if not file_name:
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (self.stats.id, starttime.year, starttime.timetuple().tm_yday, endtime.year, endtime.timetuple().tm_yday, window)
        
        out_dir = kwargs.get("out_dir", None)
        if not out_dir:
            out_dir = os.path.join(seisvo_paths["lte"])
        
        file_name_full = os.path.join(out_dir, file_name)
        if file_name_full.split('.')[-1] != 'lte':
            file_name_full += '.lte'

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)

        lte = StationLTE.new(self, file_name_full, ltebase, njobs)
        
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
            st = self.__read_file__(chan, date=item, stream=True, headonly=True)
            
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



    
