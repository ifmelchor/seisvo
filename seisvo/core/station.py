#!/usr/bin/python3
# coding=utf-8

import os
import utm
import datetime as dt
import numpy as np
from glob import glob

from seisvo import LTE_PATH
from seisvo.core import get_respfile
from seisvo.core.obspyext import UTCDateTime, read2, Stream2
from seisvo.signal import freq_bins, time_bins
from seisvo.signal.polarization import PolarAnalysis
from seisvo.plotting import pplot_control

class Station(object):
    def __init__(self, StaFile):
        self.stats = StaFile
        self.resp_ = get_respfile(self.stats.net, self.stats.code, self.stats.loc)
        self.set_dates()


    def __str__(self):
        return self.stats.__str__()


    def __get_channel__(self, channel):
        
        if isinstance(channel, type(None)):
            channel = self.stats.chan

        else: 
            if not isinstance(channel, (str, list, np.ndarray)):
                raise ValueError('channel should be list or string')

        if isinstance(channel, str):
            if channel not in self.stats.chan:
                raise ValueError('channel %s not available' % channel)
            
            else:
                true_chan = [channel]
            
        else:
            true_chan = []
            for ch in channel:
                if ch in self.stats.chan:
                    true_chan.append(ch)
                else:
                    print('warn: channel %s not available' % ch)
        
        channel = true_chan
        
        if not channel:
            raise ValueError(' no channel defined!')
        
        return channel


    def is_infrasound(self):
        # Check if the station is an infrasound channel
        if len(self.stats.chan) == 1:
            if self.stats.chan[0][-1] == 'P':
                return True
        return False


    def is_component(self, component):
        # Check if the component exist
        ans = False
        for true_chan in self.stats.chan:
            if true_chan[-1] == component:
                ans = True
        return ans


    def get_chan(self, component):
        # get all the channel with the same component
        ans = []
        for true_chan in self.stats.chan:
            if true_chan[-1] == component:
                ans +=[true_chan]
        return ans


    def is_three_component(self):
        if self.is_component('Z') and self.is_component('E') and self.is_component('N'):
            return True
        else:
            return False


    def get_latlon(self, degree=True):
        """
        Get longitude coord.
        """

        lat = self.stats.lat
        lon = self.stats.lon

        if degree:
            if self.stats.type == 'utm':
                zn = self.stats.zone_number
                zl = self.stats.zone_letter
                return utm.to_latlon(lat, lon, zn, zl)
            else:
                return (lat, lon)

        else:
            if self.stats.type == 'utm':
                zn = self.stats.zone_number
                zl = self.stats.zone_letter
                return (lat, lon, zn, zl)
            else:
                return utm.from_latlon(lat, lon)


    def set_dates(self, chan=None):
        """
        Set startdate and enddate of the station object.
        By default, first channel is used to compute.
        :param chan: optional, first channel by default
        """

        if not chan:
            chan = self.stats.chan[0]
        else:
            if chan not in self.stats.chan:
                raise TypeError('Channel not loaded')

        list = glob(os.path.join(self.stats.sdsdir, '*', self.stats.net, self.stats.code, '%s.D' % chan, '*'))
        datelist = [i[-8:] for i in list]
        datelist.sort()

        if datelist:
            startdate = dt.datetime.strptime(datelist[0], '%Y.%j')
            enddate = dt.datetime.strptime(datelist[-1], '%Y.%j')

            sd_st = self.is_file(chan, date=startdate, stream=True)
            ed_st = self.is_file(chan, date=enddate, stream=True)

            self.starttime = sd_st[0].stats.starttime.datetime
            self.endtime = ed_st[0].stats.endtime.datetime

        else:
            self.starttime = None
            self.endtime = None


    def is_file(self, chan, julian_date=(), date=None, stream=False, headonly=False):
        """
        Return True or False if exist file for a specific date of the station's channel
        :param chan: channel e.j. 'SHZ'
        :param julian_date: tuple of ints. year and julian day
        :param date: datetime object. optinal
        :param stream: if True, return stream object
        :return: file_path or False
        """

        if chan not in self.stats.chan:
            raise ValueError('Channel not loaded')

        if julian_date:
            year = julian_date[0]
            yday = julian_date[1]

        elif date:
            year = date.year
            yday = date.timetuple().tm_yday

        else:
            raise TypeError('A date must be specified')

        file_name = '%s.%s.%s.%s.D.%s.%03d' % (
            self.stats.net,
            self.stats.code,
            self.stats.loc,
            chan,
            year,
            yday
            )

        file = os.path.join(
            self.stats.sdsdir,
            str(year),
            self.stats.net,
            self.stats.code,
            '%s.D' % chan,
            file_name)

        if os.path.isfile(file):
            if stream:
                return read2(file, headonly=headonly)
            else:
                return file

        else:
            return False


    def get_filelist(self, chan, startdate=None, enddate=None):
        """
        Return a list of files for a station's channel
        :param chan: channel e.j. 'SHZ'
        :param julian_date: tuple of ints. year and julian day
        :param date: datetime object. optinal
        :return: file_path or False
        """

        if chan not in self.stats.chan:
            raise TypeError('Channel not loaded')

        if not startdate:
            startdate = self.stats.starttime

        if not enddate:
            enddate = self.stats.endtime

        day_diff = (enddate - startdate).days
        date_list = [startdate + dt.timedelta(days=i) for i in range(day_diff+1)]
        file_list = [self.is_file(chan, date=i) for i in date_list]

        return file_list


    def plot_control(self, chan, startdate=None, enddate=None):
        """
        Plot information of your mseed files
        :param chan: channel e.j. 'SHZ'
        :param starttime: datetime object, optional
        :param endtime: datetime object. optinal
        """

        if chan not in self.stats.chan:
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
            st = self.is_file(chan, date=item, stream=True, headonly=True)
            
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


    def get_stream(self, starttime, endtime, channel=None, **kwargs):
        """
        Get stream from station object
        :param starttime: datetime
        :param endtime: datetime
        :param chan: list or string e.j.:'SHZ' or 'HHN'. optional
        :param kwargs: remove_response, prefilt, fill_values
        :return: stream2 object
        """

        sample_rate = kwargs.get('sample_rate', None)
        remove_response = kwargs.get('remove_response', False)
        rrkwargs = kwargs.get('rrkwargs', {})
        prefilt = kwargs.get('prefilt', [])
        fill_value = kwargs.get('fill_value', None)
        verbose = kwargs.get('verbose', True)

        if starttime > endtime:
            raise ValueError('starttime is greater than endtime')
        
        if starttime < self.starttime:
            raise ValueError('starttime not valid')

        if endtime > self.endtime:
            raise ValueError('endtime not valid')

        # most of MSEED files do not start in 00:00:00 and (finish in 23:59:59), this is why we need to define time delta to overcome this lack.
        time_delta = dt.timedelta(minutes=10)

        if starttime - time_delta > self.starttime:
            starttime_mod = starttime - time_delta
        else:
            starttime_mod = starttime

        if endtime + time_delta < self.endtime:
            endtime_mod = endtime + time_delta
        else:
            endtime_mod = endtime

        day_diff = (endtime_mod.date() - starttime_mod.date()).days + 1
        date_list = [starttime_mod.date() + dt.timedelta(days=i) for i in range(day_diff)]

        stream = Stream2()
        sample_rate_list = []
        for day in date_list:
            for ch in self.__get_channel__(channel):
                st_day = self.is_file(ch, date=day, stream=True)
                if st_day:
                    stream += st_day
                    sample_rate_list.append(int(st_day[0].stats.sampling_rate))
                else:
                    if verbose:
                        print(' warn [STA.ID: %s.%s]. No data for %s' % (self.stats.code, self.stats.loc, day.strftime('%d %b %Y')))

        if len(list(set(sample_rate_list))) > 1:
            raise ValueError('error: stream with mixed sampling rates. Revise your data.')

        stream.merge(fill_value=fill_value)
        st = stream.slice(UTCDateTime(starttime), UTCDateTime(endtime))

        if not st:
            print(' warn:  No data found for selected dates.')
            return None

        else:
            t1 = max([tr.stats.starttime for tr in st])
            t2 = min([tr.stats.endtime for tr in st])
            st = Stream2(st.slice(t1, t2))

            if remove_response:
                if self.resp_: # can be list or RespFile

                    if isinstance(self.resp_, list):
                        for resp in self.resp_:
                            if resp.endtime > t2:
                                st = st.remove_response2(resp_dict=resp, **rrkwargs)
                            
                    else:
                        st = st.remove_response2(resp_dict=self.resp_, **rrkwargs)
                
                else:
                    print('warn: no responsefile found')

            if sample_rate:
                if sample_rate < sample_rate_list[0]:
                    st = Stream2(st.resample(sample_rate))

            if prefilt:
                st = st.filter2(fq_band=prefilt)

            return st


    def polar_analysis(self, starttime, endtime, step_sec=None, olap=0.0, fq_band=(0.5, 10), full_analysis=False, **kwargs):
        
        if not self.is_three_component():
            raise ValueError(' polarization analysis is only for three-component station')
        
        st_kwargs = dict(
            sample_rate = kwargs.get('sample_rate', None),
            remove_response = kwargs.get('remove_response', False),
            prefilt = kwargs.get('prefilt', []),
            fill_value = kwargs.get('fill_value', None),
            rrkwargs = kwargs.get('rrkwargs', {}),
            verbose = kwargs.get('verbose', True)
        )

        stream = self.get_stream(starttime, endtime, **st_kwargs)

        z_data = stream.get_component('Z')
        n_data = stream.get_component('N').get_data(detrend=True)
        e_data = stream.get_component('E').get_data(detrend=True)

        sampling_rate = z_data.stats.sampling_rate

        pa_kwargs = dict(
            taper = kwargs.get('taper', False),
            taper_p = kwargs.get('taper_p', 0.05),
            time_bandwidth = kwargs.get('time_bandwidth', 3.5)
        )

        if isinstance(step_sec, float):
            npts_mov_avg = step_sec*sampling_rate
        else:
            npts_mov_avg = None

        pa = PolarAnalysis(z_data.get_data(detrend=True), n_data, e_data, sampling_rate, npts_mov_avg=npts_mov_avg, olap=olap, fq_band=fq_band, full_analysis=full_analysis, **pa_kwargs)

        if full_analysis:
            return pa.freq, pa.polar_dgr, pa.rect, pa.azimuth, pa.elevation
        
        else:
            return pa.freq, pa.polar_dgr


    def lte(self, starttime, endtime, channel=None, interval=20, int_olap=0.0, step=1, step_olap=0.0, **kwargs):
        """ Compute LTE file

        Parameters
        ----------
        starttime : datetime
        endtime : datetime
        channel : _type_, optional
            _description_, by default None
        interval : int, optional
            interval length for LTE file, by default 20
        int_olap : float, optional
            overlap of the interval length, by default 0.0
        step : int, optional
            interval length for moving average of the interval, by default 1
        step_olap : float, optional
            overlap for moving average in the interval, by default 0.0

        Returns
        -------
        LTE
            LTE object
        """

        from seisvo.file.lte import LTE

        channel = self.__get_channel__(channel)

        if starttime > endtime:
            raise ValueError('starttime is greater than endtime')
        
        if starttime < self.starttime:
            raise ValueError('starttime not valid')
        
        if endtime > self.endtime:
            raise ValueError('endtime not valid')
        
        if not (0 <= int_olap < 1):
            raise ValueError(' time_olap should be float between 0 and 1')
    
        if step > interval:
            raise ValueError(' step is greater than interval')
        
        if not (0 <= step_olap < 1):
            raise ValueError(' step_olap should be float between 0 and 1')

        # defining kwargs
        lte_header = dict(
            id = self.stats.id,
            channel = channel,
            starttime = starttime.strftime('%Y-%m-%d %H:%M:%S'),
            endtime = endtime.strftime('%Y-%m-%d %H:%M:%S'),
            interval = interval,
            int_olap = int_olap,
            step = step,
            step_olap = step_olap,
            fq_band = kwargs.get('fq_band', (0.5, 10)),
            sample_rate = kwargs.get('sample_rate', self.stats.sampling_rate),
            remove_response = kwargs.get('remove_response', False),
            polar_degree = kwargs.get('polar_degree', False),
            polar_analysis = kwargs.get('polar_analysis', False),
            f_params = kwargs.get('f_params', False),
            opt_params = kwargs.get('opt_params', False),
            f_threshold = kwargs.get("f_threshold", False),
            PE_tau = kwargs.get("pe_tau", 2),
            PE_order = kwargs.get("pe_order", 7),
            time_bandwidth = kwargs.get('time_bandwidth', 3.5)
        )

        # other kwargs
        file_name = kwargs.get("file_name", None)
        out_dir = kwargs.get("out_dir", './')
        auto_chunk = kwargs.get("auto_chunk", True)
        ltekwargs = kwargs.get("ltekwargs", {})

        # compute time and freq bins
        nro_time_bins = time_bins(starttime, endtime, interval, int_olap)
        lte_header['nro_time_bins'] = nro_time_bins

        nro_freq_bins = freq_bins(lte_header['sample_rate']*60*step, lte_header['sample_rate'], fq_band=lte_header['fq_band'], nfft='uppest')
        lte_header['nro_freq_bins'] = nro_freq_bins

        if lte_header['polar_analysis']:
            lte_header['polar_degree'] = True
            
        if isinstance(lte_header['f_threshold'], float):
            lte_header['f_params'] = True
        
        # create hdf5 file and process data
        if not file_name:
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (lte_header['id'], starttime.year, starttime.timetuple().tm_yday, endtime.year, endtime.timetuple().tm_yday, interval)
        
        if not out_dir:
            out_dir = os.path.join(LTE_PATH)
        
        file_name_full = os.path.join(out_dir, file_name)

        if file_name_full.split('.')[-1] != 'lte':
            file_name_full += '.lte'

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)

        lte = LTE.new(self, file_name_full, lte_header, auto_chunk=auto_chunk, **ltekwargs)

        return lte


    def plot(self, starttime, sde, channel='all', delta=30, app=False, **kwargs):
        """ Plot seismogram of the station in a simple GUI. 

        Args:
            channel (str or list): string or list of string with the channels to plot
            starttime (datetime): start time to plot
            delta (int, optional): in minutes. Defaults to 15.
            return_fig (bool, optional): Return fig and axes objects. Defaults to False.
        """

        from seisvo.plotting.gui.gstation import plot_station_gui

        window = plot_station_gui(self, starttime, sde, channel=channel, delta=delta, app=app, **kwargs)
        
        if app:
            
            return window



    
