#!/usr/bin/python3
# coding=utf-8

import os
import utm
from datetime import datetime, timedelta
from glob import glob

from seisvo import __seisvo__
from seisvo.core import get_respfile
from seisvo.core.obspyext import UTCDateTime, read2, Stream2
from seisvo.signal import freq_bins
from seisvo.plotting import pplot_control

class Station(object):
    def __init__(self, StaFile):
        self.stats = StaFile
        self.resp_ = get_respfile(self.stats.net, self.stats.code, self.stats.loc)
        self.set_dates()


    def __str__(self):
        return self.stats.__str__()


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
            startdate = datetime.strptime(datelist[0], '%Y.%j')
            enddate = datetime.strptime(datelist[-1], '%Y.%j')

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
        date_list = [startdate + timedelta(days=i) for i in range(day_diff+1)]
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
        date_list = [startdate + timedelta(days=i) for i in range(day_diff+1)]

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

        starttime_2 = starttime - timedelta(minutes=2)
        endtime_2 = endtime + timedelta(minutes=2)

        if not isinstance(channel, (str, list, type(None))):
            raise ValueError('channel should be list or string')


        if isinstance(channel, type(None)):
            channel = self.stats.chan

        else:
            if isinstance(channel, str):
                if channel not in self.stats.chan:
                    raise ValueError('channel %s not available' % channel)
            
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

        day_diff = (endtime_2.date() - starttime_2.date()).days
        date_list = [starttime_2.date() + timedelta(days=i) for i in range(day_diff+1)]

        stream = Stream2()
        sample_rate_list = []
        for day in date_list:
            for ch in channel:
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
            print(' warn:  No data found for selected dates. Return None')
            return None

        else:
            t1 = max([tr.stats.starttime for tr in st])
            t2 = min([tr.stats.endtime for tr in st])
            st = Stream2(st.slice(t1, t2))

            if remove_response:
                if self.resp_:
                    st = st.remove_response2(resp_dict=self.resp_, **rrkwargs)
                else:
                    print('warn: no responsefile found')

            if isinstance(sample_rate, int):
                if sample_rate < sample_rate_list[0]:
                    st.resample(sample_rate)

            if prefilt:
                st = st.filter2(fq_band=prefilt)

            return st


    def lte(self, starttime, endtime, time_step, channel=None, avg_step=1, **kwargs):
        """Compute LTE file

        Parameters
        ----------
        starttime : datetime
        endtime : datetime
        time_step : int, time window in minutes
        time_olap : int, time window in minutes (< time_step)
        chan : list or string, optional
            channels to be analyzed. If None, computes all available channels, by default None
        avg_step : None or int, optional
            time window for moving average (< time_step). If None, avg_step = time_step, by default 1
        kwargs : sample_rate, polar_degree, polarization_attr, 

        Returns
        -------
        LTE file
        """
        from seisvo.file.lte import LTE

        # defining kwargs
        sample_rate = kwargs.get('sample_rate', None)
        remove_response = kwargs.get('remove_response', False)
        polar_degree = kwargs.get('polar_degree', False)
        polarization_attr = kwargs.get('polarization_attr', False)
        fq_band = kwargs.get('fq_band', (1, 10))
        file_name = kwargs.get("file_name", None)
        out_dir = kwargs.get("out_dir", None)
        pe_tau = kwargs.get("tau", 1)
        pe_order = kwargs.get("p_order", 5)
        threshold = kwargs.get("threshold", 0.0)
        time_bandwidth = kwargs.get('time_bandwidth', 3.5)
        th_outlier = kwargs.get("th_outlier", 3.2)
        ltekwargs = kwargs.get("ltekwargs", {})

        if not chan:
            # get vertical component as default
            chan = [c for c in self.stats.chan if c[-1]=='Z'][0]

        lte_stats = {}
        lte_stats['id'] = self.stats.id
        lte_stats['chan'] = channel
        # info = {'id': '%s.%s' % (self.stats.id, chan)}
        lte_stats['starttime'] = starttime.strftime('%Y-%m-%d %H:%M:%S')
        lte_stats['endtime'] = endtime.strftime('%Y-%m-%d %H:%M:%S')

        lte_stats['time_step'] = time_step
        if avg_step is None:
            lte_stats['avg_step'] = time_step
        else:
            lte_stats['avg_step'] = avg_step

        # computing time_bins
        nro_time_bins = 0
        start_time = starttime
        time_delta = timedelta(minutes=time_step)
        while start_time + time_delta <= endtime:
            nro_time_bins += 1
            start_time += time_delta
        lte_stats['time_bins'] = nro_time_bins

        # computing freq_bins
        if not sample_rate:
            sample_rate = self.stats.sampling_rate
        lte_stats['sampling_rate'] = sample_rate
        
        nro_freq_bins = freq_bins(sample_rate*60*info['avg_step'], sample_rate, fq_band=fq_band, nfft='uppest')
        lte_stats['fq_band'] = fq_band
        lte_stats['freq_bins'] = nro_freq_bins

        # create hdf5 file
        if not file_name:
            lte_id = '.'.join(info['id'].split('.')[1:])
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (lte_id, starttime.year, starttime.timetuple().tm_yday,
                endtime.year, endtime.timetuple().tm_yday, time_step)
        
        if not out_dir:
            out_dir = os.path.join(__seisvo__, 'lte', self.stats.net)
        
        file_name_full = os.path.join(out_dir, file_name)

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)

        # define other parameters
        info['tau'] = pe_tau
        info['p_order'] = pe_order
        info['threshold'] = threshold
        info['time_bandwidth'] = time_bandwidth
        info['remove_response'] = int(remove_response)
        info['th_outlier'] = th_outlier
        info['polar_degree'] = polar_degree
        info['matrix_return'] = polarization_attr

        # create new file and process data
        lte = LTE.new(file_name_full, info, **ltekwargs)

        return lte, file_name


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



    
