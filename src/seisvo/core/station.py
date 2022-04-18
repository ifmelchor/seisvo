#!/usr/bin/python3
# coding=utf-8

from obspy import UTCDateTime
from datetime import datetime, timedelta
from glob import glob

import os
import numpy as np

from seisvo import __seisvo__
from seisvo.signal import freq_bins
from seisvo.file.lte import LTE
from seisvo.utils.obspyext import read2, Stream2
from seisvo.utils.imports import get_respfile

class Station(object):
    def __init__(self, StaInfo):
        self.info = StaInfo
        self.set_dates()

    def __str__(self):
        return self.info.__str__()

    def is_infrasound(self):
        # Check if the station has an infrasound channel
        if len(self.info.chan) == 1:
            if self.info.chan[0][-1] == 'P':
                return True
        return False

    def is_response(self):
        # Check if the instrumental response can be removed
        if get_respfile(self.info.net, self.info.code, self.info.loc):
            return True
        else:
            return False

    def is_component(self, component):
        # Check if the component exist
        ans = False
        for true_chan in self.info.chan:
            if true_chan[-1] == component:
                ans = True
        return ans


    def get_chan(self, component):
        # get all the channel with the same component
        ans = []
        for true_chan in self.info.chan:
            if true_chan[-1] == component:
                ans +=[true_chan]
        return ans


    def is_3component(self):
        # Check if the station is three component

        if self.is_component('Z') and self.is_component('E') and self.is_component('N'):
            return True

        else:
            return False


    def get_latlon(self, degree=True):
        """
        Get longitude coord.
        """

        import utm

        lat = self.info.lat
        lon = self.info.lon

        if degree:
            if self.info.type == 'utm':
                zn = self.info.zone_number
                zl = self.info.zone_letter
                return utm.to_latlon(lat, lon, zn, zl)
            else:
                return (lat, lon)

        else:
            if self.info.type == 'utm':
                zn = self.info.zone_number
                zl = self.info.zone_letter
                return (lat, lon, zn, zl)
            else:
                return utm.from_latlon(lat, lon)


    def get_sampling_rate(self, starttime=None):
        # reed first file
        one_min = timedelta(minutes=1)

        if starttime is None:
            starttime = self.info.starttime

        st = self.get_stream(starttime, starttime + one_min)
        return st[0].stats.sampling_rate


    def set_dates(self, chan=None):
        """
        Set startdate and enddate of the station object.
        By default, first channel is used to compute.
        :param chan: optional, first channel by default
        """

        if not chan:
            chan = self.info.chan[0]
        else:
            if chan not in self.info.chan:
                raise TypeError('Channel not loaded')

        list = glob(os.path.join(self.info.sdsdir, '*', self.info.net, self.info.code, '%s.D' % chan, '*'))
        datelist = [i[-8:] for i in list]
        datelist.sort()

        if datelist:
            startdate = datetime.strptime(datelist[0], '%Y.%j')
            enddate = datetime.strptime(datelist[-1], '%Y.%j')

            sd_st = self.is_file(chan, date=startdate, stream=True)
            ed_st = self.is_file(chan, date=enddate, stream=True)

            self.info.starttime = sd_st[0].stats.starttime.datetime
            self.info.endtime = ed_st[0].stats.endtime.datetime

        else:
            self.info.starttime = None
            self.info.endtime = None


    def is_file(self, chan, julian_date=(), date=None, stream=False, status=False):
        """
        Return True or False if exist file for a specific date of the station's channel
        :param chan: channel e.j. 'SHZ'
        :param julian_date: tuple of ints. year and julian day
        :param date: datetime object. optinal
        :param stream: if True, return stream object
        :return: file_path or False
        """

        if chan not in self.info.chan:
            raise TypeError('Channel not loaded')

        if julian_date:
            year = julian_date[0]
            yday = julian_date[1]

        elif date:
            year = date.year
            yday = date.timetuple().tm_yday

        else:
            raise TypeError('A date must be specified')

        file_name = '%s.%s.%s.%s.D.%s.%03d' % (
            self.info.net,
            self.info.code,
            self.info.loc,
            chan,
            year,
            yday
            )

        file = os.path.join(
            self.info.sdsdir,
            str(year),
            self.info.net,
            self.info.code,
            '%s.D' % chan,
            file_name)

        if os.path.isfile(file):
            if stream:
                return read2(file)
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

        if chan not in self.info.chan:
            raise TypeError('Channel not loaded')

        if not startdate:
            startdate = self.info.starttime

        if not enddate:
            enddate = self.info.endtime

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

        from seisvo.utils.plotting import pplot_control

        if chan not in self.info.chan:
            raise TypeError('Channel not loaded')

        if not startdate:
            startdate = self.info.starttime

        if not enddate:
            enddate = self.info.endtime

        day_diff = (enddate - startdate).days
        date_list = [startdate + timedelta(days=i) for i in range(day_diff+1)]

        npts = []
        sample_rate = []
        nro_traces = []
        max_value = []

        for i, item in enumerate(date_list):
            print ('  reading data... (%d%%)' % (100*i/len(date_list)), end='\r')
            st = self.is_file(chan, date=item, stream=True)
            if st:
                npts += [sum([len(tr.data) for tr in st])]
                sample_rate += [st[0].stats.sampling_rate]
                nro_traces += [len(st)]
                max_value += [max([max(tr.data) for tr in st])]
            else:
                npts += [None]
                sample_rate += [None]
                nro_traces += [None]
                max_value += [None]

        print('\t\t')

        title = '%s \n %s -- %s' % (self.info.id,
            startdate.strftime('%Y.%m.%d'),
            enddate.strftime('%Y.%m.%d'))
        pplot_control(title, date_list, nro_traces, sample_rate, npts, max_value)


    def get_stream(self, starttime, endtime, chan=None, component=None, 
        remove_response=False, sample_rate=None, pad_zeros=False, **kwargs):
        """
        Get stream from station object
        :param starttime: datetime
        :param endtime: datetime
        :param component: string e.j.:'Z' or 'N'. optional
        :param chan: string e.j.:'SHZ' or 'HHN'. optional
        :param remove_response: boolean
        :param pad_zeros: boolean
        :param kwargs: taper, taper_p, prefilt, fill_values
        :return: stream2 object
        """

        prefilt = kwargs.get('prefilt', [])
        fill_value = kwargs.get('fill_value', None)
        verbose = kwargs.get('verbose', True)

        if starttime > endtime:
            raise ValueError('starttime is greater than endtime')

        starttime_2 = starttime - timedelta(minutes=2)
        endtime_2 = endtime + timedelta(minutes=2)

        if chan:
            component = chan[-1]

        if component:
            if not self.is_component(component):
                raise TypeError('component not loaded')
            else:
                chan = self.get_chan(component)

        else:
            chan = self.info.chan

        day_diff = (endtime_2.date() - starttime_2.date()).days
        date_list = [starttime_2.date() + timedelta(days=i) for i in range(day_diff+1)]

        stream = Stream2()
        for day in date_list:
            for ch in chan:
                st_day = self.is_file(ch, date=day, stream=True)
                if st_day:
                    stream += st_day
                else:
                    if verbose:
                        print(' Error [%s.%s]: No data for %s' % (self.info.code,
                    self.info.loc, day.strftime('%Y.%m.%d')))

        if pad_zeros:
            fill_value = 0

        stream.merge(fill_value=fill_value)
        st = stream.slice(UTCDateTime(starttime), UTCDateTime(endtime))

        if not st:
            print(' Warning: [get_stream] No data found for selected dates')
            return None

        else:
            if pad_zeros:
                for trace in st:
                    true_starttime = trace.stats.starttime.datetime
                    true_endtime = trace.stats.endtime.datetime
                    fs = trace.stats.sampling_rate

                    if starttime != true_starttime:
                        sec_diff = abs((starttime - true_starttime).total_seconds())
                        zeros = np.zeros(int(fs*sec_diff))
                        trace.data = np.concatenate((zeros, trace.data))
                        trace.starttime = UTCDateTime(starttime)

                    if endtime != true_endtime:
                        sec_diff = abs((true_endtime - endtime).total_seconds())
                        zeros = np.zeros(int(fs*sec_diff))
                        trace.data = np.concatenate((trace.data, zeros))
                        trace.stats.npts += int(fs*sec_diff)
            else:
                t1 = max([tr.stats.starttime for tr in st])
                t2 = min([tr.stats.endtime for tr in st])
                st = Stream2(st.slice(t1, t2))

            # removing instrument response
            if remove_response:
                resp = get_respfile(self.info.net, self.info.code, self.info.loc)
                st = st.remove_response2(resp_dict=resp, **kwargs)

            if sample_rate:
                if sample_rate < st[0].stats.sampling_rate:
                    st.resample(sample_rate)

            if prefilt:
                st = st.filter2(fq_band=prefilt)

            return st


    def lte(self, starttime, endtime, chan, time_step, polargram=False, sample_rate=None,
        avg_step=1, polarization_attr=False, **kwargs):

        info = {'id': '%s.%s' % (self.info.id, chan)}
        info['starttime'] = starttime.strftime('%Y-%m-%d %H:%M:%S')
        info['endtime'] = endtime.strftime('%Y-%m-%d %H:%M:%S')

        info['time_step'] = time_step
        if avg_step is None:
            info['avg_step'] = time_step
        else:
            info['avg_step'] = avg_step

        # computing time_bins
        nro_time_bins = 0
        start_time = starttime
        time_delta = timedelta(minutes=time_step)
        while start_time + time_delta <= endtime:
            nro_time_bins += 1
            start_time += time_delta
        info['time_bins'] = nro_time_bins

        # computing freq_bins
        if not sample_rate:
            sample_rate = int(self.get_sampling_rate(starttime))
        info['sampling_rate'] = sample_rate
        
        fq_band = kwargs.get('fq_band', (1, 10))
        nro_freq_bins = freq_bins(sample_rate*60*info['avg_step'], sample_rate, fq_band=fq_band, nfft='uppest')
        info['fq_band'] = fq_band
        info['freq_bins'] = nro_freq_bins

        # create hdf5 file
        file_name = kwargs.get("file_name", None)
        if not file_name:
            lte_id = '.'.join(info['id'].split('.')[1:])
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (lte_id, starttime.year, starttime.timetuple().tm_yday,
                endtime.year, endtime.timetuple().tm_yday, time_step)
        
        out_dir = kwargs.get("out_dir", None)
        if not out_dir:
            out_dir = os.path.join(__seisvo__, 'lte', self.info.net)
        
        file_name_full = os.path.join(out_dir, file_name)

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)

        # define other parameters
        info['tau'] = kwargs.get("tau", 1)
        info['p_order'] = kwargs.get("p_order", 5)
        info['threshold'] = kwargs.get("threshold", 0.0)
        info['time_bandwidth'] = kwargs.get('time_bandwidth', 3.5)
        info['remove_response'] = int(kwargs.get("remove_response", False))
        info['polargram'] = polargram
        info['matrix_return'] = polarization_attr

        # create new file and process data
        lte = LTE.new(file_name_full, info, **kwargs)

        return lte, file_name


    def plot(self, channel, starttime, delta=15, return_fig=False, **kwargs):
        """ Plot seismogram of the station in a simple GUI. 

        Args:
            channel (str or list): string or list of string with the channels to plot
            starttime (datetime): start time to plot
            delta (int, optional): in minutes. Defaults to 15.
            return_fig (bool, optional): Return fig and axes objects. Defaults to False.
        """

        from seisvo.gui.gstation import plot_station, get_fig
        
        if return_fig:
            fig, axes = get_fig(self, starttime, channel, delta, **kwargs)
            return fig, axes

        else:
            plot_station(self, channel, starttime, delta, **kwargs)



    
