#!/usr/bin/python3
# coding=utf-8

'''

Read write and operate air files

'''

import h5py
import numpy as np
import pandas as pd
from obspy.core.util.attribdict import AttribDict
from datetime import datetime, timedelta
from seisvo.utils.plotting import plot_out_ccorr


class AiRstats(AttribDict):
    def __init__(self, header):
        super(AiRstats, self).__init__(header)

    def __str__(self):
        priorized_keys = [
        'code',
        'sampling_rate',
        'starttime',
        'endtime',
        'time_width',
        'overlap',
        'time_bins',
        'fq_band',
        'radii',
        'vel_air', 
        'h_src', 
        'src_dgr',
        'sensors_loc']

        return self._pretty_str(priorized_keys)


class AiR(object):
    def __init__(self, air_file):
        self.air_file = air_file
        self.set_stats()
        self.set_delta_times()
        self.__groups__ = ('crosscorr', 'pressure')
        self.__attrs__ = ('crosscorr', 'p_avg', 'p_max')


    def __str__(self):
        return self.stats.__str__()


    def set_stats(self):
        with h5py.File(self.air_file, "r") as f:
            hdr = f['header']
            self.stats = AiRstats({
                'code': hdr.attrs['code'],
                'sampling_rate': hdr.attrs['sampling_rate'],
                'starttime': datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                'endtime': datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                'time_width': int(hdr.attrs['time_width']),
                'time_bins': int(hdr.attrs['time_bins']),
                'fq_band': hdr.attrs['fq_band'],
                'overlap': hdr.attrs['overlap'],
                'radii': hdr.attrs['radii'],
                'vel_air': hdr.attrs['vel_air'],
                'h_src': hdr.attrs['h_src'],
                'src_dgr': hdr.attrs['src_dgr'],
                'sensors_loc': hdr.attrs['sensors_loc']
                })


    def set_delta_times(self):
        with h5py.File(self.air_file, "r") as f:
            dt = f['delta_times']
            self.delta_times = {}
            for key in dt.keys():
                self.delta_times[key] = dt.get(key)[:]

    
    def get_group(self, attr):
        if attr in ('p_avg', 'p_max'):
            return 'pressure'

        elif attr in ('crosscorr'):
            return 'crosscorr'
        
        else:
            return None


    def print_shape(self):
        " Shape and chunk shape:"

        with h5py.File(self.air_file, "r") as f:

            dset = f.get('crosscorr')['crosscorr']
            print('  shape: crosscorr/crosscorr  >>> ', dset.shape)
            print('  chunk_shape: crosscorr/crosscorr  >>> ', dset.chunks)

            dset = f.get('pressure')['p_avg']
            print('  shape: pressure/p_avg  >>> ', dset.shape)
            print('  chunk_shape: pressure/p_avg  >>> ', dset.chunks)

            dset = f.get('pressure')['p_max']
            print('  shape: pressure/p_max  >>> ', dset.shape)
            print('  chunk_shape: pressure/p_max  >>> ', dset.chunks)


    def save_data(self, group, attr, data, item):
        """
        This function replace a dataset of the hdf5 file
        """
        with h5py.File(self.air_file, "r+") as f:
            dset = f.get(group)[attr]
            dset[item,:] = data[0,:]
            f.flush()


    def get_dataset(self, group, attr, index):
        with h5py.File(self.air_file, "r") as f:
            dset = f.get(group)[attr][index[0]:index[1],:]        
        return dset


    def get_index(self, starttime=None, endtime=None):

        if starttime:
            st_diff = int((starttime - self.stats.starttime).total_seconds())
            start_index = int(st_diff/(self.stats.time_width-self.stats.overlap))
        else:
            start_index = 0

        if start_index < 0 or start_index >= self.stats.time_bins:
            raise ValueError(' Select a correct period of time ')

        if endtime:
            ed_diff = int((endtime - self.stats.starttime).total_seconds())
            end_index = int(ed_diff/(self.stats.time_width-self.stats.overlap))+1
        else:
            end_index = self.stats.time_bins

        if end_index < 0 or end_index > self.stats.time_bins:
            raise ValueError(' Select a correct period of time ')
        
        return (start_index, end_index)


    def get_time(self, starttime=None, endtime=None):
        npts = self.get_index(starttime, endtime)
        n = int(self.stats.time_width-self.stats.overlap)
        times = [self.stats.starttime + timedelta(seconds=k*n) for k in range(npts[0],npts[1])]
        return np.array(times)


    def get_attr(self, attr=None, starttime=None, endtime=None):
        """
        Get array from file
        :param starttime: datetime
        :param endtime: datetime
        :param attr: string. see print_attr_list for details.
        :return: numpy array
        """

        gr = self.get_group(attr)
        (i, j) = self.get_index(starttime, endtime)
        return self.get_dataset(gr, attr, (i,j))


    def get_dict(self, starttime=None, endtime=None):
        out = {}
        out['time'] = self.get_time(starttime=starttime, endtime=endtime)
        out['mcorr'] = self.get_attr('crosscorr', starttime=starttime, endtime=endtime)
        out['p_avg'] = self.get_attr('p_avg', starttime=starttime, endtime=endtime)
        out['p_max'] = self.get_attr('p_max', starttime=starttime, endtime=endtime)
        out['info'] = {
            'sampling_rate': self.stats.sampling_rate,
            'starttime' : self.stats.starttime,
            'time_width': self.stats.time_width,
            'overlap': self.stats.overlap,
            'fq_band': self.stats.fq_band,
            'code': self.stats.code,
            'sensors_loc': self.stats.sensors_loc
        }
        out['model'] = {
            'radii' : self.stats.radii,
            'vel_air' : self.stats.vel_air,
            'h_src' : self.stats.h_src,
            'src_dgr' : self.stats.src_dgr
        }

        return out


    def plot_gui(self, **kwargs):
        from seisvo.gui.giarr import gplot
        gplot(air_file=self, **kwargs)


    def daily_plot(self, startday=None, endday=None, maxSMB=0.5):

        if startday is None:
            startday = self.get_time()[0]

        if endday is None:
            endday = self.get_time()[-1]

        rang = (endday-startday).days
        times = [startday + timedelta(days=x) for x in range(0, rang+1)]
        
        for t1, t2 in zip(times[:-1],times[1:]):
            print(t1, t2)
            out = self.get_dict(t1, t2)
            plot_out_ccorr(out, maxSMB=maxSMB, save=True)


    @staticmethod
    def new(air_file, headers, delta_times, sstep, chunk_weight=30):
        """
        Create new AiR (hdf5) file
        """

        f = h5py.File(air_file, "w-")
        
        # header dataset
        hdr = f.create_dataset('header',(1,))
        hdr.attrs['code'] = headers['code']
        hdr.attrs['sensors_loc'] = headers['sensors_loc']
        hdr.attrs['starttime'] = headers['starttime']
        hdr.attrs['endtime'] = headers['endtime']
        hdr.attrs['fq_band'] = headers['fq_band']

        hdr.attrs['sampling_rate'] = headers['sampling_rate']
        hdr.attrs['time_width'] = headers['time_width']
        hdr.attrs['overlap'] = headers['overlap']

        # model
        hdr.attrs['radii'] = headers['radii']
        hdr.attrs['vel_air'] = headers['vel_air']
        hdr.attrs['h_src'] = headers['h_src']
        hdr.attrs['src_dgr'] = headers['src_dgr']

        azmtbins = int(np.diff(headers['src_dgr']) + 1) # [-1] - headers['src_dgr'][0]
        
        # intervals
        interval_length = sstep.interval_.total_seconds()/60
        last_interval_length = sstep.last_interval_.total_seconds()/60

        nro_timebins = sstep.total_steps()

        if last_interval_length == interval_length:
            chunk_size = (sstep.steps_*chunk_weight, azmtbins)
        else:
            chunk_size = (max(sstep.steps_, sstep.laststeps_)*chunk_weight, azmtbins)
        
        # if chunk_size[0] > 
        
        hdr.attrs['time_bins'] = nro_timebins

        # print shape info
        print(" hdf5_memory info: %s " % air_file)
        print(' --- all datasets sizes are ', (nro_timebins, azmtbins))
        print(' --- saving in memory chunk sizes of ', chunk_size)
        print('')

        # crosscorr datasets
        crosscorr = f.create_group("crosscorr")
        crosscorr.create_dataset('crosscorr', (nro_timebins, azmtbins), chunks=chunk_size, dtype=np.float32)

        # crosscorr datasets
        press = f.create_group("pressure")
        press.create_dataset('p_avg', (nro_timebins, azmtbins), chunks=chunk_size, dtype=np.float32)
        press.create_dataset('p_max', (nro_timebins, azmtbins), chunks=chunk_size, dtype=np.float32)

        # delta_times (review)
        dts = f.create_group("delta_times")
        for dt_key in delta_times.keys():
            i = dts.create_dataset(dt_key, (len(delta_times[dt_key]),), dtype=np.int8)
            i.attrs['sta'] = dt_key
            i[:] = delta_times[dt_key]

        f.flush()
        f.close()

        return AiR(air_file)