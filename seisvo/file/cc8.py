#!/usr/bin/python3
# coding=utf-8

'''

Read write and operate .lte files

'''


import os
import h5py
import numpy as np

from obspy.core.util.attribdict import AttribDict


class cc8stats(AttribDict):
    def __init__(self, header):
        super(cc8stats, self).__init__(header)

    def __add_attr__(self, list):
        self.attributes = list

    def __str__(self):
        priorized_keys = ['id','starttime','endtime','window_length','overlap','nro_time_bins','last_time_bin','sample_rate','fq_band','pmax','pinc','np', 'nini', 'nwin', 'nadv', 'ccerr','sup_file']

        return self._pretty_str(priorized_keys)


class CC8(object):
    def __init__(self, filename, sup_path=None):

        if not os.path.isfile(filename):
            return ValueError(' file not found')

        self.file_ = filename

        try:
            self.__set_stats__()
        except:
            return ValueError(' file %s cannot be read' %filename)

        # search for supfiles
        if sup_path and self.stats.sup_file:
            
            sup_file = os.path.join(cc8_sup_path, )

            if os.path.isfile(sup_file):
            # try to read


    def __set_stats__(self):
        with h5py.File(self.file_, "r") as f:
            hdr = f['header']
            lte_stats = dict(
                id = hdr.attrs['id'],
                starttime = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                endtime = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                ...
                nro_time_bins = hdr.attrs['nro_time_bins'],
                last_time_bin = hdr.attrs['last_time_bin'],
                fq_band = hdr.attrs['fq_band'],
                sample_rate = hdr.attrs['sample_rate'],
                ...
                )
        
        self.stats = cc8stats(lte_stats)
    

    @staticmethod
    def __new__(filename, headers, sup_file=False, sup_file_path=None, kwargs**):

        cc8main = h5py.File(filename, "w-")
        
        # header dataset
        hdr = cc8main.create_dataset('header',(1,))
        hdr.attrs['id'] = headers['id']
        hdr.attrs['starttime'] = headers['starttime']
        hdr.attrs['endtime'] = headers['endtime']
        hdr.attrs[''] = headers['']
        hdr.attrs['sample_rate'] = headers['sample_rate']
        hdr.attrs['fq_band'] = headers['fq_band']
        hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
        hdr.attrs['time_bandwidth'] = headers['time_bandwidth']
        hdr.attrs['last_time_bin'] = -1
        timebins = headers['nro_time_bins']

        # set chunks 
        chunk_shape1 = chunk_shape2 = True
        chunk_info = 'auto'