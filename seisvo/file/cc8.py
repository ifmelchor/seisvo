#!/usr/bin/python3
# coding=utf-8

'''

Read write and operate .lte files

'''

import os
import h5py
import numpy as np
import datetime as dt
import functools
import julia
from tqdm import tqdm
from multiprocessing import Process, Queue, Pool
from obspy.core.util.attribdict import AttribDict

class cc8stats(AttribDict):
    def __init__(self, header):
        super(cc8stats, self).__init__(header)

    def __add_attr__(self, list):
        self.attributes = list

    def __str__(self):
        priorized_keys = ['id','starttime','endtime','window_length','overlap','nro_time_bins','last_time_bin','sample_rate','fq_band','np','pmax','pinc', 'ccerr']

        return self._pretty_str(priorized_keys)


class CC8Process(object):
    def __init__(self, h5file, args):
        self.file_ = h5file
        self.init_args = args
        self.n = -1
        self.processes = []
        self.queue = Queue()

    def _wrapper(self, *args):
        julia.Julia(compiled_modules=False)
        from julia import CC8 as cc8
        ret = cc8.run(
            args[0],
            args[1],
            args[2],
            self.init_args[0],
            self.init_args[1],
            self.init_args[2],
            self.init_args[3],
            self.init_args[4],
            self.init_args[5],
            self.init_args[6],
            self.init_args[7],
            )

        # data[self.n] = {
        data = {
            'ip': [r.ip-1 for r in ret],
            'slows':[r.slowness for r in ret],
            'bazth':[r.bazimuth for r in ret],
            'ccmax':[r.macc for r in ret],
            'ccmap':[r.sumap for r in ret]
        }

        self.queue.put((self.n, data))
    

    def run(self, *args):
        self.n += 1
        p = Process(target=self._wrapper, args=args)
        self.processes.append(p)
        p.start()
    

    def reset(self):
        self.processes = []
        self.queue = Queue()


    def wait(self):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        
        self.data = {}
        for n, data in rets:
            self.data[n] = data
        
        for p in self.processes:
            p.join()
    

    def save(self):
        available_n = list(self.data.keys())
        available_n.sort() # nro interval per job

        for n in available_n:
            data = self.data[n]

            #save into h5file
            with h5py.File(self.file_, "r+") as h5f:
                nbin = h5f['header'].attrs["last_time_bin"]
                nwin = self.init_args[4]
                
                # iterate over slowness domain
                for i in data["ip"]:
                    ip_group = h5f.get(str(i))
                    for attr in ('slows','bazth','ccmax'):
                        ip_group[attr][nbin+1:nbin+1+nwin] = data[attr][i]
                    ip_group['ccmap'][nbin+1:nbin+1+nwin,:,:] = data['ccmap'][i]

                h5f['header'].attrs.modify('last_time_bin', nbin+nwin)
                h5f.flush()

        self.reset()


class CC8(object):
    def __init__(self, filename):
        if not os.path.isfile(filename):
            raise ValueError(' file not found')
        
        self.file_ = filename
        
        try:
            self.__set_stats__()
        except:
            raise ValueError(' file %s cannot be read' %filename)


    def __set_stats__(self):
        with h5py.File(self.file_, "r") as f:
            hdr = f['header']
            self.stats = cc8stats(
                dict(
                id = hdr.attrs['id'],
                starttime = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                endtime = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                window_length = int(hdr.attrs['window_length']),
                overlap = int(hdr.attrs['overlap']),
                nro_time_bins = int(hdr.attrs['nro_time_bins']),
                last_time_bin = int(hdr.attrs['last_time_bin']),
                sample_rate = int(hdr.attrs['sample_rate']),
                fq_band = hdr.attrs['fq_band'],
                np = int(hdr.attrs['np']),
                pmax = list(hdr.attrs['pmax']),
                pinc = list(hdr.attrs['pinc']),
                ccerr = hdr.attrs['ccerr'],
                )
            )

    
    def __compute__(self, sarray, headers):
        if self.stats.last_time_bin > 0:
            raise ValueError("last_time_bin is > 0. Please check cc8 file.")

        # loop over the intervals
        interval = headers['interval']
        nro_intervals = headers['nro_intervals']
        start_time = self.stats.starttime
        end_time = self.stats.endtime
        nini_sec = headers['nini']
        nwin = int(headers['nwin'])
        lwin = int(headers['lwin'])
        fsem = int(self.stats.sample_rate)
        njobs = int(headers['njobs']) # must be < nro_intervals

        # init process
        cc8p = CC8Process(
            self.file_,[
                self.stats.pmax,
                self.stats.pinc,
                fsem, 
                lwin, 
                nwin,
                headers['nadv'],
                self.stats.ccerr,
                int(nini_sec*fsem)
                ]
            )

        start = start_time
        while start + interval <= end_time:
                end = start + interval
                
                # prepare parameters
                mdata, stats = sarray.get_mdata(start, end, toff_in_sec=nini_sec, sample_rate=fsem, fq_band=self.stats.fq_band, return_stats=True)
                xsta_utm = [st.lon for st in stats]
                ysta_utm = [st.lat for st in stats]
               
                # run the process
                cc8p.run(mdata, xsta_utm, ysta_utm)

                if len(cc8p.processes) == njobs:
                    cc8p.wait()
                    cc8p.save()
                
                # advance 
                start += interval
        
        # wait until finishes
        if len(cc8p.processes) > 0:
            cc8p.wait()
            cc8p.save()
    
    def __idx__(self, start_time=None, end_time=None):
        total_interval = (self.stats.endtime - self.stats.starttime).total_seconds()
        if start_time:
            diff = (start_time - self.stats.starttime).total_seconds()
            idxi = int(diff*self.nro_time_bins/total_interval)
            if idxi > self.stats.last_time_bin:
                raise ValueError("start_time > end_time !!")
        else:
            idxi = 0
        
        if end_time:
            diff = (end_time - self.stats.starttime).total_seconds()
            idxf = int(diff*self.nro_time_bins/total_interval) + 1

            if idxf >= self.stats.last_time_bin:
                idxf = -1
        else:
            idxf = -1
        
        # get time
        start = self.stats.starttime
        delta = dt.timedelta(seconds=self.stats.window_length)
        olap = dt.timedelta(seconds=self.stats.overlap)
        time = [start]
        while start + delta <= self.stats.endtime:
            time.append(start + delta)
            start += delta - olap
        
        return time[idxi:idxf], (idxi, idxf)

    def get(self, start_time=None, end_time=None, ip=0, ccmap=False):

        data = {}
        time, (idxi, idxf) = self.__idx__(start_time, end_time)

        with h5py.File(self.file_, "r") as h5f:
            ip_group = h5f.get(str(ip))
            for attr in ('slows','bazth','ccmax'):
                data[attr] = ip_group[attr][idxi:idxf]
            
            if ccmap:
                data['ccmap'] = ip_group['ccmap'][idxi:idxf,:,:]
        
        return time, data
        

    @staticmethod
    def new(sarray, filename, headers):
        
        # print info
        print('')
        print(' CC8 file INFO')
        print(' -------------')
        print('  file   ::  %s ' % filename)
        for key in ['id', 'starttime', 'endtime' ,'window_length' ,'overlap' ,'nro_time_bins' ,'sample_rate' ,'fq_band' , 'np', 'pmax', 'pinc', 'ccerr']:
            print(f'  {key}  ::  {headers[key]}')
        print('')

        # create file
        with h5py.File(filename, "w-") as cc8h5:
        
            # add header
            hdr = cc8h5.create_dataset('header',(1,))
            hdr.attrs['id'] = headers['id']
            hdr.attrs['starttime'] = headers['starttime']
            hdr.attrs['endtime'] = headers['endtime']
            hdr.attrs['window_length'] = headers['window_length']
            hdr.attrs['overlap'] = headers['overlap']
            hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
            hdr.attrs['last_time_bin'] = -1
            hdr.attrs['sample_rate'] = headers['sample_rate']
            hdr.attrs['fq_band'] = headers['fq_band']
            hdr.attrs['np'] = headers['np']
            hdr.attrs['pmax'] = headers['pmax']
            hdr.attrs['pinc'] = headers['pinc']
            hdr.attrs['ccerr'] = headers['ccerr']

            # add datasets
            timebins = headers['nro_time_bins']
            for n in range(headers['np']):
                nite = 2*int(headers["pmax"][n]/headers["pinc"][n]) + 1
                np_n = cc8h5.create_group(str(n))
                np_n.create_dataset('slows', (headers['nro_time_bins'],),          chunks=True, dtype=np.float32)
                np_n.create_dataset('bazth', (headers['nro_time_bins'],),          chunks=True, dtype=np.float32)
                np_n.create_dataset('ccmax', (headers['nro_time_bins'],),          chunks=True, dtype=np.float32)
                np_n.create_dataset('ccmap', (headers['nro_time_bins'],nite,nite), chunks=True, dtype=np.float32)

            cc8h5.flush()

        cc8 = CC8(filename)
        cc8.__compute__(sarray, headers)
        
        return cc8


