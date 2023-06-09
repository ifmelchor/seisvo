#!/usr/bin/python3
# coding=utf-8

'''

Read and create .cc8 files

'''

import os
import h5py
import numpy as np
import pandas as pd
import datetime as dt
from .cc8utils import CC8out, _CC8Process
from ..stats import CC8stats

# from tqdm import tqdm
ATTR_LIST = ["rms", "maac", "slow", "bazm", "slowmap", "bazmbnd", "slowbnd"]

def check_cc8_attr(attr):
    if isinstance(attr, str):
        if attr in ATTR_LIST:
            return [attr]
        else:
            return []

    if isinstance(attr, (list, tuple)):
        return [at for at in attr if at in ATTR_LIST]
   

class CC8(object):
    def __init__(self, cc8_file):
        assert os.path.isfile(cc8_file)
    
        stats_dict = {"file":cc8_file}
        with h5py.File(cc8_file, "r") as f:
            hdr = f['header']
            stats_dict["id"]              = str(hdr.attrs['id'])
            stats_dict["locs"]            = list(hdr.attrs['locs'])
            stats_dict["starttime"]       = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S')
            stats_dict["endtime"]         = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S')
            stats_dict["window"]          = int(hdr.attrs['window'])
            stats_dict["overlap"]         = float(hdr.attrs['overlap'])
            stats_dict["nro_time_bins"]   = int(hdr.attrs['nro_time_bins'])
            stats_dict["last_time_bin"]   = list(hdr.attrs['last_time_bin'])
            stats_dict["sample_rate"]     = int(hdr.attrs['sample_rate'])
            stats_dict["fq_bands"]        = [(float(fqb[0]), float(fqb[1])) for fqb in hdr.attrs['fq_bands']]
            stats_dict["slow_max"]        = [float(smax) for smax in hdr.attrs['slow_max']]
            stats_dict["slow_inc"]        = [float(sinc) for sinc in hdr.attrs['slow_inc']]
            stats_dict["cc_thres"]        = float(hdr.attrs['cc_thres'])

        self.stats  = CC8stats(stats_dict)


    def __str__(self):
        return self.stats.__str__()


    def __read__(self, attr, fq_idx, slow_idx, nt):
        n0, nf = nt
        with h5py.File(self.stats.file, "r") as f:
            if attr == 'slowmap':
                ts = f.get(str(fq_idx)).get(str(slow_idx))[attr][n0:nf,:]
            else:
                ts = f.get(str(fq_idx)).get(str(slow_idx))[attr][n0:nf]
        
        return ts


    def __compute__(self, array, headers, njobs):
        if any([ltb > 0 for ltb in self.stats.last_time_bin]):
            raise ValueError("last_time_bin is > 0. Please check cc8 file.")

        excluded_locs = array.get_excluded_locs(self.stats.locs)
        xutm = np.array([array.utm[loc]["easting"] for loc in self.stats.locs])
        yutm = np.array([array.utm[loc]["northing"] for loc in self.stats.locs])

        # init process
        cc8p = _CC8Process(self.stats, xutm, yutm,\
            int(headers['lwin']), int(headers['nwin']),\
            float(headers['nadv']), int(headers["toff_sec"]))

        start = self.stats.starttime
        interval = dt.timedelta(minutes=headers['interval'])
        nro_ints, nint = headers['nro_intervals'], 1
        toff_delta = dt.timedelta(seconds=int(headers["toff_sec"]))

        # loop over the intervals
        while start + interval <= self.stats.endtime:
            end = start + interval
            stream = array.get_stream(start, end, toff_sec=int(headers["toff_sec"]), exclude_locs=excluded_locs)
            
            # check that stream times ares equal to start and end
            if stream and stream.get_bounds() == (start-toff_delta, end+toff_delta):
                data   = stream.to_array(detrend=True)
            else:
                data = None

            for fqn, fqband in enumerate(self.stats.fq_bands):
                cc8p.run(data, np.array(fqband), int(fqn+1), start, end)

                if len(cc8p.processes) == njobs:
                    cc8p.wait(f"{nint}/{nro_ints}")
                
            # advance 
            start += interval
            nint  += 1
        
        # wait until finishes
        if len(cc8p.processes) > 0:
            cc8p.wait(f"{nint}/{nro_ints}")


    def get_time(self, starttime=None, endtime=None):

        if not starttime:
            starttime = self.stats.starttime
        
        if not endtime:
            endtime = self.stats.endtime
        
        full_interval = (self.stats.starttime, self.stats.endtime)
        interval = (starttime, endtime)
        
        return get_time(full_interval, interval, self.stats.window, self.stats.window_olap)


    def get(self, attr=None, fq_idx=None, slow_idx=None, starttime=None, endtime=None):
        if not fq_idx:
            fq_idx = self.stats.fqidx
        else:
            fq_idx = self.stats.check_idx(fq_idx, "fq")
        
        if not slow_idx:
            slow_idx = self.stats.sidx
        else:
            slow_idx = self.stats.check_idx(slow_idx, "slow")
        
        if not attr:
            attr = ATTR_LIST
        else:
            attr = check_cc8_attr(attr)
        
        dout = {}
        dout["time"], dout["dtime"], (n0,nf) = self.get_time(starttime=starttime, endtime=endtime)

        for nfq in fq_idx:
            for ns in slow_idx:
                for at in attr:
                    attr_key = '/'.join([nfq, ns, at])
                    ts  = _CC8_read(self.stats, at, nfq, ns, (n0,nf))
                    if ts[ts>0].any():
                        dout[attr_key] = ts

        return CC8out(self.stats, dout)

    
    @staticmethod
    def new(array, filename, headers, njobs):
        # print info
        print('')
        print(' CC8 file INFO')
        print(' -------------')
        print('  file   ::  %s ' % filename)
        for key in ['id', 'locs', 'starttime', 'endtime', 'window', 'overlap', 'nro_time_bins', 'sample_rate', 'fq_bands', 'slow_max', 'slow_inc', 'cc_thres']:
            print(f'  {key}  ::  {headers[key]}')
        print('')

        timebins = int(headers['nro_time_bins'])
        nites = [1 + 2*int(pmax/pinc) for pmax, pinc in zip(headers["slow_max"], headers["slow_inc"])]

        # create file
        with h5py.File(filename, "w-") as cc8h5:
            # add header
            hdr = cc8h5.create_dataset('header',(1,))
            hdr.attrs['id']            = headers['id']
            hdr.attrs['locs']          = headers['locs']
            hdr.attrs['starttime']     = headers['starttime']
            hdr.attrs['endtime']       = headers['endtime']
            hdr.attrs['window']        = headers['window']
            hdr.attrs['overlap']       = headers['overlap']
            hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
            hdr.attrs['last_time_bin'] = [-1]*len(headers['fq_bands'])
            hdr.attrs['sample_rate']   = headers['sample_rate']
            hdr.attrs['fq_bands']      = headers['fq_bands']
            hdr.attrs['slow_max']      = headers['slow_max']
            hdr.attrs['slow_inc']      = headers['slow_inc']
            hdr.attrs['cc_thres']      = headers['cc_thres']

            # add datasets
            for fq_n in range(1, len(headers['fq_bands'])+1):
                fq_n = cc8h5.create_group(str(fq_n))
                
                for s_n in range(1, len(nites)+1):
                    nite = nites[s_n-1]
                    np_n = fq_n.create_group(str(s_n))
                    
                    for attr in ("slow", "bazm", "maac", "rms"):
                        np_n.create_dataset(attr, (timebins,), chunks=True, dtype=np.float32)
                    
                    np_n.create_dataset('slowmap', (timebins, nite, nite), chunks=True, dtype=np.float32)
                    np_n.create_dataset('slowbnd', (timebins, 2), chunks=True, dtype=np.float32)
                    np_n.create_dataset('bazmbnd', (timebins, 2), chunks=True, dtype=np.float32)

            cc8h5.flush()
        
        cc8 = CC8(filename)
        cc8.__compute__(array, headers, njobs)

        return cc8

