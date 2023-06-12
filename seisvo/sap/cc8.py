#!/usr/bin/python3
# coding=utf-8

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
   

def _new_CC8(array, cc8file, headers, njobs):
    """
        Create new CC8 (hdf5) file
    """

    with h5py.File(cc8file, "w-") as cc8h5:
        # add header
        hdr = cc8h5.create_dataset('header',(1,))
        hdr.attrs['id']            = headers['id']
        hdr.attrs['locs']          = headers['locs']
        hdr.attrs['starttime']     = headers['starttime']
        hdr.attrs['endtime']       = headers['endtime']
        hdr.attrs['window']        = headers['window']
        hdr.attrs['overlap']       = headers['overlap']
        hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
        hdr.attrs['nro_slow_bins'] = headers['nro_slow_bins']
        hdr.attrs['last_time_bin'] = [-1]*len(headers['fq_bands'])
        hdr.attrs['sample_rate']   = headers['sample_rate']
        hdr.attrs['fq_bands']      = headers['fq_bands']
        hdr.attrs['slow_max']      = headers['slow_max']
        hdr.attrs['slow_inc']      = headers['slow_inc']
        hdr.attrs['cc_thres']      = headers['cc_thres']

        # print info
        print('')
        print(' CC8 file INFO')
        print(' ----------------')
        print(" hdf5_memory info: %s " % cc8file)
        print(' dataset size: ', (headers['nro_time_bins'],))
        print('')
        print('    CC8 stats   ')
        print(' ----------------')
        
        for info_key in ['id', 'locs', 'starttime', 'endtime', 'window',\
            'overlap', 'nro_time_bins', 'sample_rate', 'fq_bands',\
            'slow_max', 'slow_inc', 'cc_thres']:
            if info_key == "window":
                print(f'   {info_key} [sec]:  {headers[info_key]}')
            else:
                print(f'   {info_key}:  {headers[info_key]}')        
        print(' ----------------\n')

        # add datasets
        timebins = headers['nro_time_bins']
        for fq_n in range(1, len(headers['fq_bands'])+1):
            fq_n = cc8h5.create_group(str(fq_n))
            
            for sn, nite in enumerate(headers['slow_bins']):
                np_n = fq_n.create_group(str(sn+1))
                
                for attr in ("slow", "bazm", "maac", "rms"):
                    np_n.create_dataset(attr, (timebins,), chunks=True, dtype=np.float32)
                
                np_n.create_dataset('slowmap', (timebins, nite, nite), chunks=True, dtype=np.float32)
                np_n.create_dataset('slowbnd', (timebins, 2), chunks=True, dtype=np.float32)
                np_n.create_dataset('bazmbnd', (timebins, 2), chunks=True, dtype=np.float32)

        cc8h5.flush()
    
    cc8 = CC8(filename)
    cc8.__compute__(array, headers, njobs)

    return cc8



class CC8(object):
    def __init__(self, cc8_file):
        assert os.path.isfile(cc8_file)
    
        stats_dict = {"file":cc8_file}
        with h5py.File(cc8_file, "r") as f:
            hdr = f['header']
            stats_dict["id"]              = str(hdr.attrs['id'])
            stats_dict["locs"]            = list(hdr.attrs['locs'])
            stats_dict["starttime"]       = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S.%f')
            stats_dict["endtime"]         = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S.%f')
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
        """
        Process CC8 file
        """

        # init process
        cc8p = _CC8Process(self.stats, headers)
        excluded_locs = array.get_excluded_locs(self.stats.locs)

        # define timedelta parameters
        toff_sec = dt.timedelta(seconds=headers["toff_sec"])
        start = self.stats.starttime
        delta = dt.timedelta(seconds=(headers["interval"]*60)+headers["int_extra_sec"])

        # loop over intervals
        for nint in range(1, headers["nro_intervals"]+1):
            end = start + delta    
            stream = array.get_stream(start, end, toff_sec=headers["toff_sec"], exclude_locs=excluded_locs)

            # send to a cpu
            if nint == headers["nro_intervals"]:
                last = True
            else:
                last = False

            if stream and stream.get_bounds() == (start-toff_delta, end+toff_delta):
                data = stream.to_array(detrend=True)
            else:
                data = None
            
            for fqn, fqband in enumerate(self.stats.fq_bands):
                cc8p.run(data, start, end, int(fqn+1), last)

                if len(cc8p.processes) == njobs:
                    cc8p.wait(nint)
            
            # advance 
            start = end
        
        # check if there are any process running
        if len(cc8p.processes) > 0:
            cc8p.wait(nint)


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
