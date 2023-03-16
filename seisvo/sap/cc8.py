#!/usr/bin/python3
# coding=utf-8

'''

Read write and operate .cc8 files

'''

import os
import pickle
import h5py
import numpy as np
import pandas as pd
import datetime as dt
from .cc8utils import CC8stats, CC8Process, CC8out

# from tqdm import tqdm
ATTR_LIST = ["rms", "maac", "slow", "bazm", "slowmap", "bazmbnd", "slowbnd"]


class CC8(object):
    def __init__(self, cc8_file):
        if not os.path.isfile(cc8_file):
            raise ValueError(' file not found')
        
        self.file_ = cc8_file
        self.__set_stats__()


    def __set_stats__(self):
        stats_dict = {}
        
        with h5py.File(self.file_, "r") as f:
            hdr = f['header']
            stats_dict["id"]              = str(hdr.attrs['id'])
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
        self.nites_ = [1 + 2*int(pmax/pinc) for pmax, pinc in zip(self.stats.slow_max, self.stats.slow_inc)]
        self.fqidx_ = [str(fqix) for fqix in range(1, len(self.stats.fq_bands)+1)]
        self.sidx_  = [str(sidx) for sidx in range(1, len(self.stats.slow_max)+1)]


    def __idx__(self, start_time=None, end_time=None):
        total_interval = (self.stats.endtime - self.stats.starttime).total_seconds()
        if start_time:
            diff = (start_time - self.stats.starttime).total_seconds()
            idxi = int(diff*self.nro_time_bins/total_interval)
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


    def __checkidx__(self, idx, type):

        if isinstance(idx, (str, int)):
            cond1 = type == "slow" and str(idx) in self.sidx_
            cond2 = type == "fq" and str(idx) in self.fqidx_
            
            if cond1 or cond2:
                return [str(idx)]
            else:
                return []
                
        if isinstance(idx, (list, tuple)):
            if type == "fq":
                return [str(ix) for ix in idx if str(ix) in self.fqidx_]
            
            if type == "slow":
                return [str(ix) for ix in idx if str(ix) in self.sidx_]
            
        return []
    

    def __checkattr__(self, attr):

        if isinstance(attr, str):
            if attr in ATTR_LIST:
                return [attr]
            else:
                return []

        if isinstance(attr, (list, tuple)):
            return [at for at in attr if at in ATTR_LIST]
        
        return []


    def __str__(self):
        return self.stats.__str__()


    def get_time(self, starttime=None, endtime=None):

        if not starttime:
            starttime = self.stats.starttime
        
        if not endtime:
            endtime = self.stats.endtime
        
        assert endtime > starttime
        assert starttime >= self.stats.starttime
        assert endtime <= self.stats.endtime

        durmin = (endtime-starttime).total_seconds()/60
        time_list     = np.linspace(0, durmin, self.stats.nro_time_bins)
        datetime_list = pd.date_range(starttime, endtime, periods=self.stats.nro_time_bins).to_pydatetime()

        n0 = np.argmin(np.abs(starttime - datetime_list))
        nf = np.argmin(np.abs(endtime - datetime_list))+1
        
        return time_list, datetime_list, (n0,nf)


    def get(self, attr=None, fq_idx=None, slow_idx=None, starttime=None, endtime=None):

        if not fq_idx:
            fq_idx = self.fqidx_
        else:
            fq_idx = self.__checkidx__(fq_idx, "fq")
        
        if not slow_idx:
            slow_idx = self.sidx_
        else:
            slow_idx = self.__checkidx__(slow_idx, "slow")
        
        if not attr:
            attr = ATTR_LIST
        else:
            attr = self.__checkattr__(attr)
        
        dout = {}
        dout["time"], dout["dtime"], (n0,nf) = self.get_time(starttime=starttime, endtime=endtime)

        for nfq in fq_idx:
            for ns in slow_idx:
                for at in attr:
                    attr_key = '/'.join([nfq, ns, at])
                    ts  = self.__read__(at, nfq, ns, (n0,nf))
                    if ts[ts>0].any():
                        dout[attr_key] = ts

        return CC8out(self, dout)

        
    def __read__(self, attr, fq_idx, slow_idx, nt):

        assert attr in ATTR_LIST

        n0, nf = nt
        with h5py.File(self.file_, "r") as f:
            if attr == 'slowmap':
                ts = f.get(str(fq_idx)).get(str(slow_idx))[attr][n0:nf,:]
            else:
                ts = f.get(str(fq_idx)).get(str(slow_idx))[attr][n0:nf]
        
        return ts


    def __compute__(self, sarray, headers, njobs):
        if any([ltb > 0 for ltb in self.stats.last_time_bin]):
            raise ValueError("last_time_bin is > 0. Please check cc8 file.")

        # init process
        cc8p = CC8Process(self, int(headers['lwin']), int(headers['nwin']), float(headers['nadv']), int(headers["toff_sec"]))

        # loop over the intervals
        start = self.stats.starttime
        interval = dt.timedelta(minutes=headers['interval'])
        nro_ints = headers['nro_intervals']
        nint = 1
        while start + interval <= self.stats.endtime:
            end = start + interval

            for fqn, fqband in enumerate(self.stats.fq_bands):
                # prepare parameters
                mdata, stats = sarray.get_mdata(start, end, toff_sec=int(headers["toff_sec"]), sample_rate=self.stats.sample_rate, fq_band=(), return_stats=True)
                xsta_utm = np.array([st.lon for st in stats])
                ysta_utm = np.array([st.lat for st in stats])
            
                # run the process
                cc8p.run(mdata, xsta_utm, ysta_utm, np.array(fqband), int(fqn+1), start, end)

                if len(cc8p.processes) == njobs:
                    cc8p.wait(f"{nint}/{nro_ints}")
                
            # advance 
            start += interval
            nint  += 1
        
        # wait until finishes
        if len(cc8p.processes) > 0:
            cc8p.wait(f"{nint}/{nro_ints}")


    def __write__(self, cc8dict, fqn, nwin):
        with h5py.File(self.file_, "r+") as h5f:
            ltblist = h5f['header'].attrs["last_time_bin"]

            nbin = ltblist[fqn-1]
            for nsi in range(1, len(self.nites_)+1):
                h5f[str(fqn)][str(nsi)]["slowmap"][nbin+1:nbin+1+nwin,:,:] = cc8dict[nsi]["slowmap"]
                h5f[str(fqn)][str(nsi)]["slowbnd"][nbin+1:nbin+1+nwin,:] = cc8dict[nsi]["slowbnd"]
                h5f[str(fqn)][str(nsi)]["bazmbnd"][nbin+1:nbin+1+nwin,:] = cc8dict[nsi]["bazmbnd"]
                    
                for attr in ("slow", "bazm", "maac", "rms"):
                    h5f[str(fqn)][str(nsi)][attr][nbin+1:nbin+1+nwin] = cc8dict[nsi][attr]

            ltblist[fqn-1] = nbin + nwin
            h5f['header'].attrs.modify('last_time_bin', ltblist)
            h5f.flush()


    @staticmethod
    def __check_process__(sarray, starttime, endtime, verbose=True):
        
        ans = sarray.check_files(starttime, endtime, sarray.sta_code, loc=sarray.locs)

        if not ans:
            if verbose:
                print(" No missing files!")
            return True
        
        else:
            if verbose:
                for md in ans:
                    print(md)
                print(f" Missing files: {len(ans)}")
            return False


    @staticmethod
    def new(sarray, filename, headers, njobs):
        
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
        cc8.__compute__(sarray, headers, njobs)
        
        return cc8


    @staticmethod
    def attr_filt(attr_list, which):
        """
        Filter a list of attributes to only vector or scalar.
        which must be a string of "scalar" or "vector"
        """

        assert which in ("scalar", "vector")
        assert isinstance(attr_list, list)

        filter_list = []

        for attr in attr_list:
            if which == "scalar" and attr in ["slow", "bazm", "maac", "rms"]:
                filter_list.append(attr)
            
            if which == "vector" and attr in ["slowmap"]:
                filter_list.append(attr)
                
        return filter_list