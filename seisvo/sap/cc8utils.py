#!/usr/bin/env python3
# coding=utf-8

import time
import h5py
import numpy as np
import multiprocessing as mp 
from ..signal import get_CC8, get_Stats, get_PDF

# from .plotting import cc8_plot, slowness_map_motion, simple_slowness_plot, plot_slowbaz_tmap

def _error_of_bounds(slowbnd, azmbnd):
    bazb     = (np.pi/180)*azmbnd
    slodiff  = np.abs(slowbnd[:,1] - slowbnd[:,0])
    bazdiff  = np.abs(bazb[:,1] - bazb[:,0])
    error    = (slodiff + bazdiff) / 2
    return error


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
 

class _CC8Process(object):
    def __init__(self, cc8stats, headers):
        self.cc8stats  = cc8stats
        self.headers   = headers
        self.processes = []
        self.queue     = mp.Queue()
        self.n         = -1


    def _wrapper(self, data, start, end, fqidx, last):
        if last:
            nwin = self.headers["last_nwin"]
        else:
            nwin = self.headers["nwin"]

        fqband = self.cc8stats.fq_bands[fqidx-1]
        fs   = self.cc8stats.sample_rate
        xutm = self.headers["utm"]["x"]
        yutm = self.headers["utm"]["y"]
        
        t0 = time.time()
        cc8_ans = get_CC8(data, fs, xutm, yutm, fqband, self.cc8stats.slow_max,\
            self.cc8stats.slow_int, lwin=self.headers["lwin"], nwin=nwin,\
            nadv=self.headers["nadv"], cc_thres=self.cc8stats.cc_thres,\
            toff=self.headers["toff_sec"])
        t1 = time.time()
        
        self.queue.put((self.n, cc8_ans, fqidx, start, end, t1-t0, nwin))
    

    def run(self, data, start, end, fqidx, last):
        self.n += 1
        p = mp.Process(target=self._wrapper, args=(data, start, end, fqidx, last))
        self.processes.append(p)
        p.start()
    

    def reset(self):
        self.processes = []
        self.queue = mp.Queue()


    def wait(self, nint):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        
        self.data = {}
        for n, data, nfq, starttime, endtime, proct, nwin in rets:
            self.data[n] = {}
            self.data[n]["nfq"]   = nfq
            self.data[n]["ans"]   = data
            self.data[n]["start"] = starttime
            self.data[n]["end"]   = endtime
            self.data[n]["proct"] = proct
            self.data[n]["nwin"]  = nwin
        
        for p in self.processes:
            p.join()
        
        self.save(f"{nint}/{self.headers['nro_intervals']}")
    
    
    def get_empty_dict(self, nwin):
        cc8_ans = {}
        vector_nan = np.full([nwin,], np.nan)
        nite = self.cc8stats.nro_slow_bins
        matrix_nan = np.full([nwin, nite, nite], np.nan)
        matrixbnd_nan = np.full([nwin, 2], np.nan)
        cc8_ans = {}

        for attr in ("slow", "bazm", "maac", "rms"):
            cc8_ans[attr] = vector_nan
        
        cc8_ans["slowmap"] = matrix_nan
        cc8_ans["slowbnd"] = matrixbnd_nan
        cc8_ans["bazmbnd"] = matrixbnd_nan
        
        return cc8_ans


    def save(self, int_prct):
        available_n = list(self.data.keys())
        available_n.sort() # nro interval per job
        
        for n in available_n:
            start   = self.data[n]["start"].strftime("%d %b%y %H:%M:%S.%f")
            end     = self.data[n]["end"].strftime("%d %b%y %H:%M:%S.%f")
            cc8_ans = self.data[n]["ans"]
            nfq     = self.data[n]["nfq"]
            proct   = self.data[n]["proct"]
            nwin    = self.data[n]["nwin"]

            if not cc8_ans:
                proc_status = "FAIL"
                cc8_ans = self.get_empty_dict(nwin)
            else:
                proc_status = "OK"
            
            self.cc8stats.__write__(cc8_ans, nfq, nwin)
            print(f" [interval {int_prct} / fq band {nfq}]  >>  data wrote {proc_status} from {start} to {end}  [in {proct:.1f} sec]")

        self.reset()


