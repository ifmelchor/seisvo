#!/usr/bin/env python3
# coding=utf-8

import time
import h5py
import numpy as np
import multiprocessing as mp 
from ..signal import jl_zlcc

def attr_filt(attr_list, which):
    """
    Filter a list of attributes to only vector or scalar.
    which must be a string of "scalar" or "vector"
    """

    assert which in ("scalar", "vector")
    assert isinstance(attr_list, list)

    filter_list = []

    for attr in attr_list:
        if which == "scalar" and attr in ["slow", "baz", "maac", "rms", "error"]:
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
        
        t0 = time.time()

        cc8_ans = jl_zlcc(data, self.cc8stats.sample_rate, self.headers["utm"]["x"], 
                          self.headers["utm"]["y"], self.cc8stats.fq_bands[fqidx-1], 
                          self.cc8stats.slow_max, self.cc8stats.slow_int, 
                          self.headers["lwin"], nwin, self.headers["nadv"], 
                          self.cc8stats.ccerr_thr, self.headers["toff_sec"], 
                          self.cc8stats.slow2, self.headers["maac_thr"],
                          self.headers["slow_max2"], self.headers["slow_int2"]
                          )
        
        tp = time.time() - t0 # processing time
        
        self.queue.put((self.n, cc8_ans, fqidx, start, end, tp, nwin))
    

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

        # compute slowness space
        nite = 1 + 2*int(self.cc8stats.slow_max/self.cc8stats.slow_int)

        matrix_nan = np.full([nwin, nite, nite], np.nan)
        matrixbnd_nan = np.full([nwin, 2], np.nan)
        cc8_ans = {}

        for attr in ("slow", "baz", "maac", "rms", "error"):
            cc8_ans[attr] = vector_nan
        
        if self.cc8stats.slowmap:
            cc8_ans["slowmap"] = matrix_nan
        
        if self.cc8stats.slow2:
            cc8_ans["slow2"] = vector_nan
            cc8_ans["baz2"] = vector_nan
        
        cc8_ans["slowbnd"] = matrixbnd_nan
        cc8_ans["bazbnd"] = matrixbnd_nan
        
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


