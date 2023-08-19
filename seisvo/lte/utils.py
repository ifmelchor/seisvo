#!/usr/bin/env python3
# coding=utf-8

import time
import numpy as np
import multiprocessing as mp
from ..signal import get_LTE
# from .peaks import Peaks
# from ..plotting.lte import ltaoutsta_plot


class _LTEProcess(object):
    def __init__(self, ltestats, lteheaders):
        self.ltestats = ltestats
        self.ltehdr   = lteheaders
        self.queue    = mp.Queue()
        self.n        = -1
        self.processes = []


    def wrapper(self, data, start, end, last):
        kwargs = self.ltestats.get_kwdict()

        if last:
            kwargs["nwin"] = self.ltehdr["last_nwin"]
        else:
            kwargs["nwin"] = self.ltehdr["nwin"]

        kwargs["lwin"] = self.ltehdr["lwin"]
        kwargs["wadv"] = self.ltehdr["wadv"]
        kwargs["nswin"] = self.ltehdr["nswin"]
        kwargs["lswin"] = self.ltehdr["lswin"]
        kwargs["swadv"] = self.ltehdr["swadv"]

        t0 = time.time()
        lte_ans = get_LTE(data, self.ltestats.sample_rate,\
            self.ltestats.get_codelist(), self.ltestats.fq_band, **kwargs)
        proctime = time.time() - t0

        self.queue.put((self.n, lte_ans, start, end, proctime, kwargs["nwin"]))
    

    def run(self, data, start, end, last):
        self.n += 1
        p = mp.Process(target=self.wrapper, args=(data, start, end, last))
        self.processes.append(p)
        p.start()
    

    def get_empty_lte(self, nwin):
        nfb  = self.ltestats.nro_freq_bins

        lte_ans = {}
        if self.ltestats.type == "station":
            nfb  = self.ltestats.nro_freq_bins
            matrix_nan = np.full([nwin, nfb], np.nan)
            vector_nan = np.full([nwin,], np.nan)

            for n, chan in enumerate(self.ltestats.channel):
                lte_ans[chan] = {}
                lte_ans[chan]["specgram"]  = matrix_nan
                for attr in ("perm_entr", "vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar"):
                    lte_ans[chan][attr] = vector_nan

            if self.ltestats.polar:
                lte_ans["polar"] = {}
                for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh"):
                    lte_ans["polar"][attr] = matrix_nan
        
        if self.ltestats.type == "network":
            
            matrix_nan = np.full([nwin, nfb], np.nan)
            station_matnan = np.full([nwin, nfb, len(self.ltestats.stations)], np.nan, dtype=np.complex64)
            
            for sta in self.ltestats.stations:
                lte_ans[sta] = matrix_nan
            
            lte_ans["csw"] = matrix_nan
            lte_ans["vt"]  = station_matnan

        return lte_ans
    

    def reset(self):
        self.processes = []
        self.queue = mp.Queue()


    def wait(self, nint):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        
        self.data = {}
        for n, lte_ans, start, end, proct, nwin in rets:
            self.data[n] = {}
            self.data[n]["lte"]   = lte_ans
            self.data[n]["start"] = start
            self.data[n]["end"]   = end
            self.data[n]["proct"] = proct
            self.data[n]["nwin"]  = nwin
        
        for p in self.processes:
            p.join()
        
        self.save(f"{nint}/{self.ltehdr['nro_intervals']}")
        

    def save(self, int_prct):
        available_n = list(self.data.keys())
        available_n.sort() # nro interval per job

        for n in available_n:
            lte_ans = self.data[n]["lte"]
            start   = self.data[n]["start"].strftime("%d %b%y %H:%M")
            end     = self.data[n]["end"].strftime("%d %b%y %H:%M")
            proct   = self.data[n]["proct"]
            nwin    = self.data[n]["nwin"]

            if not lte_ans:
                proc_status = "FAIL"
                lte_ans = self.get_empty_lte(nwin)
            else:
                proc_status = "OK"

            self.ltestats.__write__(lte_ans, nwin)
            print(f" [interval {int_prct}]  >> data wrote {proc_status} from {start} to {end}  [in {proct:.1f} sec]")

        # clear process
        self.reset()