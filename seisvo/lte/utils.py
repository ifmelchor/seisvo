#!/usr/bin/env python3
# coding=utf-8

import scipy
import time as ttime
import numpy as np
import multiprocessing as mp
from sklearn.preprocessing import MinMaxScaler
from ..signal import get_stats, get_pdf_data, get_LTE
# from .peaks import Peaks
from ..plotting.lte import ltaoutsta_plot


class _LTEProcess(object):
    def __init__(self, ltestats, nwin, lwin, nswin, lswin, nadv):
        self.ltestats = ltestats
        self.nwin     = nwin
        self.lwin     = lwin
        self.nswin    = nswin
        self.lswin    = lswin
        self.nadv     = nadv
        self.queue    = mp.Queue()
        self.n        = -1
        self.processes = []


    def wrapper(self, data, start, end):
        t0 = time.time()
        lte_ans = get_LTE(data, self.ltestats.sample_rate,\
            self.ltestats.get_codelist(), self.ltestats.fq_band,\
            self.nwin, self.lwin, self.nswin, self.lswin,\
            self.nadv, self.ltestats.get_dict())
        proctime = time.time() - t0

        self.queue.put((self.n, lte_ans, start, end, proctime))
    

    def get_empty_lte(self):
        lte_ans = {}

        if self.ltestats.type == "station":
            nfb  = self.ltestats.nro_freq_bins
            matrix_nan = np.full([self.nwin, nfb], np.nan)
            vector_nan = np.full([self.nwin,], np.nan)

            for chan in self.ltestats.channel:
                lte_ans[chan] = {}
                lte_ans[chan]["specgram"] = matrix_nan

                for attr in ("perm_entr", "energy", "fq_dominant", "fq_centroid"):
                    lte_ans[chan][attr] = vector_nan
            
            if self.ltestats.opt_params:
                lte_ans["opt"] = {}
                for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar"):
                    lte_ans["opt"][attr] = vector_nan

            if self.ltestats.polar:
                lte_ans["polar"] = {}
                for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh"):
                    lte_ans["polar"][attr] = matrix_nan
        
        if self.ltestats.type == "network":
            nfb  = self.ltestats.nro_freq_bins
            matrix_nan = np.full([self.nwin, nfb], np.nan)
            station_matnan = np.full([self.nwin, nfb, len(self.ltestats.stations)], np.nan, dtype=np.complex64)
            
            for sta in self.ltestats.stations:
                lte_ans[sta] = matrix_nan
            
            lte_ans["csw"] = matrix_nan
            lte_ans["vt"]  = station_matnan

        return lte_ans


    def run(self, data, start, end):
        self.n += 1
        p = mp.Process(target=self.wrapper, args=(data, start, end))
        self.processes.append(p)
        p.start()
    

    def reset(self):
        self.processes = []
        self.queue = mp.Queue()


    def wait(self, int_prct):
        rets = []
        for p in self.processes:
            ret = self.queue.get()
            rets.append(ret)
        
        self.data = {}
        for n, lte_ans, start, end, proct in rets:
            self.data[n] = {}
            self.data[n]["lte"]   = lte_ans
            self.data[n]["start"] = start
            self.data[n]["end"]   = end
            self.data[n]["proct"] = proct
        
        for p in self.processes:
            p.join()
        
        self.save(int_prct)
        

    def save(self, int_prct):
        available_n = list(self.data.keys())
        available_n.sort() # nro interval per job

        for n in available_n:
            lte_ans = self.data[n]["lte"]
            start   = self.data[n]["start"].strftime("%d %b%y %H:%M")
            end     = self.data[n]["end"].strftime("%d %b%y %H:%M")
            proct   = self.data[n]["proct"]

            # if lte_ans is none, create NaN dictionary
            if not lte_ans:
                proc_status = "FAIL"
                lte_ans = self.get_empty_lte()
            else:
                proc_status = "OK"

            self.ltestats.__write__(lte_ans, self.nwin)
            print(f" >> {start} -- {end}  ::  data wrote  {proc_status}  ::  job {int_prct} :: process time {proct:.1f} sec")

        # clear process
        self.reset()


def _out_stats(dout, ltetype):
    chan_list = []
    attr_list = []

    if ltetype == "station":
        attr_in_chan_list = []
        for key, item in dout.items():
            if isinstance(item, dict):
                chan_list.append(key)
                for ck in item:
                    if ck not in attr_list:
                        attr_list.append(ck)
                        attr_in_chan_list.append(ck)
            else:
                if k not in ("freq", "time", "dtime"):
                    attr_list.append(k)

        return {chan_list, attr_list, attr_in_chan_list}
    

    if ltetype == "network":
        # for key, item in dout.items():
        return None


class LTEout(object):
    """
    this object allow 
    """
    def __init__(self, ltestats, dout):
        self.ltestats = ltestats
        self.starttime = dout["dtime"][0]
        self.endtime = dout["dtime"][-1]
        self._dout = dout

        # define stats

        dstats = _out_stats(dout, ltestats.type)
        self.npts = len(dout["time"])
        self.npfs = len(dout["freq"])
    

    def __str__(self):
        txt_to_return =  f'\n   >>LTE file    ::: {self.ltestats.file}'
        txt_to_return += f'\n   >channels      :  {self.chan_list}'
        txt_to_return += f'\n   >starttime     :  {self.starttime_.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >endtime       :  {self.endtime_.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >attribute     :  {self.attr_list}'
        txt_to_return +=  f'\n'
        return txt_to_return
    
    
    def check_chan(self, chan):
        # check chan list
        assert ltestats.type == ""
        assert isinstance(chan, (list, tuple, str))

        if isinstance(chan, str):
            if chan in self.chan_list:
                chan_list = [chan]
            else:
                print("channel %s not found" % chan)
                return None
        else:
            chan_list = []
            for ch in chan:
                if ch in self.chan_list:
                    chan_list.append(ch)
                else:
                    print("channel %s not found" % ch)
            
            if not chan_list:
                return None

        return chan_list
    

    def check_attr(self, attr, which=None):
        # check attr list
        assert isinstance(attr, (list, tuple, str))
        
        if isinstance(attr, str):
            if attr in self.attr_list:
                attr_list = [attr]
            else:
                print("attribute %s not found" % attr)
                return None
        
        else:
            attr_list = []
            for at in attr:
                if at in self.attr_list:
                    attr_list.append(at)
                else:
                    print("attribute %s not found" % at)
            
            if not attr_list:
                return None
        
        if which:
            attr_list = self.lte.attr_filt(attr_list, which)
            
        return attr_list


    def get_stats(self, attr=None, chan=None):
        """
        Return (min, max, mean, mode) of an specific scalar-only attribute
        """

        # check attr
        if not attr:
            attr = self.attr_list
        attr_list = self.__check_attr__(attr, "scalar")
        
        # check chan
        if not chan:
            chan = self.chan_list
        chan_list = self.__check_chan__(chan)

        dout = {}

        for attr in attr_list:
            if attr in self._chanattr:
                for chan in chan_list:
                    key = "/".join([chan, attr])
                    data = self._dout[chan][attr]

                    dout[key] = get_stats(data[np.isfinite(data)])
            
            else:
                data = self._dout[attr]
                dout[attr] = get_stats(data[np.isfinite(data)])

        return dout
    

    def get_pdf(self, vector_attr, chan=None, bandwidth=None, y_min=None, y_max=None):
        chan_list = self.__check_chan__(chan)

        assert isinstance(vector_attr, str)

        if not self.lte.any_vector([vector_attr]):
            return None
        
        freq = self._dout["freq"]

        if vector_attr in self._chanattr:
            dout = {}
            for chan in chan_list:
                key = "/".join([chan, vector_attr])
                data = self._dout[chan][vector_attr]
                dout[key] = get_pdf_data(data, bandwidth, xmin=y_min, xmax=y_max, db_scale=True)
                
        else:
            data = self._dout[vector_attr]
            dout = get_pdf_data(data, bandwidth, xmin=y_min, xmax=y_max)
        
        return freq, dout


    def plot(self, chan, attr, **kwargs):
        
        chan_list = self.__check_chan__(chan)
        attr_list = self.__check_attr__(attr)

        fig, _ = ltaoutsta_plot(self, chan_list, attr_list, plot=True, **kwargs)

        return fig