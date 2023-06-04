#!/usr/bin/env python3
# coding=utf-8

import time
import h5py
import numpy as np
import multiprocessing as mp 
from ..signal import get_CC8, get_stats, get_pdf_data

# from .plotting import cc8_plot, slowness_map_motion, simple_slowness_plot, plot_slowbaz_tmap

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
    def __init__(self, cc8stats, xutm, yutm, lwin, nwin, nadv, toff):
        self.cc8stats  = cc8stats
        self.xutm      = xutm
        self.yutm      = yutm
        self.lwin      = lwin
        self.nwin      = nwin
        self.nadv      = nadv
        self.toff      = toff
        self.n         = -1
        self.processes = []
        self.queue     = mp.Queue()


    def _wrapper(self, data, fqband, fqidx, starttime, endtime):
        t0 = time.time()
        cc8_ans = get_CC8(data, self.cc8stats.sample_rate, self.xutm, self.yutm,\
            fqband, self.cc8stats.slow_max, self.cc8stats.slow_inc,\
            lwin=self.lwin, nwin=self.nwin, nadv=self.nadv, cc_thres=self.cc8stats.cc_thres, toff=self.toff)
        t1 = time.time()
        self.queue.put((self.n, cc8_ans, fqidx, starttime, endtime, t1-t0))
    

    def run(self, *args, **kwargs):
        self.n += 1
        p = mp.Process(target=self._wrapper, args=args)
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
        for n, data, nfq, starttime, endtime, proct in rets:
            self.data[n] = {}
            self.data[n]["nfq"]   = nfq
            self.data[n]["ans"]   = data
            self.data[n]["start"] = starttime
            self.data[n]["end"]   = endtime
            self.data[n]["proct"] = proct
        
        for p in self.processes:
            p.join()
        
        self.save(int_prct)


    def save(self, int_prct):
        available_n = list(self.data.keys())
        available_n.sort() # nro interval per job
        
        for n in available_n:
            start   = self.data[n]["start"].strftime("%d %b%y %H:%M:%S")
            end     = self.data[n]["end"].strftime("%d %b%y %H:%M:%S")
            cc8_ans = self.data[n]["ans"]
            nfq     = self.data[n]["nfq"]
            proct   = self.data[n]["proct"]

            if not cc8_ans:
                proc_status = "FAIL"
                vector_nan = np.full([self.nwin,], np.nan)

                cc8_ans = {}
                for nsi in range(1, len(self.cc8stats.nites)+1):
                    nite = self.cc8stats.nites[nsi-1]
                    matrix_nan = np.full([self.nwin, nite, nite], np.nan)
                    matrixbnd_nan = np.full([self.nwin, 2], np.nan)
                    cc8_ans[nsi] = {}

                    for attr in ("slow", "bazm", "maac", "rms"):
                        cc8_ans[nsi][attr] = vector_nan
                    
                    cc8_ans[nsi]["slowmap"] = matrix_nan
                    cc8_ans[nsi]["slowbnd"] = matrixbnd_nan
                    cc8_ans[nsi]["bazmbnd"] = matrixbnd_nan

            else:
                proc_status = "OK"
            
            self.cc8stats.__write__(cc8_ans, nfq, self.nwin)

            print(f" >> {start} -- {end}  ::  data wrote  {proc_status}  ::  job {int_prct} :: fq_band {nfq}/{len(self.cc8stats.fq_bands)}  ::  process time {proct:.1f} sec")

        self.reset()


class CC8out(object):
    def __init__(self, cc8stats, dout):
        self.cc8stats = cc8stats
        self.starttime_ = dout["dtime"][0]
        self.endtime_   = dout["dtime"][-1]
        self._dout = dout
        self.__set_stats__()


    def __set_stats__(self):
        self.attr_list = []
        fqslo_    = []
        _fqidx    = []
        _slidx    = []

        for k in self._dout.keys():
            klist = k.split("/")
            if len(klist) > 1:
                fqslo_.append('/'.join(klist[:-1]))
                _fqidx.append(klist[0])
                _slidx.append(klist[1])

                if klist[2] not in self.attr_list:
                    self.attr_list.append(klist[2])

        self.fqslo_ = list(set(fqslo_))
        self.fqslo_.sort()

        self._fqidx = list(set(_fqidx))
        self._fqidx.sort()

        self._slidx = list(set(_slidx))
        self._slidx.sort()


    def __str__(self):
        txt_to_return =  f'\n   >>LTE file    ::: {self.cc8stats.file}'
        txt_to_return += f'\n   > ID           :  {self.cc8stats.id}'
        txt_to_return += f'\n   >starttime     :  {self.starttime_.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >endtime       :  {self.endtime_.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >attribute     :  {self.attr_list}'
        txt_to_return +=  f'\n'
        return txt_to_return


    def __check_attr__(self, attr, which=None):

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
            attr_list = attr_filt(attr_list, which)
            
        return attr_list


    def get_stats(self, attr_list=None, fq_idx=None, slow_idx=None, db_scale=True):
        """
        Return (min, max, mean, mode) of an specific scalar-only attribute
        """

        if not attr_list:
            attr_list = self.attr_list
        
        attr_list = self.__check_attr__(attr_list, "scalar")
        
        if not fq_idx:
            fq_idx = self._fqidx
        else:
            fq_idx = self.cc8stats.check_idx(fq_idx, "fq")

        if not slow_idx:
            slow_idx = self._slidx
        else:
            slow_idx = self.cc8stats.check_idx(slow_idx, "slow")

        dout = {}

        for nfqi in fq_idx:
            for nsi in slow_idx:
                fqslo = "/".join([nfqi, nsi])

                if fqslo in self.fqslo_:
                    for attr in attr_list:
                        key = "/".join([fqslo, attr])
                        data = self._dout[key]

                        if attr == "bazm":
                            data[data>400] = np.nan

                        if attr == 'rms' and db_scale:
                            data = 10*np.log10(data)
                    
                        dout[key] = get_stats(data[np.isfinite(data)])
        
        return dout


    def get_pdf(self, attr, fq_idx=None, slow_idx=None, bandwidth=0.02, maac_th=None, rel_prob=False, **kwargs):
        
        attr_list = self.__check_attr__(attr, "scalar")
        
        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            fq_idx = str(fq_idx)
            assert fq_idx in self._fqidx
        
        if not slow_idx:
            slow_idx = self._slidx[0]
        else:
            slow_idx = str(slow_idx)
            assert slow_idx in self._slidx

        fqslo = "/".join([fq_idx, slow_idx])
        key = "/".join([fqslo, attr_list[0]])
        data = self._dout[key]
        
        if attr == "bazm":
            data[data>400] = np.nan
        
        if attr == "rms":
            data[data<=0] = np.nan 
            data = 10*np.log10(data)
        
        if isinstance(maac_th, float):
            if "maac" in self.attr_list:
                # filter data to fullfill with maac values
                key = "/".join([fqslo, "maac"])
                maac_data = self._dout[key]
                data = data[maac_data>=maac_th]
            else:
                print("warn :: no maac data found in attr_list")

        # remove nan
        data = data[np.isfinite(data)]
        if data.any():
            y, pdf = get_pdf_data(data, bandwidth, db_scale=False, **kwargs)
            
            if rel_prob:
                pdf *= np.log(data.shape[0])
            
            return y, pdf
        
        else:
            return None, None


    def plot(self, cc_th=0.5, fq_idx=None, slow_idx=None, datetime=False, **kwargs):

        if not fq_idx:
            fq_idx = self._fqidx
        else:
            fq_idx = self.cc8stats.check_idx(fq_idx, "fq")

        if not slow_idx:
            slow_idx = self._slidx
        else:
            slow_idx = self.cc8stats.check_idx(slow_idx, "slow")

        fq_slo_idx = []
        for nfqi in fq_idx:
            for nsi in slow_idx:
                fqslo = "/".join([nfqi, nsi])
                if fqslo in self.fqslo_:
                    fq_slo_idx.append(fqslo)

        ans = cc8_plot(self, fq_slo_idx[0], cc_th=cc_th, datetime=datetime, **kwargs)

        return ans
        
    
    def plot_slowmap(self, fq_idx, slow_idx, fps=30, plot=True, save=False, starttime=None, endtime=None, filename=None):

        assert "slowmap" in self.attr_list

        if isinstance(fq_idx, int):
            fq_idx = str(fq_idx)
        
        if isinstance(slow_idx, int):
            slow_idx = str(slow_idx)

        fq_slo_idx = "/".join([fq_idx, slow_idx])
        motion = slowness_map_motion(self, fq_slo_idx, fps, starttime=starttime, endtime=endtime, plot=plot)

        if save:
            if not filename:
                filename = "_".join(self.cc8.stats.id, fq_slo_idx) + '.mp4'
            else:
                if filename.split('.')[-1] != "mp4":
                    filename += '.mp4'
            
            motion.save(filename, fps=fps, extra_args=['-vcodec', 'libx264'])
    

    def get_maac_probmap(self, fq_idx, slow_idx, starttime=None, endtime=None, **kwargs):

        assert "slowmap" in self.attr_list

        if isinstance(fq_idx, int):
            fq_idx = str(fq_idx)
        
        assert fq_idx in self._fqidx
        
        if isinstance(slow_idx, int):
            slow_idx = str(slow_idx)
        
        assert slow_idx in self._slidx

        key = "/".join([fq_idx, slow_idx, "slowmap"])
        data = self._dout[key]

        # load SAP.jl
        from juliacall import Main as jl
        jl.seval("using SAP")
        slowprob = jl.mpm(jl.Array(data))

        plot = kwargs.get("plot", False)
        fileout = kwargs.get("fileout", None)

        if plot or fileout:
            slomax = self.cc8stats.slow_max[int(slow_idx)-1]
            sloinc = self.cc8stats.slow_inc[int(slow_idx)-1]
            title  = f"\n {self.starttime_} -- {self.endtime_}"
            simple_slowness_plot(slowprob, slomax, sloinc, title=title, bar_label="most probable MAAC", **kwargs)

        return slowprob
    

    def get_slobaz_tmap(self, fq_idx, slow_idx, cc_th, starttime=None, endtime=None, savefig=True, **kwargs):

        assert "slowmap" in self.attr_list

        if isinstance(fq_idx, int):
            fq_idx = str(fq_idx)
        
        assert fq_idx in self._fqidx
        
        if isinstance(slow_idx, int):
            slow_idx = str(slow_idx)
        
        assert slow_idx in self._slidx

        key = "/".join([fq_idx, slow_idx, "slowmap"])
        data = self._dout[key]
        slomax = self.cc8stats.slow_max[int(slow_idx)-1]
        sloinc = self.cc8stats.slow_inc[int(slow_idx)-1]

        from juliacall import Main as jl
        jl.seval("using SAP")

        slowbaz_tmap = jl.slobaztmap(jl.Array(data), sloinc, slomax, cc_th)

        if savefig:
            title = f"MAAC > {cc_th}" + f"\n {self.starttime_} -- {self.endtime_}"
            plot_slowbaz_tmap(slowbaz_tmap, title, **kwargs)

        return slowbaz_tmap




