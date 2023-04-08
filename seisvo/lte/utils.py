#!/usr/bin/env python3
# coding=utf-8

import scipy
import time as ttime
import numpy as np
import multiprocessing as mp
from obspy.core.util.attribdict import AttribDict
from sklearn.preprocessing import MinMaxScaler

from seisvo.signal import get_tseries_stats, get_pdf_data
from .peaks import Peaks
from .plotting import ltaoutsta_plot


class LTEstats(AttribDict):
    def __init__(self, header):
        super(LTEstats, self).__init__(header)
        self.set_type()


    def add_attr(self, listin):
        self.attributes = listin


    def set_type(self):
        if len(self.id.split('.')) == 1:
            self.type = "network"
        else:
            self.type = "station"


    def __str__(self):
        if self.type == "station":
            priorized_keys = ['id','channel','starttime','endtime','window','subwindow','subwindow_olap','sample_rate','rm_sens','nro_time_bins','last_time_bin','fq_band','nro_freq_bins','polar','opt_params']

        else:
            # network not implemented yet!
            priorized_keys = []

        return self._pretty_str(priorized_keys)


class LTEProc(object):
    def __init__(self, lte, nwin, lwin, nswin, lswin, nadv):
        self.lte   = lte
        self.nwin  = nwin
        self.lwin  = lwin
        self.nswin = nswin
        self.lswin = lswin
        self.nadv  = nadv
        self.queue = mp.Queue()
        self.n  = -1
        self.processes = []


    def sta_wrapper(self, *args):
        start_time = ttime.time()
        lte_ans = {}

        if isinstance(args[0], np.ndarray):
            # load julia function
            from juliacall import Main as jl
            jl.seval("using LTE")

            # prepare data
            jldata = jl.Array(args[0])
            jlchan = "/".join(self.lte.stats.channel)
            jlband = jl.Tuple(self.lte.stats.fq_band)
            # jl = julia.Julia(compiled_modules=False)
            # jl._LegacyJulia__julia.using("LTE")
            # sta_run_func = jl._LegacyJulia__julia.eval("sta_run")

            # exec jl function
            jlans = jl.sta_run(jldata, jlchan, self.lte.stats.sample_rate, self.nwin, self.lwin, self.nswin, self.lswin, self.nadv, jlband, self.lte.stats.time_bandwidth, self.lte.stats.pad, self.lte.stats.opt_params, self.lte.stats.polar, self.lte.stats.PE_order, self.lte.stats.PE_tau, self.lte.stats.opt_twin, self.lte.stats.opt_th)
            try:
                jlans = dict(jlans)
                # convert the jl dict into a python dict
                for chan in self.lte.stats.channel:
                    lte_ans[chan] = {}
                    lte_ans[chan]["specgram"] = np.array(jlans[chan]["specgram"])

                    for attr in ("perm_entr", "energy", "fq_dominant", "fq_centroid"):
                        lte_ans[chan][attr] = np.array(jlans[chan][attr])
                    
                    if self.lte.stats.opt_params:
                        lte_ans["opt"] = {}
                        for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar"):
                            lte_ans["opt"][attr] = np.array(jlans["opt"][attr])
                    
                    if self.lte.stats.polar:
                        lte_ans["polar"] = {}
                        for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh"):
                            lte_ans["polar"][attr] = np.array(jlans["polar"][attr])
            except:
                # you can catch error here and do whatever
                pass

        proc_duration = ttime.time() - start_time
        
        # save data
        self.queue.put((self.n, lte_ans, args[1], args[2], proc_duration))
    

    def net_wrapper(self, *args):
        start_time = ttime.time()
        lte_ans = {}

        if isinstance(args[0], np.ndarray):
            # load julia function
            from juliacall import Main as jl
            jl.seval("using LTE")

            # prepare data
            jldata = jl.Array(args[0])
            jlstat = "/".join(self.lte.stats.stations)
            jlband = jl.Tuple(self.lte.stats.fq_band)
            # jl = julia.Julia(compiled_modules=False)
            # jl._LegacyJulia__julia.using("LTE")
            # sta_run_func = jl._LegacyJulia__julia.eval("sta_run")

            # exec jl function
            jlans = jl.net_run(jldata, jlstat, self.lte.stats.sample_rate, self.nswin, self.lswin, self.nadv, jlband, self.lte.stats.time_bandwidth, self.lte.stats.pad)
            
            try:
                jlans = dict(jlans)
                # convert the jl dict into a python dict
                for sta in self.lte.stats.stations:
                    lte_ans[sta] = {}
                    lte_ans[sta]["specgram"] = np.array(jlans[chan]["specgram"])
                
                lte_ans["csw"] = np.array(jlans["csw"])
                lte_ans["vt"]  = np.array(jlans["vt"])
                
            except:
                # you can catch error here and do whatever
                pass

        proc_duration = ttime.time() - start_time
        
        # save data
        self.queue.put((self.n, lte_ans, args[1], args[2], proc_duration))


    def get_empty_lte(self):
        lte_ans = {}

        if self.lte.type == "station":
            nfb  = self.lte.stats.nro_freq_bins
            matrix_nan = np.full([self.nwin, nfb], np.nan)
            vector_nan = np.full([self.nwin,], np.nan)

            for chan in self.lte.stats.channel:
                lte_ans[chan] = {}
                lte_ans[chan]["specgram"] = matrix_nan

                for attr in ("perm_entr", "energy", "fq_dominant", "fq_centroid"):
                    lte_ans[chan][attr] = vector_nan
            
            if self.lte.stats.opt_params:
                lte_ans["opt"] = {}
                for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar"):
                    lte_ans["opt"][attr] = vector_nan


            if self.lte.stats.polar:
                lte_ans["polar"] = {}
                for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh"):
                    lte_ans["polar"][attr] = matrix_nan
        
        else:
            nfb  = self.lte.stats.nro_freq_bins
            matrix_nan = np.full([self.nwin, nfb], np.nan)

            for sta in self.lte.stats.stations:
                lte_ans[sta] = {}
                lte_ans[sta]["specgram"] = matrix_nan
            
            lte_ans["csw"] = matrix_nan
            lte_ans["vt"]  = np.full([self.nwin, nfb, len(self.lte.stats.stations)], np.nan)

        return lte_ans


    def run(self, *args):
        self.n += 1

        if self.lte.type == "station":
            p = mp.Process(target=self.sta_wrapper, args=args)
        else:
            p = mp.Process(target=self.net_wrapper, args=args)
        
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

            self.lte.__write__(lte_ans, self.nwin)
            print(f" >> {start} -- {end}  ::  data wrote  {proc_status}  ::  job {int_prct} :: process time {proct:.1f} sec")


        # clear process
        self.reset()


class LTEoutSTA(object):
    def __init__(self, lte, dout):
        self.lte = lte
        self.starttime_ = dout["dtime"][0]
        self.endtime_ = dout["dtime"][-1]
        self._dout = dout
        self.__set_stats__()
    

    def __set_stats__(self):
        self.chan_list = []
        self.attr_list = []
        self._chanattr = []

        for k in self._dout.keys():
            if isinstance(self._dout[k], dict):
                self.chan_list.append(k)

                for ck in self._dout[k]:
                    if ck not in self.attr_list:
                        self._chanattr.append(ck)
                        self.attr_list.append(ck)
            
            else:
                if k not in ("freq", "time", "dtime"):
                    self.attr_list.append(k)
        
        self.npts_ = len(self._dout["time"])
        self.npfs_ = len(self._dout["freq"])
            
    
    def __str__(self):
        txt_to_return =  f'\n   >>LTE file    ::: {self.lte.file_}'
        txt_to_return += f'\n   >channels      :  {self.chan_list}'
        txt_to_return += f'\n   >starttime     :  {self.starttime_.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >endtime       :  {self.endtime_.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >attribute     :  {self.attr_list}'
        txt_to_return +=  f'\n'
        return txt_to_return
    
    
    def __check_chan__(self, chan):

        # check chan list
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

                    dout[key] = get_tseries_stats(data[np.isfinite(data)])
            
            else:
                data = self._dout[attr]
                dout[attr] = get_tseries_stats(data[np.isfinite(data)])

        return dout
    

    def get_pdf(self, vector_attr, chan=None, bandwidth=0.01, y_min=None, y_max=None):

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
    

    def extract(self, chan, peak_distance=5, peak_thresholds={}, shape_min=2, **kwargs):
        """

        Extract dominant peaks as describes in Melchor et al. 2022 (https://doi.org/10.1016/j.jsames.2022.103961)

        Parameters
        ----------
        fq_range : tuple, optional
            by default stats.fq_band
        peak_thresholds : dict, optional
            by default {'fq_delta': 0.05, 'specgram_th': 0.95, 'pd_th': 0.8, 'rect_th': 0.7, 'pd_std':0.1, 'r_std':0.1, 'azm_std':10, 'dip_std':10}

        Returns
        -------
        [Peak object]

        """
        
        if chan not in self.chan_list:
            return

        # check that VECTOR_PARAMS are loaded
        if not all([attr in self.attr_list for attr in ['specgram', 'degree', 'elevation', 'rect', 'azimuth']]):
            return

        # define default peak thresholds and include user's definitions
        self.peak_thresholds = {
            'fq_delta': 0.05,
            'specgram_th': 0.7,
            'degree_th': 0.7,
            'rect_th': 0.7,
            'degree_std':0.1,
            'rect_std':0.1,
            'azimuth_std':10,
            'elevation_std':10
        }
        for key in default_peak_thresholds.keys():
            if peak_thresholds.get(key, None):
                self.peak_thresholds[key] = peak_thresholds.get(key)
                
        # plot model parameters
        if kwargs.get("verbose", True):
            print(f'\n  File: {os.path.basename(self.file_)}')
            print('  ----- Model param ------')
            for key, item in self.peak_thresholds.items():
                print(f'  {key:>10}     ', item)
            print('  ------------------------\n')

        fq_band = kwargs.get("fq_range", self.lte.stats.fq_band)
        # get fq_npts
        fq0 = np.argmin(np.abs(self._dout["freq"] - fq_band[0]))
        fq1 = np.argmin(np.abs(self._dout["freq"] - fq_band[1])) + 1
        
        # get attributes
        sxx = 10*np.log10(self._dout["specgram"])
        pd = self._dout["degree"]
        rect = self._dout['rect']
        thetaH = self._dout['azimuth']
        thetaV = self._dout['elevation']

        # correct 180 ambiguity for azimuth

        peaks = {}
        nro_Wpeaks = 0
        nro_Lpeaks = 0
        
        total_time = len(self._dout["time"])
        for t in range(total_time):
            sxx_t = sxx[t, fq0:fq1]
            pd_t = pd[t, fq0:fq1]

            if np.isfinite(sxx_t).all():
                # prepare to save data
                peaks[t] = dict(fq=[], specgram=[], degree=[], rect=[], azimuth=[], elevation=[])

                # normalize sxx and find peaks
                sxx_t_norm = MinMaxScaler().fit_transform(sxx_t.reshape(-1,1))
                peaks_pos,_ = scipy.signal.find_peaks(sxx_t_norm.reshape(-1,), height=self.peak_thresholds['specgram_th'], distance=peak_distance)
                
                # iterate over the founded peaks
                for p in peaks_pos:
                    fp = freq[fq0:fq1][p]

                    # define frequency range around the peak's frequency                       
                    fp_0 = fp - self.peak_thresholds['fq_delta']
                    fp_1 = fp + self.peak_thresholds['fq_delta']
                    fq_0_pos = np.argmin(np.abs(freq-fp_0))
                    fq_1_pos = np.argmin(np.abs(freq-fp_1)) + 1

                    # compute the averages of sxx and pd and pd's standard devation
                    sp_peak_avg = sxx_t[fq_0_pos:fq_1_pos].mean()
                    pd_peak = pd_t[fq_0_pos:fq_1_pos]
                    pd_peak_avg = pd_peak.mean()
                    pd_peak_std = pd_peak.std()

                    # check whether the peak is well-polarized and save it
                    if pd_peak_avg > self.peak_thresholds['degree_th'] and pd_peak_std < self.peak_thresholds['degree_std']:
                        peaks[t]['fq'] += [fp]
                        peaks[t]['specgram'] += [sp_peak_avg]
                        peaks[t]['degree'] += [pd_peak_avg]
                        nro_Wpeaks += 1

                        # prepare polar data
                        rect_peak = rect[t, fq_0_pos:fq_1_pos]
                        tH_peak = thetaH[t, fq_0_pos:fq_1_pos]
                        tV_peak = thetaV[t, fq_0_pos:fq_1_pos]
                            
                        # remove low pd values (below threshold)
                        rect_peak = rect_peak[np.where(pd_peak > self.peak_thresholds['degree_th'])]
                        tH_peak = tH_peak[np.where(pd_peak > self.peak_thresholds['degree_th'])]
                        tV_peak = tV_peak[np.where(pd_peak > self.peak_thresholds['degree_th'])]

                        # define None as default
                        rect_peak_val = None
                        tH_peak_val = None
                        tV_peak_val = None
                        
                        # filtering
                        if rect_peak.any() and rect_peak.shape[0] >= shape_min:
                            rect_peak_avg = rect_peak.mean()
                            rect_peak_std = rect_peak.std()

                            # check whether the rectilinearity of the peak is statistically robust
                            if rect_peak_std < self.peak_thresholds['rect_std']:
                                rect_peak_val = rect_peak_avg

                            # check whether the peak is linearly-polarized and compute averages and stds
                            if rect_peak_avg > self.peak_thresholds['rect_th']:
                                tH_peak = tH_peak[np.where(rect_peak > self.peak_thresholds['rect_th'])]
                                tV_peak = tV_peak[np.where(rect_peak > self.peak_thresholds['rect_th'])]
                                nro_Lpeaks += 1

                                if tH_peak.any() and tH_peak.shape[0] >= shape_min:
                                    tH_peak_rad = tH_peak*np.pi/180
                                    tH_peak_avg = stats.circmean(tH_peak_rad, high=np.pi)*180/np.pi
                                    tH_peak_std = stats.circstd(tH_peak_rad, high=np.pi)*180/np.pi

                                    # check whether the azimuth of the peak is statistically robust
                                    if tH_peak_std <= self.peak_thresholds['azimuth_std']:
                                        tH_peak_val = tH_peak_avg

                                if tV_peak.any() and tH_peak.shape[0] >= shape_min:
                                    tV_peak_avg = tV_peak.mean()
                                    tV_peak_std = tV_peak.std()

                                    # check whether the elevation of the peak is statistically robust
                                    if tV_peak_std <= self.peak_thresholds['elevation_std']:
                                        tV_peak_val = tV_peak_avg

                        # save results
                        peaks[t]['rect'] += [rect_peak_val]
                        peaks[t]['azimuth'] += [tH_peak_val]
                        peaks[t]['elevation'] += [tV_peak_val]
        
        return Peaks(self, chan, peaks, (nro_Wpeaks, nro_Lpeaks), fq_band)


    def plot(self, chan, attr, **kwargs):
        
        chan_list = self.__check_chan__(chan)
        attr_list = self.__check_attr__(attr)

        fig, _ = ltaoutsta_plot(self, chan_list, attr_list, plot=True, **kwargs)

        return fig