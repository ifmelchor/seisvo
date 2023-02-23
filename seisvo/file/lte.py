#!/usr/bin/python3
# coding=utf-8

'''

Read write and operate .lte files

'''

import os
import h5py
import numpy as np
import pandas as pd
import datetime as dt
import multiprocessing

from itertools import chain
from functools import partial
from tqdm import tqdm
from scipy import signal, stats
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from obspy.core.util.attribdict import AttribDict
from pyentrp import entropy as ent

import time as ttime
from seisvo.signal.polarization import PolarAnalysis
from seisvo.signal.proba import get_PDF, get_KDE


SCALAR_PARAMS = ['fq_dominant', 'fq_centroid', 'energy', 'pentropy', 'rsam', 'mf', 'hf', 'vlf', 'lf', 'vlar', 'lrar', 'rmar', 'dsar']
VECTORAL_PARAMS = ['specgram', 'degree', 'elevation', 'rect', 'azimuth']
AMP_PARAMS = ['pentropy', 'rsam', 'dsar', 'mf', 'hf', 'lf', 'vlf', 'vlar', 'lrar', 'rmar']
SPEC_PARAMS = ['specgram', 'fq_dominant', 'fq_centroid', 'energy']
POLAR_PARAMS = ['degree', 'elevation', 'rect', 'azimuth']
LTE_GROUPS = ['amplitude', 'spectral', 'polar']

def get_group(attr):
    if attr in AMP_PARAMS:
        return 'amplitude'

    if attr in SPEC_PARAMS:
        return 'spectral'

    if attr in POLAR_PARAMS:
        return 'polar'


default_bandwidth = {
            'freq': 0.01,
            'specgram': 0.05,
            'rect': 0.25,
            'azimuth': 0.5,
            'elevation': 0.5
        }


default_peak_thresholds = {
    'fq_delta': 0.05,
    'specgram_th': 0.7,
    'degree_th': 0.7,
    'rect_th': 0.7,
    'degree_std':0.1,
    'rect_std':0.1,
    'azimuth_std':10,
    'elevation_std':10
}


class LTEstats(AttribDict):
    def __init__(self, header):
        super(LTEstats, self).__init__(header)

    def __add_attr__(self, list):
        self.attributes = list

    def __str__(self):
        priorized_keys = ['id','channel','starttime','endtime','interval','int_olap','step','step_olap','sample_rate','remove_response','fq_band','nro_time_bins','last_time_bin','nro_freq_bins','polar_degree','polar_analysis','f_params', 'opt_params','f_threshold','PE_tau','PE_order','time_bandwidth']

        return self._pretty_str(priorized_keys)


class LTEProcess(object):
    def __init__(self, h5file, ncomp, args):
        self.file_ = h5file
        self.read_nfs()
        self.ncomp = ncomp
        self.init_args = args
        self.n = -1
        self.processes = []
        self.queue = Queue()
    
    def read_nfs(self):
        with h5py.File(self.file_, "r+") as h5f:
            self.nfs = h5f['header'].attrs["nro_freq_bins"]

    def _wrapper(self, *args):
        data = {}
        if arg[0] != None:
            julia.Julia(compiled_modules=False)
            from julia import LTE as lte
            ret = lte.lte_run(
                args[0],
                args[1],
                self.init_args[0],
                self.init_args[1],
                self.init_args[2],
                self.init_args[3],
                self.init_args[4],
                self.init_args[5],
                self.init_args[6],
                self.init_args[7],
                self.init_args[8],
                self.init_args[9],
                self.init_args[10],
                self.init_args[11],
                self.init_args[12],
                self.init_args[13],
                self.init_args[14]
                )
            
            # fill ddict
            for n in range(self.ncomp):
                data[n+1] = {}
                data[n+1]["spec"] = ret[n+1]["spec"]
                for attr in ("erg", "dfq", "cfq", "pe"):
                        data[n+1][attr] = ret[n+1][attr] 

                if self.init_args[8]:
                    data["opt"] = {}
                    for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf"):
                        data["opt"][attr]  = ret[attr]

                    if self.init_args[9]:
                        data["opt"]["dsar"]  = ret["dsar"]        

                if self.init_args[10]:
                    data["polar"] = {}
                    for attr in ("degree", "rect", "azimuth", "elev"):
                        data["polar"][attr]  = ret[attr]
        else:
            # fill ddict with nan
            nwin = self.init_args[1]
            for n in range(args[0].shape[0]):
                data[n+1] = {}
                data[n+1]["spec"] = np.zeros((nwin,self.nfs))*np.nan
                for attr in ("erg", "dfq", "cfq", "pe"):
                        data[n+1][attr] = np.zeros((nwin,))*np.nan

                if self.init_args[8]:
                    data["opt"] = {}
                    for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf"):
                        data["opt"][attr]  = np.zeros((nwin,))*np.nan

                    if self.init_args[9]:
                        data["opt"]["dsar"]  = np.zeros((nwin,))*np.nan        

                if self.init_args[10]:
                    data["polar"] = {}
                    for attr in ("degree", "rect", "azimuth", "elev"):
                        data["polar"][attr]  = np.zeros((nwin,self.nfs))*np.nan

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
                
                # add data
                


                h5f['header'].attrs.modify('last_time_bin', nbin+nwin)
                h5f.flush()

        self.reset()


class LTE(object):
    def __init__(self, lte_file, v2=False):

        if not os.path.isfile(lte_file):
            return ValueError(' lte_file do not found')
        
        self.file_ = lte_file
        self.v2 = v2
        self.__set_stats__()
        self.__set_attrs__()


    def __str__(self):
        return self.stats.__str__()


    def is_attr(self, attr, only_scalars=False, only_vectors=False):
        if not isinstance(attr, list):
            list_attr = [attr]
        else:
            list_attr = attr

        if not all([attr in self.stats.attributes for attr in list_attr]):
            attrs = [at not in self.stats.attributes for at in list_attr]
            pos = list(filter(lambda x: attrs[x], range(len(list_attr))))
            not_availabel_attr = np.array(list_attr)[pos]
            print('warn: attributes %s not available' % not_availabel_attr)

            return_list = [attr for attr in list_attr if attr in self.stats.attributes]
            
            if not return_list:
                print('available attr: %s' % self.stats.attributes)

        else:
            return_list = list_attr
        
        if only_scalars:
            return_list = [attr for attr in return_list if attr in SCALAR_PARAMS]
        
        if only_vectors:
            return_list = [attr for attr in return_list if attr in VECTORAL_PARAMS]
            
        return return_list


    def is_matrix(self, attr):
        return attr in VECTORAL_PARAMS


    def is_chan(self, chan):
        if isinstance(chan, list):
            chan_list = [ch for ch in chan if ch in self.stats.channel]
        else:
            if isinstance(chan, type(None)):
                chan_list = list(self.stats.channel)
            else:
                if chan in self.stats.channel:
                    chan_list = [chan]
                else:
                    chan_list = []
        
        return chan_list


    def __set_stats__(self):
        with h5py.File(self.file_, "r") as f:
            hdr = f['header']

            if self.v2:
                lte_stats = dict(
                    id = hdr.attrs['id'],
                    channel = hdr.attrs['channel'],
                    starttime = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                    endtime = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                    nro_time_bins = hdr.attrs['nro_time_bins'],
                    last_time_bin = hdr.attrs['last_time_bin'],
                    step = hdr.attrs['step'],
                    step_olap = hdr.attrs['step_olap'],
                    fq_band = hdr.attrs['fq_band'],
                    pad = hdr.attrs['pad'],
                    nro_freq_bins = hdr.attrs['nro_freq_bins'],
                    sample_rate = hdr.attrs['sample_rate'],
                    remove_response = hdr.attrs['remove_response'],
                    polar = hdr.attrs['polar'],
                    opt_params = hdr.attrs['opt_params'],
                    opt_th = hdr.attrs['opt_th'],
                    opt_twin = hdr.attrs['opt_twin'],
                    PE_tau = hdr.attrs["PE_tau"],
                    PE_order = hdr.attrs["PE_order"],
                    time_bandwidth = hdr.attrs['time_bandwidth']
                    )
            else:
                lte_stats = dict(
                    id = hdr.attrs['id'],
                    channel = hdr.attrs['channel'],
                    starttime = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                    endtime = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                    interval = hdr.attrs['interval'],
                    int_olap = hdr.attrs['int_olap'],
                    nro_time_bins = hdr.attrs['nro_time_bins'],
                    last_time_bin = hdr.attrs['last_time_bin'],
                    step = hdr.attrs['step'],
                    step_olap = hdr.attrs['step_olap'],
                    fq_band = hdr.attrs['fq_band'],
                    nro_freq_bins = hdr.attrs['nro_freq_bins'],
                    sample_rate = hdr.attrs['sample_rate'],
                    remove_response = hdr.attrs['remove_response'],
                    polar_degree = hdr.attrs['polar_degree'],
                    polar_analysis = hdr.attrs['polar_analysis'],
                    f_params = hdr.attrs['f_params'],
                    opt_params = hdr.attrs['opt_params'],
                    f_threshold = hdr.attrs["f_threshold"],
                    PE_tau = hdr.attrs["PE_tau"],
                    PE_order = hdr.attrs["PE_order"],
                    time_bandwidth = hdr.attrs['time_bandwidth']
                    )
        
        self.stats = LTEstats(lte_stats)


    def __set_attrs__(self):
        attrs = SCALAR_PARAMS[0:4]
        attrs.append(VECTORAL_PARAMS[0])
        
        if self.stats.opt_params:
            attrs += SCALAR_PARAMS[4:12]

            if self.stats.remove_response:
                attrs.append(SCALAR_PARAMS_OPT[-1]) 
         
        if self.stats.polar:
            attrs += VECTORAL_PARAMS[1:]
                
        self.stats.__add_attr__(attrs)



    # def __save__(self, dict_in, tbin):
    #     with h5py.File(self.file_, "r+") as f:
    #         for chan in self.stats.channel:
    #             chan_dict = dict_in[chan]
    #             chan_gr = f.get(chan)

    #             for attr, val in chan_dict.items():
    #                 group = get_group(attr)
    #                 dset = chan_gr.get(group)[attr]

    #                 if isinstance(val, type(None)):
    #                     val = np.nan
                    
    #                 if self.is_matrix(attr):
    #                     dset[tbin, :] = val
                    
    #                 else:
    #                     dset[tbin,] = val
            
    #         if self.stats.polar_degree:
    #             for attr, val in dict_in['polar'].items():
    #                 dset = f.get('polar')[attr]

    #                 if isinstance(val, type(None)):
    #                     val = np.nan
                    
    #                 dset[tbin, :] = val
            
    #         f['header'].attrs.modify('last_time_bin', tbin+1)
    #         f.flush()


    def __compute__(self, station, headers, njobs):
        
        with tqdm(total=total_bins) as pbar:
            start = self.stats.starttime
            interval = dt.timedelta(minutes=headers["interval"])

            ltep = LTEProcess(
                self.file_,[
                int(self.stats.sample_rate),
                int(headers["nwin"]),
                int(headers["lwin"]),
                int(headers["nadv"]),
                0,
                self.stats.fq_band,
                self.stats.time_bandwidth,
                self.stats.pad,
                bool(self.stats.remove_response*self.stats.opt_params),
                self.stats.opt_params,
                self.stats.polar,
                self.stats.PE_order,
                self.stats.PE_tau,
                self.stats.opt_twin,
                self.stats.opt_th
                ]
            )

            while start + interval <= self.stats.endtime:
                end = start + interval
                
                mdata = station.get_mdata(start, end, channel=self.stats.channel, remove_response=self.stats.remove_response, sample_rate=self.stats.sample_rate, sort="ZNE")

                if mdata and self.stats.remove_response and self.stats.opt_params:
                    zchan = station.get_chan("Z")[0]
                    sdata = station.get_mdata(start, end, channel=zchan, remove_response=self.stats.remove_response, sample_rate=self.stats.sample_rate, rrkwargs={"output":'DISP'})
                else:
                    sdata = None
                
                ltep.run(mdata, sdata) # send to a cpu
                
                if len(ltep.processes) == njobs:
                    ltep.wait()
                    ltep.save()

                # advance 
                start += interval
                pbar.update()
            
            # save last precesse



    def get(self, attr, chan=None, starttime=None, endtime=None, **kwargs):
        
        attr_list = self.is_attr(attr)

        if not attr_list:
            print(' attribute unknown')
            return None
        
        chan_list = self.is_chan(chan)

        if not list(chan_list):
            print(' channel unknown')
            return None

        time, (n_0, n_f) = self.get_time(starttime=starttime, endtime=endtime, **kwargs)

        dout = {}

        if len(chan_list) > 1:
            multichanel = True
            for chan in chan_list:
                dout[chan] = {}
        
        else:
            multichanel = False
        
        rout = kwargs.get('remove_outlier')
        
        for attr in attr_list:
            group = get_group(attr)

            with h5py.File(self.file_, "r") as f:
                if group == 'polar':
                    data = f.get('polar')[attr][n_0:n_f,:]
                    freq = f.get('freq')[:]

                    # remove outlines
                    if rout:
                        for chan in chan_list:
                            if isinstance(rout, (float, int)):
                                spikes = self.get_outlier(chan, rout, n_0, n_f)
                            else:
                                spikes = None
                            
                        if isinstance(spikes, np.ndarray):
                            data[spikes] = np.nan    

                    dout[attr] = (data, freq)
                
                else:
                    for chan in chan_list:

                        if isinstance(rout, (float, int)):
                            spikes = self.get_outlier(chan, rout, n_0, n_f)
                        else:
                            spikes = None
                                
                        if self.is_matrix(attr):
                            data = f.get(chan)[group][attr][n_0:n_f,:]
                            freq = f.get('freq')[:]

                            if isinstance(spikes, np.ndarray):
                                data[spikes] = np.nan

                            to_return = (data, freq)
                        
                        else:
                            data = f.get(chan)[group][attr][n_0:n_f]
                            if isinstance(spikes, np.ndarray):
                                data[spikes] = np.nan
                            
                            to_return = data

                        if multichanel:
                            dout[chan][attr] = to_return
                        
                        else:
                            dout[attr] = to_return
    
        return time, dout    


    def get_outlier(self, chan, threshold, n0, nf):
        with h5py.File(self.file_, "r") as f:
            erg = 10*np.log10(f.get(chan)['spectral']['energy'][n0:nf])
            outliers = np.where(erg > threshold)[0]
        
        return outliers


    def get_time(self, starttime=None, endtime=None, **kwargs):
        
        if not starttime or starttime < self.stats.starttime:
            starttime = self.stats.starttime
            n_0 = 0
        else:
            delta = (starttime - self.stats.starttime).total_seconds()/60
            n_0 = ((delta/self.stats.interval) - self.stats.int_olap) / (1 - self.stats.int_olap)
            n_0 = int(np.ceil(n_0))

        if not endtime or endtime > self.stats.endtime:
            endtime = self.stats.endtime
            n_f = self.stats.nro_time_bins

        else:
            delta = (endtime - self.stats.starttime).total_seconds()/60
            n_f = ((delta/self.stats.interval) - self.stats.int_olap) / (1 - self.stats.int_olap)
            n_f = int(np.ceil(n_f))
        
        # return datetime or int in hours

        if kwargs.get('datetime', False):
            time1 = [starttime]
            time2 = [starttime + dt.timedelta(minutes=k*self.stats.interval-(k-1)*self.stats.interval*self.stats.int_olap) for k in range(n_0,n_f)]
            true_time = time1 + time2[1:]
        
        else:
            total_duration = (endtime - starttime).total_seconds()/3600
            true_time = np.linspace(0, total_duration, n_f-n_0)

        return true_time, (n_0, n_f)


    def get_stats(self, attr, chan=None, starttime=None, endtime=None):
        """
        Return (min,max,mean,mode)
        """

        def attr_stats(x):
            v_min = x.min()
            v_max = x.max()
            v_mean = x.mean()

            x_range = np.linspace(v_min, v_max, 500)
            gkde = stats.gaussian_kde(x)
            kde = gkde(x_range)
            v_mode = x_range[np.argmax(kde)]

            return [v_min, v_max, v_mean, v_mode]

        attr_list = self.is_attr(attr, only_scalars=True)
        chan_list = self.is_chan(chan)

        if len(chan_list) > 1:
            multichanel = True
        else:
            multichanel = False

        if list(attr_list) and list(chan_list):
            ans = self.get(attr_list, chan=chan_list, starttime=starttime, endtime=endtime)

            if isinstance(ans, type(None)):
                return None
            
            din = ans[1] # time: 0

            dout = {}
            if multichanel:
                for chan, ch_dict in din.items():
                    dout[chan] = {}

                    for attr, data in ch_dict.items():
                        x_filt = data[np.isfinite(data)]

                        if attr == 'energy':
                            x_filt = 10*np.log10(x_filt)

                        dout[chan][attr] = attr_stats(x_filt)
            
            else:
                for attr, data in din.items():
                    x_filt = data[np.isfinite(data)]
                    
                    if attr == 'energy':
                        x_filt = 10*np.log10(x_filt)
                    
                    dout[attr] = attr_stats(x_filt)
            
            return dout


    def pdf(self, attr, chan=None, starttime=None, endtime=None, bandwidth=0.01, **kwargs):
        
        if not chan:
            chan = self.stats.channel[0]
        
        if chan not in self.stats.channel:
            raise ValueError(' chan not found')
        
        if attr not in self.stats.attributes:
            raise ValueError(' attr not availabel')
        

        _, dout = self.get(attr, chan=chan, starttime=starttime, endtime=endtime)
        y_size = kwargs.get('y_space', 1000)

        if self.is_matrix(attr):
            data = dout[attr][0]
            freq = dout[attr][1]

            masked_data = np.ma.masked_invalid(data)
            data = np.ma.compress_rowcols(masked_data, axis=0)

            if attr == 'specgram':
                data = 10*np.log10(data)
            
            dmin = kwargs.get('y_min', data.min())
            dmax = kwargs.get('y_max', data.max())

            y_space = np.linspace(dmin, dmax, y_size).reshape(y_size, 1)
            pdf = get_PDF(data, y_space, bandwidth, **kwargs)
        
            return freq, y_space.reshape(y_size,), pdf
        
        else:
            data = dout[attr]
            dmin = kwargs.get('y_min', data.min())
            dmax = kwargs.get('y_max', data.max())

            if attr == 'energy':
                data = 10*np.log10(data)

            y_space = np.linspace(dmin, dmax, y_size).reshape(y_size, 1)
            pdf = get_PDF(data.reshape(-1,1), y_space, bandwidth, **kwargs)

            return y_space.reshape(y_size,), pdf


    def __get_peaks__(self, chan, starttime=None, endtime=None, fq_range=(), peak_distance=5, peak_thresholds={}, **kwargs):
        """

        Return peaks as describes in Melchor et al 2021

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
        
        if chan not in self.stats.channel:
            raise ValueError(' chan not found')

        if not self.stats.polar_degree:
            raise ValueError(' No polarization data in LTE')

        self.peak_thresholds = default_peak_thresholds
        
        # define parameters
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

        if not fq_range:
            fq_range = self.stats.fq_band

        attr_list = ['specgram', 'degree']
        if self.stats.polar_analysis:
            attr_list += ['elevation', 'rect', 'azimuth']
        
        time, dout = self.get(['energy']+attr_list, chan=chan, starttime=starttime, endtime=endtime, **kwargs)
        erg = dout['energy']
        sxx, freq = dout['specgram']
        pd, _ = dout['degree']

        # get fq_npts
        freq = np.array(freq)
        fq0 = np.argmin(np.abs(freq-fq_range[0]))
        fq1 = np.argmin(np.abs(freq-fq_range[1]))

        if self.stats.polar_analysis:
            rect, _ = dout['rect']
            thetaH, _ = dout['azimuth']
            thetaV, _ = dout['elevation']

        sxx = 10*np.log10(sxx)
        total_time = len(time)
        peaks = {}
        nro_Wpeaks = 0
        nro_Lpeaks = 0
        with tqdm(total=total_time) as pbar:
            for t in range(total_time):
                sxx_t = sxx[t, fq0:fq1+1]
                pd_t = pd[t, fq0:fq1+1]
                if np.isfinite(sxx_t).all():
                    sxx_t_norm = MinMaxScaler().fit_transform(sxx_t.reshape(-1,1))
                    peaks_pos,_ = signal.find_peaks(sxx_t_norm.reshape(-1,), height=self.peak_thresholds['specgram_th'], distance=peak_distance)

                    peaks[t] = dict(fq=[], specgram=[], degree=[])
                    if self.stats.polar_analysis:
                        peaks[t]['rect'] = []
                        peaks[t]['azimuth'] = []
                        peaks[t]['elevation'] = []

                    for p in peaks_pos:
                        fp = freq[fq0:fq1][p]

                        # define freq range                        
                        fp_0 = fp - self.peak_thresholds['fq_delta']
                        fp_1 = fp + self.peak_thresholds['fq_delta']
                        fq_0_pos = np.argmin(np.abs(freq-fp_0))
                        fq_1_pos = np.argmin(np.abs(freq-fp_1))


                        # compute PSD mean, PD mean and PD std
                        sp_peak = sxx_t[fq_0_pos:fq_1_pos].mean()
                        pd_peak = pd_t[fq_0_pos:fq_1_pos]
                        pd_peak_avg, pd_peak_std = pd_peak.mean(), pd_peak.std()

                        if pd_peak_avg > self.peak_thresholds['degree_th'] and pd_peak_std < self.peak_thresholds['degree_std']:
                            peaks[t]['fq'] += [fp]
                            nro_Wpeaks += 1

                            # only when peak is well polairzed we add peak
                            peaks[t]['specgram'] += [sp_peak]
                            peaks[t]['degree'] += [pd_peak_avg]

                            if self.stats.polar_analysis:
                                rect_peak = rect[t, fq_0_pos:fq_1_pos]
                                tH_peak = thetaH[t, fq_0_pos:fq_1_pos]
                                tV_peak = thetaV[t, fq_0_pos:fq_1_pos]
                                
                                rect_peak_val = None
                                tH_peak_val = None
                                tV_peak_val = None

                                # remove low PD values
                                rect_peak = rect_peak[np.where(pd_peak > self.peak_thresholds['degree_th'])]
                                tH_peak = tH_peak[np.where(pd_peak > self.peak_thresholds['degree_th'])]
                                tV_peak = tV_peak[np.where(pd_peak > self.peak_thresholds['degree_th'])]
                                
                                if rect_peak.any():
                                    rect_peak_avg, rect_peak_std = rect_peak.mean(), rect_peak.std()
                                    
                                    if rect_peak_std < self.peak_thresholds['rect_std']:
                                        rect_peak_val = rect_peak_avg

                                        if rect_peak_avg > self.peak_thresholds['rect_th']:
                                            
                                            tH_peak = tH_peak[np.where(rect_peak > self.peak_thresholds['rect_th'])]
                                            tV_peak = tV_peak[np.where(rect_peak > self.peak_thresholds['rect_th'])]
                                            nro_Lpeaks += 1
                                            
                                            if tH_peak.any():
                                                tH_peak_rad = tH_peak*np.pi/180
                                                tH_peak_avg = stats.circmean(tH_peak_rad, high=np.pi)*180/np.pi
                                                tH_peak_std = stats.circstd(tH_peak_rad, high=np.pi)*180/np.pi
                                                
                                                if tH_peak_std <= self.peak_thresholds['azimuth_std']:
                                                    tH_peak_val = tH_peak_avg
                                            
                                            if tV_peak.any():
                                                tV_peak_avg, tV_peak_std = tV_peak.mean(), tV_peak.std()

                                                if tV_peak_std <= self.peak_thresholds['elevation_std']:
                                                    tV_peak_val = tV_peak_avg
                                
                                peaks[t]['rect'] += [rect_peak_val]
                                peaks[t]['azimuth'] += [tH_peak_val]
                                peaks[t]['elevation'] += [tV_peak_val]

                pbar.update()
        
        return (nro_Wpeaks, nro_Lpeaks), peaks


    def get_Peaks(self, chan, starttime=None, endtime=None, fq_range=(), peak_thresholds={}, **kwargs):
        return Peaks(self, chan, starttime=starttime, endtime=endtime, fq_range=fq_range, peak_thresholds=peak_thresholds, **kwargs)
    

    def get_WLpeaks(self, time_length, file_out=None, chan_list=None, starttime=None, endtime=None, fq_range=(), peak_thresholds={}):

        # check starttime and endtime
        if starttime:
            start_time = starttime
            assert start_time >= self.stats.starttime, "start_time error: check lte file"
        else:
            start_time = self.stats.starttime
        
        if endtime:
            end_time = endtime
            assert end_time <= self.stats.endtime, "end_time error: check lte file"
        else:
            end_time = self.stats.endtime

        time_window = dt.timedelta(minutes=time_length)
        nro_intervals = int(((end_time-start_time).total_seconds()/60) / time_length)

        nro_Wpeaks = np.empty(nro_intervals)
        nro_Lpeaks = np.empty(nro_intervals)
        time = []

        if not chan_list:
            chan_list = self.stats.channel

        for n in range(0, nro_intervals):
            start = start_time + dt.timedelta(minutes=n*time_length)
            end = start + time_window

            nro_wpeak_n = 0
            nro_lpeak_n = 0
            for chan in chan_list:
                pks, _ = self.__get_peaks__(chan, starttime=start, endtime=end, fq_range=fq_range, peak_thresholds=peak_thresholds, verbose=False)
                nro_wpeak_n += pks[0]
                nro_lpeak_n += pks[1]
            
            time += [start]
            nro_Wpeaks[n] = nro_wpeak_n
            nro_Lpeaks[n] = nro_lpeak_n

        data = {
            "time":time,
            "NW": nro_Wpeaks,
            "NL": nro_Lpeaks
            }

        if file_out:
            df = pd.DataFrame(data)
            df.to_csv(f'{file_out}.csv')
        
        return data


    def plot(self, chan, list_attr, starttime=None, endtime=None, interval=None, init_gui=False, lde=None, **kwargs):
        # check times
        if not starttime:
            starttime = self.stats.starttime 
        
        if starttime < self.stats.starttime:
            raise ValueError(' start < lte.stats.starttime')
        
        if not endtime:
            endtime = self.stats.endtime
        
        if endtime > self.stats.endtime:
            raise ValueError(' end > lte.stats.endtime')

        # check interval
        max_int = (endtime - starttime).days

        if not interval or max_int < interval:
            interval = max_int

        # check attributes
        list_attr = self.is_attr(list_attr)
        
        if not list_attr:
            raise ValueError('not attr available')

        # check channel 
        chan_list = self.is_chan(chan)

        if not chan_list:
            raise ValueError('not chan available')
            
        if init_gui:
            from seisvo import LDE, default_LDE_dir
            from seisvo.plotting.gui.glte import plot_gui

            if not lde:
                lde_file = os.path.join(default_LDE_dir, self.lte.stats.id + '.lde')

                if not os.path.isdir(default_LDE_dir):
                    os.makedirs(default_LDE_dir)

                lde = LDE(lde_file)
        
            else:
                if isinstance(lde, LDE):
                    lde = lde
            
                elif isinstance(lde, str):
                    if lde.split('.')[-1] != 'lde':
                        lde += '.lde'
                    lde = LDE(lde)
            
                else:
                    raise ValueError('lde should be LDE, string or None')

            plot_gui(chan_list, self, starttime, endtime, interval, list_attr, lde=lde, **kwargs)

        else:
            from seisvo.plotting.base.lte import plotLTE, plt
            
            fig = plt.figure(figsize=kwargs.get('figsize',(12,9)))
            plte = plotLTE(chan_list, fig, self, starttime, endtime, interval, list_attr, **kwargs)
            
            return plte
    

    def plot_peak_tevo(self, chan, list_attr, fq_range_dict, marker_dict={}, color_dict={}, starttime=None, endtime=None, fig=None, **kwargs):
        
        # example of dictionary inputs:
        #  > fq_range_dict = dict(fq1:{'width': 
        #                            ...

        from seisvo.plotting.base.lte import plotDPeakTEVO
            
        attr_list = self.is_attr(list_attr, only_vectors=True)

        if attr_list and self.stats.polar_degree:
            chan_list = self.is_chan(chan)

            if not starttime:
                starttime = self.stats.starttime
            
            if not endtime:
                endtime = self.stats.endtime

            tevo = plotDPeakTEVO(self, starttime, endtime, chan_list, attr_list, fq_range_dict, fig, **kwargs)

            return tevo.fig


    @staticmethod
    def new(station, lte_file, headers, njobs):
        """
        Create new LTE (hdf5) file
        """

        with h5py.File(lte_file, "w-") as f:
            # header dataset
            hdr = f.create_dataset('header',(1,))
            hdr.attrs['id'] = headers['id']
            hdr.attrs['channel'] = headers['channel']
            hdr.attrs['starttime'] = headers['starttime']
            hdr.attrs['endtime'] = headers['endtime']
            hdr.attrs['step'] = headers['step']
            hdr.attrs['step_olap'] = headers['step_olap']
            
            hdr.attrs['sample_rate'] = headers['sample_rate']
            hdr.attrs['pad'] = headers['pad']
            hdr.attrs['fq_band'] = headers['fq_band']
            hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
            hdr.attrs['last_time_bin'] = -1
            
            # optional param
            hdr.attrs['polar'] = headers['polar']
            hdr.attrs['opt_params'] = headers['opt_params']
            hdr.attrs['remove_response'] = headers['remove_response']
            hdr.attrs['opt_twin'] = headers['opt_twin']
            hdr.attrs['opt_th'] = headers['opt_th']
            
            hdr.attrs['PE_tau'] = headers['pe_tau']
            hdr.attrs['PE_order'] = headers['pe_order']
            hdr.attrs['time_bandwidth'] = headers['time_bandwidth']

            freqbins = headers['nro_freq_bins']
            timebins = headers['nro_time_bins']

            # print info
            print('')
            print(' LTE file INFO')
            print(' -------------')
            print(" hdf5_memory info: %s " % lte_file)
            print(' --- dataset size: ', (timebins, freqbins))
            print('')
            print(' LTE stats:')
            for info_key in ['id', 'channel', 'starttime', 'endtime', 'step' ,'step_olap' ,'sample_rate' ,'remove_response' ,'fq_band' ,'polar' ,'opt_params']:
                print(f' {info_key}:  {headers[info_key]}')

            # create datasets structure
            # f.create_dataset('freq', (freqbins,), dtype=np.float32)

            for chan in headers['channel']:
                chgr = f.create_group(chan)
                amp = chgr.create_group("amplitude")
                spec = chgr.create_group("spectral")
                amp.create_dataset('pentropy', (timebins,), chunks=True, dtype=np.float32)
                spec.create_dataset('specgram', (timebins, freqbins), chunks=True, dtype=np.float32)
                spec.create_dataset('fq_dominant', (timebins,), chunks=True, dtype=np.float32)
                spec.create_dataset('fq_centroid', (timebins,), chunks=True, dtype=np.float32)
                spec.create_dataset('energy', (timebins,), chunks=True, dtype=np.float32)


                if headers['opt_params']:
                    amp.create_dataset('vlf', (timebins,), chunks=True, dtype=np.float32)
                    amp.create_dataset('lf', (timebins,), chunks=True, dtype=np.float32)
                    amp.create_dataset('mf', (timebins,), chunks=True, dtype=np.float32)
                    amp.create_dataset('hf', (timebins,), chunks=True, dtype=np.float32)
                    amp.create_dataset('rsam', (timebins,), chunks=True, dtype=np.float32)
                    amp.create_dataset('vlar', (timebins,), chunks=True, dtype=np.float32)
                    amp.create_dataset('lrar', (timebins,), chunks=True, dtype=np.float32)
                    amp.create_dataset('rmar', (timebins,), chunks=True, dtype=np.float32)

                    if headers['remove_response'] and chan[-1] == 'Z':
                        amp.create_dataset('dsar', (timebins,), chunks=True, dtype=np.float32)

            if headers['polar']:
                polar = f.create_group("polar")
                polar.create_dataset('degree', (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('rect', (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('elevation', (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('azimuth', (timebins, freqbins), chunks=True, dtype=np.float32)
                
            f.flush()

        lte = LTE(lte_file, v2=True)
        lte.__compute__(station, headers, njobs)

        return lte


class Peaks(LTE):
    def __init__(self, lte, chan, starttime=None, endtime=None, fq_range=(), peak_thresholds={}, **kwargs):

        if isinstance(lte, str):
            super().__init__(lte)
        
        elif isinstance(lte, LTE):
            super().__init__(lte.lte_file)
        else:
            raise ValueError ('lte should be string with lte file path or lte object')

        ans = super().__get_peaks__(chan, starttime=starttime, endtime=endtime, fq_range=fq_range, peak_thresholds=peak_thresholds, **kwargs)

        self.total_peaks_, self.dominant_peaks_ = ans
        self.chan = chan
        # self.peak_thresholds_ = self.peak_thresholds

        if not fq_range:
            self.fq_range = self.stats.fq_band
        else:
            self.fq_range = fq_range
        
        if not starttime:
            self.starttime = self.stats.starttime
        else:
            self.starttime = starttime

        if not endtime:
            self.endtime = self.stats.endtime
        else:
            self.endtime = endtime

        self.peaks_ = {}
        self.nro_ = 0
        self.df_ = None
    
    
    def __str__(self):
        txt_to_return =  f'\n   >>LTE file    ::: {self.file_}'
        txt_to_return += f'\n   >channel       : {self.chan}'
        txt_to_return += f'\n   >starttime     : {self.stats.starttime.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >endtime       : {self.stats.endtime.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >fq_range      : {self.fq_range}'
        txt_to_return += f'\n   >polar_analysis: {self.stats.polar_analysis}'
        txt_to_return += f'\n   >Nro dom. peaks: {self.total_peaks_[0]}'
        txt_to_return += f'\n   >Nro chr. peaks: {self.nro_}'
        txt_to_return +=  f'\n'
        return txt_to_return


    def fit(self, threshold=0.0, peak_width='auto', **kwargs):
        """
        Computes the dominant peaks following Melchor et al 2021?
        peak_width can be 'auto' or float. If 'auto', for each peak, the width is computing by signal.peak_widths.
        """

        if not isinstance(peak_width, float) and not peak_width=='auto':
            raise ValueError (" peak_width must be float or 'auto'")
        
        if isinstance(threshold, float) and not 0 <= threshold <= 1:
            raise ValueError (" threshold must bounded by (0,1)")

        self.threshold = threshold
        self.peak_width = peak_width
        
        # fit kernel distribution
        self.__fit_fq_kde__(**kwargs)

        # search for dominant frequencies
        fq_pdf_norm = self.fq_pdf_/self.fq_pdf_.max()
        peaks, _ = signal.find_peaks(fq_pdf_norm, height=threshold, distance=5)
        self.nro_ = len(peaks)

        # search for width of dominant frequencies
        if isinstance(peak_width, float):
            half_width = peak_width/2
        elif peak_width == 'auto':
            delta = np.diff(self.fq_space_.reshape(-1,))[0]
            peak_width = signal.peak_widths(self.fq_pdf_, peaks)[0]
            half_width = None
        else:
            raise ValueError('peak_width must be float or "auto"')

        for ip, p in enumerate(peaks):

            if half_width:
                pwidth = peak_width/2
            
            else:
                p_half_width = (delta*peak_width[ip])/2
                pwidth = max(p_half_width, 0.05)
                pwidth = min(p_half_width, 0.50)

            self.peaks_[ip+1] = {
                'fq':self.fq_space_[p],
                'width':pwidth,
                'fq_prob':self.fq_pdf_[p],
                'cp':None,
                'cl':None,
                'sp':None,
                'rect':None,
                'thH':None,
                'thV':None
                }
        
        if self.nro_ > 0:
            self.__fit_polar__(**kwargs)
        
        self.__dataframe__()


    def __fit_fq_kde__(self, **kwargs):
        # compute KDE for peak distribution weighted by spectral energy
        all_dominant_peaks = [nd['fq'] for _, nd in self.dominant_peaks_.items()]
        all_dominant_peaks = np.array(list(chain(*all_dominant_peaks)))

        bandwidth = kwargs.get('bandwidth', {})
        self.bandwidth_fc = bandwidth.get('freq', default_bandwidth['freq'])

        v_min = kwargs.get('v_min', None)
        if not v_min:
            v_min = all_dominant_peaks.min()
            
        v_max = kwargs.get('v_max', None)
        if not v_max:
            v_max = all_dominant_peaks.max()
        
        v_space = kwargs.get('v_space', 1000) # fqa.shape[0]*5
        self.fq_space_ = np.linspace(v_min, v_max, v_space).reshape(v_space, 1)

        if kwargs.get('weighted', False):
            sp_l = [item['sp'] for item in self.all_peaks if item]
            spa = np.array(list(chain(*sp_l)))
            self.weights = spa/spa.max()

        else:
            self.weights = None

        kde = get_KDE(all_dominant_peaks.reshape(-1, 1), self.bandwidth_fc, weight=self.weights)

        self.fq_pdf_ = np.exp(kde.score_samples(self.fq_space_))


    def __fit_polar__(self, **kwargs):
        self.dist_throld = kwargs.get('dist_throld', 0.5)
        self.n_sample_min = kwargs.get('n_sample_min', 5)
        
        # default bandwidth
        usr_bandwidth = kwargs.get('bandwidth', {})
        self.bandwidth = default_bandwidth
        for key in self.bandwidth.keys():
            if isinstance(usr_bandwidth.get(key, None), float) or isinstance(usr_bandwidth.get(key, None), np.ndarray):
                self.bandwidth[key] = usr_bandwidth.get(key)
        
        print('\n  ---- Bandwidth init model ---')
        for key, item in self.bandwidth.items():
            if isinstance(item, np.ndarray):
                print(f'     {key:>10}    {item.min():.2f}--{item.max():.2f}')
            else:
                print(f'     {key:>10}    {item:.2f}')

        print('  ---- --------- ---- ---- ---\n')

        print('\n     #peak   attr   bandwidth')

        with tqdm(total=self.nro_) as pbar:
            for f in range(1, self.nro_+1):
                dfq = self.peaks_[f]['fq'] # dominant frequency
                sp = []
                pd = []
                
                if self.stats.polar_analysis:
                    rect = [] 
                    th_H = [] 
                    th_V = []

                for _, tdict in self.dominant_peaks_.items():
                    for n, fq in enumerate(tdict['fq']):
                        if dfq-self.peaks_[f]['width'] <= fq <= dfq+self.peaks_[f]['width']:
                            sp += [tdict['specgram'][n]]
                            pd += [tdict['degree'][n]]

                            if self.stats.polar_analysis:
                                rect += [tdict['rect'][n]]
                                th_H += [tdict['azimuth'][n]]
                                th_V += [tdict['elevation'][n]]
                
                # list sp, pd, etc. may contain None
                sp = np.array([sp_k for sp_k in sp if sp_k])
                pd = np.array([pd_k for pd_k in pd if pd_k])

                if self.stats.polar_analysis:
                    rect = np.array([rect_k for rect_k in rect if rect_k])
                    th_H = np.array([th_H_k for th_H_k in th_H if th_H_k])
                    th_V = np.array([th_V_k for th_V_k in th_V if th_V_k])
                
                # if the number of samples is low remove frequency.
                if sp.shape[0] > self.n_sample_min:
                
                    sp_kde = get_KDE(sp.reshape(-1,1), self.bandwidth['specgram'])

                    print(f'     {f:^5}   spec   {sp_kde.bandwidth:.2f}')

                    sp_x = np.linspace(sp.min(), sp.max(), 500).reshape(-1,1)

                    sp_dist = np.exp(sp_kde.score_samples(sp_x))
                    sp_dist /= sp_dist.max()
                    args = np.where(sp_dist > self.dist_throld)[0]
                    sp_range = (sp_x[args[0]], sp_x[args[-1]])
                    self.peaks_[f]['sp'] = {
                                'val':sp_x[np.argmax(sp_dist)],
                                'range':sp_range,
                                'n':len(sp),
                                'kde':sp_kde
                            }
                    # second condition: number of dominant peaks with high polarization degree.
                    if rect.shape[0] > self.n_sample_min:
                        rect_kde = get_KDE(rect.reshape(-1,1), self.bandwidth['rect'])

                        print(f'     {f:^5}   rect   {rect_kde.bandwidth:.2f}')

                        rect_x = np.linspace(0, 1, 500).reshape(-1,1)
                        rect_dist = np.exp(rect_kde.score_samples(rect_x))
                        rect_dist /= rect_dist.max()
                        args = np.where(rect_dist > self.dist_throld)[0]
                        rect_range = (rect_x[args[0]], rect_x[args[-1]])

                        # to save
                        self.peaks_[f]['rect'] = {
                            'val':rect_x[np.argmax(rect_dist)],
                            'range':rect_range,
                            'n':len(rect),
                            'kde':rect_kde
                        }
                        self.peaks_[f]['cp'] = len(rect)/len(sp)

                        # third condition: number of dominant peaks with high rectilinearity.
                        if th_H.shape[0] > self.n_sample_min and th_V.shape[0] > self.n_sample_min:
                            
                            # compute Cl
                            self.peaks_[f]['cl'] = (len(th_H) + len(th_V)) / (2*len(rect))
                        
                            thH_kde = get_KDE(th_H.reshape(-1,1), self.bandwidth['azimuth'])
                            thH_x = np.linspace(0, 180, 500).reshape(-1,1)

                            print(f'     {f:^5}   thH    {thH_kde.bandwidth:.2f}')

                            thH_dist = np.exp(thH_kde.score_samples(thH_x))
                            thH_dist /= thH_dist.max()
                            args = np.where(thH_dist > self.dist_throld)[0]

                            # save azimuth
                            self.peaks_[f]['thH'] = {
                                'val':thH_x[np.argmax(thH_dist)],
                                'range':(thH_x[args[0]], thH_x[args[-1]]),
                                'n':len(th_H),
                                'kde':thH_kde
                            }

                            thV_kde = get_KDE(th_V.reshape(-1,1), self.bandwidth['elevation'])
                            thV_x = np.linspace(0, 90, 500).reshape(-1,1)

                            print(f'     {f:^5}   thV    {thH_kde.bandwidth:.2f}')

                            thV_dist = np.exp(thV_kde.score_samples(thV_x))
                            thV_dist /= thV_dist.max()
                            args = np.where(thV_dist > self.dist_throld)[0]

                            # save elevation
                            self.peaks_[f]['thV'] = {
                                'val':thV_x[np.argmax(thV_dist)],
                                'range':(thV_x[args[0]], thV_x[args[-1]]),
                                'n':len(th_V),
                                'kde':thV_kde
                            }
                else:
                    del self.peaks_[f]

                    pbar.update()


    def __dataframe__(self):
        if self.nro_ == 0:
            raise ValueError ('no dominant frequencies to save!')
        
        index = list(self.peaks_.keys())
        
        data = {
            'fq':[info['fq'][0] for _, info in self.peaks_.items()],
            'fq_prob':[info['fq_prob'] for _, info in self.peaks_.items()],
            'width':[info['width'] for _, info in self.peaks_.items()]
        }

        cl = []
        cp = []
        for  _, info in self.peaks_.items():
            if isinstance(info['cp'], type(None)):
                cp.append(np.nan)
            else:
                cp.append(info['cp'])

            if isinstance(info['cl'], type(None)):
                cl.append(np.nan)
            else:
                cl.append(info['cl'])
        data['cp'] = cp
        data['cl'] = cl

        for key_str in ('sp', 'rect', 'thH', 'thV'):
            k_data = []
            k_range = []
            k_n = []
            for  _, info in self.peaks_.items():
                if isinstance(info[key_str], type(None)):
                    k_data.append(np.nan)
                    k_range.append(np.nan)
                    k_n.append(np.nan)
                else:
                    k_data.append(info[key_str]['val'][0])
                    r = (info[key_str]['range'][0][0], info[key_str]['range'][1][0])
                    k_range.append(r)
                    k_n.append(info[key_str]['n'])
            
            if key_str == 'sp':
                data['S'] = k_data
                data['S_range'] = k_range
                data['N_T'] = k_n
            
            if key_str == 'rect':
                data['R'] = k_data
                data['R_range'] = k_range
                data['N_R'] = k_n
            
            if key_str == 'thH':
                data['H'] = k_data
                data['H_range'] = k_range
                data['N_H'] = k_n
            
            if key_str == 'thV':
                data['V'] = k_data
                data['V_range'] = k_range
                data['N_V'] = k_n

        self.df_ = pd.DataFrame(data, index=index)
    
    
    def to_json(self, fout=None, out_path='./'):
        if isinstance(self.df_, pd.DataFrame):
            if not fout:
                lte_file = os.path.basename(self.file_).split('.')
                lte_file.insert(-1, self.chan)
                
                if lte_file[-1] == 'lte':
                    lte_file[-1] = 'json'
                else:
                    lte_file += ['json']
                
                out = '.'.join(lte_file)
            
            json_file = os.path.join(out_path, out)

            if os.path.isfile(json_file):
                os.remove(json_file)

            self.df_.to_json(json_file)


    def get_dominant_peaks(self, fq_range):     
        dout = {}
        for t, tdict in self.dominant_peaks_.items():
            if tdict['fq']:
                dout[t] = {'fq':[], 'specgram':[]} 
                
                if self.stats.polar_analysis:
                    dout[t]['rect'] = []
                    dout[t]['azimuth'] = []
                    dout[t]['elevation'] = []
                
                for n, fq in enumerate(tdict['fq']):
                    if fq_range[0] <= fq <= fq_range[1]:
                        dout[t]['fq'].append(fq)
                        sp = tdict['specgram'][n]
                        dout[t]['specgram'].append(sp)

                        if self.stats.polar_analysis:
                            rect = tdict['rect'][n]
                            tH = tdict['azimuth'][n]
                            tV = tdict['elevation'][n]

                            if rect:
                                dout[t]['rect'].append(rect)
                            
                            if tH:
                                dout[t]['azimuth'].append(tH)
                            
                            if tV:
                                dout[t]['elevation'].append(tV)

        return dout


    def get_fq_area(self, fq, fqtol=0.5):

        aW = np.zeros((self.stats.nro_time_bins,))
        aL = np.zeros((self.stats.nro_time_bins,))
        
        for i in range(self.stats.nro_time_bins):
            if i in self.dominant_peaks_:
                d = self.dominant_peaks_[i]
                for f, pd, r in zip(d['fq'], d['degree'], d['rect']):
                    if (fq>f-fqtol) and (fq<f+fqtol):
                        if pd and pd > self.peak_thresholds['degree_th']:
                            aW[i] = 1
                            if r and r > self.peak_thresholds['rect_th']:
                                aL[i] = 1
        return aW, aL

    # PLOTTING    

    def plot_spec_pdf(self, show=True, **kwargs):
        if self.nro_ == 0:
            raise ValueError ('no dominant frequencies to plot!')

        from seisvo.plotting.base.lte import plotPeaksSpecPDF

        fig = plotPeaksSpecPDF(self, show=show, **kwargs)

        return fig


    def plot_peak_pdf(self, n, plot=True, **kwargs):
        if self.nro_ == 0:
            raise ValueError ('no dominant frequencies to plot!')

        from seisvo.plotting.base.lte import plotPeaksPDF
        
        fig = plotPeaksPDF(self, n, plot=plot, **kwargs)
        
        return fig