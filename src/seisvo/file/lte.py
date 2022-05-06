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

from itertools import chain
from tqdm import tqdm
from scipy import signal
from scipy.stats import gaussian_kde
from sklearn.preprocessing import MinMaxScaler, StandardScaler

from obspy.core.util.attribdict import AttribDict
from antropy import perm_entropy

from seisvo.signal.polarization import PolarAnalysis
# import seisvo.utils.plotting as sup
# from seisvo.signal.pdf import get_PDF, get_KDE

SCALAR_PARAMS = ['fq_dominant', 'fq_centroid', 'energy', 'pentropy', 'rsam', 'dsar']
SCALAR_PARAMS_OPT = ['mf', 'hf', 'vlf', 'lf', 'vlar', 'lrar', 'rmar']
SCALAR_PARAMS_F = ['rsam_f', 'mf_f', 'hf_f', 'vlf_f', 'lf_f', 'vlar_f', 'lrar_f', 'rmar_f', 'dsar_f']

VECTORAL_PARAMS = ['specgram', 'degree']
VECTORAL_PARAMS_OPT = ['elevation', 'rectlinearity', 'azimuth']

AMP_PARAMS = ['pentropy', 'rsam', 'dsar', 'mf', 'hf', 'lf', 'vlf', 'vlar', 'lrar', 'rmar', 'rsam_f', 'dsar_f', 'mf_f', 'hf_f', 'lf_f', 'vlf_f', 'vlar_f', 'lrar_f', 'rmar_f']
SPEC_PARAMS = ['specgram', 'fq_dominant', 'fq_centroid', 'energy']
POLAR_PARAMS = ['degree', 'elevation', 'rectlinearity', 'azimuth']

LTE_GROUPS = ['amplitude', 'spectral', 'polar']


def get_group(attr):
    if attr in AMP_PARAMS:
        if attr in SCALAR_PARAMS_F:
            return 'amplitude_f'
        else:
            return 'amplitude'

    if attr in SPEC_PARAMS:
        return 'spectral'

    if attr in POLAR_PARAMS:
        return 'polar'


default_bandwidth = {
            'frec': 0.01,
            'sp': 0.05,
            'rect': 0.25,
            'angle': 0.5
        }


default_peak_thresholds = {
    'fq_delta': 0.05,
    'spec_th': 0.7,
    'pd_th': 0.8,
    'rect_th': 0.7,
    'pd_std':0.1,
    'r_std':0.1,
    'azm_std':10,
    'dip_std':10
}


class LTEstats(AttribDict):
    def __init__(self, header):
        super(LTEstats, self).__init__(header)

    def __add_attr__(self, list):
        self.attributes = list

    def __str__(self):
        priorized_keys = ['id','channel','starttime','endtime','interval','int_olap','step','step_olap','sample_rate','remove_response','fq_band','nro_time_bins','nro_freq_bins','polar_degree','polar_analysis','f_params', 'opt_params','f_threshold','PE_tau','PE_order','time_bandwidth']

        return self._pretty_str(priorized_keys)


def remove_outlier(data, seg_twindow, th):
    scaler = StandardScaler()
    data_f = np.copy(data)
    data_standard = scaler.fit_transform(data.reshape(-1,1)).reshape(-1,)
    # searching outliers
    outs_position = np.where(data_standard > th)
    if outs_position:
        for i in list(map(int, outs_position[0])):
            f_0 = i - int(seg_twindow/10)
            if f_0 >= 0: # el outlier esta dentro de la ventana
                if f_0 + seg_twindow >= len(data_f):
                    # el final de la ventana esta fuera
                    data_f[f_0:] = np.nan
                else:
                    # el final de la ventana esta dentro
                    data_f[f_0:f_0+seg_twindow] = np.nan
            else:
                # el outlier esta al comienza de la ventana
                data_f[0:seg_twindow] = np.nan
    
    if data_f[~np.isnan(data_f)].size > 0:
        return np.nanmean(data_f)
    else:
        return np.nan


class LTE(object):
    def __init__(self, lte_file):

        if not os.path.isfile(lte_file):
            return ValueError(' lte_file do not found')
        
        self.lte_file = lte_file
        self.__set_stats__()
        self.__set_attrs__()


    def __str__(self):
        return self.stats.__str__()


    def is_attr(self, attr):
        return attr in self.stats.attributes


    def is_matrix(self, attr):
        if attr in SCALAR_PARAMS + SCALAR_PARAMS_OPT + SCALAR_PARAMS_F:
            return False

        if attr in VECTORAL_PARAMS + VECTORAL_PARAMS_OPT:
            return True


    def __set_stats__(self):
        with h5py.File(self.lte_file, "r") as f:
            hdr = f['header']
            lte_stats = dict(
                id = hdr.attrs['id'],
                channel = hdr.attrs['channel'],
                starttime = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                endtime = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                interval = hdr.attrs['interval'],
                int_olap = hdr.attrs['int_olap'],
                nro_time_bins = hdr.attrs['nro_time_bins'],
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
                PE_tau = hdr.attrs["pe_tau"],
                PE_order = hdr.attrs["pe_order"],
                time_bandwidth = hdr.attrs['time_bandwidth']
                )
        
        self.stats = LTEstats(lte_stats)


    def __set_attrs__(self):
        if self.stats.remove_response:
            attrs = SCALAR_PARAMS
        else:
            attrs = SCALAR_PARAMS[:-1]
        
        if self.stats.opt_params:
            attrs += SCALAR_PARAMS_OPT
                 
        if self.stats.f_params:
            attrs += SCALAR_PARAMS_F[0]

            if self.stats.remove_response:
                attrs += SCALAR_PARAMS_F[-1]

            if self.stats.opt_params:
                attrs += SCALAR_PARAMS_F[1:-1]

        if self.stats.polar_degree:
            attrs += VECTORAL_PARAMS
            
            if self.stats.polar_analysis:
                attrs += VECTORAL_PARAMS_OPT
        
        else:
            attrs += [VECTORAL_PARAMS[0]]
                
        self.stats.__add_attr__(attrs)


    def __empty_dict__(self):
        chan_dict = {
            'specgram':None,
            'fq_dominant':None,
            'fq_centroid':None,
            'energy':None,
            'pentropy':None,
            'rsam':None}
        
        if self.stats.remove_response:
            chan_dict['dsar']=None
        
        if self.stats.opt_params:
            empty_dict['vlf']=None
            empty_dict['lf']=None
            empty_dict['mf']=None
            empty_dict['hf']=None
            empty_dict['vlar']=None
            empty_dict['lrar']=None
            empty_dict['rmar']=None
        
        if self.stats.f_params:
            chan_dict['rsam_f']=None

            if self.stats.remove_response:
                chan_dict['dsar_f']=None
            
            if self.stats.opt_params:
                chan_dict['vlf_f']=None
                chan_dict['lf_f']=None
                chan_dict['mf_f']=None
                chan_dict['hf_f']=None
                chan_dict['vlar_f']=None
                chan_dict['lrar_f']=None
                chan_dict['rmar_f']=None
        
        empty_dict = {}
        for chan in self.stats.channel:
            empty_dict[chan] = chan_dict

        if self.stats.polar_degree:
            empty_dict['polar'] = {'degree': None}

            if self.stats.polar_analysis:
                empty_dict['polar']['rectlinearity'] = None
                empty_dict['polar']['azimuth'] = None
                empty_dict['polar']['elevation'] = None
        
        return empty_dict


    def __save__(self, dict_in, tbin):
        with h5py.File(self.lte_file, "r+") as f:
            for chan in self.stats.channel:
                chan_dict = dict_in[chan]
                chan_gr = f.get(chan)

                for attr, val in chan_dict.items():
                    group = get_group(attr)
                    dset = chan_gr.get(group)[attr]

                    if isinstance(val, type(None)):
                        val = np.nan
                    
                    if self.is_matrix(attr):
                        dset[tbin, :] = val
                    
                    else:
                        dset[tbin,] = val
            
            if self.stats.polar_degree:
                for attr, val in dict_in['polar'].items():
                    dset = f.get('polar')[attr]

                    if isinstance(val, type(None)):
                        val = np.nan
                    
                    dset[tbin, :] = val
            
            f.flush()


    def __compute__(self, station, days_in_memory=1, **kwargs):

        total_day = (self.stats.starttime - self.stats.endtime).total_seconds()/3600/24
        if total_day < days_in_memory:
            days_in_memory = -1
        
        if days_in_memory > 0:
            memory_delta = dt.timedelta(days=days_in_memory)
        else:
            memory_delta = dt.timedelta(days=total_day)
        
        total_bins = self.stats.nro_time_bins
        interval_delta = dt.timedelta(minutes=self.stats.interval)
        olap_delta = dt.timedelta(minutes=self.stats.interval*self.stats.int_olap)

        time = self.stats.starttime
        m_etime = time + memory_delta
        m_stream = station.get_stream(time, time + memory_delta, channel=self.stats.channel, remove_response=self.stats.remove_response, sample_rate=self.stats.sample_rate)

        psdkwargs = {
            'taper':kwargs.get('taper', False),
            'taper_p':kwargs.get('taper_p', 0.05),
            'time_bandwidth':self.stats.time_bandwidth
        }
        
        with tqdm(total=total_bins) as pbar:
            tbin = 0
            while time + interval_delta <= self.stats.endtime:
                time_text = time.strftime('%d %B %Y | %H:%M UTC  ::  ')
                int_endtime = time + interval_delta
                
                if int_endtime > m_etime:
                    m_etime = time + memory_delta
                    try:
                        m_stream = station.get_stream(time, time + memory_delta, channel=self.stats.channel, remove_response=self.stats.remove_response)
                    except:
                        m_stream = None
                
                if not m_stream:
                    pbar.write('\x1b[0;31;40m' + time_text + 'no data (stream error)' + '\x1b[0m')
                    # save nan

                else:
                    save_dict = self.__empty_dict__()

                    for tr in m_stream:
                        chan = tr.stats.channel

                        try:
                            chan_spec_params = tr.psd(starttime=time, endtime=int_endtime, fq_band=self.stats.fq_band, mov_avg_step=self.stats.step, olap_step=self.stats.step_olap, drm_params=True, **psdkwargs)
                        except:
                            chan_spec_params = None

                        if not chan_spec_params:
                            pbar.write('\x1b[0;31;40m' + time_text + tr.stats.channel + '(spectral error)' + '\x1b[0m')
                        
                        else:
                            data = tr.get_data(starttime=time, endtime=int_endtime, fq_band=self.stats.fq_band)
                            h = perm_entropy(data, order=self.stats.PE_order, delay=self.stats.PE_tau, normalize=True)
                            rsam = tr.get_data(starttime=time, endtime=int_endtime, fq_band=(2,4.5), abs=True)
                            
                            save_dict[chan]['specgram'] = chan_spec_params[0]
                            save_dict[chan]['fq_dominant'] = chan_spec_params[3]
                            save_dict[chan]['fq_centroid'] = chan_spec_params[2]
                            save_dict[chan]['energy'] = chan_spec_params[4]
                            save_dict[chan]['pentropy'] = h
                            save_dict[chan]['rsam'] = rsam.mean()

                            if self.stats.opt_params:
                                vlf = tr.get_data(starttime=time, endtime=int_endtime, fq_band=(.01,.1), abs=True)
                                lf = tr.get_data(starttime=time, endtime=int_endtime, fq_band=(.1,2), abs=True)
                                mf = tr.get_data(starttime=time, endtime=int_endtime, fq_band=(4,8), abs=True)
                                hf = tr.get_data(starttime=time, endtime=int_endtime, fq_band=(8,16), abs=True)
                                vlar = vlf/lf
                                lrar = lf/rsam
                                rmar = rsam/mf
                                
                                save_dict[chan]['vlf'] = vlf.mean()
                                save_dict[chan]['lf'] = lf.mean()
                                save_dict[chan]['mf'] = mf.mean()
                                save_dict[chan]['hf'] = hf.mean()
                                save_dict[chan]['vlar'] = vlar.mean()
                                save_dict[chan]['lrar'] = lrar.mean()
                                save_dict[chan]['rmar'] = rmar.mean()
                            
                            if self.stats.remove_response:
                                rr_tr = tr.remove_response2(output='DISP')
                                dst_mf = rr_tr.get_data(starttime=time, endtime=int_endtime, fq_band=(4,8), abs=True)
                                dst_hf = rr_tr.get_data(starttime=time, endtime=int_endtime, fq_band=(8,16), abs=True)
                                dsar = dst_mf/dst_hf
                                save_dict[chan]['dsar'] = dsar.mean()

                            if self.stats.f_params:
                                seg_time_window = kwargs.get('seg_tw', 150)
                                seg_twindow = int(seg_time_window*self.stats.sampling_rate)
                                rsam_f = np.array([remove_outlier(rsam, seg_twindow, self.stats.f_threshold)])
                                save_dict[chan]['rsam_f'] = rsam_f.mean()

                                if self.stats.opt_params:
                                    mf_f = np.array([remove_outlier(mf, seg_twindow, self.stats.f_threshold)])
                                    hf_f = np.array([remove_outlier(hf, seg_twindow, self.stats.f_threshold)])
                                    lf_f = np.array([remove_outlier(lf, seg_twindow, self.stats.f_threshold)])
                                    vlf_f = np.array([remove_outlier(vlf, seg_twindow, self.stats.f_threshold)])
                                    vlar_f = np.array([remove_outlier(vlar, seg_twindow, self.stats.f_threshold)])
                                    lrar_f = np.array([remove_outlier(lrar, seg_twindow, self.stats.f_threshold)])
                                    rmar_f = np.array([remove_outlier(rmar, seg_twindow, self.stats.f_threshold)])
                                    
                                    save_dict[chan]['vlf_f'] = vlf_f.mean()
                                    save_dict[chan]['lf_f'] = lf_f.mean()
                                    save_dict[chan]['mf_f'] = mf_f.mean()
                                    save_dict[chan]['hf_f'] = hf_f.mean()
                                    save_dict[chan]['vlar_f'] = vlar_f.mean()
                                    save_dict[chan]['lrar_f'] = lrar_f.mean()
                                    save_dict[chan]['rmar_f'] = rmar_f.mean()
                                
                                if self.stats.remove_response:
                                    dsar_f = np.array([remove_outlier(dsar, seg_twindow, self.stats.f_threshold)])
                                    save_dict[chan]['dsar_f'] = dsar_f.mean()
                    
                    # polarization parameters
                    if self.stats.polar_degree:
                        z_data = m_stream.get_component('Z').get_data(starttime=time, endtime=int_endtime, detrend=True)
                        n_data = m_stream.get_component('N').get_data(starttime=time, endtime=int_endtime, detrend=True)
                        e_data = m_stream.get_component('E').get_data(starttime=time, endtime=int_endtime, detrend=True)

                        npts_mov_avg = self.stats.sample_rate*self.stats.step*60
                        pa_kwargs = dict(
                            taper = kwargs.get('taper', False),
                            taper_p = kwargs.get('taper_p', 0.05),
                            time_bandwidth = self.stats.time_bandwidth
                            )
                        
                        try:
                            pa = PolarAnalysis(z_data, n_data, e_data, self.stats.sample_rate, npts_mov_avg=npts_mov_avg, olap=self.stats.step_olap, fq_band=self.stats.fq_band, full_analysis=self.stats.polar_analysis, **pa_kwargs)
                        
                        except ValueError:
                            pa = None
                        
                        if not pa:
                            pbar.write('\x1b[0;31;40m' + time_text + '(polar error)' + '\x1b[0m')

                        else:
                            save_dict['polar']['degree'] = pa.polar_dgr

                            if self.stats.polar_analysis:
                                save_dict['polar']['rectilinearity'] = pa.rect
                                save_dict['polar']['azimuth'] = pa.azimuth
                                save_dict['polar']['elevation'] = pa.elevation
                    
                self.__save__(save_dict, tbin)
                time = time + interval_delta - olap_delta
                tbin += 1
                pbar.update()


    def get(self, attr, starttime=None, endtime=None):
        
        if not starttime:
            n_0 = 0
        else:
            delta = (starttime - self.stats.starttime).total_seconds()/60
            n_0 = ((delta/self.stats.interval) - self.stats.int_olap) / (1 - self.stats.int_olap)
            n_0 = int(np.ceil(n_0))

        if not endtime:
            n_f = self.stats.nro_time_bins
        else:
            delta = (endtime - self.stats.starttime).total_seconds()/60
            n_f = ((delta/self.stats.interval) - self.stats.int_olap) / (1 - self.stats.int_olap)
            n_f = int(np.ceil(n_f))
        
        if self.is_attr(attr):
            group = get_group(attr)


        else:
            return None


    def check_list_attr(self, list_attr):
        if not all([attr in self.stats.attributes for attr in list_attr]):
            attrs = [at not in self.stats.attributes for at in list_attr]
            pos = list(filter(lambda x: attrs[x], range(len(list_attr))))
            not_availabel_attr = np.array(list_attr)[pos]
            print('warn: attributes %s not available' %not_availabel_attr)

            new_list_attr = [attr for attr in list_attr if attr in self.stats.attributes]
            if not new_list_attr:
                print('available attr: %s' % self.stats.attributes)

        else:
            new_list_attr = list_attr

        return new_list_attr


    def get_dataset(self, group, attr, index):
        with h5py.File(self.lte_file, "r") as f:
            if index[1]:
                if attr in ('specgram', 'degree'):
                    dset1 = f.get(group)[attr][index[0]:index[1],:]
                    dset2 = f.get('spectral')['freq'][:]
                    dset = (dset1, dset2)

                elif attr in ('dip', 'rect', 'azm') and self.stats.matrix_return:
                    dset1 = f.get('polar')[attr][index[0]:index[1],:]
                    dset2 = f.get('spectral')['freq'][:]
                    dset = (dset1, dset2)

                else:
                    dset = f.get(group)[attr][index[0]:index[1],]
            
            else:
                if attr in ('specgram', 'degree'):
                    dset1 = f.get(group)[attr][index[0],:]
                    dset2 = f.get('spectral')['freq'][:]
                    dset = (dset1, dset2)

                elif attr in ('dip', 'rect', 'azm') and self.stats.matrix_return:
                    dset1 = f.get('polar')[attr][index[0],:]
                    dset2 = f.get('spectral')['freq'][:]
                    dset = (dset1, dset2)

                else:
                    dset = f.get(group)[attr][index[0],]

        return dset


    def get_index(self, starttime=None, endtime=None):
        if starttime and endtime:
            st_diff = int((starttime - self.stats.starttime).total_seconds()/60)
            start_index = int(st_diff/self.stats.time_step)
            cond1 = start_index < 0 or start_index >= self.stats.time_bins
                
            ed_diff = int((endtime - self.stats.starttime).total_seconds()/60)
            end_index = int(ed_diff/self.stats.time_step)
            cond2 = end_index < 0 or end_index > self.stats.time_bins
            
            if cond1 or cond2:
                raise ValueError(' [0] Select a correct period of time ')

        if starttime and not endtime:
            st_diff = int((starttime - self.stats.starttime).total_seconds()/60)
            start_index = int(st_diff/self.stats.time_step)
            cond1 = start_index < 0 or start_index >= self.stats.time_bins
            end_index = None

            if cond1:
                raise ValueError(' [1] Select a correct period of time ')
        
        if not starttime and not endtime:
            start_index = 0
            end_index = self.stats.time_bins
        
        return (start_index, end_index)


    def get_attr(self, attr, starttime=None, endtime=None):
        """
        Get array from file
        :param starttime: datetime
        :param endtime: datetime
        :param attr: string.
        :return: numpy array
        """

        if attr not in self.__attrs__:
            raise ValueError(' attribute nor in %s' %self.__attrs__)

        gr = get_group(attr)
        (i, j) = self.get_index(starttime, endtime)
        dset = self.get_dataset(gr, attr, (i,j))

        return dset


    def get_time(self, starttime=None, endtime=None):
        i = self.get_index(starttime, endtime)
        times = [self.stats.starttime + dt.timedelta(minutes=k*self.stats.time_step) for k in range(i[0],i[1])]
        return np.array(times)


    def get_stats(self, attr, starttime=None, endtime=None):
        """
        Return (min,max,mean,mode)
        """

        if attr not in SCALAR_PARAMS:
            raise ValueError ('attr must be a scalar parameter')

        data = self.get_attr(attr, starttime=starttime, endtime=endtime)
        true_data = data[np.isfinite(data)]

        if attr == 'energy':
            data[np.where(data == 0)] = np.nan
            true_data = 10*np.log10(true_data)

        if true_data.shape[0] > 2:
            v_min = true_data.min()
            v_max = true_data.max()
            v_mean = true_data.mean()

            x_range = np.linspace(v_min, v_max, 500)
            gkde = gaussian_kde(true_data)
            kde = gkde(x_range)
            v_mode = x_range[np.argmax(kde)]

            return(v_min, v_max, v_mean, v_mode)

        else:
            print(' sample size is less than 2')
            return None


    def get_dict_stats(self, list_attr, starttime=None, endtime=None):
        dout = {}
        true_list = self.check_list_attr(list_attr, return_list=True)

        for att in true_list:
            if att in SCALAR_PARAMS:
                dout[att] = self.get_stats(att, starttime=starttime, endtime=endtime)

        return dout


    def get_peaks(self, fq_range=(), peak_thresholds={}):
        """

        Return peaks as describes in Melchor et al 2021

        Parameters
        ----------
        fq_range : tuple, optional
            by default stats.fq_band
        peak_thresholds : dict, optional
            by default {'fq_delta': 0.05, 'spec_th': 0.95, 'pd_th': 0.8, 'rect_th': 0.7, 'pd_std':0.1, 'r_std':0.1, 'azm_std':10, 'dip_std':10}

        Returns
        -------
        [Peak object]
        """

        if not self.stats.polargram or not self.stats.matrix_return:
            raise ValueError(' err: No polarization data in LTE')

        self.peak_thresholds = default_peak_thresholds
        
        # define parameters
        for key in default_peak_thresholds.keys():
            if peak_thresholds.get(key, None):
                self.peak_thresholds[key] = peak_thresholds.get(key)
                
        # plot model parameters
        print(f'\n  File: {os.path.basename(self.lte_file)}')
        print('  ----- Model param ------')
        for key, item in self.peak_thresholds.items():
            print(f'  {key:>10}     ', item)
        print('  ------------------------\n')

        if not fq_range:
            fq_range = self.stats.fq_band

        # erg = self.get_attr('energy')
        erg = self.get_dataset('energy', 'energy', index=(0, self.stats.time_bins))
        sxx, fq = self.get_attr('specgram')

        # do nan when energy in zero
        sxx[np.where(erg == 0), :] = np.nan

        # polar atributes
        sxx = 10*np.log10(sxx)
        pd, _ = self.get_attr('degree')
        rect, _ = self.get_attr('rect')
        thetaH, _ = self.get_attr('azm')
        thetaV, _ = self.get_attr('dip')

        # correct 180 degrees in azimuth
        thetaH[np.where(thetaH < 0)] = thetaH[np.where(thetaH < 0)] + 180

        # get fq_npts
        fq0 = np.argmin(np.abs(np.array(fq)-fq_range[0]))
        fq1 = np.argmin(np.abs(np.array(fq)-fq_range[1]))

        peaks = []
        for t in range(sxx.shape[0]):
            print(f'     getting peaks {int((t+1)*100/sxx.shape[0]):3}%', end='\r')

            sxx_t = sxx[t, fq0:fq1]

            if not sxx_t.all() or np.isnan(sxx_t).any() or not np.isfinite(sxx_t).all():
                peaks.append(None)
            
            else:
                sxx_t_norm = MinMaxScaler().fit_transform(sxx_t.reshape(-1,1))
                peaks_pos,_ = signal.find_peaks(sxx_t_norm.reshape(-1,), height=self.peak_thresholds['spec_th'], distance=5)
                peaks_info = dict(fq=[], sp=[], pd=[], r=[], azm=[], dip=[])

                for p in peaks_pos:
                    fp = fq[fq0:fq1][p]
                    sp = sxx_t[p]
                    peaks_info['fq'] += [fp]
                    peaks_info['sp'] += [sp]

                    pd_ans = np.nan
                    r_ans = np.nan
                    azm_ans = np.nan
                    dip_ans = np.nan

                    fp_0 = fp - self.peak_thresholds['fq_delta']
                    fp_1 = fp + self.peak_thresholds['fq_delta']
                    fq_npts = [
                        np.argmin(np.abs(np.array(fq)-fp_0)), 
                        np.argmin(np.abs(np.array(fq)-fp_1))
                        ]
                    
                    pd_peak = pd[t, fq_npts[0]:fq_npts[1]]
                    r_peak = rect[t, fq_npts[0]:fq_npts[1]]
                    azm_peak = thetaH[t, fq_npts[0]:fq_npts[1]]
                    dip_peak = thetaV[t, fq_npts[0]:fq_npts[1]]

                    pp_avg, pp_std = pd_peak.mean(), pd_peak.std()

                    if pp_avg > self.peak_thresholds['pd_th'] and pp_std < self.peak_thresholds['pd_std']:
                        pd_ans = pp_avg
                        r = r_peak[np.where(pd_peak > self.peak_thresholds['pd_th'])]
                        azm = azm_peak[np.where(pd_peak > self.peak_thresholds['pd_th'])]
                        dip = dip_peak[np.where(pd_peak > self.peak_thresholds['pd_th'])]
                        r_avg, r_std = r.mean(), r.std()

                        if r_std < self.peak_thresholds['r_std']:
                            r_ans = r_avg

                            if r_avg > self.peak_thresholds['rect_th'] and any(r[r > self.peak_thresholds['rect_th']]):
                                azm_avg = azm[np.where(r > self.peak_thresholds['rect_th'])].mean()
                                azm_std = azm[np.where(r > self.peak_thresholds['rect_th'])].std()

                                if azm_std <= self.peak_thresholds['azm_std']:
                                    azm_ans = azm_avg

                                dip_avg = dip[np.where(r > self.peak_thresholds['rect_th'])].mean()
                                dip_std = dip[np.where(r > self.peak_thresholds['rect_th'])].std()

                                if dip_std <= self.peak_thresholds['dip_std']:
                                    dip_ans = dip_avg
                    
                    peaks_info['pd'] += [pd_ans]
                    peaks_info['r'] += [r_ans]
                    peaks_info['azm'] += [azm_ans]
                    peaks_info['dip'] += [dip_ans]

                if not peaks_info['fq']:
                    peaks.append(None)
                else:
                    peaks.append(peaks_info)

        n_pks = len(list(filter(lambda a: a != None, peaks)))
        print(f'\n     Bins with peaks info: {int(100*n_pks/len(peaks)):3}%')
        
        return peaks


    def plot(self, starttime=None, endtime=None, interval=None, plot_attr=[], return_fig=False, **kwargs):
        from seisvo.gui.glte import plot_gui, get_fig

        if not starttime:
            starttime = self.stats.starttime

        else:
            if starttime < self.stats.starttime:
                starttime = self.stats.starttime

        if not endtime:
            if interval:
                endtime = starttime + timedelta(days=interval)

            else:
                endtime = self.stats.endtime
                interval = (endtime - starttime).days
        else:
            if endtime > self.stats.endtime:
                endtime = self.stats.endtime
            interval = (endtime - starttime).days


        if not plot_attr:
            list_attr = ['energy', 'pentropy', 'fq_dominant', 'fq_centroid', 'specgram']

            if self.stats.polargram:
                list_attr += ['fq_polar', 'degree_max', 'degree_wavg', 'degree']

        else:
            if set(plot_attr).issubset(set(LIST_ATTR)):
                list_attr = plot_attr
            else:
                raise ValueError(' Error reading attributes.')

        if return_fig:
            fig, axes = get_fig(self, starttime, endtime, interval, list_attr, **kwargs)

            return fig, axes
        
        else:
            plot_gui(self, starttime, endtime, interval, list_attr, **kwargs)


    def get_pdf(self, attr, starttime=None, endtime=None, bandwidth=0.01, **kwargs):
        """
        Compute the PDF
        """

        y_size = kwargs.get('y_space', 1000)

        if self.is_matrix(attr):
            matrix, x_space = self.get_attr(attr, starttime=starttime, endtime=endtime)
            masked_matrix = np.ma.masked_invalid(matrix.T)
            matrix = np.ma.compress_cols(masked_matrix)
            
            if attr == 'specgram':
                matrix = 10*np.log10(matrix.T)
            
            y_min = kwargs.get('y_min', matrix.min())
            y_max = kwargs.get('y_max', matrix.max())
            y_space = np.linspace(y_min, y_max, y_size).reshape(y_size, 1)
            pdf = get_PDF(matrix, y_space, bandwidth, **kwargs)

            return x_space, y_space.reshape(y_size,), pdf
        
        else:
            data = self.get_attr(attr, starttime=starttime, endtime=endtime)
            y_min = kwargs.get('y_min', data.min())
            y_max = kwargs.get('y_max', data.max())
            y_space = np.linspace(y_min, y_max, y_size).reshape(y_size, 1)
            pdf = get_PDF(data.reshape(-1,1), y_space, bandwidth, **kwargs)

            return y_space.reshape(y_size,), pdf


    def get_Peaks(self, fq_range=(), peak_thresholds={}):
        return Peaks(self, fq_range=fq_range, peak_thresholds=peak_thresholds)


    @staticmethod
    def new(station, lte_file, headers, auto_chunk=True, **ltekwargs):
        """
        Create new LTE (hdf5) file
        """

        f = h5py.File(lte_file, "w-")
        
        # header dataset
        hdr = f.create_dataset('header',(1,))
        
        hdr.attrs['id'] = headers['id']
        hdr.attrs['channel'] = headers['channel']
        hdr.attrs['starttime'] = headers['starttime']
        hdr.attrs['endtime'] = headers['endtime']
        hdr.attrs['interval'] = headers['interval']
        hdr.attrs['int_olap'] = headers['int_olap']
        hdr.attrs['step'] = headers['step']
        hdr.attrs['step_olap'] = headers['step_olap']
        hdr.attrs['sample_rate'] = headers['sample_rate']
        hdr.attrs['remove_response'] = headers['remove_response']
        hdr.attrs['fq_band'] = headers['fq_band']
        hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
        hdr.attrs['nro_freq_bins'] = headers['nro_freq_bins']
        hdr.attrs['polar_degree'] = headers['polar_degree']
        hdr.attrs['polar_analysis'] = headers['polar_analysis']
        hdr.attrs['f_params'] = headers['f_params']
        hdr.attrs['opt_params'] = headers['opt_params']
        hdr.attrs['f_threshold'] = headers['f_threshold']
        hdr.attrs['PE_tau'] = headers['PE_tau']
        hdr.attrs['PE_order'] = headers['PE_order']
        hdr.attrs['time_bandwidth'] = headers['time_bandwidth']

        freqbins = headers['nro_freq_bins']
        timebins = headers['nro_time_bins']

        # set chunks 
        if auto_chunk:
            chunk_shape1 = chunk_shape2 = True
            chunk_info = 'auto'
        
        else:
            if timebins > 1000:
                chunk_info = 200
                chunk_shape1 = (200, freqbins)
                chunk_shape2 = (200,)
            
            elif timebins > 10000:
                chunk_info = 500
                chunk_shape1 = (500, freqbins)
                chunk_shape2 = (500,)
            
            elif timebins > 100000:
                chunk_info = 1000
                chunk_shape1 = (1000, freqbins)
                chunk_shape2 = (1000,)

            else:
                chunk_shape1 = chunk_shape2 = None
                chunk_info = 'none'

        # print shape info
        print(' LTE file INFO')
        print(' -------------')
        print(" hdf5_memory info: %s " % lte_file)
        print(' --- dataset size: ', (timebins, freqbins))
        print(' --- chunk size: ', chunk_info)
        print('')
        
        # print info
        print(' LTE stats:')
        for info_key in ['id', 'channel', 'starttime', 'endtime' ,'interval' ,'int_olap' ,'step' ,'step_olap' ,'sample_rate' ,'remove_response' ,'fq_band' ,'nro_time_bins' ,'nro_freq_bins' ,'polar_degree' ,'polar_analysis' ,'f_params', 'opt_params', 'f_threshold' ,'PE_tau' ,'PE_order', 'time_bandwidth']:
            print(f' {info_key}:  {headers[info_key]}')

        for chan in headers['channel']:
            chgr = f.create_group(chan)

            spec = chgr.create_group("spectral")
            spec.create_dataset('freq', (freqbins,), dtype=np.float32)
            spec.create_dataset('specgram', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)
            spec.create_dataset('fq_dominant', (timebins,), chunks=chunk_shape2, dtype=np.float32)
            spec.create_dataset('fq_centroid', (timebins,), chunks=chunk_shape2, dtype=np.float32)
            spec.create_dataset('energy', (timebins,), chunks=chunk_shape2, dtype=np.float32)

            amp = chgr.create_group("amplitude")
            amp.create_dataset('pentropy', (timebins,), chunks=chunk_shape2, dtype=np.float32)
            amp.create_dataset('rsam', (timebins,), chunks=chunk_shape2, dtype=np.float32)

            if headers['remove_response']:
                amp.create_dataset('dsar', (timebins,), chunks=chunk_shape2, dtype=np.float32)

            if headers['opt_params']:
                amp.create_dataset('vlf', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                amp.create_dataset('lf', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                amp.create_dataset('mf', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                amp.create_dataset('hf', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                amp.create_dataset('vlar', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                amp.create_dataset('lrar', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                amp.create_dataset('rmar', (timebins,), chunks=chunk_shape2, dtype=np.float32)

            if headers['f_params']:
                ampF = chgr.create_group("amplitude_f")
                ampF.create_dataset('rsam_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                
                if headers['remove_response']:
                    ampF.create_dataset('dsar_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                
                if headers['opt_params']:
                    ampF.create_dataset('mf_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                    ampF.create_dataset('hf_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                    ampF.create_dataset('lf_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                    ampF.create_dataset('vlf_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                    ampF.create_dataset('vlar_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                    ampF.create_dataset('lrar_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
                    ampF.create_dataset('rmar_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)

        if headers['polar_degree']:
            polar = f.create_group("polar")
            polar.create_dataset('degree', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)

            if headers['polar_analysis']:
                polar.create_dataset('rectlinearity', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)
                polar.create_dataset('elevation', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)
                polar.create_dataset('azimuth', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)
                
        f.flush()
        f.close()

        lte = LTE(lte_file)
        lte.__compute__(station, ltekwargs)

        return lte


class Peaks(LTE):
    def __init__(self, lte, fq_range=(), peak_thresholds={}):

        if isinstance(lte, str):
            super().__init__(lte)
        
        elif isinstance(lte, LTE):
            super().__init__(lte.lte_file)
        else:
            raise ValueError ('lte should be string with lte file path or lte object')

        self.all_peaks = super().get_peaks(fq_range, peak_thresholds)
        self.total_dominant_peaks_ = sum([len(i['fq']) for i in self.all_peaks if i])
        
        self.peaks_ = {}
        self.nro_ = 0
        self.df_ = None


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

        return self


    def __fit_fq_kde__(self, **kwargs):
        # compute KDE for peak distribution weighted by spectral energy
        
        fq_l = [item['fq'] for item in self.all_peaks if item]
        fqa = np.array(list(chain(*fq_l)))

        bandwidth = kwargs.get('bandwidth', {})
        self.bandwidth_fc = bandwidth.get('frec', default_bandwidth['frec'])

        v_min = kwargs.get('v_min', None)
        if not v_min:
            v_min = fqa.min()
            
        v_max = kwargs.get('v_max', None)
        if not v_max:
            v_max = fqa.max()
        
        v_space = kwargs.get('v_space', 1000) # fqa.shape[0]*5
        self.fq_space_ = np.linspace(v_min, v_max, v_space).reshape(v_space, 1)

        if kwargs.get('weighted', False):
            sp_l = [item['sp'] for item in self.all_peaks if item]
            spa = np.array(list(chain(*sp_l)))
            self.weights = spa/spa.max()

        else:
            self.weights = None

        kde = get_KDE(fqa.reshape(-1, 1), self.bandwidth_fc, weight=self.weights)

        self.fq_pdf_ = np.exp(kde.score_samples(self.fq_space_))


    def __fit_polar__(self, **kwargs):
        self.dist_throld = kwargs.get('dist_throld', 0.5)
        self.n_sample_min = kwargs.get('n_sample_min', 10)
        
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

        for f in range(1, self.nro_+1):
            dfq = self.peaks_[f]['fq']
            sp, pd, rect, th_H, th_V = [], [], [], [], []
            for item in self.all_peaks:
                if item:
                    for n, fq in enumerate(item['fq']):
                        if dfq-self.peaks_[f]['width'] <= fq <= dfq+self.peaks_[f]['width']:
                            sp += [item['sp'][n]]
                            pd += [item['pd'][n]]
                            rect += [item['r'][n]]
                            th_H += [item['azm'][n]]
                            th_V += [item['dip'][n]]
            
            sp = np.array(sp)[np.where(~np.isnan(sp))]
            
            if len(sp) < 10*self.n_sample_min:
                continue
            
            sp_kde = get_KDE(sp.reshape(-1,1), self.bandwidth['sp'])

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

            pd = np.array(pd)[np.where(~np.isnan(pd))]
            if pd.shape[0] > self.n_sample_min:
                rect = np.array(rect)[np.where(~np.isnan(rect))]
                
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

                    if rect_range[1] >= self.peak_thresholds['rect_th']:
                        th_H = np.array(th_H)[np.where(~np.isnan(th_H))]
                        th_V = np.array(th_V)[np.where(~np.isnan(th_V))]
                        if th_H.shape[0] > self.n_sample_min and th_V.shape[0] > self.n_sample_min:
                            
                            # save Cl
                            self.peaks_[f]['cl'] = (len(th_H) + len(th_V)) / (2*len(rect))
                        
                            thH_kde = get_KDE(th_H.reshape(-1,1), self.bandwidth['angle'])
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

                            thV_kde = get_KDE(th_V.reshape(-1,1), self.bandwidth['angle'])
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


    def plot_time_evo(self, fq_range=(), out='time_evo', format='pdf', spj=5, rj=3, pj=1):
        sup.plot_lte_peaks_evo(self, fq_range=fq_range, out=out, format=format, spj=spj, rj=rj, pj=pj)


    def plot_peak_evo(self, nro, fq_off=0.1, pd_throld=0.8, r_throld=0.75, out=None, out_dir='./', format='pdf', show=False):
        sup.plot_lte_peak_evo(self, self.peaks[nro]['fq'], fq_off=fq_off, pd_throld=pd_throld, r_throld=r_throld, out=out, out_dir=out_dir, format=format, show=show)
        

    def plot_spec_pdf(self, plot=True, **kwargs):
        if self.nro_ == 0:
            raise ValueError ('no dominant frequencies to plot!')

        fig = sup.plot_peaks_spec_pdf(self, plot=plot, **kwargs)

        return fig


    def plot_peak_pdf(self, n, plot=True, **kwargs):
        if self.nro_ == 0:
            raise ValueError ('no dominant frequencies to plot!')
        
        fig = sup.plot_peaks_pdf(self, n, plot=plot, **kwargs)
        
        return fig


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
    
    
    def to_json(self, out=None, out_path='./'):
        if isinstance(self.df_, pd.DataFrame):
            if not out:
                lte_file = os.path.basename(self.lte_file).split('.')
                if lte_file[-1] == 'lte':
                    lte_file[-1] = 'json'
                else:
                    lte_file += ['json']
                out = '.'.join(lte_file)
            
            json_file = os.path.join(out_path, out)

            if os.path.isfile(json_file):
                os.remove(json_file)

            self.df_.to_json(json_file)
