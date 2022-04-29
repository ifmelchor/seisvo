#!/usr/bin/python3
# coding=utf-8

'''

Read write and operate .lte files

'''

import os
import h5py
import numpy as np
import pandas as pd
from scipy import signal
from scipy.stats import gaussian_kde
from antropy import perm_entropy
from pytictoc import TicToc
from seisvo.plotting import plot_gram
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from obspy.core.util.attribdict import AttribDict
from datetime import datetime, timedelta
from itertools import chain

# import seisvo.utils.plotting as sup
from seisvo.core.network import Network
from seisvo.signal.pdf import get_PDF, get_KDE

SCALAR_PARAMS = [
    'fq_dominant',
    'fq_centroid',
    'energy',
    'pentropy',
    'rsam',
    'mf',
    'hf',
    'dsar',
    'rsam_f',
    'mf_f',
    'hf_f',
    'dsar_f'
]

SCALAR_OPT_PARAMS = [
    'dsar',
    'dsar_f',
]

VECTORAL_PARAMS = [
    'specgram',
    'degree'
    ]

VECTORAL_OPT_PARAMS = [
    'dip',
    'rect',
    'azm'
    ]

AMP_PARAMS = [
    'rsam',
    'mf',
    'hf',
    'pentropy',
    'dsar',
    'rsam_f',
    'mf_f',
    'hf_f',
    'dsar_f'
]

SPEC_PARAMS = [
    'specgram',
    'fq_dominant',
    'fq_centroid',
    'energy'
]

POLAR_PARAMS = [
    'degree',
    'dip',
    'azm',
    'rect',
    'fq_polar',
    'degree_max',
    'degree_wavg'
]


def get_group(attr):
        if attr in AMP_PARAMS:
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

    def add_atributtes(self, list):
        self.attributes = list

    def __str__(self):
        priorized_keys = [
        'id',
        'sampling_rate',
        'starttime',
        'endtime',
        'time_step',
        'time_bins',
        'freq_bins',
        'remove_response',
        'fq_band',
        'avg_step',
        'threshold',
        ]
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
        self.lte_file = lte_file
        self.__set_stats__()

        self.__groups__ = (
            'amplitude', 
            'spectral', 
            'polar'
            )

        self.__attrs__ = SCALAR_PARAMS
        
        if self.stats.polargram:
            self.__attrs__ += VECTORAL_PARAMS
            
            if self.stats.matrix_return:
                self.__attrs__ += VECTORAL_OPT_PARAMS
        
        else:
            self.__attrs__ += [VECTORAL_PARAMS[0]]
        
        if self.stats.remove_response:
            self.__attrs__ += SCALAR_OPT_PARAMS
        
        self.stats.add_atributtes(self.__attrs__)


    def __str__(self):
        return self.stats.__str__()


    def is_matrix(self, attr):
        if attr in SCALAR_PARAMS:
            return False

        if attr in VECTORAL_PARAMS:
            return True


    def __set_stats__(self):
        with h5py.File(self.lte_file, "r") as f:
            hdr = f['header']
            self.stats = LTEstats({
                'id': hdr.attrs['id'],
                'sampling_rate' : hdr.attrs['sampling_rate'],
                'starttime': datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                'endtime': datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                'time_step': int(hdr.attrs['time_step']),
                'time_bins': int(hdr.attrs['time_bins']),
                'freq_bins': int(hdr.attrs['freq_bins']),
                'remove_response': bool(hdr.attrs['remove_response']),
                # 'th_outlier': float(hdr.attrs['th_outlier']),
                'polargram': bool(hdr.attrs['polargram']),
                'matrix_return': hdr.attrs['matrix_return'],
                'fq_band': hdr.attrs['fq_band'],
                'avg_step': hdr.attrs['avg_step'],
                'p_order': hdr.attrs['p_order'],
                'tau':hdr.attrs['tau'],
                'threshold':hdr.attrs['threshold']
                })

        net = self.stats.id.split('.')[0]
        sta = self.stats.id.split('.')[1]
        loc = self.stats.id.split('.')[2]

        self.station = Network(net).get_sta(sta, loc=loc)


    def check_list_attr(self, list_attr, return_list=True):
        if not all([attr in self.__attrs__ for attr in list_attr]):
            attrs = [at not in self.__attrs__ for at in list_attr]
            pos = list(filter(lambda x: attrs[x], range(len(list_attr))))
            not_availabel_attr = np.array(list_attr)[pos]
            print('warn: attributes %s not available' %not_availabel_attr)

            new_list_attr = [attr for attr in list_attr if attr in self.__attrs__]
            if not new_list_attr:
                print('available attr: %s' % self.__attrs__)

        else:
            new_list_attr = list_attr

        if return_list:
            return new_list_attr
        
        else:
            return bool(new_list_attr)


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
        times = [self.stats.starttime + timedelta(minutes=k*self.stats.time_step) for k in range(i[0],i[1])]
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


    def __save_data__(self, group, attr, data, item):
        """
        Save data to a hdf5 file
        """

        with h5py.File(self.lte_file, "r+") as f:

            dset = f.get(group)[attr]

            if group in ('spectral', 'polar'):
                if attr == 'freq':
                    dset[:,] = data[:,]

                elif attr in ('degree', 'specgram'):
                    dset[item, :] = data[0, :]

                elif attr in ('dip', 'azm', 'rect') and self.stats.matrix_return:
                    dset[item, :] = data[0, :]

                else:
                    dset[item,] = data[0,]

            else:
                dset[item,] = data[0,]

            f.flush()


    def __save_dict__(self, dict_in, tbin, nan=False):
        if nan:
            nan_data = np.array([np.nan]).reshape(1,)
            nan_array = np.full((1, self.stats.freq_bins), np.nan)

            dict_in = {}
            for attr in self.__attrs__:
                if attr in SPEC_PARAMS:
                    dict_in[attr] = nan_data
                if attr in VECTORAL_PARAMS:
                    dict_in[attr] = nan_array

        for attr in dict_in.keys():
            gr = self.get_group(attr)
            self.__save_data__(gr, attr, dict_in[attr], tbin)


    def __compute__(self, **kwargs):
        """ 
        Compute parameters.
        """

        start_time = self.stats.starttime
        time_delta = timedelta(minutes=self.stats.time_step)
        t = TicToc()

        if self.stats.avg_step == self.stats.time_step:
            avg_step = None
        else:
            avg_step = self.stats.avg_step

        for tbin in range(self.stats.time_bins):
            t.tic()
            save_dict = {}
            step_pcent = int(100*tbin/self.stats.time_bins)
            str_time = start_time.strftime('%Y-%m-%d %H:%M:%S')
            status_text = None
            freq = None

            # channel
            chan = self.stats.id.split('.')[-1]
            
            try:
                st3 = self.station.get_stream(start_time, start_time + time_delta, 
                        remove_response=self.stats.remove_response, sample_rate=self.stats.sampling_rate)
                trace = st3.get(component=chan[-1])

                # if self.stats.polargram:
                #     st3 = self.station.get_stream(start_time, start_time + time_delta, 
                #         remove_response=self.stats.remove_response, sample_rate=self.stats.sampling_rate)
                #     trace = st3.get_component(chan[-1])
                #     # start = max([tr.stats.starttime for tr in st3])
                #     # end = min([tr.stats.endtime for tr in st3])
                
                # else:
                #     st = self.station.get_stream(start_time, start_time + time_delta, 
                #         remove_response=self.stats.remove_response, sample_rate=self.stats.sampling_rate,
                #         component=chan[-1])
                #     # start = start_time
                #     # end = start_time + time_delta
                #     trace = st[0]
                # stats = trace.stats

            except:
                status_text = '\x1b[0;31;40m' + 'Fail 1 (stream)' + '\x1b[0m'
                self.__save_dict__(save_dict, tbin, nan=True)
                t.toc(' %d%% ... Time: %s --- Status: : %s' % (step_pcent, str_time, status_text))
                start_time += time_delta
                continue

            try:
                spec = trace.psd(avg_step=avg_step, plot=False, return_fig=False, opt_return=True, **kwargs)
                
                save_dict['fq_dominant'] = np.array([spec[2][0]]).reshape(1,)
                save_dict['fq_centroid'] = np.array([spec[2][1]]).reshape(1,)
                save_dict['specgram'] = spec[0].reshape(1, self.stats.freq_bins)

                rsam = np.abs(trace.get_data(detrend=True, fq_band=(2,4.5)))
                mf = np.abs(trace.get_data(detrend=True, fq_band=(4,8)))
                hf = np.abs(trace.get_data(detrend=True, fq_band=(8,16)))
                seg_twindow = int(150*self.stats.sampling_rate)
                
                save_dict['rsam'] = rsam.mean().reshape(1,)
                save_dict['rsam_f'] = np.array([remove_outlier(rsam, seg_twindow, self.stats.th_outlier)]).reshape(1,)

                save_dict['mf'] = mf.mean().reshape(1,)
                save_dict['mf_f'] = np.array([remove_outlier(mf, seg_twindow, self.stats.th_outlier)]).reshape(1,)

                save_dict['hf'] = hf.mean().reshape(1,)
                save_dict['hf_f'] = np.array([remove_outlier(hf, seg_twindow, self.stats.th_outlier)]).reshape(1,)

                if self.stats.remove_response:
                    dst_mf = self.station.get_stream(start_time, start_time + time_delta, component='Z',
                        remove_response='DISP', sample_rate=self.stats.sampling_rate, prefilt=(4,8))
                    
                    dst_hf = self.station.get_stream(start_time, start_time + time_delta, component='Z',
                        remove_response='DISP', sample_rate=self.stats.sampling_rate, prefilt=(8,16))
                    
                    dsar = np.abs(dst_mf[0].data).mean()/np.abs(dst_hf[0].data).mean()
                    save_dict['dsar'] = dsar.mean().reshape(1,)
                    save_dict['dsar_f'] = np.array([remove_outlier(dsar, seg_twindow, self.stats.th_outlier)]).reshape(1,)

                rms = trace.rms(avg_step=avg_step, fq_bands=self.stats.fq_band, **kwargs)
                save_dict['energy'] = np.array([rms[0]]).reshape(1,)
            
                tr_data = trace.get_data(detrend=True, fq_band=self.stats.fq_band)
                h = [perm_entropy(tr_data.data, order=self.stats.p_order, delay=self.stats.tau, normalize=True)]
                save_dict['pentropy'] = np.array(h).reshape(1,)

            except:
                status_text = '\x1b[0;31;40m' + 'Fail 2 (spec)' + '\x1b[0m'
                self.__save_dict__(save_dict, tbin, nan=True)
                t.toc(' %d%% ... Time: %s --- Status: : %s' % (step_pcent, str_time, status_text))
                start_time += time_delta
                continue

            if np.isnan([save_dict[attr] for attr in ('fq_dominant', 'energy', 'pentropy')]).any():
                status_text = '\x1b[0;31;40m' + 'Fail 3 (spec)' + '\x1b[0m'
                self.__save_dict__(save_dict, tbin, nan=True)
                t.toc(' %d%% ... Time: %s --- Status: : %s' % (step_pcent, str_time, status_text))
                start_time += time_delta
                continue

            if self.stats.polargram:
                try:
                    with np.errstate(all='raise'):
                        polar_ans = st3.polardegree(avg_step=avg_step, olap=0, opt_return=True, 
                                         plot=False, threshold=self.stats.threshold, 
                                         matrix_return=self.stats.matrix_return, **kwargs)

                        save_dict['degree'] = polar_ans[0].reshape(1, self.stats.freq_bins)

                        if self.stats.matrix_return:
                            save_dict['azm'] = polar_ans[3].reshape(1, self.stats.freq_bins)
                            save_dict['dip'] = polar_ans[4].reshape(1, self.stats.freq_bins)
                            save_dict['rect'] = polar_ans[5].reshape(1, self.stats.freq_bins)
                
                except:
                    status_text = '\x1b[0;31;40m' + 'Fail 4 (polr)' + '\x1b[0m'
                    nan_array = np.full((1, self.stats.freq_bins), np.nan)

                    save_dict['degree'] = nan_array
                    if self.stats.matrix_return:
                        save_dict['azm'] = save_dict['dip'] = save_dict['rect'] = nan_array

            if not status_text:
                status_text = '\x1b[0;32;40m' + '{0:^13}'.format('Ok') + '\x1b[0m'

            if not freq:
                freq = spec[1].reshape(self.stats.freq_bins,)
                self.__save_data__("spectral", 'freq', freq, item=0)

            self.__save_dict__(save_dict, tbin)
            t.toc(' %d%% ... Time: %s --- Status: : %s' % (step_pcent, str_time, status_text))
            start_time += time_delta

        print('')
        print('\n Data processed successfully!')
        print('')


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
    def new(lte_file, headers, **kwargs):
        """
        Create new LTE (hdf5) file
        """

        f = h5py.File(lte_file, "w-")
        # header dataset
        hdr = f.create_dataset('header',(1,))
        hdr.attrs['id'] = headers['id']
        hdr.attrs['sampling_rate'] = headers['sampling_rate']
        hdr.attrs['starttime'] = headers['starttime']
        hdr.attrs['endtime'] = headers['endtime']
        hdr.attrs['time_step'] = headers['time_step']
        hdr.attrs['time_bins'] = headers['time_bins']
        hdr.attrs['time_bandwidth'] = headers['time_bandwidth']
        hdr.attrs['avg_step'] = headers['avg_step']
        hdr.attrs['remove_response'] = headers['remove_response']
        hdr.attrs['th_outlier'] = headers['th_outlier']
        hdr.attrs['fq_band'] = headers['fq_band']
        hdr.attrs['freq_bins'] = headers['freq_bins']
        hdr.attrs['p_order'] = headers['p_order']
        hdr.attrs['tau'] = headers['tau']
        hdr.attrs['threshold'] = headers['threshold']
        hdr.attrs['polargram'] = headers['polargram']
        hdr.attrs['matrix_return'] = headers['matrix_return']

        freqbins = headers['freq_bins']
        timebins = headers['time_bins']

        # set chunks
        time0 = datetime.strptime(headers['starttime'], '%Y-%m-%d %H:%M:%S')
        time1 = datetime.strptime(headers['endtime'], '%Y-%m-%d %H:%M:%S')
        nro_days = (time1 - time0).days
        day_bin = int(60/(headers['time_step'])*24)

        if nro_days > 10:
            chunk_shape1 = (5*day_bin, freqbins)
            chunk_shape2 = (5*day_bin,)

        elif 2 < nro_days <= 10:
            chunk_shape1 = (day_bin, freqbins)
            chunk_shape2 = None

        else:
            chunk_shape1 = chunk_shape2 = None

        # print shape info
        print(" hdf5_memory info: %s " % lte_file)
        print(' --- all datasets sizes are ', (timebins, freqbins))
        print(' --- saving in memory chunk sizes of ', chunk_shape1)
        print('')

        # spectral datasets
        spec = f.create_group("spectral")
        spec.create_dataset('freq', (freqbins,), dtype=np.float32)
        spec.create_dataset('specgram', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)
        spec.create_dataset('fq_dominant', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        spec.create_dataset('fq_centroid', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        spec.create_dataset('energy', (timebins,), chunks=chunk_shape2, dtype=np.float32)

        # polargram datasets
        if headers['polargram']:
            polar = f.create_group("polar")
            polar.create_dataset('degree', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)

            if headers['matrix_return']:
                polar.create_dataset('dip', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)
                polar.create_dataset('azm', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)
                polar.create_dataset('rect', (timebins, freqbins), chunks=chunk_shape1, dtype=np.float32)

        # amplitude
        amp = f.create_group("amplitude")
        amp.create_dataset('pentropy', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        amp.create_dataset('rsam', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        amp.create_dataset('rsam_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        amp.create_dataset('mf', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        amp.create_dataset('mf_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        amp.create_dataset('hf', (timebins,), chunks=chunk_shape2, dtype=np.float32)
        amp.create_dataset('hf_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)

        if headers['remove_response']:
            amp.create_dataset('dsar', (timebins,), chunks=chunk_shape2, dtype=np.float32)
            amp.create_dataset('dsar_f', (timebins,), chunks=chunk_shape2, dtype=np.float32)

        f.flush()
        f.close()

        lte = LTE(lte_file)
        lte.__compute__(**kwargs)

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
