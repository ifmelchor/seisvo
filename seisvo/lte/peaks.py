#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np
import pandas as pd
import scipy

from seisvo.signal.proba import get_KDE
from .plotting import plotPeaksSpecPDF, plotPeaksPDF


class Peaks(object):
    def __init__(self, lteout, chan, peaks_dict, total_peaks, fq_band):
        self.chan = chan
        self.starttime = starttime 
        self.endtime = endtime
        self.fq_band = fq_band

        self.lteout_ = lteout
        self.total_peaks_ = total_peaks
        self.dominant_peaks_ = peaks_dict

        self.peaks_ = {}
        self.nro_ = 0
        self.df_ = None


    def __str__(self):
        txt_to_return =  f'\n   >>LTE file    ::: {self.lteout.lte.file_}'
        txt_to_return += f'\n   >channel       :  {self.chan}'
        txt_to_return += f'\n   >starttime     :  {self.starttime.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >endtime       :  {self.endtime.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   >fq_range      :  {self.fq_band}'
        txt_to_return += f'\n   >Nro dom. peaks:  {self.total_peaks_[0]}'
        txt_to_return += f'\n   >Nro dom. freq :  {self.nro_}'
        txt_to_return +=  f'\n'
        return txt_to_return


    def fit(self, threshold=0.0, peak_width='auto', **kwargs):
        """
        Computes the dominant peaks following Melchor et al. 2022 (https://doi.org/10.1016/j.jsames.2022.103961)
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
            peak_width = scipy.signal.peak_widths(self.fq_pdf_, peaks)[0]
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
                rect = [] 
                th_H = [] 
                th_V = []

                for _, tdict in self.dominant_peaks_.items():
                    for n, fq in enumerate(tdict['fq']):
                        if dfq-self.peaks_[f]['width'] <= fq <= dfq+self.peaks_[f]['width']:
                            sp += [tdict['specgram'][n]]
                            pd += [tdict['degree'][n]]
                            rect += [tdict['rect'][n]]
                            th_H += [tdict['azimuth'][n]]
                            th_V += [tdict['elevation'][n]]
                
                # list sp, pd, etc. may contain None
                sp = np.array([sp_k for sp_k in sp if sp_k])
                pd = np.array([pd_k for pd_k in pd if pd_k])
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
                dout[t] = {'fq':[], 'specgram':[], 'rect':[], 'azimuth':[], 'elevation':[]}
                
                for n, fq in enumerate(tdict['fq']):
                    if fq_range[0] <= fq <= fq_range[1]:
                        dout[t]['fq'].append(fq)

                        sp = tdict['specgram'][n]
                        dout[t]['specgram'].append(sp)

                        rect = tdict['rect'][n]
                        if rect:
                            dout[t]['rect'].append(rect)
                        
                        tH = tdict['azimuth'][n]
                        if tH:
                            dout[t]['azimuth'].append(tH)
                        
                        tV = tdict['elevation'][n]
                        if tV:
                            dout[t]['elevation'].append(tV)

        return dout


    def get_fq_area(self, fq, fqtol=0.5):
        aW = np.zeros((self.lteout_.lte.stats.nro_time_bins,))
        aL = np.zeros((self.lteout_.lte.stats.nro_time_bins,))
        
        for i in range(self.lteout_.lte.stats.nro_time_bins):
            if i in self.dominant_peaks_:
                d = self.dominant_peaks_[i]
                for f, pd, r in zip(d['fq'], d['degree'], d['rect']):
                    if (fq>f-fqtol) and (fq<f+fqtol):
                        if pd and pd > self.peak_thresholds['degree_th']:
                            aW[i] = 1
                            if r and r > self.peak_thresholds['rect_th']:
                                aL[i] = 1
        return aW, aL


    def plot_spec_pdf(self, show=True, **kwargs):
        if self.nro_ == 0:
            raise ValueError ('no dominant frequencies to plot!')

        fig = plotPeaksSpecPDF(self, show=show, **kwargs)

        return fig


    def plot_peak_pdf(self, n, plot=True, **kwargs):
        if self.nro_ == 0:
            raise ValueError ('no dominant frequencies to plot!')
        
        fig = plotPeaksPDF(self, n, plot=plot, **kwargs)
        
        return fig

