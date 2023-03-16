#!/usr/bin/env python3
# coding=utf-8

import os
import pickle
import numpy as np
import pandas as pd
import scipy

from itertools import chain
from seisvo.signal.proba import get_KDE
from .plotting import plotPeaksSpecPDF, plotPeaksPDF


class Peaks(object):
    def __init__(self, peaks_dict, lte_file, starttime, endtime, fq_band, peaks_threshold, **kwargs):
        self._dout     = peaks_dict
        self.starttime = starttime 
        self.endtime   = endtime
        self.fq_band   = fq_band
        self.ltefile_  = lte_file
        self.model_    = peaks_threshold
        self.chan      = list(peaks_dict.keys())

        self.file_  = kwargs.get("pksfile", None)
    

    def write(self, pks_file):
        # write pickle object
        with open(pks_file+".pks", 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return None


    @staticmethod
    def load(pks_file):
        with open(pks_file, 'rb') as handle:
            pks = pickle.load(handle)

        kwargs = {
            "pksfile":pks_file,
            "char_peaks":pks.peaks_,
            "nro_char_peaks":pks.nro_,
            "dataframes":pks.df_,
            }

        return Peaks(pks._dout, pks.ltefile_, pks.starttime, pks.endtime, pks.fq_band, pks.model_, **kwargs)


    def __str__(self):
        txt_to_return = f'\n  {"LTE File":>10}     :   {self.ltefile_}\n'
        # txt_to_return += f'\n  File    : {self.file_}'
        
        dprint = {
            "starttime":self.starttime.strftime("%d %B %Y %H:%M"),
            "endtime":self.endtime.strftime("%d %B %Y %H:%M"),
            "fq_band": self.fq_band,
            "nW": [i["nW"] for _,i in self._dout.items()],
            "nL": [i["nL"] for _,i in self._dout.items()]
        }

        for key, item in dprint.items():
            txt_to_return += f'\n  {key:>10}     :  {item}'

        txt_to_return += '\n\n  ----- Model param ------ \n'
        for key, item in self.model_.items():
            txt_to_return += f'\n  {key:>10}     :  {item}'

        txt_to_return +=  f'\n'
        return txt_to_return


    def fit(self, chan=None, ns_min=100, threshold=0.0, peak_width='auto', to_dataframe=True, **kwargs):
        """
        Computes the dominant peaks following Melchor et al. 2022 (https://doi.org/10.1016/j.jsames.2022.103961)
        peak_width can be 'auto' or float. If 'auto', for each peak, the width is computing by signal.peak_widths.
        """

        assert 0 <= threshold <= 1

        if not isinstance(peak_width, float) and not peak_width=='auto':
            raise ValueError (" peak_width must be float or 'auto'")

        if chan and chan not in self.chan:
            raise ValueError("chan not found")

        if chan:
            all_dpks = [nd for _, nd in self._dout[chan]["pks"].items()]
        else:
            all_dpks = [nd for chan in self.chan for _, nd in self._dout[chan]["pks"].items()]
        
        # join all frequencies
        all_fq_pks = [nd['fq'] for nd in all_dpks]
        all_pks = np.array(list(chain(*all_fq_pks)))
        
        # compute PDF
        fq_space = np.linspace(all_pks.min(), all_pks.max(), 1000).reshape(-1, 1)
        fq_bwd  = kwargs.get("fq_bandwidth", 0.01)
        pks_kde = get_KDE(all_pks.reshape(-1, 1), fq_bwd)
        fq_pdf  = np.exp(pks_kde.score_samples(fq_space))

        # search for dominant frequencies
        fq_pdf_norm = fq_pdf/fq_pdf.max()
        peaks, _ = scipy.signal.find_peaks(fq_pdf_norm, height=threshold, distance=5)
        nro_dfq  = len(peaks)

        # search for width of dominant frequencies
        if isinstance(peak_width, float):
            half_width = peak_width/2
        else:
            delta = np.diff(fq_space.reshape(-1,))[0]
            peak_width = scipy.signal.peak_widths(fq_pdf, peaks)[0]
            half_width = None
        
        to_return = [(fq_space, fq_pdf)]
        
        # define dominant frequencies
        d_peaks = {}
        for ip, p in enumerate(peaks):
            if half_width:
                pwidth = peak_width/2
            else:
                p_half_width = (delta*peak_width[ip])/2
                pwidth = max(p_half_width, 0.05)
                pwidth = min(p_half_width, 0.50)

            d_peaks[ip+1] = {
                'fq':fq_space[p],
                'width':pwidth,
                'fq_prob':fq_pdf[p],
                'cp':None,
                'cl':None,
                'sp':None,
                'rect':None,
                'thH':None,
                'thV':None
                }
        
        if nro_dfq > 0:
            for f in range(1, nro_dfq + 1):
                dfq = d_peaks[f]['fq'] # dominant frequency
                sp, rect, th_H, th_V = [], [], [], []

                for tdict in all_dpks:
                    for n, fq in enumerate(tdict['fq']):
                        if dfq - d_peaks[f]['width'] <= fq <= dfq + d_peaks[f]['width']:
                            sp   += [tdict['sxx'][n]]
                            # pd   += [tdict['dgr'][n]]
                            rect += [tdict['rect'][n]]
                            th_H += [tdict['azim'][n]]
                            th_V += [tdict['elev'][n]]
                
                sp   = np.array(sp, dtype=np.float32)
                if np.isnan(sp).any():
                    sp = sp[np.isfinite(sp)]
                
                if sp.shape[0] > ns_min:
                    sxx_bwd  = kwargs.get("sxx_bandwidth", 0.05)
                    sp_kde   = get_KDE(sp.reshape(-1,1), sxx_bwd)
                    sp_x = np.linspace(sp.min(), sp.max(), 500).reshape(-1,1)
                    sp_dist  = np.exp(sp_kde.score_samples(sp_x))
                    sp_dist /= sp_dist.max()
                    sp_args = np.where(sp_dist > 0.5)[0]
                    sp_range = (sp_x[sp_args[0]], sp_x[sp_args[-1]])

                    d_peaks[f]['sp'] = {
                        'val':sp_x[np.argmax(sp_dist)],
                        'range':sp_range,
                        'n':len(sp),
                        'kde':sp_kde
                    }

                    rect = np.array(rect, dtype=np.float32)
                    if np.isnan(rect).any():
                        rect = rect[np.isfinite(rect)]
                    
                    if rect.shape[0] > ns_min:
                        rect_bwd = kwargs.get("rect_bandwidth", 0.1)
                        rect_kde = get_KDE(rect.reshape(-1,1), rect_bwd)
                        rect_space = np.linspace(0, 1, 500).reshape(-1,1)
                        rect_dist = np.exp(rect_kde.score_samples(rect_space))
                        rect_dist /= rect_dist.max()
                        args = np.where(rect_dist > 0.5)[0]
                        rect_range = (rect_space[args[0]], rect_space[args[-1]])

                        d_peaks[f]['rect'] = {
                            'val':rect_space[np.argmax(rect_dist)],
                            'range':rect_range,
                            'n':len(rect),
                            'kde':rect_kde
                        }
                        d_peaks[f]['cp'] = len(rect)/len(sp)

                        th_H = np.array(th_H, dtype=np.float32)
                        if np.isnan(th_H).any():
                            th_H = th_H[np.isfinite(th_H)]

                        th_V = np.array(th_V, dtype=np.float32)
                        if np.isnan(th_V).any():
                            th_V = th_V[np.isfinite(th_V)]
                        
                        if th_H.shape[0] > ns_min and th_V.shape[0] > ns_min:
                            d_peaks[f]['cl'] = (len(th_H) + len(th_V)) / (2*len(rect))

                            azim_bwd = kwargs.get("azim_bandwidth", 5)
                            thH_kde = get_KDE(th_H.reshape(-1,1), azim_bwd)
                            thH_x = np.linspace(0, 180, 500).reshape(-1,1)
                            thH_dist = np.exp(thH_kde.score_samples(thH_x))
                            thH_dist /= thH_dist.max()
                            args = np.where(thH_dist > 0.5)[0]

                            d_peaks[f]['thH'] = {
                                'val':thH_x[np.argmax(thH_dist)],
                                'range':(thH_x[args[0]], thH_x[args[-1]]),
                                'n':len(th_H),
                                'kde':thH_kde
                            }

                            elev_bwd = kwargs.get("elev_bandwidth", 5)
                            thV_kde = get_KDE(th_V.reshape(-1,1), elev_bwd)
                            thV_x = np.linspace(0, 90, 500).reshape(-1,1)
                            thV_dist = np.exp(thV_kde.score_samples(thV_x))
                            thV_dist /= thV_dist.max()
                            args = np.where(thV_dist > 0.5)[0]

                            # save elevation
                            d_peaks[f]['thV'] = {
                                'val':thV_x[np.argmax(thV_dist)],
                                'range':(thV_x[args[0]], thV_x[args[-1]]),
                                'n':len(th_V),
                                'kde':thV_kde
                            }

                else:
                    del d_peaks[f]

            # convert to dataframe
            if to_dataframe:
                index = list(d_peaks.keys())
                data = {
                    'fq':[info['fq'][0] for _, info in d_peaks.items()],
                    'fq_prob':[info['fq_prob'] for _, info in d_peaks.items()],
                    'width':[info['width'] for _, info in d_peaks.items()]
                }
                cl = []
                cp = []
                for  _, info in d_peaks.items():
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
                    for  _, info in d_peaks.items():
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

                to_return.append(pd.DataFrame(data, index=index))

            else:
                to_return.append(d_peaks)
        
        else:
            to_return.append(None)

        return to_return


    def write_json(self, df, fout=None):
        assert isinstance(df, pd.DataFrame)

        if not fout:

            if self.file_:
                fout = '.'.join(os.path.basename(self.file_).split('.')[:-1])
            else:
                fout = '.'.join(os.path.basename(self.ltefile_).split('.')[:-1])
        
        if fout.split('.')[-1] != "json":
            fout += ".json"

        if os.path.isfile(fout):
            os.remove(fout)

        df.to_json(fout)


    def get_dpeaks(self, chan, fq_range):
        
        assert chan in self.chan

        dout = {}
        for t, tdict in self._dout[chan]["pks"].items():
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


    # def get_fq_area(self, fq, fqtol=0.5):
    #     aW = np.zeros((self.lteout_.lte.stats.nro_time_bins,))
    #     aL = np.zeros((self.lteout_.lte.stats.nro_time_bins,))
        
    #     for i in range(self.lteout_.lte.stats.nro_time_bins):
    #         if i in self.dominant_peaks_:
    #             d = self.dominant_peaks_[i]
    #             for f, pd, r in zip(d['fq'], d['degree'], d['rect']):
    #                 if (fq>f-fqtol) and (fq<f+fqtol):
    #                     if pd and pd > self.peak_thresholds['degree_th']:
    #                         aW[i] = 1
    #                         if r and r > self.peak_thresholds['rect_th']:
    #                             aL[i] = 1
    #     return aW, aL


    @staticmethod
    def plot(self, peaks_dict, n=None, fq_out=[], plot=True):

        if n:
            assert n in list(peaks_dict.keys())
            fig = plotPeaksPDF(peaks_dict[n], plot=plot)

        elif fq_out:
            fig = plotPeaksSpecPDF(peaks_dict, fq_out, plot=plot)

        else:
            fig = None
            print("you should define 'n' or 'fq_out' to plot")
        
        return fig

