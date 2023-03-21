#!/usr/bin/env python3
# coding=utf-8

import os
import scipy
import pickle
import numpy as np
import pandas as pd
import datetime as dt
from itertools import chain
from seisvo.signal.proba import get_KDE
from .plotting import plotPeaksSpecPDF, plotPeaksPDF, plotPeakTimeEvo


class Peaks(object):
    def __init__(self, peaks_dict, lte_file, starttime, endtime, delta, fq_band, peaks_threshold, **kwargs):
        self._dout     = peaks_dict
        self.starttime = starttime 
        self.endtime   = endtime
        self.delta     = delta
        self.fq_band   = fq_band
        self.ltefile_  = lte_file
        self.model_    = peaks_threshold
        self.chan      = list(peaks_dict.keys())
        self.df_       = None
        self.file_     = kwargs.get("pksfile", None)
    

    def write(self, pks_file=None):
        if not pks_file:
            pks_file = os.path.basename(self.ltefile_)
            pks_file = pks_file.replace(".lte", ".pks")

        # write pickle object
        with open(pks_file, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        print(f"  >>>  {pks_file}  created!")


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


    def fit(self, starttime=None, endtime=None, chan=None, ns_min=100, threshold=0.0, peak_width='auto', to_json=True, dpks=None, **kwargs):
        """
        Computes the dominant peaks following Melchor et al. 2022 (https://doi.org/10.1016/j.jsames.2022.103961)
        peak_width can be 'auto' or float. If 'auto', for each peak, the width is computing by signal.peak_widths.
        """

        assert 0 <= threshold <= 1

        if not isinstance(peak_width, float) and not peak_width=='auto':
            raise ValueError (" peak_width must be float or 'auto'")

        if not dpks:
            dpks = self.get(starttime=starttime, endtime=endtime, chan=chan, return_stats=False)

        if not dpks:
            print("warn :: no dominant peaks found")
            return
        
        true_chan = list(dpks.keys())

        # join all frequencies
        all_fq_pks = [nd for chan in true_chan for nd in dpks[chan]["pks"]["fq"]]
        all_fq_pks = np.array(list(chain(*all_fq_pks)))
        
        # compute PDF
        fq_space = np.linspace(all_fq_pks.min(), all_fq_pks.max(), 1000).reshape(-1, 1)
        fq_bwd  = kwargs.get("fq_bandwidth", 0.01)
        pks_kde = get_KDE(all_fq_pks.reshape(-1, 1), fq_bwd)
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
        
        # characterize dominant frequencies
        if nro_dfq > 0:
            # sxx
            all_sxx_pks = [nd for chan in true_chan for nd in dpks[chan]["pks"]["sxx"]]
            all_sxx_pks = np.array(list(chain(*all_sxx_pks)))
            # rect
            all_rect_pks = [nd for chan in true_chan for nd in dpks[chan]["pks"]["rect"]]
            all_rect_pks = np.array(list(chain(*all_rect_pks)))
            # azimuth
            all_azim_pks = [nd for chan in true_chan for nd in dpks[chan]["pks"]["azim"]]
            all_azim_pks = np.array(list(chain(*all_azim_pks)))
            # elevation
            all_elev_pks = [nd for chan in true_chan for nd in dpks[chan]["pks"]["elev"]]
            all_elev_pks = np.array(list(chain(*all_elev_pks)))

            for f in range(1, nro_dfq + 1):
                dfq = d_peaks[f]['fq'] # dominant frequency
                fhigh = dfq + d_peaks[f]['width']
                flow  = dfq - d_peaks[f]['width']
                fqn = np.where((all_fq_pks>=flow)&(all_fq_pks<=fhigh))

                sxx_dfq  = all_sxx_pks[fqn]
                rect_dfq = all_rect_pks[fqn]
                azim_dfq = all_azim_pks[fqn]
                elev_dfq = all_elev_pks[fqn]

                if np.isnan(sxx_dfq).any():
                    sxx_dfq = sxx_dfq[np.isfinite(sxx_dfq)]
                
                if sxx_dfq.shape[0] >= ns_min:
                    sxx_bwd  = kwargs.get("sxx_bandwidth", 0.05)
                    sp_kde   = get_KDE(sxx_dfq.reshape(-1,1), sxx_bwd)
                    sp_space = np.linspace(sxx_dfq.min(), sxx_dfq.max(), 500).reshape(-1,1)
                    sp_dist  = np.exp(sp_kde.score_samples(sp_space))
                    sp_dist /= sp_dist.max()
                    sp_args  = np.where(sp_dist > 0.5)[0]
                    sp_range = (sp_space[sp_args[0]], sp_space[sp_args[-1]])

                    d_peaks[f]['sp'] = {
                        'val':sp_space[np.argmax(sp_dist)],
                        'range':sp_range,
                        'n':len(sxx_dfq),
                        'kde':sp_kde,
                        'min':sxx_dfq.min(),
                        'max':sxx_dfq.max()
                    }

                    if np.isnan(rect_dfq).any():
                        rect_dfq = rect_dfq[np.isfinite(rect_dfq)]
                    
                    if rect_dfq.shape[0] >= ns_min:
                        rect_bwd = kwargs.get("rect_bandwidth", 0.1)
                        rect_kde = get_KDE(rect_dfq.reshape(-1,1), rect_bwd)
                        rect_space = np.linspace(0, 1, 500).reshape(-1,1)
                        rect_dist = np.exp(rect_kde.score_samples(rect_space))
                        rect_dist /= rect_dist.max()
                        args = np.where(rect_dist > 0.5)[0]
                        rect_range = (rect_space[args[0]], rect_space[args[-1]])

                        d_peaks[f]['rect'] = {
                            'val':rect_space[np.argmax(rect_dist)],
                            'range':rect_range,
                            'n':len(rect_dfq),
                            'kde':rect_kde
                        }
                        d_peaks[f]['cp'] = len(rect_dfq)/len(sxx_dfq)

                        azim_dfq = azim_dfq[np.isfinite(azim_dfq)]
                        elev_dfq = elev_dfq[np.isfinite(elev_dfq)]
                        
                        if azim_dfq.shape[0] >= ns_min and elev_dfq.shape[0] > ns_min:
                            d_peaks[f]['cl'] = (len(azim_dfq) + len(elev_dfq)) / (2*len(rect_dfq))

                            azim_bwd  = kwargs.get("azim_bandwidth", 5)
                            thH_kde   = get_KDE(azim_dfq.reshape(-1,1), azim_bwd)
                            thH_space = np.linspace(0, 180, 500).reshape(-1,1)
                            thH_dist  = np.exp(thH_kde.score_samples(thH_space))
                            thH_dist /= thH_dist.max()
                            args      = np.where(thH_dist > 0.5)[0]

                            d_peaks[f]['thH'] = {
                                'val':thH_space[np.argmax(thH_dist)],
                                'range':(thH_space[args[0]], thH_space[args[-1]]),
                                'n':len(azim_dfq),
                                'kde':thH_kde
                            }

                            elev_bwd  = kwargs.get("elev_bandwidth", 5)
                            thV_kde   = get_KDE(elev_dfq.reshape(-1,1), elev_bwd)
                            thV_space = np.linspace(0, 90, 500).reshape(-1,1)
                            thV_dist  = np.exp(thV_kde.score_samples(thV_space))
                            thV_dist /= thV_dist.max()
                            args      = np.where(thV_dist > 0.5)[0]

                            # save elevation
                            d_peaks[f]['thV'] = {
                                'val':thV_space[np.argmax(thV_dist)],
                                'range':(thV_space[args[0]], thV_space[args[-1]]),
                                'n':len(elev_dfq),
                                'kde':thV_kde
                            }

                else:
                    del d_peaks[f]

            # convert to dataframe
            if to_json:
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

                self.df_ = pd.DataFrame(data, index=index)
                self.__to_json__(kwargs.get("fout", None))
        
        to_return.append(d_peaks)

        return to_return


    def __to_json__(self, fout):

        if not fout:
            if self.file_:
                fout = '.'.join(os.path.basename(self.file_).split('.')[:-1])
            else:
                fout = '.'.join(os.path.basename(self.ltefile_).split('.')[:-1])
        
        if fout.split('.')[-1] != "json":
            fout += ".json"

        if os.path.isfile(fout):
            os.remove(fout)

        self.df_.to_json(fout)


    def get(self, starttime=None, endtime=None, chan=None, return_stats=True):
                
        if isinstance(chan, str):
            chan = [chan]
        
        elif isinstance(chan, list):
            true_chan = [ch for ch in chan if ch in self.chan]
        
        else:
            chan = self.chan
        
        if not starttime:
            starttime = self.starttime
        else:
            assert starttime >= self.starttime
        
        if not endtime:
            endtime = self.endtime
        else:
            assert endtime <= self.endtime
        
        time = self.get_time(to_array=True)

        dout = {}
        
        for ch in chan:
            p = self._dout[ch]["pks"]
            t = list(p.keys())
            t.sort()
            tt = time[t]
            tn = np.where((starttime <= tt)&(endtime >= tt))
            t  = np.array(t)[tn]
            pks = list(map(p.get, t))
            
            if len(pks) > 0:
                pks_dict = {}
                for attr in ("fq", "sxx", "rect", "azim", "elev"):
                    pks_dict[attr] = [pk[attr] for pk in pks]

                dout[ch] = {"time":tt[tn], "pks":pks_dict}
        
        # get stats
        if dout and return_stats:
            for ch in chan:
                dout[ch]["stats"] = {}
                for attr in ("fq", "sxx", "rect", "azim", "elev"):
                    ts = np.array(list(chain(*dout[ch]["pks"][attr])))
                    ts = ts[np.isfinite(ts)]
                    ts_kde = scipy.stats.gaussian_kde(ts)
                    dout[ch]["stats"][attr] = (ts.min(), ts.max(), ts_kde)

        return dout


    def get_time(self, starttime=None, endtime=None, to_array=False):

        if not starttime:
            starttime = self.starttime
        
        if not endtime:
            endtime = self.endtime
        
        start_diff = (starttime - self.starttime).total_seconds() # sp
        n0 =  int(np.floor(start_diff/(self.delta*60)))
    
        end_diff = (endtime - self.starttime).total_seconds() # sp
        nf =  int(np.ceil(end_diff/(self.delta*60))) + 1

        ts = [starttime + dt.timedelta(minutes=float(self.delta/2)) + dt.timedelta(minutes=float(k*self.delta)) for k in range(n0,nf)]

        if to_array:
            ts = np.array(ts)

        return ts



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
    def load(pks_file):
        with open(pks_file, 'rb') as handle:
            pks = pickle.load(handle)

        return Peaks(pks._dout, pks.ltefile_, pks.starttime, pks.endtime, pks.delta, pks.fq_band, pks.model_, pksfile=pks_file)


    @staticmethod
    def plot(peaks_dict, n=None, fq_out=[], plot=True):

        if n:
            assert n in list(peaks_dict.keys())
            fig = plotPeaksPDF(peaks_dict[n], plot=plot)

        elif fq_out:
            fig = plotPeaksSpecPDF(peaks_dict, fq_out, plot=plot)

        else:
            fig = None
            print("you should define 'n' or 'fq_out' to plot")
        
        return fig


    def plot_tevo(self, starttime=None, endtime=None, chan=None, n=None, plot=True, **kwargs):

        if not n:
            pkt = self.get(starttime=starttime, endtime=endtime, chan=chan)
            fig = plotPeakTimeEvo(pkt, plot=plot, **kwargs)

            return fig

