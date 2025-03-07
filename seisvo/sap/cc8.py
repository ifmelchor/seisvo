#!/usr/bin/python3
# coding=utf-8

import os
import csv
import h5py
import pickle
import numpy as np
import pandas as pd
import datetime as dt
import time as pytime

from itertools import groupby
from operator import itemgetter

from .cc8utils import _CC8Process
from ..stats import CC8stats
from ..signal import get_Stats, get_PDF, slowness_vector
from ..plotting.array import simple_cc8_plot, simple_slowmap, _detections

# from tqdm import tqdm
ATTR_LIST = ["rms", "maac", "slow", "baz", "error", "slowmap", "bazbnd", "slowbnd"]


def check_cc8_attr(attr):
    if isinstance(attr, str):
        if attr in ATTR_LIST:
            return [attr]
        else:
            return []

    if isinstance(attr, (list, tuple)):
        return [at for at in attr if at in ATTR_LIST]
   

def _new_CC8(array, cc8file, headers, njobs):
    """
        Create new CC8 (hdf5) file
    """

    with h5py.File(cc8file, "w-") as cc8h5:
        # add header
        hdr = cc8h5.create_dataset('header',(1,))
        hdr.attrs['id']            = headers['id']
        hdr.attrs['locs']          = headers['locs']
        hdr.attrs['starttime']     = headers['starttime']
        hdr.attrs['endtime']       = headers['endtime']
        hdr.attrs['window']        = headers['window']
        hdr.attrs['overlap']       = headers['overlap']
        hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
        hdr.attrs['nro_slow_bins'] = headers['nro_slow_bins']
        hdr.attrs['last_time_bin'] = [-1]*len(headers['fq_bands'])
        hdr.attrs['sample_rate']   = headers['sample_rate']
        hdr.attrs['fq_bands']      = headers['fq_bands']
        hdr.attrs['slow_max']      = headers['slow_max']
        hdr.attrs['slow_int']      = headers['slow_int']
        hdr.attrs['cc_thres']      = headers['cc_thres']
        hdr.attrs['toff_sec']      = headers['toff_sec']
        hdr.attrs['slowmap']       = headers['slowmap']

        # print info
        print('')
        print(' CC8 file INFO')
        print(' ----------------')
        print(" hdf5_memory info: %s " % cc8file)
        print(' dataset size: ', (headers['nro_time_bins'],))
        print('')
        print('    CC8 stats   ')
        print(' ----------------')
        
        for info_key in ['id', 'locs', 'starttime', 'endtime', 'window',\
            'overlap', 'nro_time_bins', 'sample_rate', 'fq_bands',\
            'slow_max', 'slow_int', 'cc_thres']:
            if info_key == "window":
                print(f'   {info_key} [sec]:  {headers[info_key]}')
            else:
                print(f'   {info_key}:  {headers[info_key]}')        
        print(' ----------------\n')

        # add datasets
        timebins = headers['nro_time_bins']
        slowbins = headers['nro_slow_bins']

        for fq_n in range(1, len(headers['fq_bands'])+1):
            fq_n = cc8h5.create_group(str(fq_n))

            for attr in ("slow", "baz", "maac", "rms", "error"):
                fq_n.create_dataset(attr, (timebins,), chunks=True, dtype=np.float32)
            
            if headers['slowmap']:
                fq_n.create_dataset('slowmap', (timebins, slowbins, slowbins), chunks=True, dtype=np.float32)
            
            fq_n.create_dataset('slowbnd', (timebins, 2), chunks=True, dtype=np.float32)
            fq_n.create_dataset('bazbnd', (timebins, 2), chunks=True, dtype=np.float32)

        cc8h5.flush()
    
    cc8 = CC8(cc8file)
    cc8.__compute__(array, headers, njobs)

    return cc8


class CC8(object):
    def __init__(self, cc8file):

        assert os.path.isfile(cc8file)
    
        stats_dict = {"file":cc8file}
        with h5py.File(cc8file, "r") as f:
            hdr = f['header']
            stats_dict["id"]              = str(hdr.attrs['id'])
            stats_dict["locs"]            = list(hdr.attrs['locs'])
            stats_dict["starttime"]       = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S.%f')
            stats_dict["endtime"]         = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S.%f')
            stats_dict["window"]          = float(hdr.attrs['window'])
            stats_dict["overlap"]         = float(hdr.attrs['overlap'])
            stats_dict["nro_time_bins"]   = int(hdr.attrs['nro_time_bins'])
            stats_dict["last_time_bin"]   = list(hdr.attrs['last_time_bin'])
            stats_dict["sample_rate"]     = int(hdr.attrs['sample_rate'])
            stats_dict["fq_bands"]        = [(float(fqb[0]), float(fqb[1])) for fqb in hdr.attrs['fq_bands']]
            stats_dict["slow_max"]        = float(hdr.attrs['slow_max'])
            stats_dict["slow_int"]        = float(hdr.attrs['slow_int'])
            stats_dict['nro_slow_bins']   = int(hdr.attrs['nro_slow_bins'])
            stats_dict["cc_thres"]        = float(hdr.attrs['cc_thres'])
            try:
                stats_dict["slowmap"] = bool(hdr.attrs['slowmap'])
            except:
                stats_dict["slowmap"] = True

        # compute the true time series
        # i.e., the times point to the middle of the window
        self.stats  = CC8stats(stats_dict)
        advance     = dt.timedelta(seconds=float(self.stats.window*(1-self.stats.overlap)))
        half_delta  = dt.timedelta(seconds=float(self.stats.window/2))
        t0 = self.stats.starttime + half_delta
        timelist  = [t0]
        
        for n in range(self.stats.nro_time_bins-1):
            t0 += advance
            timelist.append(t0)

        self.time_ = np.array(timelist, dtype=dt.datetime)


    def __str__(self):
        return self.stats.__str__()


    def __read__(self, attr, fq_idx, nt):
        n0, nf = nt
        with h5py.File(self.stats.file, "r") as f:
            if attr == 'slowmap':
                if self.stats.slowmap:
                    ts = f.get(str(fq_idx))[attr][n0:nf,:]
                else:
                    ts = None
            else:
                ts = f.get(str(fq_idx))[attr][n0:nf]
        
        return ts


    def __compute__(self, array, headers, njobs):
        """
        Process CC8 file
        """

        # init process
        cc8p = _CC8Process(self.stats, headers)
        excluded_locs = array.get_excluded_locs(self.stats.locs)

        # define timedelta parameters
        delta = dt.timedelta(seconds=float(headers["interval"]*60))
        toff  = dt.timedelta(seconds=headers["int_extra_sec"])
        olap  = dt.timedelta(seconds=float(headers["overlap"]*headers["window"]))
        toff_sec = dt.timedelta(seconds=headers["toff_sec"])
        
        for nint in range(1, headers["nro_intervals"]+1):
            if nint >= 2:
                start_int = end_int - olap
            else:
                start_int = self.stats.starttime
            
            end_int   = start_int + delta + toff
            
            if nint == headers["nro_intervals"]:
                last = True
                end_int = self.stats.endtime
            
            else:
                last = False

            stream = array.get_stream(start_int, end_int,\
                toff_sec=headers["toff_sec"], exclude_locs=excluded_locs)
            
            if stream:
                tsi, tei = stream.get_bounds()
                check_bounds = tsi <= start_int and tei >= end_int
                if check_bounds:
                    data = stream.to_array(detrend=True)
                else:
                    data = None
            else:
                data = None
            
            for fqn, fqband in enumerate(self.stats.fq_bands):
                cc8p.run(data, start_int, end_int, int(fqn+1), last)

                if len(cc8p.processes) == njobs:
                    cc8p.wait(nint)
        
        # check if there are any process running
        if len(cc8p.processes) > 0:
            cc8p.wait(nint)


    def get_time(self, starttime, endtime):
        real_time = np.array(self.time_)
        n0 = np.argmin(np.abs(real_time - starttime))
        nf = np.argmin(np.abs(real_time - endtime)) + 1
        return real_time[n0:nf], (n0, nf)


    def get(self, starttime=None, endtime=None, interval=None, attr=None, slowmap=True, fq_idx=None):
        """
        Return a CC8out object
        interval in minutes result in endtime = starttime + interval
        """

        if not starttime or starttime < self.time_[0]:
            starttime = self.time_[0]

        if interval:
            endtime = starttime + dt.timedelta(minutes=interval)
        else:
            if not endtime or endtime > self.time_[-1]:
                endtime = self.time_[-1]

        if not fq_idx:
            fq_idx = self.stats.fqidx
        else:
            fq_idx = self.stats.check_idx(fq_idx)[0]
        
        if not attr:
            attr = []
            for a in ATTR_LIST:
                if a == "slowmap":
                    if slowmap and self.stats.slowmap:
                        attr.append(a)
                else:
                    attr.append(a)
        else:
            attr = check_cc8_attr(attr)
        
        dout = {}
        dout["dtime"], (n0,nf) = self.get_time(starttime, endtime)
        dout["fqidx"]    = []
        dout["attlist"]  = attr
        
        for nfq in fq_idx:
            dout["fqidx"].append(nfq)

            for at in attr:
                attr_key = '/'.join([nfq, at])
                ts  = self.__read__(at, nfq, (n0,nf))
                
                if ts[ts>0].any():
                    dout[attr_key] = ts
        
        return CC8out(self.stats, dout)


    def gui(self, interval, starttime=None, db=None, **kwargs):
        """
        GUI for navigate the cc8 file.
        interval in minutes
        """

        from ..gui import load_cc8widget

        if not starttime:
            starttime = self.time_[0]

        fq_idx = kwargs.get("fq_idx", self.stats.fqidx[0])
        widget = load_cc8widget(self, starttime, interval, fq_idx, db, **kwargs)

        return


    def auto_detect(self, n_nearest, n_min=1, fq_idx=None, **nidx_kwargs):
        """
        An event is defined as a successive detection (controlled by n_min).
        This algorithm search for events and write its characteristics into a csv.
        Return a dataframe object
        """

        cco  = self.get(fq_idx=fq_idx)
        maac = cco.get_data('maac')
        rms  = cco.get_data('rms')
        baz  = cco.get_data('baz')
        slow = cco.get_data('slow')
        nidx = cco.get_nidx(return_full=False, **nidx_kwargs)

        maxwt   = 0
        
        # columns
        ndList        = []
        starttimeList = []
        durationList  = []
        maac_avgList  = []
        rms_avgList   = []
        d_maacList    = []
        d_rmsList     = []
        slow_avgList  = []
        slow_stdList  = []
        baz_avgList   = []
        baz_stdList   = []

        for k, g in groupby(enumerate(nidx), lambda ix: ix[0]-ix[1]):
            ans = list(map(itemgetter(1), g))
            nd  = len(ans)

            starttime = self.time_[ans[0]]
            duration  = nd*self.stats.window - (nd-1)*(self.stats.window*self.stats.overlap)

            if nd >= n_min:
                if nd > 1:
                    maac_avg = maac[ans].mean()
                    rms_avg  = rms[ans].mean()
                    baz_avg  = baz[ans].mean()
                    baz_std  = baz[ans].std()
                    slow_avg = slow[ans].mean()
                    slow_std = slow[ans].std()
                    bi = ans[0]
                    bf = ans[-1]
                else:
                    maac_avg = maac[bi]
                    rms_avg  = rms[bi]
                    baz_avg  = baz[bi]
                    slow_avg = slow[bi]
                    bi = bf = ans[0]
                    baz_std = slow_std = np.nan
            else:
                continue

            d_maac = 0
            d_rms  = 0
            n = 0
            for i in range(1, n_nearest):
                if bi-i >= 0:
                    d_maac += maac[bi-i] - maac_avg
                    d_rms  += rms[bi-i] - rms_avg
                    n += 2
                if bf+i < len(maac):
                    d_maac += maac[bf+i] - maac_avg
                    d_rms  += rms[bf+i] - rms_avg
                    n += 2
            d_maac = abs(d_maac) / n
            d_rms  = abs(d_rms) / n

            # save into list
            ndList.append(nd)
            starttimeList.append(starttime)
            durationList.append(duration)
            maac_avgList.append(maac_avg)
            rms_avgList.append(rms_avg)
            d_maacList.append(d_maac)
            d_rmsList.append(d_rms)
            slow_avgList.append(slow_avg)
            slow_stdList.append(slow_std)
            baz_avgList.append(baz_avg)
            baz_stdList.append(baz_std)

        ndList        = pd.Series(ndList)
        starttimeList = pd.Series(starttimeList)
        durationList  = pd.Series(durationList)
        maac_avgList  = pd.Series(maac_avgList)
        rms_avgList   = pd.Series(rms_avgList)
        d_maacList    = pd.Series(d_maacList)
        d_rmsList     = pd.Series(d_rmsList)
        slow_avgList  = pd.Series(slow_avgList)
        slow_stdList  = pd.Series(slow_stdList)
        baz_avgList   = pd.Series(baz_avgList)
        baz_stdList   = pd.Series(baz_stdList)

        sers = [starttimeList, ndList, durationList, maac_avgList, rms_avgList, d_maacList, d_rmsList, slow_avgList, slow_stdList, baz_avgList, baz_stdList]
        df   = pd.concat(sers, axis=1)
        df.columns = ['time', 'n', 'duration', 'maac', 'rms', 'd_maac', 'd_rms', 'slow', 'slow_u', 'baz', 'baz_u']
        return df


class CC8out(object):
    def __init__(self, cc8stats, dout):
        self.cc8stats = cc8stats
        self.starttime = dout["dtime"][0]
        self.endtime   = dout["dtime"][-1]
        self._dout = dout

        # set stats for checks
        self.attr_list = dout["attlist"]
        self._fqidx = dout["fqidx"]
        self._fqidx.sort()


    def __str__(self):
        txt_to_return =  f'\n    LTE file     ::: {self.cc8stats.file}'
        txt_to_return += f'\n    ID            :  {self.cc8stats.id}'
        txt_to_return += f'\n    starttime     :  {self.starttime.strftime("%d %B %Y %H:%M:%S.%f")}'
        txt_to_return += f'\n    endtime       :  {self.endtime.strftime("%d %B %Y %H:%M:%S.%f")}'
        txt_to_return += f'\n    attributes    :  {self.attr_list}'
        txt_to_return += f'\n    freq idx      :  {self._fqidx}'
        txt_to_return +=  f'\n'
        return txt_to_return

    @property
    def time(self):
        return self._dout["dtime"]


    def check_attr(self, attr):

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

        return attr_list


    def get_nidx(self, fq_idx=None, return_full=True, **kwargs):
        """
        Return time bins only all attributes are available
        """

        # get kwargs
        maac_th  = kwargs.get("maac_th", 0)
        maac_rv  = kwargs.get("maac_rv", False)
        max_err  = kwargs.get("max_err", 0.0)
        rms_th   = kwargs.get("rms_th", 0)
        rms_rv   = kwargs.get("rms_rv", False)
        baz_int  = kwargs.get("baz_int", [])
        rms_int  = kwargs.get("rms_int", [])

        if "maac" not in self.attr_list:
            return None

        if not fq_idx:
            fq_idx   = self._fqidx[0]

        fq_idx = self.cc8stats.check_idx(fq_idx)[0]
        error = self._dout["/".join([fq_idx, "error"])]
        maac  = self._dout["/".join([fq_idx, "maac"])]
        baz   = self._dout["/".join([fq_idx, "baz"])]
        slow   = self._dout["/".join([fq_idx, "slow"])]
        rms   = self._dout["/".join([fq_idx, "rms"])]

        # init nidx
        nidx = np.arange(len(maac), dtype=np.int32)

        # slow max range
        nidx = np.where(slow<self.cc8stats.slow_max, nidx, np.nan)

        # apply maac threshold
        if maac_rv:
            nidx = np.where(maac<maac_th, nidx, np.nan)
        else:
            nidx = np.where(maac>maac_th, nidx, np.nan)

        # apply azimuth interval
        if baz_int:
            bazmin, bazmax = baz_int
            nidx = np.where(((baz<bazmax) & (baz>bazmin)), nidx, np.nan)

        # apply rms threshold
        nidx = np.where(rms>0, nidx, np.nan)
        rms  = 10 * np.log10(rms)
        if rms_int:
            rmsmin, rmsmax = rms_int
            nidx = np.where(((rms<=rmsmax) & (rms>=rmsmin)), nidx, np.nan)
        else:
            if rms_rv:
                nidx = np.where(rms<rms_th, nidx, np.nan)
            else:
                nidx = np.where(rms>rms_th, nidx, np.nan)

        # apply max_error
        if max_err > 0:
            nidx = np.where(error<max_err, nidx, np.nan)
            
        if return_full:
            return nidx
        else:
            return nidx[~np.isnan(nidx)].astype(dtype=np.int32)


    def get_data(self, attr, fq_idx=None, rms_in_db=True, **nidx_kwargs):
        """
        Return the time series of specific attribute.
        attr can be list or string, if is a list it returns a dict
        """

        if not fq_idx:
            fq_idx   = self._fqidx[0]

        attr_list = self.check_attr(attr)
        fq_idx    = self.cc8stats.check_idx(fq_idx)[0]
        
        # always return full array
        nidx_kwargs["return_full"] = True
        nidx  = self.get_nidx(fq_idx=fq_idx, **nidx_kwargs)
        nidx  = np.isnan(nidx)

        data_dict = {}
        for attr in attr_list:
            # get time series
            key = "/".join([fq_idx, attr])
            data = np.copy(self._dout[key])

            if attr == "rms" and rms_in_db:
                data = 10*np.log10(data)

            if attr not in ("slowmap", "slowbnd", "bazbnd"):
                data = np.where(nidx, np.nan, data)
                data_dict[attr] = data
            
            else:
                if attr == "slowmap":
                    print(" slowmap not implemented yet")
                    continue

                for n, idx_nan in enumerate(nidx):
                    if idx_nan:
                        data[n,:] = [np.nan, np.nan]

                data_dict[attr] = data

        if len(attr_list) > 1:
            return data_dict
        
        else:
            return data_dict[attr_list[0]]


    def get_stats(self, attr, fq_idx=None, rms_in_db=True, **nidx_kwargs):
        """
        Return (min, max, mean, mode) of an specific scalar-only attribute
        """
        
        attr  = self.check_attr(attr)[0]
        data  = self.get_data(attr, fq_idx=fq_idx, rms_in_db=rms_in_db, **nidx_kwargs)
        stats = get_Stats(data[np.isfinite(data)])

        return stats


    def get_pdf(self, attr, fq_idx=None, rms_in_db=True, vmin=None, vmax=None, data=None, n_min=10, **nidx_kwargs):
        """
        Return the PDF of a scalar attribute
        """

        if not isinstance(data, np.ndarray):
            attr = self.check_attr(attr)[0]
            data = self.get_data(attr, fq_idx=fq_idx, rms_in_db=rms_in_db, **nidx_kwargs)
        
        data = data[np.isfinite(data)]
        data = data.reshape(-1,1)
        
        if vmin == None:
            vmin = data.min()

        if vmax == None:
            vmax = data.max()

        space = np.linspace(vmin, vmax, 1000)

        if data.shape[0] > n_min:
            pdf = get_PDF(data, space)
        else:
            pdf = np.full((space.shape[0], data.shape[1]), np.nan)

        return pdf, space


    def plot(self, fq_idx=None, show_title=True, maac_th=0.75, max_err=1.0, rms_th=0.0, rms_lim=[], **fig_kwargs):

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx

        fq_band = self.cc8stats.fq_bands[int(fq_idx)-1]

        attr_list = ["rms", "maac", "slow", "baz", "slowbnd", "bazbnd"]
        datattr   = self.get_data(attr_list, fq_idx=fq_idx, max_err=max_err, maac_th=maac_th, rms_th=rms_th)
        
        datapdf = {
            "slow":self.get_pdf("slow", vmin=0, vmax=slomax, data=datattr["slow"]),
            "baz":self.get_pdf("baz", vmin=0, vmax=360, data=datattr["baz"])
        }

        fig_kwargs["rms_rv"]  = [
            self.get_data("rms", fq_idx=fq_idx, max_err=max_err, maac_th=maac_th, maac_rv=True),
            self.get_data("maac", fq_idx=fq_idx, max_err=max_err, maac_th=maac_th, rms_th=rms_th, rms_rv=True)
        ]
        
        fig_kwargs["maac_rv"] = [
            self.get_data("maac", fq_idx=fq_idx, max_err=max_err, maac_th=maac_th, maac_rv=True),
            self.get_data("maac", fq_idx=fq_idx, max_err=max_err, maac_th=maac_th, rms_th=rms_th, rms_rv=True)
        ]
        
        if show_title:
            title = f"{self.cc8stats.id} \n Fq {fq_band} :: Slomax/Sloint [{slomax}/{sloint}] \n {self._dout['dtime'][0]}"
            fig_kwargs["title"] = title
        
        if not rms_lim:
            fig_kwargs["rms_min"] = np.nanmin(datattr["rms"])
            fig_kwargs["rms_max"] = np.nanmax(datattr["rms"])
        else:
            fig_kwargs["rms_min"] = rms_lim[0]
            fig_kwargs["rms_max"] = rms_lim[1]


        ans = simple_cc8_plot(self._dout["dtime"], datattr, datapdf, **fig_kwargs)

        return ans
    

    def plot_beamform(self, ntime=None, off_sec=0, fq_idx=None, show_title=True, fq_band=[], **fig_kwargs):
        """
        Plot shifted traces
        """
        
        arr, exclude_locs = self.cc8stats.get_array()

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx
        
        if not fq_band:
            fq_band = self.cc8stats.fq_bands[int(fq_idx)-1]
        
        slowarg = {
            "slomax":self.cc8stats.slow_max,
            "sloint":self.cc8stats.slow_int,
            "exclude_locs":exclude_locs,
            "fq_band":fq_band,
        }

        data_dict = self.get_data(["maac", "slow", "baz", "error"], fq_idx=fq_idx)
        maac  = data_dict["maac"]
        slow  = data_dict["slow"]
        baz   = data_dict["baz"]
        error = data_dict["error"]
        
        if not ntime:
            ntime = np.nanargmax(maac)
            print(f" Best MAAC [{maac[ntime]:.1f}] is at position {ntime} with slow {slow[ntime]:.2f} and baz {baz[ntime]:.1f} and ERROR {error[ntime]:.2f}")

        bazt     = baz[ntime]
        slowt    = slow[ntime]
        half_w   = dt.timedelta(seconds=float(self.cc8stats.window/2))
        dtoff    = dt.timedelta(seconds=off_sec)
        start    = self._dout["dtime"][ntime] - half_w - dtoff
        end      = self._dout["dtime"][ntime] + half_w + dtoff
        duration = (end - start).total_seconds()
        startw   = duration/2 - self.cc8stats.window/2
        endw     = startw + self.cc8stats.window

        fig      = arr.beamform(start, end, slowt, bazt, slowarg=slowarg,
            shadow_times=(startw, endw), taper=False, plot=True, **fig_kwargs)

        return fig


    def plot_smap(self, ntime=None, fq_idx=None, show_title=True, **fig_kwargs):

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx

        fq_band = self.cc8stats.fq_bands[int(fq_idx)-1]

        data_dict  = self.get_data(["maac", "slow", "baz", "slowmap"], fq_idx=fq_idx)
        # data_dict = self.get_data(["maac", "slow", "baz"], fq_idx=fq_idx)
        data  = data_dict["slowmap"]
        maac  = data_dict["maac"]
        slow  = data_dict["slow"]
        baz   = data_dict["baz"]

        if not ntime:
            ntime = np.argmax(maac)
            print(f" Best MAAC [{maac[ntime]:.1f}] is at position {ntime} with slow {slow[ntime]:.2f} and baz {baz[ntime]:.1f}")

        bazt    = baz[ntime]
        slowt   = slow[ntime]
        maact   = maac[ntime]
        timet   = self._dout["dtime"][ntime]

        # make fig
        if show_title:
            title = f" nidx :: {ntime}  >> {timet} \n Fq {fq_band} :: Slomax/Sloint [{slomax}/{sloint}] \n MAAC {maact:.1f} :: Slow {slowt:.2f} [s/km] :: Baz {bazt:.1f}"
            fig_kwargs["title"] = title

        fig = simple_slowmap(data[ntime,:,:], sloint, slomax, **fig_kwargs)

        return fig
    

    def prob_slowmap(self, fq_idx=None, nidx=None, plot=True, **fig_kwargs):

        if not fq_idx:
            fq_idx   = self._fqidx[0]
        
        fq_band   = self.cc8stats.fq_bands[int(fq_idx)-1]
        slowmap  = self._dout["/".join([fq_idx, "slowmap"])]

        # get filtered slowmap
        if isinstance(nidx, np.ndarray):
            slowmap = slowmap[np.isfinite(nidx)]

        # do pdf map
        nites = self.cc8stats.nro_slow_bins
        pdfmap = np.zeros((nites,nites))
        for ii in range(nites):
            for jj in range(nites):
                data = slowmap[:,ii,jj]
                data = data[np.isfinite(data)]
                if data.shape[0] > 1:
                    pdfmap[ii,jj] = get_Stats(data)[3]

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        if plot:
            if not fig_kwargs:
                fig_kwargs = {}

            fig_kwargs["cmap"] = fig_kwargs.get("cmap", "gist_earth_r")
            fig = simple_slowmap(pdfmap, sloint, slomax, cc_th=0, **fig_kwargs)

        return pdfmap, (sloint, slomax)
    

    def compute_smap(self, ntime, fq_idx=None, slowarg={}, show_title=True, **fig_kwargs):
        """
        This function computes the smap for a specific window (ntime)
        """

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx

        arr, exloc = self.cc8stats.get_array()

        if not slowarg:
            slowarg = {
                "slomax":self.cc8stats.slow_max,
                "sloint":self.cc8stats.slow_int,
                "fq_band":self.cc8stats.fq_bands[int(fq_idx)-1],
                "exclude_locs":exloc,
                "cc_thres":self.cc8stats.cc_thres,
                "slow0":[0,0]
            }
        
        timet  = self._dout["dtime"][ntime]
        timet -= dt.timedelta(seconds=self.cc8stats.window/2)

        fig = arr.slowmap(timet, self.cc8stats.window, slowarg=slowarg, plot=True, show_title=show_title, **fig_kwargs) 
    
        return fig


    def slowness_vector(self, slow, baz, slow0=[0.,0.]):
        """
        Give information about slowness vector. 
        return_dtimes :: also give delta times for each stations of the array
        """

        pxy, pij = slowness_vector(slow, baz, self.cc8stats.slow_max, self.cc8stats.slow_int, slow0=slow0)

        return (pxy, pij)


    # def recompute_cc8(self, slomax, sloint, fq_idx=None, **nidx_kwargs):
        
    #     # este codigo tiene que seleccionar los periodos de tiempo en donde se cumplan las condiciones de nidx_kwargs, y en cada uno de los tiempos, recalcular el vector lentitud aparente utilizando un slomax y sloint adecuado.

    #     if not fq_idx:
    #         fq_idx = self._fqidx[0]
    #     else:
    #         assert fq_idx in self._fqidx

    #     arr, exloc = self.cc8stats.get_array()

    #     data = self.get_data(["slow", "baz"], fq_idx=fq_idx, **nidx_kwargs)
        
    #     slowarg0 = {
    #         "slomax":slomax,
    #         "sloint":sloint,
    #         "fq_band":self.cc8stats.fq_bands[int(fq_idx)-1],
    #         "exclude_locs":exloc,
    #         "cc_thres":self.cc8stats.cc_thres
    #     }
#      for each nidx to recompute
#      search slow0
#      fill slowarg --> slowarg["slow0"] = [0,0]



class CC8detect(object):
    def __init__(self, cc8, pkfile=None):
        self.cc8 = cc8
        if pkfile:
            self._load(pkfile)
    
    def _fit(self, pkfile='df.pkl', **fit_kwargs):
        dout = {}
        dout["n_nearest"] = fit_kwargs.get("n_nearest", 5)
        dout["n_min"]     = fit_kwargs.get("n_min", 3)
        dout["fq_idx"]    = fit_kwargs.get("fq_idx", 1)
        dout["maac_th"]   = fit_kwargs.get("maac_th", 0.7)
        dout["max_err"]   = fit_kwargs.get("max_err", 0.02)
        dout["df"] = self.cc8.auto_detect(
            dout["n_nearest"], n_min=dout["n_min"], 
            fq_idx=dout["fq_idx"], maac_th=dout["maac_th"], 
            max_err=dout["max_err"])

        # write pickle
        with open(pkfile, 'wb') as handle:
            pickle.dump(dout, handle, protocol=pickle.HIGHEST_PROTOCOL)

        pytime.sleep(1)
        self._load(pkfile)
    

    def _load(self, pkfile):
        with open(pkfile, 'rb') as f:
            data = pickle.load(f)

        self.n_nearest = data["n_nearest"]
        self.n_min     = data["n_min"]
        self.fq_idx    = data["fq_idx"]
        self.maac_th   = data["maac_th"]
        self.max_err   = data["max_err"]
        self.df_       = data["df"]
        self.times_    = self.df_.time.dt.to_pydatetime()
    

    def filter(self, df=None, **kwargs):
        d_maac_min   = kwargs.get("d_maac_min", None)
        d_maac_max   = kwargs.get("d_maac_max", None)
        d_rms_min    = kwargs.get("d_rms_min", None)
        d_rms_max    = kwargs.get("d_rms_max", None)
        baz_std_max  = kwargs.get("baz_std_max", None)
        baz_std_min  = kwargs.get("baz_std_min", None)
        slow_std_max = kwargs.get("slow_std_max", None)
        slow_std_min = kwargs.get("slow_std_min", None)
        time_range   = kwargs.get("time_range", [])
        dur_range    = kwargs.get("dur_range", [])
        baz_range    = kwargs.get("baz_range", [])
        slo_range    = kwargs.get("slow_range", [])

        if not isinstance(df, pd.DataFrame):
            df = self.df_

        cond = df.maac > 0.0

        if d_maac_min:
            dm = df.d_maac
            cond &= dm > d_maac_min

        if d_maac_max:
            dm = df.d_maac
            cond &= dm < d_maac_max

        if d_rms_min:
            dr = df.d_rms
            cond &= dr > d_rms_min

        if d_rms_max:
            dr = df.d_rms
            cond &= dr < d_rms_max

        if baz_std_max:
            baz_std = df.baz_u
            cond &= baz_std < baz_std_max

        if baz_std_min:        
            baz_std = df.baz_u
            cond &= baz_std > baz_std_min

        if slow_std_max:
            slow_std = df.slow_u
            cond &= slow_std < slow_std_max

        if slow_std_min:
            slow_std = df.slow_u
            cond &= slow_std > slow_std_min

        if time_range:
            ti, tf = time_range
            cond &= pd.Series(np.where((self.times_>=ti) & (self.times_<=tf), True, False))

        if dur_range:
            d1, d2 = dur_range
            dur = df.duration
            cond &= (dur >= d1) & (dur <= d2)

        if baz_range:
            baz1, baz2 = baz_range
            baz = df.baz
            cond &= (baz >= baz1) & (baz <= baz2)

        if slo_range:
            slo1, slo2 = slo_range
            slo = df.slow
            cond &= (slo >= slo1) & (slo <= slo2)

        return df[cond]


    def beamform(self, n, **bm_kwargs):

        row  = self.df_.iloc[n]
        time = self.times_[n]

        ans = self._beamform(self.cc8, time, row.duration, self.fq_idx, **bm_kwargs)
        
        return ans
    

    def to_bmm(self, off_sec=1, taper=False, fout='events-1.pkl', **filt_kwargs):
        
        """
        Save data into pickle file
        """

        df   = self.filter(**filt_kwargs)
        idx  = df.index.to_list()
        nidx = len(idx)

        max_duration = int(np.ceil(df.duration.max())) + 1

        if isinstance(off_sec, (list,tuple)):
            osL = off_sec[0]
            osR = off_sec[1]
            ofs = osL + osR
        else:
            ofs = 2*off_sec
        
        size = int(max_duration + ofs) * self.cc8.stats.sample_rate
        bmm_dict = {
            "time" : df.time.dt.to_pydatetime(),
            "dur"  : df['duration'].to_numpy(),
            "eid"  : np.array(idx),
            "bmm"  : np.zeros((nidx, size)),
            "baz"  : np.zeros(nidx),
            "bazu" : np.zeros((nidx, 2)),
            "slow" : np.zeros(nidx),
            "slowu": np.zeros((nidx, 2)),
            "maac" : np.zeros(nidx),
            "error": np.zeros(nidx),
            "rms"  : np.zeros(nidx)
        }

        for k, n in enumerate(idx):
            full, bm = self.beamform(n, off_sec=off_sec, sloint=0.01, taper=taper, return_full=True)

            bmm_dict["baz"][k]   = full["baz"][0] 
            bmm_dict["slow"][k]  = full["slow"][0]
            bmm_dict["maac"][k]  = full["maac"][0]
            bmm_dict["error"][k] = full["error"][0]
            bmm_dict["rms"][k]   = full["rms"][0]
            bmm_dict["bazu"][k,:]  = full["bazbnd"][0]
            bmm_dict["slowu"][k,:] = full["slowbnd"][0]

            off = (size - len(bm)) / 2

            if off.is_integer():
                off = int(off)
                bmm_dict["bmm"][k,off:-off] = bm
            else:
                off = int(off)
                bmm_dict["bmm"][k,off+1:-off] = bm

            print(f" .... {k} [{int(100 * (k/nidx))}%]", end="\r")

        with open(fout, 'wb') as handle:
            pickle.dump(bmm_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return
    

    def get_ratio(self, interval_in_hr=1, return_best=False, **filt_kwargs):
        df    = self.filter(**filt_kwargs)
        dft   = df.time.dt.to_pydatetime() # convert to datetime
        
        # create interval range
        start = self.cc8.stats.starttime
        end   = self.cc8.stats.endtime
        delta = dt.timedelta(hours=interval_in_hr)
        time  = np.arange(start, end, delta, dtype=dt.datetime)

        if return_best:
            best = {}
            n = 1
        else:
            count = []

        for t1 in time:
            t2  = t1 + delta
            idx = np.where((dft>=t1) & (dft<=t2))[0]
            nt = len(idx)

            if return_best:
                if nt >= 1:
                    i0 = min(idx)
                    i1 = i0 + nt
                    sdf = df.iloc[i0:i1]

                    # get clear+stable
                    sdf = self.filter(df=sdf, d_maac_min=0.1, d_rms_min=1.2, baz_std_max=2.0, slow_std_max=0.05)
                    
                    # get best d_rms
                    sdf = sdf[sdf["d_rms"]==sdf["d_rms"].max()]

                    if len(sdf) >= 1:
                        best[n] = {
                            "time":sdf.time.dt.to_pydatetime()[0],
                            "duration":sdf.duration.to_numpy()[0],
                            "n":sdf.n.to_numpy()[0],
                            "baz":sdf.baz.to_numpy()[0],
                            "baz_u":sdf.baz_u.to_numpy()[0],
                            "slow":sdf.slow.to_numpy()[0],
                            "slow_u":sdf.slow_u.to_numpy()[0]

                        }

                        n += 1
            else:
                count.append(nt)

        if return_best:
            return best

        else:
            count = np.array(count, dtype=np.float32)
            count[count==0] = np.nan

            return time, count
    

    @staticmethod
    def _beamform(cc8, time, duration, fq_idx, off_sec=0.75, sloint=0.01, taper=False, plot=True, return_full=True, return_st=False, fq_band=None):

        ar, exclude_locs = cc8.stats.get_array()
        winTime   = dt.timedelta(seconds=cc8.stats.window/2)
        time -= winTime

        if isinstance(off_sec, (list, tuple)):
            offDeltaL = dt.timedelta(seconds=off_sec[0])
            nptL = int(off_sec[0]*100)
            offDeltaR = dt.timedelta(seconds=off_sec[1])
            nptR = int(off_sec[1]*100)
        else:
            offDeltaL = offDeltaR = dt.timedelta(seconds=off_sec)
            nptL = nptR = int(off_sec*100)
        
        starttime = time - offDeltaL
        endtime   = time + dt.timedelta(seconds=duration) + offDeltaR

        if not fq_band:
            fq_band = cc8.stats.fq_bands[fq_idx-1]

        slowarg = {
                "slomax":cc8.stats.slow_max, 
                "sloint":sloint,
                "fq_band":fq_band,
                "exclude_locs":exclude_locs
            }

        full, bm = ar.beamform(starttime, endtime, 0, 0, slowarg=slowarg, taper=taper, plot=plot, return_full=True)

        # compute the SNR
        signal = np.abs(bm[nptL:-nptR]).mean()
        noiseL = np.abs(bm[:nptL]).mean()
        noiseR = np.abs(bm[-nptR:]).mean()
        snr = 2*signal/(noiseR+noiseL)
        full["snr"] = snr

        return full, bm