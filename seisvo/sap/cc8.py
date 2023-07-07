#!/usr/bin/python3
# coding=utf-8

import os
import h5py
import numpy as np
import pandas as pd
import datetime as dt

from .cc8utils import _CC8Process
from ..stats import CC8stats
from ..signal import get_Stats, get_PDF
from ..plotting.array import simple_cc8_plot, simple_slowmap, window_wvfm
from ..gui import load_cc8widget

# from tqdm import tqdm
ATTR_LIST = ["rms", "maac", "slow", "bazm", "slowmap", "bazmbnd", "slowbnd"]


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

            for attr in ("slow", "bazm", "maac", "rms"):
                fq_n.create_dataset(attr, (timebins,), chunks=True, dtype=np.float32)
            
            fq_n.create_dataset('slowmap', (timebins, slowbins, slowbins), chunks=True, dtype=np.float32)
            fq_n.create_dataset('slowbnd', (timebins, 2), chunks=True, dtype=np.float32)
            fq_n.create_dataset('bazmbnd', (timebins, 2), chunks=True, dtype=np.float32)

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

        self.stats  = CC8stats(stats_dict)

        # compute the true time series
        start       = self.stats.starttime
        advance     = dt.timedelta(seconds=float(self.stats.window*(1-self.stats.overlap)))
        half_delta  = dt.timedelta(seconds=float(self.stats.window/2))
        self.time_ = [start+half_delta]

        for n in range(self.stats.nro_time_bins-1):
            start += advance
            self.time_.append(start + half_delta)


    def __str__(self):
        return self.stats.__str__()


    def __read__(self, attr, fq_idx, nt):
        n0, nf = nt
        with h5py.File(self.stats.file, "r") as f:
            if attr == 'slowmap':
                ts = f.get(str(fq_idx))[attr][n0:nf,:]
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
            tsi, tei = stream.get_bounds()
            check_bounds = tsi <= start_int and tei >= end_int

            if stream and check_bounds:
                data = stream.to_array(detrend=True)
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

        n0 = np.argmin(np.abs(np.array(self.time_) - starttime))
        nf = np.argmin(np.abs(np.array(self.time_) - endtime)) + 1
        
        return self.time_[n0:nf], (n0, nf)


    def get(self, starttime=None, endtime=None, attr=None, slowmap=True, fq_idx=None):
        
        if not starttime:
            starttime = self.stats.starttime
        else:
            assert starttime >= self.stats.starttime

        if not endtime:
            endtime = self.stats.endtime
        else:
            assert endtime <= self.stats.endtime

        if not fq_idx:
            fq_idx = self.stats.fqidx
        else:
            fq_idx = self.stats.check_idx(fq_idx)[0]
        
        if not attr:
            attr = []
            for a in ATTR_LIST:
                if a == "slowmap" and slowmap:
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


    def gui(self, interval, starttime=None, olap=0.1, **kwargs):
        """
        GUI for navigate the cc8 file.
        interval in minutes
        """

        if not starttime:
            starttime = self.stats.starttime

        fq_idx   = kwargs.get("fq_idx", self.stats.fqidx[0])
        maac_th  = kwargs.get("maac_th", 0.6)
        max_err  = kwargs.get("max_err", 0.7)
        baz_int  = kwargs.get("baz_int", [])

        widget = load_cc8widget(self, starttime, interval, fq_idx, olap=olap, maac_th=maac_th, baz_int=baz_int, max_err=max_err)


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
        max_err  = kwargs.get("max_err", 0)
        rms_th   = kwargs.get("rms_th", 0)
        rms_rv   = kwargs.get("rms_rv", False)
        baz_int  = kwargs.get("baz_int", [])

        if "maac" not in self.attr_list:
            return None

        if not fq_idx:
            fq_idx   = self._fqidx[0]

        fq_idx   = self.cc8stats.check_idx(fq_idx)[0]
        mackey = "/".join([fq_idx, "maac"])
        rmskey = "/".join([fq_idx, "rms"])
        bazkey = "/".join([fq_idx, "bazm"])
        maac = np.copy(self._dout[mackey])
        baz  = np.copy(self._dout[bazkey])
        rms  = np.copy(self._dout[rmskey])

        # init nidx
        nidx = np.arange(len(maac))

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
        if rms_rv:
            nidx = np.where(rms<rms_th, nidx, np.nan)
        else:
            nidx = np.where(rms>rms_th, nidx, np.nan)

        # apply max_error
        if max_err > 0:
            slokey = "/".join([fq_idx, "slow"])
            slowb  = np.copy(self._dout[slokey+"bnd"])
            bazb   = (np.pi/180)*np.copy(self._dout[bazkey+"bnd"])
            slodiff  = np.abs(slowb[:,1] - slowb[:,0])
            bazdiff  = np.abs(bazb[:,1] - bazb[:,0])
            error  = (slodiff+bazdiff)/2
            nidx = np.where(error<max_err, nidx, np.nan)

        if return_full:
            return nidx
        else:
            return nidx[np.isfinite(nidx)].astype(dtype=np.int32)


    def get_data(self, attr, fq_idx=None, rms_in_db=True, **nidx_kwargs):
        """
        Return the time series of specific attribute.
        attr can be list or string, if is a list it returns a dict
        """

        if not fq_idx:
            fq_idx   = self._fqidx[0]

        fq_idx   = self.cc8stats.check_idx(fq_idx)[0]

        attr_list = self.check_attr(attr)
        data_dict = {}
        nidx      = self.get_nidx(fq_idx=fq_idx, **nidx_kwargs)
        slow_max  = self.cc8stats.slow_max

        if "slowbnd" in attr_list:
            slokey = "/".join([fq_idx, "slow"])
            slow   = self._dout[slokey]
        
        if "bazmbnd" in attr_list:
            bazkey   = "/".join([fq_idx, "bazm"])
            baz      = self._dout[bazkey]

        for attr in attr_list:
            # get time series
            key = "/".join([fq_idx, attr])
            data = np.copy(self._dout[key])

            if attr == "slowmap":
                data_dict[attr] = data
                continue

            if attr == "slow":
                data[data>slow_max] = np.nan

            if attr == "rms":
                data[data<=0] = np.nan
                data = 10*np.log10(data)

            if isinstance(nidx, np.ndarray):
                for n, idx_nan in enumerate(np.isnan(nidx)):
                    if idx_nan:
                        if attr in ("slowbnd", "bazmbnd"):
                            data[n,:] = [np.nan, np.nan]
                        else:
                            data[n] = np.nan
                    else:
                        if attr == "slowbnd":
                            x = slow[n]
                            x0, x1 = data[n,:]
                            d1, d2 = x-x0, x1-x
                            data[n,:] = [d1, d2]

                        if attr == "bazmbnd":
                            x = baz[n]
                            x0, x1 = data[n,:]
                            d1, d2 = x-x0, x1-x
                            data[n,:] = [d1, d2]

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


    def get_pdf(self, attr, fq_idx=None, rms_in_db=True, vmin=None, vmax=None, data=None, **nidx_kwargs):
        """
        Return the PDF of a scalar attribute
        """

        if not isinstance(data, np.ndarray):
            attr = self.check_attr(attr)[0]
            data = self.get_data(attr, fq_idx=fq_idx, rms_in_db=rms_in_db, **nidx_kwargs)
        
        data = data[np.isfinite(data)]
        data = data.reshape(-1,1)
        
        if not vmin:
            vmin = data.min()

        if not vmax:
            vmax = data.max()

        space = np.linspace(vmin, vmax, 1000)
        pdf = get_PDF(data, space)

        return pdf, space


    def plot(self, fq_idx=None, show_title=True, maac_th=0.5, max_err=0.0, baz_int=[], rms_lim=[], **fig_kwargs):

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx

        fq_band = self.cc8stats.fq_bands[int(fq_idx)-1]

        attr_list = ["rms", "maac", "slow", "bazm", "slowbnd", "bazmbnd"]
        datattr   = self.get_data(attr_list, fq_idx=fq_idx, max_err=max_err, maac_th=maac_th, baz_int=baz_int)
        
        datapdf = {
            "slow":self.get_pdf("slow", vmin=0, vmax=slomax, data=datattr["slow"]),
            "bazm":self.get_pdf("bazm", vmin=0, vmax=360, data=datattr["bazm"])
        }

        fig_kwargs["maac_rv"] = self.get_data("maac", fq_idx=fq_idx, max_err=max_err, maac_th=maac_th, maac_rv=True, baz_int=baz_int)
        
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
    

    def get_beamform(self, starttime, endtime, slow, baz, fq_idx=None, fq_band=[], **fig_kwargs):
        """
        Return a waveform shifted between starttime and endtime for a specific slowness and back-azimuth
        """

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx

        if not fq_band:
            fq_band = self.cc8stats.fq_bands[int(fq_idx)-1]

        arr, exclude_locs = self.cc8stats.get_array()
        deltas, _ = arr.get_deltatimes(slow, baz, slomax, sloint, exclude_locs=exclude_locs)
        stream = arr.get_stream(starttime, endtime, prefilt=fq_band, toff_sec=600, exlcude_locs=exclude_locs)

        # shift stream
        wvfm_dict = {}
        interval  = endtime - starttime
        for delta, tr in zip(deltas, stream):
            of_npts = int(600*self.cc8stats.sample_rate)
            d_npts  = int(delta*self.cc8stats.sample_rate)
            data    = tr.get_data()
            data_sh = data[of_npts+d_npts:-of_npts+d_npts]
            wvfm_dict[tr.stats.location] = data_sh

        duration = interval.total_seconds()
        time = np.linspace(0, duration, len(data_sh))
        fig = window_wvfm(wvfm_dict, time, None, None, **fig_kwargs)

        return fig


    def plot_wvfm(self, ntime=None, off_sec=0, fq_idx=None, show_title=True, **fig_kwargs):
        """
        Plot shifted traces
        """
        
        arr, exclude_locs = self.cc8stats.get_array()

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx

        fq_band = self.cc8stats.fq_bands[int(fq_idx)-1]

        data_dict = self.get_data(["maac", "slow", "bazm"], fq_idx=fq_idx)
        maac  = data_dict["maac"]
        slow  = data_dict["slow"]
        baz   = data_dict["bazm"]
        
        if not ntime:
            ntime = np.argmax(maac)
            print(f" Best MAAC [{maac[ntime]:.1f}] is at position {ntime} with slow {slow[ntime]:.2f} and baz {baz[ntime]:.1f}")

        bazt    = baz[ntime]
        slowt   = slow[ntime]
        maact   = maac[ntime]
        half_w  = dt.timedelta(seconds=float(self.cc8stats.window/2))
        start   = self._dout["dtime"][ntime] - half_w
        end     = start + half_w + half_w

        # get stream and delta array
        stream = arr.get_stream(start, end, prefilt=fq_band, toff_sec=600, exlcude_locs=exclude_locs)
        deltas, _ = arr.get_deltatimes(slowt, bazt, slomax, sloint, exclude_locs=exclude_locs)

        # shift stream
        wvfm_dict = {}
        dtoff = dt.timedelta(seconds=off_sec)
        for delta, tr in zip(deltas, stream):
            starttime = start + dt.timedelta(seconds=delta) - dtoff
            endtime = starttime + half_w + half_w + dtoff + dtoff
            wvfm_dict[tr.stats.location] = tr.get_data(starttime=starttime, endtime=endtime)
    
        # create time array
        duration = self.cc8stats.window + 2*off_sec
        time = np.linspace(0, duration, int(duration*self.cc8stats.sample_rate)+1)
        
        # make fig
        if show_title:
            title = f"Fq {fq_band} :: Slomax/Sloint [{slomax}/{sloint}] \n MAAC {maact:.1f} :: Slow {slowt:.2f} [s/km] :: Baz {bazt:.1f}"
            fig_kwargs["title"] = title

        startw = duration/2 - self.cc8stats.window/2
        endw   = startw + self.cc8stats.window

        fig = window_wvfm(wvfm_dict, time, startw, endw, **fig_kwargs)

        return fig


    def plot_smap(self, ntime=None, fq_idx=None, show_title=True, **fig_kwargs):

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        if not fq_idx:
            fq_idx = self._fqidx[0]
        else:
            assert fq_idx in self._fqidx

        fq_band = self.cc8stats.fq_bands[int(fq_idx)-1]

        data_dict = self.get_data(["maac", "slow", "bazm", "slowmap"], fq_idx=fq_idx)
        data  = data_dict["slowmap"]
        maac  = data_dict["maac"]
        slow  = data_dict["slow"]
        baz   = data_dict["bazm"]

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
    

    def prob_slowmap(self, fq_idx=None, nidx=None, **nidx_kwargs):

        if not fq_idx:
            fq_idx   = self._fqidx[0]
        
        fq_band   = self.cc8stats.fq_bands[int(fq_idx)-1]
        maac_th  = nidx_kwargs.get("maac_th", 0)
        maac_rv  = nidx_kwargs.get("maac_rv", False)

        # get full slowmap
        smapkey  = "/".join([fq_idx, "slowmap"])
        slowmap  = self._dout[smapkey]

        # get filtered slowmap
        if not isinstance(nidx, np.ndarray):
            nidx    = self.get_nidx(**nidx_kwargs)

        fsmap   = slowmap[np.isfinite(nidx)]

        # filter each slowmap
        for nix in range(fsmap.shape[0]):
            if maac_rv:
                fsmap[nix,:,:] = np.where(fsmap[nix,:,:]<maac_th,fsmap[nix,:,:],np.nan)
            else:
                fsmap[nix,:,:] = np.where(fsmap[nix,:,:]>maac_th,fsmap[nix,:,:],np.nan)

        # do pdf map
        nites = self.cc8stats.nro_slow_bins
        pdfmap = np.zeros((nites,nites))
        for ii in range(nites):
            for jj in range(nites):
                data = fsmap[:,ii,jj]
                data = data[np.isfinite(data)]
                if data.shape[0] > 1:
                    pdfmap[ii,jj] = get_Stats(data)[3]

        sloint = self.cc8stats.slow_int
        slomax = self.cc8stats.slow_max

        fig_kwargs = {}
        fig_kwargs["cmap"] = "gist_earth_r"
        fig_kwargs["vlim"] = [np.nanmin(pdfmap), np.nanmax(pdfmap)]
        fig_kwargs["bar_label"] = "PDF MAAC > 0.7"

        starttime = self._dout["dtime"][0]
        endtime   = self._dout["dtime"][-1]
        fig_kwargs["title"] = f"{starttime}  -- {endtime} \n Fq {fq_band} :: Slomax/Sloint [{slomax}/{sloint}]"

        fig = simple_slowmap(pdfmap, sloint, slomax, **fig_kwargs)

        return pdfmap