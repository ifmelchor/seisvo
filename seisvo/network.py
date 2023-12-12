#!/usr/bin/python3
# coding=utf-8

import os
import glob
import json
import numpy as np
import datetime as dt

from  seisvo import seisvo_paths
from .stats import NetworkStats
from .obspyext import Stream2
from .station import Station
from .sap import _new_CC8
from .signal import SSteps, get_freq, array_response, get_CSW, get_CC8, get_PSD, array_delta_times
from .lte.base import _new_LTE
from .utils import nCPU
from .plotting.array import location_map, traces_psd, simple_slowmap, beamform_wvfm, plot_slowmap

default_slowarg = {
    "slomax":4.0, 
    "sloint":0.1, 
    "fq_band":[1., 3.], 
    "slow0":[0.,0.], 
    "cc_thres":0.05, 
    "exclude_locs":[]
}

def _slowarg(slowarg):
    if slowarg:
        slowarg0 = default_slowarg
        for key in list(slowarg0.keys()):
            if key in list(slowarg.keys()):
                slowarg0[key] = slowarg[key]
        slowarg = slowarg0
    else:
        slowarg = default_slowarg

    return slowarg


def get_network(net_code):
    availabel_network_codes = []
    
    for netfile in glob.glob(seisvo_paths["networks"]+"/*.json"):
        netcode = os.path.basename(netfile).split('.')[0]
        availabel_network_codes.append(netcode)
        
        if netcode == net_code:
            with open(netfile, 'r') as f:
                network_dict = json.load(f)

            network_dict["file"] = netfile
            return NetworkStats(netcode, network_dict)
    
    print(" warn :: no network file found ")
    print(f"    loaded networks are: {availabel_network_codes}")


class Network(object):
    def __init__(self, net_code):
        self.stats = get_network(net_code)


    def __str__(self):
        return self.stats.__str__()


    def __len__(self):
        return len(self.stats.stations)


    def __getitem__(self, item):
        if isinstance(item, int):
            return Station(self.stats.stations[item])
        
        if isinstance(item, str):
            sta_code = item.split('.')[0]
            
            try:
                loc = item.split('.')[1]
            except IndexError:
                loc = ''
            
            return self.get_sta(sta_code, loc=loc)
    

    def check_stalist(self, stalist, time_interval, return_sta=False):

        assert isinstance(stalist, (tuple, list))

        true_sta_list = []
        station_obj_list = []
        starttime = time_interval[0]
        endtime = time_interval[1]

        for sta in stalist:
            sta_sp = sta.split(".")
            sta_code = sta_sp[0]
            if len(sta_sp) > 1:
                sta_loc  = sta_sp[1]
            else:
                sta_loc  = ""
            sta_id = '.'.join([sta_code, sta_loc])

            if sta_id in ['.'.join([s.stats.code, s.stats.location]) for s in self]:
                sta_ob = self.get_sta(sta_code, loc=sta_loc)

                if time_interval:
                    if sta_ob.stats.starttime <= starttime and sta_ob.stats.endtime >= endtime:
                        station_obj_list.append(sta_ob)
                        true_sta_list.append(sta)
                else:
                    station_obj_list.append(sta_ob)
                    true_sta_list.append(sta)

        if return_sta:
            return true_sta_list, station_obj_list

        else:
            return true_sta_list


    def get_sta(self, sta_code, loc=''):
        """
        Get a station object
        """
        for sta in self:
            if sta_code == sta.stats.code and loc == sta.stats.location:
                return sta

        print(" station not found")
        return


    def get_stream(self, starttime, endtime, sta_code=[], component="Z", toff_sec=0, avoid_exception=False, return_stats=False, **st_kwargs):
        """
        This code get Stream object for the network.
        sta_code can be a list or a string
        """

        if toff_sec:
            starttime = starttime - dt.timedelta(seconds=toff_sec)
            endtime = endtime + dt.timedelta(seconds=toff_sec)

        if isinstance(sta_code, str):
            sta_code = [sta_code]

        stream = Stream2()
        stats = []
        for sta in self:
            if (sta_code and sta.stats.code in sta_code) or not sta_code:
                if component:
                    st_kwargs["channel"] = sta.get_chan(component)     
                
                try:
                    stream += sta.get_stream(starttime, endtime, **st_kwargs)
                    stats.append(sta.stats)
                except:
                    if avoid_exception:
                        return None
                    else:
                        print(f" error reading stream in {sta.stats.id}")
                        continue

        if return_stats:
            return stream, stats
        else:
            return stream


    def check_files(self, startday, endday, sta_code, loc=None, plot=False):
        """
        This code plot the availability of the stations of the network.
        :param startdate: datetime
        :param enddate: datetime
        :param sta_code: string
        :param loc: string, optional
        :param plot: boolean, by default False
        """
        
        # if plot:
        #     from seisvo.utils.plotting import plot_check

        sta_list = []
        for sta in self:
            if sta_code == sta.stats.code:
                if loc:
                    if isinstance(loc, str):
                        if sta.stats.location == loc:
                            for chan in sta.stats.channels:
                                sta_list += ['%s.%s' % (sta.stats.id, chan)]
                    
                    if isinstance(loc, (list, tuple)):
                        for l in loc:
                            if sta.stats.location == l:
                                for chan in sta.stats.channels:
                                    sta_list += ['%s.%s' % (sta.stats.id, chan)]

        if not sta_list:
            print(' "checking ERROR :: Station list is empty')
            return False

        day_list = [startday + dt.timedelta(days=k) for k in range(0,(endday-startday).days+1)]
        
        list_missing_days = []
        for i, day in enumerate(day_list):
            for j, sta_code_chan in enumerate(sta_list):
                sta_code = sta_code_chan.split('.')[1]
                chan     = sta_code_chan.split('.')[-1]
                sta_loc  = sta_code_chan.split('.')[2]
                sta      = self.get_sta(sta_code, loc=sta_loc)
                sta_loc_day_file = sta.__read_file__(chan, date=day)
                
                if not sta_loc_day_file:
                    missing_str = ".".join([sta_code,sta_loc,chan])
                    missing_str += f"/{day}"
                    list_missing_days += [missing_str]
        
        return list_missing_days


    def get_csw(self, starttime, endtime, sta_list, window, win_freq=0.33, olap=0.5, return_vt=False, **kwargs):
        """
        Compute the cross spectral width (CSW) between start and end time.
        window [float, in seconds] > 0
        if return_pdf is True, return PDF of all window segments of the trace
        """

        assert starttime < endtime

        sta_list = self.check_stalist(sta_list, (starttime,endtime))
        assert sta_list

        # load kwargs
        rm_sens     = kwargs.get('rm_sens', True)
        rm_resp     = kwargs.get('rm_resp', False)
        sample_rate = kwargs.get('sample_rate', 40)
        njobs       = kwargs.get('njobs', 1)
        pad         = kwargs.get("pad", 1.0)
        time_bandw  = kwargs.get('time_bandwidth', 3.5)
        fq_band     = kwargs.get('fq_band', (0.5, 10))

        stream = self.get_stream(starttime, endtime, sta_code=sta_list,\
            avoid_exception=False, remove_sensitivity=rm_sens, sample_rate=sample_rate)

        ans = None
        if stream and stream.get_bounds() == (starttime,endtime):
            data = stream.to_array()
            ans = get_CSW(data, sample_rate, window*sample_rate,\
                olap=olap, fq_band=fq_band, NW=time_bandw, pad=pad,\
                win_freq=win_freq, return_vt=return_vt)

        return ans


    def lte(self, starttime, endtime, sta_list, window, subwindow, win_olap=0.5,\
        subw_olap=0.75, interval=None, **kwargs):
        """ Compute Network LTE file

        Parameters
        ----------
        starttime : datetime
        
        endtime : datetime
        
        sta_list : list of strings
            list of stations id

        window : int [min]
            length of the time window (in min) to reduce
        
        subwindow : int [min]
            length of the time window (in min) for moving average over the window. 
            If ``subwindow=0`` no moving average is applied.
        
        subwindow_olap : float
            overlap percent for moving average over the interval, by default 0.75

        Returns
        -------
        LTE
            LTE object
        """

        assert starttime < endtime
        assert window > subwindow

        sta_list = check_stalist(sta_list, time_interval=(starttime,endtime))
        assert sta_list
        
        # load kwargs
        rm_sens     = kwargs.get('rm_sens', True)
        rm_resp     = kwargs.get('rm_resp', False)
        sample_rate = kwargs.get('sample_rate', 40)
        njobs       = kwargs.get('njobs', 1)
        pad         = kwargs.get("pad", 1.0)
        time_bandw  = kwargs.get('time_bandwidth', 3.5)
        fq_band     = kwargs.get('fq_band', (0.5, 10))
        file_name   = kwargs.get("file_name", None)
        out_dir     = kwargs.get("out_dir", './')

        # defining base params
        ltebase = dict(
            id              = self.stats.code,
            type            = "network",
            stations        = tuple(sta_list),
            starttime       = starttime.strftime('%Y-%m-%d %H:%M:%S'),
            endtime         = endtime.strftime('%Y-%m-%d %H:%M:%S'),
            interval        = interval,
            window          = window,
            window_olap     = win_olap,
            subwindow       = subwindow,
            subwindow_olap  = subw_olap,
            sample_rate     = int(sample_rate),
            pad             = pad,
            time_bandwidth  = time_bandw,
            fq_band         = fq_band,
            rm_sens         = rm_sens
        )

        ss = SSteps(starttime, endtime, window, interval=interval, win_olap=window_olap, subwindow=subwindow, subw_olap=subwindow_olap)

        # check for integer values
        assert ss.nro_intervals.is_integer() and ss.nro_intervals >= 1
        assert ss.total_nwin.is_integer() and ss.total_nwin >= 1
        assert ss.int_nwin.is_integer() and ss.int_nwin >= 1
        assert ss.subwindow >= 1

        ltebase["lwin"] = int(ss.window*sample_rate)
        ltebase["nwin"] = int(ss.int_nwin)
        ltebase["wadv"] = float(ss.win_adv)
        ltebase["nro_intervals"] = int(ss.nro_intervals)
        ltebase["nro_time_bins"] = int(ss.total_nwin)
        ltebase["lswin"] = int(ss.subwindow*sample_rate)
        ltebase["nswin"] = int(ss.nsubwin)
        ltebase["swadv"] = float(ss.subw_adv)
        nfs, _ = get_freq(ltebase["lswin"], sample_rate, fq_band=fq_band, pad=pad)
        ltebase["nro_freq_bins"] = len(nfs)

        if njobs >= nCPU or njobs == -1:
            njobs = nCPU - 2
        
        if njobs > nro_intervals:
            njobs = nro_intervals
            print(f"warn  ::  njobs set to {njobs}")

        if njobs >= nCPU or njobs == -1:
            njobs = nCPU - 2

        # create hdf5 file and process data
        if not file_name:
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (self.stats.code, starttime.year, starttime.timetuple().tm_yday, endtime.year, endtime.timetuple().tm_yday, window)
        
        file_name_full = os.path.join(out_dir, file_name)
        if file_name_full.split('.')[-1] != 'lte':
            file_name_full += '.lte'

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)
        
        lte = _new_LTE(self, file_name_full, ltebase, njobs)
        return lte


class Array(Network):
    def __init__(self, net_code, sta_code):
        super().__init__(net_code)

        assert sta_code in self.stats.stations_code
        self.sta_code = sta_code

        # change self
        to_remove = [sta.stats for sta in self if sta.stats.code != sta_code]
        for sta_stats in to_remove:
            self.stats.stations.remove(sta_stats)

        # compute UTM position
        self.utm = self.stats.get_latlon(utm=True, in_km=True, hide_sta=True)
        self.locs = list(self.utm.keys())

        # check sample rate
        fslist = list(set([sta.stats.sample_rate for sta in self]))
        if len(fslist) > 1:
            fslist = [None]
            print(" warn :: multiple sampling rates defined in config file")
        self.sample_rate = fslist[0]

    
    def get_sta(self, loc):
        return super().get_sta(self.sta_code, loc=loc)


    def get_utm(self, exclude_locs=[]):

        if exclude_locs:
            locs = [loc for loc in self.locs if loc not in exclude_locs]
        else:
            locs = self.locs

        utmloc = {"x":[],"y":[]}
        for loc, utm in self.utm.items():
            if loc in locs:
                utmloc["x"].append(utm["easting"])
                utmloc["y"].append(utm["northing"])

        return locs, utmloc


    def aperture(self, exclude_locs=[]):
        a = 0

        for i in range(len(self)):
            iloc = self.locs[i]

            if iloc in exclude_locs:
                continue
            
            for j in range(i,len(self)):
                jloc = self.locs[j]
                
                if jloc in exclude_locs:
                    continue

                if i != j:
                    x1 = self.utm[iloc]["easting"]
                    x2 = self.utm[jloc]["easting"]
                    y1 = self.utm[iloc]["northing"]
                    y2 = self.utm[jloc]["northing"]
                    d = np.sqrt(np.abs(x1-x2)**2 + np.abs(y2-y1)**2)
                    
                    if d > a:
                        a = d
        return a

    
    def response(self, slomax, sloint, fq_band=(1.,10.), fq_int=0.1, exclude_locs=[], plot=True):
        """
        Return array response dictionary
        """
        xutm, yutm = [],[]

        for loc, utm in self.utm.items():
            if loc not in exclude_locs:
                xutm.append(float(utm["easting"]))
                yutm.append(float(utm["northing"]))

        ans = array_response(xutm, yutm, slomax, sloint, fq_band=fq_band, fq_int=fq_int)
        power = ans["power"]
        
        if plot:
            simple_slowmap(power/power.max(), slomax, sloint, title=fq_band)
        
        return ans


    def get_stream(self,starttime, endtime, component="Z", toff_sec=0, return_stats=False, exclude_locs=[], **st_kwargs):
        ans = super().get_stream(starttime, endtime, component=component, toff_sec=toff_sec, return_stats=return_stats, **st_kwargs)

        if return_stats:
            stream, stats = ans
        else:
            stream = ans

        if exclude_locs:
            new_stream = Stream2()
            for tr in stream:
                if tr.stats.location not in exclude_locs:
                    new_stream.append(tr)
            stream = new_stream

        if return_stats:
            return stream, stats
        else:
            return stream


    def get_cc8(self, starttime, endtime, window, overlap, slowarg={}, toff_sec=10):
        """
        compute CC8 algorithm
        window in seconds (float) and overlap between 0 an 1.
        """

        slowarg = _slowarg(slowarg)
        nite    = 1 + 2*int(slowarg["slomax"]/slowarg["sloint"])

        # locations and positions
        _, utmloc = self.get_utm(exclude_locs=slowarg["exclude_locs"])
        
        # compute steps
        ss     = SSteps(starttime, endtime, window, win_olap=overlap, logfile=True)
        ssdict = ss.to_dict()

        # get data array
        stream = self.get_stream(starttime, endtime, toff_sec=toff_sec,\
            exclude_locs=slowarg["exclude_locs"], sample_rate=self.sample_rate, 
            avoid_exception=True)
        data = stream.to_array(detrend=True)
        lwin = int(ss.window * self.sample_rate)

        # wrapper into julia
        ans = get_CC8(data, self.sample_rate, utmloc["x"], utmloc["y"], slowarg["fq_band"], 
            slowarg["slomax"], slowarg["sloint"], lwin=lwin, nwin=ssdict["nwin"], 
            nadv=ssdict["wadv"], toff=toff_sec, cc_thres=slowarg["cc_thres"])
        
        return ans


    def cc8(self, starttime, endtime, window, overlap, interval=30, slow_max=3., slow_int=0.05, fq_bands=[(1.,3.)], cc_thres=0.05, exclude_locs=[], slowmap=True, **kwargs):
        
        """ Compute CC8 file

        Parameters
        ----------
        starttime : datetime
        
        endtime : datetime
        
        window : int [sec]
            length of the time window (in sec)
        
        overlap : float [0--1)
            overlap percent
        
        interval : int (minutes)
            length of the interval to store in memory
        
        exclude_locs : list (optional)
            list of strings to exclude locations
        
        slow_max : list 
            maximum slowness in km/s, by default [3.]

        slow_int : list
            slowness intervals in km/s, by default [0.1]
        
        fq_bands : list of tuples
            frequency bands (in Hz) to analyze, by default [(1.,5.)]
        
        cc_thres : float
            therhold level to keep computing slowness, by default 0.75

        Returns
        -------
        CC8
            CC8 object
        """

        # do simple checks
        assert starttime < endtime
        assert slow_max > slow_int
        assert window/60 < interval

        # load parameters
        sample_rate = self.sample_rate
        njobs       = kwargs.get('njobs', 1)
        toff_sec    = kwargs.get('toff_sec', 3)
        fileout     = kwargs.get("fileout", None)
        
        # compute slowness invervals
        nites = 1 + 2*int(slow_max/slow_int)
        
        # define locations and positions
        locs, utmloc = self.get_utm(exclude_locs=exclude_locs)
        
        # define the header of the cc8 file
        cc8base = dict(
            id              = '.'.join([self.stats.code, self.sta_code]),
            locs            = locs,
            utm             = utmloc,
            starttime       = starttime.strftime('%Y-%m-%d %H:%M:%S.%f'),
            endtime         = endtime.strftime('%Y-%m-%d %H:%M:%S.%f'),
            interval        = interval,
            window          = window,
            overlap         = overlap,
            sample_rate     = sample_rate,
            fq_bands        = fq_bands,
            slow_max        = slow_max,
            slow_int        = slow_int,
            nro_slow_bins   = nites,
            cc_thres        = cc_thres,
            slowmap         = slowmap,
            toff_sec        = toff_sec
        )

        # compute the steps and save info
        ssdict = SSteps(starttime, endtime, window, interval=interval,\
            win_olap=overlap, logfile=True).to_dict()
        
        cc8base["lwin"] = int(window*sample_rate)
        cc8base["nwin"] = ssdict["nwin"]
        cc8base["nadv"] = ssdict["wadv"]
        cc8base["last_nwin"] = ssdict["last_nwin"]
        cc8base["nro_intervals"] = ssdict["nro_intervals"]
        cc8base["nro_time_bins"] = ssdict["total_nwin"]
        cc8base["int_extra_sec"] = ssdict["int_extra_sec"]

        # set a valid njob
        if njobs >= nCPU or njobs == -1:
            njobs = nCPU - 2
        
        if njobs > cc8base["nro_intervals"]*len(fq_bands):
            njobs = cc8base["nro_intervals"]*len(fq_bands)
            print(f"warn  ::  njobs set to {njobs}")
        
        # define name
        if not fileout:
            fileout = "./%s.%s%03d-%s%03d.cc8" % (cc8base["id"],\
                starttime.year, starttime.timetuple().tm_yday,\
                endtime.year, endtime.timetuple().tm_yday)
        
        # check if filename ends with .cc8
        if fileout.split('.')[-1] != "cc8":
            fileout += ".cc8"

        # check if file exist and remove
        if os.path.isfile(fileout):
            os.remove(fileout)
            print(' file %s removed.' % fileout)
        
        # create a cc8 file
        _new_CC8(self, fileout, cc8base, njobs)
        
        return None


    def plot_map(self, exclude_locs=[], show=True):
        return location_map(self, exclude_locs=exclude_locs, show=show)


    def get_excluded_locs(self, loc_list):
        exclude_loc = [loc for loc in self.locs if loc not in loc_list]
        return exclude_loc


    def psd(self, starttime, endtime, window, olap=0.25, fq_band=(0.5,10.), exclude_locs=[], plot=True, show=True, **st_kwargs):

        stream = self.get_stream(starttime, endtime, exclude_locs=exclude_locs, **st_kwargs)
        
        if not stream:
            return None

        psd_dict = {}
        freq     = None
        for tr in stream:
            try:
                fs = int(tr.stats.sampling_rate)
                psd, freq = get_PSD(tr.data, fs, lwin=int(fs*window), olap=olap,\
                    fq_band=fq_band)
                psd_dict[tr.stats.location] = psd
            except:
                pass

        if isinstance(freq, np.ndarray):
            if plot:
                fig = traces_psd(psd_dict, freq, db_scale=True, title=f"{starttime} -- {endtime}", show=show)
                return fig

            else:
                return psd_dict, freq

        else:
            return None


    def deltatimes(self, slow, baz, slowarg={}, return_xy=False):
        """
        compute the waveforms displaced by given a slowness vector (slow, baz)
        """

        slowarg  = _slowarg(slowarg)
        _, utmloc   = self.get_utm(exclude_locs=slowarg["exclude_locs"])
        ans, tol = array_delta_times(slow, baz, slowarg["slomax"], slowarg["sloint"],\
            self.sample_rate, utmloc["x"], utmloc["y"], pxy0=[0.,0.], return_xy=return_xy)

        print(f" [deltatimes] error = {tol}")

        if return_xy:
            return ans, tol
        
        else:
            return ans/self.sample_rate, tol


    def beamform(self, starttime, endtime, slow, baz, slowarg={}, return_full=False, taper=False, plot=True, **fig_kwargs):
        """
        Return a waveform shifted between starttime and endtime for a specific slowness and back-azimuth
        """

        slowarg   = _slowarg(slowarg)
        # deltas, tol = self.deltatimes(slow, baz, slowarg=slowarg)
        # print(f" [deltatimes] error = {tol}")
        window    = (endtime-starttime).total_seconds()
        full, deltas = self.slowmap(starttime, window, slowarg=slowarg, return_dt=True, plot=False)
        stream    = self.get_stream(starttime, endtime, prefilt=slowarg["fq_band"], toff_sec=10, exclude_locs=slowarg["exclude_locs"])

        # shift stream
        wvfm_dict = {}
        for delta, tr in zip(deltas, stream):
            of_npts  = int(10*self.sample_rate)
            d_npts   = int(delta*self.sample_rate)
            data     = tr.get_data()
            data_sh  = data[of_npts+d_npts:-of_npts+d_npts]
            wvfm_dict[tr.stats.location] = data_sh
        
        # beamform
        bf = np.empty((len(data_sh), len(deltas)))
        for n, (_, data) in enumerate(wvfm_dict.items()):
            bf[:, n] = data/np.abs(data).max()
        
        suma = np.sum(bf, axis=1) / bf.shape[1]
        if taper:
            suma *= np.hanning(len(suma))

        duration = (endtime - starttime).total_seconds()
        time     = np.linspace(0, duration, len(data_sh))

        if return_full:
            full["wvfm"] = wvfm_dict
            return full, suma

        if plot:
            fig = beamform_wvfm(wvfm_dict, suma, time, **fig_kwargs)
            return fig
        else:
            return wvfm_dict, suma, time


    def slowmap(self, starttime, window, slowarg={}, return_dt=False, plot=True, show_title=True, **fig_kwargs):
        """
        Compute the slowness map

        return_dt should be True/False or "xy"
        """

        slowarg = _slowarg(slowarg)

        # init some parameters
        toff_sec = 10
        nite     = 1 + 2*int(slowarg["slomax"]/slowarg["sloint"])
        endtime  = starttime + dt.timedelta(seconds=window)
        lwin     = int(window * self.sample_rate)
        nwin     = 1
        nadv     = 0
        
        # get data
        stream = self.get_stream(starttime, endtime, toff_sec=toff_sec,\
            exclude_locs=slowarg["exclude_locs"], sample_rate=self.sample_rate, avoid_exception=True)
        data   = stream.to_array(detrend=True)

        # locations and positions
        _, utmloc = self.get_utm(exclude_locs=slowarg["exclude_locs"])
        
        # compute CC8 algorithm
        ans    = get_CC8(data, self.sample_rate, np.array(utmloc["x"]), np.array(utmloc["y"]),\
            slowarg["fq_band"], slowarg["slomax"], slowarg["sloint"], lwin=lwin, nwin=nwin, \
            nadv=nadv, cc_thres=slowarg["cc_thres"], toff=toff_sec, slow0=slowarg["slow0"])

        if return_dt:
            slow = ans["slow"][0]
            baz  = ans["baz"][0]

            if return_dt == "xy":
                deltas, tol = self.deltatimes(slow, baz, slowarg=slowarg, return_xy=True)
            else:
                deltas, tol = self.deltatimes(slow, baz, slowarg=slowarg, return_xy=False)
            
            return ans, deltas
        
        if plot:
            if not fig_kwargs:
                fig_kwargs = {}
            
            if show_title:    
                rms    = 10*np.log10(ans["rms"][0])
                error  = ans["error"][0]
                maact  = ans["maac"][0]
                slowt  = ans["slow"][0]
                bazt   = ans["baz"][0]

                fig_kwargs["title"] = f' {starttime}  [{window} sec] :: Fq {slowarg["fq_band"]} :: Slomax/Sloint [{slowarg["slomax"]}/{slowarg["sloint"]} s/km] \n RMS {rms:.1f} [dB] :: MAAC {maact:.2f} :: SLOW {slowt:.2f} [s/km] :: BAZ {bazt:.1f} :: ERR {error:.1f} \n'
        
            slomap   = ans["slowmap"][0,:,:]

            fig = simple_slowmap(slomap, slowarg["sloint"], slowarg["slomax"], **fig_kwargs)
            return fig, ans
        else:
            return ans


    def plot(self, starttime, window, taper=False, offsec=3, slowarg={}):
        return plot_slowmap(self, starttime, window, taper=taper, offsec=offsec, slowarg=slowarg)



class SoundArray(Network):
    def __init__(self, net_code, sta_code, chan='P', **kwargs):
        super().__init__(net_code)
        self.sta_code = sta_code
        self.central_sta = None
        self.chan = chan

        self.__setstacode__(sta_code, chan)
        self.__check__()
        self.set_central(sta_code=kwargs.get('central_sta'))
        self.set_model(model=kwargs.get('model', None))


    def __check__(self):
        # check sample rate
        sample_rate = list(set([sta.stats.sampling_rate for sta in self.stations]))
        if len(sample_rate) == 1:
            self.sample_rate = int(sample_rate[0])
        else:
            print('warn: different sampling rate in the array!')
            self.sample_rate = int(min(sample_rate))


    def set_central(self, sta_code=None):
        """By default, the central station of the array is defined in the config file. But it can be stated manually specifying the sta_code.

        Parameters
        ----------
        sta_code : [str], optional
        """

        self.lat0 = None
        self.lon0 = None
        for sta in self.station:
            if not sta_code:
                if sta.stats.central:
                    self.lat0 = sta.get_latlon(degree=False)[0]
                    self.lon0 = sta.get_latlon(degree=False)[1]
                    self.central_sta = sta.stats.id
            else:
                if sta.stats.code == sta_code:
                    self.lat0 = sta.get_latlon(degree=False)[0]
                    self.lon0 = sta.get_latlon(degree=False)[1]
                    self.central_sta = sta.stats.id
        
        if not self.lat0:
            print("warn: central station is not defined")
            self.central_sta = None


    def set_model(self, model=None):
        """By default, the model is a dictionary with keys: radii, vel_air, h_src, src_dgr.

        Parameters
        ----------
        model : [dict], optional
        """
        self.model = infrasound_model_default

        if model:
            if isinstance(model, dict):
                for k in model.keys():
                    self.model[k] = model[k]
            else:
                print("model must be a dictionary. Default parameters are used.")
        
        print('\n   -- Input Model --')
        for key, item in self.model.items():
            print(f'  {key:>8s} : ', item)

        h_mean = np.array([sta.stats.elev for sta in self.station]).mean()
        self.h_mean_ = h_mean
        h_diff = self.model['h_src'] - h_mean

        self.model_azimuths = np.array(np.arange(self.model['src_dgr'][0], self.model['src_dgr'][1]+1, 1), dtype='int')
        th = self.model_azimuths * np.pi / 180
        self.nro_srcs = len(th)
        x_src = self.lat0 + self.model['radii'] * np.sin(th)
        y_src = self.lon0 + self.model['radii'] * np.cos(th)

        self.dt_times = {}
        for sta in self:
            lat_lon = sta.get_latlon(degree=False)
            xx = lat_lon[0] - x_src
            yy = lat_lon[1] - y_src
            zz = sta.stats.elev - h_diff
            dist_src = np.sqrt(xx ** 2 + yy ** 2 + zz ** 2)
            self.dt_times[sta.stats.id] = dist_src / self.model['vel_air']

        dt_central = self.dt_times[self.central_sta]
        for sta in self:
            self.dt_times[sta.stats.id] = np.array(np.around((self.dt_times[sta.stats.id] - dt_central) * self.sample_rate), dtype='int')
    

    def get_stream(self, starttime, endtime, azm=None, time_pad=30, model=None, **kwargs):
        
        pad_delta = dt.timedelta(seconds=time_pad)
        delta = 1/self.sample_rate

        if model:
            self.set_model(model=model)
        
        if azm:
            if azm not in list(self.model_azimuths):
                raise ValueError('azm not found')
            azm_idx = list(self.model_azimuths).index(azm)

        stream = Stream2()
        for sta in self:
            if azm:
                azm_delta = dt.timedelta(seconds=int(self.dt_times[sta.stats.id][azm_idx]) * delta)
            else:
                azm_delta = dt.timedelta(seconds=0)

            stream += sta.get_stream(starttime-pad_delta+azm_delta, endtime+pad_delta+azm_delta, component='P', remove_response=True, **kwargs)

        return stream

    
    def ripepe_ccorr(self, time_width, overlap, starttime, interval=None, endtime=None, **kwargs):
        """This code compute the cross correlation between infrasound stations following the Ripepe algorithm.

        Parameters
        ----------
        starttime : [datetime]
        interval : [int], optional
            [in minutes], by default None
        endtime : [datetime], optional

        kwargs : fq_band (default: (0.5, 5)),
                 air_file: True/False or string (name of the file)
                 out_dir: string with path of the air_file
                 


        Raises
        ------
        ValueError
            [if interval and endtime are undefined]
        """


        if not interval and not endtime:
            raise ValueError('interval (in minutes) or endtime must be defined')

        if not endtime:
            time_delta = dt.timedelta(minutes=interval)
            endtime = starttime + time_delta
        
        if endtime < starttime:
            raise ValueError('endtime must be greater than starttime')
        
        fq_band_ = kwargs.get('fq_band', (0.5, 5))
        air_file = kwargs.get('air', False)
        sample_rate = kwargs.get('sample_rate', self[0].stats.sample_rate)

        info = {}
        info['code'] = '%s.%s' % (self.stats.code, self.sta_code)
        info['sampling_rate'] = sample_rate
        info['starttime'] = starttime.strftime('%Y-%m-%d %H:%M:%S')
        info['endtime'] = endtime.strftime('%Y-%m-%d %H:%M:%S')
        info['time_width'] = time_width
        info['overlap'] = overlap
        info['fq_band'] = fq_band_
        info['radii'] = self.model['radii']
        info['vel_air'] = self.model['vel_air']
        info['h_src'] = self.model['h_src']
        info['src_dgr'] = self.model['src_dgr']
        info['sensors_loc'] = self.stations_info

        ss = SSteps(starttime, endtime, time_width, overlap)
        ss.fit(best_interval=True, interval_range=(15,30))
        ss.print()
        interval = ss.interval_.total_seconds()/60

        if ss.nro_intervals%2 > 0:
            last_interval = ss.last_interval_
        else:
            last_interval = False
        
        info['interval'] = interval
        info['nro_intervals'] = ss.nro_intervals
        info['last_interval'] = ss.last_interval_.total_seconds()/60
        info['steps_per_interval'] = ss.steps_
        
        if ss.rest_ > 0:
            raise ValueError ('steprest_per_interval is not zero!')

        if air_file:
            if isinstance(air_file, str):
                file_name = air_file
            else:
                file_name = '%s.%s%03d-%s%03d.air' % (
                    info['code'],
                    starttime.year,
                    starttime.timetuple().tm_yday,
                    endtime.year,
                    endtime.timetuple().tm_yday
                    )
            out_dir = kwargs.get('out_dir', './')
            file_name_path = os.path.join(out_dir, file_name)

            if not os.path.isdir(out_dir):
                os.makedirs(out_dir)

            if os.path.isfile(file_name_path):
                os.remove(file_name_path)
            
            air_file = AiR.new(file_name_path, info, self.dt_times, ss)
            iter_items = iter(range(ss.total_steps()))

        else:
            dout = {
                'time':[],
                'mcorr':None,
                'p_avg':None,
                'p_max':None,
                'info':info,
                'model':self.model,
            }
            iter_items = None

        nint = int(np.floor(ss.nro_intervals))
        # digits = len(str(abs(nint)))
        start = starttime
        for i in range(nint):
            txt_to_print = f' {i+1:>5}/{nint:<5}  ::'
            
            end = start + ss.interval_
            
            stream = self.get_stream(
                start - dt.timedelta(seconds=ss.window_length),
                end + dt.timedelta(seconds=ss.window_length),
                prefilt=fq_band_,
                remove_response=True,
                sample_rate=sample_rate
                )
            
            ans = cross_corr(
                stream,
                self.dt_times,
                self.nro_srcs,
                ss,
                iter_bin=iter_items,
                air_file=air_file,
                last_interval=False,
                txt_print=txt_to_print
                )

            if not air_file:
                dout['time'] += ans[0]

                if dout['mcorr'] is None:
                    dout['mcorr'] = ans[1]
                    dout['p_avg'] = ans[2]
                    dout['p_max'] = ans[3]

                else:
                    dout['mcorr'] = np.concatenate((dout['mcorr'], ans[1]), axis=0)
                    dout['p_avg'] = np.concatenate((dout['p_avg'], ans[2]), axis=0)
                    dout['p_max'] = np.concatenate((dout['p_max'], ans[3]), axis=0)

            start += ss.interval_
        
        if last_interval:
            txt_to_print = ' "last_interval" ::'
            end = start + last_interval
            
            stream = self.get_stream(
                start - dt.timedelta(seconds=ss.window_length),
                end + dt.timedelta(seconds=ss.window_length),
                prefilt=fq_band_,
                remove_response=True
                )
            
            ans = cross_corr(
                stream,
                self.dt_times,
                self.nro_srcs,
                ss,
                iter_bin=iter_items,
                air_file=air_file,
                last_interval=True,
                txt_print=txt_to_print
                )

            if not air_file:
                dout['time'] += ans[0]

                if dout['mcorr'] is None:
                    dout['mcorr'] = ans[1]
                    dout['p_avg'] = ans[2]
                    dout['p_max'] = ans[3]

                else:
                    dout['mcorr'] = np.concatenate((dout['mcorr'], ans[1]), axis=0)
                    dout['p_avg'] = np.concatenate((dout['p_avg'], ans[2]), axis=0)
                    dout['p_max'] = np.concatenate((dout['p_max'], ans[3]), axis=0)
        
        if air_file:
            print('\n done.')
            return air_file
        
        return dout


    def get_available(self, start_time, endt_time):
        """El objetivo de esta función es que te devuelva aquellas fechas en las que hay datos. Pero no esta implementado aún y no se recomienda su uso.
        """

        missing_files = self.check_files(start_time, endt_time, self.sta_code)
        dmising = col.Counter(missing_files)
        missing_dates = list(set([dt.datetime.strptime(x.split('-')[1], '%Y%j') for x in missing_files]))
        missing_dates.sort()

        one_day = dt.timedelta(days=1)
        starttime_list = [start_time]
        endtime_list = [endt_time]
        recursive_days = []
        for i,j in zip(missing_dates[:-1], missing_dates[1:]):
            if j != i+one_day:
                starttime_list.append(i+one_day)
                endtime_list.append(i)
            else:
                recursive_days.append(i)

        if recursive_days:
            starttime_list.append(max(recursive_days)+2*one_day)
            endtime_list.append(min(recursive_days))

        starttime_list.sort()
        endtime_list.sort()

        return (starttime_list, endtime_list)

    
    def get_waveform(self, out, time_bin, azm=None, time_pad=30, fq_band=()):

        if not fq_band:
            fq_band = out['info']['fq_band']

        delta = 1/out['info']['sampling_rate']
        time = out['time'][time_bin]
        time_width = dt.timedelta(seconds=out['info']['time_width'])

        if not azm:
            ans = self.get_max_values(out)
            azm_bin = ans[1][time_bin]
        else:
            azm_bin = self.get_azimuth(out, azm)

        st = Stream2()
        for sta in self.station:
            if time_pad:
                pad_time = dt.timedelta(seconds=time_pad)
                starttime = time - pad_time
                endtime = time + pad_time

            else:
                sec = float((int(self.dt_times[sta.stats.id][azm_bin]) - 1) * delta)
                starttime = time + dt.timedelta(seconds=sec)
                endtime = time + dt.timedelta(seconds=sec) + time_width

            st += sta.get_stream(starttime, endtime, remove_response=True, prefilt=fq_band)

        return st
    

    def get_starttime(self, out, time_bin, azm=None):
        v_bars = {}
        delta = 1/out['info']['sampling_rate']
        time = out['time'][time_bin]
        time_width = dt.timedelta(seconds=out['info']['time_width'])

        if not azm:
            ans = self.get_max_values(out)
            azm_bin = ans[1][time_bin]
        else:
            azm_bin = self.get_azimuth(out, azm)

        for sta in self.station:
            dt_delta = dt.timedelta(seconds=float((int(self.dt_times[sta.stats.id][azm_bin]) - 1) * delta))
            v_bars[sta.stats.id] = (time + dt_delta, time + dt_delta + time_width)

        return v_bars
    

    def plot_gui(self, **kwargs):
        from seisvo.gui.giarr import gplot
        gplot(array=self, **kwargs)


    def get_azimuth(self, azm_bin):
        azm_list = list(self.model_azimuths)
        return  azm_list.index(azm_bin)


    @staticmethod
    def get_max_values(out, time_bin=None):
        mcorr = out['mcorr']
        p_max = out['p_max']
        p_avg = out['p_avg']
        azms = range(out['model'].get('src_dgr')[0], out['model'].get('src_dgr')[1])

        mcorr_max = []
        azm_max = []
        pmax = []
        pavg = []

        if not time_bin:
            for x in range(mcorr.shape[0]):
                r = mcorr[x, :]
                mcorr_max += [r.max()]
                src = np.argmax(r)
                azm_max += [azms[src]]
                pmax += [p_max[x, src]]
                pavg += [p_avg[x, src]]

        else:
            x = time_bin
            r = mcorr[x, :]
            mcorr_max += [r.max()]
            src = np.argmax(r)
            azm_max += [azms[src]]
            pmax += [p_max[x, src]]
            pavg += [p_avg[x, src]]

        return [np.array(mcorr_max), np.array(azm_max), np.array(pmax), np.array(pavg)]

