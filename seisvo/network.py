#!/usr/bin/python3
# coding=utf-8

import os
import glob
import json
import numpy as np
import datetime as dt

from seisvo import seisvo_paths
from .stats import NetworkStats
from .obspyext import Stream2
from .station import Station
from .sap import CC8
from .signal import SSteps, get_freq, array_response
from .lte.base import _new_LTE
from .utils import nCPU


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
    

    def check_stalist(self, stalist, time_interval=None, return_sta=False):

        assert isinstance(stalist, (tuple, list))

        true_sta_list = []
        station_obj_list = []

        if time_interval:
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
                    if sta_ob.starttime <= starttime and sta_ob.endtime >= endtime:
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
                    kwargs["channel"] = sta.get_chan(component)     
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


    # def gui(self, starttime, station_list, component="Z", delta=30, sde_file=None, **kwargs):

    #     from seisvo.plotting.gui.gnetwork import init_network_gui

    #     # check station_list
    #     true_station_list = []
    #     for sta in self.station:
    #         sta_id = '.'.join([sta.stats.code, sta.stats.loc])
    #         if sta_id in station_list:
    #             true_station_list += [sta_id]

    #     if not true_station_list:
    #         print(" no stations loaded!. Revise network file!")
    #         return
        
    #     if component not in ["Z", "N", "E"]:
    #         print(" componente must be Z, N, or E")
    #         return
        
    #     if sde_file:
    #         if sde_file.split('.')[-1] != '.db':
    #             sde_file += '.db'

    #     else:
    #         sde_file = os.path.join(seisvo_paths["database"], self.stats.code + '_sde.db')
        
    #     if isinstance(kwargs.get("specgram"), int):
    #         specgram = station_list[kwargs.get("specgram")]
    #         kwargs["specgram"] = specgram
        
    #     if specgram not in true_station_list:
    #         kwargs["specgram"] = None

    #     init_network_gui(self, true_station_list, starttime, delta, component, sde_file, **kwargs)


    def lte(self, starttime, endtime, sta_list, window, subwindow, win_olap=0.5, subw_olap=0.75, **kwargs):
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

        sta_list, sta_obj = check_stalist(sta_list, time_interval=(starttime,endtime), return_sta=True)

        if not sta_ob_list:
            print("error:: no data to proces")
            return None
        
        # load kwargs
        sample_rate = kwargs.get('sample_rate', 40)
        njobs       = kwargs.get('njobs', 1)
        pad         = kwargs.get("pad", 1.0)
        time_bandw  = kwargs.get('time_bandwidth', 3.5)
        fq_band     = kwargs.get('fq_band', (0.5, 10))
        validate    = kwargs.get("validate", True) # if True, ssteps will ask for confirmation
        file_name   = kwargs.get("file_name", None)
        out_dir     = kwargs.get("out_dir", seisvo_paths["lte"])

        # defining base params
        ltebase = dict(
            id              = self.stats.code,
            type            = "network",
            stations        = sta_list,
            starttime       = starttime.strftime('%Y-%m-%d %H:%M:%S'),
            endtime         = endtime.strftime('%Y-%m-%d %H:%M:%S'),
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

        ss = SSteps(starttime, endtime, window*60, win_olap=win_olap, subwindow=subwindow*60, subw_olap=subw_olap, validate=validate)

        lwin = int(ss.subwindow*sample_rate)
        ltebase["lswin"] = lwin
        ltebase["nswin"] = int(ss.nsubwin)
        ltebase["nadv"]  = float(ss.subw_adv)
        ltebase["nro_freq_bins"] = len(get_freq(lwin, sample_rate, fq_band=fq_band, pad=pad)[0])
        ltebase["nro_time_bins"] = int(ss.int_nwin*ss.nro_intervals)

        if njobs >= nCPU or njobs == -1:
            njobs = nCPU - 2

        # create hdf5 file and process data
        if not file_name:
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (self.stats.code, starttime.year, starttime.timetuple().tm_yday, endtime.year, endtime.timetuple().tm_yday, window)
        
        if not out_dir:
            out_dir = os.path.join(seisvo_paths["lte"])
        
        file_name_full = os.path.join(out_dir, file_name)
        if file_name_full.split('.')[-1] != 'lte':
            file_name_full += '.lte'

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)
        
        lte = _new_LTE(self, file_name_full, ltebase, njobs, "network")
        return lte


class Array(Network):
    def __init__(self, net_code, sta_code):
        super().__init__(net_code)
        self.sta_code = sta_code

        # change self
        to_remove = [sta.stats for sta in self if sta.stats.code != sta_code]
        for sta_stats in to_remove:
            self.stats.stations.remove(sta_stats)

        # compute UTM position
        self.utm = {}
        for sta in self:
            x, y, n, l = sta.get_latlon(return_utm=True)
            self.utm[sta.stats.location] = {
                "x":x,
                "y":y,
                "zone":(n,l)
            }
        
        self.locs = list(self.utm.keys())


    def get_aperture(self, exclude_locs=[]):
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
                    x1 = self.utm[iloc]["x"]
                    x2 = self.utm[jloc]["x"]
                    y1 = self.utm[iloc]["y"]
                    y2 = self.utm[jloc]["y"]
                    d = np.sqrt(np.abs(x1-x2)**2 + np.abs(y2-y1)**2)
                    
                    if d > a:
                        a = d
        return a

    
    def get_response(self, slow_max, slow_inc, fq_band=(1.,10.), fq_int=0.1, exclude_locs=[]):
        """
        Return array response dictionary
        """
        xutm, yutm = [],[]

        for loc, utm in self.utm.item():
            if loc not in exclude_locs:
                xutm.append(float(utm["x"]))
                yutm.append(float(utm["y"]))

        ans = array_response(xutm, yutm, slow_max, slow_inc, fq_band=fq_band, fq_int=fq_int)

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


    def cc8(self, starttime, endtime, window, overlap, interval=60, slow_max=[3.,0.5], slow_inc=[0.1,0.01], fq_bands=[(1,5)], cc_thres=0.05, filename=None, outdir=None, exclude_locs=[], **kwargs):
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

        slow_inc : list
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

        assert starttime < endtime
        assert window/60 < interval
        assert len(slow_max) == len(slow_inc)

        # load parameters
        njobs       = kwargs.get('njobs', -1)
        toff_sec    = kwargs.get('toff_sec', 10)
        sample_rate = kwargs.get("sample_rate", self[0].stats.sample_rate)
        validate    = kwargs.get("validate", True) # if False, ss do not ask for confirmation

        # defining base params
        if exclude_locs:
            locs = [loc for loc in self.locs if loc not in exclude_locs]
        else:
            locs = self.locs

        cc8base = dict(
            id              = '.'.join([self.stats.code, self.sta_code, self.comp]),
            locs            = locs,
            starttime       = starttime.strftime('%Y-%m-%d %H:%M:%S'),
            endtime         = endtime.strftime('%Y-%m-%d %H:%M:%S'),
            interval        = interval,
            window          = window,
            overlap         = overlap,
            sample_rate     = sample_rate,
            fq_bands        = fq_bands,
            slow_max        = slow_max,
            slow_inc        = slow_inc,
            cc_thres        = cc_thres
        )

        ss = SSteps(starttime, endtime, window, interval=interval, win_olap=overlap)
        cc8base["interval"] = float(interval)
        cc8base["nro_intervals"] = int(ss.nro_intervals)
        cc8base["nro_time_bins"] = int(ss.total_nwin)
        cc8base["nwin"] = int(ss.int_nwin)
        cc8base["lwin"] = int(window*sample_rate)
        cc8base["nadv"] = float(ss.win_adv)
        cc8base["toff_sec"] = toff_sec

        if njobs >= nCPU or njobs == -1:
            njobs = nCPU - 2
        
        if njobs > int(ss.nro_intervals)*len(fq_bands):
            njobs = int(ss.nro_intervals)*len(fq_bands)
            print(f"warn  ::  njobs set to {njobs}")

        if validate:
            print('')
            print(' CC8base INFO')
            print(' -------------')
            print(f'{"       ID     ":^5s} : ', cc8base.get("id"))
            print(f'{"      LOCS    ":^5s} : ', cc8base.get("locs"))
            print(f'{"    Start time":^5s} : ', cc8base.get("starttime"))
            print(f'{"      End time":^5s} : ', cc8base.get("endtime"))
            print(f'{"interval [min]":^5s} : ', cc8base.get("interval"))
            print(f'{"  window [sec]":^5s} : ', cc8base.get("window"))
            print(f'{"   window olap":^5s} : ', cc8base.get("overlap"))
            print(f'{"   nro windows":^5s} : ', cc8base.get("nwin")*cc8base.get("nro_intervals"))
            print(f'{"    |_ per int":^5s} : ', cc8base.get("nwin"))
            print('')
            name = input("\n Please, confirm that you are agree (Y/n): ")
            if name not in ('', 'Y', 'y'):
                return
        
        # define name
        if not filename:
            filename = "%s.%s%03d-%s%03d.cc8" % (cc8base["id"], starttime.year, starttime.timetuple().tm_yday, endtime.year, endtime.timetuple().tm_yday)
        else:
            if filename.split('.')[-1] != "cc8":
                filename += ".cc8"

        if not outdir:
            outdir = os.path.join(seisvo_paths["cc8"])
        
        filenamef = os.path.join(outdir, filename)
        if os.path.isfile(filenamef):
            os.remove(filenamef)
            print(' file %s removed.' % filenamef)
        
        cc8 = CC8.new(self, filenamef, cc8base, njobs)
        return cc8


    def plot(self, save=False, filename="./array_map.png"):
        fig = plot_array(self)

        if save:
            fig.savefig(filename)


    def get_excluded_locs(self, loc_list):
        exclude_loc = []
        for loc in self.locs:
            if loc not in loc_list:
                exclude_loc.append(loc)
        
        return exclude_loc


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

