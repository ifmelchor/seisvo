#!/usr/bin/python3
# coding=utf-8

import os
import numpy as np
import collections as col
import datetime as dt
import multiprocessing as mp

from seisvo import LTE_PATH, DB_PATH, CC8_PATH
from seisvo.core import get_network
from .obspyext import Stream2
from .station import Station
from seisvo.sap import CC8, Arfr, plot_array
from seisvo.signal import SSteps, get_freq
from seisvo.file.air import AiR
from seisvo.signal.infrasound import infrasound_model_default, cross_corr
from seisvo.lte import netLTE


class Network(object):
    def __init__(self, net_code):
        self.stats = get_network(net_code)
        self.station = [Station(x) for x in self.stats.stations]


    def __str__(self):
        return self.stats.__str__()


    def __len__(self):
        return len(self.station)


    def __getitem__(self, item):
        if isinstance(item, int):
            return self.station[item]
        
        if isinstance(item, str):
            sta_code = item.split('.')[0]
            
            try:
                loc = item.split('.')[1]
            except IndexError:
                loc = ''
            
            return self.get_sta(sta_code, loc=loc)
    

    def __setstacode__(self, sta_code, component):
        sta_to_remove = []
        
        for sta in self.station:
            if sta.stats.code == sta_code:
                chan_to_remove = []
                
                for ch in sta.stats.chan:
                    if ch[-1] != component:
                        chan_to_remove.append(ch)
                
                if chan_to_remove:
                    for ch in chan_to_remove:
                        sta.stats.chan.remove(ch)

                if not sta.stats.chan:
                    sta_to_remove += [sta]
           
            else:
                sta_to_remove += [sta]
        
        for sta in sta_to_remove:
            self.station.remove(sta)
        
        self.stations_info = [sta.stats.id for sta in self.station]


    def get_datebound(self, sta_code=None):
        """
        This code get the bound dates of the network
        """

        start = None
        end = None
        
        n = 0
        for sta in self.station:
            if not sta_code or sta.stats.code == sta_code:
                st = sta.stats.starttime
                et = sta.stats.endtime

                if n >= 1:
                    if st < start:
                        start = st

                    if et > end:
                        end = et
                else:
                    start = st
                    end = et

                n += 1

        return (start, end)


    def get_iarray(self, sta_code, model=None):
        """
        Get a iArray object.
        """
        return iArray(self.stats.code, sta_code, model)
    

    def get_sarray(self, sta_code, comp="Z", **kwargs):
        """
        Get a sArray object.
        """
        return sArray(self.stats.code, sta_code, comp=comp, **kwargs)


    def get_sta(self, sta_code, loc=''):
        """
        Get a station object
        """
        for sta in self.station:
            if sta_code == sta.stats.code and loc ==sta.stats.loc:
                return sta
        return None


    def get_stream(self, starttime, endtime, sta_code=None, **kwargs):
        """
        This code get Stream object for the network.
        """
        
        stream = Stream2()
        for sta in self.station:
            if not sta_code or sta_code == sta.stats.code:
                try:
                    stream += sta.get_stream(starttime, endtime, **kwargs)
                except:
                    continue
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
        for sta in self.station:
            if sta_code == sta.stats.code:
                if loc:
                    if isinstance(loc, str):
                        if sta.stats.loc == loc:
                            for chan in sta.stats.chan:
                                sta_list += ['%s.%s' % (sta.stats.id, chan)]
                    
                    if isinstance(loc, (list, tuple)):
                        for l in loc:
                            if sta.stats.loc == l:
                                for chan in sta.stats.chan:
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

                
        #         for sta in self.station:
        #             if sta.stats.id == '.'.join(y.split('.')[0:-1]):
        #                 #print ('  reading data... (%d%%)' % (100*r/nro_reads), end='\r')
        #                 status = sta.__read_file__(y.split('.')[-1], date=day)
                        
        #                 if plot:
        #                     if isinstance(status, str):
        #                         status = 1
        #                     else:
        #                         status = 0
        #                     availability[i,j] = status
        #                 else:
        #                     if not status:
        #                         missing_day = x.strftime('%Y%j')
        #                         missing_day_str = '%s-%s' %(sta.stats.id, missing_day)
        #                         list_missing_days += [missing_day_str]

        #                 r += 1
            
        # if plot:
        #     title = "Availability for Network %s" % self.stats.code
        #     plot_check(title, sta_list, availability, day_list)

        # else:
        #     return list_missing_days


    def remove(self, sta_code, loc):
        """
        Remove station from the network class
        """
        sta = self.get_sta(sta_code, loc)
        if sta:
            sta_id = sta.id
            self.station.remove(sta)
            self.stats.stations_info.remove(sta_id)


    def gui(self, starttime, station_list, component="Z", delta=30, sde_file=None, **kwargs):

        from seisvo.plotting.gui.gnetwork import init_network_gui

        # check station_list
        true_station_list = []
        for sta in self.station:
            sta_id = '.'.join([sta.stats.code, sta.stats.loc])
            if sta_id in station_list:
                true_station_list += [sta_id]

        if not true_station_list:
            print(" no stations loaded!. Revise network file!")
            return
        
        if component not in ["Z", "N", "E"]:
            print(" componente must be Z, N, or E")
            return
        
        if sde_file:
            if sde_file.split('.')[-1] != '.db':
                sde_file += '.db'

        else:
            sde_file = os.path.join(DB_PATH, self.stats.code + '.db')
        
        if isinstance(kwargs.get("specgram"), int):
            specgram = station_list[kwargs.get("specgram")]
            kwargs["specgram"] = specgram
        
        if specgram not in true_station_list:
            kwargs["specgram"] = None

        init_network_gui(self, true_station_list, starttime, delta, component, sde_file, **kwargs)


    def lte(self, starttime, endtime, sta_list, window, subwindow, win_olap=0.5, subw_olap=0.75, **kwargs):
        """ Compute (station) LTE file

        Parameters
        ----------
        starttime : datetime
        
        endtime : datetime
        
        sta_list : list of strings
            list of stations id

        window : int [hr]
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
        assert window*60 > subwindow

        # check sta_list
        sta_ob_list = []
        true_sta_list = []
        for sta in sta_list:
            sta_sp = sta.split(".")
            sta_name = sta_sp[0]
            if len(sta_sp) > 1:
                sta_loc  = sta_sp[1]
            else:
                sta_loc  = ""
            
            sta_ob = self.get_sta(sta_name, loc=sta_loc)
            
            if sta_ob.starttime <= starttime and sta_ob.endtime >= endtime:
                sta_ob_list.append(sta_ob)
                true_sta_list.append(sta)
            else:
                print("warn:: station {sta} out of temporal bounds")
        
        if not sta_ob_list:
            print("error:: no data to proceed")
            return None
        
        # load kwargs
        sample_rate = kwargs.get('sample_rate', 40)
        njobs       = kwargs.get('njobs', 1)
        pad         = kwargs.get("pad", 1.0)
        time_bandw  = kwargs.get('time_bandwidth', 3.5)
        fq_band     = kwargs.get('fq_band', (0.5, 10))
        validate    = kwargs.get("validate", True) # if True, ssteps will ask for confirmation
        file_name   = kwargs.get("file_name", None)
        out_dir     = kwargs.get("out_dir", "./")
        resp_dict   = kwargs.get("resp_dict", {})

        # check resp dict
        if resp_dict:
            for key in list(resp_dict.keys()):
                if key not in true_sta_list:
                    print("warn:: {key} in resp_dict not valid")
                    del resp_dict[key]
            rm_sens = True
        else:
            rm_sens = False

        # defining base params
        ltebase = dict(
            id              = self.stats.code,
            type            = "network",
            stations        = true_sta_list,
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

        ss = SSteps(starttime, endtime, -1, window*3600, win_olap=win_olap, subwindow=subwindow*60, subw_olap=subw_olap, validate=validate)

        lwin = int(ss.subwindow*sample_rate)
        ltebase["lswin"] = lwin
        ltebase["nswin"] = int(ss.nsubwin)
        ltebase["nadv"]  = float(ss.subw_adv)
        ltebase["nro_freq_bins"] = len(get_freq(lwin, sample_rate, fq_band=fq_band, pad=pad)[0])
        ltebase["nro_time_bins"] = int(ss.int_nwin*ss.nro_intervals)

        if njobs >= mp.cpu_count() or njobs == -1:
            njobs = mp.cpu_count() - 2

        # create hdf5 file and process data
        if not file_name:
            file_name = '%s.%s%03d-%s%03d_%s.lte' % (self.stats.code, starttime.year, starttime.timetuple().tm_yday, endtime.year, endtime.timetuple().tm_yday, window)
        
        if not out_dir:
            out_dir = os.path.join(LTE_PATH)
        
        file_name_full = os.path.join(out_dir, file_name)
        if file_name_full.split('.')[-1] != 'lte':
            file_name_full += '.lte'

        if os.path.isfile(file_name_full):
            os.remove(file_name_full)
            print(' file %s removed.' % file_name_full)
        
        lte = netLTE.new(sta_ob_list, file_name_full, ltebase, njobs, resp_dict)
        return lte


class sArray(Network):
    def __init__(self, net_code, sta_code, comp='Z', locs=None):
        super().__init__(net_code)
        self.sta_code = sta_code
        self.comp = comp
        self.__setstacode__(sta_code, comp)

        if locs:
            if isinstance(locs, (str, tuple)):
                locs = list(locs)
            self.__filterlocs__(locs)

        self.locs = [sta.stats.loc for sta in self.station]
        self.xUTM = np.array([float(sta.stats.lon) for sta in self.station])
        self.yUTM = np.array([float(sta.stats.lat) for sta in self.station])
        self.resp = Arfr(self)
    
    
    def __filterlocs__(self, loc_list):

        locs_to_remove = []
        for sta in self.station:
            if sta.stats.loc not in loc_list:
                locs_to_remove.append(sta)
        
        for sta in locs_to_remove:
            self.station.remove(sta)
        
        self.stations_info = [sta.stats.id for sta in self.station]


    def get_aperture(self):
        # only for utm data
        a = 0
        for i in range(len(self.locs)):
            for j in range(i,len(self.locs)):
                if i != j:
                    x1 = self.xUTM[i]
                    x2 = self.xUTM[j]
                    y1 = self.yUTM[i]
                    y2 = self.yUTM[j]
                    d = np.sqrt(np.abs(x1-x2)**2 + np.abs(y2-y1)**2)
                    if d > a:
                        a = d
        return a


    def get_mdata(self, start_time, end_time, toff_sec=0, sample_rate=None, fq_band=(), return_stats=False):
        if toff_sec > 0:
            toff = dt.timedelta(seconds=toff_sec)
            start_time -= toff
            end_time += toff

        st = self.get_stream(start_time, end_time, sample_rate=sample_rate)
        
        if st:
            fs = st[0].stats.sampling_rate
            lwin = int((end_time-start_time).total_seconds()*fs) + 1
            mdata = np.empty((len(st), lwin))
            stats = []

            n = 0
            for tr in st:
                tr_dat = tr.get_data(detrend=True, fq_band=fq_band)
                
                if len(tr_dat) == lwin:
                    mdata[n,:] = tr_dat
                    if return_stats:
                        sta = self.get_sta(tr.stats.station, loc=tr.stats.location)
                        stats.append(sta.stats)
                    n += 1

            if np.isnan(mdata).any():
                print("Warning: data containing NaN values")

        else:
            mdata = None
            stats = []

        if return_stats:
            return mdata, stats

        else:
            return mdata


    def cc8(self, starttime, endtime, window, overlap, slow_max=[3.,0.5], slow_inc=[0.1,0.01], fq_bands=[(1.,5.)], cc_thres=0.05, filename=None, outdir=None, interval=1, **kwargs):
        """ Compute CC8 file

        Parameters
        ----------
        starttime : datetime
        
        endtime : datetime
        
        window : int [sec]
            length of the time window (in sec)
        
        overlap : float [0--1)
            overlap percent
        
        interval : int [min]
            length of time window (in min) to store data in memory. By default is 15 min.
        
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
        sample_rate = kwargs.get("sample_rate", 40)
        validate    = kwargs.get("validate", True) # if False, ss do not ask for confirmation

        # defining base params
        cc8base = dict(
            id              = '.'.join([self.stats.code, self.sta_code, self.comp]),
            locs            = self.locs,
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

        ss = SSteps(starttime, endtime, interval, window, validate=validate, win_olap=overlap)

        if ss.nwin_eff < ss.total_nwin:
            diff = ss.total_nwin - ss.nwin_eff
            nwext = round(diff/ss.nro_intervals)
            if toff_sec <= nwext*window:
                toff_sec += int(nwext*window)
        else:
            nwext = 0

        cc8base["nwin"] = int(ss.int_nwin)+nwext
        cc8base["lwin"] = int(window*sample_rate)
        cc8base["nadv"] = float(ss.win_adv)
        cc8base["toff_sec"] = toff_sec
        cc8base["nro_intervals"] = int(ss.nro_intervals)
        cc8base["nro_time_bins"] = int(cc8base["nwin"]*cc8base["nro_intervals"])

        if njobs >= mp.cpu_count() or njobs == -1:
            njobs = mp.cpu_count() - 2
        
        if njobs > int(ss.nro_intervals)*len(fq_bands):
            njobs = int(ss.nro_intervals)*len(fq_bands)

            print(f"warn  ::  njobs set to {njobs}")

        if validate:
            self.print_cc8_validation(cc8base)
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
            outdir = os.path.join(CC8_PATH)
        
        filenamef = os.path.join(outdir, filename)

        if os.path.isfile(filenamef):
            os.remove(filenamef)
            print(' file %s removed.' % filenamef)
        
        if CC8.__check_process__(self, starttime, endtime, verbose=True):
            cc8 = CC8.new(self, filenamef, cc8base, njobs)
        
        return 


    def plot(self, save=False, filename="./array_map.png"):
        fig = plot_array(self)

        if save:
            fig.savefig(filename)

    @staticmethod
    def print_cc8_validation(cc8dict):
        print('')
        print(' CC8base INFO')
        print(' -------------')
        print(f'{"      ID     ":^5s} : ', cc8dict.get("id"))
        print(f'{"     LOCS    ":^5s} : ', cc8dict.get("locs"))
        print(f'{"   Start time":^5s} : ', cc8dict.get("starttime"))
        print(f'{"     End time":^5s} : ', cc8dict.get("endtime"))
        print(f'{" window [sec]":^5s} : ', cc8dict.get("window"))
        print(f'{"  window olap":^5s} : ', cc8dict.get("overlap"))
        print(f'{"Total windows":^5s} : ', cc8dict.get("nwin")*cc8dict.get("nro_intervals"))
        print(f'{"     interval":^5s} : ', cc8dict.get("interval"))
        print(f'{"     toff_sec":^5s} : ', cc8dict.get("toff_sec"))
        print(f'{" nwin per int":^5s} : ', cc8dict.get("nwin"))
        print('')


class iArray(Network):
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
        sample_rate = list(set([sta.stats.sampling_rate for sta in self.station]))
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

