#!/usr/bin/python3
# coding=utf-8

from seisvo import __seisvo__
from seisvo.core import get_network
from seisvo.core.obspyext import Stream2
from seisvo.core.station import Station

from seisvo.signal import SSteps
from seisvo.signal.infrasound import infrasound_model_default, cross_corr

import numpy as np
import collections as col
import datetime as dt
import os

class Network(object):
    def __init__(self, net_code):
        self.info = get_network(net_code)
        self.station = [Station(x) for x in self.info.stations]


    def __str__(self):
        return self.info.__str__()


    def __len__(self):
        return len(self.station)


    def __getitem__(self, item):
        return self.station[item]


    def get_datebound(self, sta_code=None):
        """
        This code get the bound dates of the network
        """

        start = None
        end = None
        
        n = 0
        for sta in self.station:
            if not sta_code or sta.info.code == sta_code:
                st = sta.info.starttime
                et = sta.info.endtime

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
        Get a iArray object. Only for infrasound stations
        :param loc_list: Specify the list of locations to exclude
        """
        return iArray(self.info.code, sta_code, model)


    def get_sta(self, sta_code, loc=''):
        """
        Get a station object
        """
        for sta in self.station:
            if sta_code == sta.info.code and loc ==sta.info.loc:
                return sta
        return None


    def get_stream(self, starttime, endtime, sta_code=None, **kwargs):
        """
        This code get Stream object for the network.
        """
        
        stream = Stream2()
        for sta in self.station:
            if not sta_code or sta_code == sta.info.code:
                try:
                    stream += sta.get_stream(starttime, endtime, **kwargs)
                except:
                    continue
        return stream


    def plot_map(self, zoom_scale=1, arcgis_map='World_Shaded_Relief', epsg=4326, pixel=1500, dpi=100, save=False):
        """
        Plot a map for geograpich network
        :param zoom_scale: how much to zoom from coordinates (in degrees)
        :param map: http://server.arcgisonline.com/arcgis/rest/services
        """
        
        from seisvo.utils.maps import get_map

        # desired coordinates
        coord = [self.info.latOrig, self.info.lonOrig]
        title = "Network: %s" % self.info.code
        
        fig, ax = get_map(coord, arcgis_map, zoom_scale, epsg, pixel=pixel, dpi=dpi, title=title)

        # draw station
        for station in self.station:
            (lat, lon) = station.get_latlon(degree=True)
            ax.scatter(lon, lat, marker='^', label=station.info.id, transform=proj)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        if save:
            fig.savefig('map.png', format='png', dpi=500)
        else:    
            fig.show()


    def check_files(self, startday, endday, sta_code, loc=None, plot=False):
        """
        This code plot the availability of the stations of the network.
        :param startdate: datetime
        :param enddate: datetime
        :param sta_code: string
        :param loc: string, optional
        :param plot: boolean, by default False
        """

        from datetime import timedelta
        
        if plot:
            from seisvo.utils.plotting import plot_check

        sta_list = []
        for sta in self.station:
            if sta_code == sta.info.code:
                if not loc or sta.info.loc == loc:
                    for chan in sta.info.chan:
                        sta_list += ['%s.%s' % (sta.info.id, chan)]

        if not sta_list:
            raise TypeError(' Station list is empty.')

        day_list = [startday + timedelta(days=k) for k in range(0,(endday-startday).days+1)]
        
        if plot:
            availability = np.empty(shape=(len(day_list), len(sta_list)))
        list_missing_days = []

        r = 0
        nro_reads = len(day_list)*len(sta_list)
        for i, x in enumerate(day_list):
            for j, y in enumerate(sta_list):
                for sta in self.station:
                    if sta.info.id == '.'.join(y.split('.')[0:-1]):
                        #print ('  reading data... (%d%%)' % (100*r/nro_reads), end='\r')
                        status = sta.is_file(y.split('.')[-1], date=x)
                        
                        if plot:
                            if isinstance(status, str):
                                status = 1
                            else:
                                status = 0
                            availability[i,j] = status
                        else:
                            if not status:
                                missing_day = x.strftime('%Y%j')
                                missing_day_str = '%s-%s' %(sta.info.id, missing_day)
                                list_missing_days += [missing_day_str]

                        r += 1
            
        if plot:
            title = "Availability for Network %s" % self.info.code
            plot_check(title, sta_list, availability, day_list)

        else:
            return list_missing_days


    def remove(self, sta_code, loc):
        """
        Remove station from the network class
        """
        sta = self.get_sta(sta_code, loc)
        if sta:
            sta_id = sta.id
            self.station.remove(sta)
            self.info.stations_info.remove(sta_id)


class iArray(Network):
    def __init__(self, net_code, sta_code, chan='P', **kwargs):
        super().__init__(net_code)
        self.sta_code = sta_code
        self.central_sta = None
        self.chan = chan
        self.__setstacode__()
        self.__check__()
        self.set_central(sta_code=kwargs.get('central_sta'))
        self.set_model(model=kwargs.get('model', None))


    def __setstacode__(self):
        sta_to_remove = []
        for sta in self.station:
            if sta.info.code == self.sta_code:
                for ch in sta.info.chan:
                    if ch[-1] != self.chan:
                        sta.info.chan.remove(ch)
                if not sta.info.chan:
                    sta_to_remove += [sta]
            else:
                sta_to_remove += [sta]
        
        id_list = [sta.info.id for sta in sta_to_remove]
        for sta in self.station:
            if sta.info.id in id_list:
                self.station.remove(sta)
        
        self.stations_info = [sta.info.id for sta in self.station]
    

    def __check__(self):
        # check sample rate
        sample_rate = list(set([sta.info.sampling_rate for sta in self.station]))
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
                if sta.info.central:
                    self.lat0 = sta.get_latlon(degree=False)[0]
                    self.lon0 = sta.get_latlon(degree=False)[1]
                    self.central_sta = sta.info.id
            else:
                if sta.info.code == sta_code:
                    self.lat0 = sta.get_latlon(degree=False)[0]
                    self.lon0 = sta.get_latlon(degree=False)[1]
                    self.central_sta = sta.info.id
        
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

        h_mean = np.array([sta.info.elev for sta in self.station]).mean()
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
            zz = sta.info.elev - h_diff
            dist_src = np.sqrt(xx ** 2 + yy ** 2 + zz ** 2)
            self.dt_times[sta.info.id] = dist_src / self.model['vel_air']

        dt_central = self.dt_times[self.central_sta]
        for sta in self:
            self.dt_times[sta.info.id] = np.array(np.around((self.dt_times[sta.info.id] - dt_central) * self.sample_rate), dtype='int')
    

    def __add_event__(self, isde, label, starttime, duration, air_file, pmax, pavg, channels):
        """
        Add a new episode in iSDE database, for adding a row visit database/__init__ info
        Check database atributes por kwargs
        """

        event_to_save = {}
        event_to_save['network'] = self.info.code
        event_to_save['station'] = self.sta_code
        event_to_save['label'] = label
        event_to_save['channels'] = channels
        event_to_save['starttime'] = starttime #datetime
        event_to_save['duration'] = duration
        event_to_save['air_file'] = air_file
        event_to_save['pmax'] = pmax
        event_to_save['pavg'] = pavg

        isde.add_row(event_to_save)


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
                azm_delta = dt.timedelta(seconds=int(self.dt_times[sta.info.id][azm_idx]) * delta)
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

        info = {}
        info['code'] = '%s.%s' % (self.info.code, self.sta_code)
        info['sampling_rate'] = self.sample_rate_
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
                remove_response=True
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
                sec = float((int(self.dt_times[sta.info.id][azm_bin]) - 1) * delta)
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
            dt_delta = dt.timedelta(seconds=float((int(self.dt_times[sta.info.id][azm_bin]) - 1) * delta))
            v_bars[sta.info.id] = (time + dt_delta, time + dt_delta + time_width)

        return v_bars
    

    def plot_gui(self, **kwargs):
        from seisvo.gui.giarr import gplot
        gplot(array=self, **kwargs)


    @staticmethod
    def get_azimuth(out, azm):
        azm_list = list(range(out['model'].get('src_dgr')[0], out['model'].get('src_dgr')[1]+1))
        return  azm_list.index(azm)


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




