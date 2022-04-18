
import numpy as np
from datetime import timedelta
from seisvo.file.air import AiR
from seisvo.signal.infrarray import cross_corr
from seisvo.signal import time_bins
from seisvo.utils.dates import tic, toc

class iArray(object):

    def __init__(self, network, sta_code, loc_list=None, model=None):

        self.network = self._set(network, sta_code)
        
        if loc_list:
            for loc in loc_list:
                self.network.remove(sta_code, loc)
        self.id = sta_code

        if not self.network:
            raise ValueError(' no infrasound array in network!')

        self._set_central()

        if model:
            self.model = model
            self.set_model()
        else:
            self.model = False


    def __str__(self):
        return self.network.__str__()


    def _set(self, network, sta_code):
        """
        This code creates the network related to the infrasound array
        """
        chan = 'P'

        sta_to_remove = []
        for sta in network.station:
            if sta.info.code == sta_code:
                for ch in sta.info.chan:
                    if ch[-1] != chan:
                        sta.info.chan.remove(ch)
                if not sta.info.chan:
                    sta_to_remove += [sta]
            else:
                sta_to_remove += [sta]
                
        id_list = [sta.info.id for sta in sta_to_remove]
        for sta in network.station:
            if sta.info.id in id_list:
                network.station.remove(sta)

        if not network.station:
            return False

        else:
            network.stations_info = [sta.info.id for sta in network.station]
            return network


    def get_central_station(self):
        for sta in self.network.station:
            if sta.info.central:
                return sta


    def _set_central(self):
        self.lat0 = None
        self.lon0 = None
        for sta in self.network.station:
            if sta.info.central:
                self.lat0 = sta.get_latlon(degree=False)[0]
                self.lon0 = sta.get_latlon(degree=False)[1]

        if not self.lat0 or not self.lon0:
            raise ValueError('Not central station is defined')


    def set_model(self, model=None):
        """
        This code compute the model.
        :param model: dictionary with radii, vel_air, h_src, src_dgr (tuple of min, max degree)
        """

        if model:
            self.model = model

        if not self.model:
            raise ValueError(' no model defined.')

        # load model
        radii = self.model.get('radii')
        vel_air = self.model.get('vel_air')
        h_src = self.model.get('h_src')
        src_dgr = self.model.get('src_dgr')

        # compute model parameters
        h_mean = np.array([sta.info.elev for sta in self.network.station]).mean()
        h_diff = h_src - h_mean
        th = np.arange(src_dgr[0], src_dgr[1], 1) * np.pi / 180
        self.nro_srcs = len(th)
        x_src = self.lat0 + radii * np.sin(th)
        y_src = self.lon0 + radii * np.cos(th)

        self.dt_times = {}
        for sta in self.network:
            lat_lon = sta.get_latlon(degree=False)
            xx = lat_lon[0] - x_src
            yy = lat_lon[1] - y_src
            zz = sta.info.elev - h_diff
            dist_src = np.sqrt(xx ** 2 + yy ** 2 + zz ** 2)
            self.dt_times[sta.info.id] = dist_src / vel_air

        str_arg = 'zip(%s)' % (','.join(["self.dt_times['%s']" % (sta.info.id) for sta in self.network.station]))
        t_min = [min(xyz) for xyz in eval(str_arg)]

        self.fs = self.network[0].get_sampling_rate()
        for sta in self.network:
            self.dt_times[sta.info.id] = np.array(np.around((self.dt_times[sta.info.id] - t_min) * self.fs), dtype='int')


    def ccorr(self, starttime, endtime, interval=60, time_width=6, overlap=1, fq_band=(0.5, 5)):
        """
        This code compute the cross correlation between stations following the Ripepe algorithm
        :param starttime: datetime object
        :param endtime: optional. datetime object > starttime
        :param interval: interval in min. by default 60 min
        :param time_width: time window in sec. by default 6 sec
        :param overlap: time overlap in sec. by default 1 sec.
        :param prefilt: bandwidth filter. by default 0.5-5 Hz.
        :param arr_file: boolean, if True save onto a hdf5 file
        :return: dictionary
        """

        start_time = starttime
        time_delta = timedelta(minutes=interval)

        if endtime:
            if endtime - starttime < time_delta:
                time_delta = endtime - starttime

            if endtime < starttime:
                raise ValueError('endtime is less than starttime!')

        else:
            endtime = start_time + time_delta

        # check files exist
        missing_files = self.network.check_files(start_time, endtime, self.id)

        if missing_files:
            # create out.txt with the informacion of the dates and stations
            missing_dates = list(set([dt.datetime.strptime(x.split('-')[1], '%Y%j') for x in missing_files]))
            missing_dates.sort()

            one_day = dt.timedelta(days=1)
            starttime_list = [start_time]
            endtime_list = [endtime]
            recursive_days = []
            for n, (i,j) in enumerate(zip(missing_dates[:-1], missing_dates[1:])):
                if j != i+one_day:
                    starttime_list += [i+one_day]
                    endtime_list += [i]
                else:
                    recursive_days += [i]

            if recursive_days:
                starttime_list += [max(recursive_days)+2*one_day]
                endtime_list += [min(recursive_days)]

            starttime_list.sort()
            endtime_list.sort()

            for i,j in zip(starttime_list, endtime_list):
                self.get_ccorr(i,j,interval=interval, time_width=time_width, overlap=overlap,
                    fq_band=fq_band, air_file=True)

        else:
            self.get_ccorr(start_time, endtime, interval=interval, time_width=time_width, overlap=overlap,
                    fq_band=fq_band, air_file=True)


    def get_ccorr(self, starttime, endtime=None, interval=15, time_width=6, overlap=1, fq_band=(0.5, 5), sps=100, air_file=False, counter=None):
        """
        This code compute the cross correlation between stations following the Ripepe algorithm
        :param starttime: datetime object
        :param endtime: optional. datetime object > starttime
        :param interval: interval in min. by default 60 min
        :param time_width: time window in sec. by default 6 sec
        :param overlap: time overlap in sec. by default 1 sec.
        :param prefilt: bandwidth filter. by default 0.5-5 Hz.
        :param arr_file: boolean, if True save onto a hdf5 file
        :return: dictionary
        """

        start_time = starttime
        time_delta = timedelta(minutes=interval)

        if endtime:
            if endtime - starttime < time_delta:
                time_delta = endtime - starttime

            if endtime < starttime:
                raise ValueError('endtime is less than starttime!')

        else:
            endtime = start_time + time_delta

        info = {}
        info['code'] = '%s.%s' % (self.network.info.code, self.id)
        if not sps:
            sps = self.network[0].get_sampling_rate(starttime)
        else:
            info['sampling_rate'] = sps
        info['starttime'] = starttime.strftime('%Y-%m-%d %H:%M:%S')
        info['endtime'] = endtime.strftime('%Y-%m-%d %H:%M:%S')
        info['time_width'] = time_width
        info['overlap'] = overlap
        info['fq_band'] = fq_band
        info['radii'] = self.model['radii']
        info['vel_air'] = self.model['vel_air']
        info['h_src'] = self.model['h_src']
        info['src_dgr'] = self.model['src_dgr']
        info['sensors_loc'] = self.network.stations_info

        interval = int(time_delta.total_seconds())

        if air_file:
            ext = 'air'
            if counter:
                ext += str(counter)
            
            file_name = '%s.%s%03d-%s%03d.%s' % (self.id, starttime.year, starttime.timetuple().tm_yday,
                endtime.year, endtime.timetuple().tm_yday, ext)
            air_file = AiR.new(file_name, info, self.dt_times, interval)

        else:
            out = {}
            out['time'] = []
            out['mcorr'] = None
            out['p_avg'] = None
            out['p_max'] = None
            out['info'] = info
            out['model'] = self.model

        item = 0
        step = 1
        nro_steps = time_bins(starttime, endtime, interval)

        list_elapse_time = []
        while start_time < endtime:
            tic()
            array = self.network.get_stream(start_time, start_time+time_delta, prefilt=fq_band, remove_response=True)

            kwargs = {'air_file':air_file, 'item':item}
            ans = cross_corr(array, self.dt_times, self.nro_srcs, time_width=time_width, overlap=overlap, **kwargs)
            step_pcent = int(100*(step/nro_steps))
            text = ' %d%% ... Time: %s --- ' % (step_pcent,
                                                start_time.strftime('%Y-%m-%d %H:%M:%S'))
            list_elapse_time += [toc(text)]
            print('')

            if not air_file:
                out['time'] += ans[0]

                if out['mcorr'] is None:
                    out['mcorr'] = ans[1]
                    out['p_avg'] = ans[2]
                    out['p_max'] = ans[3]

                else:
                    out['mcorr'] = np.concatenate((out['mcorr'], ans[1]), axis=0)
                    out['p_avg'] = np.concatenate((out['p_avg'], ans[2]), axis=0)
                    out['p_max'] = np.concatenate((out['p_max'], ans[3]), axis=0)

            else:
                item = ans

            # print status

            step += 1
            start_time += time_delta


        if air_file:
            return air_file, list_elapse_time

        else:
            return out, list_elapse_time


    def get_waveform(self, out, time_bin, azm=None, time_pad=30, fq_band=()):

        from seisvo import Stream2

        if not fq_band:
            fq_band = out['info']['fq_band']

        delta = 1/self.network[0].get_sampling_rate()
        time = out['time'][time_bin]
        time_width = timedelta(seconds=out['info']['time_width'])

        if not azm:
            ans = self.get_max_values(out)
            azm_bin = ans[1][time_bin]
        else:
            azm_bin = self.get_azimuth(out, azm)

        st = Stream2()
        for sta in self.network.station:
            if time_pad:
                pad_time = timedelta(seconds=time_pad)
                starttime = time - pad_time
                endtime = time + pad_time

            else:
                sec = float((int(self.dt_times[sta.info.id][azm_bin]) - 1) * delta)
                dt = timedelta(seconds=sec)
                starttime = time + dt
                endtime = time + dt + time_width

            st += sta.get_stream(starttime, endtime, remove_response=True, prefilt=fq_band)

        return st


    def get_starttime(self, out, time_bin, azm=None):

        v_bars = {}
        delta = 1 / self.network[0].get_sampling_rate()
        time = out['time'][time_bin]
        time_width = timedelta(seconds=out['info']['time_width'])

        if not azm:
            ans = self.get_max_values(out)
            azm_bin = ans[1][time_bin]
        else:
            azm_bin = self.get_azimuth(out, azm)

        for sta in self.network.station:
            dt = timedelta(seconds=float((int(self.dt_times[sta.info.id][azm_bin]) - 1) * delta))
            v_bars[sta.info.id] = (time + dt, time + dt + time_width)

        return v_bars


    def plot_waveform(self, out, time_bin, azm=None, sort_dt=False, pad_time=30, prefilt=()):

        fig = plt.figure()
        axspec = fig.add_subplot(511)
        ax0 = fig.add_subplot(512)
        ax1 = fig.add_subplot(513, sharex=ax0)
        ax2 = fig.add_subplot(514, sharex=ax1)
        ax3 = fig.add_subplot(515, sharex=ax2)

        if not azm:
            ans = self.get_max_values(out)
            azm_bin = ans[1][time_bin]
            azm = azm_bin
        else:
            azm_bin = self.get_azimuth(out, azm)

        time = out['time'][time_bin]
        r = out['mcorr'][azm_bin, time_bin]
        p_max = out['p_max'][azm_bin, time_bin]
        p_avg = out['p_avg'][azm_bin, time_bin]

        fig.suptitle('Azimuth: %s  |  R: %0.2f \n P_max: %1.5f  |  P_avg: %1.5f' % (azm, r, p_max, p_avg))

        # sort by arrival time
        if sort_dt:
            station_names = [(self.dt_times[sta.info.id][azm], sta.info.id) for sta in self.network.station]
            station_names.sort()

        st = self.get_waveform(out, time_bin, azm=azm, time_pad=pad_time, prefilt=prefilt)

        v_bars = self.get_starttime(out, time_bin, azm=azm)

        dt = st[0].stats.delta
        time0 = st[0].stats.starttime.datetime
        npts = st[0].stats.npts
        ttime = [time0 + timedelta(seconds=dt * (x + 1)) for x in range(npts)]
        time_pad = timedelta(seconds=pad_time)

        for n, station in enumerate(self.network.station):
            id = station.info.id
            sta_code = station.code
            sta_loc = station.loc
            trace = Trace2(st.select(station=sta_code, location=sta_loc)[0])
            tr = trace.data

            if n == 0:
                trace.specgram(starttime=time-time_pad, endtime=time+time_pad,
                               fq_band=out['info']['fq_band'], axes=axspec, plot=False)
                axspec.xaxis.set_major_formatter(mtick.NullFormatter())

            eval('ax%s' % n).plot(ttime, tr, 'k')
            eval('ax%s' % n).axvspan(v_bars[id][0], v_bars[id][1], alpha=0.5, color='red')
            eval('ax%s' % n).grid(True)
            eval('ax%s' % n).set_ylabel(id)
            eval('ax%s' % n).set_xlim(ttime[0], ttime[-1])

            if n != 3:
                eval('ax%s' % n).xaxis.set_minor_formatter(mtick.NullFormatter())
                eval('ax%s' % n).xaxis.set_major_formatter(mtick.NullFormatter())

        plt.show()


    def plot_sta(self, starttime, delta=5, specgram=True, **kwargs):

        from seisvo.gui.gstationarray import plot_station
        loc_specgram = kwargs.get('loc_specgram', None)
        plot_all_loc = kwargs.get('plot_all_loc', True)
        fq_band = kwargs.get('fq_band', (0.5, 5))
        remove_response = kwargs.get('remove_response', True)
        sampling_rate = kwargs.get('sampling_rate', None)

        sta_list = self.network.station
        loc_sta_list = [s.info.loc for s in sta_list]
        if loc_specgram:
            if not loc_specgram in loc_sta_list:
                if specgram:
                    print('warning: loc not found, spectrogram will be for loc %s' % loc_sta_list[0])
                loc_specgram = loc_sta_list[0]
        else:
            loc_specgram = loc_sta_list[0]

        if not sampling_rate:
            if not fq_band:
                sampling_rate = None
            
            else:
                sampling_rate = 4*fq_band[1]

        title = '%s.%s.[%s]' % (self.network.info.code, sta_list[0].info.code, ','.join(loc_sta_list))

        plot_dict = dict(starttime=starttime, delta=delta, specgram=specgram, sta_list=sta_list,
                         loc_specgram=loc_specgram, plot_all_loc=plot_all_loc, fq_band=fq_band, 
                         remove_response=remove_response, sampling_rate=sampling_rate, title=title)

        plot_station(self, plot_dict)



    def plot_gui(self, **kwargs):
        from seisvo.gui.giarr import gplot
        gplot(iarr=self, **kwargs)


    @staticmethod
    def get_azimuth(out, azm):
        azm_list = list(range(out['model'].get('src_dgr')[0], out['model'].get('src_dgr')[1]))
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


    @staticmethod
    def plot(out, maxSMB=0.5, save=False):

        import matplotlib.pyplot as plt
        import matplotlib.ticker as mtick
        import matplotlib.colors as mcolor

        plt.rc('axes', labelsize=10)
        plt.rc('axes', labelpad=4.0)
        plt.rc('axes', titlepad=6.0)
        plt.rc('axes', titlesize=10)
        plt.rc('xtick', labelsize=10)
        plt.rc('xtick', labelsize=10)
        plt.rc('ytick', labelsize=10)
        plt.rc('ytick', labelsize=10)
        plt.rc('lines', linewidth=0.5)
        plt.rc('lines', linewidth=0.5)

        def azimuthFmtt(x, pos):
            azimuths = range(out['model'].get('src_dgr')[0], out['model'].get('src_dgr')[1])
            try:
                return str(azimuths[int(x)])
            except:
                pass

        ans = iArray.get_max_values(out)
        time_bin = out['time']
        time_bin_max = len(time_bin)


        fig = plt.figure(figsize=(8, 6))
        fig.suptitle('%s   --   %s' % (out['time'][0], out['time'][-1]))

        ax1 = fig.add_axes([0.1, 0.75, 0.8, 0.15])
        ax2 = fig.add_axes([0.1, 0.55, 0.8, 0.15])
        ax3 = fig.add_axes([0.1, 0.35, 0.8, 0.15])
        ax4 = fig.add_axes([0.1, 0.15, 0.8, 0.15])
        cbar_ax = fig.add_axes([0.91, 0.15, 0.025, 0.15])

        # pressure plot
        p_max_r = ans[2][np.where(ans[0] >= maxSMB)]
        t_az = time_bin[np.where(ans[0] >= maxSMB)]
        ax1.scatter(t_az, p_max_r, color='darkred', label=r'$P_{max}$')
        ax1.plot(time_bin, ans[3], color='k', label=r'$P_{avg}$')
        ax1.legend(bbox_to_anchor=(1.1, 1.05))
        ax1.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        ax1.set_ylabel(r'P [Pa]')
        ax1.set_xlim(time_bin[0], time_bin[-1])
        ax1.grid(True)

        # corr plot
        ax2.plot(time_bin, ans[0], color='darkgreen')
        ax2.axhline(y=maxSMB, color='red', linestyle='--')
        ax2.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        ax2.set_ylabel(r'$r_{max}$')
        ax2.set_xlim(time_bin[0], time_bin[-1])
        ax2.set_xticklabels([' ']*len(ax1.get_xticks()))
        ax2.grid(True)

        # azimuth plot
        az_r = ans[1][np.where(ans[0]>=maxSMB)]
        ax3.scatter(t_az, az_r)
        ax3.set_ylabel(r'$azm_{r_{max}}$ [º]')
        ax3.set_xlim(time_bin[0], time_bin[-1])
        ax3.set_ylim(0, 359)
        ax3.yaxis.set_major_locator(mtick.FixedLocator([0,90,180,270]))
        ax3.set_xticklabels([' ']*len(ax1.get_xticks()))
        ax3.grid(True)

        # crosscorrelation plot
        mcorr = out['mcorr']
        mcorr = np.flipud(mcorr.T)
        azmbin = range(mcorr.shape[0])
        halfbin_time = (time_bin[1] - time_bin[0]) / 2.0
        halfbin_azmbin = (azmbin[1] - azmbin[0]) / 2.0
        extent = (mdates.date2num(time_bin[0] - halfbin_time),
                  mdates.date2num(time_bin[-1] + halfbin_time),
                  azmbin[0] - halfbin_azmbin,
                  azmbin[-1] + halfbin_azmbin)
        v_max = mcorr.max()
        v_min = mcorr.min()
        levels = mtick.MaxNLocator(nbins=360).tick_values(v_min, v_max)
        cmap = plt.get_cmap('inferno')
        norm = mcolor.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        im = ax4.imshow(mcorr, cmap=cmap, norm=norm, interpolation='gaussian', extent=extent, aspect='auto')
        ax4.yaxis.set_major_locator(mtick.MaxNLocator(nbins=4))
        ax4.yaxis.set_major_formatter(mtick.FuncFormatter(azimuthFmtt))

        fig.colorbar(im, cax=cbar_ax, ticks=[v_min, (v_min+v_max)/2, v_max], orientation='vertical',
                             format='%.1f')
        lims = ax4.axis('tight')
        ax4.xaxis_date()
        ax1.set_xticklabels([' ']*len(ax1.get_xticks()))
        ax4.set_xlabel('Time')
        ax4.set_ylabel('azm [º]')
        
        if save:
            date_fmt = out['time'][0].strftime('%Y%m%d')
            name = '%s.%s.png' % (out['info']['code'], date_fmt)
            fig.savefig(name, dpi=200)
            del fig
        else:
            plt.show()
