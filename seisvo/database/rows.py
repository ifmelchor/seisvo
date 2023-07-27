
import os
from scipy import signal as ss
import numpy as np
import datetime as dt
import sqlalchemy as sql
from ..network import Network
from ..sap import CC8
from ..plotting.array import window_wvfm

def get_sderow(SQLbase):

    class SDErow(SQLbase):
        __tablename__ = 'SDE'
        __table_args__ = {'extend_existing': True}

        # base attributes
        id = sql.Column(sql.Integer, primary_key=True)
        
        event_id   = sql.Column(sql.Integer, nullable=False)
        event_type = sql.Column(sql.String, nullable=False)
        
        network  = sql.Column(sql.String, nullable=False)
        station  = sql.Column(sql.String, nullable=False)
        location = sql.Column(sql.String, nullable=True)
        
        label     = sql.Column(sql.String, nullable=False)
        sublabel  = sql.Column(sql.String, nullable=True)
        starttime = sql.Column(sql.DateTime(timezone=False), nullable=False)
        duration  = sql.Column(sql.Float, nullable=False)

        # event attributes
        time_P   = sql.Column(sql.DateTime(timezone=False), nullable=True)
        weight_P = sql.Column(sql.Integer, nullable=True)
        onset_P  = sql.Column(sql.String, nullable=True)
        time_S   = sql.Column(sql.DateTime(timezone=False), nullable=True)
        weight_S = sql.Column(sql.Integer, nullable=True)
        time_F   = sql.Column(sql.DateTime(timezone=False), nullable=True)    

        # additional floats
        value_1 = sql.Column(sql.Float, nullable=True) # peak to peak
        value_2 = sql.Column(sql.Float, nullable=True) # max amplitude
        value_3 = sql.Column(sql.Float, nullable=True)
        value_4 = sql.Column(sql.Float, nullable=True)
        value_5 = sql.Column(sql.Float, nullable=True)

        # additional strings
        string_1 = sql.Column(sql.String, nullable=True)
        string_2 = sql.Column(sql.String, nullable=True)
        string_3 = sql.Column(sql.String, nullable=True)
        string_4 = sql.Column(sql.String, nullable=True)
        string_5 = sql.Column(sql.String, nullable=True)


        def __str__(self):
            text_info = " EID %i || ID %i | %s | %s.%s.%s " % (self.event_id,\
                self.id, self.event_type, self.network, self.station, self.location)
            return text_info
        

        def get_station_id(self):
            return '.'.join([self.network, self.station, self.location])


        def get_network(self):
            return Network(self.network)


        def get_station(self):
            return self.get_network().get_sta(self.station, self.location)
        

        def get_stream(self, off_seconds=0, **kwargs):
            sta = self.get_station()
            start = self.starttime - dt.timedelta(seconds=off_seconds)
            end   = self.starttime + dt.timedelta(seconds=self.duration+off_seconds)
            st = sta.get_stream(start, end, **kwargs)
            return st
        

        def get_endtime(self):
            return self.starttime + dt.timedelta(seconds=self.duration)

    return SDErow


def get_lderow(SQLbase):

    class LDErow(SQLbase):
        __tablename__ = 'LDE'
        __table_args__ = {'extend_existing': True}

        id = sql.Column(sql.Integer, primary_key=True)
        network = sql.Column(sql.String, nullable=False)
        station = sql.Column(sql.String, nullable=False)
        location = sql.Column(sql.String, nullable=True)
        label     = sql.Column(sql.String, nullable=False)
        starttime = sql.Column(sql.DateTime(timezone=False), nullable=False)
        duration  = sql.Column(sql.Float, nullable=False)
        lte_file  = sql.Column(sql.String, nullable=False) # lte file containing the event info
        # additional floats
        value_1 = sql.Column(sql.Float, nullable=True)
        value_2 = sql.Column(sql.Float, nullable=True)
        value_3 = sql.Column(sql.Float, nullable=True)
        value_4 = sql.Column(sql.Float, nullable=True)
        value_5 = sql.Column(sql.Float, nullable=True)
        # additional strings
        string_1 = sql.Column(sql.String, nullable=True)
        string_2 = sql.Column(sql.String, nullable=True)
        string_3 = sql.Column(sql.String, nullable=True)
        string_4 = sql.Column(sql.String, nullable=True)
        string_5 = sql.Column(sql.String, nullable=True)

        @property
        def endtime(self):
            return self.starttime + dt.timedelta(minutes=self.duration)
        
        @property
        def station_id(self):
            return '.'.join([self.network, self.station, self.location])

        def get_network(self):
            return Network(self.network)
        
        def get_station(self):
            net = self.get_network()
            return net.get_sta(self.station, self.location)

    return LDErow


def get_ccerow(SQLbase):

    class CCErow(SQLbase):
        __tablename__ = 'CCE'
        __table_args__ = {'extend_existing': True}

        id       = sql.Column(sql.Integer, primary_key=True)
        network  = sql.Column(sql.String, nullable=False)
        station  = sql.Column(sql.String, nullable=False)
        label    = sql.Column(sql.String, nullable=True)
        time     = sql.Column(sql.DateTime(timezone=False), nullable=False)
        slow     = sql.Column(sql.Float, nullable=False)
        baz      = sql.Column(sql.Float, nullable=False)
        maac     = sql.Column(sql.Float, nullable=False)
        rms      = sql.Column(sql.Float, nullable=False)
        cc8_file = sql.Column(sql.String, nullable=False)
        fqidx    = sql.Column(sql.String, nullable=False)
        
        # additional floats
        value_1 = sql.Column(sql.Float, nullable=True)
        value_2 = sql.Column(sql.Float, nullable=True)
        value_3 = sql.Column(sql.Float, nullable=True)
        value_4 = sql.Column(sql.Float, nullable=True)
        value_5 = sql.Column(sql.Float, nullable=True)
        
        # additional strings
        string_1 = sql.Column(sql.String, nullable=True)
        string_2 = sql.Column(sql.String, nullable=True)
        string_3 = sql.Column(sql.String, nullable=True)
        string_4 = sql.Column(sql.String, nullable=True)
        string_5 = sql.Column(sql.String, nullable=True)


        @property
        def array_id(self):
            return '.'.join([self.network, self.station])

        def get_cc8(self, path_to_cc8file="./"):
            cc8f = os.path.join(path_to_cc8file, self.cc8_file)
            return CC8(cc8f)


        def get_windowtimes(self, window, off_sec=0):
            halfw = dt.timedelta(seconds=float(window/2))
            start = self.time - halfw - dt.timedelta(seconds=off_sec)
            end   = self.time + halfw + dt.timedelta(seconds=off_sec)
            return (start, end)


        def get_network(self):
            return Network(self.network)

        
        def get_array(self):
            net = self.get_network()
            return net.get_array(self.station)


        def get_beamform(self, taper=True, off_sec=2, path_to_cc8file="./", exclude_locs=[], fq_band=[], plot=False, **fig_kwargs):

            cc8 = self.get_cc8(path_to_cc8file=path_to_cc8file)

            sloint = cc8.stats.slow_int
            slomax = cc8.stats.slow_max
            fs     = cc8.stats.sample_rate
            if not fq_band:
                fq_band = cc8.stats.fq_bands[int(self.fqidx)-1]

            arr, ex_locs = cc8.stats.get_array()
            exclude_locs += ex_locs
            deltas, _ = arr.get_deltatimes(self.slow, self.baz, slomax, sloint, fs, exclude_locs=exclude_locs)

            starttime, endtime = self.get_windowtimes(cc8.stats.window, off_sec=off_sec)
            stream    = arr.get_stream(starttime, endtime, prefilt=fq_band, toff_sec=600, exclude_locs=exclude_locs)

            # shift stream
            wvfm_dict = {}
            interval  = endtime - starttime
            for delta, tr in zip(deltas, stream):
                of_npts = int(600*fs)
                d_npts  = int(delta*fs)
                data    = tr.get_data()
                data_sh = data[of_npts+d_npts:-of_npts+d_npts]
                wvfm_dict[tr.stats.location] = data_sh

            duration = interval.total_seconds()
            time = np.linspace(0, duration, len(data_sh))

            if plot:
                fig_kwargs["title"] = f"EID : {self.id} TIME: {self.time.strftime('%Y %b %m %H:%M:%S')} \n BAZ [{self.baz:.2f}] SLOW [{self.slow:.2f}]"
                fig = window_wvfm(wvfm_dict, time, None, None, **fig_kwargs)
                
                return fig
            
            else:
                
                bf = np.empty((len(time), len(wvfm_dict)))
                for n, (_, wvfm) in enumerate(wvfm_dict.items()):
                    bf[:, n] = wvfm
                suma = np.sum(bf, axis=1) / bf.shape[1]

                if taper:
                    suma *= ss.windows.hann(len(suma))

                return time, wvfm_dict, suma

        
    

    return CCErow