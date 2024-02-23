
import os
from scipy import signal as ss
import numpy as np
import datetime as dt
import sqlalchemy as sql
from ..sap import CC8
from ..signal import get_PSD

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
        value_6 = sql.Column(sql.Float, nullable=True)
        value_7 = sql.Column(sql.Float, nullable=True)
        value_8 = sql.Column(sql.Float, nullable=True)
        value_9 = sql.Column(sql.Float, nullable=True)
        value_10 = sql.Column(sql.Float, nullable=True)
        value_11 = sql.Column(sql.Float, nullable=True)
        value_12 = sql.Column(sql.Float, nullable=True)
        value_13 = sql.Column(sql.Float, nullable=True)
        value_14 = sql.Column(sql.Float, nullable=True)
        value_15 = sql.Column(sql.Float, nullable=True)


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
            import seisvo as sv
            return sv.Network(self.network)


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
        slow_u   = sql.Column(sql.Float, nullable=False)
        baz      = sql.Column(sql.Float, nullable=False)
        baz_u    = sql.Column(sql.Float, nullable=False)
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
        
        @property
        def error(self):
            return np.sqrt(self.baz_u*self.baz_u + self.slow_u*self.slow_u)


        def get_cc8(self, path_to_cc8file="./"):
            cc8f = os.path.join(path_to_cc8file, self.cc8_file)
            return CC8(cc8f)


        def get_array(self):
            import seisvo as sv
            return sv.Array(self.network, self.station)


        def slowmap(self, slomax=None, sloint=None, path_to_cc8file="./", exclude_locs=[], fq_band=[], tol=1e-4, **fig_kwargs):

            cc8 = self.get_cc8(path_to_cc8file=path_to_cc8file)

            if not sloint:
                sloint = cc8.stats.slow_int

            if not slomax:
                slomax = cc8.stats.slow_max

            if not fq_band:
                fq_band = cc8.stats.fq_bands[int(self.fqidx)-1]

            arr, ex_locs = cc8.stats.get_array()
            exclude_locs += ex_locs
            
            slowarg = {
                "slomax":slomax,
                "sloint":sloint,
                "exclude_locs":exclude_locs,
                "fq_band":fq_band,
            }
            
            # slowarg["slow0"], _ = arr.deltatimes(self.slow, self.baz, slowarg=slowarg, tol=tol, return_xy=True)
            
            halfw     = dt.timedelta(seconds=cc8.stats.window/2)
            starttime = self.time - halfw

            fig = arr.slowmap(starttime, cc8.stats.window, slowarg=slowarg, plot=True, show_title=True, **fig_kwargs)

            return fig


        def beamform(self, taper=True, off_sec=5, slomax=None, sloint=None, path_to_cc8file="./", exclude_locs=[], fq_band=[], plot=False, **fig_kwargs):

            cc8 = self.get_cc8(path_to_cc8file=path_to_cc8file)

            if not sloint:
                sloint = cc8.stats.slow_int
            
            if not slomax:
                slomax = cc8.stats.slow_max
        
            if not fq_band:
                fq_band = cc8.stats.fq_bands[int(self.fqidx)-1]

            arr, ex_locs = cc8.stats.get_array()
            exclude_locs += ex_locs

            slowarg = {
                "slomax":slomax,
                "sloint":sloint,
                "exclude_locs":exclude_locs,
                "fq_band":fq_band
            }

            halfw     = dt.timedelta(seconds=cc8.stats.window/2)
            starttime = self.time - halfw
            endtime   = starttime + dt.timedelta(seconds=cc8.stats.window)
            
            if off_sec > 0:
                starttime -= dt.timedelta(seconds=off_sec)
                endtime   += dt.timedelta(seconds=off_sec)
                duration = (endtime - starttime).total_seconds()
                startw   = duration/2 - cc8.stats.window/2
                endw     = startw + cc8.stats.window
                shadow_times = (startw, endw)
            else:
                shadow_times = ()

            ans = arr.beamform(starttime, endtime, self.slow, self.baz, slowarg=slowarg,\
                shadow_times=shadow_times, taper=taper, plot=plot, **fig_kwargs)

            return ans
            

        def get_psd(self, fs, slomax=None, sloint=None, path_to_cc8file="./", exclude_locs=[]):

            cc8 = self.get_cc8(path_to_cc8file=path_to_cc8file)
            fq_band = cc8.stats.fq_bands[int(self.fqidx)-1]
            
            _, bmf, time = self.beamform(taper=True, off_sec=2, slomax=slomax, sloint=sloint, path_to_cc8file=path_to_cc8file, exclude_locs=exclude_locs, fq_band=fq_band, plot=False)

            return get_PSD(bmf, fs, fq_band=fq_band)



        def plot(self, off_sec=7, taper=True, path_to_cc8file="./", exclude_locs=[], fq_band=[]):

            cc8 = self.get_cc8(path_to_cc8file=path_to_cc8file)
            arr, ex_locs = cc8.stats.get_array()
            exclude_locs += ex_locs

            if not fq_band:
                fq_band = cc8.stats.fq_bands[int(self.fqidx)-1]

            slowarg = {
                "slomax":cc8.stats.slow_max,
                "sloint":cc8.stats.slow_int,
                "exclude_locs":exclude_locs,
                "fq_band":fq_band
            }

            halfw     = dt.timedelta(seconds=cc8.stats.window/2)
            starttime = self.time - halfw

            fig = arr.plot(starttime, cc8.stats.window, offsec=off_sec, taper=taper, slowarg=slowarg)

            return fig


    return CCErow