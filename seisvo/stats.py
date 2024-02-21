#!/usr/bin/python3
# coding=utf-8


import os
import utm
import pickle
import h5py
import datetime as dt
from glob import glob
from obspy import UTCDateTime
from obspy.core.util.attribdict import AttribDict
from obspy.core.inventory.response import Response
from .obspyext import read2

class _Stats(AttribDict):
    def __init__(self, header):
        super().__init__(header)
        self.keys_ = list(header.keys())


    def __add_attr__(self, listin):
        self.attributes = listin
        

class StationStats(_Stats):
    def __init__(self, network, header):
        super().__init__(header)
        self.network = network
        
        # load location
        if "location" in self.keys_:
            if not self.location:
                self.location = ""
        
        # set id and dates
        self.id = '.'.join([network.code, self.code, self.location])
        self.__set_dates__()


    def __str__(self):
        priorized_keys = [
            'id',
            'code',
            'location',
            'channels',
            'starttime',
            'endtime',
            'sample_rate'
        ]
        return self._pretty_str(priorized_keys)
        

    def check_channel(self, channel=None):
        """
        This function return a list of available channels
        >> channel can be a list or a string
        """

        if isinstance(channel, type(None)):
            true_chan = self.channels

        elif isinstance(channel, str):
            if channel not in self.channels:
                print(' [warn] channel %s not available' % channel)
                true_chan = None
            else:
                true_chan = [channel]
            
        else:
            # channel is a list/tuple/array
            true_chan = [ch for ch in channel if ch in self.channels]

        return true_chan


    def get_latlon(self, **kwargs):
        """
        Get longitude coord. in 'degree' or 'utm'.
        in 'degree' return (lat,lon)
        in 'utm' return (eastern, northern, nro, zone)
        """

        if 'lat' not in self.keys_ or 'lon' not in self.keys_:
            print(" lat/lon not defined in network JSON file")
            return None
        
        lat = self.lat
        lon = self.lon

        if not lat or not lon:
            print(f" lat/lon not defined in network JSON file for station {self.id}")
            return None

        return_utm = kwargs.get("utm", False)
        in_km = kwargs.get("in_km", False)
        if return_utm:
            (x, y, code, nro) = utm.from_latlon(lat, lon)
            if in_km:
                x /= 1000
                y /= 1000
            return (x,y,code,nro)
        
        else:
            return (lat, lon)


    def get_response(self):
        if "resp_file_path" in self.keys_ :
            filepath = os.path.dirname(self.network.file)

            if isinstance(self.resp_file_path, str):
                resp_file = os.path.join(filepath, self.resp_file_path)
                if os.path.isfile(resp_file):
                    with open(resp_file, 'rb') as f:
                        resp = pickle.load(f)
                    
                    if isinstance(resp, Response): # check if resp is a NRL object
                        return resp
                    else:
                        print(f" resp object [type: {type(resp)} is not a Response object!")
                
            if isinstance(self.resp_file_path, dict):
                resp = {}
                for chan in self.channels:
                    if chan in list(self.resp_file_path.keys()):
                        resp_file = os.path.join(filepath, self.resp_file_path[chan])
                        if os.path.isfile(resp_file):
                            with open(resp_file, 'rb') as f:
                                resp = pickle.load(f)
                            
                            if isinstance(resp, Response):
                                resp[chan] = resp
                            else:
                                print(f" resp object [type: {type(resp)} is not a Response object!")
                            
        return None
    

    def get_factor(self, channel=None):
        if "resp_factor" in self.keys_ :

            if isinstance(self.resp_factor, float):
                return self.resp_factor
            
            if isinstance(self.resp_factor, dict):
                resp_fact = {}
                for chan in self.channels:
                    if chan in list(self.resp_factor.keys()):
                        resp_fact[chan] = self.resp_factor[chan]
                
                return resp_fact

        return None


    def __read_file__(self, chan, julian_date=(), date=None, stream=False, **readkwargs):
        """
        Read file if exists for a specific date of the station's channel

        :param chan: channel e.j. 'SHZ'
        :param julian_date: tuple of ints. year and julian day
        :param date: datetime object. optional
        :param stream: if True, return stream object
        :return: file_path or stream if exist, or None if not exist

        """

        if chan not in self.channels:
            raise ValueError('Channel not loaded')

        if julian_date:
            year = julian_date[0]
            yday = julian_date[1]

        elif date:
            year = date.year
            yday = date.timetuple().tm_yday

        else:
            raise TypeError("A 'date' or 'julian_date' not specified")

        file_name = '.'.join([self.id, f'{chan}.D', str(year), f"{yday:03}"])
        file = os.path.join(self.network.sds_path, str(year), self.network.code, self.code, f'{chan}.D', file_name)
        
        if os.path.isfile(file):
            if stream:
                starttime = readkwargs.get("starttime", None)
                endtime   = readkwargs.get("endtime", None)
                headonly  = readkwargs.get("headonly", False)

                if isinstance(starttime, dt.datetime):
                    starttime = UTCDateTime(starttime)
                
                if isinstance(endtime, dt.datetime):
                    endtime = UTCDateTime(endtime)
                
                return read2(file, starttime=starttime, endtime=endtime, headonly=headonly)
            
            else:
                return file
        else:
            return False


    def __set_dates__(self):
        start = dt.datetime(2969,1,1)
        end   = dt.datetime(1969,1,1)

        sds_path = self.network.sds_path
        net_code = self.network.code
        sta_code = self.code

        for chan in self.channels:
            flist = glob(os.path.join(sds_path, '*', net_code, sta_code, f'{chan}.D', '*'))
            flist = list(filter(lambda x : len(x.split('.')[-1]) == 3, flist))
            datelist = [i[-8:] for i in flist]
            datelist.sort()

            if datelist:
                startdate = dt.datetime.strptime(datelist[0], '%Y.%j')
                sd_st = self.__read_file__(chan, date=startdate, stream=True, headonly=True)
                if sd_st:
                    starttime = sd_st[0].stats.starttime.datetime
                    if starttime < start:
                        start = starttime
                
                enddate = dt.datetime.strptime(datelist[-1], '%Y.%j')
                ed_st = self.__read_file__(chan, date=enddate, stream=True, headonly=True)
                if ed_st:
                    endtime   = ed_st[0].stats.endtime.datetime                    
                    if endtime > end:
                        end = endtime

        self.starttime = start
        self.endtime   = end


    def get_offtimes(self, starttime, endtime, offtime):
        
        if starttime - offtime > self.starttime:
            start = starttime - offtime
        else:
            start = starttime

        if endtime + offtime < self.endtime:
            end = endtime + offtime
        else:
            end = endtime
        
        return (start, end)


    def get_filelist(self, chan, startdate=None, enddate=None):
        """
        Return a list of files for a station's channel
        :param chan: channel e.j. 'SHZ'
        :param julian_date: tuple of ints. year and julian day
        :param date: datetime object. optinal
        :return: file_path or False
        """

        if chan not in self.channels:
            raise TypeError('Channel not loaded')

        if not startdate:
            startdate = self.starttime

        if not enddate:
            enddate = self.endtime

        day_diff = (enddate - startdate).days
        date_list = [startdate + dt.timedelta(days=i) for i in range(day_diff+1)]
        file_list = [self.__read_file__(chan, date=i) for i in date_list]

        return file_list


class NetworkStats(_Stats):
    def __init__(self, code, header):
        super().__init__(header)
        self.code = code

        # buid station list
        station_stats = []
        for _, sta in self.stations.items():
            station_stats.append(StationStats(self, sta))        
        self.stations = station_stats

        self.stations_id = [sta.id for sta in self.stations]
        self.stations_code = list(set([sta.code for sta in self.stations]))


    def __str__(self):
        text =  f"\n   {self.name} network  |  code: {self.code}  |  path: {self.sds_path} ---"
        text += f"\n   Stations  >>  {len(self.stations)} {self.stations_id}"

        return text
    

    def get_latlon(self, **kwargs):
        dout = {}

        hide_sta = kwargs.get("hide_sta", False)
        for sta in self.stations:

            if hide_sta:
                sta_id = sta.location
            else:
                sta_id = '.'.join([sta.code,sta.location])
            
            ans = sta.get_latlon(**kwargs)

            if ans:
                dout[sta_id] = {}

                if kwargs.get("utm", False):
                    x,y = ans[0], ans[1]
                    dout[sta_id]["easting"] = x
                    dout[sta_id]["northing"] = y
                
                else:
                    dout[sta_id]["lat"] = ans[0]
                    dout[sta_id]["lon"] = ans[1]
        
        return dout


class CC8stats(_Stats):
    def __init__(self, header):
        super().__init__(header)
        self.fqidx = [str(fqix) for fqix in range(1, len(self.fq_bands)+1)]


    def check_idx(self, idx):
        if isinstance(idx, (str, int)):
            if isinstance(idx, int):
                idx = str(idx)
            
            if idx in self.fqidx:
                return [idx]
            
            else:
                return []
                
        if isinstance(idx, (list, tuple)):
            return [str(ix) for ix in idx if str(ix) in self.fqidx]


    def __str__(self):
        priorized_keys = ['id', 'locs', 'starttime', 'endtime',\
            'window','overlap','nro_time_bins','last_time_bin',\
            'sample_rate','fq_bands', 'slow_max','slow_int', 'cc_thres', 'slowmap'
            ]
        return self._pretty_str(priorized_keys)
    

    def __write__(self, wdict, fqn, nwin):
        with h5py.File(self.file, "r+") as h5f:
            ltblist = h5f['header'].attrs["last_time_bin"]

            nbin = ltblist[fqn-1]
            
            if self.slowmap:
                h5f[str(fqn)]["slowmap"][nbin+1:nbin+1+nwin,:,:] = wdict["slowmap"]
            
            h5f[str(fqn)]["slowbnd"][nbin+1:nbin+1+nwin,:]   = wdict["slowbnd"]
            h5f[str(fqn)]["bazbnd"][nbin+1:nbin+1+nwin,:]   = wdict["bazbnd"]
                
            for attr in ("slow", "baz", "maac", "rms", "error"):
                h5f[str(fqn)][attr][nbin+1:nbin+1+nwin] = wdict[attr]

            ltblist[fqn-1] = nbin + nwin
            h5f['header'].attrs.modify('last_time_bin', ltblist)
            h5f.flush()


    def get_array(self, return_exclude_locs=True):

        net, sta = self.id.split(".")
        
        try:
            from .network import Array
            arr = Array(net, sta)

            if return_exclude_locs:
                excluded_locs = arr.get_excluded_locs(self.locs)

                return arr, excluded_locs
            
            else:
                return arr
        
        except:
            print(" error reading network file on seisvo.")


    def get_utm(self):
        array, excluded_locs = self.get_array()
        east, north = [], []
        
        for loc, utm in array.utm.items():
            if loc not in excluded_locs:
                east.append(utm["easting"])
                north.append(utm["northing"])
        
        return east, north


class LTEstats(_Stats):
    def __init__(self, header):
        super().__init__(header)

    def get_codelist(self):
        if self.type == "station":
            return tuple(self.channel)

        if self.type == "network":
            return tuple(self.stations)


    def get_kwdict(self):
        kwdict = {"type":self.type, "time_bandwidth":self.time_bandwidth, "pad":self.pad}
        
        if self.type == "station":
            kwdict["polar"] = self.polar
            kwdict["PE_order"] = self.PE_order
            kwdict["PE_tau"] = self.PE_tau
            kwdict["opt_twin"] = self.opt_twin
            kwdict["opt_th"] = self.opt_th
            return kwdict

        if self.type == "network":
            return kwdict


    def __str__(self):
        if self.type == "station":
            priorized_keys = ['type', 'id','channel','starttime',\
            'endtime','window', 'window_olap', 'subwindow','subwindow_olap','sample_rate',\
            'rm_sens','nro_time_bins','last_time_bin','fq_band',\
            'nro_freq_bins','polar']

        else:
            # network not implemented yet!
            priorized_keys = []

        return self._pretty_str(priorized_keys)
    


    def __write__(self, wdict, nwin):
        with h5py.File(self.file, "r+") as h5f:
            nbin = h5f['header'].attrs["last_time_bin"]
            
            if self.type == "station":
                for n, chan in enumerate(self.channel):
                    h5f[chan]["specgram"][nbin+1:nbin+1+nwin,:] = wdict[chan]["specgram"]
                    for attr in ("perm_entr", "dsar", "vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf"):
                        h5f[chan][attr][nbin+1:nbin+1+nwin] = wdict[chan][attr]

                # polar params
                if self.polar:
                    for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh"):
                        h5f["polar"][attr][nbin+1:nbin+1+nwin,:] = wdict["polar"][attr]

                h5f['header'].attrs.modify('last_time_bin', nbin+nwin)
                h5f.flush()

            if self.type == "network":
                # base params
                for sta in self.stations:
                    h5f[sta][nbin+1:nbin+1+nwin,:] = wdict[sta]
                
                h5f["csw"][nbin+1:nbin+1+nwin,:]  = wdict["csw"]
                h5f["vt"][nbin+1:nbin+1+nwin,:,:] = wdict["vt"]

                h5f['header'].attrs.modify('last_time_bin', nbin+nwin)
                h5f.flush()

