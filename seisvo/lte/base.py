#!/usr/bin/env python3
# coding=utf-8

import os
import h5py
import datetime as dt
import numpy as np

# from ..station import stalist_data
from ..signal.utils import get_freq
from ..signal.ltejl import get_Peaks
from ..stats import LTEstats
from .utils import _LTEProcess, LTEout
from .peaks import Peaks

STA_SCALAR_PARAMS = ['fq_dominant', 'fq_centroid', 'energy', 'perm_entr', 'rsam', 'mf', 'hf', 'vlf', 'lf', 'vlar', 'lrar', 'rmar', 'dsar']

STA_VECTOR_PARAMS = ['specgram', 'degree', 'elev', 'rect', 'azimuth', 'phyhh', 'phyvh']


def _new_LTE(base, lte_file, headers, njobs, ltetype):
    """
        Create new LTE (hdf5) file.
        base can be station or network object
    """

    if ltetype == "station":
        with h5py.File(lte_file, "w-") as f:
            # header dataset
            hdr = f.create_dataset('header',(1,))
            hdr.attrs['id']              = headers['id']
            hdr.attrs['type']            = headers['type']
            hdr.attrs['channel']         = headers['channel']
            hdr.attrs['starttime']       = headers['starttime']
            hdr.attrs['endtime']         = headers['endtime']
            hdr.attrs['window']          = headers['window']
            hdr.attrs['subwindow']       = headers['subwindow']
            hdr.attrs['subwindow_olap']  = headers['subwindow_olap']
            hdr.attrs['sample_rate']     = headers['sample_rate']
            hdr.attrs['pad']             = headers['pad']
            hdr.attrs['fq_band']         = headers['fq_band']
            hdr.attrs['nro_time_bins']   = headers['nro_time_bins']
            hdr.attrs['nro_freq_bins']   = headers['nro_freq_bins']
            hdr.attrs['last_time_bin']   = -1
            hdr.attrs['PE_tau']          = headers['pe_tau']
            hdr.attrs['PE_order']        = headers['pe_order']
            hdr.attrs['time_bandwidth']  = headers['time_bandwidth']
            # optional param
            hdr.attrs['polar']           = headers['polar']
            hdr.attrs['opt_params']      = headers['opt_params']
            hdr.attrs['rm_sens']         = headers['rm_sens']
            hdr.attrs['opt_twin']        = headers['opt_twin']
            hdr.attrs['opt_th']          = headers['opt_th']
            
            # print info
            freqbins = headers['nro_freq_bins']
            timebins = headers['nro_time_bins']

            print('')
            print(' LTE file INFO')
            print(' -------------')
            print(" hdf5_memory info: %s " % lte_file)
            print(' --- dataset size: ', (timebins, freqbins))
            print('')
            print(' LTE stats:')
            
            for info_key in ['id', 'channel', 'starttime', 'endtime', 'window', 'subwindow' ,'subwindow_olap' ,'sample_rate' ,'rm_sens' ,'fq_band' ,'polar' ,'opt_params']:

                if info_key == "window":
                    print(f' {info_key} [min]:  {headers[info_key]}')
                
                elif info_key == "subwindow":
                    print(f' {info_key} [sec]:  {headers[info_key]}')
                
                else:
                    print(f' {info_key}:  {headers[info_key]}')

            for chan in headers['channel']:
                chgr = f.create_group(chan)
                chgr.create_dataset('specgram',    (timebins, freqbins), chunks=True, dtype=np.float32)
                chgr.create_dataset('perm_entr',   (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('fq_dominant', (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('fq_centroid', (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('energy',      (timebins,), chunks=True, dtype=np.float32)
            
            if headers['opt_params']:
                opt  = f.create_group("opt")
                opt.create_dataset('lf',   (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('mf',   (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('hf',   (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('vlf',  (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('rsam', (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('vlar', (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('lrar', (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('rmar', (timebins,), chunks=True, dtype=np.float32)
                opt.create_dataset('dsar', (timebins,), chunks=True, dtype=np.float32)

            if headers['polar']:
                polar = f.create_group("polar")
                polar.create_dataset('degree',  (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('rect',    (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('azimuth', (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('elev',    (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('phyhh',   (timebins, freqbins), chunks=True, dtype=np.float32)
                polar.create_dataset('phyvh',   (timebins, freqbins), chunks=True, dtype=np.float32)

            f.flush()

        lte = StationLTE(lte_file)

    if ltetype == "network":
        with h5py.File(lte_file, "w-") as f:
            # header dataset
            hdr = f.create_dataset('header',(1,))
            hdr.attrs['id']              = headers['id']
            hdr.attrs['type']            = headers['type']
            hdr.attrs['stations']        = headers['stations']
            hdr.attrs['starttime']       = headers['starttime']
            hdr.attrs['endtime']         = headers['endtime']
            hdr.attrs['window']          = headers['window']
            hdr.attrs['window_olap']     = headers['window_olap']
            hdr.attrs['subwindow']       = headers['subwindow']
            hdr.attrs['subwindow_olap']  = headers['subwindow_olap']
            hdr.attrs['nro_time_bins']   = headers['nro_time_bins']
            hdr.attrs['last_time_bin']   = -1
            hdr.attrs['fq_band']         = headers['fq_band']
            hdr.attrs['nro_freq_bins']   = headers['nro_freq_bins']
            hdr.attrs['pad']             = headers['pad']
            hdr.attrs['sample_rate']     = headers['sample_rate']
            hdr.attrs['time_bandwidth']  = headers['time_bandwidth']
            # optional param
            hdr.attrs['rm_sens']         = headers['rm_sens']
 
            # print info
            freqbins = headers['nro_freq_bins']
            timebins = headers['nro_time_bins']

            print('')
            print(' LTE file INFO')
            print(' -------------')
            print(" hdf5_memory info: %s " % lte_file)
            print(' --- dataset size: ', (timebins, freqbins))
            print('')
            print(' LTE stats:')
            
            for info_key in ['id', 'stations', 'starttime', 'endtime', 'window', 'window_olap', 'subwindow' ,'subwindow_olap' ,'sample_rate' ,'rm_sens' ,'fq_band']:

                if info_key == "window":
                    print(f' {info_key} [min]:  {headers[info_key]}')
                
                elif info_key == "subwindow":
                    print(f' {info_key} [min]:  {headers[info_key]}')
                
                else:
                    print(f' {info_key}:  {headers[info_key]}')

            for sta in headers['stations']:
                stgr = f.create_group(sta)
                stgr.create_dataset('specgram', (timebins, freqbins), chunks=True, dtype=np.float32)
            
            f.create_dataset('csw', (timebins, freqbins), chunks=True, dtype=np.float32)
            f.create_dataset('vt', (timebins, freqbins, len(headers['stations'])), chunks=True, dtype=np.complex64)

            f.flush()

        lte = NetworkLTE(lte_file)
    
    lte.__compute__(base, headers, njobs)
    

class _LTE(object):
    def __init__(self, lte_file):
        assert os.path.isfile(lte_file)
    

    def __set_stats__(self, dstats, lattrs):
        self.stats = LTEstats(dstats)
        self.stats.__add_attr__(lattrs)


    def __compute__(self, base, headers, njobs):
        """
        Process LTE file for station and network
        base can be station or network object
        """

        # init process and starttime
        ltep = _LTEProcess(self.stats, headers["nwin"], headers["lwin"], headers["nswin"], headers["lswin"], headers["nadv"])

        start = self.stats.starttime
        interval = dt.timedelta(hours=headers["interval"])
        nro_ints = headers['nro_intervals']
        nint = 1
        while start + interval <= self.stats.endtime:
            end = start + interval

            data = None

            if self.stats.type == "station":
                stream = base.get_stream(start, end, channel=ltestats.channel,\
                    remove_sensitivity=self.stats.rm_sens, sample_rate=self.stats.sample_rate)

                if stream and stream.get_bounds() == (start,end):
                    data = stream.to_array(sort="ZNE")

            if self.stats.type == "network":
                stream = base.get_stream(start, end, sta_code=self.stats.stations, avoid_exception=True,
                    remove_sensitivity=self.stats.rm_sens, sample_rate=self.stats.sample_rate)

                if stream and stream.get_bounds() == (start,end):
                    data = stream.to_array()

            # send to a cpu
            ltep.run(data, start, end)

            # stack jobs until njobs
            if len(ltep.processes) == njobs:
                ltep.wait(f"{nint}/{nro_ints}")
            
            # advance 
            start += interval
            nint += 1
        
        # check if there are any process running
        if len(ltep.processes) > 0:
            ltep.wait(f"{nint}/{nro_ints}")


    def __read__(self, attr, n0, nf, stakwargs={}, netkwargs={}):
        """
        function for read dataseries (ts) in lte files
        """
        ts = None
        if self.stats.type == "station":
            chan = stakwargs["chan"]
            db_scale = stakwargs["db_scale"] 
            azimuth_ambiguity = stakwargs["azimuth_ambiguity"]

            with h5py.File(self.stats.file, "r") as f:
                if attr == 'specgram':
                    ts = f.get(chan)[attr][n0:nf,:]

                    if db_scale:
                        ts = 10*np.log10(ts)
                
                if attr in ('fq_dominant', 'fq_centroid', 'energy', 'perm_entr'):
                    ts = f.get(chan)[attr][n0:nf]

                    if attr == "energy" and db_scale:
                        ts = 10*np.log10(ts)
            
                if attr in ('degree', 'elev', 'rect', 'azimuth', 'phyhh', 'phyvh'):
                    ts = f.get("polar")[attr][n0:nf,:]

                    if attr == "azimuth" and azimuth_ambiguity:
                        ts = np.where(ts>180, ts-180, ts)
                
                if attr in ('rsam', 'mf', 'hf', 'vlf', 'lf', 'vlar', 'lrar', 'rmar', 'dsar'):
                    ts = f.get("opt")[attr][n0:nf]

        # if self.stats.type == "network":

        return ts


    def __str__(self):
        return self.stats.__str__()


    def get_time(self, starttime=None, endtime=None, olap=0):

        if not starttime:
            starttime = self.stats.starttime
        
        if not endtime:
            endtime = self.stats.endtime
        
        assert endtime > starttime
        assert starttime >= self.stats.starttime
        assert endtime <= self.stats.endtime

        start_diff = (starttime - self.stats.starttime).total_seconds()*self.stats.sample_rate # sp
        n0 =  int(np.floor(start_diff/(self.stats.window*60*self.stats.sample_rate)))
    
        end_diff = (endtime - self.stats.starttime).total_seconds()*self.stats.sample_rate # sp
        nf =  int(np.ceil(end_diff/(self.stats.window*60*self.stats.sample_rate)))

        datetime_list = [starttime + dt.timedelta(minutes=int(k*self.stats.window)) for k in range(n0,nf)]
        duration  = (endtime-starttime).total_seconds()/60
        time_list = np.linspace(0, duration, nf-n0)
        
        return time_list, datetime_list, (n0,nf)


class StationLTE(_LTE):
    def __init__(self, lte_file):
        super().__init__(lte_file)

        # set stats
        with h5py.File(lte_file, "r") as f:
            hdr = f['header']
            dstats = dict(
                file            = lte_file,
                type            = hdr.attrs["type"],
                id              = hdr.attrs['id'],
                channel         = list(hdr.attrs['channel']),
                starttime       = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                endtime         = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                sample_rate     = int(hdr.attrs['sample_rate']),
                nro_time_bins   = int(hdr.attrs['nro_time_bins']),
                last_time_bin   = int(hdr.attrs['last_time_bin']),
                fq_band         = [float(fq) for fq in hdr.attrs['fq_band']],
                nro_freq_bins   = int(hdr.attrs['nro_freq_bins']),
                window          = hdr.attrs['window'],
                subwindow       = hdr.attrs['subwindow'],
                subwindow_olap  = hdr.attrs['subwindow_olap'],
                rm_sens         = bool(hdr.attrs['rm_sens']),
                polar           = bool(hdr.attrs['polar']),
                opt_params      = bool(hdr.attrs['opt_params']),
                opt_th          = float(hdr.attrs['opt_th']),
                opt_twin        = float(hdr.attrs['opt_twin']),
                PE_tau          = int(hdr.attrs["PE_tau"]),
                PE_order        = int(hdr.attrs["PE_order"]),
                time_bandwidth  = float(hdr.attrs['time_bandwidth']),
                pad             = float(hdr.attrs['pad']),
            )

        # set attr
        attrs = [STA_VECTOR_PARAMS[0]] + STA_SCALAR_PARAMS[:4]

        if dstats["opt_params"]:
            attrs += STA_SCALAR_PARAMS[4:]
        
        if dstats["polar"]:
            attrs += STA_VECTOR_PARAMS[1:]
        
        self.__set_stats__(dstats, attrs)


    def check_attr(self, attr, only_scalars=False, only_vectors=False):
        """
        Check if attribute list or string is available and return available attributes
        """

        if not isinstance(attr, list):
            list_attr = [attr]
        else:
            list_attr = attr

        if not all([attr in self.stats.attributes for attr in list_attr]):
            attrs = [at not in self.stats.attributes for at in list_attr]
            pos = list(filter(lambda x: attrs[x], range(len(list_attr))))
            not_availabel_attr = np.array(list_attr)[pos]
            print('warn: attributes %s not available' % not_availabel_attr)

            return_list = [attr for attr in list_attr if attr in self.stats.attributes]
            
            if not return_list:
                print('available attr: %s' % self.stats.attributes)
                return

        else:
            return_list = list_attr
        
        if only_scalars:
            return_list = [attr for attr in return_list if attr in STA_SCALAR_PARAMS]
        
        if only_vectors:
            return_list = [attr for attr in return_list if attr in STA_VECTOR_PARAMS]
            
        return return_list
    

    def check_chan(self, chan):
        """
        Check if chan list or string is available and return available channels (only for station LTE mode)
        """

        if isinstance(chan, list):
            chan_list = [ch for ch in chan if ch in self.stats.channel]
        
        else:
            if isinstance(chan, type(None)):
                chan_list = list(self.stats.channel)
            else:
                if chan in self.stats.channel:
                    chan_list = [chan]
                else:
                    chan_list = []
        
        return chan_list


    def get(self, attr=None, chan=None, starttime=None, endtime=None, db_scale=True, azimuth_ambiguity=True):
        
        stakwargs = {"chan":None, "db_scale":db_scale, "azimuth_ambiguity":azimuth_ambiguity}
        
        # if attr is None, return all attributes
        if not attr:
            attr_list = self.stats.attributes
        else:
            attr_list = self.check_attr(attr)
        
        # if chan is None, return all channels
        if not chan:
            chan_list = self.stats.channel
        else:
            chan_list = self.check_chan(chan)
        
        # create dictionary
        dout = {}
        for chan in chan_list:
            dout[chan] = {}
        
        # define time series
        dout["time"], dout["dtime"], (n0,nf) = self.get_time(starttime=starttime, endtime=endtime)

        if any([attr in STA_VECTOR_PARAMS for attr in attr_list]):
            if self.stats.subwindow:
                lwin = self.stats.subwindow*self.stats.sample_rate
            else:
                lwin = self.stats.window*60*self.stats.sample_rate

            dout["freq"] = get_freq(lwin, self.stats.sample_rate, fq_band=self.stats.fq_band, pad=self.stats.pad)[0]

        for attr in attr_list:
            if attr in ('fq_dominant', 'fq_centroid', 'energy', 'perm_entr', 'specgram'):
                for chan in chan_list:
                    stakwargs["chan"] = chan
                    dout[chan][attr] = self.__read__(attr, n0, nf, stakwargs=stakwargs)
            else:
                stakwargs["chan"] = None
                dout[attr] = self.__read__(attr, n0, nf, stakwargs=stakwargs)

        return LTEoutSTA(self, dout)


    def get_peaks(self, chan=None, starttime=None, endtime=None, fq_range=(), peak_thresholds={}):
        """

        Extract dominant peaks as described in Melchor et al. 2022 (https://doi.org/10.1016/j.jsames.2022.103961)

        Parameters
        ----------
        fq_range : tuple, optional
            by default stats.fq_band
        peak_thresholds : dict, optional
            by default {'fq_dt': 0.05, 'sxx_th': 0.95, 'pdg_th': 0.8, 'rect_th': 0.7, 'pdg_std':0.1, 'rect_std':0.1, 'azim_std':10, 'elev_std':10}

        Returns
        -------
        [Peak object]
        """

        if not self.type == "station":
            return None 
        
        if not self.stats.polar:
            return None

        from juliacall import Main as jl

        default_peak_thresholds = {
                'fq_dt': 0.05,
                'sxx_th': 0.7,
                'pdg_th': 0.7,
                'pdg_std':0.1,
                'rect_th': 0.7,
                'rect_std':0.1,
                'azim_std':10,
                'elev_std':10,
                'npts_min':4
        }

        # build the peak_thresholds dict
        for key in default_peak_thresholds.keys():
            if not peak_thresholds.get(key, None):
                peak_thresholds[key] = default_peak_thresholds.get(key)

        # plot model parameters
        print(f'\n  File: {os.path.basename(self.file_)}')
        print('  ----- Model param ------')
        for key, item in peak_thresholds.items():
            print(f'  {key:>10}     ', item)
        print('  ------------------------\n')

        attr = ['specgram', 'degree', 'elev', 'rect', 'azimuth']
        gout = self.get(attr=attr, chan=chan, starttime=starttime, endtime=endtime)

        freq = gout._dout["freq"]
        if fq_range:
            assert fq_range[0] >= freq[0]
            assert fq_range[1] <= freq[-1]
            fleft  = np.argmin(np.abs(freq - fq_range[0]))
            fright = np.argmin(np.abs(freq - fq_range[1])) + 1
            freq = freq[fleft:fright]
        else:
            fq_range = self.fq_band
            fleft = 0 
            fright = len(freq)-1
        
        nchan = len(gout.chan_list)
        sxx = np.empty((nchan, gout.npts_, len(freq)), dtype=np.float32)
        for c, chan in enumerate(gout.chan_list):
            sxx[c,:,:] = gout._dout[chan]["specgram"][:,fleft:fright]

        freq = np.array(freq, dtype=np.float32)
        degr = gout._dout["degree"][:,fleft:fright]
        rect = gout._dout["rect"][:,fleft:fright]
        azim = gout._dout["azimuth"][:,fleft:fright]
        elev = gout._dout["elev"][:,fleft:fright]

        dpeaks = get_Peaks(chan_list, sxx, degree, rect, azimuth, elev, freq, peak_thresholds)

        return Peaks(dpeaks, self.file_, gout.starttime_, gout.endtime_, self.stats.window, fq_range, peak_thresholds)
        

    def plot(self, attr, chan, day_interval, starttime=None, lde=None):
        attr_list = self.check_attr(attr)
        chan_list = self.check_chan(chan)

        print("warn :: glte.py not upgraded yet")

        return None


    @staticmethod
    def any_vector(attr_list):
        """
        Check if any attribute of attr_list is a VECTOR_PARAM
        """

        check_list = [attr in STA_VECTOR_PARAMS for attr in attr_list]
        return any(check_list)
    

    @staticmethod
    def attr_filt(attr_list, which):
        """
        Filter a list of attributes to only vector or scalar.
        which must be a string of "scalar" or "vector"
        """

        assert which in ("scalar", "vector")
        assert isinstance(attr_list, list)

        filter_list = []

        for attr in attr_list:
            if which == "scalar" and attr in STA_SCALAR_PARAMS:
                filter_list.append(attr)
            
            if which == "vector" and attr in STA_VECTOR_PARAMS:
                filter_list.append(attr)
                
        return filter_list


class NetworkLTE(_LTE):
    def __init__(self, lte_file):
        super().__init__(lte_file)
        
        # set stats
        with h5py.File(lte_file, "r") as f:
            hdr = f['header']
            dstats = dict(
                file            = lte_file,
                type            = hdr.attrs["type"],
                id              = hdr.attrs['id'],
                stations        = list(hdr.attrs['stations']),
                starttime       = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                endtime         = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                sample_rate     = int(hdr.attrs['sample_rate']),
                nro_time_bins   = int(hdr.attrs['nro_time_bins']),
                last_time_bin   = int(hdr.attrs['last_time_bin']),
                fq_band         = [float(fq) for fq in hdr.attrs['fq_band']],
                nro_freq_bins   = int(hdr.attrs['nro_freq_bins']),
                window          = hdr.attrs['window'],
                window_olap     = hdr.attrs['window_olap'],
                subwindow       = hdr.attrs['subwindow'],
                subwindow_olap  = hdr.attrs['subwindow_olap'],
                rm_sens         = bool(hdr.attrs['rm_sens']),
                time_bandwidth  = float(hdr.attrs['time_bandwidth']),
                pad             = float(hdr.attrs['pad']),
            )

        self.__set_stats__(dstats, ["specgram", "csw", "vt"])
    

