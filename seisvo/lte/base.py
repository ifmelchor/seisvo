#!/usr/bin/env python3
# coding=utf-8

import os
import h5py
import datetime as dt
import numpy as np

# from seisvo
from ..signal import get_Peaks, get_Stats, get_PDF, get_freq
from ..stats import LTEstats
from ..plotting import plotPDF, LTESTAplot
from .utils import _LTEProcess
from .peaks import Peaks

STA_SCALAR_PARAMS = ['perm_entr', 'rsam', 'mf', 'hf', 'vlf', 'lf', 'vlar', 'lrar', 'rmar', 'dsar']

STA_VECTOR_PARAMS = ['specgram', 'degree', 'elev', 'rect', 'azimuth', 'phyhh', 'phyvh']


def any_vector(attr_list):
    """
    Check if any attribute of attr_list is a VECTOR_PARAM
    """

    check_list = [attr in STA_VECTOR_PARAMS for attr in attr_list]
    return any(check_list)


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


def _new_LTE(base, lte_file, headers, njobs):
    """
        Create new LTE (hdf5) file.
        base can be station or network object
    """

    if headers["type"] == "station":
        with h5py.File(lte_file, "w-") as f:
            # header dataset
            hdr = f.create_dataset('header',(1,))
            hdr.attrs['id']              = headers['id']
            hdr.attrs['type']            = headers['type']
            hdr.attrs['channel']         = headers['channel']
            hdr.attrs['starttime']       = headers['starttime']
            hdr.attrs['endtime']         = headers['endtime']
            hdr.attrs['window']          = headers['window']
            hdr.attrs['window_olap']     = headers['window_olap']
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
            hdr.attrs['rm_sens']         = headers['rm_sens']
            hdr.attrs['opt_twin']        = headers['opt_twin']
            hdr.attrs['opt_th']          = headers['opt_th']
            
            # print info
            freqbins = headers['nro_freq_bins']
            timebins = headers['nro_time_bins']

            print('')
            print(' LTE file INFO')
            print(' ----------------')
            print(" hdf5_memory info: %s " % lte_file)
            print(' dataset size: ', (timebins, freqbins))
            print('')
            print('    LTE stats   ')
            print(' ----------------')
            
            for info_key in ['id', 'type', 'channel', 'starttime', 'endtime',\
                'window', 'window_olap', 'subwindow', 'subwindow_olap',\
                'sample_rate', 'rm_sens' ,'fq_band' ,'polar']:

                if info_key == "window":
                    print(f'   {info_key} [sec]:  {headers[info_key]}')
                
                elif info_key == "subwindow":
                    print(f'   {info_key} [sec]:  {headers[info_key]}')
                
                else:
                    print(f'   {info_key}:  {headers[info_key]}')
            
            print(' ----------------\n')

            for n, chan in enumerate(headers['channel']):
                chgr = f.create_group(chan)
                chgr.create_dataset('specgram',  (timebins, freqbins), chunks=True, dtype=np.float32)
                chgr.create_dataset('perm_entr', (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('lf',        (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('mf',        (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('hf',        (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('vlf',       (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('rsam',      (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('vlar',      (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('lrar',      (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('rmar',      (timebins,), chunks=True, dtype=np.float32)
                chgr.create_dataset('dsar',      (timebins,), chunks=True, dtype=np.float32)

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

    if headers["type"] == "network":
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
            
            for info_key in ['id', 'type', 'stations', 'starttime', 'endtime',\
                'window', 'window_olap', 'subwindow', 'subwindow_olap',\
                'sample_rate' ,'rm_sens' ,'fq_band']:

                if info_key == "window":
                    print(f' {info_key} [min]:  {headers[info_key]}')
                
                elif info_key == "subwindow":
                    print(f' {info_key} [min]:  {headers[info_key]}')
                
                else:
                    print(f' {info_key}:  {headers[info_key]}')

            for sta in headers['stations']:
                f.create_dataset(sta, (timebins, freqbins), chunks=True, dtype=np.float32)
            
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

        # compute the true time series
        npts  = self.stats.nro_time_bins
        delta = dt.timedelta(seconds=self.stats.window*(1-self.stats.window_olap))
        half_delta = dt.timedelta(seconds=float(self.stats.window/2))
        self._dtime = [self.stats.starttime+half_delta]
        start = self.stats.starttime
        for n in range(self.stats.last_time_bin):
            start += delta
            self._dtime.append(start + half_delta)
        
        self.starttime = self._dtime[0] 
        self.endtime   = self._dtime[-1] 


    def __compute__(self, base, headers, njobs):
        """
        Process LTE file for station and network
        base can be station or network object
        """

        # init process and starttime
        ltep = _LTEProcess(self.stats, headers)
        delta = dt.timedelta(seconds=float(headers["interval"]*60))
        toff  = dt.timedelta(seconds=headers["int_extra_sec"])
        start = self.stats.starttime + toff 
            
        for nint in range(1, headers["nro_intervals"]+1):
            start -= toff
            end    = start + delta + toff
            
            if nint == headers["nro_intervals"]:
                last = True
                end = self.stats.endtime
            else:
                last = False

            data = None
            sort = None
            if self.stats.type == "station":
                stream = base.get_stream(start, end, channel=self.stats.channel,\
                    rm_sens=self.stats.rm_sens, sample_rate=self.stats.sample_rate)
                if self.stats.polar and len(stream) == 3:
                    sort = "ZNE"

            if self.stats.type == "network":
                stream = base.get_stream(start, end, sta_code=self.stats.stations,\
                    avoid_exception=True, rm_sens=self.stats.rm_sens,\
                    sample_rate=self.stats.sample_rate)

            if stream and stream.get_bounds() == (start,end):
                data = stream.to_array(sort=sort)

            ltep.run(data, start, end, last)

            # stack jobs until njobs
            if len(ltep.processes) == njobs:
                ltep.wait(nint)
            
            # advance 
            start = end
        
        # check if there are any process running
        if len(ltep.processes) > 0:
            ltep.wait(nint)


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

                if attr in ('perm_entr', 'rsam', 'mf', 'hf', 'vlf', 'lf', 'vlar', 'lrar', 'rmar', 'dsar'):
                    ts = f.get(chan)[attr][n0:nf]

                if attr in ('degree', 'elev', 'rect', 'azimuth', 'phyhh', 'phyvh'):
                    ts = f.get("polar")[attr][n0:nf,:]

                    if attr == "azimuth" and azimuth_ambiguity:
                        ts = np.where(ts>180, ts-180, ts)
                
        if self.stats.type == "network":
            with h5py.File(self.stats.file, "r") as f:
                if attr in self.stats.stations or attr == "csw":
                    ts = f[attr][n0:nf,:]
                else:
                    ts = f[attr][n0:nf,:,:]

        return ts


    def __str__(self):
        return self.stats.__str__()


    def get_time(self, starttime, endtime):

        if not starttime:
            starttime = self.starttime
        
        if not endtime:
            endtime = self.endtime
        
        n0 = np.argmin(np.abs(np.array(self._dtime) - starttime))
        nf = np.argmin(np.abs(np.array(self._dtime) - endtime)) + 1
        
        return self._dtime[n0:nf], (n0, nf) 
        

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
                starttime       = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S.%f'),
                endtime         = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S.%f'),
                sample_rate     = int(hdr.attrs['sample_rate']),
                nro_time_bins   = int(hdr.attrs['nro_time_bins']),
                last_time_bin   = int(hdr.attrs['last_time_bin']),
                fq_band         = [float(fq) for fq in hdr.attrs['fq_band']],
                nro_freq_bins   = int(hdr.attrs['nro_freq_bins']),
                window          = float(hdr.attrs['window']),
                window_olap     = float(hdr.attrs['window_olap']),
                subwindow       = float(hdr.attrs['subwindow']),
                subwindow_olap  = float(hdr.attrs['subwindow_olap']),
                rm_sens         = bool(hdr.attrs['rm_sens']),
                polar           = bool(hdr.attrs['polar']),
                opt_th          = float(hdr.attrs['opt_th']),
                opt_twin        = float(hdr.attrs['opt_twin']),
                PE_tau          = int(hdr.attrs["PE_tau"]),
                PE_order        = int(hdr.attrs["PE_order"]),
                time_bandwidth  = float(hdr.attrs['time_bandwidth']),
                pad             = float(hdr.attrs['pad']),
            )

        # set attr
        attrs = STA_SCALAR_PARAMS + [STA_VECTOR_PARAMS[0]]
        if dstats["polar"]:
            attrs += STA_VECTOR_PARAMS[1:]
        
        self.__set_stats__(dstats, attrs)


    def check_attr(self, attr, chan=None, only_scalars=False, only_vectors=False):
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


    def get(self, starttime=None, endtime=None, attr=None, chan=None, db_scale=True, azimuth_ambiguity=True):
        
        if not starttime:
            starttime = self.stats.starttime
        else:
            assert starttime >= self.stats.starttime
        
        if not endtime:
            endtime = self.stats.endtime
        else:
            assert endtime <= self.stats.endtime

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
        dout["dtime"], (n0,nf) = self.get_time(starttime, endtime)

        if any([attr in STA_VECTOR_PARAMS for attr in attr_list]):
            if self.stats.subwindow > 0:
                lwin = self.stats.subwindow*self.stats.sample_rate
            else:
                lwin = self.stats.window*self.stats.sample_rate

            dout["freq"] = get_freq(lwin, self.stats.sample_rate, fq_band=self.stats.fq_band, pad=self.stats.pad)[0]

        for attr in attr_list:
            if attr in STA_SCALAR_PARAMS + [STA_VECTOR_PARAMS[0]]:
                for n, chan in enumerate(chan_list):
                    stakwargs["chan"] = chan
                    dout[chan][attr] = self.__read__(attr, n0, nf, stakwargs=stakwargs)
            else:
                stakwargs["chan"] = None
                dout[attr] = self.__read__(attr, n0, nf, stakwargs=stakwargs)

        return LTEout(self.stats, dout)


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
        

    def plot(self, interval, attr, chan, starttime=None, olap=0.1, db=None, **fig_kwargs):
        """
        Inizialice the LTE GUI. Be default starttime is lte.starttime.
        the fig_kwars controls the cmap of vector parameters and the limits of the attributes
        """
        
        from ..gui import load_ltewidget

        # check attr and channels
        attr_list = self.check_attr(attr)
        chan_list = self.check_chan(chan)

        # init widget
        widget = load_ltewidget(self, db, starttime, interval, attr_list,\
            chan_list, olap=olap, **fig_kwargs)

        return None


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
    

class LTEout(object):
    """
    this object allow interact with station LTEout
    """

    def __init__(self, ltestats, dout):
        self.ltestats = ltestats
        self.starttime = dout["dtime"][0]
        self.endtime = dout["dtime"][-1]
        self._dout = dout

        # define some stats
        self.npts = len(dout["dtime"])
        self.npfs = len(dout["freq"])
        self.chan_list = []
        self.chan_attr_list = []
        self.polar_attr_list = []

        for key, item in dout.items():
            if key in ltestats.channel:
                self.chan_list.append(key)
                for chankey in list(item.keys()):
                    if chankey not in self.chan_attr_list:
                        self.chan_attr_list.append(chankey)
            else:
                if key not in ("freq", "dtime"):
                    self.polar_attr_list.append(key)
                    

    def __str__(self):
        txt_to_return =  f'\n   LTE file  : {self.ltestats.file}'
        txt_to_return += f'\n   starttime :  {self.starttime.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   endtime   :  {self.endtime.strftime("%d %B %Y %H:%M")}'
        txt_to_return += f'\n   channels  :  {self.chan_list}'
        txt_to_return += f'\n   attribute :  {self.chan_attr_list+self.polar_attr_list}'
        txt_to_return +=  f'\n'
        return txt_to_return
    
    
    def check_chan(self, chan):

        if isinstance(chan, str):
            if chan in self.chan_list:
                chan_list = [chan]
            else:
                chan_list = []
        else:
            chan_list = [ch for ch in chan if ch in self.chan_list]

        if not chan_list:
            print(" warn :: no channel found")

        return chan_list
    

    def check_attr(self, attr, which=None):
        # check attr list
        all_attr = self.polar_attr_list + self.chan_attr_list
        
        if isinstance(attr, str):
            if attr in all_attr:
                attr_list = [attr]
            else:
                attr_list = []
        else:
            attr_list = [at for at in attr if at in all_attr]
                    
        if which:
            attr_list = attr_filt(attr_list, which)
        
        if not attr_list:
            print(" warn :: no attributes found")
            
        return attr_list
    

    def any_vector(self, attr_list=None):

        if not attr_list:
            attr_list = self.polar_attr_list + self.chan_attr_list

        return any_vector(attr_list)


    def get_stats(self, attr, chan, bw_method=None):
        """
        Return (min, max, mean, mode) of an specific scalar-only attribute
        """

        # check attr and channel
        attr = self.check_attr(attr, "scalar")[0]

        if chan:
            chan = self.check_chan(chan)[0]
            data = self._dout[chan][attr]
        else:
            data = self._dout[attr]

        stats = get_Stats(data[np.isfinite(data)], bw_method=bw_method)

        return stats
    

    def get_pdf(self, vector_attr, chan=None, db_scale=True, ymin=None, ymax=None, plot=True, show=True, title=None, **kde_kwargs):
        """
        Return the PDF of a vector attribute,
        if "specgram", channel must be specified
        """

        attr = self.check_attr(vector_attr, "vector")[0]
        
        if chan:
            chan_list = self.check_chan(chan)[0]
            data = self._dout[chan][attr]
            if db_scale:
                data = 10*np.log(data)
        else:
            data = self._dout[attr]
        
        if not ymin:
            ymin = np.floor(np.nanmin(data))
        
        if not ymax:
            ymax  = np.ceil(np.nanmax(data))
        
        space = np.linspace(ymin, ymax, num=1000).reshape(-1,1)

        pdf = get_PDF(data, space, **kde_kwargs)

        if plot:
            return plotPDF(pdf, space.reshape(-1,), self._dout["freq"], show=show, title=title)

        return pdf, space.reshape(-1,), self._dout["freq"]


    def plot(self, chan=None, attr=["specgram", "rsam", "lf", "hf"], plot=True, **plot_kw):
        
        if not chan:
            chan_list = [self.chan_list[0]]
        else:
            chan_list = self.check_chan(chan)
        
        attr_list = self.check_attr(attr)
        
        fig, _ = LTESTAplot(self, chan_list, attr_list, plot=plot, **plot_kw)

        return fig

