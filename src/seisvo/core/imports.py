#!/usr/bin/env python
# coding=utf-8

from obspy.core.util.attribdict import AttribDict
from importlib.machinery import SourceFileLoader
from os import path
from shutil import rmtree
from seisvo import __seisvo__
import datetime as dt

NET_PATH = '%s/networks' % __seisvo__


class NetInfo(AttribDict):
    def __init__(self, header):
        super(NetInfo, self).__init__(header)

        for s, sta in enumerate(self.stations):
            self.stations[s] = StaInfo(sta)
            self.stations[s].sdsdir = self.sdsdir
            self.stations[s].net = self.code
            self.stations[s].set_id()

        self.stations_info = [sta.id for sta in self.stations]

    def __str__(self):
        priorized_keys = [
        'name',
        'code',
        'latOrig',
        'lonOrig',
        'sdsdir',
        'stations_info'
        ]
        return self._pretty_str(priorized_keys)


class StaInfo(AttribDict):

    def __init__(self, header):
        super(StaInfo, self).__init__(header)
        self.id = None
        self.starttime = None
        self.endtime = None

    def __str__(self):
        priorized_keys = [
        'id',
        'net',
        'code',
        'loc',
        'chan',
        'starttime',
        'endtime',
        'lat',
        'lon',
        'elev'
        ]
        return self._pretty_str(priorized_keys)

    def set_id(self):
        self.id = '%s.%s.%s' % (self.net, self.code, self.loc)


class RespInfo(AttribDict):
    def __init__(self, header):
        super(RespInfo, self).__init__(header)    

    def __str__(self):
        priorized_keys = [
        'code',
        'loc',
        'sensor_keys',
        'datalogger_keys'
        ]
        return self._pretty_str(priorized_keys)


def get_network(net_code):
    
    net_dir = path.join(NET_PATH, net_code)
    net_file = '%s.net' % net_code
    net_conf = path.join(net_dir, net_file)

    if path.isdir(net_dir):
        if path.isfile(net_conf):
            try:
                loader = SourceFileLoader(net_code, net_conf)
                net = loader.load_module()
                remove_cache(net_dir)
            except:
                raise ValueError('importing network file')
        else:
            raise ValueError('reading network directory')    
    else:
        raise ValueError('finding network directory')

    return NetInfo(net.network)

def get_respfile(net_code, sta, loc, starttime=None, endtime=None):
    net_dir = path.join(NET_PATH, net_code)
    resp_file = '%s.resp' % net_code
    resp_path = path.join(net_dir, resp_file)

    if path.isdir(net_dir):
        if path.isfile(resp_path):
            # try:
                loader = SourceFileLoader(net_code, resp_path)
                resp = loader.load_module()
                remove_cache(net_dir)
                for dic in resp.response:
                    if dic['code'] == sta and dic['loc'] == loc:
                        if starttime and endtime:
                            st_resp = dt.datetime(*dic['starttime'])
                            et_resp = dt.datetime(*dic['endtime'])
                            cond1 = st_resp <= starttime <= et_resp
                            cond2 = st_resp <= endtime <= et_resp
                            if cond1 and cond2:
                                return RespInfo(dic)
                        else:
                            return RespInfo(dic)
                return None
            # except:
                # raise ValueError('importing network respfile')
        else:
            raise ValueError('reading network directory')    
    else:
        raise ValueError('finding network directory')


def remove_cache(pathdir):
    cache = '%s/__pycache__' % pathdir
    
    if path.isdir(cache):
            rmtree(cache)
