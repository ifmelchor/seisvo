#!/usr/bin/python3
# coding=utf-8

import os
import pickle
import datetime as dt
from obspy.core.util.attribdict import AttribDict
from obspy.clients.nrl import NRL
from seisvo import __seisvo__
from shutil import rmtree
from importlib.machinery import SourceFileLoader

NET_PATH = '%s/networks' % __seisvo__
RESP_PATH = '%s/respfiles' % __seisvo__

class StaFile(AttribDict):
    def __init__(self, header):
        super(StaFile, self).__init__(header)
        self.id = None

    def __str__(self):
        priorized_keys = [
        'id',
        'net',
        'code',
        'loc',
        'chan',
        'sampling_rate',
        'central',
        'lat',
        'lon',
        'type',
        'elev'
        ]
        return self._pretty_str(priorized_keys)

    def set_id(self):
        self.id = '%s.%s.%s' % (self.net, self.code, self.loc)


class NetFile(AttribDict):
    def __init__(self, header):
        super(NetFile, self).__init__(header)

        for s, sta in enumerate(self.stations):
            self.stations[s] = StaFile(sta)
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


class RespFile(AttribDict):
    def __init__(self, header):
        super(RespFile, self).__init__(header)
        self.__download__()

    def __str__(self):
        priorized_keys = [
        'code',
        'loc',
        'sensor_keys',
        'datalogger_keys'
        ]
        return self._pretty_str(priorized_keys)
    
    def __download__(self):

        if self.sensor_keys and self.datalogger_keys and self.file:
            pass
        
        else:
            print(' STA:%s.%s :: no information in .resp file for download respfile of this station' % (self.code, self.loc))
            return

        path_out = RESP_PATH
        if not os.path.isdir(path_out):
            os.makedirs(path_out)
        
        file_out_name = os.path.join(path_out, self.file)

        if not os.path.isfile(file_out_name):
            nrl = NRL()
            resp = nrl.get_response(
                sensor_keys=self.sensor_keys,
                datalogger_keys=self.datalogger_keys
            )

            file_out = open(file_out_name, 'wb')
            pickle.dump(resp, file_out, protocol=pickle.HIGHEST_PROTOCOL)
            file_out.close()


    def load(self):
        file_out_name = os.path.join(RESP_PATH, self.file)
        if os.path.isfile(file_out_name):
            with open(file_out_name, 'rb') as f:
                resp = pickle.load(f)
        else:
            resp = None

        return resp


def get_network(net_code):
    net_dir = os.path.join(NET_PATH, net_code)
    net_file = '%s.net' % net_code
    net_conf = os.path.join(net_dir, net_file)

    if os.path.isdir(net_dir):
        if os.path.isfile(net_conf):
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

    return NetFile(net.network)


def get_respfile(net_code, sta, loc, starttime=None, endtime=None):
    net_dir = os.path.join(NET_PATH, net_code)
    resp_file = '%s.resp' % net_code
    resp_path = os.path.join(net_dir, resp_file)

    if os.path.isdir(net_dir):
        if os.path.isfile(resp_path):
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
                                return RespFile(dic)
                        else:
                            return RespFile(dic)
                return None
            # except:
                # raise ValueError('importing network respfile')
        else:
            raise ValueError('reading network directory')    
    else:
        raise ValueError('finding network directory')


def remove_cache(pathdir):
    cache = '%s/__pycache__' % pathdir
    
    if os.path.isdir(cache):
        rmtree(cache)