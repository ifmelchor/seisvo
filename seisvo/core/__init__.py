#!/usr/bin/python3
# coding=utf-8

import os
import pickle
import datetime as dt
from obspy.core.util.attribdict import AttribDict
from obspy.clients.nrl import NRL

from seisvo import NET_PATH, RESP_PATH
from shutil import rmtree
from importlib.machinery import SourceFileLoader


class StaFile(AttribDict):
    def __init__(self, header):
        super(StaFile, self).__init__(header)
        self.keys_ = list(header.keys())
        self.id = None

    def __str__(self):
        priorized_keys = [
        'id',
        'net',
        'code',
        'loc',
        'chan',
        'sampling_rate']

        for key in self.keys_:
            if key not in priorized_keys:
                priorized_keys += [key]

        return self._pretty_str(priorized_keys)
    
    def __setid__(self):
        self.id = '%s.%s.%s' % (self.net, self.code, self.loc)


class NetFile(AttribDict):
    def __init__(self, header):
        super(NetFile, self).__init__(header)
        self.keys_ = list(header.keys())

        for s, sta in enumerate(self.stations):
            self.stations[s] = StaFile(sta)
            self.stations[s].sdsdir = self.sdsdir
            self.stations[s].net = self.code
            self.stations[s].__setid__()

        self.stations_info = [sta.id for sta in self.stations]

    def __str__(self):
        priorized_keys = [
        'name',
        'code',
        'sdsdir',
        'stations_info']

        for key in self.keys_:
            if key not in priorized_keys:
                priorized_keys += [key]

        return self._pretty_str(priorized_keys)


class RespFile(AttribDict):
    def __init__(self, header):
        super(RespFile, self).__init__(header)
        self.keys_ = list(header.keys())

        if hasattr(self, 'sensor_keys') and hasattr(self, 'datalogger_keys'):
            self.download()
        
        self.check_times()
        self.load_factor()


    def check_times(self):
        if 'starttime' in self.keys_:
            if self.starttime:
                try:
                    year = int(self.starttime[0])
                    month = int(self.starttime[1])
                    day = int(self.starttime[2])
                    self.starttime = dt.datetime(year, month, day)
                except:
                    print(' warning: fail in defining starttime in RespFile for %s' % self.code)
        
        if 'endtime' in self.keys_:
            if self.endtime:
                try:
                    year = int(self.endtime[0])
                    month = int(self.endtime[1])
                    day = int(self.endtime[2])
                    self.endtime = dt.datetime(year, month, day)    
                except:
                    print(' warning: fail in defining starttime in RespFile for %s' % self.code)
        

    def __str__(self):
        priorized_keys = [
        'code',
        'loc'
        ]

        for key in self.keys_:
            if key not in priorized_keys:
                priorized_keys += [key]
        
        return self._pretty_str(priorized_keys)
    

    def download(self):
        if not os.path.isdir(RESP_PATH):
            os.makedirs(RESP_PATH)
        
        if self.file:
            file_out_name = os.path.join(RESP_PATH, self.file)

            if not os.path.isfile(file_out_name):
                if not all([self.sensor_keys, self.datalogger_keys]):
                    print(' STA : %s.%s >> no sensor/datalogger info in .resp file' % (self.code, self.loc))
                    return

                nrl = NRL()
                resp = nrl.get_response(
                    sensor_keys=self.sensor_keys,
                    datalogger_keys=self.datalogger_keys
                )

                file_out = open(file_out_name, 'wb')
                pickle.dump(resp, file_out, protocol=pickle.HIGHEST_PROTOCOL)
                file_out.close()
        
        else:
            print(' STA : %s.%s >> no file info in .resp file' % (self.code, self.loc))


    def load_factor(self):

        if hasattr(self, 'file'):
            file_out_name = os.path.join(RESP_PATH, self.file)
            if os.path.isfile(file_out_name):
                with open(file_out_name, 'rb') as f:
                    resp = pickle.load(f)

                self.factor = resp.instrument_sensitivity.value
                self.resp  = resp

            return

        self.resp = None
        
        if not hasattr(self, 'factor'):
            self.factor = 1



def get_network(net_code):
    net_dir = os.path.join(NET_PATH, net_code)
    net_file = '%s.net' % net_code
    net_conf = os.path.join(net_dir, net_file)

    if os.path.isdir(net_dir) and os.path.isfile(net_conf):
        loader = SourceFileLoader(net_code, net_conf)
        net = loader.load_module()
        remove_cache(net_dir)
        return NetFile(net.network)
    
    else:
        return None


def get_respfile(net_code, sta, loc):
    net_dir = os.path.join(NET_PATH, net_code)
    resp_file = '%s.resp' % net_code
    resp_path = os.path.join(net_dir, resp_file)

    to_return = []
    if os.path.isdir(net_dir) and os.path.isfile(resp_path):
        loader = SourceFileLoader(net_code, resp_path)
        resp = loader.load_module()
        remove_cache(net_dir)

        for dic in resp.response:
            if dic['code'] == sta and dic['loc'] == loc:
                to_return += [RespFile(dic)]

    if to_return:
        if len(to_return) == 1:
            return to_return[0]
        else:
            return to_return
    else:
        return None


def remove_cache(pathdir):
    cache = '%s/__pycache__' % pathdir
    
    if os.path.isdir(cache):
        rmtree(cache)