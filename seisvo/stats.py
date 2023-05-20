#!/usr/bin/python3
# coding=utf-8


import os
import pickle
from obspy.core.util.attribdict import AttribDict
from obspy.core.inventory.response import Response


class StationStats(AttribDict):
    def __init__(self, network, header):
        super(StationStats, self).__init__(header)
        self.keys_ = list(header.keys())
        self.__load_location__()
        self.id = '.'.join([network.code, self.code, self.location])
        self.network = network


    def __str__(self):
        priorized_keys = [
            'id',
            'code',
            'location',
            'channels',
            'sample_rate'
        ]

        return self._pretty_str(priorized_keys)
    

    def __load_location__(self):
        if "location" in self.keys_:
            if not self.location:
                self.location = ""


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



class NetworkStats(AttribDict):
    def __init__(self, code, header):
        super(NetworkStats, self).__init__(header)
        self.keys_ = list(header.keys())
        self.code = code

        # buid station list
        station_stats = []
        for _, sta in self.stations.items():
            station_stats.append(StationStats(self, sta))        
        self.stations = station_stats

        self.stations_id = [sta.id for sta in self.stations]


    def __str__(self):

        text =  f"\n   {self.name} network  |  code: {self.code}  |  path: {self.sds_path} ---"
        text += f"\n   Stations  >>  {len(self.stations)} {self.stations_id}"

        return text

