#!/usr/bin/python3
# coding=utf-8

__all__ = ['get', 'read', 'Network', 'Station', 'LTE', 'CC8', 'SDE', 'LDE', 'seisvo_paths']

import os

try:
    __seisvo__ = os.environ["SEISVO_PATH"]
except KeyError:
    print(' SEISVO_PATH not defined in bashrc. Please, see documentation')
    exit()

seisvo_paths = {
    "networks":os.path.join(__seisvo__, 'networks'),
    "database":os.path.join(__seisvo__, 'database'),
    "cc8":os.path.join(__seisvo__, 'files', 'cc8'),
    "lte":os.path.join(__seisvo__, 'files', 'lte')
}

for _, path in seisvo_paths.items():
    if not os.path.isdir(path):
        os.makedirs(path)
        print(f" >>> folder {path} created")

seisvo_paths["main"] = __seisvo__

from .obspyext import read2 as read
from .network import Network, Station, SeismicArray #,iArray
from .lte import LTE
from .sap.cc8 import CC8
from .database import SDE, LDE
# from .file.air import AiR
from .signal import SSteps

def get(sarg):
    argin = sarg.split("/")
    net_code = argin[0]

    try:
        net = Network(net_code)
    except Exception as ex:
        print("error [reading network] :: ", str(ex))
        return

    if len(argin) == 2:
        try:
            stacode = argin[1].split(".")
            sta = stacode[0]

            if len(stacode) == 2:
                loc = stacode[1]
            else:
                loc = ''
                
            sta = net.get_sta(sta, loc=loc)
        
        except Exception as ex:
            print("error [reading station] :: ", str(ex))
            return net

        return sta
    
    else:
        return net
