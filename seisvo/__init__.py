#!/usr/bin/python3
# coding=utf-8

import os

try:
    __seisvo__ = os.environ["SEISVO_PATH"]
    DB_PATH = os.path.join(__seisvo__, 'database')
    CC8_PATH = os.path.join(__seisvo__, 'cc8')
    LTE_PATH = os.path.join(__seisvo__, 'lte')
    NET_PATH = os.path.join(__seisvo__, 'networks')
    RESP_PATH = os.path.join(__seisvo__, 'respfiles')
    
except KeyError:
    print('"seisvo_path" not defined in bashrc. Please, see documentation')
    exit()

from .core.network import Network, iArray
from .core.station import Station
from .lte import LTE
from .sap.cc8 import CC8
from .file.air import AiR
from .database import SDE, LDE
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

