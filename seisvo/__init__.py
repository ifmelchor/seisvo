#!/usr/bin/python3
# coding=utf-8

import os

try:
    __seisvo__ = os.environ["SEISVO_PATH"]
except KeyError:
    print(' SEISVO_PATH not defined in bashrc. Please, see documentation')
    exit()

seisvo_paths = {
    "networks":os.path.join(__seisvo__, 'networks'),
    "database":os.path.join(__seisvo__, 'database')
}

for name, path in seisvo_paths.items():
    if not os.path.isdir(path):
        os.makedirs(path)
        print(f" info  ::  folder {path} created")
    
    if not os.listdir(path):
        print(f" warning  ::  {name} directory is empty! ")

seisvo_paths["main"] = __seisvo__

from .obspyext import read2 as read
from .network import Network, Station, Array
from .lte import LTE
from .sap.cc8 import CC8
from .database import SDE, LDE


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
