#!/usr/bin/python3
# coding=utf-8

import os

try:
    __seisvo__ = os.environ["SEISVO_PATH"]
    LDE_PATH = os.path.join(__seisvo__, 'database', 'lde')
    LTE_PATH = os.path.join(__seisvo__, 'lte')
    NET_PATH = os.path.join(__seisvo__, 'networks')
    RESP_PATH = os.path.join(__seisvo__, 'respfiles')
    
except KeyError:
    print('"seisvo_path" not defined in bashrc. Please, see documentation')
    exit()

from .core.network import Network, iArray
from .core.station import Station
from .file.lte import LTE
from .file.air import AiR
from .database import SDE, LDE
from .signal import SSteps


