#!/usr/bin/python3
# coding=utf-8

import os

try:
    __seisvo__ = os.environ["seisvo_path"]
    default_LDE_dir = os.path.join(__seisvo__, 'database', 'lde')
    default_LTE_dir = os.path.join(__seisvo__, 'lte')
except KeyError:
    print('"seisvo_path" not defined in bashrc. Please, see documentation')
    exit()

from .core.network import Network, iArray
from .core.station import Station
from .file.lte import LTE
from .file.air import AiR
from .database import SDE, LDE


