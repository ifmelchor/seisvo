#!/usr/bin/env python3
# coding=utf-8

import os
import h5py
from .base import StationLTE, NetworkLTE

def LTE(lte_file):
    if not os.path.isfile(lte_file):
        return ValueError(' lte_file do not found')

    with h5py.File(lte_file, "r") as f:
        hdr = f['header']
        ltetype = hdr.attrs["type"]

    if ltetype == "station":
        return StationLTE(lte_file)

    elif ltetype == "network":
        return NetworkLTE(lte_file)

    else:
        print("error reading file")
        return None






