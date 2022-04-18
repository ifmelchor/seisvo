#!/usr/bin/python3
# coding=utf-8

import os

try:
    __seisvo__ = os.environ["seisvo_path"]

except KeyError:
    print('"seisvo_path" not defined in bashrc. Please, see documentation')
    exit()
