#!/usr/bin/env python
# coding=utf-8

# import os
from matplotlib.backends.qt_compat import QtWidgets
from .events import EventWidget
from .lte import LTEWidget
from .cc8 import CC8Widget
from .cce import CCEWidget
from .station import StationMainWidget

def load_ltewidget(lte, db, starttime, interval, attr_list, chan_list, olap=0.1, init_app=True, **kwargs):

    if init_app:
        app = QtWidgets.QApplication([])
    
    LTE_widget = LTEWidget(lte, db, starttime, interval, attr_list, chan_list, olap=0.1, **kwargs)
    LTE_widget.show()

    if init_app:
        app.exec_()
    
    return LTE_widget


def load_eventwidget(sde, event_id, station_id, init_app=True):

    if init_app:
        app = QtWidgets.QApplication([])
    
    SDE_widget = EventWidget(sde, event_id=event_id, station_id=station_id)
    SDE_widget.show()

    if init_app:
        app.exec_()
    
    return SDE_widget


def load_cc8widget(cc8, starttime, interval, fq_idx, db, init_app=True, **kwargs):

    if init_app:
        app = QtWidgets.QApplication([])
    
    cc_widget = CC8Widget(cc8, starttime, interval, fq_idx, db, **kwargs)
    cc_widget.show()

    if init_app:
        app.exec_()
    
    return cc_widget


def load_cceWidget(cce, eid_list, path_to_cc8file="./"):

    app = QtWidgets.QApplication([])
    
    CCE_widget = CCEWidget(cce, eid_list, path_to_cc8file=path_to_cc8file)
    CCE_widget.show()
    
    app.exec_()

    return CCE_widget


def load_stationwidget(network, station, starttime, interval=60, olap=0.1, sde_db=None, spec_v=[10,45], init_app=True, **sta_kwargs):

    if init_app:
        app = QtWidgets.QApplication([])
    
    widget = StationMainWidget(network, station, starttime, interval, olap, sde_db, spec_v, **sta_kwargs)
    widget.show()

    if init_app:
        app.exec_()
    
    return widget