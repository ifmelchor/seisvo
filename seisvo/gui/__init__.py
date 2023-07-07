#!/usr/bin/env python
# coding=utf-8

# import os
from matplotlib.backends.qt_compat import QtWidgets
from .events import EventWidget
from .lte import LTEWidget
from .cc8 import CC8Widget

def load_ltewidget(lte, lde, starttime, interval, attr_list, chan_list, olap=0.1, init_app=True, **kwargs):

    if init_app:
        app = QtWidgets.QApplication([])
    
    LTE_widget = LTEWidget(lte, lde, starttime, interval, attr_list, chan_list, olap=0.1, **kwargs)
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


def load_cc8widget(cc8, starttime, interval, fq_idx,\
        olap=0.1, maac_th=0.6, max_err=1, baz_int=[], init_app=True):

    if init_app:
        app = QtWidgets.QApplication([])
    
    cc_widget = CC8Widget(cc8, starttime, interval, fq_idx,\
        olap=olap, maac_th=maac_th, max_err=max_err, baz_int=baz_int)
    cc_widget.show()

    if init_app:
        app.exec_()
    
    return cc_widget


