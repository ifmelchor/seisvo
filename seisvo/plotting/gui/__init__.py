#!/usr/bin/env python
# coding=utf-8

# import os
from matplotlib.backends.qt_compat import QtWidgets
from .events import EventWidget
from .lte import LTEWidget

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



