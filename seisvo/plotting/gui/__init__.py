#!/usr/bin/env python
# coding=utf-8

import os
from pynotifier import Notification
from matplotlib.backends.qt_compat import QtWidgets
from .events import EventWidget

path = os.path.dirname(os.path.realpath(__file__))
icons_path = os.path.join(path, 'icons')


def plot_event(sde, event_id, station_id, init_app=True):

    if init_app:
        app = QtWidgets.QApplication([])
    
    SDE_widget = EventWidget(sde, event_id=event_id, station_id=station_id)
    SDE_widget.show()

    if init_app:
        app.exec_()
    
    return SDE_widget


def notify(title, description, status='info', duration=2):
    """
    status: error, info, warn
    """
    n = Notification(
        title=title, 
        description=description, 
        icon_path='%s/%s.png' %(icons_path, status),
        duration=duration
        )
    n.send()


