#!/usr/bin/env python3
# coding=utf-8

import datetime as dt
from sqlalchemy.ext.mutable import MutableList
# from seisvo.network import iArray
from seisvo.utils import in_interval
from .events import iEvent
from .base import sql, DataBase, SQLbase

class iSDErow(SQLbase):
    __tablename__ = 'iSDE'

    id = sql.Column(sql.Integer, primary_key=True)
    network = sql.Column(sql.String, nullable=False)
    station = sql.Column(sql.String, nullable=False)
    
    label = sql.Column(sql.String, nullable=False)
    label2 = sql.Column(sql.String, nullable=True)
    
    starttime = sql.Column(sql.DateTime(timezone=False), nullable=False)
    duration = sql.Column(sql.Float, nullable=False)

    pmax = sql.Column(sql.Float, nullable=False)
    pavg = sql.Column(sql.Float, nullable=False)
    azimuth = sql.Column(sql.Float, nullable=False)
    ccorr = sql.Column(sql.Float, nullable=False)
    channels = sql.Column(MutableList.as_mutable(sql.PickleType), nullable=False)

    air_file = sql.Column(sql.String, nullable=False)

    # additional floats
    value_1 = sql.Column(sql.Float, nullable=True)
    value_2 = sql.Column(sql.Float, nullable=True)
    value_3 = sql.Column(sql.Float, nullable=True)
    value_4 = sql.Column(sql.Float, nullable=True)
    value_5 = sql.Column(sql.Float, nullable=True)

    # additional strings
    string_1 = sql.Column(sql.String, nullable=True)
    string_2 = sql.Column(sql.String, nullable=True)
    string_3 = sql.Column(sql.String, nullable=True)
    string_4 = sql.Column(sql.String, nullable=True)
    string_5 = sql.Column(sql.String, nullable=True)

    def get_array(self, **kwargs):
        return iArray(self.network, self.station, **kwargs)


class iSDE(DataBase):
    def __init__(self, sql_path):
        super().__init__(sql_path, 'iSDE')
    
    def __add_event__(self, net_code, sta_code, label, starttime, duration, air_file, pmax, pavg, channels):
        """
        Add a new event in iSDE database
        """

        event_to_save = {}
        event_to_save['network'] = net_code
        event_to_save['station'] = sta_code
        event_to_save['label'] = label
        event_to_save['channels'] = channels
        event_to_save['starttime'] = starttime #datetime
        event_to_save['duration'] = duration
        event_to_save['air_file'] = air_file
        event_to_save['pmax'] = pmax
        event_to_save['pavg'] = pavg

        self.add_row(iSDErow, event_to_save)


    def __getitem__(self, eid):
        return iEvent(eid, self)
    

    def relabel_event(self, id, new_label):
        if self.is_id(id):
            info = dict(label=new_label)
            self.update_row(iSDErow, id, info)
    

    def get_events_id(self, label=None, time_interval=()):
        row_list = self.get_id()

        if label:
            if isinstance(label, str):
                row_list = list(filter(lambda e: e.label==label, row_list))
            
            elif isinstance(label, list):
                row_list_copy = []
                for lbl in label:
                    row_list_copy += list(filter(lambda e: e.label==lbl, row_list))
                
                row_list = row_list_copy
        
        if time_interval:
            row_list = list(filter(
                lambda e: in_interval(
                    e.starttime,
                    e.starttime+dt.timedelta(seconds=e.duration),
                    time_interval
                    ),row_list))

        row_list.sort(key=lambda x: x.starttime)
        row_list = [e.id for e in row_list]

        return row_list

