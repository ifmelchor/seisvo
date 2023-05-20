#!/usr/bin/env python3
# coding=utf-8

import os
import datetime as dt
from seisvo.network import Network
from seisvo.utils import in_interval
from .events import Episode
from .base import sql, DataBase, SQLbase


class LDErow(SQLbase):
    __tablename__ = 'LDE'
    id = sql.Column(sql.Integer, primary_key=True)

    network = sql.Column(sql.String, nullable=False)
    station = sql.Column(sql.String, nullable=False)
    location = sql.Column(sql.String, nullable=True)
    
    label = sql.Column(sql.String, nullable=False)
    starttime = sql.Column(sql.DateTime(timezone=False), nullable=False)
    duration = sql.Column(sql.Float, nullable=False)
    lte_file = sql.Column(sql.String, nullable=False) # lte file containing the event info
    
    # lte file supplementary
    lte_file_sup = sql.Column(sql.String, nullable=True)
    
    def get_network(self):
        return Network(self.network)

    def get_station(self):
        net = self.get_network()
        return net.get_sta(self.station, self.location)



class LDE(DataBase):
    def __init__(self, sql_path):
        super().__init__(sql_path, 'LDE')
        self.id = os.path.basename(sql_path).split('.lde')[0]

    
    def __len__(self):
        return len(self.get_id_list())


    def __getitem__(self, eid):
        return self.get(eid)
    

    def get(self, eid):
        return Episode(eid, self)


    def __add_episode__(self, net_code, sta_code, location, label, starttime, duration, lte_file, lte_file_sup=None):
        """
        Add a new episode in LDE database, for adding a row visit database/__init__ info
        Check database atributes por kwargs
        """

        event_to_save = {}
        event_to_save['network'] = net_code
        event_to_save['station'] = sta_code
        event_to_save['location'] = location
        
        event_to_save['label'] = label
        event_to_save['starttime'] = starttime #datetime
        event_to_save['duration'] = duration
        event_to_save['lte_file'] = lte_file
        event_to_save['lte_file_sup'] = lte_file_sup

        self.add_row(LDErow, event_to_save)


    def relabel_row(self, id, new_label):
        if self.is_id(id):
            info = dict(label=new_label)
            self.update_row(LDErow, id, info)
    

    def get_episodes_id(self, label=None, time_interval=()):
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
                    e.starttime+dt.timedelta(hours=e.duration),
                    time_interval
                    ),row_list))

        row_list.sort(key=lambda x: x.starttime)
        row_list = [e.id for e in row_list]

        return row_list
    

    def __update_lte_sup__(self, id, lte_sup_path):
        if self.is_id(id):
            info = dict(lte_file_sup=lte_sup_path)
            self.update_row(LDErow, id, info)
        
