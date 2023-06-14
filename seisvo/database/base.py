#!/usr/bin/env python3
# coding=utf-8

import os
import datetime as dt
import sqlalchemy as sql
import subprocess as sp
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

from .rows import SDErow
from .events import Event
from ..plotting.gui import plot_event
# from seisvo.network import Network # iArray

SQLbase = declarative_base()

class _DataBase(object):
    def __init__(self, sql_path, dbtype):

        if dbtype in ('SDE', 'LDE', 'iSDE'):
            self.type = dbtype
        else:
            raise ValueError(' DataBase must be SDE or LDE')

        if not os.path.isfile(sql_path):
            engine = sql.create_engine('sqlite:///%s' % sql_path)
            SQLbase.metadata.create_all(engine)
            print(' Database %s created' % sql_path)
        else:
            # check is database coincides with type
            print(' Database %s ready' % sql_path)

        self.sql_path = sql_path
        

    def open(self):
        call_file = ["sqlitebrowser", self.sql_path]
        sp.Popen(call_file, stdout=sp.PIPE, stderr=sp.STDOUT)


    def _is_id(self, id, row):
        return int(id) in self._get_id_list(row)
        

    def _get_id_list(self, row):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        event_list = session.query(row).all()
        id_list = [x.id for x in event_list]
        session.close()
        return id_list
        

    def _get_id(self, row, id=None, with_info=True):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        if id:
            if self._is_id(id):
                event = session.query(row).filter(row.id == id).all()[0]
            else:
                event = None
        else:
            if with_info:
                event = session.query(row).filter(row.starttime != None).all()
            else:
                event = session.query(row).all()
        session.close()
        return event
    

    def _get_label_dict(self, row):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        events = session.query(row).all()
        labels = [event.label for event in events if event.label]
        session.close()
        return dict((l,labels.count(l)) for l in set(labels))


    def _add_row(self, row, dict_info):
        """
        Save a new row using a dict and return the ID of the saved row
        """

        # open db an add new event
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()

        if self.type == 'SDE':
            new_evnt = row(
                network=dict_info.get('network'),
                station=dict_info.get('station'),
                location=dict_info.get('location'),
                label=dict_info.get('label'),
                starttime=dict_info.get('starttime'),
                duration=dict_info.get('duration'), # in seconds
                event_type=dict_info.get('event_type'),
                event_id=dict_info.get('event_id')
                )
        
        # if self.type == 'LDE':
        #     new_evnt = row(
        #         network=dict_info.get('network'),
        #         station=dict_info.get('station'),
        #         location=dict_info.get('location'),
        #         lte_file=dict_info.get('lte_file'),
        #         lte_file_sup=dict_info.get('lte_file_sup', None),
        #         label=dict_info.get('label'),
        #         starttime=dict_info.get('starttime'),
        #         duration=dict_info.get('duration') # in hours
        #         )

        # if self.type == 'iSDE':
        #     new_evnt = row(
        #         network=dict_info.get('network'),
        #         station=dict_info.get('station'),
        #         air_file=dict_info.get('air_file'),
        #         channels=dict_info.get('channels'),
        #         label=dict_info.get('label'),
        #         label2=dict_info.get('label2'),
        #         starttime=dict_info.get('starttime'),
        #         duration=dict_info.get('duration'), # in seconds
        #         pmax=dict_info.get('pmax'),
        #         pavg=dict_info.get('pavg'),
        #         azimuth=dict_info.get('azimuth'),
        #         ccorr=dict_info.get('ccorr')
        #         )
        
        session.add(new_evnt)
        session.commit()
        session.refresh(new_evnt)
        session.close()

        return new_evnt.id


    def _remove_row(self, row, id):
        if self._is_id(id):
            engine = sql.create_engine('sqlite:///%s' % self.sql_path)
            SQLbase.metadata.bind = engine
            DBSession = sessionmaker(bind=engine)
            session = DBSession()
            session.query(row).filter(row.id == id).delete()
            session.commit()
            session.close()
            return True


    def _update_row(self, row, id, info_dict):
        if self._is_id(id):
            engine = sql.create_engine('sqlite:///%s' % self.sql_path)
            SQLbase.metadata.bind = engine
            DBSession = sessionmaker(bind=engine)
            session = DBSession()

            for item, key in info_dict.items():
                session.query(row).filter(row.id == id).update(
                {item : key}, synchronize_session=False)
                session.commit()
            
            session.close()
            return True


class SDE(_DataBase):
    def __init__(self, sql_path):
        super().__init__(sql_path, 'SDE')
        self.__set_attr__()
        print(self)


    def __set_attr__(self):
        eid_list = self.get_event_list() 
        
        if eid_list:
            self.first_eid = eid_list[0]
            self.last_eid  = eid_list[-1]
        else:
            print(" [info]  Empty database")
            self.first_eid = None
            self.last_eid = 0

    
    def __len__(self):
        return len(self.get_event_list())


    def __getitem__(self, eid):
        if eid in self.get_event_list():
            return Event(eid, self)
        else:
            print(f" eid {eid} not exist")
            return None


    def __str__(self):
        text = f"  SQL Database  [{self.type}] ::  {self.sql_path} "
        text += f"\n     nro. events: {len(self)}"
        return text
        

    def get_id(self, id=None, with_info=True):
        return self._get_id(SDErow, id=id, with_info=with_info)


    def add_event(self, station_id_list, label, starttime, duration, etype="S"):
        """
        Add a new event in SDE database, for adding a row visit database/__init__ info
        Check database atributes por kwargs
        """

        eid = self.last_eid + 1

        for station_id in list(set(station_id_list)):
            event_to_save = {
                'network'   : station_id.split(".")[0],
                'station'   : station_id.split(".")[1],
                'location'  : station_id.split(".")[2],
                'event_type': etype,
                'label'     : label,
                'starttime' : starttime,
                'duration'  : duration,
                'event_id'  : eid,
            }
            id = self._add_row(SDErow, event_to_save)
        
        self.last_eid += 1
    

    def get_event_list(self, label=None, time_interval=(), nro_station=None):
        row_list = self.get_id()

        if label:
            if isinstance(label, str):
                row_list = list(filter(lambda e: e.label==label, row_list))
            
            if isinstance(label, list):
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
        eid_list = [e.event_id for e in row_list]

        if nro_station:
            out_list = []
            for eid in eid_list:
                if len(self[eid]) >= nro_station:
                    out_list.append(eid)
            return out_list
        
        return list(set(eid_list))
    

    def is_eid(self, eid):
        return eid in self.get_event_list()


    def get_event(self, eid=None):
        """
        Rreturn a list of events associated to event_id, 
        else return a dict of lists.
        """

        eid_list = self.get_event_list()
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        
        if eid:
            if isinstance(eid, int) and self.is_eid(eid):
                event = session.query(SDErow).filter(SDErow.event_id == eid).all()

            if isinstance(eid, (list, tuple)):
                event = {}
                for eid in eid_list:
                    if self.is_eid(eid):
                        event[eid] = session.query(SDErow).filter(SDErow.event_id == eid).all()
        
        else:
            event = {}
            for eid in eid_list:
                event[eid] = session.query(SDErow).filter(SDErow.event_id == eid).all()
        
        session.close()

        return event


    def remove_event(self, eid):
        if isinstance(eid, int) and self.is_eid(eid):
            for row in self.get_event(eid):
                self._remove_row(SDErow, row.id)
        
        if isinstance(eid, (list, tuple)):
            for eidi, row_list in self.get_event(eid).items():
                if self.is_eid(eidi):
                    for row in row_list:
                        self._remove_row(SDErow, row.id)
        
        self.__set_attr__()
    

    def relabel_event(self, eid, new_label):
        assert self.is_eid(eid)
        assert isinstance(new_label, str)

        event = self.get_event(eid)
        for row in event:
            self._update_row(SDErow, row.id, dict(label=new_label))


    def plot_gui(self, event_id=None, station_id=None, init_app=True, **kwargs):
        eid_list = self.get_event_list(**kwargs)

        if event_id:
            assert event_id in eid_list
        
        event_id = eid_list[0]
        widget = plot_event(self, event_id, station_id, init_app=init_app)

        return widget


    def write_hypo71(self, min_phase=3, out_file='./h71.phase'):
        """
        Create a phase file in HYPO71 format
        """

        f = open(out_file, 'w')
        try:
            for event in self:
                if event:
                    if event.nro_phase >= min_phase:
                        event.to_hypo71(phsfile=f)
                        f.write("\n")
            f.close()
            print(f" OK :: file {out_file} created!")
        
        except Exception as e:
            f.close()
            print(e)
            print(" FAIL :: error found")





