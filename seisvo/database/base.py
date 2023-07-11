#!/usr/bin/env python3
# coding=utf-8

import os
import datetime as dt
import sqlalchemy as sql
import subprocess as sp
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

from .rows import get_sderow, get_lderow, get_ccerow
from .events import Event
from ..utils import in_interval

SQLbase = declarative_base()

class _DataBase(object):
    def __init__(self, sql_path, dbtype, new):

        if dbtype in ('SDE', 'LDE', 'CCE'):
            self.type = dbtype
            
            if self.type == "SDE":
                self.sql_row = get_sderow(SQLbase)
            
            if self.type == "LDE":
                self.sql_row = get_lderow(SQLbase)

            if self.type == "CCE":
                self.sql_row = get_ccerow(SQLbase)
        
        else:
            raise ValueError(' DataBase must be SDE or LDE')

        if new:
            engine = sql.create_engine('sqlite:///%s' % sql_path)
            SQLbase.metadata.create_all(engine)
            print(' Database %s created' % sql_path)
        
        else:
            print(' Database %s ready' % sql_path)

        self.sql_path = sql_path
        

    def open(self):
        call_file = ["sqlitebrowser", self.sql_path]
        sp.Popen(call_file, stdout=sp.PIPE, stderr=sp.STDOUT)


    def _is_id(self, id):
        return int(id) in self._get_id_list()
        

    def _get_id_list(self):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        event_list = session.query(self.sql_row).all()
        id_list = [x.id for x in event_list]
        session.close()
        return id_list
        

    def _get_id(self, id=None, with_info=True):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        if id:
            if self._is_id(id):
                event = session.query(self.sql_row).filter(self.sql_row.id == id).all()[0]
            else:
                event = None
        else:
            if with_info:
                event = session.query(self.sql_row).filter(self.sql_row.starttime != None).all()
            else:
                event = session.query(self.sql_row).all()
        session.close()
        return event
    

    def _get_label_dict(self):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        events = session.query(self.sql_row).all()
        labels = [event.label for event in events if event.label]
        session.close()
        return dict((l,labels.count(l)) for l in set(labels))


    def _add_row(self, dict_info):
        """
        Save a new row using a dict and return the ID of the saved row
        """

        # open db an add new event
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()

        if self.type == 'SDE':
            new_evnt = self.sql_row(
                network=dict_info.get('network'),
                station=dict_info.get('station'),
                location=dict_info.get('location'),
                label=dict_info.get('label'),
                starttime=dict_info.get('starttime'),
                duration=dict_info.get('duration'), # in seconds
                event_type=dict_info.get('event_type'),
                event_id=dict_info.get('event_id')
                )
        
        if self.type == 'LDE':
            new_evnt = self.sql_row(
                network=dict_info.get('network'),
                station=dict_info.get('station'),
                location=dict_info.get('location'),
                lte_file=dict_info.get('lte_file'),
                label=dict_info.get('label'),
                starttime=dict_info.get('starttime'),
                duration=dict_info.get('duration') # in minutes
                )

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


    def _remove_row(self, id):
        if self._is_id(id):
            engine = sql.create_engine('sqlite:///%s' % self.sql_path)
            SQLbase.metadata.bind = engine
            DBSession = sessionmaker(bind=engine)
            session = DBSession()
            session.query(self.sql_row).filter(self.sql_row.id == id).delete()
            session.commit()
            session.close()
            return True


    def _update_row(self, id, info_dict):
        if self._is_id(id):
            engine = sql.create_engine('sqlite:///%s' % self.sql_path)
            SQLbase.metadata.bind = engine
            DBSession = sessionmaker(bind=engine)
            session = DBSession()

            for item, key in info_dict.items():
                session.query(self.sql_row).filter(self.sql_row.id == id).update(
                {item : key}, synchronize_session=False)
                session.commit()
            
            session.close()
            return True


class SDE(_DataBase):
    def __init__(self, sql_path):
        if not os.path.isfile(sql_path):
            create = True
        else:
            create = False

        super().__init__(sql_path, 'SDE', new=create)
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
        return self._get_id(id=id, with_info=with_info)


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
            id = self._add_row(event_to_save)
        
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
                event = session.query(self.sql_row).filter(self.sql_row.event_id == eid).all()

            if isinstance(eid, (list, tuple)):
                event = {}
                for eid in eid_list:
                    if self.is_eid(eid):
                        event[eid] = session.query(self.sql_row).filter(self.sql_row.event_id == eid).all()
        
        else:
            event = {}
            for eid in eid_list:
                event[eid] = session.query(self.sql_row).filter(self.sql_row.event_id == eid).all()
        
        session.close()

        return event


    def remove_event(self, eid):
        if isinstance(eid, int) and self.is_eid(eid):
            for row in self.get_event(eid):
                self._remove_row(row.id)
        
        if isinstance(eid, (list, tuple)):
            for eidi, row_list in self.get_event(eid).items():
                if self.is_eid(eidi):
                    for row in row_list:
                        self._remove_row(row.id)
        
        self.__set_attr__()
    

    def relabel_event(self, eid, new_label):
        assert self.is_eid(eid)
        assert isinstance(new_label, str)

        event = self.get_event(eid)
        for row in event:
            self._update_row(row.id, dict(label=new_label))


    def plot_gui(self, event_id=None, station_id=None, init_app=True, **kwargs):
        eid_list = self.get_event_list(**kwargs)

        if event_id:
            assert event_id in eid_list
        
        event_id = eid_list[0]

        from ..gui import load_eventwidget
        
        widget = load_eventwidget(self, event_id, station_id, init_app=init_app)

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


class LDE(_DataBase):
    def __init__(self, sql_path):
        if not os.path.isfile(sql_path):
            create = True
        else:
            create = False

        super().__init__(sql_path, 'LDE', new=create)


    def __len__(self):
        return len(self.get_id_list())


    def __getitem__(self, eid):
        return self.get_id(eid)


    def get_id(self, id=None):
        return self._get_id(id=id)


    def add_episode(self, row_dict):
        episode_to_save = {}
        episode_to_save['network']   = row_dict.get("network")
        episode_to_save['station']   = row_dict.get("station")
        episode_to_save['location']  = row_dict.get("location")
        episode_to_save['starttime'] = row_dict.get("starttime") #datetime
        episode_to_save['duration']  = row_dict.get("duration")
        episode_to_save['label']     = row_dict.get("label")
        episode_to_save['lte_file']  = row_dict.get("lte_file")
        return self._add_row(episode_to_save)
        

    def relabel_event(self, eid, new_label):
        assert self.is_eid(eid)
        assert isinstance(new_label, str)
        self._update_row(eid, dict(label=new_label))


    def get_episode_list(self, label=None, time_interval=()):
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
    
    
    def is_eid(self, eid):
        return eid in self.get_episode_list()
    

    def remove_event(self, eid):
        if isinstance(eid, int) and self.is_eid(eid):
            self._remove_row(eid)
        
        if isinstance(eid, (list, tuple)):
            for eidi in eid:
                if self.is_eid(eidi):
                    self._remove_row(eid)
    

class CCE(_DataBase):
    def __init__(self, sql_path):
        
        if not os.path.isfile(sql_path):
            create = True
        else:
            create = False

        super().__init__(sql_path, 'CCE', new=create)


