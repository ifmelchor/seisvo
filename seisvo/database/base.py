#!/usr/bin/env python3
# coding=utf-8

import os
import sqlalchemy as sql
import datetime as dt
from subprocess import Popen, PIPE, STDOUT

from seisvo.network import Network # iArray
from seisvo.plotting.gui import plot_event
from .events import Event

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

SQLbase = declarative_base()

class DataBase(object):
    def __init__(self, sql_path, type):

        if type in ('SDE', 'LDE', 'iSDE'):
            self.type = type
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
        Popen(call_file, stdout=PIPE, stderr=STDOUT)


    def is_id(self, id):
        return int(id) in self.get_id_list()
        

    def get_id_list(self):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        
        if self.type == 'SDE':
            event_list = session.query(SDErow).all()
        
        if self.type == 'LDE':
            event_list = session.query(LDErow).all()
        
        if self.type == 'iSDE':
            event_list = session.query(iSDErow).all()
        
        id_list = [x.id for x in event_list]
        session.close()

        return id_list
    

    def get_id(self, id=None, with_info=True):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()

        if self.type == 'SDE':
            if id:
                if self.is_id(id):
                    event = session.query(SDErow).filter(SDErow.id == id).all()[0]
                else:
                    event = None
            else:
                if with_info:
                    event = session.query(SDErow).filter(SDErow.starttime != None).all()
                else:
                    event = session.query(SDErow).all()
        
        if self.type == 'LDE':
            if id:
                if self.is_id(id):
                    event = session.query(LDErow).filter(LDErow.id == id).all()[0]
                else:
                    event = None
            else:
                event = session.query(LDErow).all()
        
        if self.type == 'iSDE':
            if id:
                if self.is_id(id):
                    event = session.query(iSDErow).filter(iSDErow.id == id).all()[0]
                else:
                    event = None
            else:
                event = session.query(iSDErow).all()
        
        session.close()
        return event
    

    def get_label_dict(self):
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        
        if self.type == 'SDE':
            events = session.query(SDErow).all()
        
        if self.type == 'LDE':
            events = session.query(LDErow).all()
        
        if self.type == 'iSDE':
            events = session.query(iSDErow).all()
        
        labels = [event.label for event in events if event.label]
        session.close()
        return dict((l,labels.count(l)) for l in set(labels))


    def add_row(self, row, dict_info):
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
        
        if self.type == 'LDE':
            new_evnt = row(
                network=dict_info.get('network'),
                station=dict_info.get('station'),
                location=dict_info.get('location'),
                lte_file=dict_info.get('lte_file'),
                lte_file_sup=dict_info.get('lte_file_sup', None),
                label=dict_info.get('label'),
                starttime=dict_info.get('starttime'),
                duration=dict_info.get('duration') # in hours
                )

        if self.type == 'iSDE':
            new_evnt = row(
                network=dict_info.get('network'),
                station=dict_info.get('station'),
                air_file=dict_info.get('air_file'),
                channels=dict_info.get('channels'),
                label=dict_info.get('label'),
                label2=dict_info.get('label2'),
                starttime=dict_info.get('starttime'),
                duration=dict_info.get('duration'), # in seconds
                pmax=dict_info.get('pmax'),
                pavg=dict_info.get('pavg'),
                azimuth=dict_info.get('azimuth'),
                ccorr=dict_info.get('ccorr')
                )
        
        session.add(new_evnt)
        session.commit()
        session.refresh(new_evnt)
        session.close()

        return new_evnt.id


    def remove_row(self, row, id):
        if self.is_id(id):
            engine = sql.create_engine('sqlite:///%s' % self.sql_path)
            SQLbase.metadata.bind = engine
            DBSession = sessionmaker(bind=engine)
            session = DBSession()
            session.query(row).filter(row.id == id).delete()
            session.commit()
            session.close()
            return True


    def update_row(self, row, id, info_dict):
        if self.is_id(id):
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





