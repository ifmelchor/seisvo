#!/usr/bin/env python3
# coding=utf-8

"""
    This code creates the structure of the database
    and define the SDE and LDE databases

    warn!! LDE is still unlisted!
"""

from seisvo import __seisvo__
from seisvo.core.network import Network, iArray
from events import Event, iEvent, Episode

import sqlalchemy as sql
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.mutable import MutableList
from sqlalchemy.orm import sessionmaker
from subprocess import Popen, PIPE, STDOUT
import datetime as dt
import os

def in_interval(st, et, time_interval):
    cond1 = st >= time_interval[0] and et <= time_interval[1] # start and end in
    cond2 = time_interval[1] > st > time_interval[0] and et > time_interval[1] # start in, end out
    cond3 = time_interval[0] > st and time_interval[0] < et < time_interval[1] # end out, start in
    cond4 = st < time_interval[0] and time_interval[1] < et # start and end out
    return cond1 or cond2 or cond3 or cond4

SQLbase = declarative_base()

class SDErow(SQLbase):
    __tablename__ = 'SDE'

    # base attributes
    id = sql.Column(sql.Integer, primary_key=True)
    event_id = sql.Column(sql.Integer, nullable=False)
    event_type = sql.Column(sql.String, nullable=False)
    
    network = sql.Column(sql.String, nullable=False)
    station = sql.Column(sql.String, nullable=False)
    location = sql.Column(sql.String, nullable=True)
    
    label = sql.Column(sql.String, nullable=True)
    sublabel = sql.Column(sql.String, nullable=True)
    
    starttime = sql.Column(sql.DateTime(timezone=False), nullable=True)
    duration = sql.Column(sql.Float, nullable=True)

    # event attributes
    onset_P = sql.Column(sql.String, nullable=True)
    first_motion_P = sql.Column(sql.String, nullable=True)
    time_P = sql.Column(sql.Float, nullable=True)
    weight_P = sql.Column(sql.Integer, nullable=True)
    time_S = sql.Column(sql.Float, nullable=True)
    weight_S = sql.Column(sql.Integer, nullable=True)
    peak_to_peak = sql.Column(sql.Float, nullable=True)
    max_period = sql.Column(sql.Float, nullable=True)

    # others
    event_amplitude = sql.Column(sql.Float, nullable=True)
    event_duration = sql.Column(sql.Float, nullable=True)
    event_frequency = sql.Column(sql.Float, nullable=True)
    azimuth = sql.Column(sql.Float, nullable=True)
    incidence = sql.Column(sql.Float, nullable=True)

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

    def __str__(self):
        text_info = " EID %i || ID %i | %s | %s.%s.%s " % (
            self.event_id, self.id, self.event_type, 
            self.network, self.station, self.location)
        return text_info


    def get_network(self):
        return Network(self.network)


    def get_station(self):
        net = self.get_network()
        return net.get_sta(self.station, self.location)


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


class _DataBase(object):
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
        if isinstance(id, int):
            return id in self.get_id_list()
        else:
            print(' id must be an integer')
            

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


    def add_row(self, dict_info):
        """
        Save a new row using a dict and return the ID of the saved row
        """

        # open db an add new event
        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()

        if self.type == 'SDE':
            new_evnt = SDErow(
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
            new_evnt = LDErow(
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
            new_evnt = iSDErow(
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


    def remove_row(self, id):
        if self.is_id(id):
            engine = sql.create_engine('sqlite:///%s' % self.sql_path)
            SQLbase.metadata.bind = engine
            DBSession = sessionmaker(bind=engine)
            session = DBSession()
            
            if self.type == 'SDE':
                session.query(SDErow).filter(SDErow.id == id).delete()
                session.commit()
            
            if self.type == 'LDE':
                session.query(LDErow).filter(LDErow.id == id).delete()
                session.commit()
            
            if self.type == 'iSDE':
                session.query(iSDErow).filter(iSDErow.id == id).delete()
                session.commit()

            session.close()
            return True


    def update_row(self, id, info_dict):
        if self.is_id(id):
            engine = sql.create_engine('sqlite:///%s' % self.sql_path)
            SQLbase.metadata.bind = engine
            DBSession = sessionmaker(bind=engine)
            session = DBSession()

            if self.type == 'SDE':
                for item, key in info_dict.items():
                    session.query(SDErow).filter(SDErow.id == id).update(
                    {item : key}, synchronize_session=False)
                    session.commit()
            
            if self.type == 'LDE':
                for item, key in info_dict.items():
                    session.query(LDErow).filter(LDErow.id == id).update(
                    {item : key}, synchronize_session=False)
                    session.commit()
            
            if self.type == 'iSDE':
                for item, key in info_dict.items():
                    session.query(iSDErow).filter(iSDErow.id == id).update(
                    {item : key}, synchronize_session=False)
                    session.commit()

            session.close()


class SDE(_DataBase):
    def __init__(self, sql_path):
        super().__init__(sql_path, 'SDE')
        self.__check__()

    
    def __len__(self):
        eid_len = len(self.get_eid_list())
        return eid_len


    def __getitem__(self, eid):
        return Event(eid, self)


    def __check__(self):
        # this function checks the database is ready for operate
        for eid in self.get_eid_list():
            stations = Event(eid, self).stations
            if len(stations) != len(set(stations)):
                # it exist repeated items
                raise ValueError(' Duplicate stations with same Event ID found.\n SDE database should be revised manually.')
            # correct automatically


    def __add_event__(self, net_code, sta_code, location, label, starttime, duration, **kwargs):
        """
        Add a new event in SDE database, for adding a row visit database/__init__ info
        Check database atributes por kwargs
        """

        event_to_save = {}
        event_to_save['network'] = net_code
        event_to_save['station'] = sta_code
        event_to_save['location'] = location

        if self.is_infrasound():
            event_to_save['event_type'] = 'P'
        else:
            event_to_save['event_type'] = 'S'
        
        event_to_save['label'] = label
        event_to_save['starttime'] = starttime #datetime
        event_to_save['duration'] = duration
        event_to_save['event_id'] = self.last_eid() + 1

        id = self.add_row(event_to_save)
        self.update_row(id, kwargs)
    

    def get_eid_list(self, label=None, time_interval=(), nro_station=None):
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
        row_list = [e.event_id for e in row_list]

        if nro_station:
            out_list = []
            for eid in row_list:
                if len(self[eid]) >= nro_station:
                    out_list.append(eid)
            return out_list
        
        else:
            return row_list
    

    def last_eid(self):
        if self.get_eid_list():
            return max(self.get_eid_list())
        else:
            return 0


    def is_eid(self, eid):
        if isinstance(eid, int):
            return eid in self.get_eid_list()
        
        else:
            print(' event_id must be an integer')
    

    def get_row(self, eid, network, station, location):
        rows = self.get_id(with_info=False)
        events_filt = list(filter(lambda e: e.event_id == eid and e.network == network and e.station == station and e.location == location, rows))
        if events_filt:
            events_filt = [e.id for e in events_filt]
            if len(events_filt) > 1:
                print( 'warn: more than one id found')
            return events_filt[0]
        else:
            return None


    def get_event(self, eid=None):
        """
        Rreturn a list of events associated to event_id, 
        else return a dict of lists.
        """

        eid_list = self.get_eid_list()

        engine = sql.create_engine('sqlite:///%s' % self.sql_path)
        SQLbase.metadata.bind = engine
        DBSession = sessionmaker(bind=engine)
        session = DBSession()
        if eid:
            if self.is_eid(eid):
                event = session.query(SDErow).filter(SDErow.event_id == eid).all()
        else:
            event = {}
            for eid in eid_list:
                event[eid] = session.query(SDErow).filter(SDErow.event_id == eid).all()
        session.close()
        return event


    def relabel_event(self, eid, new_label):
        if self.is_eid(eid):
            for e in self.get_event(eid):
                if e.label:
                    info = dict(label=new_label)
                    self.update_row(e.id, info)


    def append_event(self, eid, info):
        if self.is_eid(eid):
            id = self.add_row(info)
            return id
    

    def remove_event(self, eid):
        if self.is_eid(eid):
            for e in self.get_event(eid):
                self.remove_row(e.id)


    def plot_gui(self, eid=None, label=None, time_interval=(), app=False, **gui_kwargs):
        from seisvo.gui.gsde import plot_SDE_event

        eid_list = self.get_eid_list(label=label, time_interval=time_interval)

        if eid_list:
            # gui creation
            window = plot_SDE_event(self, eid_list, eid=eid, app=app, **gui_kwargs)

            if app:
                return window
        
        else:
            print(' Empty EID list')


class LDE(_DataBase):
    def __init__(self, sql_path):
        super().__init__(sql_path, 'LDE')
        self.id = os.path.basename(sql_path).split('.lde')[0]

    
    def __len__(self):
        return len(self.get_id_list())


    def __getitem__(self, eid):
        self.get(eid)
    

    def get(self, eid):
        return Episode(eid, self)


    def __add_episode__(self, net_code, sta_code, location, label, starttime, duration, lte_file):
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

        self.add_row(event_to_save)


    def relabel_row(self, id, new_label):
        if self.is_id(id):
            info = dict(label=new_label)
            self.update_row(id, info)
    

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
            self.update_row(id, info)
        

class iSDE(_DataBase):
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

        self.add_row(event_to_save)


    def __getitem__(self, eid):
        return iEvent(eid, self)
    

    def relabel_event(self, id, new_label):
        if self.is_id(id):
            info = dict(label=new_label)
            self.update_row(id, info)
    

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
