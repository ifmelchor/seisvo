#!/usr/bin/env python3
# coding=utf-8

import datetime as dt
from seisvo.network import Network
from .base import sql, DataBase, SQLbase
from .events import Event


class SDErow(SQLbase):
    __tablename__ = 'SDE'

    # base attributes
    id = sql.Column(sql.Integer, primary_key=True)
    event_id = sql.Column(sql.Integer, nullable=False)
    event_type = sql.Column(sql.String, nullable=False)
    
    network = sql.Column(sql.String, nullable=False)
    station = sql.Column(sql.String, nullable=False)
    location = sql.Column(sql.String, nullable=True)
    
    label = sql.Column(sql.String, nullable=False)
    sublabel = sql.Column(sql.String, nullable=True)
    
    starttime = sql.Column(sql.DateTime(timezone=False), nullable=False)
    duration = sql.Column(sql.Float, nullable=False)

    # event attributes
    time_P = sql.Column(sql.DateTime(timezone=False), nullable=True)
    weight_P = sql.Column(sql.Integer, nullable=True)
    onset_P = sql.Column(sql.String, nullable=True)

    time_S = sql.Column(sql.DateTime(timezone=False), nullable=True)
    weight_S = sql.Column(sql.Integer, nullable=True)

    time_F = sql.Column(sql.DateTime(timezone=False), nullable=True)
    
    peak_to_peak = sql.Column(sql.Float, nullable=True)

    # others
    # event_amplitude = sql.Column(sql.Float, nullable=True)
    # event_duration = sql.Column(sql.Float, nullable=True)
    # event_frequency = sql.Column(sql.Float, nullable=True)
    # azimuth = sql.Column(sql.Float, nullable=True)
    # incidence = sql.Column(sql.Float, nullable=True)

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
    

    def get_station_id(self):
        return '.'.join([self.network, self.station, self.location])


    def get_network(self):
        return Network(self.network)


    def get_station(self):
        return self.get_network().get_sta(self.station, self.location)
    

    def get_stream(self, off_seconds=0, **kwargs):
        sta = self.get_station()
        start = self.starttime - dt.timedelta(seconds=off_seconds)
        end   = self.starttime + dt.timedelta(seconds=self.duration+off_seconds)
        st = sta.get_stream(start, end, **kwargs)
        return st
    

    def get_endtime(self):
        return self.starttime + dt.timedelta(seconds=self.duration)


class SDE(DataBase):
    def __init__(self, sql_path):
        super().__init__(sql_path, 'SDE')
        self.__set_attr__()


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
            id = self.add_row(event_to_save)
        
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
                self.remove_row(row.id)
        
        if isinstance(eid, (list, tuple)):
            for eidi, row_list in self.get_event(eid).items():
                if self.is_eid(eidi):
                    for row in row_list:
                        self.remove_row(row.id)
        
        self.__set_attr__()
    

    def clone_event(self, eid):
        if self.is_eid(eid):
            row_list = get_event(eid)
            station_id_list = [row.get_station_id() for row in row_list]
            label = row_list[0].label
            starttime = row_list[0].starttime
            duration = row_list[0].duration
            etype = row_list[0].event_type
            self.add_event(station_id_list, label, starttime, duration, etype)


    def plot_gui(self, event_id=None, station_id=None, init_app=True, **kwargs):
        eid_list = self.get_id_list(**kwargs)

        if event_id:
            assert event_id in eid_list
        
        event_id = eid_list[0]
        widget = plot_event(self, event_id, station_id, init_app=init_app)

        return widget
