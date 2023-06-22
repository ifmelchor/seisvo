
import datetime as dt
import sqlalchemy as sql
from ..network import Network

def get_sderow(SQLbase):

    class SDErow(SQLbase):
        __tablename__ = 'SDE'
        __table_args__ = {'extend_existing': True}

        # base attributes
        id = sql.Column(sql.Integer, primary_key=True)
        
        event_id   = sql.Column(sql.Integer, nullable=False)
        event_type = sql.Column(sql.String, nullable=False)
        
        network  = sql.Column(sql.String, nullable=False)
        station  = sql.Column(sql.String, nullable=False)
        location = sql.Column(sql.String, nullable=True)
        
        label     = sql.Column(sql.String, nullable=False)
        sublabel  = sql.Column(sql.String, nullable=True)
        starttime = sql.Column(sql.DateTime(timezone=False), nullable=False)
        duration  = sql.Column(sql.Float, nullable=False)

        # event attributes
        time_P   = sql.Column(sql.DateTime(timezone=False), nullable=True)
        weight_P = sql.Column(sql.Integer, nullable=True)
        onset_P  = sql.Column(sql.String, nullable=True)
        time_S   = sql.Column(sql.DateTime(timezone=False), nullable=True)
        weight_S = sql.Column(sql.Integer, nullable=True)
        time_F   = sql.Column(sql.DateTime(timezone=False), nullable=True)    

        # additional floats
        value_1 = sql.Column(sql.Float, nullable=True) # peak to peak
        value_2 = sql.Column(sql.Float, nullable=True) # max amplitude
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
            text_info = " EID %i || ID %i | %s | %s.%s.%s " % (self.event_id,\
                self.id, self.event_type, self.network, self.station, self.location)
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

    return SDErow


def get_lderow(SQLbase):

    class LDErow(SQLbase):
        __tablename__ = 'LDE'
        __table_args__ = {'extend_existing': True}

        id = sql.Column(sql.Integer, primary_key=True)
        network = sql.Column(sql.String, nullable=False)
        station = sql.Column(sql.String, nullable=False)
        location = sql.Column(sql.String, nullable=True)
        label     = sql.Column(sql.String, nullable=False)
        starttime = sql.Column(sql.DateTime(timezone=False), nullable=False)
        duration  = sql.Column(sql.Float, nullable=False)
        lte_file  = sql.Column(sql.String, nullable=False) # lte file containing the event info
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

        @property
        def endtime(self):
            return self.starttime + dt.timedelta(minutes=self.duration)
        
        @property
        def station_id(self):
            return '.'.join([self.network, self.station, self.location])

        def get_network(self):
            return Network(self.network)
        
        def get_station(self):
            net = self.get_network()
            return net.get_sta(self.station, self.location)

    return LDErow