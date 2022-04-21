#!/usr/bin/env python3
# coding=utf-8

import os
import datetime as dt
from seisvo import __seisvo__
from seisvo.core.network import Network
from seisvo.file.lte import LTE, Peaks

class Event(object):
    def __init__(self, event_id, sde):
        self.id = event_id
        self.sde = sde
        self.label = None
        self.starttime = None
        self.endtime = None
        self.duration = None
        self.stations = []
        self.event_duration = None
        self._setattr()


    def __len__(self):
        return len(self.stations)


    def __str__(self):
        text_info = " Event ID: %s\n" % self.id
        text_info += "    Database       : %s\n" % self.sde.sql_path
        text_info += "    Label          : %s\n" % self.label
        text_info += "    Starttime      : %s\n" % self.starttime.strftime('%Y-%m-%d %H:%M')
        text_info += "    Duration [min] : %.2f\n" % self.duration
        text_info += "    Stations       : %s\n" % (self.stations)
        return text_info


    def _setattr(self):
        self.events = self.sde.get_event(self.id)

        max_duration = []
        for e in self.events:
            if e.label:
                self.label = e.label
                self.starttime = e.starttime
                self.duration = e.duration
                self.endtime = e.starttime + dt.timedelta(seconds=e.duration)

            if e.event_duration:
                max_duration.append(e.event_duration)
            
            station_id = '.'.join([e.network, e.station, e.location])
            self.stations.append(station_id)
        
        if max_duration:
            self.event_duration = max(max_duration)
            

    def __getitem__(self, i):
        return self.events[i]


    def get_row(self, station_id):
        network = station_id.split('.')[0]
        station = station_id.split('.')[1]
        location = station_id.split('.')[2]
        n_row = self.sde.get_row(self.id, network, station, location)
        if n_row:
            return self.sde.get_id(n_row)
        else:
            return None


    def plot(self, **kwargs):
        from seisvo.utils.plotting import plot_sde_event
        
        ans = plot_sde_event(self, **kwargs)

        return ans


    def to_hypo71(self, out_dir=None):
        """
        Write a hypo71 phase file based on event info.

        Parameters
        ----------
        out_dir : [str]
            Output directory. If output is None, return a string

        """

        stime = self.starttime

        if not out_dir:
            return_string = []
        else:
            return_string = None
            sql_name = os.path.basename(self.sde.sql_path)

            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)
                print( 'Directory %s created' %out_dir)
            
            file_name = os.path.join(out_dir, sql_name + '_' + str(self.id) + '.obs')

        for staid in self.stations:
            row = self.get_row(staid)
            line = f'{row.station:<4}'
            if row.onset_P:
                line += f'{row.onset_P.lower()}P'
            else:
                line += ' P'
            if row.first_motion_P:
                line += '{}'.format(row.first_motion_P)
            else:
                line += ' '
            
            stime + dt.timedelta(seconds=row.time_P)
            timestr = stime.strftime('%y%m%d%H%M')

            if isinstance(row.weight_P, int):
                line += '{} {}'.format(row.weight_P, timestr)
            else:
                line += '  {}'.format(timestr)

            if row.time_P:
                line += f'{float(row.time_P):>5.2f}'
            else:
                line += '     '

            line += '       '

            if row.time_S:
                line += f'{float(row.time_S):>5.2f} S '

                if isinstance(row.weight_S, int):
                    line += '{}'.format(row.weight_S)
                else:
                    line += ' '
            else:
                line += ' '*5 + ' '*3 + ' '
            
            line += '   '

            if row.peak_to_peak:
                line += f'{float(row.peak_to_peak):>1.2f}'
            else:
                line += '    '
            
            if row.max_period:
                line += f'{float(row.max_period):>3.1f}'
            else:
                line += '   '
            
            line += ' '*20
            if row.event_duration:
                line += f'{float(row.event_duration):>5.0f}'
            
            if isinstance(return_string, list):
                return_string.append(line)
            
            else:
                if not os.path.isfile(file_name):
                    print('  >> %s  added' % file_name)
                with open(file_name, 'a') as fle:
                    fle.write(line+'\n')
                
            
        return return_string


    def get_stream(self, remove_response=True, **kwargs):
        stream = None
        
        for sta in self.stations:
            net_code = sta.split('.')[0]
            sta_code = sta.split('.')[1]
            loc = sta.split('.')[2]
            net = Network(net_code)
            sta = net.get_sta(sta_code, loc=loc)
            st1 = sta.get_stream(self.starttime, self.starttime+dt.timedelta(seconds=self.duration), remove_response=remove_response, **kwargs)
            if not stream:
                stream = st1
            else:
                stream += st1

        return stream


class Episode(object):
    def __init__(self, event_id, lde):
        self.id = event_id
        self.lde = lde
        self.main_lte = None
        self.sup_lte = None
        self.label = None
        self.starttime = None
        self.endtime = None
        self.duration = None
        self.station = []
        self._setattr()
    

    def _setattr(self):
        self.row_ = self.lde.get_id(self.id)
        self.station_ = self.row_.get_station()

        self.label = self.row_.label
        self.starttime = self.row_.starttime
        self.duration = self.row_.duration
        self.endtime = self.row_.starttime + dt.timedelta(hours=self.row_.duration)

        if os.path.isfile(self.row_.lte_file):
            self.main_lte = LTE(self.row_.lte_file)
        else:
            print('warn: main lte file not found!')
        
        if os.path.isfile(self.row_.lte_file_sup):
            self.sup_lte = LTE(self.row_.lte_file_sup)
    

    def getPeaks(self, fq_range=(), peak_thresholds={}):
        if self.sup_lte:
            return Peaks(self.sup_lte, fq_range=fq_range, peak_thresholds=peak_thresholds)


    def compute_sup_lte(self, dir_out=None, time_step=1, sample_rate=None, avg_step=None, **ltekwargs):

        if not dir_out:
            dir_out = os.path.join(__seisvo__, 'database', self.row_.network, 'sup')

        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)
            
        file_name = '%s_%i.lte' % (self.lde.id, self.id)
        file_name_path = os.path.join(dir_out, file_name)

        if os.path.isfile(file_name_path):
            os.remove(file_name_path)
            self.lde.__update_lte_sup__(self.id, None)

        ltekwargs['file_name'] = file_name
        ltekwargs['out_dir'] = dir_out

        self.station_.lte(
            self.starttime, 
            self.endtime, 
            time_step,
            polargram=True, 
            sample_rate=sample_rate, 
            avg_step=avg_step,
            polarization_attr=True, 
            **ltekwargs)

        self.lde.__update_lte_sup__(self.id, file_name_path)
