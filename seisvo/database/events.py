#!/usr/bin/env python3
# coding=utf-8

import os
import numpy as np
import datetime as dt
from scipy.stats import gaussian_kde
from seisvo import seisvo_paths

# import seisvo.file.lte as sfl
# from seisvo.file.air import AiR


class Event(object):
    def __init__(self, event_id, sde):
        self.id = event_id
        self.sde = sde
        self.rows_ = self.sde.get_event(self.id)
        self.label     = self.rows_[0].label
        self.starttime = self.rows_[0].starttime
        self.duration  = self.rows_[0].duration
        self.endtime   = self.rows_[0].get_endtime()
        self.stations = [row.get_station_id() for row in self.rows_]
        

    def __len__(self):
        return len(self.stations)


    def __str__(self):
        text_info = " Event ID: %s\n" % self.id
        text_info += "    Database       : %s\n" % self.sde.sql_path
        text_info += "    Label          : %s\n" % self.label
        text_info += "    Starttime      : %s\n" % self.starttime.strftime('%Y-%m-%d %H:%M')
        text_info += "    Duration [sec] : %.2f\n" % self.duration
        text_info += "    Stations       : %s\n" % (self.stations)
        return text_info
     

    def __getitem__(self, i):
        return self.rows_[i]


    def get_stream(self, toff_sec=0, **kwargs):
        stream = None
        toff = dt.timedelta(seconds=toff_sec)
        for row in self.rows_:
            sta = row.get_station()
            st1 = sta.get_stream(self.starttime-toff, self.endtime+toff, **kwargs)
            if st1:
                if not stream:
                    stream = st1
                else:
                    stream += st1
        return stream


    def append_station(self, station_id, etype="S"):
        if station_id not in self.stations:
            dinfo = {
                'network'   : station_id.split(".")[0],
                'station'   : station_id.split(".")[1],
                'location'  : station_id.split(".")[2],
                'event_type': etype,
                'label'     : self.label,
                'starttime' : self.starttime,
                'duration'  : self.duration,
                'event_id'  : self.id,
            }
            self.sde.add_row(dinfo)
            self.stations.append(station_id)
    

    def remove_station(self, station_id):
        row = self.get_row(station_id)
        if row and len(self.stations)>1:
            self.sde.remove_row(row.id)
            self.stations.remove(station_id)


    def relabel(self, new_label):
        for row in self.rows_:
            self.sde.update_row(row.id, dict(label=new_label))
            self.label = new_label


    def get_row(self, station_id):
        for row in self.rows_:
            staid = '.'.join([row.network, row.station, row.location])
            if staid == station_id:
                return row


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


    def get_phases(self, nro_phases=False):
        phase = {}
        nphase = 0
        for sta in self.stations:
            row = self.get_row(sta)
            
            if row.time_P or row.time_S or row.time_F:
                phase[sta] = {}
                
                if row.time_P:
                    phase[sta]["P"] = row.time_P
                    nphase += 1
                
                if row.time_S:
                    phase[sta]["S"] = row.time_S
                    nphase += 1
                
                if row.time_F:
                    phase[sta]["F"] = row.time_F
                    nphase += 1

        if nro_phases:
            return nphase

        return phase



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
    
    
    def __str__(self):
        text_info = " Event ID: %s ::%s (%s) \n" % (self.lde.id, self.label, self.id)
        text_info += "    LTE_file      : %s\n" % self.row_.lte_file
        text_info += "    LTE_file_sup  : %s\n" % self.row_.lte_file_sup
        text_info += "    Starttime     : %s\n" % self.starttime.strftime('%Y-%m-%d %H:%M')
        text_info += "    Duration [hr] : %.2f\n" % self.duration
        return text_info


    def _setattr(self):
        self.row_ = self.lde.get_id(self.id)
        self.station_ = self.row_.get_station()

        self.label = self.row_.label
        self.starttime = self.row_.starttime
        self.duration = self.row_.duration
        self.endtime = self.row_.starttime + dt.timedelta(hours=self.row_.duration)

        if self.row_.lte_file and os.path.isfile(self.row_.lte_file):
            self.main_lte = sfl.LTE(self.row_.lte_file)
        else:
            print('warn [ID %i]: main lte file not found' %self.id)
        
        if self.row_.lte_file_sup and os.path.isfile(self.row_.lte_file_sup):
            self.sup_lte = sfl.LTE(self.row_.lte_file_sup)


    def compute_sup_lte(self, channel, interval, int_olap=0.0, step=1, step_olap=0.0, out_dir=None, replace=True, **ltekwargs):

        if not out_dir:
            out_dir = os.path.join(seisvo_paths["main"], 'database', self.row_.network, 'sup')

        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
            
        file_name = '%s_%i.lte' % (self.station_.stats.id, self.id)
        file_name_path = os.path.join(out_dir, file_name)

        if os.path.isfile(file_name_path):

            if replace:
                os.remove(file_name_path)
                self.lde.__update_lte_sup__(self.id, None)
            else:
                print('File exist. Do nothing.')
                return None

        ltekwargs['polar_analysis'] = True
        ltekwargs['file_name'] = file_name
        ltekwargs['out_dir'] = out_dir

        self.station_.lte(
            self.starttime,
            self.endtime,
            channel=channel,
            interval=interval,
            int_olap=int_olap,
            step=step,
            step_olap=step_olap,
            **ltekwargs)

        self.lde.__update_lte_sup__(self.id, file_name_path)


    def get_stream(self, **kwargs):
        return self.station_.get_stream(self.starttime, self.endtime, **kwargs)
    

    def get_values(self, which_lte='sup', attrs=[]):
        """Return min, max, mean, mode of SCALARs parameters

        Parameters
        ----------
        which_lte : str, optional
            'main' or 'sup', by default 'sup'
        attrs : list or str, optional
            the attributes related to SCALAR parameters, by default []

        Returns
        -------
        dict
        """

        if which_lte == 'sup':
            lte = self.sup_lte
        
        elif which_lte == 'main':
            lte = self.main_lte
        
        else:
            print(' which_lte should be "main" or "sup"')
            return
        
        if isinstance(attrs, (list, tuple, str)):

            if isinstance(attrs, str):
                if attrs not in lte.stats.attributes and sfl.SCALAR_PARAMS + sfl.SCALAR_OPT_PARAMS:
                    print(' attr not found or is not scalar')
                    return
                else:
                    attrs = [attrs]
            
            else:
                checked_attrs = []
                for atr in attrs:
                    if atr in lte.stats.attributes and sfl.SCALAR_PARAMS + sfl.SCALAR_OPT_PARAMS:
                        checked_attrs.append(atr)
                    else:
                        print('warn: %s not found or is not scalar' %atr)
                attrs = checked_attrs
        
        else:
            print(' attr should be string or list')
            return
        
        if not attrs:
            print(' attrs not found')
            return

        dict_out = {}

        for attr in attrs:
            data = lte.get_attr(attr)
            data = data[np.isfinite(data)]

            if attr == 'energy':
                data = 10*np.log10(data)
            
            v_min = data.min()
            v_max = data.max()     
            x_range = np.linspace(v_min, v_max, 500)
            gkde = gaussian_kde(data)
            kde = gkde(x_range)
            mode = x_range[np.argmax(kde)]

            dict_out[attr] = [v_min, v_max, data.mean(), mode]
        
        return dict_out


    def get_matrix(self, attr, which_lte='sup'):

        if which_lte == 'sup':
            lte = self.sup_lte
        
        elif which_lte == 'main':
            lte = self.main_lte
        
        else:
            print(' which_lte should be "main" or "sup"')
            return
        
        if attr not in lte.stats.attributes and sfl.VECTORAL_PARAMS + sfl.VECTORAL_OPT_PARAMS:
            print(' attr not found or is not vectorial')
        
        return lte.get_attr(attr)


    def get_peaks(self, fq_range=(), peak_thresholds={}):
        if self.sup_lte:
            return sfl.Peaks(self.sup_lte, fq_range=fq_range, peak_thresholds=peak_thresholds)
    

    def plot(self, fig=None, return_fig=False, **kwargs):
        from seisvo.gui.glde import plot_event
        return plot_event(self, fig=fig, return_fig=True, **kwargs)


    def plot_polar(self, **kwargs):
        from seisvo.gui.glde import plot_event_polar
        if self.lte_file_sup:
            plot_event_polar(self, **kwargs)



class iEvent(object):
    def __init__(self, event_id, isde):
        self.id = event_id
        self.isde = isde
        self.label = None
        self.starttime = None
        self.endtime = None
        self.duration = None
        self.air = None
        self.model_ = None
        self.channels = []
        self._setattr()
    

    def _setattr(self):
        self.row_ = self.lde.get_id(self.id)
        self.label = self.row_.label
        self.starttime = self.row_.starttime
        self.duration = self.row_.duration
        self.endtime = self.row_.starttime + dt.timedelta(seconds=self.row_.duration)

        if os.path.isfile(self.row_.air_file):
            self.air = AiR(self.row_.lte_file)
            self.model_ = dict(
                radii=self.air.stats.radii,
                vel_air=self.air.stats.vel_air,
                h_src=self.air.h_src,
                src_dgr=self.air.stats.src_dgr
                )
        else:
            print('warn: air file not found!')
    
        self.array_ = self.row_.get_array(model=self.model_)
    

    def get_stream(self, azm, model=None, time_pad=0, **kwargs):
        return self.array_.get_stream(self.starttime, self.endtime, azm, model=model, time_pad=time_pad, **kwargs)
    

    

