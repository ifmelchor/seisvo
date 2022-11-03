#!/usr/bin/python3
# coding=utf-8

'''

Read write and operate .lte files

'''


import os
import h5py
import numpy as np

from obspy.core.util.attribdict import AttribDict


class cc8stats(AttribDict):
    def __init__(self, header):
        super(cc8stats, self).__init__(header)

    def __add_attr__(self, list):
        self.attributes = list

    def __str__(self):
        priorized_keys = ['id','starttime','endtime','window_length','overlap','nro_time_bins','last_time_bin','sample_rate','fq_band','np','pmax','pinc', 'ccerr','sup_file']

        return self._pretty_str(priorized_keys)

# class CC8s(object):
    # class for supfile


class CC8(object):
    def __init__(self, filename, sup_path=None):
        if not os.path.isfile(filename):
            return ValueError(' file not found')
        
        self.file_ = filename
        
        try:
            self.__set_stats__()
        except:
            return ValueError(' file %s cannot be read' %filename)
        
        # search for supfiles
        if self.stats.sup_file:
            self.supfiles_ = []
            for n in range(self.stats.np):
                supfile_n = filename + f'.{n}'
                
                if sup_path and os.path.isdir(sup_path):
                    supfile_n = os.path.join(sup_path, supfile_n)
                
                if os.path.isfile(supfile_n):
                    # read supfile
                    sfilen = CC8s(supfile_n)
                    self.supfiles_.append(sfilen)


    def __set_stats__(self):
        with h5py.File(self.file_, "r") as f:
            hdr = f['header']
            self.stats = cc8stats(
                dict(
                id = hdr.attrs['id'],
                starttime = dt.datetime.strptime(hdr.attrs['starttime'], '%Y-%m-%d %H:%M:%S'),
                endtime = dt.datetime.strptime(hdr.attrs['endtime'], '%Y-%m-%d %H:%M:%S'),
                window_length = hdr.attrs['window_length'],
                overlap = hdr.attrs['overlap'],
                nro_time_bins = hdr.attrs['nro_time_bins'],
                last_time_bin = hdr.attrs['last_time_bin'],
                sample_rate = hdr.attrs['sample_rate'],
                fq_band = hdr.attrs['fq_band'],
                np = hdr.attrs['np'],
                pmax = hdr.attrs['pmax'],
                pinc = hdr.attrs['pinc'],
                ccerr = hdr.attrs['ccerr'],
                sup_file = hdr.attrs['sup_file']
                )
            )
    
    def __compute__(self, sarray, headers):

        if self.stats.last_time_bin > 0:


        # loop over the intervals
        interval = headers["interval"]
        start_time = self.stats.starttime
        end_time = self.stats.endtime
        nini = int(headers['nini'])
        nwin = int(headers['nwin'])
        lwin = int(headers['lwin'])
        fsem = int(self.stats.sample_rate)

        start = start_time
        while start + interval <= end_time:
            end = start + interval
            mdata, stats = sarray.get_mdata(start, end, toff=nini*self.stats.sample_rate, sample_rate=self.stats.sample_rate, fq_band=self.stats.fq_band, return_stats=True)

            # prepare parameters
            xsta_utm = [st.lon for st in stats]
            ysta_utm = [st.lat for st in stats]
            
            # run cc8
            cc8run_jl(mdata, xsta_utm, ysta_utm, self.stats.pmax, self.stats.pinc, fsem, lwin, nwin, headers['nadv'], self.stats.ccerr, nini, self.file_, headers["supfile"])

            start += interval

    

    @staticmethod
    def new(sarray, filename, headers, sup_path=None):
        
        # print info
        print('')
        print(' CC8 file INFO')
        print(' -------------')
        print('  file   ::  %s ' % filename)
        for key in ['id', 'starttime', 'endtime' ,'window_length' ,'overlap' ,'nro_time_bins' ,'sample_rate' ,'fq_band' , 'np', 'pmax', 'pinc', 'ccer', 'sup_file']:
            print(f'  {key}  ::  {headers[key]}')
        print('')

        # create file
        cc8main = h5py.File(filename, "w-")
        
        # add header
        hdr = cc8main.create_dataset('header',(1,))
        hdr.attrs['id'] = headers['id']
        hdr.attrs['starttime'] = headers['starttime']
        hdr.attrs['endtime'] = headers['endtime']
        hdr.attrs['window_length'] = headers['window_length']
        hdr.attrs['overlap'] = headers['overlap']
        hdr.attrs['nro_time_bins'] = headers['nro_time_bins']
        hdr.attrs['last_time_bin'] = -1
        hdr.attrs['sample_rate'] = headers['sample_rate']
        hdr.attrs['fq_band'] = headers['fq_band']
        hdr.attrs['np'] = headers['np']
        hdr.attrs['pmax'] = headers['pmax']
        hdr.attrs['pinc'] = headers['pinc']
        hdr.attrs['ccerr'] = headers['ccerr']
        hdr.attrs['sup_file'] = headers['sup_file']

        # add datasets
        timebins = headers['nro_time_bins']

        for n in range(headers['np']):
            np_n = cc8main.create_group(str(n))
            np_n.create_dataset('slow', (timebins,), chunks='auto', dtype=np.float32)
            np_n.create_dataset('baz', (timebins,), chunks='auto', dtype=np.float32)
            np_n.create_dataset('cc_max', (timebins,), chunks='auto', dtype=np.float32)
            
            bds_n = np_n.create_group('bounds')
            bds_n.create_dataset('baz_min', (timebins,), chunks='auto', dtype=np.float32)
            bds_n.create_dataset('baz_max', (timebins,), chunks='auto', dtype=np.float32)
            bds_n.create_dataset('slo_min', (timebins,), chunks='auto', dtype=np.float32)
            bds_n.create_dataset('slo_max', (timebins,), chunks='auto', dtype=np.float32)
        
        cc8main.flush()
        cc8main.close()

        if headers['sup_file']:
            if sup_path:
                if not os.path.isdir(sup_path):
                    os.makedirs(sup_path)
                    print(f" >>> dir {sup_path} created")
                headers["supfile"] = os.path.join(sup_path, filename)
            else:
                headers["supfile"] = filename

            for n in range(headers['np']):
                # create supfiles
                supfile_n = filename + f'.{n}'
                if sup_path:
                    supfile_n = os.path.join(sup_path, supfile_n)
                cc8s_n = h5py.File(supfile_n, "w-")
                
                nite = 2*int(header["pmax"][n]/header["pinc"][n]) + 1
                dset = cc8main.create_dataset('sumap',(timebins,nite,nite), chunks='auto', dtype=np.float32)
                dset.attrs['last_time_bin'] = -1

                cc8s_n.flush()
                cc8s_n.close()

        else:
            headers["supfile"] = None

        cc8 = CC8(filename, sup_path=sup_path)
        cc8.__compute__(sarray, headers)

        return cc8
