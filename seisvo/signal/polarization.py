#!/usr/bin/python3
# coding=utf-8

import functools
from lib2to3.pgen2.token import N_TOKENS
import numpy as np
import cmath as cm
import multiprocessing
from seisvo.signal import freq_bins
from seisvo.signal.spectrum import cosine_taper, MTSpec
from mtspec import mtspec
import time as ttime

class PolarAnalysis(object):
    def __init__(self, Zdata, Ndata, Edata, sample_rate, npts_mov_avg=None, olap=0.0, fq_band=(0.5, 10),full_analysis=False, **kwargs):

        self.data = [Zdata, Ndata, Edata]
        npts = list(set(list(map(len, self.data))))

        if len(npts) > 1:
            raise ValueError(' data with different lengths')
        
        self.npts = npts[0]
        nfft = kwargs.get("nfft", 'uppest')
        
        if npts_mov_avg:
            npts_mov_avg = int(npts_mov_avg)

            if not isinstance(olap, float) or not 0 <= olap < 1:
                raise ValueError(' olap must be a float between 0 and 1')
            
            if npts_mov_avg > self.npts:
                raise ValueError(' npts_mov_avg must be an integer lower than npts')
        
            self.freq, self.fq_pos, self.nfft = freq_bins(npts_mov_avg, sample_rate, fq_band=fq_band, nfft=nfft, get_freq=True)
        
        else:
            self.freq, self.fq_pos, self.nfft = freq_bins(self.npts, sample_rate, fq_band=fq_band, nfft=nfft, get_freq=True)

        self.sample_rate = sample_rate
        self.npts_mov_avg = npts_mov_avg
        self.olap_step = olap
        self.fq_band = fq_band
        self.taper = kwargs.get('taper', False)
        self.taper_p = kwargs.get('taper_p', 0.05)
        self.time_bandwidth = kwargs.get('time_bandwidth', 3.5)
        self.njobs = kwargs.get('njobs', multiprocessing.cpu_count())
        self.full_analysis = full_analysis

        self.fit()


    def fit(self):
        # msize = (len(self.freq),)
        # self.degree = np.zeros(msize, dtype='float64')
        
        # if self.full_analysis:
        #     self.azimuth = np.zeros(msize, dtype='float64')
        #     self.elevation = np.zeros(msize, dtype='float64')
        #     self.rect = np.zeros(msize, dtype='float64')

        if self.npts_mov_avg:
            npts_olap = int(self.npts_mov_avg*self.olap_step)
        else:
            # one step
            self.npts_mov_avg = self.npts
            npts_olap = 0

        npts_start = 0

        data_split = []
        while npts_start + self.npts_mov_avg <= self.npts:
            Zdata_n = self.data[0][npts_start:npts_start + self.npts_mov_avg]
            Ndata_n = self.data[1][npts_start:npts_start + self.npts_mov_avg]
            Edata_n = self.data[2][npts_start:npts_start + self.npts_mov_avg]
            data_split += [[Zdata_n, Ndata_n, Edata_n]]
            npts_start += self.npts_mov_avg - npts_olap
        self.n = len(data_split)

        if self.n > 1:
            self.degree = np.empty(shape=(len(self.freq), self.n))
            if self.full_analysis:
                self.rect = np.empty(shape=(len(self.freq), self.n))
                self.azimuth = np.empty(shape=(len(self.freq), self.n))
                self.elevation = np.empty(shape=(len(self.freq), self.n))

            ans = list(map(self.__process__, data_split))
            for n in range(self.n):
                self.degree[:,n] = ans[n][0]

                if self.full_analysis:
                    self.rect[:,n] = ans[n][1]
                    self.azimuth[:,n] = ans[n][2]
                    self.elevation[:,n] = ans[n][3]

        else:
            # for data in data_split:
            ans = self.__process__(data_split)

            # for n_ans in ans:
            self.degree = np.array(ans[0])
            if self.full_analysis:
                self.rect = np.array(ans[1])
                self.azimuth = np.array(ans[2])
                self.elevation = np.array(ans[3])


    def __process__(self, data):
        """
        This code computes the averaged cross-spectral hermitian matrix and compute polarization attributes.

        data is a list of lists, on which each list contains Z, N, E data.
        """

        if len(data) > 1:
            data = [data]

        # compute cross-spectral elements
        cross_spec_func = functools.partial(self.__cross_spec__, data)
        index_list = [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
        with multiprocessing.Pool(self.njobs) as p:
            cross_spec_ans = list(p.map(cross_spec_func, index_list))

        # build hermitian matrix
        cmatrix = np.zeros((3, 3, len(self.freq)), dtype='complex128')
        for csa, (i, j) in zip(cross_spec_ans, index_list):
            cmatrix[i,j,:] = csa
            if i != j:
                cmatrix[j,i,:] = np.conj(csa)

        get_polar_degree = functools.partial(self.__polar_dgr__, cmatrix)
        index_list = range(len(self.freq))

        with multiprocessing.Pool(self.njobs) as p:
            polar_degree_ans = list(p.map(get_polar_degree, list(index_list)))

        if self.full_analysis:
            polar_dgr = [ans[0] for ans in polar_degree_ans]
            z_list = [ans[1] for ans in polar_degree_ans]
        else:
            polar_dgr = polar_degree_ans
        
        if self.full_analysis:
            rect = np.array(list(map(get_rectiliniarity, z_list)))
            azimuth = np.array(list(map(get_azimuth, z_list)))
            elevation = np.array(list(map(get_elevation, z_list)))
            return (polar_dgr, rect, azimuth, elevation)
        
        else:
            return polar_dgr


    def __polar_dgr__(self, cmatrix, fq_idx):
        _, s, vh = np.linalg.svd(cmatrix[:,:,fq_idx], full_matrices=False)
        polar_dgr = (3 * np.sum(s**2) - np.sum(s) ** 2) / (2 * np.sum(s) ** 2)
        
        if self.full_analysis:
            z = rotate(vh[0,:])
            return (polar_dgr, z)
        
        else:
            return polar_dgr


    def __cross_spec_xy__(self, data_t, n, k):
            xdata = data_t[n]
            ydata = data_t[k]

            if self.taper:
                tap = cosine_taper(self.npts, p=self.taper_p)
                xdata = xdata * tap
                ydata = ydata * tap

            # using mtspec is more rapid that MTSpec because it is a fortran wrapper

            # MTSx = MTSpec(xdata, dt=float(1/self.sample_rate), nfft=self.nfft, nw=self.time_bandwidth)
            # xTrue = MTSx.yk
            x = mtspec(data=xdata, delta=float(1/self.sample_rate), nfft=self.nfft, time_bandwidth=self.time_bandwidth, optional_output=True)[3]
            # MTSy = MTSpec(ydata, dt=float(1/self.sample_rate), nfft=self.nfft, nw=self.time_bandwidth)
            # y = MTSy.yk
            y = mtspec(data=ydata, delta=float(1/self.sample_rate), nfft=self.nfft, time_bandwidth=self.time_bandwidth, optional_output=True)[3]
        
            nro_tapers = int(2* self.time_bandwidth) - 1
            psd_xy_t = np.zeros((len(self.freq),), dtype='complex128')
            
            for nt in range(nro_tapers):
                psd_xy_t += x[self.fq_pos[0]:self.fq_pos[1], nt] * np.conj(y[self.fq_pos[0]:self.fq_pos[1], nt])
        
            return psd_xy_t
        

    def __cross_spec__(self, data, idx):
        # data is a list os lists
        
        fxy = functools.partial(self.__cross_spec_xy__, n=idx[0], k=idx[1])
        csans = list(map(fxy, data))
        
        if self.n > 1:
            psd_xy = np.zeros((len(self.freq),), dtype='complex128')
            for ans in csans:
                psd_xy += ans
            
            return psd_xy/self.n # mov average cross spectral element
        
        else:
            return csans[0] # cross spectral element


def rotate(z):
    def rot_dgr(alpha):
        z_rot = z * (np.cos(alpha)+1j*np.sin(alpha))
        a = np.array([x.real for x in z_rot])
        b = np.array([x.imag for x in z_rot])
        return np.abs(np.vdot(a,b))

    all_alpha = np.linspace(0, np.pi*2, 100)
    ans = np.array(list(map(rot_dgr, all_alpha)))
    phi_min = all_alpha[np.argmin(ans)]

    z_rot = z * (np.cos(phi_min)+1j*np.sin(phi_min))
    return z_rot


def get_rectiliniarity(z):
    a = np.array([x.real for x in z])
    b = np.array([x.imag for x in z])
    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)

    if a_norm > b_norm:
        minor, major = b_norm, a_norm
    
    else:
        minor, major = a_norm, b_norm

    rect = 1 - (minor/major)

    return rect


def get_azimuth(z):
    
    # decompose horizontal comp.
    abs_zN, phyN = cm.polar(z[1])
    abs_zE, phyE = cm.polar(z[2])
    
    _ , phyH = cm.polar(z[2]*z[2] + z[1]*z[1])
    th_l = -0.5*phyH + np.arange(0,5) * np.pi / 2
    
    func = abs_zN*abs_zN * np.cos(th_l + phyN)*np.cos(th_l + phyN) + abs_zE*abs_zE * np.cos(th_l+phyE)*np.cos(th_l+phyE)

    th_H = th_l[np.argmax(func)]
    zN_rot = z[1] * np.exp(-1j*th_H)
    zE_rot = z[2] * np.exp(-1j*th_H)
    tH = np.arctan(zE_rot.real/zN_rot.real)
    
    arg = (z[0]*z[2].conjugate()).real
    if arg < 0:
        if tH < 0:
            tH += np.pi
    else:
        if tH > 0:
            tH -= np.pi
    
    # measure clockwise from north
    tH = (np.pi/2) - tH
    if tH < 0:
      tH += 2 * np.pi

    # we only seek for the direction (no orientation)
    if tH > np.pi:
        tH -= np.pi

    # ----------------------------------------------------------
    
    ## get phiHH, which is the phase difference between horizontals
    # phiHH = (phyN - phyE) * 180 /np.pi
    
    # if(phiHH > 180.0):
    #     phiHH -= 360.0
    
    # elif (phiHH < -180.0):
    #     phiHH += 360.0

    ## get phiVH, between -90 and 90
    # phiVH = (th_H - cm.polar(z[0])[1]) * 180 / np.pi
    
    # if phiVH > 90:
    #     phiVH -= 180.0
    
    # if phiVH < -90:
    #     phiVH += 180.0
    
    return tH * 180 / np.pi


def get_elevation(z):

    # decompose first component
    abs_zV, phyV = cm.polar(z[0])
    
    # compute horizontal comp
    zH = np.sqrt(z[1]*z[1] + z[2]*z[2])
    abs_zH, phyH = cm.polar(zH)
    _ , phyVH = cm.polar(z[0]*z[0] + zH*zH)

    th_m = -0.5*phyVH + np.arange(0,5)*np.pi/2
    func = abs_zV*abs_zV*np.cos(th_m + phyV)*np.cos(th_m + phyV) + abs_zH*abs_zH*np.cos(th_m + phyH)*np.cos(th_m + phyH)
    th_V = th_m[np.argmax(func)]

    if(zH.imag < 0):
      zH = complex(-1*zH.real,-1*zH.imag)
    
    ztmp = complex(np.cos(th_V), -1*np.sin(th_V))
    tV = np.arctan(np.abs( (z[0]*ztmp).real / (zH*ztmp).real))
    tV  = np.pi/2 - tV

    # this gives the elevation angle, 0 for vertical, 90 for horizontal
    return tV * 180/np.pi


