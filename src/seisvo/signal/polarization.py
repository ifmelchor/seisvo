#!/usr/bin/python3
# coding=utf-8

import numpy as np
import cmath as cm
from itertools import product
from seisvo.signal import freq_bins
from seisvo.signal.spectrum import cosine_taper, mtspec


class PolarAnalysis(object):
    def __init__(self, Zdata, Ndata, Edata, sample_rate, npts_mov_avg=None, olap=0.0, fq_band=(0.5, 10),full_analysis=False, **kwargs):

        self.data = [Zdata, Ndata, Edata]
        npts = list(set(list(map(len, self.data))))

        if len(npts) > 1:
            raise ValueError(' data with different lengths')
        
        self.npts = npts[0]
        nfft = kwargs.get("nfft", 'uppest')
        
        if npts_mov_avg:
            if not isinstance(olap, float) or not 0 <= olap < 1:
                raise ValueError(' olap must be a float between 0 and 1')
            
            if not isinstance(npts_mov_avg, int) or not npts_mov_avg < self.npts:
                raise ValueError(' npts_mov_avg must be an integer lower than npts')
        
            self.freq, _, self.nfft = freq_bins(npts_mov_avg, sample_rate, fq_band=fq_band, nfft=nfft, get_freq=True)
        
        else:
            self.freq, self.fq_pos, self.nfft = freq_bins(self.npts, sample_rate, fq_band=fq_band, nfft=nfft, get_freq=True)

        self.sample_rate = sample_rate
        self.npts_mov_avg = npts_mov_avg
        self.olap_step = olap
        self.fq_band = fq_band
        self.taper = kwargs.get('taper', False)
        self.taper_p = kwargs.get('taper_p', 0.05)
        self.time_bandwidth = kwargs.get('time_bandwidth', 3.5)
        self.full_analysis = full_analysis

        self.fit()

    def fit(self):
        msize = (len(self.freq),)
        self.polar_dgr = np.zeros(msize, dtype='float64')
        
        if self.full_analysis:
            self.azimuth = np.zeros(msize, dtype='float64')
            self.elevation = np.zeros(msize, dtype='float64')
            self.rect = np.zeros(msize, dtype='float64')

        if self.npts_mov_avg:
            npts_olap = int(self.npts_mov_avg*self.olap_step)
        else:
            # one step
            self.npts_mov_avg = self.npts
            npts_olap = 0

        npts_start = 0
        step = 0
        while npts_start + self.npts_mov_avg <= self.npts:
            Zdata_n = self.data[0][npts_start:npts_start + self.npts_mov_avg]
            Ndata_n = self.data[1][npts_start:npts_start + self.npts_mov_avg]
            Edata_n = self.data[2][npts_start:npts_start + self.npts_mov_avg]
            ans = self.get_polar_attributes(Zdata_n, Ndata_n, Edata_n)

            self.polar_dgr += ans[0]

            if self.full_analysis:
                self.rect += ans[1]
                self.azimuth += ans[2]
                self.elevation += ans[3]

            npts_start += self.npts_mov_avg - npts_olap
            step += 1
        
        self.polar_dgr /= step
        self.rect /= step
        self.azimuth /= step
        self.elevation /= step


    def get_polar_attributes(self, Zdata, Ndata, Edata):
        data = [Zdata, Ndata, Edata]
        cmatrix = np.zeros((3, 3, len(self.freq)), dtype='complex128')

        index = set(product(set(range(3)), repeat=2))
        for i, j in list(index):
            cmatrix[i,j] = self.__cross_spec__(data[i], data[j])
        
        polar_dgr = np.empty((len(self.freq),), dtype=np.float32)

        z_list = []
        for fq in range(len(self.freq)):
            _, s, vh = np.linalg.svd(cmatrix[:,:,fq], full_matrices=False)
            polar_dgr[fq] = (3 * np.sum(s**2) - np.sum(s) ** 2) / (2 * np.sum(s) ** 2)

            if self.full_analysis:
                z = vh[0,:]
                z_list += [rotate(z)]
        
        if self.full_analysis:
            rect = np.array(list(map(get_rectiliniarity, z_list)))
            azimuth = np.array(list(map(get_azimuth, z_list)))
            elevation = np.array(list(map(get_elevation, z_list)))
            return polar_dgr, rect, azimuth, elevation
        
        else:
            return polar_dgr 


    def __cross_spec__(self, xdata, ydata):
        delta = float(1/self.sample_rate)
        
        xdata_mean = np.nanmean(xdata[np.isfinite(xdata)])
        xdata = xdata - xdata_mean

        ydata_mean = np.nanmean(ydata[np.isfinite(ydata)])
        ydata = ydata - ydata_mean

        if self.taper:
            tap = cosine_taper(self.npts, p=self.taper_p)
            xdata = xdata * tap
            ydata = ydata * tap

        x = mtspec(data=xdata, delta=delta, nfft=self.nfft, time_bandwidth=self.time_bandwidth, optional_output=True)
        y = mtspec(data=ydata, delta=delta, nfft=self.nfft, time_bandwidth=self.time_bandwidth, optional_output=True)

        psd_xy = np.zeros((len(self.freq),), dtype='complex128')
        
        nro_tapers = int(2* self.time_bandwidth) - 1
        for k in range(nro_tapers):
            psd_xy += x[3][self.fq_pos[0]:self.fq_pos[1], k] * np.conj(y[3][self.fq_pos[0]:self.fq_pos[1], k])
        
        return psd_xy


def rotate(z):
    def rot_dgr(alpha):
        z_rot = z * (np.cos(alpha)+1j*np.sin(alpha))
        a = np.array([x.real for x in z_rot])
        b = np.array([x.imag for x in z_rot])
        return np.abs(np.vdot(a,b))
    ans = np.array(list(map(rot_dgr, np.arange(0, np.pi/2, 1e-2))))
    ans_min = int(np.argmin(ans))
    phi_min = np.arange(0, np.pi/2, 1e-2)[ans_min]
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
    abs_zN, phyN = cm.polar(z[1])
    abs_zE, phyE = cm.polar(z[2])
    _ , phyH = cm.polar(z[2]**2+z[1]**2)
    th_l = -0.5*phyH + np.arange(1,10)*np.pi/2
    func = abs_zN**2*np.cos(th_l+phyN)**2 + abs_zE**2*np.cos(th_l+phyE)**2
    th_H = th_l[np.argmin(func)]
    zN_rot = z[1] * np.exp(-1j*th_H)
    zE_rot = z[2] * np.exp(-1j*th_H)
    tH = np.arctan(np.real(zE_rot)/np.real(zN_rot))*180/np.pi
    arg = np.real(z[0]*np.conjugate(z[2]))
    if arg < 0:
        tH += 90
    else:
        tH -= 90
    if tH < 0:
        tH += 180
    return tH


def get_elevation(z):
    abs_zV, phyV = cm.polar(z[0])
    zH = np.sqrt(z[1]**2 + z[2]**2)
    abs_zH, phyH = cm.polar(zH)
    _ , phyVH = cm.polar(z[0]**2 + zH**2)
    th_m = -0.5*phyVH + np.arange(1,10)*np.pi/2
    func = abs_zV**2*np.cos(th_m + phyV)**2 + abs_zH**2*np.cos(th_m + phyH)**2
    th_V = th_m[np.argmax(func)]
    zV_rot = np.real(z[0] * np.exp(-1j*th_V))
    zH_rot = np.real(zH * np.exp(-1j*th_V))
    tV = np.arctan(np.abs(zV_rot/zH_rot))
    return tV*180/np.pi




