#!/usr/bin/python3
# coding=utf-8

from __future__ import division
import numpy as np
import cmath as cm
from numba import njit, jit, float32, complex128


def pplot(a, b, text=None):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    print(a, b)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X, Y, Z = ([0,0], [0,0], [0,0])
    U, V, W = ([a[0],b[0]], [a[1],b[1]], [a[2],b[2]])
    ax.quiver(X, Y, Z, U, V, W, color=['k', 'r'], length=0.2)
    ax.set_title(text)
    plt.show()


def polar_degree(idata, jdata, kdata, sample_rate, per_lap=0.7, avg_step=None, 
    fq_band=[], matrix_return=False, **kwargs):
    """
    Compute the polargram...
    :param avg_step: int, in minutes
    :param opt_return: boolean, if True return azimuth, dip and rectiliniarity
    """

    from seisvo.signal import freq_bins
    from seisvo.signal.spectrum import cross_spectrum_matrix

    npts = len(idata)
    
    if not avg_step:
        npts_step = npts
    
    else:
        npts_step = int(avg_step*60*sample_rate-1)

    nfft = kwargs.get("nfft", 'uppest')
    freq_n, (fnptlo, fnpthi), nfft_n = freq_bins(npts_step, sample_rate,
            fq_band=fq_band, nfft=nfft, get_freq=True)

    npts_olap = int(npts_step*per_lap)

    if npts_olap >= npts_step:
        raise ValueError('overlap must be shorter than step')

    polar_dgr = np.zeros((len(freq_n),), dtype='float64')
    
    if matrix_return:
        azm = np.zeros((len(freq_n),), dtype='float64')
        dip = np.zeros((len(freq_n),), dtype='float64')
        rect = np.zeros((len(freq_n),), dtype='float64')

    threshold = kwargs.get('threshold', 0)
    npts_start = 0
    k = 0
    while npts_start + npts_step <= npts:
        idata_n = idata[npts_start:npts_start + npts_step]
        jdata_n = jdata[npts_start:npts_start + npts_step]
        kdata_n = kdata[npts_start:npts_start + npts_step]
        spec_mat, _ = cross_spectrum_matrix(idata_n, jdata_n, kdata_n, sample_rate, fq_band=fq_band, **kwargs)

        # polarization analysis
        ans = polarization_content(spec_mat, matrix_return=matrix_return, threshold=threshold)

        polar_dgr += ans[0]
        
        if matrix_return:
            azm += ans[1]
            dip += ans[2]
            rect += ans[3]

        k += 1
        npts_start += npts_step - npts_olap

    if matrix_return:
        return freq_n, polar_dgr/k, azm/k, dip/k, rect/k
    
    else:
        return freq_n, polar_dgr/k


@jit(float32(complex128[:]))
def get_azimuth(z):
    """
    Code to give the azimuth and the phy2-phy3 acoording to Park et al (1987).
    """

    abs_z2, phy2 = cm.polar(z[1])
    abs_z3, phy3 = cm.polar(z[2])

    _ , phy23 = cm.polar(z[2]**2+z[1]**2)

    th_l = -0.5*phy23 + np.arange(1,10)*np.pi/2
    func = abs_z2**2*np.cos(th_l+phy2)**2 + abs_z3**2*np.cos(th_l+phy3)**2
    th_H = th_l[np.argmin(func)]

    z2_rot = z[1] * np.exp(-1j*th_H)
    z3_rot = z[2] * np.exp(-1j*th_H)

    # z2 shuld be NORTH and z3 shuld de EAST
    azm = np.real(np.arctan(z2_rot/z3_rot))*180/np.pi 
    
    # range -90 -- +90
    # arg = np.real(z[0]*np.conjugate(z[1]))
    # if arg < 0:
    #     azimuth = azm + 90
    # else:
    #     azimuth = azm - 90

    return azm


@jit(float32(complex128[:]))
def get_dip(z):
    """
    Code to give the incidence acoording to Park et al (1987).
    """
    abs_z1, phy1 = cm.polar(z[0])
    zH = np.sqrt(z[1]**2 + z[2]**2)
    abs_zH, phyH = cm.polar(zH)
    _ , phy1H = cm.polar(z[0]**2 + zH**2)

    th_m = -0.5*phy1H + np.arange(1,10)*np.pi/2
    func = abs_z1**2*np.cos(th_m + phy1)**2 + abs_zH**2*np.cos(th_m + phyH)**2
    th_V = th_m[np.argmax(func)]

    z1_rot = np.real(z[0] * np.exp(-1j*th_V))
    zH_rot = np.real(zH * np.exp(-1j*th_V))

    dip = np.arctan(np.abs(z1_rot/zH_rot))
    return dip*180/np.pi


@jit
def rotate(z, dgr):
    ans = np.empty((dgr.shape[0],), dtype=float32)
    for k in range(dgr.shape[0]):
        z_rot = z * np.complex(np.cos(dgr[k]), np.sin(dgr[k]))
        #z_rot = z * np.exp(np.complex(0, dgr[k]), dtype='complex128')
        a = np.array([x.real for x in z_rot])
        b = np.array([x.imag for x in z_rot])
        ans[k] = np.abs(np.vdot(a,b))
    return ans


def get_rectiliniarity(z):
    """
    This code returns the rectilinearity acoording to Schimmel and Gallard 2004.

    --exlpain--:
    First, we rotate the first eigenvector in the complex domain by a phase "phy" which
    produces a orthogonal relation between real vectors a and b. Then semimajor and
    semiminor axes of polarization ellipse are identify. Finally the rectiliniarity is
    computed as: 1 - |b|/|a|

    :param z: first eigenvector
    :type z: complex array
    :return: float
    """

    degree = np.arange(0, np.pi/2, 1e-2)
    ans = rotate(z, degree)
    ans_min = int(np.argmin(ans))
    phi_min = np.arange(0, np.pi/2, 1e-2)[ans_min]

    z_rot = z * np.exp(np.complex(0, phi_min))
    a = np.array([x.real for x in z_rot])
    b = np.array([x.imag for x in z_rot])
    a_norm = np.linalg.norm(a)
    b_norm = np.linalg.norm(b)

    if a_norm > b_norm:
        minor, major = b_norm, a_norm
    else:
        minor, major = a_norm, b_norm

    rectiliniarity = 1 - (minor/major)

    return rectiliniarity


def decompose(spec_matrix, full_return=False):
    vh_vector = np.empty((spec_matrix.shape[2],3), dtype=np.complex128)
    polar_dgr = np.empty((spec_matrix.shape[2],), dtype=np.float32)
    if full_return:
        s_matrix = np.empty((spec_matrix.shape[2],3), dtype=np.float32)
    
    def decompose_k(k):
        _, s, vh = np.linalg.svd(spec_matrix[:,:,k], full_matrices=False)
        beta2 = (3 * np.sum(s**2) - np.sum(s) ** 2) / (2 * np.sum(s) ** 2)
        polar_dgr[k] = beta2
        vh_vector[k,:] = vh[0,:]
        if full_return:
            s_matrix[k,:] = s

    # map(decompose_k, list(range(spec_matrix.shape[2])))

    for k in range(spec_matrix.shape[2]):
        _, s, vh = np.linalg.svd(spec_matrix[:,:,k], full_matrices=False)
        beta2 = (3 * np.sum(s**2) - np.sum(s) ** 2) / (2 * np.sum(s) ** 2)
        polar_dgr[k] = beta2
        vh_vector[k,:] = vh[0,:]
        if full_return:
            s_matrix[k,:] = s
    
    if full_return:
        return (polar_dgr, vh_vector, s_matrix)
    else:
        return (polar_dgr, vh_vector)


def polarization_content(spec_matrix, matrix_return=False, threshold=0.75):
    """
    Compute polarization degree from the cross-spectral matrix.

    :param spec_matrix: cross-correlation spectral hermitian complex matrix
    :type spec_matrix: 3x3xN complex numpy array
    :param threshold: minimum point to compute optional return
    :type threshold: float. default 0.75
    :param opt_return: if True, return azimuth, dip, and rectiliniarity
    :type opt_return: boolean
    :return: array with a float (polarization degree).
    """

    ans = decompose(spec_matrix)
    polar_dgr = ans[0]

    if matrix_return:
        # return a 2D array
        azm = np.zeros((polar_dgr.shape[0],), dtype='float32')
        dip = np.zeros((polar_dgr.shape[0],), dtype='float32')
        rect = np.zeros((polar_dgr.shape[0],), dtype='float32')
        for k in range(polar_dgr.shape[0]):
            if polar_dgr[k] >= threshold:
                z = ans[1][k,:]
                azm[k] = get_azimuth(z)
                dip[k] = get_dip(z)
                rect[k] = get_rectiliniarity(z)
            else:
                azm[k] = np.nan
                dip[k] = np.nan
                rect[k] = np.nan

    else:
        # return 1D array with max values
        if polar_dgr.max() >= threshold:
            k = np.argmax(polar_dgr)
            z = ans[1][k,:]
            azm = get_azimuth(z)
            dip = get_dip(z)
            rect = get_rectiliniarity(z)
        else:
            azm = np.nan
            dip = np.nan
            rect = np.nan

    return [polar_dgr, azm, dip, rect]



