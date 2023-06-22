#!/usr/bin/env python3
# coding=utf-8

import numpy as np
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from .utils import get_dominant, get_centroid


def azimuth_mean(x, ambiguity180=False):
    # convert to radian an computes sin(x)
    xrad = np.array(x)*np.pi/180
    
    if ambiguity180:
        xrad_mean = np.arctan(np.sin(xrad).sum()/np.abs(np.cos(xrad)).sum())
        xrad_std = np.sin(xrad).std()

    else:
        # the ambiguity is for 360 degree
        xrad_mean = np.arctan(np.sin(xrad).sum()/np.cos(xrad).sum())
        xrad_std = np.sin(xrad).std()

    return xrad_mean*180/np.pi, xrad_std*180/np.pi


def get_Stats(x, bw_method=None):
    """
    Return min, max, mean, and mode of a time series
    """
    v_min = x.min()
    v_max = x.max()
    v_mean = x.mean()

    x_range = np.linspace(v_min, v_max, 500)
    gkde = gaussian_kde(x, bw_method=bw_method)
    pdf = gkde(x_range)

    v_mode = get_dominant(x_range, pdf)

    return [v_min, v_max, v_mean, v_mode]


def get_PDF(array, space, **kde_kwargs):
    """
    Compute the PDF of a matrix of shape (M,N), where M is the nro of function of N length.

    Parameters
    ----------
    array : arrray of shape (M,N)
    space : arrray of shape (L,1)
    
    kde_kwargs ::
        >> bandwidth : (optional) float or ndarrray (optional), by default None
        >> bw_method : see gaussian_kde of scipy doc.
    
    help :: if bandwidth is an array, it search the best score for all items using GridSearch, 
            if None, it use gaussian_kde instead.

    Returns
    -------
    PDF array of shape (L, N)

    """

    # empty pdf
    pdf = np.empty((space.shape[0], array.shape[1]))

    for i in range(array.shape[1]):
        x = array[:,i].reshape(-1,1)
        kde = get_KDE(x, **kde_kwargs)

        if isinstance(kde, KernelDensity):
            pdf[:, i] = np.exp(kde.score_samples(space))
        else:
            pdf[:, i] = kde(space.reshape(-1,))
    
    return pdf


def get_KDE(x, bw_method=None, bandwidth=None, kernel="gaussian", weight=None, print_bestbw=False):
    """
    Compute KDE. IF bandwidth is specified, uses KernelDensity of scikit-learn, else it uses gaussian_kde of scipy.

    returns kde object

    """
    if bw_method or bandwidth == None:
        kde = gaussian_kde(x.reshape(-1,), bw_method=bw_method)
    else:
        if isinstance(bandwidth, np.ndarray):
            kde = KernelDensity(kernel=kernel)
            grid = GridSearchCV(kde, {'bandwidth': bandwidth}, n_jobs=-1)
            grid.fit(x, sample_weight=weight)
            kde = grid.best_estimator_
            if print_bestbw:
                print(f' Best bandwidth: {kde.bandwidth}')
        else:
            kde = KernelDensity(kernel=kernel, bandwidth=bandwidth)
    
        if isinstance(weight, np.ndarray):
            kde.fit(x, sample_weight=weight)
        else:
            kde.fit(x)
    
    return kde


# revisar esta function
def get_pPSD(matrix, weights=None):

    dominant_gdk = np.empty(matrix.shape[1])
    centroid_gdk = np.empty(matrix.shape[1])

    nro_timebins = matrix.shape[0]
    
    for f in range(matrix.shape[1]):
        m_f = matrix[:, f]
        m_f[np.isneginf(m_f)] = np.nan
        m_f_avg = np.nanmean(m_f)

        if weights is not None:
            wm_f = weights[:, f]

        # avoid spikes (with 5db > average)
        m_f[m_f > np.abs(5*m_f_avg)] = np.nan

        if weights is not None:
            wm_f[np.isnan(m_f)] = np.nan

        # remove nans
        if any(np.isnan(m_f)):
            m_f = m_f[np.isfinite(m_f)]

            if weights is not None:
                wm_f = wm_f[np.isfinite(wm_f)]

        if len(m_f) > nro_timebins/2:
            time_range = np.linspace(m_f.min(), m_f.max(), len(m_f)*5)
            
            if weights is not None:
                gkde = gaussian_kde(m_f, weights=wm_f)
            
            else:
                gkde = gaussian_kde(m_f)
            
            kde = gkde(time_range)
            dominant_gdk[f] = get_dominant(time_range, kde)
            centroid_gdk[f] = get_centroid(time_range, kde)

        else:
            dominant_gdk[f] = centroid_gdk[f] = np.nan

    return [dominant_gdk, centroid_gdk]