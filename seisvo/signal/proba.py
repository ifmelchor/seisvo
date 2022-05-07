#!/usr/bin/env python3
# coding=utf-8

import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV


def get_PDF(array, x_array, bandwidth, **kwargs):
    """
    Compute the PDF

    Parameters
    ----------
    array : [ndarrray of shape (M,N)]
    x_array : [ndarrray of shape (L,1)]
    bandwidth : float or ndarrray
        [if bandwidth is ndarray it search the best score for all items]
    kernel : string
        [ default "gaussian"]

    Returns
    -------
    [ndarray of shape (L, N)]
        [description]
    """
    
    if not isinstance(bandwidth, float) and not isinstance(bandwidth, np.ndarray):
        raise ValueError ('bandwidth must be float or ndarray')

    pdf = np.empty((x_array.shape[0], array.shape[1]))

    for i in range(array.shape[1]):
        x_train = array[:,i].reshape(array.shape[0],1)
        kde = get_KDE(x_train, bandwidth, **kwargs)
        pdf[:, i] = np.exp(kde.score_samples(x_array))
    
    return pdf


def get_KDE(x, bandwidth, **kwargs):

    weight = kwargs.get('weight', None)
    kernel = kwargs.get('kernel', "gaussian")

    if isinstance(bandwidth, np.ndarray):
        kde = KernelDensity(kernel=kernel)
        grid = GridSearchCV(kde, {'bandwidth': bandwidth}, n_jobs=-1)
        grid.fit(x, sample_weight=weight)
        kde = grid.best_estimator_

        if kwargs.get('verbose', False):
            print(f' Best bandwidth: {kde.bandwidth}')
    
    else:
        kde = KernelDensity(kernel=kernel, bandwidth=bandwidth)
    

    if isinstance(weight, np.ndarray):
        kde.fit(x, sample_weight=weight)
    else:
        kde.fit(x)
    
    return kde


def get_pPSD(matrix, weights=None):
    from seisvo.signal.spectrum import get_centroid
    from scipy.stats import gaussian_kde
    
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
            dominant_gdk[f] = time_range[np.argmax(kde)]
            centroid_gdk[f] = get_centroid(time_range, kde)

        else:
            dominant_gdk[f] = centroid_gdk[f] = np.nan

    return [dominant_gdk, centroid_gdk]