#!/usr/bin/env python3
# coding=utf-8

import functools
import numpy as np
import datetime as dt

infrasound_model_default = {
    'radii': 10000,
    'vel_air': 342,
    'h_src': 1000,
    'src_dgr': (0,359)
}


def corr_coeff_index(nro):
    """
    Esta function me devuelve los indices de correlacion cruzada entre pares de sensores
    :return: lista de pares de tuplas
    """
    i = 0
    list1 = []
    list2 = []
    nro_sensors = range(nro)
    while i < nro:
        list1 += [nro_sensors[i]] * (nro - (i + 1))
        list2 += [nro_sensors[x - 1] for x in range(i + 2, nro + 1)]
        i += 1

    return zip(list1, list2)


def cross_corr(stream, dtimes, nro_srcs, ss, **kwargs):
    '''
    Compute cross correlation between pair of stations
    :param stream:
    :param dtimes:
    :param nro_srcs:
    :param time_width:
    :param overlap:
    :return:
    '''

    air_file = kwargs.get('air_file', None)
    last_interval = kwargs.get('last_interval', None)
    iter_bin = kwargs.get('iter_bin', None)
    pre_txt = kwargs.get('txt_print', '')

    if len(set([tr.stats.npts for tr in stream])) > 1:
        print('\n\n')
        print(stream)
        raise ValueError ('\n error: npts not coincide')

    nro_sensors = len(stream)
    stats = stream[0].stats
    npt_delta = int(stats.sampling_rate * ss.window_length)
    npt_olap = int(stats.sampling_rate * ss.overlap)

    cross_time_matrix = False
    times = []

    if last_interval:
        nro_steps = ss.laststeps_
    
    else:
        nro_steps = ss.steps_

    npt_start = npt_delta # esto es así porque hemos añadido una ventana adicional para el solapamiento!
    for kv in range(nro_steps):
        control_txt = pre_txt + f'{int(100*(kv/nro_steps)):>7}%'
        print(control_txt, end="\r")

        cross_matrix = np.ndarray((1, nro_srcs))
        press_avg = np.ndarray((1, nro_srcs))
        press_max = np.ndarray((1, nro_srcs))

        for k in range(nro_srcs):
            matrix = np.ndarray((1, nro_sensors, npt_delta))

            for c, trace in enumerate(stream):
                sta_name = '%s.%s.%s' % (trace.stats.network, trace.stats.station, trace.stats.location)
                npts_c_in = npt_start + dtimes[sta_name][k]
                npts_c_fi = npts_c_in + npt_delta
                matrix[0, c] = trace.data[npts_c_in:npts_c_fi]
            
            signal_avg = matrix.mean(axis=1)
            signal_avg -= signal_avg.mean()
            press_avg[0, k] = np.abs(signal_avg).mean()
            press_max[0, k] = np.abs(signal_avg).max()

            list_corr = []
            
            # esto se puede paralelizar con multiprocessing si vemos que tarda mucho
            crosscorr_func = functools.partial(__croscormat__, matrix)
            list_corr = list(map(crosscorr_func, corr_coeff_index(nro_sensors)))
            array_corr = np.array(list_corr)
            
            # for j, i in corr_coeff_index(nro_sensors):
            #     mx = np.corrcoef(matrix[0, j], matrix[0, i])[0, 1]
            #     if mx < 0:
            #         mx = 0
            #     list_corr += [mx]

            if np.count_nonzero(np.isnan(array_corr)) > 1:
                print(' Warning: NaN in the correlation approach')

            cross_matrix[0, k] = np.nanmean(array_corr)

        if air_file:
            item_bin = next(iter_bin)
            air_file.save_data("crosscorr", 'crosscorr', cross_matrix, item=item_bin)
            air_file.save_data("pressure", 'p_avg', press_avg, item=item_bin)
            air_file.save_data("pressure", 'p_max', press_max, item=item_bin)

        else:
            if cross_time_matrix is False:
                cross_time_matrix = cross_matrix
                press_time_avg = press_avg
                press_time_max = press_max
            else:
                cross_time_matrix = np.concatenate((cross_time_matrix, cross_matrix), axis=0)
                press_time_avg = np.concatenate((press_time_avg, press_avg), axis=0)
                press_time_max = np.concatenate((press_time_max, press_max), axis=0)

            t = stats.starttime + dt.timedelta(seconds=npt_start/stats.sampling_rate)
            times.append(t.datetime)
        
        npt_start += npt_delta - npt_olap
    
    if not air_file:
        return (times, cross_time_matrix, press_time_avg, press_time_max)


def __croscormat__(matrix, idx):
    i, j = idx
    mx = np.corrcoef(matrix[0, j], matrix[0, i])[0, 1]
    if mx < 0:
        mx = 0
    return mx