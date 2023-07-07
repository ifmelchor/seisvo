#!/usr/bin/env python3
# coding=utf-8

import numpy as np


def get_CC8(data, fs, xutm, yutm, fq_band, slow_max, slow_inc, **kwargs):
    """
    julia wrapper for the processing of the CC8 algorithm
    """

    dictans = {}
    
    if isinstance(data, np.ndarray):
        #assert len(xutm) == len(yutm) == data.shape[0]
        #assert len(fq_band) == 2
        #assert len(slow_max) == len(slow_inc)
        data = data.astype(dtype=np.float64)

        # convert to array
        from juliacall import Main as jl
        data = jl.Array(np.array(data))
        xutm = jl.Array(np.array(xutm))
        yutm = jl.Array(np.array(yutm))
        fq_band = jl.Array(np.array(fq_band))
        slow_max = jl.Array(np.array([slow_max]))
        slow_inc = jl.Array(np.array([slow_inc]))

        # get kwargs
        lwin = int(kwargs["lwin"])
        nwin = int(kwargs["nwin"])
        nadv = float(kwargs["nadv"])
        cc_thres = float(kwargs["cc_thres"])
        toff = int(kwargs["toff"])

        # load julia
        jl.seval("using SAP")
        
        try:
            # run and save
            jlans = jl.CC8(data, xutm, yutm, slow_max, slow_inc, fq_band, int(fs), lwin, nwin, nadv, cc_thres, toff)
            pyans = dict(jlans)
            dictans = {}
            for attr in ("slow", "bazm", "maac", "rms", "slowmap", "slowbnd", "bazmbnd"):
                dictans[attr] = np.array(pyans[1][attr])

        except Exception as exc:
            print("\n ------ ERROR INFO ------")
            print(exc)
            print(" ------------------------\n")

    return dictans


def array_delta_times(slowness, bazimuth, slomax, slomint, xUTM, yUTM, etol=1e-2):

    from juliacall import Main as jl
    jl.seval("using SAP")

    assert len(xUTM) == len(yUTM)

    etol      = float(etol)
    slow_max_ = float(slomax)
    slow_inc_ = float(slomint)
    slowness  = float(slowness)
    bazimuth  = float(bazimuth)
    utm_east  = jl.Array(np.array(xUTM))
    utm_north = jl.Array(np.array(yUTM))
    deltas, tol  = jl.get_dtimes(slowness, bazimuth, slow_max_, slow_inc_,utm_east, utm_north, etol)

    return np.array(deltas), np.array(tol)


def array_response(xUTM, yUTM, slow_max, slow_inc, fq_band=(1.,10.), fq_int=0.1):
    
    from juliacall import Main as jl
    jl.seval("using SAP")

    assert len(xUTM) == len(yUTM)

    slow_max_ = float(slow_max)
    slow_inc_ = float(slow_inc)

    x = jl.Array(np.array(xUTM))
    y = jl.Array(np.array(yUTM))
    jlarf  = jl.array_response(x, y, slow_max_, slow_inc_, float(fq_band[0]), float(fq_band[1]), float(fq_int))
    
    ans = {
	    "freq":np.array(jlarf.freqs),
	    "power":np.array(jlarf.power),
	    "sx":np.array(jlarf.sx),
	    "sy":np.array(jlarf.sy)
    }

    return ans