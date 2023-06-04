#!/usr/bin/env python3
# coding=utf-8

import numpy as np


def get_CC8(data, fs, xutm, yutm, fq_band, slow_max, slow_inc, **kwargs):

    ans = {}
    if isinstance(data, np.ndarray):
        #assert len(xutm) == len(yutm) == data.shape[0]
        #assert len(fq_band) == 2
        #assert len(slow_max) == len(slow_inc)
        nites = len(slow_max)

        # convert to array
        data = jl.Array(np.array(data))
        xutm = jl.Array(np.array(xutm))
        yutm = jl.Array(np.array(yutm))
        fq_band = jl.Array(np.array(fq_band))
        slow_max = jl.Array(np.array(slow_max))
        slow_inc = jl.Array(np.array(slow_inc))

        # get kwargs
        lwin = int(kwargs["lwin"])
        nwin = int(kwargs["nwin"])
        nadv = float(kwargs["nadv"])
        cc_thres = float(kwargs["cc_thres"])
        toff = int(kwargs["toff"])

        # run in julia
        from juliacall import Main as jl
        jl.seval("using SAP")
        jlans = jl.CC8(data, xutm, yutm, slow_max, slow_inc, fq_band, int(fs), lwin, nwin, nadv, cc_thres, toff)

        # convert to python dict
        pyans = dict(jlans)
        for nsi in range(1, nites+1):
            ans[nsi] = {}
            for attr in ("slow", "bazm", "maac", "rms", "slowmap", "slowbnd", "bazmbnd"):
                    ans[nsi][attr] = np.array(pyans[nsi][attr])

    return ans


def array_response(xUTM, yUTM, slow_max, slow_inc, fq_band=(1.,10.), fq_int=0.1):
    
    from juliacall import Main as jl
    jl.seval("using SAP")

    assert len(xUTM) == len(yUTM)

    slow_max_ = float(slow_max)
    slow_inc_ = float(slow_inc)

    x = jl.Array(xUTM)
    y = jl.Array(yUTM)
    jlarf  = jl.array_response(x, y, slow_max_, slow_inc_, float(fq_band[0]), float(fq_band[1]), float(fq_int))
    
    ans = {
	    "freq":np.array(jlarf.freqs),
	    "power":np.array(jlarf.power),
	    "sx":np.array(jlarf.sx),
	    "sy":np.array(jlarf.sy)
    }

    return ans