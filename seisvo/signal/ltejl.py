#!/usr/bin/python3
# coding=utf-8

import numpy as np
from .utils import get_freq, SSteps


def get_PSD(data, fs, lwin=None, olap=0., fq_band=(0.5,15), NW=3.5, pad=1.0, full_return=False):

    """
    Compute the PSD using multitaper algorithm.
    Requires julia packatge LTE.jl

    for apply moving average, use lwin and olap:
        >> lwin in samples
        >> olap between 0 and 1 is the percent of overlap
    
    """
    
    if lwin:
        assert lwin < len(data)
        nsteps = np.floor(SSteps.nsteps(int(len(data)/fs), int(lwin/fs), olap))
        npts = lwin = int(lwin)
        nwin = int(nsteps)
        nadv = float(1-olap)
        assert nadv > 0
    else:
        npts = len(data)
        nwin = None
        lwin = None
        full_return = False

    freq, fqr = get_freq(npts, fs, fq_band=fq_band, pad=pad)
    fqr = np.array(fqr)

    # import LTE in julia
    from juliacall import Main as jl
    jl.seval("using LTE")

    if full_return:
        psd, _ = jl.LTE._full_psd(jl.Array(data), int(fs),\
            lwin, nwin, nadv, jl.Array(fqr), float(NW), float(pad))
    else:
        psd, _ = jl.LTE._psd(jl.Array(data), int(fs),\
            lwin, nwin, nadv, jl.Array(fqr), float(NW), float(pad))

    return np.array(psd), freq


def get_CSM(data, fs, lwin=None, olap=0., fq_band=(0.5,15), NW=3.5, pad=1.0,):
    """
    Compute the CrossSpectralMatrix (CSM) using multitaper algorithm.
    Requires julia packatge LTE.jl

    for apply moving average, use lwin and olap:
        >> lwin in samples
        >> olap between 0 and 1 is the percent of overlap
    
    """

    if lwin:
        assert lwin < len(data)
        npts = lwin
        nsteps = np.floor(SSteps.nsteps(int(len(data)/fs), int(lwin/fs), olap))
        nwin = int(nsteps)
        nadv = float(1-olap)
        assert nadv > 0
    else:
        npts = len(data)
        nwin = None
        lwin = None

    freq, fqr = get_freq(npts, fs, fq_band=fq_band, pad=pad)

    # import LTE in julia
    from juliacall import Main as jl
    jl.seval("using LTE")

    # compute SVD of CrossSpectralMatrix (CSM)
    csm = jl.LTE._csm(jl.Array(data), int(fs), lwin, nwin, nadv, jl.Array(fqr), float(NW), float(pad))

    # for CSW, you should apply spectral normalization before CSM
    # and then compute the CSW from CSM

    return freq, csm


def get_Polar(data, fs, lwin=None, olap=0, fq_band=(1., 5.), NW=3.5, pad=1.0, return_all=False, full_return=True):
    """
    Compute the polarization analysis using multitaper algorithm.
    Requires julia packatge LTE.jl

    for apply moving average, use lwin and olap:
        >> lwin in samples
        >> olap between 0 and 1 is the percent of overlap
    
    """

    ncomp, tnpts = data.shape
    assert ncomp == 3

    if lwin:
        assert lwin < tnpts
        nsteps = np.floor(SSteps.nsteps(int(tnpts/fs), int(lwin/fs), olap))
        npts = lwin = int(lwin)
        nwin = int(nsteps)
        nadv = float(1-olap)
        assert nadv > 0
    else:
        npts = len(data)
        nwin = None
        lwin = None

    freq, _ = get_freq(npts, fs, fq_band=fq_band, pad=pad)

    # import LTE in julia
    from juliacall import Main as jl
    jl.seval("using LTE")

    _, polar = jl.polar_run(jl.Array(data), jl.Tuple(fq_band),\
         int(fs), NW, pad, nwin, lwin, nadv, return_all, full_return)

    if full_return:
        polar = {
            "degree" :np.array(polar.degree),
            "rect"   :np.array(polar.rect),
            "azimuth":np.array(polar.azimuth),
            "elev"   :np.array(polar.elev),
            "phyHH"  :np.array(polar.phyhh),
            "phyVH"  :np.array(polar.phyvh)
        }

    else:
        # polar is polarization degree
        polar = np.array(polar)

    return freq, polar


def get_Peaks(chan_list, sxx, degree, rect, azimuth, elev, freq, peak_thresholds):

    assert len(chan_list) == sxx.shape[0]
    assert len(freq) == sxx.shape[2]
    
    # import LTE
    from juliacall import Main as jl
    jl.seval("using LTE")

    # define PeakThresholds type
    pksth = jl.PeakThresholds(
        float(peak_thresholds["fq_dt"]),
        float(peak_thresholds["sxx_th"]),
        float(peak_thresholds["pdg_th"]),
        float(peak_thresholds["pdg_std"]),
        float(peak_thresholds["rect_th"]),
        float(peak_thresholds["rect_std"]),
        float(peak_thresholds["azim_std"]),
        float(peak_thresholds["elev_std"]),
        int(peak_thresholds["npts_min"]))
    
    # extract the peaks
    sxx = jl.Array(sxx)
    freq = jl.Array(freq)
    degr = jl.Array(degree)
    rect = jl.Array(rect)
    azim = jl.Array(azimuth)
    elev = jl.Array(elev)
    jlout = jl.extract(sxx, degr, rect, azim, elev, freq, pksth)

    # convert to a python object
    dpeaks = {}
    for c, chan in enumerate(chan_list):
        dpeaks[chan] = {"nW":jlout[1][0][c], "nL":jlout[1][1][c], "pks":{}}
        
        for t, pks in dict(jlout[0][c]).items():
            dpeaks[chan]["pks"][t] = {
                "fq"  :np.array(dict(pks)["fq"],dtype=np.float32),
                "sxx" :np.array(dict(pks)["specgram"],dtype=np.float32),
                "dgr" :np.array(dict(pks)["degree"],dtype=np.float32),
                "rect":np.array(dict(pks)["rect"],dtype=np.float32),
                "azim":np.array(dict(pks)["azimuth"],dtype=np.float32),
                "elev":np.array(dict(pks)["elev"],dtype=np.float32)
            }
    
    return dpeaks


def get_LTE(data, fs, chan_list, fq_band, **lte_dict):

    lte_ans = {}

    if isinstance(data, np.ndarray):
        # load julia function
        from juliacall import Main as jl
        jl.seval("using LTE")
        chan = "/".join(chan_list)
        data = jl.Array(data)
        band = jl.Tuple(fq_band)
        fs = int(fs)

        if lte_dict["type"] == "station":
            nwin  = int(lte_dict["nwin"])
            lwin  = int(lte_dict["lwin"])
            nswin = int(lte_dict["nswin"])
            lswin = int(lte_dict["lswin"])
            nadv  = float(lte_dict["nadv"])
            NW    = float(lte_dict["time_bandwidth"])
            pad   = float(lte_dict["pad"])
            optp  = lte_dict["opt_params"]
            polar = lte_dict["polar"]
            peo   = int(lte_dict["PE_order"])
            pet   = int(lte_dict["PE_tau"])
            optw  = float(lte_dict["opt_twin"])
            opth  = float(lte_dict["opt_th"])

            # run in julia
            ans = jl.sta_run(data, chan, fs, nwin, lwin, nswin, lswin, nadv, band, NW, pad, optp, polar, peo, pet, optw, opth)
            
            try:
                jlans = dict(ans)
                # convert the jl dict into a python dict
                for chan in chan_list:
                    lte_ans[chan] = {}
                    lte_ans[chan]["specgram"] = np.array(jlans[chan]["specgram"])

                    for attr in ("perm_entr", "energy", "fq_dominant", "fq_centroid"):
                        lte_ans[chan][attr] = np.array(jlans[chan][attr])
                    
                    if optp:
                        lte_ans["opt"] = {}
                        for attr in ("vlf", "lf", "vlar", "rsam", "lrar", "mf", "rmar", "hf", "dsar"):
                            lte_ans["opt"][attr] = np.array(jlans["opt"][attr])
                    
                    if polar:
                        lte_ans["polar"] = {}
                        for attr in ("degree", "rect", "azimuth", "elev", "phyhh", "phyvh"):
                            lte_ans["polar"][attr] = np.array(jlans["polar"][attr])
            except:
                # you can catch error here and do whatever
                pass
        
        if lte_dict["type"] == "network":
            # nwin  = int(lte_dict["nwin"])
            # lwin  = int(lte_dict["lwin"])
            nswin = int(lte_dict["nswin"])
            lswin = int(lte_dict["lswin"])
            nadv  = float(lte_dict["nadv"])
            NW    = float(lte_dict["time_bandwidth"])
            pad   = float(lte_dict["pad"])
            ans = jl.net_run(data, chan, fs, nswin, lswin, nadv, band, NW, pad)

            try:
                jlans = dict(ans)
                for sta in chan_list:
                    lte_ans[sta] = np.array(jlans[sta]).reshape(1,-1)
            
                lte_ans["csw"] = np.array(jlans["csw"]).reshape(1,-1)
                lte_ans["vt"]  = np.array(jlans["vt"]).reshape(1, -1, len(chan_list))
            except:
                # you can catch error here and do whatever
                pass

    return lte_ans

