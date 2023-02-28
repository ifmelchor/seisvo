#!/usr/bin/env python3
# coding=utf-8

from .cc8 import np, CC8
from .plotting import plot_array_response, plot_array


class Arfr(object):
    def __init__(self, sarray):
        self.sar_  = sarray
        self.arf_  = None
    

    def fit(self, slow_max, slow_inc, fq_band=(1., 10.), fq_int=0.1):
        from juliacall import Main as jl
        jl.seval("using SAP")

        self.fq_band_  = fq_band
        self.slow_max_ = float(slow_max)
        self.slow_inc_ = float(slow_inc)

        x = jl.Array(self.sar_.xUTM)
        y = jl.Array(self.sar_.yUTM)
        jlarf  = jl.array_response(x, y, self.slow_max_, self.slow_inc_, float(fq_band[0]), float(fq_band[1]), fq_int)
        
        self.freqs_ = np.array(jlarf.freqs) 
        self.power_ = np.array(jlarf.power)
        self.sx_    = np.array(jlarf.sx)
        self.sy_    = np.array(jlarf.sy)
        self.arf_   = True

        return self


    def plot(self, save=False, filename="./arf.png"):
        if self.arf_:
            fig = plot_array_response(self)

            if save:
                fig.savefig(filename)



