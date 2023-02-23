#!/usr/bin/env python3
# coding=utf-8

import julia
from .cc8 import CC8
from .plotting import plot_array_response, plot_array

class Arfr(object):
    def __init__(self, sarray):
        self.sar_  = sarray
        self.arf_  = None
    

    def fit(self, slow_max, slow_inc, fq_band=(1., 10.), fq_int=0.1):
        julia.Julia(compiled_modules=False)
        from julia import SAP

        self.fq_band_  = fq_band
        self.slow_max_ = float(slow_max)
        self.slow_inc_ = float(slow_inc)
        self.arf_      = SAP.array_response(self.sar_.xUTM, self.sar_.yUTM, self.slow_max_, self.slow_inc_, float(fq_band[0]), float(fq_band[1]), float(fq_int))

        return self


    def plot(self, save=False, filename="./arf.png"):
        if self.arf_:
            fig = plot_array_response(self)

            if save:
                fig.savefig(filename)



