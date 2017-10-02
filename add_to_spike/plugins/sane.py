#!/usr/bin/env python 
# encoding: utf-8

"plugin for Sane and PG_Sane "

from __future__ import print_function

import unittest
import numpy as np
from numpy.fft import fft, ifft, rfft, irfft

from spike.NPKData import NPKData_plugin,  as_cpx, as_float, _base_fft,\
            _base_ifft, _base_rfft, _base_irfft
from spike.Algo.sane import sane
from spike.util.signal_tools import filtering

import sys #
if sys.version_info[0] < 3:
    pass
else:
    xrange = range

def sane_plugin(npkd, rank, orda=None, iterations=1, axis=0, trick=True, optk=False, ktrick=False):
    """
    Apply "sane" denoising to data
    rank is about 2 x number_of_expected_lines
    Manages real and complex cases.
    Handles the case of hypercomplex for denoising of 2D FTICR for example.
    
    sane algorithm. Name stands for Support Selection for Noise Elimination.
    From a data series return a denoised series denoised
    data : the series to be denoised - a (normally complex) numpy buffer
    rank : the rank of the analysis
    orda : is the order of the analysis
        internally, a Hankel matrix (M,N) is constructed, with M = orda and N = len(data)-orda+1
        if None (default) orda = (len(data)+1)/2
    iterations : the number of time the operation should be repeated
    optk : if set to True will calculate the rank giving the best recovery for an automatic estimated noise level. 
    trick : permits to enhanced the denoising by using a cleaned signal as the projective space. "Support Selection"
    ktrick : if a value is given, it permits to change the rank on the second pass.
             The idea is that for the first pass a rank large enough as to be used to compensate for the noise while
             for the second pass a lower rank can be used. 
    
    """
    if npkd.dim == 1:
        if npkd.axis1.itype == 0:   # real
            buff = as_cpx(_base_ifft(_base_rfft(npkd.buffer)))       # real case, go to analytical signal
        else:   #complex
            buff = npkd.get_buffer()                       # complex case, makes complex
        sane_result = sane( buff, rank, orda = orda, trick = trick, iterations = iterations) # performs denoising
        if npkd.axis1.itype == 0:   # real
            buff = _base_irfft(_base_fft(as_float(sane_result)))      # real case, comes back to real
            npkd.set_buffer(buff)
        else:
            npkd.buffer = as_float(sane_result)             # complex case, makes real
    elif npkd.dim == 2:
         todo = npkd.test_axis(axis)
         if todo == 2:
             for i in xrange(npkd.size1):
                 r = npkd.row(i).sane(rank=rank, orda=orda, iterations=iterations)
                 npkd.set_row(i,r)
         elif todo == 1:
             for i in xrange(npkd.size2):
                 r = npkd.col(i).sane(rank=rank, orda=orda, iterations=iterations)
                 npkd.set_col(i,r)
    elif npkd.dim == 3:
         raise Exception("not implemented yet")
    return npkd



NPKData_plugin("sane", sane_plugin)
