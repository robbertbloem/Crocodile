from __future__ import print_function
from __future__ import division

import inspect
import re

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.Debug as DEBUG
import Crocodile.Resources.Constants as CONST

def import_data_LV_A(path, base_filename):
    """
    Imports data for 'LV_file_format.1'
    
    """

    try:
        # data
        path_and_filename = path + base_filename + ".csv"
        data = numpy.loadtxt(path_and_filename, dtype = "float", delimiter = ",")
        data = data.T

        # time axis in fringes
        path_and_filename = path + base_filename + "_t1.csv"
        t1_axis = numpy.loadtxt(path_and_filename, dtype = "float", delimiter = ",") 

        # frequency axis
        path_and_filename = path + base_filename + "_w3.csv"
        w3_axis = numpy.loadtxt(path_and_filename, dtype = "float", delimiter = ",") 

        # phase
        path_and_filename = path + base_filename + "_phase.txt"
        f = open(path_and_filename)
        for line in f:
            temp = line
        f.close()
        if temp == "NaN":
            DEBUG.printWarning("Phase is Not-a-Number, will be set to 0 ", inspect.stack())
            phase = 0
        else:
            phase = float(temp)

        # last pump
        path_and_filename = path + base_filename + "_lastpump.txt"
        f = open(path_and_filename)
        for line in f:
            lastpump = line
        f.close()

        # determine number of fringes
        n_fringes = int((len(t1_axis)+1)/2)
        n_pixels = len(w3_axis)

        # convert NaN to zeros
        data = numpy.nan_to_num(data)

        # labview measures 4000-N to 4000+N, we want the data split into 4000-N to 4000 (non-rephasing) and 4000 to 4000+N (rephasing)
        R = data[n_fringes-1:, :]
        NR = numpy.flipud(data[:n_fringes, :])

        # for the FFT, we don't want 4000 to be zero. The axes run from 0 to N
        # also: calculate the axis in fs        
        t1fr_axis = numpy.arange(n_fringes)
        t1fs_axis = numpy.arange(n_fringes) * CONST.hene_fringe_fs

        # return everything
        return R, NR, t1fs_axis, t1fr_axis, w3_axis, phase, lastpump, n_fringes, n_pixels

    except IOError:
        DEBUG.printError("Unable to import LabView data from file " + path + base_filename, inspect.stack())
        # raise
        return False
        





 















 