from __future__ import print_function
from __future__ import division

import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe as PE


class pe_tw(PE.pe):
    """
    pe_tw: implementation class for photo echo time-frequency. This includes all importing functions.
    
    """

    def __init__(self, objectname, base_filename, time_stamp, population_time, undersampling, flag_verbose = False):

        self.verbose("New Pe.pe.pe_tw class", flag_verbose)

        PE.pe.__init__(self, objectname, flag_verbose = flag_verbose)

        self.base_filename = base_filename
        self.r_axis[1] = population_time
        self.time_stamp = time_stamp 
        self.undersampling = undersampling
        
        # not the best place...
        self.path = self.base_filename + "_" + str(self.time_stamp) + "_T" + str(self.r_axis[1]) + "/"   
        self.n_pixels = 32  


    def absorptive(self, window_function = "none", window_length = 0):

        self.super_absorptive(axes = [1,0], window_function = "none", window_length = 0)
        




