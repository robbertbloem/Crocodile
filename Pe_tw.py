from __future__ import print_function
from __future__ import division

import inspect
import os

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe as PE
import Crocodile.Resources.IOMethods as IOM
import Crocodile.Resources.Constants as CONST

reload(PE)
reload(IOM)

class pe_tw(PE.pe):
    """
    pe_tw: implementation class for photo echo time-frequency. This includes all importing functions.
    
    the VB6-stuff is found in resources/Pe_tw_vb6.py. Good riddance! 
    
    """

    def __init__(self, objectname, root_filename, time_stamp, population_time, undersampling, flag_verbose = False):

        self.verbose("New Pe.pe.pe_tw class", flag_verbose)

        PE.pe.__init__(self, objectname, flag_verbose = flag_verbose)


        self.r_axis[1] = population_time
        self.time_stamp = time_stamp 
        self.base_filename = root_filename + "_" + str(self.time_stamp) + "_T" + str(self.r_axis[1])
        self.sub_type = root_filename
        
        self.undersampling = undersampling
        self.path = self.base_filename + "/"  
        
        # somewhere else? 
        self.n_pixels = 32
        self.reference = numpy.ones(self.n_pixels)
        self.data_type_version = 0
        
        self.n_shots = 0
        self.n_steps = 0
        



    def absorptive(self, window_function = "none", window_length = 0, flag_verbose = False):

        self.super_absorptive(axes = [1,0], window_function = window_function, window_length = window_length, flag_verbose = flag_verbose)
    
    



class pe_LV(pe_tw):
    """
    Class to process results from LabView.
    

    
    
    
    """
    

    def __init__(self, objectname, root_filename, time_stamp, population_time, flag_verbose = False):
    
        self.verbose("New Pe.pe.pe_LV class", flag_verbose)
    
        pe_tw.__init__(self, objectname = objectname, root_filename = root_filename, time_stamp = time_stamp, population_time = population_time, undersampling = 0, flag_verbose = flag_verbose)  
        
        self.n_fringes = 0 
        self.undersampling = 0
    

    def import_data(self, flag_verbose = False):
        """
        There can be different file formats, indicated by the file 'LV_file_format.N', where N is the format. The file format should be changed when a new or adapted import method is needed. 
        For example, the format '1' writes data as .csv, while '2' writes data as .bin. Format '1' needs import method 'IOM.import_data_LV_A'. Format 2 can use a modified version of that (add an extra flag indicating .bin files, it should default to .csv to maintain compatibility with older scripts) or make a completly new method 'IOM.import_data_LV_B'. 
        Note that the 'LV_file_format.N' are numbered by numbers, while the versions of 'IOM.import_data_LV_X' are indicated by letters.
        'LV_file_format.666' is reserved for testing purposes.
        
        """
        
        self.verbose("Import data", flag_verbose)
        
        try:
            filelist = os.listdir(self.path)
        except OSError:
            self.printError("Directory can not be found.", inspect.stack())
            return False
        
        if "LV_file_format.1" in filelist:
            self.verbose("  LabView file format 1", flag_verbose)
            res = IOM.import_data_LV_A(self.path, self.base_filename)
            if res:
                R, NR, t1fs_axis, t1fr_axis, w3_axis, phase, lastpump, n_fringes, n_pixels = res
            else:
                # res == False means error in importing
                return False
            
        else:
            self.printError("Unknown file format. Please check 'LV_file_format.N' is present.", inspect.stack())
            return False
        
        self.r[0] = R
        self.r[1] = NR
        self.r_axis[0] = t1fs_axis
        self.r_axis[2] = w3_axis
        self.s_axis[2] = w3_axis
        self.phase_degrees = phase
        
        self.comment = lastpump
        
        if self.n_fringes == 0:
            self.n_fringes = n_fringes
        if self.n_pixels == 0:
            self.n_pixels = n_pixels
        
        self.r_units = ["fs", "fs", "cm-1"]
        
        self.n_scans += 1
        
        return True
        







  



