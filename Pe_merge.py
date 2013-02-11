from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe as PE
import PythonTools.Debug as DEBUG

class pe_merge(PE.pe):
    
    def __init__(self, objectname, class_plus = [], class_min = [], flag_verbose = False):
        
        PE.pe.__init__(self, objectname, flag_verbose = flag_verbose)

        # self.obj_id = objectname
        # self.objectname = objectname

        # first check if it makes sense to continue
        # at least r or s should have the same dimensions
        flag_r = True
        flag_s = True
        r_size = False
        s_size = False
        for obj in class_plus + class_min:
            
            if hasattr(obj, "s"):
                temp = numpy.shape(obj.s)
                if temp == (1,):
                    flag_s = False
                elif s_size:
                    if s_size != temp:
                        flag_s = False
                else:
                    s_size = temp
            else:
                flag_s = False

            if hasattr(obj, "r"):
                temp = numpy.shape(obj.r)
                if temp == (2,):
                    flag_r = False
                elif r_size:
                    if r_size != temp:
                        flag_r = False
                else:
                    r_size = temp
            else:
                flag_r = False            

        if flag_r == False and flag_s == False:
            return False

        n_p = 0
        n_m = 0   
        for obj in class_plus:
            n_p += obj.n_scans

        for obj in class_min:
            n_m += obj.n_scans
        
        
        # add or subtract the data
        # frankly, it is stupid to set the axis all the time
        for obj in class_plus:
            if flag_r:
                for i in range(len(obj.r)):
                    self.r[i] += obj.r[i] * obj.n_scans / n_p
                for i in range(len(obj.r_axis)):
                    self.r_axis[i] = obj.r_axis[i]
            if flag_s:   
                self.s += obj.s * obj.n_scans / n_p
                for i in range(len(obj.s_axis)):
                    self.s_axis[i] = obj.s_axis[i]

        for obj in class_min:       
            if flag_r:
                for i in range(len(obj.r)):
                    self.r[i] -= obj.r[i] * obj.n_scans / n_m
                for i in range(len(obj.r_axis)):
                    self.r_axis[i] = obj.r_axis[i]
            if flag_s:   
                self.s -= obj.s * obj.n_scans / n_m
                for i in range(len(obj.s_axis)):
                    self.s_axis[i] = obj.s_axis[i]

        
        # SKIPPED VARIABLES:
        # phase_rad: set by phase_degrees
        # zeropad_by: set by zeropad_to
        # b, b_axis, b_count: would be stupid
        # base_filename: no meaning without measurement
        # dimensions: already set by Pe class
        # f, f_axis: why would you?
        # measurements: already set by Pe class
        # obj_id, objectname: already set by init
        # r, r_axis: set above
        # source_path: no meaning without measurement
        # sub_type: set manually
        # _time_stamp: set by init
    
        # set some variable are initialized to an array. Set to False so that     getattr(obj_to, key) == False
        self.r_units = False
        self.s_resolution = False
        self.s_units = False
        
        index = 0
        for obj in class_plus + class_min:
            # for stuff that independent for adding or subtracting
            
            for key in sorted(self.__dict__):   

                # three variables with setter methods. Ie. the _ has to be removed. Comments are appended.
                if key == "_comment":
                    self.comment = getattr(obj, key)
                elif key == "_phase_degrees":
                    check_value_set_key(self, obj, key = "phase_degrees", index = index, flag_verbose = flag_verbose)
                elif key == "_zeropad_to":
                    if flag_r:
                        check_value_set_key(self, obj, key = "zeropad_to", index  = index, flag_verbose = flag_verbose)

                # normal variables, related to r
                elif key in ["r_correction", "r_correction_applied", "r_resolution", "r_units", "undersampling"]:
                    if flag_r:
                        check_value_set_key(self, obj, key, index, flag_verbose = flag_verbose) 
                # normal variables, related to s
                elif key in ["s_resolution", "s_units"]:
                    if flag_s:
                        check_value_set_key(self, obj, key, index, flag_verbose = flag_verbose)  
                # normal variables, other
                elif key in ["date"]:
                    check_value_set_key(self, obj, key, index, flag_verbose = flag_verbose)                 
                # skipped variables
                elif key in ["s", "s_axis", "r", "r_axis", "_phase_rad", "_zeropad_by", "b", "b_axis", "b_count", "base_filename", "date", "dimensions", "f", "f_axis", "measurements", "n_scans", "obj_id", "objectname", "r", "r_axis", "source_path", "sub_type", "_time_stamp"]:
                    pass
                
                # unknown variables
                else:
                    self.printWarning("Unknown key: " + key)
            
            index += 1



def check_value_set_key(obj_to, obj_from, key, index, flag_verbose = False):
    """
    checks if a key is set. If it is set, it checks if the values are the same. If not, it gives an error and gives the variable a value numpy.nan. If the value is not set and index == 0, it will set the value. 
    
    CHANGELOG:
    20130208/RB: started the function
    
    """  
    # if key is already set, check if it is the same
    if getattr(obj_to, key):
        if getattr(obj_to, key) != getattr(obj_from, key):
            # if value is different, set to numpy.nan
            # make sure getters/setters can handle this
            setattr(obj_to, key, numpy.nan)
            DEBUG.verbose(key + " does not match. Set to NaN", flag_verbose = flag_verbose)
    elif index == 0:
        # if key was not set, set it now
        setattr(obj_to, key, getattr(obj_from, key))
    
    
    
    


        
        
