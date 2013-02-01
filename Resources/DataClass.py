from __future__ import print_function
from __future__ import division

import inspect
import time

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile
import Crocodile.Resources.ClassTools as CT

class dataclass(CT.ClassTools):
    """
    This class stores all the data. 
    Most of the variables are lists which point to ndarrays with the real data.
    """

    def __init__(self, objectname, dimensions, measurements, flag_verbose = False):
        """
        croc.DataClasses.messdata

        INPUT:
        - objectname: the name
        - dimensions: the number of axes. Linear/2D/3D = 1/2/3. If something is measured in frequency and time (certain pump-probe) the dimension is also 2.
        - measurements: number of diagrams. Is 1 for everything except 2D-PE (=2) or 3D (=4?)

        CHANGELOG:
        20130131/RB: cleaned up version from croc. The specific experimental details have been moved to their respective classes. 
        
        """
        
        self.verbose("New dataclass", flag_verbose)
        
        self.obj_id = objectname
        self.sub_type = ""
        
        self.objectname = objectname

        # file stuff
        self.source_path = ""
        self.base_filename = "" 
        self.time_stamp = ""
        self.date = ""

        # organizational stuff
        self.dimensions = dimensions    # t1, t2, t3 etc. so 2D -> 3 dimensions
        self.measurements = measurements        # number of diagrams or measurements etc

        # data
        # b: bins
        # r: response, time domain
        # f: fourier transformed
        # s: spectra     
        # _count is the bin count
        # _axis is the axes   
        
        self.b = [0] * measurements
        self.b_axis = [0] * dimensions
        self.b_count = [0] * measurements
        
        self.r = [0] * measurements       
        self.r_axis = [0] * dimensions    

        self.f = [0] * measurements 
        self.f_axis = [0] * dimensions  

        self.s = [0]                    
        self.s_axis = [0] * dimensions    

        # _units: cm-1, fs etc.
        # _resolution: _axis[i] - _axis[i-1]
        self.r_units = [0] * dimensions
        self.r_resolution = [0] * measurements 

        self.s_units = [0] * dimensions
        self.s_resolution = [0] * dimensions   
        
        # correct spectrometer inaccuracies
        self.r_correction = [0] * dimensions  
        self.r_correction_applied = [0] * dimensions

        # zeropadding
        # _zeropad_to: a specific number of samples
        # _zeropad_by: how many times it should be zeropadded
        self._zeropad_to = None
        self._zeropad_by = 1.0 

        # other experimental stuff
        self._phase_degrees = False         
        self._phase_rad = False
        self.undersampling = False
        self._comment = ""

        self.debug = False

    # comments
    @property
    def comment(self):  
        return self._comment
    @comment.setter
    def comment(self, text):
        self._comment = self._comment + time.strftime("%d/%m/%Y %H:%M:%S: ", time.localtime()) + text + "\n"

    # phase   
    @property
    def phase_degrees(self):
        return self._phase_degrees
    @phase_degrees.setter   
    def phase_degrees(self, phase):
        self._phase_degrees = phase
        self._phase_rad = phase*numpy.pi/180

    @property
    def phase_rad(self):
        return self._phase_rad
    @phase_rad.setter
    def phase_rad(self, phase):
        self._phase_rad = phase
        self._phase_degrees = phase * 180/numpy.pi

    # zeropadding
    @property
    def zeropad_to(self):
        return self._zeropad_to
    @zeropad_to.setter
    def zeropad_to(self, zpt):
        if len(self.s) != 1:
            self.printWarning("Zeropadding has changed after spectrum was calculated.", inspect.stack())
        self._zeropad_to = int(zpt)
        self._zeropad_by = zpt / numpy.shape(self.r[0])[0]

    @property
    def zeropad_by(self):
        return self._zeropad_by
    @zeropad_by.setter
    def zeropad_by(self, zp_by):
        if len(self.s) != 1:
            self.printWarning("Zeropadding has changed after spectrum was calculated.", inspect.stack())
        self._zeropad_to = int(zp_by * numpy.shape(self.r[0])[0])
        self._zeropad_by = zp_by 




if __name__ == "__main__": 
    
    pass

    
    
    
    
    
    
    
    
    
    
    
    

