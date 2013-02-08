from __future__ import print_function
from __future__ import division

import inspect
import time

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile
import PythonTools.ClassTools as CT

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
        # self._zeropad_by = 1.0 # now only a getter/setter method

        # other experimental stuff
        self._phase_degrees = False         
        # self._phase_rad = False
        self.undersampling = False
        self._comment = ""


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
        if type(phase) == bool:
            self.printWarning("phase_degrees is boolean type, will be set to " +  str(phase * 1), inspect.stack())
            phase *= 1
        if numpy.isnan(phase):
            self.printError("phase_degrees can not be numpy.nan. It is not set.", inspect.stack())
        else:
            self._phase_degrees = phase

    @property
    def phase_rad(self):
        return self.phase_degrees * numpy.pi / 180
    @phase_rad.setter
    def phase_rad(self, phase):
        if type(phase) == bool:
            self.printWarning("phase_degrees is boolean type, will be set to " +  str(phase * 1), inspect.stack())
            phase *= 1
        if numpy.isnan(phase):
            self.printError("phase_degrees can not be numpy.nan. It is not set.", inspect.stack())
        else:
            self._phase_degrees = phase * 180 / numpy.pi

    # zeropadding
    @property
    def zeropad_to(self):
        return self._zeropad_to
    @zeropad_to.setter
    def zeropad_to(self, zpt):
        if len(self.s) != 1:
            self.printWarning("Zeropadding has changed after spectrum was calculated.", inspect.stack())        
        if numpy.isnan(zpt):
            self.printError("zeropad_to can not be numpy.nan. It is not set.", inspect.stack())
        elif type(self.r[0]) == int:
            self.printError("zeropad_by can not be determined because r[0] has no length. zeropad_to is not set", inspect.stack())
        else:
            self._zeropad_to = int(zpt)
            
    
    @property
    def zeropad_by(self):
        if self._zeropad_to == None:
            return 1.0
        else:
            return self._zeropad_to / numpy.shape(self.r[0])[0]
    @zeropad_by.setter
    def zeropad_by(self, zp_by):
        if len(self.s) != 1:
            self.printWarning("Zeropadding has changed after spectrum was calculated.", inspect.stack())
        if numpy.isnan(zp_by):
            self.printError("zeropad_by can not be numpy.nan. It is not set.", inspect.stack())
        elif type(self.r[0]) == int:
            self.printError("zeropad_to can not be determined because r[0] has no length. zeropad_to is not set", inspect.stack())
        else:
            self._zeropad_to = int(zp_by * numpy.shape(self.r[0])[0])
            




if __name__ == "__main__": 
    
    pass

    
    
    
    
    
    
    
    
    
    
    
    

