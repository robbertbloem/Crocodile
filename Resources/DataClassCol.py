from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

# import imp
import enum

import inspect
import time
import os

import numpy

import Crocodile
import PythonTools.ClassTools as CT

# imp.reload(CT)

class MeasurementMethod(enum.Enum):
    show_shots = 0
    show_spectrum = 1
    scan_spectrum = 2
    find_t0_interference_intensity = 3
    find_t0_interference = 4
    find_t0_crystal = 5
    find_t0_find_phase = 6
    find_t0_interferogram = 7
    find_t0_fast = 8
#     scan_spectrum = 9
    point_cloud = 10
    ft_2d_ir = 11
    pump_probe = 12
    vcd = 13


class dataclass(CT.ClassTools):
    """
    This class stores all the data. 
    Most of the variables are lists which point to ndarrays with the real data.
    """

    def __init__(self, objectname, measurement_method, flag_verbose = 0):
        """
        croc.DataClasses.messdata

        INPUT:
        - objectname: the name
        - dimensions: the number of axes. Linear/2D/3D = 1/2/3. If something is measured in frequency and time (certain pump-probe) the dimension is also 2.
        - measurements: number of diagrams. Is 1 for everything except 2D-PE (=2) or 3D (=4?)
        
        
        REMARKS:
        
        There are 4 stages for the data:
        b: bins, the raw measurement. 
            - The signal may or may not have been calculated
            - This may or may not be divided by the count. 
            - The extra bins (t1<0) are not yet removed.
        r: response. Divided by the count and only t1>=0
        f: fourier transformed
        s: spectra    
        
        As in LabVIEW, there are 7 dimensions:
        0. pi: pixels
        1. sh/bi: shots or bins, for 2D measurements always bins
        2. ds: datastates
        3. sp: spectra
        4. sm: slow modulation
        5. de: delays
        6. du: dummies
        7. sc: scans
        
        Every stage has a some data:
        no suffix: the data, 7D 
        _count is the bin count, 6D (pixels are missing)
        _axes are the axes, 7D
        _intf is the interferogram, 7D
        
        CHANGELOG:
        20130131/RB: cleaned up version from croc. The specific experimental details have been moved to their respective classes. 
        
        """
        
        self.verbose("New dataclass", flag_verbose)
        self.flag_verbose = flag_verbose
        
        self.obj_id = objectname
        self.sub_type = ""
        
        self.objectname = objectname

        # file stuff
        self._file_dict = {
            "data_folder": "",
            "date": "",
            "basename": "",
            "timestamp": "",
            "base_folder": "",
            "base_filename": "",
            "extension": "",
        }

        # this allows us to use getter/setter methods
        self._data_folder = ""
        self._date = "" 
        self._basename = ""
        self._timestamp = ""
        self._extension = ""
        
        self.file_format = -1
        self.measurement_method = measurement_method
        self.measurement_type = ""

        # organizational stuff
        dimensions = 8
        measurements = 1
        
        self.dimensions = dimensions
        self.measurements = measurements    
        self.add_datastates = True
        self.n_sig = 0 # shots to calculate signal

# 
#         self.b = [0] * 8
#         self.b_noise = [0] * 8
#         self.b_axes = [0] * 8
#         self.b_count = [0] * 8
#         self.b_units = [""] * 8
#         self.b_n = [0] * 8
#           
#         self.r = [0] * 8
#         self.r_noise = [0] * 8
#         self.r_axes = [0] * 8  
#         self.r_units = [""] * 8
#         self.f_n = [0] * 8
#         
#         self.f = [0] * 8
#         self.f_noise = [0] * 8
#         self.f_axes = [0] * 8
#         self.f_units = [""] * 8
#         self.f_n = [0] * 8
#         
#         self.s = [0] * 8
#         self.s_noise = [0] * 8
#         self.s_axes = [0] * 8
#         self.s_units = [""] * 8
#         self.s_n = [0] * 8
#         
# 
#         self.b_intf = [0] * 8
#         self.b_intf_axes = [0] * 8
#         self.b_intf_n = [0] * 8
#         self.b_intf_units = [0] * 8
#         
#         self.r_intf = [0] * 8
#         self.r_intf_axes = [0] * 8
#         self.r_intf_n = [0] * 8
#         self.r_intf_units = [0] * 8
#         
#         self.f_intf = [0] * 8
#         self.f_intf_axes = [0] * 8
#         self.f_intf_n = [0] * 8
#         self.f_intf_units = [0] * 8

        self.n_scans = 0

        # zeropadding
        # _zeropad_to: a specific number of samples
        self._zeropad_to = None
        self._phase_degrees = None         
        self.undersampling = False

        self._comment = ""

    # file stuff
    @property
    def file_dict(self):
        return self._file_dict
    @file_dict.setter
    def file_dict(self, dict):
        self._file_dict = dict
        
    @property
    def data_folder(self):
        return self._data_folder
    @data_folder.setter
    def data_folder(self, text):
        if os.name == "posix":
            if text[-1] != "/":
                text += "/"
        elif os.name == "nt":
            if text[-2:] != "\\\\":
                if text[-1] != "\\":
                    text += r"\\\\"  
                else:
                    text += "\\"  
        else:
            print("unknown operating system")
            
        self._data_folder = text
        self._file_dict["data_folder"] = text
        self.construct_file_paths()

    @property
    def date(self):
        return self._date   
    @date.setter
    def date(self, text):
        self._date  = text
        self._file_dict["date"] = text
        self.construct_file_paths()
        
    @property
    def basename(self):
        return self._basename
    @basename.setter
    def basename(self, text):
        self._basename = text
        self._file_dict["basename"] = text
        self.construct_file_paths()

    @property
    def timestamp(self):
        return self._timestamp
    @timestamp.setter
    def timestamp(self, text):
        self._timestamp = text
        self._file_dict["timestamp"] = text
        self.construct_file_paths()

    @property
    def extension(self):
        return self._extension   
    @date.setter
    def extension(self, text):
        if text[0] != ".":
            text = "." + text
        self._extension = text
        self._file_dict["extension"] = text
        self.construct_file_paths()

    def set_file_dict(self, data_folder, date, basename, timestamp, extension = "csv"):
    
        self.data_folder = data_folder
        self.date = date
        self.basename = basename
        self.timestamp = timestamp
        self.extension = extension

        

    def construct_file_paths(self):

        if os.name == "posix":
            sep = "/"
        elif os.name == "nt":
            sep = r"\\"
        else:
            print("unknown operating system")
        
        self._file_dict["base_folder"] = self._file_dict["data_folder"] + self._file_dict["date"] + sep + self._file_dict["basename"] + "_" + self._file_dict["timestamp"] + sep
        
        self._file_dict["base_filename"] = self._file_dict["base_folder"] + self._file_dict["basename"] + "_" + self._file_dict["timestamp"]
      


    # comments
    @property
    def comment(self):  
        return self._comment
    @comment.setter
    def comment(self, text):
        self._comment = self._comment + time.strftime("%d/%m/%Y %H:%M:%S: ", time.localtime()) + text + "\n"

# 
#     # time_stamp
#     @property
#     def time_stamp(self):
#         return self._time_stamp
#     
#     @time_stamp.setter   
#     def time_stamp(self, ts):
#         self._time_stamp = str(ts)
#         if len(self._time_stamp) == 3:
#             self._time_stamp = "0" + self._time_stamp
#         try:
#             res = int(self._time_stamp)
#             if len(self._time_stamp) == 4 and res >= 2400:
#                 self.printWarning("time_stamp >= 2400. That's not an acceptable time.", inspect.stack())
#             elif len(self._time_stamp) == 6 and res >= 240000:
#                 self.printWarning("time_stamp >= 2400. That's not an acceptable time.", inspect.stack())
#         except ValueError:
#             self.printWarning("time_stamp is not numerical. Please use 0000 or something.", inspect.stack())
#         

    # phase   
    @property
    def phase_degrees(self):
        return self._phase_degrees
    
    @phase_degrees.setter   
    def phase_degrees(self, phase):
        if type(phase) == bool:
            self.printWarning("phase_degrees is boolean type, will be set to " +  str(phase * 1), inspect.stack())
            phase *= 1
        if phase == None:
            # None and numpy.isnan don't work together
            self.printWarning("phase_degrees is set to None", inspect.stack()) 
            self._phase_degrees = None 
        elif numpy.isnan(phase):
            self.printError("phase_degrees can not be numpy.nan. Will be set to None", inspect.stack())
            self._phase_degrees = None 
        else:
            self._phase_degrees = phase

    @property
    def phase_rad(self):
        if self._phase_degrees == None:
            return None
        else:
            return self.phase_degrees * numpy.pi / 180
    
    @phase_rad.setter
    def phase_rad(self, phase):
        if type(phase) == bool:
            self.printWarning("phase_rad is boolean type, will be set to " +  str(phase * 1), inspect.stack())
            phase *= 1
        if phase == None:
            # None and numpy.isnan don't work together
            self.printWarning("phase_rad is set to 1.0", inspect.stack()) 
            self._phase_degrees = None 
        elif numpy.isnan(phase):
            self.printError("phase_rad can not be numpy.nan. Will be set to None", inspect.stack())
            self._phase_degrees = None 
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

    
    
    
    
    
    
    
    
    
    
    
    

