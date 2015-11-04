from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect
import os
import imp

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.DataClassCol as DCC
import Crocodile.Resources.IOMethods as IOM
import Crocodile.Resources.Constants as CONST

imp.reload(DCC)
imp.reload(IOM)

class pe_col(DCC.dataclass):
    """
    pe_col: implementation class for colinear 2D-IR setup. This includes all importing functions.
    
    
    
    """

    def __init__(self, objectname, flag_verbose = False):

        self.verbose("New pe_col class", flag_verbose)
    
        DCC.dataclass.__init__(self, objectname = objectname, flag_verbose = flag_verbose)


    def set_file_info(self, data_folder, date, basename, timestamp):
    
        self.set_file_dict(data_folder, date, basename, timestamp)    


    def import_data(self):
        """
        This method splits up file importing for supporting files and measurement files. It first finds the file format file. It will use this to test if file_dict is correct. 
        """
    
        if self._file_dict["base_folder"] == "" or self._file_dict["base_filename"] == "":
            self.printError("No file information set.", inspect.stack()) 
            return False   

        try: 
            self.file_format = IOM.find_LV_fileformat(
                base_folder = self._file_dict["base_folder"], 
                flag_verbose = self.flag_verbose
            )   
            if self.file_format == -1:  
                self.printError("File_format file found, but was not able to parse it.", inspect.stack()) 
                return False  
        except FileNotFoundError:
            self.printError("File_format file not found.", inspect.stack()) 
            return False 
             

   
        self.import_supporting_data()
        
        return True  



    def import_supporting_data(self):
        
        w3_axis_wn, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
  