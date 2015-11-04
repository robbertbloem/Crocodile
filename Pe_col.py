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




#         PE.pe.__init__(self, objectname, measurements = 1, flag_verbose = flag_verbose)

#         self.time_stamp = timestamp 
# 
#         if MacOSX:
#             self.basename = base_folder + date + "/" + base_filename + "_" + timestamp + "/" + base_filename + "_" + timestamp
#         else:
#             self.basename = base_folder + date + "\\" + base_filename + "_" + timestamp + "\\" + base_filename + "_" + timestamp
#         
#         
#         
# 
#     
#     
        