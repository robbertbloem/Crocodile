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
import Crocodile.Resources.Mathematics as MATH
import Crocodile.Resources.Plotting as PL

imp.reload(DCC)
imp.reload(IOM)

class show_shots(DCC.dataclass):
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
        # check if _file_dict is set
        if self._file_dict["base_folder"] == "" or self._file_dict["base_filename"] == "":
            self.printError("No file information set.", inspect.stack()) 
            return False   
        
        # find LV_file_format, also a check if _file_dict is correct
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

        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        sh = numpy.arange(n_sh)

        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        spds, dump, dump = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        self.sp = spds[:,0]
        self.ds = spds[:,1]
        for i in range(n_ds):
            if self.ds[i] == "-1":
                self.ds[i] = 1
            else:
                self.ds[i] = 0
        self.ds = numpy.array(self.ds, dtype = "int")

        n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
        sc = numpy.arange(n_sc)
         
        self.b_n = [n_w3, n_sh, 2*n_ds, n_sp, n_sc]
        self.b = numpy.empty(self.b_n)
        self.b_axes = [w3_axis, sh, self.ds, self.sp, sc]
        self.b_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "Scans"]
        
#         self.b_n = [n_w3, n_sh, n_ds, n_sp, n_sc]


        self.b_choppers = numpy.zeros([3, n_sh, n_sc])
        self.b_specials = numpy.zeros([15, n_sh, n_sc])

        for sc in range(self.b_n[4]): 
        
            suffix = "specials_" + str(sc)
            self.b_specials[:,:,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T
            
            suffix = "choppers_" + str(sc)
            self.b_choppers[:,:,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)
            
            for ds in range(self.b_n[2]): 
                for sp in range(self.b_n[3]): 
                    # import probe
                    suffix = "sp" + str(sp) + "_ds" + str(ds) + "_pixels_" + str(sc)
                    self.b[:,:,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)  
  

    def make_plot(self):  

        axes = numpy.arange(self.b_n[1])
        
        fig = plt.figure()
        ax = fig.add_subplot(111)    
        
        for ds in range(self.b_n[2]):
            for sc in range(self.b_n[4]):
                data = self.b[15,:,ds,0,sc]
                mask = numpy.isfinite(data)            
                ax.plot(axes[mask], data[mask], marker = ".", linestyle = "none")
            
        plt.show()




if __name__ == "__main__": 
    pass
