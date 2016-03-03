from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Resources.MosquitoHelper as MH

class show_shots(MH.MosquitoHelperMethods):
    """
    show_shots
    """

    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New show_shots class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, flag_verbose = flag_verbose)

# 
#     def set_file_info(self, data_folder, date, basename, timestamp):
#         self.set_file_dict(data_folder, date, basename, timestamp)    

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



class show_spectrum(MH.MosquitoHelperMethods):
    """
    show_spectrum
    """

    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New show_spectrum class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, flag_verbose = flag_verbose)


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

#         n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
#         sh = numpy.arange(n_sh)
        
        n_sh = 1
        sh = numpy.array([0])
        
        n_ds = IOM.find_number_of_datastates(self._file_dict["base_folder"], flag_verbose = self.flag_verbose)
        n_ds = int(n_ds/2)
        
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
        self.b_noise = numpy.empty(self.b_n)
        self.b_axes = [w3_axis, sh, self.ds, self.sp, sc]
        self.b_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "Scans"]
        
        self.s_n = numpy.copy(self.b_n)
        self.s_n[2] = 1
        self.s = numpy.zeros(self.s_n)
        self.s_noise = numpy.zeros(self.s_n)

        sh = 0        
        for sp in range(self.b_n[3]): 
            for sc in range(self.b_n[4]): 
                
                ds = 0
                suffix = "sp%i_ds%i_signal_%i" % (sp, ds, sc)
                self.s[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose) 
                suffix = "sp%i_ds%i_signal_noise_%i" % (sp, ds, sc)
                self.s_noise[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)
              
                for ds in range(self.b_n[2]): 
                    suffix = "sp%i_ds%i_intensity_%i" % (sp, ds, sc)
                    self.b[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose) 
                    suffix = "sp%i_ds%i_intensity_noise_%i" % (sp, ds, sc)
                    self.b_noise[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)



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



class find_t0_fast(MH.MosquitoHelperMethods):
    """
    find_t0_fast
    """

    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New find_t0_fast class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, flag_verbose = flag_verbose)


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

        t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        self.bin_sign = bin_sign
        self.t1_zero_index = t1_zero_index

        n_ds = IOM.find_number_of_datastates(self._file_dict["base_folder"], flag_verbose = self.flag_verbose)
        
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        self.ds = spds[:,1]
        for i in range(n_ds):
            if self.ds[i] == "-1":
                self.ds[i] = 1
            else:
                self.ds[i] = 0
        self.ds = numpy.array(self.ds, dtype = "int")
        
        self.b_intf_n = [n_t1_bins, n_ds, n_sp]
        
        self.b_intf = numpy.empty(self.b_intf_n)
        self.b_intf_axes = [t1_bins, spds[:,1], self.ds]
        self.b_intf_units = ["T1 (bins)", "Datastates", "Spectra"]




   



if __name__ == "__main__": 
    pass