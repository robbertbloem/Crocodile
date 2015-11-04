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


    def import_data(self, import_temp_scans = False):
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
             

   
        self.import_supporting_data(import_temp_scans)
        self.import_measurement_data()
        
        return True  



    def import_supporting_data(self, import_temp_scans = False):
        
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        self.bin_sign = bin_sign
        
        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        spds, n_sp_2, n_ds_2 = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        if n_ds == n_ds_2:
            self.measurement_type = "intensity"
        else:
            self.measurement_type = "signal"
        
        sm, sm_names, n_sm = IOM.import_slow_modulation(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        de, n_de = IOM.import_delays(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        du = numpy.array([0])
        n_du = 1
        
        if import_temp_scans:
            n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
        else:
            n_sc = 1
        sc = numpy.arange(n_sc)

        # bin axis
        # signal may not have been calculated
        # maybe not yet divided by count
        if self.measurement_type == "signal":
            _n_ds = 1
        else:
            _n_ds = 2 * n_ds
        
        self.b_n = [n_w3, n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
        self.b_intf_n = [n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]

        self.b = numpy.empty(self.b_n)
        self.b_count = numpy.empty(self.b_n[1:]) 
        self.b_axes = [w3_axis, t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
        self.b_units = ["Wavenumbers (cm-1)", "T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.b_intf = numpy.empty(self.b_intf_n)
        self.b_intf_axes = [t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
        self.b_intf_units = ["T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.r_n = [n_w3, n_t1_fs, 1, n_sp, n_sm, n_de, n_du, n_sc]
        self.r = numpy.empty(self.r_n)
        self.r_axes = [w3_axis, t1_fs, [0], spds[:,0], sm, de, du, sc]
        self.r_units = ["Wavenumbers (cm-1)", "T1 (fs)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        # probe and reference are saved separately
        if self.measurement_type == "signal":
            _n_ds = 1
        else:
            _n_ds = n_ds
        
        self.b_n = [n_w3, n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
        self.b_intf_n = [n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]


#         self.r_intf = numpy.empty(self.b_intf_n)
#         r_t1_bins = 
#         self.r_intf_axes = [t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
#         self.r_intf_units = ["Wavenumbers (cm-1)", "T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]        

        
    def import_measurement_data(self):
        
        if self.measurement_type == "signal" and self.b_n[7] == 1:  
            self.verbose("Average scan, signal", self.flag_verbose)
                
            for sp in range(self.b_n[3]):
                for sm in range(self.b_n[4]):
                    for de in range(self.b_n[5]): 
                        for du in range(self.b_n[6]):
                    
                            # import files
                            suffix = "signal_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                            self.b[:,:,0,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                            # import interferogram
                            suffix = "interferogram_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                            self.b_intf[:,0,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T
                            

        elif self.measurement_type == "signal" and self.b_n[7] > 1:
            self.verbose("All scans, signal", self.flag_verbose)
            for sp in range(self.b_n[3]):
                for sm in range(self.b_n[4]):
                    for de in range(self.b_n[5]): 
                        for du in range(self.b_n[6]):
                            for sc in range(self.b_n[7]):
                    
                                # import files
                                suffix = "signal_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du) + "_" + str(sc)
                                self.b[:,:,0,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T


                                # import interferogram
                                suffix = "interferogram_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du) + "_" + str(sc)
                                self.b_intf[:,0,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T                           

        elif self.measurement_type == "intensity" and self.b_n[7] == 1:
            self.verbose("Average scan, intensity", self.flag_verbose)
            for ds in range(self.b_n[2]):    
                for sp in range(self.b_n[3]):
                    for sm in range(self.b_n[4]):
                        for de in range(self.b_n[5]): 
                            for du in range(self.b_n[6]):
                    
                                # import probe
                                suffix = "probe_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b[:,:,ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                # import reference
                                suffix = "ref_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b[:,:,2*ds+1,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                # import count
                                suffix = "count_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b_count[:,:,ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                # import interferogram
                                suffix = "interferogram_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b_intf[:,ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

        elif self.measurement_type == "intensity" and self.b_n[7] > 1:
            self.verbose("All scans, intensity", self.flag_verbose)
            for ds in range(self.b_n[2]):    
                for sp in range(self.b_n[3]):
                    for sm in range(self.b_n[4]):
                        for de in range(self.b_n[5]): 
                            for du in range(self.b_n[6]):
                                for sc in range(self.b_n[7]):
                    
                                    # import probe
                                    suffix = "probe_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b[:,:,ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                    # import reference
                                    suffix = "ref_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b[:,:,2*ds+1,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                    # import count
                                    suffix = "count_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b_count[:,:,ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                    # import interferogram
                                    suffix = "interferogram_ds" + str(ds) + "sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b_intf[:,ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

        else:
            self.printError("Failed to import measurement files", inspect.stack())                    


        if self.bin_sign:
            self.b = self.b[:,::-1,:,:,:,:,:,:]
            self.b_intf = self.intf[::-1,:,:,:,:,:,:]
            self.b_count = self.b_count[:,::-1,:,:,:,:,:,:]





















