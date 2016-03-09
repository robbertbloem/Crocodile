from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect
import os
import imp

import enum

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
    scan_spectrum = 9
    point_cloud = 10
    2d_ir = 11
    pp = 12
    vcd = 13
    
    

class MosquitoHelperMethods(DCC.dataclass):
    
    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New Mosquito class", flag_verbose)
        DCC.dataclass.__init__(self, objectname = objectname, flag_verbose = flag_verbose)  
        
        
    def import_data(self, measurement_method, import_temp_scans = False):
        """
        This method splits up file importing for supporting files and measurement files. It first finds the file format file. It will use this to test if file_dict is correct. 
        """
        
        # all methods
        # find the file format
        # quit on error
        if find_file_format(self) == False:
            return False

        # all methods
        # for scan spectrum this is the array with all the wavelengths -- it is not necessarily the same as the number of pixels
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        # shots or bins
        if measurement_method in [
            MeasurementMethod["find_t0_fast"], 
            MeasurementMethod["2d_ir"],
        ]:
            # fast scan
            t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
            self.bin_sign = bin_sign
            self.t1_zero_index = t1_zero_index

        else:
            # non-fast scan
            n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
            sh = numpy.arange(n_sh)

        if measurement_method in [
            MeasurementMethod["pp"], 
            MeasurementMethod["2d_ir"],
        ]:
            # datastates, as written to files
            # 0: signal
            # >0: intensity
            n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
            # spectra
            n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

            # datastates, originally in the measurement
            # if the signal is saved, this is different from n_ds
            spds, n_sp_2, n_ds_2 = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
            self.ds = spds[:,1]
            for i in range(n_ds):
                if self.ds[i] == "-1":
                    self.ds[i] = 1
                else:
                    self.ds[i] = 0
            self.ds = numpy.array(self.ds, dtype = "int")
        
            if n_ds == n_ds_2:
                self.measurement_type = "intensity"
            else:
                self.measurement_type = "signal"

            if self.measurement_type == "signal":
                # data already divided by count
                # 1 data state: the signal
                _n_ds = 1
            else:
                # data not yet divided by count
                # 2*n_datastates
                _n_ds = 2 * n_ds

        else:
            self.n_ds = 1
            self.ds = numpy.array([0])

        # slow modulation
        if measurement_method in [
            MeasurementMethod["pp"], 
            MeasurementMethod["2d_ir"],
        ]:
            sm, sm_names, n_sm = IOM.import_slow_modulation(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        else:
            n_sm = 1
            sm_names = ["none"]
            sm = numpy.zeros((1,1)) 
        
        if measurement_method in [
            MeasurementMethod["pp"], 
            MeasurementMethod["2d_ir"], 
            MeasurementMethod["find_t0"], 
            MeasurementMethod["find_t0_fast"],
        ]:
            de, n_de = IOM.import_delays(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
            
        else:
            de = [0]
            n_de = 1
        
        # dummy dimension
        du = numpy.array([0])
        n_du = 1

        # number of scans
        # 
        if import_temp_scans or measurement_method in [
            MeasurementMethod["show_shots"], 
            MeasurementMethod["show_spectrum"],
        ]:
            n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
        else:
            n_sc = 1
        sc = numpy.arange(n_sc)



#         n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
#         sc = numpy.arange(n_sc)
        
        # show shots
        self.b_n = [n_w3, n_sh, 2*n_ds, n_sp, n_sc]
        self.b = numpy.empty(self.b_n)
        self.b_axes = [w3_axis, sh, self.ds, self.sp, sc]
        self.b_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "Scans"]
        
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


        # 2d
        self.b_n = [n_w3, n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
        self.b_intf_n = [n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]

        self.b = numpy.empty(self.b_n)
        self.b_count = numpy.empty(self.b_n[1:]) 
        self.b_axes = [w3_axis, t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
        self.b_units = ["w3 (cm-1)", "T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.b_intf = numpy.empty(self.b_intf_n)
        self.b_intf_axes = [t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
        self.b_intf_units = ["T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.r_n = [n_w3, n_t1_fs, 1, n_sp, n_sm, n_de, n_du, n_sc]
        self.r = numpy.empty(self.r_n)
        self.r_axes = [w3_axis, t1_fs, [0], spds[:,0], sm, de, du, sc]
        self.r_units = ["w3 (cm-1)", "T1 (fs)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.f_n = [n_w3, -1, 1, n_sp, n_sm, n_de, n_du, n_sc]
        self.f_axes = [w3_axis, [0], [0], spds[:,0], sm, de, du, sc]
        self.f_units = ["w3 (cm-1)", "w1 (cm-1)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.s_n = [n_w3, -1, 1, n_sp, n_sm, n_de, n_du, n_sc]
        self.s_axes = [w3_axis, [0], [0], spds[:,0], sm, de, du, sc]
        self.s_units = ["w3 (cm-1)", "w1 (cm-1)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]


        # probe and reference are saved separately
        if self.measurement_type == "signal":
            _n_ds = 1
        else:
            _n_ds = n_ds
        
        self.b_n = [n_w3, n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
        self.b_intf_n = [n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
        
        return True



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
                                suffix = "probe_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b[:,:,2*ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                # import reference
                                suffix = "reference_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b[:,:,2*ds+1,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                # import count
                                suffix = "count_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b_count[:,ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)

                                # import interferogram
                                suffix = "interferogram_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b_intf[:,ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)

        elif self.measurement_type == "intensity" and self.b_n[7] > 1:
            self.verbose("All scans, intensity", self.flag_verbose)
            for ds in range(self.b_n[2]):    
                for sp in range(self.b_n[3]):
                    for sm in range(self.b_n[4]):
                        for de in range(self.b_n[5]): 
                            for du in range(self.b_n[6]):
                                for sc in range(self.b_n[7]):
                    
                                    # import probe
                                    suffix = "probe_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b[:,:,2*ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                    # import reference
                                    suffix = "reference_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b[:,:,2*ds+1,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                    # import count
                                    suffix = "count_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b_count[:,:,ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                                    # import interferogram
                                    suffix = "interferogram_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b_intf[:,ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

        else:
            self.printError("Failed to import measurement files", inspect.stack())                    


        if self.bin_sign:
            self.b = self.b[:,::-1,:,:,:,:,:,:]
            self.b_intf = self.intf[::-1,:,:,:,:,:,:]
            self.b_count = self.b_count[:,::-1,:,:,:,:,:,:]


    def find_file_format(self):
    
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
            
        return True





if __name__ == "__main__": 
    pass
