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


    
    

class MosquitoHelperMethods(DCC.dataclass):
    
    def __init__(self, objectname, measurement_method, flag_verbose = False):
        self.verbose("New Mosquito class", flag_verbose)
                
        DCC.dataclass.__init__(self, objectname = objectname, measurement_method = measurement_method, flag_verbose = flag_verbose)  



    def assign_axes(self, dim, axes, units, size):

        self.b_axes[dim] = self.r_axes[dim] = self.f_axes[dim] = self.s_axes[dim] = axis
        self.b_units[dim] = self.r_units[dim] = self.f_units[dim] = self.s_units[dim] = units
        self.b_n[dim] = self.r_n[dim] = self.f_n[dim] = self.s_n[dim] = size
        
#     def assign_b_axes(self, dim, axes, units, size):
# 
#         self.b_axes[dim] = axis
#         self.b_units[dim] = units
#         self.b_n[dim] = size
# 
#     def assign_r_axes(self, dim, axes, units, size):
# 
#         self.r_axes[dim] = axis
#         self.r_units[dim] = units
#         self.r_n[dim] = size
#    
#     def assign_f_axes(self, dim, axes, units, size):
# 
#         self.f_axes[dim] = axis
#         self.f_units[dim] = units
#         self.f_n[dim] = size
#         
#     def assign_s_axes(self, dim, axes, units, size):
# 
#         self.s_axes[dim] = axis
#         self.s_units[dim] = units
#         self.s_n[dim] = size
# 
#     def import_data(self, import_temp_scans = False, ssh_reload_data = False):
#         """
#         This method splits up file importing for supporting files and measurement files. It first finds the file format file. It will use this to test if file_dict is correct. 
#         
#         INPUT:
#         import_temp_scans (Bool, False): 
#         ssh_reload_data (Bool, False): show shots. Data imported from csv will be saved as a numpy array. If False, if the numpy array file exists, this will be used. Otherwise, or when True, the data will be imported from csv (again). 
#         """
#         
#         MM = DCC.MeasurementMethod
#         
#         # all methods
#         # find the file format
#         # quit on error
#         if self.find_file_format() == False:
#             return False
#         
#         # all methods
#         # for scan spectrum this is the array with all the wavelengths -- it is not necessarily the same as the number of pixels
#         w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)     
#         
#         dim = 0
#         self.assign_axes(dim, w3_axes, units = "w3 (cm-1)", size = n_w3)
# 
# 
#         # number of scans
#         if import_temp_scans or self.measurement_method in [
#             MM["show_shots"], 
#             MM["show_spectrum"],
#         ]:
#             n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
#         else:
#             n_sc = 1
#         dim = 7
#         self.assign_axes(dim, numpy.arange(n_sc), units = "scans", size = n_sc)
# 
# 
# #         n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
# #         sc = numpy.arange(n_sc)        
#         
#         # shots or bins
#         if self.measurement_method in [
#             MM["find_t0_fast"], 
#             MM["ft_2d_ir"],
#         ]:
#             # fast scan
#             t1_bins, t1_fs, self.bin_sign, n_t1_bins, n_t1_fs, self.t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
# 
#             dim = 1
#             self.assign_b_axes(dim, t1_bins, units = "t1 (bins)", size = n_t1_bins)    
#             self.assign_r_axes(dim, t1_fs, units = "t1 (fs)", size = n_t1_fs)    
#             
# 
#         elif self.measurement_method in [
#             MM["show_shots"],
#         ]:
# 
#             n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
#             sh = numpy.arange(n_sh) 
#             
#             dim = 1
#             self.assign_r_axes(dim, numpy.arange(n_sc * n_sh), units = "shots", size = n_sc * n_sh)  
#             
#             n_sc = 1
#             dim = 7
#             self.assign_axes(dim, numpy.arange(n_sc), units = "scans", size = n_sc)
#             
#         else:
#             n_sh = 1
#             sh = numpy.array([0])
        
        
    def import_data(self, import_temp_scans = False, ssh_reload_data = False):
        """
        This method splits up file importing for supporting files and measurement files. It first finds the file format file. It will use this to test if file_dict is correct. 
        
        INPUT:
        import_temp_scans (Bool, False): 
        ssh_reload_data (Bool, False): show shots. Data imported from csv will be saved as a numpy array. If False, if the numpy array file exists, this will be used. Otherwise, or when True, the data will be imported from csv (again). 
        """
        
        MM = DCC.MeasurementMethod
        
        # all methods
        # find the file format
        # quit on error
        if self.find_file_format() == False:
            return False
        
        # all methods
        # for scan spectrum this is the array with all the wavelengths -- it is not necessarily the same as the number of pixels
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        
        # shots or bins
        if self.measurement_method in [
            MM["find_t0_fast"], 
            MM["ft_2d_ir"],
        ]:
            # fast scan
            t1_bins, t1_fs, self.bin_sign, n_t1_bins, n_t1_fs, self.t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        elif self.measurement_method in [
            MM["show_shots"],
        ]:
            # non-fast scan
            n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
#             sh = numpy.arange(n_sh)        
        else:
            n_sh = 1
#             sh = numpy.array([0])

        # datastates, as written to files
        # 0: signal
        # >0: intensity
        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        if n_ds == 0:
            self.measurement_type = "signal"
        else:
            self.measurement_type = "intensity"
        ds = numpy.arange(n_ds)
    
        # spectra
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        sp = numpy.arange(n_sp)

        # the spds file is about how the signal is calculated
        # self.n_sig and self.add_ds were introduced in LV_fileformat 4. For pp and 2D-IR this is handled in the importer method
        # for show shots/spectrum etc the spds file was saved from LV_fileformat 4 onward. 
        
        if (self.measurement_method in [
            MM["pump_probe"], 
            MM["ft_2d_ir"],
        ]) or (
            self.measurement_method in [
                MM["show_shots"],
                MM["show_spectrum"],
                MM["scan_spectrum"],
            ] and self.file_format >= 4
        ):


            # datastates, originally in the measurement
            # if the signal is saved, this is different from n_ds
            spds, n_sp_2, n_ds_2, self.n_sig, self.add_ds = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
            ds = spds[:,1]
            for i in range(n_ds):
                if ds[i] == "-1":
                    ds[i] = 1
                else:
                    ds[i] = 0
            ds = numpy.array(ds, dtype = "int")
        
            if n_ds == n_ds_2:
                self.measurement_type = "intensity"
            else:
                self.measurement_type = "signal"

#             if self.measurement_type == "signal":
#                 # data already divided by count
#                 # 1 data state: the signal
#                 _n_ds = 1
#             else:
#                 # data not yet divided by count
#                 # 2*n_datastates
#                 _n_ds = 2 * n_ds

        elif self.measurement_method in [
                MM["show_shots"],
                MM["show_spectrum"],
                MM["scan_spectrum"],
            ] and self.file_format < 4:
            
            # for compatibility with some measurements with the VCD setup
            self.n_sig = 6
            self.add_ds = False
            

        # slow modulation
        if self.measurement_method in [
            MM["pump_probe"], 
            MM["ft_2d_ir"],
        ]:
            sm, sm_names, n_sm = IOM.import_slow_modulation(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        else:
            n_sm = 1
            sm_names = ["none"]
            sm = numpy.zeros((1,1)) 
        
        # delays
        if self.measurement_method in [
            MM["pump_probe"], 
            MM["ft_2d_ir"], 
            MM["find_t0_fast"],
        ]:
            de, n_de = IOM.import_delays(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
            
        else:
            de = numpy.array([0])
            n_de = 1
        
        # dummy dimension
        du = numpy.array([0])
        n_du = 1

        # number of scans
        if import_temp_scans or self.measurement_method in [
            MM["show_shots"], 
            MM["show_spectrum"],
        ]:
            n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
        else:
            n_sc = 1
        sc = numpy.arange(n_sc)



#         n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
#         sc = numpy.arange(n_sc)
        
        #### IMPORT THE DATA ####
        
        if self.measurement_method in [
            MM["show_shots"],
        ]:
#             n_sc = 1

            total_shots = n_sh * n_sc
            n_signals = int(n_sc * n_sh / self.n_sig)
            
            self.r_n = [n_w3, total_shots, 2*n_ds, n_sp, n_sm, n_de, n_du, 1]
            self.r = numpy.empty(self.r_n)
            self.r_axes = [w3_axis, numpy.arange(total_shots), ds, sp, sm, de, du, numpy.array([0])]
            self.r_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "x", "x", "x", "Scans"]
            
            self.f_n = [n_w3, n_signals, 2*n_ds, n_sp, n_sm, n_de, n_du, 1]
            self.f = numpy.empty(self.r_n)
            self.f_axes = [w3_axis, numpy.arange(n_signals), ds, sp, sm, de, du, numpy.array([0])]
            self.f_units = ["w3 (cm-1)", "Signals", "Datastates", "Spectra", "x", "x", "x", "Scans"]
            
            self.s_n = [n_w3, n_signals, 1, n_sp, n_sm, n_de, n_du, 1]
            self.s = numpy.empty(self.r_n)
            self.s_axes = [w3_axis, numpy.arange(n_signals), ds, sp, sm, de, du, numpy.array([0])]
            self.s_units = ["w3 (cm-1)", "Signals", "x", "Spectra", "x", "x", "x", "Scans"]
            
            paf = self.file_dict["base_filename"] + "_f.npy"
            
            if os.path.isfile(paf) and ssh_reload_data == False: 

                self.f = numpy.load(paf, mmap_mode = "r")

                paf = self.file_dict["base_filename"] + "_r_choppers.npy"
                self.r_choppers = numpy.load(paf, mmap_mode = "r")

                paf = self.file_dict["base_filename"] + "_r_specials.npy"
                self.r_specials = numpy.load(paf, mmap_mode = "r")    
        
                paf = self.file_dict["base_filename"] + "_r.npy"
                self.r = numpy.load(paf, mmap_mode = "r")    

            else:

                for _sc in range(n_sc): 
        
                    s = _sc * n_sh
                    e = (_sc+1) * n_sh
        
                    suffix = "specials_" + str(_sc)
                    temp_s = IOM.import_file(self._file_dict, suffix, self.flag_verbose).T

                    suffix = "choppers_" + str(_sc)
                    temp_c = IOM.import_file(self._file_dict, suffix, self.flag_verbose)
 
                    if sc == 0:
                        n_specials, dump = numpy.shape(temp_s)
                        n_choppers, dump = numpy.shape(temp_c)
                        self.r_choppers = numpy.zeros([n_choppers, n_sh * n_sc])
                        self.r_specials = numpy.zeros([n_specials, n_sh * n_sc])                  
                    
                    self.r_specials[:,s:e] = temp_s
                    self.r_choppers[:,s:e] = temp_c
                
                    for _ds in range(self.r_n[2]): 
                        for _sp in range(self.r_n[3]): 
                            # import probe
                            suffix = "sp" + str(_sp) + "_ds" + str(_ds) + "_data_" + str(_sc)
                            self.b[:, s:e, _ds, _sp, 0,0,0,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)  


                x = 1                
                for _ch in range(n_choppers):        
                    temp = numpy.amax(self.r_choppers[_ch,:])
                    if temp > 0.5:
                        self.r_choppers[_ch,:] /= temp
                        self.r_choppers[_ch,:] *= x
                    x *= 2
                
                ds_list = numpy.sum(self.r_choppers, axis = 0).astype(int)
                
                n_signals = int(n_sc * n_sh / self.n_sig)
                
                ds_count = numpy.zeros(n_ds)
                for _ds in range(n_ds):
                    ds_count[_ds] = len(numpy.where(ds_list == _ds)[0]) / n_signals
                
                
                self.r = numpy.zeros((n_w3, n_signals, n_ds)) 
                count = numpy.zeros((n_signals, n_ds)) 
                
                for sig in range(n_signals):
                    s = sig * self.n_sig
                    e = (sig+1) * self.n_sig
    
                    for _sh in range(self.n_sig):        
                        self.r[:, sig, ds_list[s+_sh]] += (self.b[:, s+_sh, 2*ds_list[s+_sh], 0, 0, 0, 0, 0] / (self.b[:, s+_sh, 2*ds_list[s+_sh]+1, 0, 0, 0, 0, 0] * ds_count[ds_list[s+_sh]]))
                        
                paf = "{base_filename}_r.npy".format(base_filename = self.file_dict["base_filename"])
                numpy.save(paf, self.r)

                paf = "{base_filename}_f.npy".format(base_filename = self.file_dict["base_filename"])
                numpy.save(paf, self.f)
                
                paf = "{base_filename}_r_choppers.npy".format(base_filename = self.file_dict["base_filename"])
                numpy.save(paf, self.r_choppers)

                paf = "{base_filename}_r_specials.npy".format(base_filename = self.file_dict["base_filename"])
                numpy.save(paf, self.r_specials)
                


        elif self.measurement_method in [
            MM["show_spectrum"],
        ]:
            
            self.r_n = [n_w3, 1, 2*n_ds, n_sp, n_sm, n_de, n_du, n_sc]
            self.r = numpy.empty(self.r_n)
            self.r_axes = [w3_axis, numpy.arange(total_shots), ds, sp, sm, de, du, numpy.array([0])]
            self.r_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "x", "x", "x", "Scans"]
            
            self.s_n = [n_w3, 1, 1, n_sp, n_sm, n_de, n_du, n_sc]
            self.s = numpy.empty(self.r_n)
            self.s_axes = [w3_axis, numpy.arange(n_signals), ds, sp, sm, de, du, numpy.array([0])]
            self.s_units = ["w3 (cm-1)", "Signals", "x", "Spectra", "x", "x", "x", "Scans"]            
            
            
            
            for _ds in range(self.r_n[2]): 
                for _sp in range(self.r_n[3]): 
                    if self.r_n[7] == 0:
            
                        suffix = "sp{sp}_ds{ds}_intensity".format(sp = _sp, ds = _ds)
                        self.r[:, 0, _ds, _sp, 0,0,0,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)    
                        
                    else:
                        for _sc in range(self.r_n[7]): 
                            suffix = "sp{sp}_ds{ds}_intensity_{sc}".format(sp = _sp, ds = _ds, sc = _sc)                
                            self.r[:, 0, _ds, _sp, 0,0,0, _sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)    
             


        elif self.measurement_method in [
            MM["ft_2d_ir"],
        ]:

            # probe and reference are saved separately
            if self.measurement_type == "signal":
                _n_ds = 1
                
                self.b_n = [n_w3, n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
                self.b_intf_n = [n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]

                self.b = numpy.empty(self.b_n)
                self.b_count = numpy.empty(self.b_n[1:]) 
                self.b_axes = [w3_axis, t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
                self.b_units = ["w3 (cm-1)", "T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

                self.b_intf = numpy.empty(self.b_intf_n)
                self.b_intf_axes = [t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
                self.b_intf_units = ["T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

            else:
                _n_ds = 2 * n_ds

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


        elif self.measurement_method in [
            MM["pump_probe"],
        ]:

            # probe and reference are saved separately
            if self.measurement_type == "signal":
                _n_ds = 1
                _ds = numpy.array([0])
                
                self.r_n = [n_w3, n_sh, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
                self.r = numpy.empty(self.r_n)
                self.r_axes = [w3_axis, sh, _ds, spds[:,0], sm, de, du, sc]
                self.r_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

            else:
                _n_ds = n_ds
                _ds = ds

                self.s_n = [n_w3, n_sh, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
                self.s_axes = [w3_axis, sh, _ds, spds[:,0], sm, de, du, sc]
                self.s_units = ["w3 (cm-1)", "w1 (cm-1)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]


        elif self.measurement_method in [
            MM["vcd"],
        ]:

            # probe and reference are saved separately
            if self.measurement_type == "signal":
                _n_ds = 1
                _ds = numpy.array([0])
            else:
                _n_ds = n_ds
                _ds = ds

            self.r_n = [n_w3, n_sh, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
            self.r = numpy.empty(self.r_n)
            self.r_axes = [w3_axis, sh, _ds, spds[:,0], sm, de, du, sc]
            self.r_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

            self.s_n = [n_w3, n_sh, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
            self.s_axes = [w3_axis, sh, _ds, spds[:,0], sm, de, du, sc]
            self.s_units = ["w3 (cm-1)", "w1 (cm-1)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

# 

#         
#             self.b_n = [n_w3, n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]
#             self.b_intf_n = [n_t1_bins, _n_ds, n_sp, n_sm, n_de, n_du, n_sc]







  
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


    


    def b_to_r(self):
        """
        Calculate the signal from the probe and reference intensity.
        """ 
        print()
        print("hi")
        
        if self.measurement_type == "intensity":
            
            b_size = numpy.array(numpy.shape(self.b))
            b_size[2] = 1
            b_temp = numpy.zeros(b_size)
            b_intf_temp = numpy.zeros(b_size[1:])

            for bi in range(self.b_n[1]):  
                for sp in range(self.b_n[3]):
                    for sm in range(self.b_n[4]):
                        for de in range(self.b_n[5]): 
                            for du in range(self.b_n[6]):
                                for sc in range(self.b_n[7]):
                                    
                                    temp = numpy.zeros((self.b_n[0], 2))
                                    temp_intf = 0
                                    count = numpy.zeros((2))
                                
                                    for ds in range(self.b_n[2]):  
                                    
                                        if self.b_count[bi, ds, sp, sm, de, du, sc] != 0:
                                            temp[:,self.ds[ds]] = self.b[:, bi, 2*ds, sp, sm, de, du, sc] / self.b[:, bi, 2*ds+1, sp, sm, de, du, sc] 
                                            
                                            temp_intf = self.b_intf[bi, ds, sp, sm, de, du, sc] / self.b_count[bi, ds, sp, sm, de, du, sc]
                                            
                                            count[self.ds[ds]] += 1
            
            
                                b_intf_temp[bi, 0, sp, sm, de, du, sc]  = temp_intf / (count[0] + count[1])
            
                                if count[0] > 0 and count[1] > 0:
                                    b_temp[:, bi, 0, sp, sm, de, du, sc] = (temp[:, 0] / count[0]) / (temp[:, 1] / count[1])     
                                elif count[0] == 0 and count[1] > 0:
                                    b_temp[:, bi, 0, sp, sm, de, du, sc] = 1 / (temp[:, 1] / count[1]) 
                                elif count[0] > 0 and count[1] == 0:
                                    b_temp[:, bi, 0, sp, sm, de, du, sc] = (temp[:, 0] / count[0])  
                                else:  
                                    pass                                                

            self.r = b_temp[:,self.t1_zero_index:,:,:,:,:,:,:]
            self.r_intf = b_intf_temp[:,:,:,:,:,:,:]
            

        else:
            self.r = self.b[:,self.t1_zero_index:,:,:,:,:,:,:]
            self.r_intf = self.b_intf[:,:,:,:,:,:,:]
        








    def calculate_phase(self, n_points = 5, w_range = [0,-1]):
        """
        This methods does the FFT of the interferometer. It looks for the peak in the real part. It then uses n_points, centered around the peak, to calculate the phase.
        
        
        INPUT:
        - n_points (int, default 5): number of points used to calculate the phase. If even, it looks for the highest neigbor of the peak and centers around that. 
        - range (list with length 2, default: [0,-1]). The wavenumber range to look for the phase. 
        
        EXAMPLES:
        For a real part with values: [3,1,2,5,4,3]
        n_points = 2: [5,4]
        n_points = 3: [2,5,4]
        n_points = 4: [2,5,4,3]
        
        """

        N_bins = self.b_n[1]
        N_bins_half = int(N_bins/2)
        dt = self.r_axes[1][1] - self.r_axes[1][0]

        # not the same as w1 axis, that one is truncated to exclude t1 < 0
        w_axis = MATH.make_ft_axis(N_bins, dt = dt, undersampling = 0, normalized_to_period = 0, zero_in_middle = False, flag_verbose = False)
        w_axis = w_axis[:N_bins_half]
        
        i_range = [0,-1]
        if w_range != [0,-1]:
            i_range[0] = numpy.where(w_axis > w_range[0])[0][0]
            i_range[1] = numpy.where(w_axis > w_range[1])[0][0]      
            self.verbose("calculate_phase: peak searching between indices %i and %i (%.1f and %.1f cm-1)" % (i_range[0], i_range[1], w_axis[i_range[0]], w_axis[i_range[1]]), self.flag_verbose)

        r_intf = numpy.zeros(self.r_n[1])
        for bi in range(self.r_n[1]): 
            r_intf[bi] = numpy.mean(self.r_intf[bi,:,:,:,:,:,:])
        r_intf = numpy.roll(r_intf, -self.t1_zero_index)
        f = numpy.fft.fft(r_intf)
        f = f[:N_bins_half] 
        
        self.find_phase(f, w_axis, n_points, i_range)   
        


    def find_phase(self, f, w_axis, n_points, i_range = [0,-1]):
        """
        This function finds the phase. It uses two tricks.
        
        The user can set the number of points used to calculate the phase. 
        
        1. Finding peaks
        I find the index of the maximum of the real values. If the measurement has a lot of problems there may be false peaks. 
        
        - For an odd n_points, I simply take the values to the left and right. So for n_points = 5 and a maximum at index 100, I use indices [98,99,100,101,102]. 
        - For an even n_points, I look at the value of the neighboring datapoints. I place the center between the maximum and the highest of the two neighbors. So for y=[3,1,2,5,4,3] and n_points = 2, I use indices [3,4].
        
        Comparison with LabView: LabView can find peaks in between indices, at 'index' 111.3 for example (it doesn't consider the x-axis indices).
        
        2. Handling angle around 180 degrees. 
        The angle in numpy is limited from -180 to 180 degrees. That is a problem when the angles are [-175, -179, 181]. I first check if all individual elements of the selection are between -135 and +135 degrees. If so, there shouldn't be a problem. If not, I phase shift the fourier-array, calculate  the angles again and then subtract the phase shift again. 
        
        Comparison with LabView: LabView can unpack the angle so that it is not limited to +/- 180 degrees. 
        
        """
        
        y = numpy.real(f)
        i_max = numpy.argmax(y[i_range[0]:i_range[1]])
        idx = numpy.arange(i_max, i_max + n_points)
        
        i_max += i_range[0]
        
        if n_points % 2 == 0:        
            if y[i_max-1] > y[i_max+1]:
                idx -= (n_points/2)
            else:
                idx -= (n_points/2 - 1)           
        else:
            idx -= int(n_points/2)
        
        idx += i_range[0]

        angle = numpy.angle(f)
        check_angle = numpy.angle(f * numpy.exp(1j*numpy.pi/2))

        if numpy.all(angle[idx] > -3*numpy.pi/4) and numpy.all(angle[idx] < 3*numpy.pi/4):
            self.phase_rad = numpy.mean(angle[idx])
        else:
            self.phase_rad = numpy.mean(check_angle[idx]) - numpy.pi/2

        print("Phase in degrees: %.1f" % (self.phase_rad * 180 / numpy.pi))




    def make_fft(self, phase_cheat_deg = 0):

        phase_rad = self.phase_rad + phase_cheat_deg * numpy.pi / 180

        if phase_cheat_deg != 0:
            self.printWarning("Phase cheating with %.1f degrees. Used phase is %.1f" % (phase_cheat_deg, phase_rad * 180 / numpy.pi))

        N_FFT_bins = self.r_n[1]
        N_FFT_half = int(N_FFT_bins / 2)
        N_FFT_bins = 2 * N_FFT_half

        w1_axis = MATH.make_ft_axis(N_FFT_bins, dt = self.r_axes[1][1] - self.r_axes[1][0], undersampling = 0, normalized_to_period = 0, zero_in_middle = False, flag_verbose = self.flag_verbose)
        self.f_axes[1] = w1_axis
        
        w1_axis = w1_axis[:N_FFT_half]
        self.s_axes[1] = w1_axis

        self.f_n[1] = N_FFT_bins
        self.s_n[1] = N_FFT_half
        
        self.f = numpy.zeros(self.f_n, dtype = "complex")
        self.s = numpy.zeros(self.s_n)

        ds = 0
        for pi in range(self.r_n[0]):  
            for sp in range(self.r_n[3]):
                for sm in range(self.r_n[4]):
                    for de in range(self.r_n[5]): 
                        for du in range(self.r_n[6]):
                            for sc in range(self.r_n[7]):
                                r = self.r[pi, :, ds, sp, sm, de, du, sc]
                                f = numpy.fft.fft(r)                      
                                self.f[pi, :, ds, sp, sm, de, du, sc] = f
                                s = numpy.real(numpy.exp(-1j * phase_rad) * f[:N_FFT_half])
                                self.s[pi, :, ds, sp, sm, de, du, sc] = s
                     




    def check_axis(self, ax):
        
        new_axis = False
        if ax == False:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            new_axis = True
            
        return ax, new_axis


if __name__ == "__main__": 
    pass
