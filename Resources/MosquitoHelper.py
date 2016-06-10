from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect
import os
import imp

# import enum

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.Debug as DEBUG

import Crocodile.Resources.DataClassCol as DCC
import Crocodile.Resources.IOMethods as IOM
import Crocodile.Resources.Constants as CONST
import Crocodile.Resources.Mathematics as MATH
import Crocodile.Resources.Plotting as PL
import Crocodile.Resources.Functions as FU
import Crocodile.Resources.Equations as EQ

imp.reload(DCC)
imp.reload(IOM)


    
    

class MosquitoHelperMethods(DCC.dataclass):
    """
    Class with supporting functions for Mosquito. 
    
    """

    
    def __init__(self, objectname, measurement_method, flag_verbose = 0):
        self.verbose("New Mosquito class", flag_verbose)
                
        DCC.dataclass.__init__(self, objectname = objectname, measurement_method = measurement_method, flag_verbose = flag_verbose)  


    def import_data_show_shots(self, reload_data = False):
        """
        Import data from show shots. 
    
        INPUT:
        - reload_data (Bool, False): if False, try to import the numpy binary file first. If that is not present, or if True, it will import the original csv files. 
    
        DESCRIPTION:
        - 
    
        CHANGELOG:
        201604-RB: started function
    
        """
        if self.find_file_format() == False:
            return False

        # 0: pixels
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # 1: number of shots
        n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # 1: number of scans
        n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)

        # 2: datastates
        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
    
        # 3: spectra
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        if self.file_format < 4:
            self.n_sig = 6
            self.add_ds = False
        else:
            spds, n_sp_2, n_ds_2, self.n_sig, self.add_ds = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        total_shots = n_sh * n_sc
        n_signals = int(n_sc * n_sh / self.n_sig)
        
        empty_list = numpy.array([0])
        
        self.r_n = [n_w3, total_shots, 2*n_ds, n_sp, 1, 1, 1, 1]
        self.r = numpy.empty(self.r_n)
        self.r_axes = [w3_axis, numpy.arange(total_shots), spds[:,1], numpy.arange(n_sp), empty_list, empty_list, empty_list, empty_list]
        self.r_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "x", "x", "x", "Scans"]
        
        self.f_n = [n_w3, n_signals, n_ds, n_sp, 1, 1, 1, 1]
        self.f = numpy.zeros(self.f_n)
        self.f_axes = [w3_axis, numpy.arange(n_signals), numpy.arange(n_ds), numpy.arange(n_sp), empty_list, empty_list, empty_list, empty_list]
        self.f_units = ["w3 (cm-1)", "Signals", "Datastates", "Spectra", "x", "x", "x", "Scans"]
        
        self.s_n = [n_w3, n_signals, 1, n_sp, 1, 1, 1, 1]
        self.s = numpy.empty(self.s_n)
        self.s_axes = [w3_axis, numpy.arange(n_signals), empty_list, numpy.arange(n_sp), empty_list, empty_list, empty_list, empty_list]
        self.s_units = ["w3 (cm-1)", "Signals", "x", "Spectra", "x", "x", "x", "Scans"]
        
#         print(self.f[0,0,0,0,0,0,0,0])
        
        # import data
        
        paf = self.file_dict["base_filename"] + "_f.npy"
        
        if os.path.isfile(paf) and reload_data == False: 

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

                if _sc == 0:
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
                        self.r[:, s:e, _ds, _sp, 0,0,0,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)  


            x = 1                
            for _ch in range(n_choppers):      
                # if chop_max is not >0.5, the chopper was not on
                # if chop_min is not <0.5, it probably is the temperature alert. 
                chop_max = numpy.amax(self.r_choppers[_ch,:])
                chop_min = numpy.amin(self.r_choppers[_ch,:])
                
                if chop_max > 0.5 and chop_min < 0.5:
                    self.r_choppers[_ch,:] /= temp
                    self.r_choppers[_ch,:] *= x
                else: 
                    self.r_choppers[_ch,:] = 0
                x *= 2
            
            ds_list = numpy.sum(self.r_choppers, axis = 0).astype(int)
            
            n_signals = int(n_sc * n_sh / self.n_sig)
            
            ds_count = numpy.zeros(n_ds)
            for _ds in range(n_ds):
                ds_count[_ds] = len(numpy.where(ds_list == _ds)[0]) / n_signals

            count = numpy.zeros((n_signals, n_ds)) 

            for sig in range(n_signals):
                s = sig * self.n_sig
                e = (sig+1) * self.n_sig

                for _sh in range(self.n_sig):        
                    self.f[:, sig, ds_list[s+_sh],0,0,0,0,0] += (self.r[:, s+_sh, 2*ds_list[s+_sh], 0, 0, 0, 0, 0] / (self.r[:, s+_sh, 2*ds_list[s+_sh]+1, 0, 0, 0, 0, 0] * ds_count[ds_list[s+_sh]]))
                    
            paf = "{base_filename}_r.npy".format(base_filename = self.file_dict["base_filename"])
            numpy.save(paf, self.r)

            paf = "{base_filename}_f.npy".format(base_filename = self.file_dict["base_filename"])
            numpy.save(paf, self.f)
            
            paf = "{base_filename}_r_choppers.npy".format(base_filename = self.file_dict["base_filename"])
            numpy.save(paf, self.r_choppers)

            paf = "{base_filename}_r_specials.npy".format(base_filename = self.file_dict["base_filename"])
            numpy.save(paf, self.r_specials)
        
#         print(self.f[0,0,0,0,0,0,0,0])
        
        if self.add_ds:
        
            N = numpy.zeros(self.s_n)
            D = numpy.zeros(self.s_n)
            N_count = 0
            D_count = 0
            
            for ds in range(n_ds):
                if self.r_axes[2] == 1:
                    N[:,0,:,:,:,:,:,:] += self.f[:,ds,:,:,:,:,:,:]
                    N_count += 1
                else:
                    D[:,0,:,:,:,:,:,:] += self.f[:,ds,:,:,:,:,:,:]
                    D_count += 1
                
                if N_count > 0 and D_count > 0:
                    self.s = (N / N_count) / (D / D_count)
                elif N_count == 0 and D_count > 0:
                    self.s = 1 / (D / D_count)
                elif N_count > 0 and D_count == 0:
                    self.s = (N / N_count)
                else:
                    pass

        else:
        
            N = numpy.ones(self.s_n)
            D = numpy.ones(self.s_n)
            N_count = 0
            D_count = 0
            
            for ds in range(n_ds):
                if self.r_axes[2][ds] == 1:
                    N[:,:,0,:,:,:,:,:] *= self.f[:,:,ds,:,:,:,:,:]
                    N_count += 1
                else:
                    D[:,:,0,:,:,:,:,:] *= self.f[:,:,ds,:,:,:,:,:]
                    D_count += 1
                
                if N_count > 0 and D_count > 0:
                    self.s = (N / N_count) / (D / D_count)
                elif N_count == 0 and D_count > 0:
                    self.s = 1 / (D / D_count)
                elif N_count > 0 and D_count == 0:
                    self.s = (N / N_count)
                else:
                    pass


    def import_data_show_spectrum(self, **kwargs):
        """
        Import data from show spectrum. 
    
        INPUT:
        - -
    
        DESCRIPTION:
        Show spectrum has r (intensity) and s (signal). The individual scans are imported. The shots are averaged for every scan. There is no slow modulation or delays. 
        There are at least 2 datastates: probe and reference. Datastates are only saved when the checkbox "separate datastates" is checked in Messquito. 
        WORKAROUND NEEDED
        Contrary to most other methods the user can change settings during the measurement. This will lead to problems when importing the data. 
    
        CHANGELOG:
        201604-RB: started function
    
        """
        if self.find_file_format() == False:
            return False

        # 0: pixels
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # 1: number of shots
        n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # 1: number of scans
        n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)

        # 2: datastates
        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        #### FIX ####
        n_ds = 1
    
        # 3: spectra
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        if self.file_format < 4:
            self.n_sig = 6
            self.add_ds = False
        else:
            spds, n_sp_2, n_ds_2, self.n_sig, self.add_ds = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        empty_list = numpy.array([0])
        
        self.r_n = [n_w3, 1, 2*n_ds, n_sp, 1, 1, 1, n_sc]
        self.r = numpy.empty(self.r_n)
        self.r_noise = numpy.empty(self.r_n)
        self.r_axes = [w3_axis, empty_list, spds[:,1], numpy.arange(n_sp), empty_list, empty_list, empty_list, numpy.arange(n_sc)]
        self.r_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "x", "x", "x", "Scans"]
        
        self.s_n = [n_w3, 1, 1, n_sp, 1, 1, 1, n_sc]
        self.s = numpy.empty(self.s_n)
        self.s_noise = numpy.empty(self.s_n)
        self.s_axes = [w3_axis, empty_list, empty_list, numpy.arange(n_sp), empty_list, empty_list, empty_list, numpy.arange(n_sc)]
        self.s_units = ["w3 (cm-1)", "Signals", "x", "Spectra", "x", "x", "x", "Scans"]
        
        for _sc in range(n_sc):                  
            for _sp in range(self.r_n[3]): 

                # import signal
                suffix = "sp" + str(_sp) + "_ds0_signal_" + str(_sc)
                self.s[:,0,0,_sp, 0,0,0,_sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)  

                # import signal noise
                suffix = "sp" + str(_sp) + "_ds0_signal_" + str(_sc)
                self.s_noise[:,0,0,_sp, 0,0,0,_sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose) 

                for _ds in range(self.r_n[2]): 
                
                    # import probe and reference
                    suffix = "sp" + str(_sp) + "_ds" + str(_ds) + "_intensity_" + str(_sc)
                    self.r[:,0,_ds,_sp, 0,0,0,_sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)  

                    # import probe and reference noise
                    suffix = "sp" + str(_sp) + "_ds" + str(_ds) + "_intensity_" + str(_sc)
                    self.r_noise[:,0,_ds,_sp, 0,0,0,_sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)  
                

    def import_data_2dir(self, **kwargs):
#     import_temp_scans = False, t1_offset = 0):
        """
        Imports data from 2D-IR measurements. 
        
        INPUT:
        - import_temp_scans (Bool, False): If True, import individual scans instead of the averaged one. 
        - t1_offset (int, 0): if t0 is not where it should be, adjust it. The value is in integer and is the number of bins. 
        
        OUTPUT:
        
        DESCRIPTION:
        1. The function will first import the metadata. 
        2. The correct datastructures will be created. 
        3. The data itself will be imported.
        4. If needed the response will be calculated. 
        The data can be saved/imported in two ways:
        1. Save as signal or as intensity. This is determined by the user in LabVIEW. The value is n_datastates will determine this. 
        2. Import final, averaged data; or individual scans. The latter is only possible when that is saved. 
        
        CHANGELOG:
        201604/RB: extracted the important parts from a universal function
        
        """
        
        if "import_temp_scans" in kwargs and kwargs["import_temp_scans"]:
            import_temp_scans = True
        else:
            import_temp_scans = False
        
        if "t1_offset" in kwargs:
            t1_offset = kwargs["t1_offset"]
        else:
            t1_offset = 0
        
        if self.find_file_format() == False:
            return False

        # 0: pixels
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # 1: number of bins
        t1_bins, t1_fs, self.bin_sign, n_t1_bins, n_t1_fs, self.t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, t1_offset = t1_offset, flag_verbose = self.flag_verbose)
        
        # 2: datastates
        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
    
        # 3: spectra
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # spds is how the signal is/should be calculated
        spds, n_sp_2, n_ds_2, self.n_sig, self.add_ds = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        if n_ds == 0:
            self.measurement_type = "signal"
            n_ds = n_ds_2
        else:
            self.measurement_type = "intensity"
        
        # 4: slow modulation
        sm, sm_names, n_sm = IOM.import_slow_modulation(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        # 5: delays
        de, n_de = IOM.import_delays(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        # 6: dummy dimension
        n_du = 1
        du = numpy.arange(n_du)

        # 7: number of scans
        if import_temp_scans:
            n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
            sc = numpy.arange(n_sc)
        else:
            n_sc = 1
            sc = numpy.arange(n_sc)

        # DEFINE THE DATASTRUCTURES

        # probe and reference are saved separately
        if self.measurement_type == "intensity":

            self.b_n = [n_w3, n_t1_bins, 2*n_ds, n_sp, n_sm, n_de, n_du, n_sc]
            self.b_count_n = [n_t1_bins, n_ds, n_sp, n_sm, n_de, n_du, n_sc]
            self.b_intf_n = [n_t1_bins, n_ds, n_sp, n_sm, n_de, n_du, n_sc]

            self.b = numpy.empty(self.b_n)
            self.b_count = numpy.empty(self.b_count_n) 
            self.b_axes = [w3_axis, t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
            self.b_units = ["w3 (cm-1)", "T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

            self.b_intf = numpy.zeros(self.b_intf_n)
            self.b_intf_axes = [t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
            self.b_intf_units = ["T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]
        
        else:
            n_ds = 1

        self.r_n = [n_w3, n_t1_fs, n_ds, n_sp, n_sm, n_de, n_du, n_sc]
        self.r = numpy.zeros(self.r_n)
        self.r_axes = [w3_axis, t1_fs, [0], spds[:,0], sm, de, du, sc]
        self.r_units = ["w3 (cm-1)", "T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]
        
        self.r_intf_n = [n_t1_bins, n_ds, n_sp, n_sm, n_de, n_du, n_sc]
        self.r_intf = numpy.zeros(self.r_intf_n)
        self.r_intf_axes = [t1_bins, spds[:,1], spds[:,0], sm, de, du, sc]
        self.r_intf_units = ["T1 (bins)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.f_n = [n_w3, n_t1_fs, 2, n_sp, n_sm, n_de, n_du, n_sc]
        self.f = numpy.zeros(self.f_n)
        self.f_axes = [w3_axis, [0], [0], spds[:,0], sm, de, du, sc]
        self.f_units = ["w3 (cm-1)", "w1 (cm-1)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        self.s_n = [n_w3, n_t1_fs, 1, n_sp, n_sm, n_de, n_du, n_sc]
        self.s = numpy.zeros(self.s_n)
        self.s_axes = [w3_axis, [0], [0], spds[:,0], sm, de, du, sc]
        self.s_units = ["w3 (cm-1)", "w1 (cm-1)", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]
        
        # IMPORT THE DATA
        
        # signal, no temp scans
        if self.measurement_type == "signal" and self.r_n[7] == 1:  
            self.verbose("Average scan, signal", self.flag_verbose)
                
            for sp in range(self.r_n[3]):
                for sm in range(self.r_n[4]):
                    for de in range(self.r_n[5]): 
                        for du in range(self.r_n[6]):
                    
                            # import files
                            suffix = "signal_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                            self.r[:,:,0,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T[:,self.t1_zero_index:self.t1_zero_index + self.r_n[1]]
                            # the portion of the data that is relevant is imported, i.e. t1>0 and the number of bins we actually need (is an even number)

                            # import interferogram
                            suffix = "interferogram_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                            self.r_intf[:,0,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T
                            
        # signal, temp scans
        elif self.measurement_type == "signal" and self.r_n[7] > 1:
            self.verbose("All scans, signal", self.flag_verbose)
            for sp in range(self.r_n[3]):
                for sm in range(self.r_n[4]):
                    for de in range(self.r_n[5]): 
                        for du in range(self.r_n[6]):
                            for sc in range(self.r_n[7]):
                    
                                # import files
                                suffix = "signal_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du) + "_" + str(sc)
                                self.r[:,:,0,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T[:,self.t1_zero_index:self.t1_zero_index + self.r_n[1]]

                                # import interferogram
                                suffix = "interferogram_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du) + "_" + str(sc)
                                self.r_intf[:,0,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T                           

        # intensity, no temp scans
        elif self.measurement_type == "intensity" and self.b_n[7] == 1:
            self.verbose("Average scan, intensity", self.flag_verbose)
            
            # probe and reference are two datastates but imported in one go
            n_ds = int(self.b_n[2]/2)
            
            for ds in range(n_ds):    
                for sp in range(self.b_n[3]):
                    for sm in range(self.b_n[4]):
                        for de in range(self.b_n[5]): 
                            for du in range(self.b_n[6]):
                    
                                # import probe
                                suffix = "probe_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b[:,:,2*ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

                                # import reference
                                suffix = "reference_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b[:,:,2*ds+1,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

                                # import count
                                suffix = "count_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b_count[:,ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1)

                                # import interferogram
                                suffix = "interferogram_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                self.b_intf[:,ds,sp,sm,de,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1)

        # intensity, temp scans
        elif self.measurement_type == "intensity" and self.b_n[7] > 1:
            self.verbose("All scans, intensity", self.flag_verbose)
            
            # probe and reference are two datastates but imported in one go
            n_ds = int(self.b_n[2]/2)
        
            for ds in range(n_ds):    
                for sp in range(self.b_n[3]):
                    for sm in range(self.b_n[4]):
                        for de in range(self.b_n[5]): 
                            for du in range(self.b_n[6]):
                                for sc in range(self.b_n[7]):
                    
                                    # import probe
                                    suffix = "probe_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b[:,:,2*ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

                                    # import reference
                                    suffix = "reference_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b[:,:,2*ds+1,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

                                    # import count
                                    suffix = "count_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b_count[:,:,ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

                                    # import interferogram
                                    suffix = "interferogram_ds" + str(ds) + "_sp" + str(sp) + "_sm" + str(sm) + "_de" + str(de) + "_du" + str(du)
                                    self.b_intf[:,ds,sp,sm,de,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

        else:
            self.printError("Failed to import measurement files", inspect.stack())   
            
        # calculate r
        if self.measurement_type == "intensity":
            self.b_to_r()

        # if needed, reverse the bin axis
        if self.bin_sign and self.measurement_type == "intensity":
            self.b = self.b[:,::-1,:,:,:,:,:,:]
            self.b_intf = self.b_intf[::-1,:,:,:,:,:,:]
            self.b_count = self.b_count[:,::-1,:,:,:,:,:,:]            


    def import_data_scan_spectrum(self):
        """
        Import data from Scan Spectrum
    
        INPUT:
        - 
    
        OUTPUT:
        - 
    
        DESCRIPTION:
        The pixel axis is now the wavelengths scanned
        
    
        CHANGELOG:
        201604-RB: started function
    
        """
    
        wl, n_wl = IOM.import_wavelengths(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        wl = 1e7 / wl
        
        empty = numpy.arange(1)
        self.r_n = [n_wl, 1, 2, 1, 1, 1, 1, 1]
        self.r = numpy.empty(self.r_n)
        self.r_axes = [wl, empty, numpy.arange(2), empty, empty, empty, empty, empty]
        self.r_units = ["wavenumber (cm-1)", "x", "Datastates", "x", "x", "x", "x", "x"] 
        
        suffix = "data_0"
        self.r[:,0,:,0,0,0,0,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T
        

    def import_data_pump_probe(self, **kwargs):
        
        if "import_temp_scans" in kwargs and kwargs["import_temp_scans"]:
            import_temp_scans = True
        else:
            import_temp_scans = False
        
        if "t1_offset" in kwargs:
            t1_offset = kwargs["t1_offset"]
        else:
            t1_offset = 0

        if self.find_file_format() == False:
            return False

        # 0: pixels
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # 1: number of bins
#         t1_bins, t1_fs, self.bin_sign, n_t1_bins, n_t1_fs, self.t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, t1_offset = t1_offset, flag_verbose = self.flag_verbose)
        
        # 2: datastates
        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
    
        # 3: spectra
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        # spds is how the signal is/should be calculated
        spds, n_sp_2, n_ds_2, self.n_sig, self.add_ds = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        if n_ds == 0:
            self.measurement_type = "signal"
            n_ds = n_ds_2
        else:
            self.measurement_type = "intensity"
        
        # 4: slow modulation
        sm, sm_names, n_sm = IOM.import_slow_modulation(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        print(sm_names)

        # 5: delays
        de, n_de = IOM.import_delays(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

        # 6: dummy dimension
        n_du = 1
        du = numpy.arange(n_du)

        # 7: number of scans
        if import_temp_scans:
            n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
            sc = numpy.arange(n_sc)
        else:
            n_sc = 1
            sc = numpy.arange(n_sc)


        # DEFINE THE DATASTRUCTURES

        # probe and reference are saved separately
        if self.measurement_type == "intensity":

            self.r_n = [n_w3, 1, n_ds, n_sp, n_sm, n_de, n_du, n_sc]
            self.r = numpy.zeros(self.r_n)
            self.r_noise = numpy.zeros(self.r_n)
            self.r_axes = [w3_axis, [0], [0], spds[:,0], sm, de, du, sc]
            self.r_units = ["w_probe (cm-1)", "x", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]

        else:
            n_ds = 1

        self.s_n = [n_w3, 1, 1, n_sp, n_sm, n_de, n_du, n_sc]
        self.s = numpy.zeros(self.s_n)
        self.s_noise = numpy.zeros(self.s_n)
        self.s_axes = [w3_axis, [0], [0], spds[:,0], sm, de, du, sc]
        self.s_units = ["w_probe (cm-1)", "x", "Datastates", "Spectra", "Slow modulation", "Delays (fs)", "Dummies", "Scans"]



        # IMPORT THE DATA
        
        # signal, no temp scans
        if self.measurement_type == "signal" and self.s_n[7] == 1:  
            self.verbose("Average scan, signal", self.flag_verbose)
                
            for sp in range(self.s_n[3]):
                for sm in range(self.s_n[4]):
                    for du in range(self.s_n[6]):
                
                        # import files
                        suffix = "signal_sp" + str(sp) + "_sm" + str(sm) + "_du" + str(du)
                        if self.s_n[5] == 1:
                            self.s[:,0,0,sp,sm,0,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1)

                        else:
                            self.s[:,0,0,sp,sm,:,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

                        # import noise
                        suffix = "signal_noise_sp" + str(sp) + "_sm" + str(sm) + "_du" + str(du)
                        if self.s_n[5] == 1:
                            self.s_noise[:,0,0,sp,sm,0,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1)

                        else:
                            self.s_noise[:,0,0,sp,sm,:,du,0] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T


        # signal, temp scans
        if self.measurement_type == "signal" and self.s_n[7] == 1:  
            self.verbose("Average scan, signal", self.flag_verbose)
                
            for sp in range(self.s_n[3]):
                for sm in range(self.s_n[4]):
                    for du in range(self.s_n[6]):
                        for sc in range(self.s_n[7]):
                
                            # import files
                            suffix = "signal_sp" + str(sp) + "_sm" + str(sm) + "_du" + str(du) + "_" + str(sc)
                            if self.s_n[5] == 1:
                                self.s[:,0,0,sp,sm,0,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1)

                            else:
                                self.s[:,0,0,sp,sm,:,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T

                            # import noise
                            suffix = "signal_noise_sp" + str(sp) + "_sm" + str(sm) + "_du" + str(du)
                            if self.s_n[5] == 1:
                                self.s_noise[:,0,0,sp,sm,0,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1)

                            else:
                                self.s_noise[:,0,0,sp,sm,:,du,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose - 1).T






        
    def find_file_format(self, flag_verbose = False):
        """
        Find file format.
    
        INPUT:
        - flag_verbose: deprecated
    
        OUTPUT:
        - True or False: If False, the file format file was not found.
    
        DESCRIPTION:
        Changes to the file format affect the importing scripts. Using the file format the scripts can be made to work for both the current and older formats. 
    
        CHANGELOG:
        201604-RB: started function
    
        """
    
        if self.flag_verbose:
            print("MosquitoHelper.find_file_format")
    
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

        if self.flag_verbose:
            print("LV fileformat: {a}".format(a = self.file_format))
            print("DONE: MosquitoHelper.find_file_format\n")
  
        return True


    def no_count_cheat(self):
        
        ix = numpy.where(self.b_count == 0)[0]
        self.b_count[ix] = 1
        self.b_to_r()
        
        DEBUG.printWarning("NO COUNT CHEAT: count 0 set to 1. USE FOR DEBUGGING ONLY!", inspect.stack())
        
    
    def b_to_r(self):
        """
        B to R. B is the raw data and is not yet divided by the count. 
    
        INPUT:
        - if self.measurement_type == 'signal' 
    
        OUTPUT:
        - 
    
        DESCRIPTION:
    
    
        CHANGELOG:
        201604-RB: started function
    
        """
        
        if self.measurement_type == "intensity":
            
            b_size = numpy.array(numpy.shape(self.b))
            b_size[2] = 1
            b_temp = numpy.zeros(b_size)
            b_intf_temp = numpy.zeros(b_size[1:])
            
            n_ds = int(self.b_n[2]/2)
            
            # select the relevant bins (i.e. t1 > 0 and an even length)
            s = self.t1_zero_index
            e = self.t1_zero_index + self.r_n[1]

            if numpy.all(self.b_count > 0):

                for ds in range(n_ds):
                    self.r_intf[:,0,:,:, :,:,:,:] += self.b_intf[:,ds,:,:, :,:,:,:] / self.b_count[:,ds,:,:, :,:,:,:]
                    
                self.r_intf /= n_ds
            
                if self.add_ds:
            
                    N = numpy.zeros(self.r_n)
                    D = numpy.zeros(self.r_n)

                    for pi in range(self.b_n[0]): 
                        N_count = 0
                        D_count = 0
                        for ds in range(self.r_n[2]):
                            if self.b_axes[2][ds] == 1:
                                N[pi,:,0,:,:,:,:,:] += self.b[pi,s:e,2*ds,:,:,:,:,:] / self.b[pi,s:e,2*ds+1,:, :,:,:,:]
                                N_count += 1
                            else:
                                D[pi,:,0,:,:,:,:,:] += self.b[pi,s:e,2*ds,:,:,:,:,:] / self.b[pi,s:e,2*ds+1,:, :,:,:,:]
                                D_count += 1
                        
                    if N_count > 0 and D_count > 0:
                        self.r = (N / N_count) / (D / D_count)
                    elif N_count == 0 and D_count > 0:
                        self.r = 1 / (D / D_count)
                    elif N_count > 0 and D_count == 0:
                        self.r = (N / N_count)
                    else:
                        pass
                        

                else:
            
                    N = numpy.ones(self.r_n)
                    D = numpy.ones(self.r_n)

                    for pi in range(self.b_n[0]): 
                        for ds in range(self.b_n[2]):
                            if self.b_axes[2] == 1:
                                N[pi,:,0,:,:,:,:,:] *= self.b[pi,:,2*ds,:,:,:,:,:,:] / self.b[pi,:,2*ds+1,:, :,:,:,:,:]
                            else:
                                D[pi,:,0,:,:,:,:,:] *= self.b[pi,:,2*ds,:,:,:,:,:,:] / self.b[pi,:,2*ds+1,:, :,:,:,:,:]

        else:
#             DEBUG.
            pass


        nan_warning = 0
        for bi in range(self.b_count_n[0]):
            
            if numpy.all(self.b_count[bi,:,:, :,:,:,:] > 0):
                self.r_intf[bi,:,:, :,:,:,:] = self.b_intf[bi,:,:, :,:,:,:] / self.b_count[bi,:,:, :,:,:,:]
            else:
                self.r_intf[bi,:,:, :,:,:,:] = numpy.nan
                nan_warning += 1
                
        if nan_warning:
            mean = numpy.nanmean(self.r_intf)
            ix = numpy.where(numpy.isnan(self.r_intf[:,0,0, 0,0,0,0]))[0]
            self.r_intf[ix] = mean
                
        print("nan warning:", nan_warning)
        

    def calculate_phase(self, **kwargs):
        """
        This methods does the FFT of the interferometer. It looks for the peak in the real part. It then uses n_points, centered around the peak, to calculate the phase.
        
        
        INPUT:
        - n_points (int, default 5): number of points used to calculate the phase. If even, it looks for the highest neigbor of the peak and centers around that. 
        - range (list with length 2, default: [0,-1]). The wavenumber range to look for the phase. 
        - flag_plot (Bool, False): Makes a plot of (1) the interferograms of each separate measurement, (2) the FFT of the average of those and (3) the angle.   
        
        EXAMPLES:
        For a real part with values: [3,1,2,5,4,3]
        n_points = 2: [5,4]
        n_points = 3: [2,5,4]
        n_points = 4: [2,5,4,3]
        
        """
        if "n_points" in kwargs:
            n_points = kwargs["n_points"]
        else:
            n_points = 5
            
        if "w_range" in kwargs:
            w_range = kwargs["w_range"]
        else:
            w_range = [0,-1]
            
        if "flag_plot" in kwargs and kwargs["flag_plot"]:
            flag_plot = True
        else:
            flag_plot = False
            
        # an even number of points
        N_bins = self.r_intf_n[0] 
        N_bins_half = int(N_bins/2)
        N_bins = 2 * N_bins_half
        
        # not the same as w1 axis, that one is truncated to exclude t1 < 0
        dt = self.r_axes[1][1] - self.r_axes[1][0]
        w_axis = MATH.make_ft_axis(N_bins, dt = dt, undersampling = 0, normalized_to_period = 0, zero_in_middle = False, flag_verbose = False)
        w_axis = w_axis[:N_bins_half]
        
        # if w_range is given, find the range of interest
        # needed when there are overtones
        i_range = [0,-1]
        if w_range != [0,-1]:
            i_range[0] = numpy.where(w_axis > w_range[0])[0][0]
            i_range[1] = numpy.where(w_axis > w_range[1])[0][0]      
            self.verbose("calculate_phase: peak searching between indices %i and %i (%.1f and %.1f cm-1)" % (i_range[0], i_range[1], w_axis[i_range[0]], w_axis[i_range[1]]), self.flag_verbose)
        
        # average the interferograms for all the measurements
        # it shouldn't change from delay to delay
        r_intf = numpy.zeros(self.r_intf_n[0])
        for bi in range(self.r_intf_n[0]): 
            r_intf[bi] = numpy.nanmean(self.r_intf[bi,:,:,:,:,:,:])
        r_intf -= numpy.nanmean(r_intf)
        # replace NaN by zero
        r_intf = numpy.nan_to_num(r_intf)
        
        # do the FFT
        r_intf_roll = numpy.roll(r_intf, -self.t1_zero_index)
        f = numpy.fft.fft(r_intf_roll)
        f = f[:N_bins_half] 
        
        # do the actual phase calculation
        angle = self.find_phase(f, w_axis, n_points, i_range)   

        if flag_plot:
    
            fig = plt.figure()
            ax_n = 3
            ax = [0] * ax_n
        
            for ax_i in range(ax_n):
                ax[ax_i] = fig.add_subplot(3,1,ax_i+1)   

            ax_i = 0
            ax[ax_i].plot(r_intf)        

            for sp in range(self.s_n[3]):
                for sm in range(self.s_n[4]):
                    for de in range(self.s_n[5]): 
                        for du in range(self.s_n[6]):
                            for sc in range(self.s_n[7]):  
                                ax[ax_i].plot(self.r_intf[:, 0, sp, sm, de, du, sc])

            ax_i = 1
            ax[ax_i].plot(w_axis, numpy.real(f))
            ax[ax_i].plot(w_axis, numpy.imag(f))
        
            ax_i = 2
            ax[ax_i].plot(w_axis, angle)

            ax_i = 0
            ax[ax_i].set_xlabel("Bins")
            ax[ax_i].set_ylabel("V")

            ax_i = 1
            ax[ax_i].set_xlabel("w1")
            ax[ax_i].set_ylabel("FT(V)")
            
            ax_i = 2
            ax[ax_i].set_xlabel("w1")
            ax[ax_i].set_ylabel("Angle")


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
        
        # find the maximum value
        y = numpy.abs(f)
        i_max = numpy.argmax(y[i_range[0]:i_range[1]])
        idx = numpy.arange(i_max, i_max + n_points)
        i_max += i_range[0]
        
        # take some points on both sides of the maximum
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
            
        self.f_intf_angle = angle

        # make a linear fit of the angle
        # the slope gives the bin-offset
        A_start = [numpy.mean(angle[idx]), 1]
        A_out = MATH.fit(w_axis[idx], angle[idx], EQ.linear, A_start)

        shift = A_out[1] / (2 * numpy.pi * 1e-15 * CONST.c_cms * CONST.hene_fringe_fs)
        
        if shift > 1 or shift < -1:
            self.printWarning("=== CHANGE ZERO BIN!!! ===")
            print("Current zero bin: {bin}".format(bin =  self.t1_zero_index))
            print("Best zero bin: {bin:4.1f}".format(bin =  self.t1_zero_index - shift))
            print("Recommended shift is: {bin:4.0f} bins".format(bin = shift))
            print("Phase below is for the current zero bin.")

        print("Phase in degrees: %.1f" % (self.phase_rad * 180 / numpy.pi))
        
        return angle


    def make_fft(self, **kwargs): 
        """
    
        INPUT:
        definition: n_bins_t1: the number of bins where t1>=0
        - phase_cheat_deg (num, 0): correct problems with the phase. The phase will only be adjusted for this particular function call. 
        - zeropad_to (int, -1): zeropad to a number of points. This can be shorter than n_bins_t1. If value is <0 (default), n_bins_t1 will be used 
        - zeropad_by (num, -1): zeropad to a length zeropad_by * n_bins_t1. 1 would be n_bins_t1. If value is <0 (default), n_bins_t1 will be used. If zeropad_by is given, this argument will be ignored. 
        - window_function (str, False): choose between none, ones, triangular, gaussian. 

        OUTPUT:
        - none
    
        DESCRIPTION:
    
    
        CHANGELOG:
        201604-RB: started function
    

        """

        if "phase_cheat_deg" in kwargs:
            phase_cheat_deg = kwargs["phase_cheat_deg"]
        else:
            phase_cheat_deg = 0

        if "zeropad_to" in kwargs:
            zeropad_to = kwargs["zeropad_to"]
        else:
            zeropad_to = -1
            
        if "zeropad_by" in kwargs:
            zeropad_by = kwargs["zeropad_by"]
        else:
            zeropad_by = -1
            
        if "window_function" in kwargs:
            window_function = kwargs["window_function"]
        else:
            window_function = False
        
        # cheat the phase. 
        phase_rad = self.phase_rad + phase_cheat_deg * numpy.pi / 180
        if phase_cheat_deg != 0:
            self.printWarning("Phase cheating with %.1f degrees. Used phase is %.1f" % (phase_cheat_deg, phase_rad * 180 / numpy.pi))


        if zeropad_to > 0:
            self.zeropad_to = zeropad_to
        elif zeropad_by > 0:
            self.zeropad_by = zeropad_by
        else:
            self.zeropad_to = self.f_n[1]
        
        N_FFT_bins = self.zeropad_to #self.r_n[1]
        N_FFT_half = int(N_FFT_bins / 2)
        N_FFT_bins = 2 * N_FFT_half
        self.zeropad_to = N_FFT_bins

        w1_axis = MATH.make_ft_axis(N_FFT_bins, dt = self.r_axes[1][1] - self.r_axes[1][0], undersampling = 0, normalized_to_period = 0, zero_in_middle = False, flag_verbose = self.flag_verbose)
        
        self.f_axes[1] = w1_axis
        self.s_axes[1] = w1_axis[:N_FFT_half]

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
                                # for zeropadding and windowfunctions, the response should be around zero, so subtract the mean. 
                                r -= numpy.mean(r)
                                if window_function:
                                    r = MATH.window_functions(r, window_function, N_FFT_bins)
                                
                                f = numpy.fft.fft(r, n = self.zeropad_to)                      
                                self.f[pi, :, ds, sp, sm, de, du, sc] = f
                                s = numpy.real(numpy.exp(-1j * phase_rad) * f[:N_FFT_half])
                                self.s[pi, :, ds, sp, sm, de, du, sc] = s
                     

    def check_axis(self, ax):
        """
    
    
        INPUT:
        - 
    
        OUTPUT:
        - 
    
        DESCRIPTION:
    
    
        CHANGELOG:
        201604-RB: started function
    
        """       
        new_axis = False
        if ax == False:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            new_axis = True
            
        return ax, new_axis


    def multiplot_ranges_helper(self, key, index, name, **kwargs):
        """
        Checks if KEY is in KWARGS. 
    
        INPUT:
        - 
    
        OUTPUT:
        - 
    
        DESCRIPTION:
        - should be improved with checks.
    
        CHANGELOG:
        20160422-RB: started function        
        
        """
        if key in kwargs:
            if type(kwargs[key]) == int:
                res = numpy.arange(kwargs[key])
            else:
                res = kwargs[key]
            
            if len(res) >= self.s_n[index]:  
                self.printWarning("Index {index} ({name}) is out of bounds with size {a} (max is {b})".format(index = index, name = name, a = len(res), b = self.s_n[index]))
                res = res[:self.s_n[index]]
                
        else:
            res = numpy.arange(self.s_n[index])       
            
        return res


    def multiplot_ranges(self, **kwargs):
        """
        Checks for keywords pi (pixels), bish (bins/shots), sp (spectra), ds (datastates), sm (slow modulation), de (delays), du (dummies), sc (scans). The value can be a list or an integer with indices. If the indices are out of range, they will be ignored. 
    
        INPUT:
        - kwargs. 
    
        OUTPUT:
        - ranges. 
    
        CHANGELOG:
        201604-RB: started function
        20160422-RB: wrote helper function
    
        """

        pi = self.multiplot_ranges_helper("pi", 0, "pixels", **kwargs)
        
        bish = self.multiplot_ranges_helper("bish", 1, "bins/shots", **kwargs)
        
        sp = self.multiplot_ranges_helper("sp", 2, "spectra", **kwargs)
        
        ds = self.multiplot_ranges_helper("ds", 3, "datastates", **kwargs)
        
        sm = self.multiplot_ranges_helper("sm", 4, "slow modulation", **kwargs)
        
        de = self.multiplot_ranges_helper("de", 5, "delays", **kwargs)
        
        du = self.multiplot_ranges_helper("du", 6, "dummies", **kwargs)
        
        sc = self.multiplot_ranges_helper("sc", 7, "scans", **kwargs)

        return pi, bish, sp, ds, sm, de, du, sc


    def sum_t1(self, x_axis, y_axis, data, **kwargs):

        if "x_range" in kwargs:
            x_range = kwargs["x_range"]
        else:
            x_range = [0,0]

        if "y_range" in kwargs:
            y_range = kwargs["y_range"]
        else:
            y_range = [0,-1]

        # determine the range to be plotted
        x_min, x_max, y_min, y_max = FU.find_axes(x_axis, y_axis, x_range, y_range, self.flag_verbose)
    
        # find the area to be plotted
        x_min_i, x_max_i = FU.find_axes_indices(x_axis, x_min, x_max)
        y_min_i, y_max_i = FU.find_axes_indices(y_axis, y_min, y_max)
          
        # truncate the data
        data, x_axis, y_axis = FU.truncate_data(data, y_axis, x_axis, y_min_i, y_max_i, x_min_i, x_max_i) 
        
        data = numpy.mean(data, axis = 1)
        
        return data
              

    def export_2d_data_for_gnuplot_helper(self, x_axis, y_axis, data, path_and_filename, **kwargs):
        """
        
    
        INPUT:
        - data:
        - x_range, y_range (optional)
    
        OUTPUT:

    
        CHANGELOG:
        20160502-RB: started function

        """
        
        if "x_range" in kwargs:
            x_range = kwargs["x_range"]
        else:
            x_range = [0,0]

        if "y_range" in kwargs:
            y_range = kwargs["y_range"]
        else:
            y_range = [0,-1]

        # determine the range to be plotted
        x_min, x_max, y_min, y_max = FU.find_axes(x_axis, y_axis, x_range, y_range, self.flag_verbose)
    
        # find the area to be plotted
        x_min_i, x_max_i = FU.find_axes_indices(x_axis, x_min, x_max)
        y_min_i, y_max_i = FU.find_axes_indices(y_axis, y_min, y_max)
          
        # truncate the data
        data, x_axis, y_axis = FU.truncate_data(data, y_axis, x_axis, y_min_i, y_max_i, x_min_i, x_max_i) 
        
        n_x = len(x_axis)
        n_y = len(y_axis)
            
        with open(path_and_filename, 'w') as f:
        
            for x_i in range(n_x):
                for y_i in range(n_y):
                    f.write("{y}, {x}, {z}\n".format(y = y_axis[y_i], x = x_axis[x_i], z = data[y_i, x_i]))
            
                f.write("\n")
            
        f.close()
        
        
    def check_results_folder(self):
        """
        Check if the results folder exists, if not, make it. 
        """
        
        try:
            os.stat(self._file_dict["result_folder"])
        except:
            os.mkdir(self._file_dict["result_folder"])  
        

    def remove_pixel(self, remove_idx):
        """
         
        INPUT:
        - remove_idx (list with int): indices of the pixels that have to be removed. Go from low to high index.
    
        OUTPUT:

        COMMENTS:

        
        CHANGELOG:
        20160525-RB: started function

        """   
        # when we delete pixels, the indices change. Start at the highest index.
        remove_idx = remove_idx[::-1]
    
        self.printWarning("=== REMOVING PIXELS ===")
        if self.flag_verbose:
            for i in remove_idx:
                print("Removing pixel {index} ({wn:4.2f} cm-1)".format(index = i, wn = self.s_axes[0][i]))
        else:   
            s = "Removing pixels:"
            for i in remove_idx:
                s = ("{string} {index},".format(string = s, index = i))
            print(s)

        # make numpy array    
        if type(remove_idx) == list:
            remove_idx = numpy.array(remove_idx)
        
        # delete the pixels 
        if hasattr(self, 's'):
            if hasattr(self, 's_noise'):
                self.s, self.s_axes[0], self.s_n[0], self.s_noise = self.remove_pixel_helper(data = self.s, axis = self.s_axes[0], n = self.s_n[0], remove_idx = remove_idx, noise = self.s_noise)
            else:
                self.s, self.s_axes[0], self.s_n[0], dump = self.remove_pixel_helper(data = self.s, axis = self.s_axes[0], n = self.s_n[0], remove_idx = remove_idx, noise = -1)

        if hasattr(self, 'r'):
            if hasattr(self, 'r_noise'):
                self.r, self.r_axes[0], self.r_n[0], self.r_noise = self.remove_pixel_helper(data = self.r, axis = self.r_axes[0], n = self.r_n[0], remove_idx = remove_idx, noise = self.r_noise)
            else:
                self.r, self.r_axes[0], self.r_n[0], dump = self.remove_pixel_helper(data = self.r, axis = self.r_axes[0], n = self.r_n[0], remove_idx = remove_idx, noise = -1)

        if hasattr(self, 'f'):
            if hasattr(self, 'f_noise'):
                self.f, self.f_axes[0], self.f_n[0], self.f_noise = self.remove_pixel_helper(data = self.f, axis = self.f_axes[0], n = self.f_n[0], remove_idx = remove_idx, noise = self.f_noise)
            else:
                self.f, self.f_axes[0], self.f_n[0], dump = self.remove_pixel_helper(data = self.f, axis = self.f_axes[0], n = self.f_n[0], remove_idx = remove_idx, noise = -1)

            
    def remove_pixel_helper(self, data, axis, n, remove_idx, noise = -1):
    
        data = numpy.delete(data, remove_idx, 0)
        if type(noise) == numpy.ndarray:
            noise = numpy.delete(noise, remove_idx, 0) 
        axis = numpy.delete(axis, remove_idx, 0)
        # the axis is shorter
        n -= len(remove_idx)  
        
        return data, axis, n, noise
            
            
        


    def average_pixels(self, average_idx):
        """
           
        INPUT:
        - average_idx (list with list with 2 ints): list with lists of the 2 pixels that have to be averaged. Go from low to high index.
    
        OUTPUT:

    
        CHANGELOG:
        20160525-RB: started function

        """    
        # when we delete pixels, the indices change. Start at the highest index.
        average_idx = average_idx[::-1]

        self.printWarning("=== AVERAGING PIXELS ===")
        if self.flag_verbose:
            for combi in average_idx:
                print("Averaging pixels {index1} ({wn1:4.2f} cm-1) and {index2} ({wn2:4.2f} cm-1)".format(index1 = combi[0], wn1 = self.s_axes[0][combi[0]], index2 = combi[1], wn2 = self.s_axes[0][combi[1]]))
        else:   
            s = "Averaging pixels:"
            for combi in average_idx:
                s = ("{string} {index1} and {index2};".format(string = s, index1 = combi[0], index2 = combi[1]))
            print(s)


        if hasattr(self, 's'):
            if hasattr(self, 's_noise'):
                self.s, self.s_axes[0], self.s_n[0], self.s_noise = self.average_pixels_helper(data = self.s, axis = self.s_axes[0], n = self.s_n[0], average_idx = average_idx, noise = self.s_noise)
            else:
                self.s, self.s_axes[0], self.s_n[0], dump = self.average_pixels_helper(data = self.s, axis = self.s_axes[0], n = self.s_n[0], average_idx = average_idx, noise = -1)       

        if hasattr(self, 'r'):
            if hasattr(self, 'r_noise'):
                self.r, self.r_axes[0], self.r_n[0], self.r_noise = self.average_pixels_helper(data = self.r, axis = self.r_axes[0], n = self.r_n[0], average_idx = average_idx, noise = self.r_noise)
            else:
                self.r, self.r_axes[0], self.r_n[0], dump = self.average_pixels_helper(data = self.r, axis = self.r_axes[0], n = self.r_n[0], average_idx = average_idx, noise = -1)

        if hasattr(self, 'f'):
            if hasattr(self, 'f_noise'):
                self.f, self.f_axes[0], self.f_n[0], self.f_noise = self.average_pixels_helper(data = self.f, axis = self.f_axes[0], n = self.f_n[0], remove_idx = remove_idx, noise = self.f_noise)
            else:
                self.f, self.f_axes[0], self.f_n[0], dump = self.average_pixels_helper(data = self.f, axis = self.f_axes[0], n = self.f_n[0], average_idx = average_idx, noise = -1)

        if hasattr(self, 'b'):
            if hasattr(self, 'b_noise'):
                self.b, self.b_axes[0], self.b_n[0], self.b_noise = self.average_pixels_helper(data = self.b, axis = self.b_axes[0], n = self.b_n[0], remove_idx = remove_idx, noise = self.b_noise)
            else:
                self.b, self.b_axes[0], self.b_n[0], dump = self.average_pixels_helper(data = self.b, axis = self.b_axes[0], n = self.b_n[0], average_idx = average_idx, noise = -1)


    def average_pixels_helper(self, data, axis, n, average_idx, noise = -1):

        for combi in average_idx:
            # average the pixels and delete one of them
            data[combi[0]] = (data[combi[0]] + data[combi[1]]) / 2
            data = numpy.delete(data, combi[1], 0)
            # the same but for the axes
            axis[combi[0]] = (axis[combi[0]] + axis[combi[1]]) / 2
            axis = numpy.delete(axis, combi[1], 0)
            
            if type(noise) == numpy.ndarray:
                noise[combi[0]] = (noise[combi[0]] + noise[combi[1]]) / 2
                noise = numpy.delete(noise, combi[1], 0)   
        
        # the pixel axis is shorter
        n -= len(average_idx)
        
        return data, axis, n, noise


    def merge_spx_modulation(self, axes, sm_spx, idx_for_scaling = False):
        """
           
        INPUT:
        - axes (list with ndarray): list with the axes to combine
        - sm_spx (list): where the slow modulation states go. 
        idx_for_scaling (list): indices for the scaling. NOT IMPLEMENTED YET.
    
        OUTPUT:

    
        CHANGELOG:
        20160525-RB: started function

        """           

        if hasattr(self, 's'):
            if hasattr(self, 's_noise'):
                self.s, self.s_axes, self.s_n, self.s_noise = self.merge_spx_modulation_helper(data = self.s, axis = self.s_axes, n = self.s_n, axes = axes, sm_spx = sm_spx, noise = self.s_noise)
            else:
                self.s, self.s_axes, self.s_n, dump = self.merge_spx_modulation_helper(data = self.s, axis = self.s_axes, n = self.s_n, axes = axes, sm_spx = sm_spx, noise = -1)
                
        if hasattr(self, 'r'):
            if hasattr(self, 'r_noise'):
                self.r, self.r_axes, self.r_n, self.r_noise = self.merge_spx_modulation_helper(data = self.r, axis = self.r_axes, n = self.r_n, axes = axes, sm_spx = sm_spx, noise = self.r_noise)
            else:
                self.r, self.r_axes, self.r_n, dump = self.merge_spx_modulation_helper(data = self.r, axis = self.r_axes, n = self.r_n, axes = axes, sm_spx = sm_spx, noise = -1)

        if hasattr(self, 'f'):
            if hasattr(self, 'f_noise'):
                self.f, self.f_axes, self.f_n, self.f_noise = self.merge_spx_modulation_helper(data = self.f, axis = self.f_axes, n = self.f_n, axes = axes, sm_spx = sm_spx, noise = self.f_noise)
            else:
                self.f, self.f_axes, self.f_n, dump = self.merge_spx_modulation_helper(data = self.f, axis = self.f_axes, n = self.f_n, axes = axes, sm_spx = sm_spx, noise = -1)
                
        if hasattr(self, 'b'):
            if hasattr(self, 'b_noise'):
                self.b, self.b_axes, self.b_n, self.b_noise = self.merge_spx_modulation_helper(data = self.b, axis = self.b_axes, n = self.b_n, axes = axes, sm_spx = sm_spx, noise = self.b_noise)
            else:
                self.b, self.b_axes, self.b_n, dump = self.merge_spx_modulation_helper(data = self.b, axis = self.b_axes, n = self.b_n, axes = axes, sm_spx = sm_spx, noise = -1)


    def merge_spx_modulation_helper(self, data, axis, n, axes, sm_spx, noise = -1):

        # make the new spectrometer axis
        axis_new, axis_sort_idx = self.merge_spx_axes(axes[0], axes[1])

        # make new data structures
        _n = n[:]
        _n[0] *= 2
        _n[4] = int(_n[4]/2)
        t = numpy.result_type(data)
        _data = numpy.zeros(_n, dtype = t)             
    
        for de in range(n[5]):  
            for sm in range(n[4]):
                if sm < int(len(sm_spx)/2):
                    _data[:n[0],:,:,:, sm_spx[sm],de,:,:] = data[:,:,:,:, sm,de,:,:]
                else:
                    _data[n[0]:,:,:,:, sm_spx[sm],de,:,:] = data[:,:,:,:, sm,de,:,:]
        # reorder the pixels
        _data[:,:,:,:, :,:,:,:] = _data[axis_sort_idx,:,:,:, :,:,:,:]

        if type(noise) == numpy.ndarray:
            _noise = numpy.zeros(_n, dtype = t)
            for de in range(n[5]):  
                for sm in range(n[4]):
                    if sm < int(len(sm_spx)/2):
                        _noise[:n[0],:,:,:, sm_spx[sm],de,:,:] = noise[:,:,:,:, sm,de,:,:]
                    else:
                        _noise[n[0]:,:,:,:, sm_spx[sm],de,:,:] = noise[:,:,:,:, sm,de,:,:]
            _noise[:,:,:,:, :,:,:,:] = _noise[axis_sort_idx,:,:,:, :,:,:,:]
        else:
            _noise = -1
  
        # update the new axes
        _axis = axis[:]
        _axis[0] = axis_new[axis_sort_idx]
        _axis[4] = numpy.array([[0,1]])

        return _data, _axis, _n, _noise 

        

    def merge_spx_axes(self, axis_1, axis_2):
    
        axis_new = numpy.concatenate((axis_1, axis_2))

        axis_sort_idx = numpy.argsort(axis_new)
    
        return axis_new, axis_sort_idx        


    def average_scans(self, scans):
        
        if len(scans) > 1 and self.s_n[7] > 1:
            self.s_n[7] = 1
            self.r_n[7] = 1
            
            self.s_axes[7] = numpy.arange(1)
            self.r_axes[7] = numpy.arange(1)
            
            s = numpy.zeros(self.s_n)
            r = numpy.zeros(self.r_n)
            
            for sc in scans:        
                s[:,:,:,:, :,:,:,0] += self.s[:,:,:,:, :,:,:,sc]
                r[:,:,:,:, :,:,:,0] += self.r[:,:,:,:, :,:,:,sc]

            s /= len(scans)
            r /= len(scans)

            self.s = s
            self.r = r

            if self.measurement_method in ["Show Spectrum"]:
                s_noise = numpy.zeros(self.s_n)
                r_noise = numpy.zeros(self.r_n)

                for sc in scans:        
                    s_noise[:,:,:,:, :,:,:,0] += self.s_noise[:,:,:,:, :,:,:,sc]
                    r_noise[:,:,:,:, :,:,:,0] += self.r_noise[:,:,:,:, :,:,:,sc]

                s_noise /= len(scans)
                r_noise /= len(scans)
                
                self.s_noise = s_noise
                self.r_noise = r_noise


        elif len(scans) > 1:
            self.printWarning("No scans to average.")
            
            
    def find_indices(self, wn): 
        
        w3_i = self.find_indices_helper(self.s_axes[0], wn)
        w1_i = self.find_indices_helper(self.s_axes[1], wn)
        
        return w3_i, w1_i         
        
        
    def find_indices_helper(self, list, value):
        
        if value < list[0]:
            if (list[1] - list[0]) > (list[0] - value):
                i = 0
            else:
                i = -1            

        elif value > list[-1]:  
            if (list[-1] - list[-2]) > (value - list[-1]):
                i = len(list) - 1
            else:
                i = -1
        
        else:
            i = numpy.where(list > value)[0][0]

            if (list[i] - value) > (value - list[i-1]):
                i -= 1
                
        return i

if __name__ == "__main__": 
    pass
