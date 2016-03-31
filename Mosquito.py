from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.DataClassCol as DCC
import Crocodile.Resources.MosquitoHelper as MH
import Crocodile.Resources.IOMethods as IOM
import Crocodile.Resources.Functions as FU

import imp

imp.reload(MH)
imp.reload(IOM)

class VCD(MH.MosquitoHelperMethods):

    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New VCD class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = DCC.MeasurementMethod["vcd"], flag_verbose = flag_verbose)
        


class show_shots(MH.MosquitoHelperMethods):
    """
    show_shots
    """

    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New show_shots class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = DCC.MeasurementMethod["show_shots"], flag_verbose = flag_verbose)
        
        self.ss_colors = ["k", "r", "green", "blue", "yellow", "orange", "gold", "purple", "brown", "pink", "darkgreen", "lightblue", "grey"]


    def make_plot(self, ax = False, pixels = 15, scans = -1, concat = True):  
        """
        
        INPUT:
        ax (matplotlib axes object, default: False): If False, a new figure will be made.
        pixels (int, list or ndarray, default: 15): the pixel or pixels you want to view.
        scans (int, list or ndarray, default: -1): the scan or scans you want to view. If <0, all scans will be shown. 
        
        """

        shot_axes = numpy.arange(self.b_n[1])
        
        if scans < 0:
            scan_range = numpy.arange(self.b_n[7])
        elif type(scans) == "int":
            scan_range = numpy.array[scans]
        elif type(scans) == numpy.ndarray:
            scan_range = scans
        elif type(scans) == "list":
            scan_range = numpy.array(scans)
        else:
            scan_range = numpy.arange(self.b_n[7])
        
        if type(pixels) == "int":
            pixel_range = numpy.array([pixels])
        elif type(pixels) == numpy.ndarray:
            pixel_range = pixels
        elif type(pixels) == "list":
            pixel_range = numpy.array(pixels)
        else:
            pixel_range = numpy.array([15])
        
        ax, new_axis = self.check_axis(ax)
            
#         if ax == False:
#             fig = plt.figure()
#             ax = fig.add_subplot(111)    
        
        
        
        for ds in range(self.b_n[2]):
            color = self.ss_colors[ds]
            for pi in pixel_range:
                for sc in scan_range:
                    data = self.b[pi,:,ds,0,0,0,0,sc]
                    mask = numpy.isfinite(data)   
                
                    if concat:
                        temp = sc * self.b_n[1]
                    else:
                        temp = 0

                    ax.plot(shot_axes[mask] + temp, data[mask], marker = ".", linestyle = "none", color = color)
            
        plt.show()


    def plot_signal(self, ax = False):
        
        pixels = numpy.arange(32)
        
        shots_per_signal = 6
        
        total_shots = self.b_n[1] * self.b_n[7]
        
        total_signals = int(total_shots / shots_per_signal)
        signals = numpy.zeros((32, total_signals))
        
        half_signals = numpy.zeros((32, 4, total_signals))
        
        signals_per_scan = int(self.b_n[1] / shots_per_signal)
        
        pi = 12
        
        ax, new_axis = self.check_axis(ax)
        
        sig_counter = 0
        for sc in range(self.b_n[7]):   
            for sh in range(signals_per_scan):
                sh_s = sh * shots_per_signal
                sh_e = sh_s + shots_per_signal
                
                if sc == 0 and sh == 0:
                    print(self.b[pi,sh_s:sh_e, :, 0,0,0,0,sc])
                
                for pi in pixels:
                    L0p = numpy.nanmean(self.b[pi,sh_s:sh_e, 0, 0,0,0,0,sc])
                    L0r = numpy.nanmean(self.b[pi,sh_s:sh_e, 1, 0,0,0,0,sc])
                    R0p = numpy.nanmean(self.b[pi,sh_s:sh_e, 2, 0,0,0,0,sc])
                    R0r = numpy.nanmean(self.b[pi,sh_s:sh_e, 3, 0,0,0,0,sc])
                    L1p = numpy.nanmean(self.b[pi,sh_s:sh_e, 4, 0,0,0,0,sc])
                    L1r = numpy.nanmean(self.b[pi,sh_s:sh_e, 5, 0,0,0,0,sc])
                    R1p = numpy.nanmean(self.b[pi,sh_s:sh_e, 6, 0,0,0,0,sc])
                    R1r = numpy.nanmean(self.b[pi,sh_s:sh_e, 7, 0,0,0,0,sc])
                    
                    L0 = L0p / L0r
                    R0 = R0p / R0r
                    L1 = L1p / L1r
                    R1 = R1p / R1r
                    
                    half_signals[pi, :, sig_counter] = numpy.array([L0, R0, L1, R1])
                    
                    signals[pi, sig_counter] = -numpy.log10((L1 * R0) / (R1 * L0))

                    
                sig_counter += 1
            

class pump_probe(MH.MosquitoHelperMethods):
    """
    pump probe
    """

    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New pump_probe class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = DCC.MeasurementMethod["pump_probe"], flag_verbose = flag_verbose)



# class show_spectrum(MH.MosquitoHelperMethods):
#     """
#     show_spectrum
#     """
# 
#     def __init__(self, objectname, flag_verbose = False):
#         self.verbose("New show_spectrum class", flag_verbose)
#         MH.MosquitoHelperMethods.__init__(self, objectname = objectname, flag_verbose = flag_verbose)
# 
# 
# #     def import_data(self):
# #         """
# #         This method splits up file importing for supporting files and measurement files. It first finds the file format file. It will use this to test if file_dict is correct. 
# #         """
# #         # check if _file_dict is set
# #         if self._file_dict["base_folder"] == "" or self._file_dict["base_filename"] == "":
# #             self.printError("No file information set.", inspect.stack()) 
# #             return False   
# #         
# #         # find LV_file_format, also a check if _file_dict is correct
# #         try: 
# #             self.file_format = IOM.find_LV_fileformat(
# #                 base_folder = self._file_dict["base_folder"], 
# #                 flag_verbose = self.flag_verbose
# #             )   
# #             if self.file_format == -1:  
# #                 self.printError("File_format file found, but was not able to parse it.", inspect.stack()) 
# #                 return False  
# #         except FileNotFoundError:
# #             self.printError("File_format file not found.", inspect.stack()) 
# #             return False 
# # 
# #         w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
# # 
# # #         n_sh = IOM.import_nshots(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
# # #         sh = numpy.arange(n_sh)
# #         
# #         n_sh = 1
# #         sh = numpy.array([0])
# #         
# #         n_ds = IOM.find_number_of_datastates(self._file_dict["base_folder"], flag_verbose = self.flag_verbose)
# #         n_ds = int(n_ds/2)
# #         
# #         n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
# # 
# #         spds, dump, dump = IOM.import_spectraAndDatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
# #         
# #         self.sp = spds[:,0]
# #         self.ds = spds[:,1]
# #         for i in range(n_ds):
# #             if self.ds[i] == "-1":
# #                 self.ds[i] = 1
# #             else:
# #                 self.ds[i] = 0
# #         self.ds = numpy.array(self.ds, dtype = "int")
# # 
# #         n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
# #         sc = numpy.arange(n_sc)
# #         
# #         self.b_n = [n_w3, n_sh, 2*n_ds, n_sp, n_sc]
# #         self.b = numpy.empty(self.b_n)
# #         self.b_noise = numpy.empty(self.b_n)
# #         self.b_axes = [w3_axis, sh, self.ds, self.sp, sc]
# #         self.b_units = ["w3 (cm-1)", "Shots", "Datastates", "Spectra", "Scans"]
# #         
# #         self.s_n = numpy.copy(self.b_n)
# #         self.s_n[2] = 1
# #         self.s = numpy.zeros(self.s_n)
# #         self.s_noise = numpy.zeros(self.s_n)
# # 
# #         sh = 0        
# #         for sp in range(self.b_n[3]): 
# #             for sc in range(self.b_n[4]): 
# #                 
# #                 ds = 0
# #                 suffix = "sp%i_ds%i_signal_%i" % (sp, ds, sc)
# #                 self.s[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose) 
# #                 suffix = "sp%i_ds%i_signal_noise_%i" % (sp, ds, sc)
# #                 self.s_noise[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)
# #               
# #                 for ds in range(self.b_n[2]): 
# #                     suffix = "sp%i_ds%i_intensity_%i" % (sp, ds, sc)
# #                     self.b[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose) 
# #                     suffix = "sp%i_ds%i_intensity_noise_%i" % (sp, ds, sc)
# #                     self.b_noise[:,sh,ds,sp,sc] = IOM.import_file(self._file_dict, suffix, self.flag_verbose)
# 
# 
# 
#     def make_plot(self):  
# 
#         axes = numpy.arange(self.b_n[1])
#         
#         fig = plt.figure()
#         ax = fig.add_subplot(111)    
#         
#         for ds in range(self.b_n[2]):
#             for sc in range(self.b_n[4]):
#                 data = self.b[15,:,ds,0,sc]
#                 mask = numpy.isfinite(data)            
#                 ax.plot(axes[mask], data[mask], marker = ".", linestyle = "none")
#             
#         plt.show()
# 
# 
# 
# class find_t0_fast(MH.MosquitoHelperMethods):
#     """
#     find_t0_fast
#     """
# 
#     def __init__(self, objectname, flag_verbose = False):
#         self.verbose("New find_t0_fast class", flag_verbose)
#         MH.MosquitoHelperMethods.__init__(self, objectname = objectname, flag_verbose = flag_verbose)
# 
# 
# #     def import_data(self):
# #         """
# #         This method splits up file importing for supporting files and measurement files. It first finds the file format file. It will use this to test if file_dict is correct. 
# #         """
# #         # check if _file_dict is set
# #         if self._file_dict["base_folder"] == "" or self._file_dict["base_filename"] == "":
# #             self.printError("No file information set.", inspect.stack()) 
# #             return False   
# #         
# #         # find LV_file_format, also a check if _file_dict is correct
# #         try: 
# #             self.file_format = IOM.find_LV_fileformat(
# #                 base_folder = self._file_dict["base_folder"], 
# #                 flag_verbose = self.flag_verbose
# #             )   
# #             if self.file_format == -1:  
# #                 self.printError("File_format file found, but was not able to parse it.", inspect.stack()) 
# #                 return False  
# #         except FileNotFoundError:
# #             self.printError("File_format file not found.", inspect.stack()) 
# #             return False 
# # 
# #         t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
# #         
# #         self.bin_sign = bin_sign
# #         self.t1_zero_index = t1_zero_index
# # 
# #         n_ds = IOM.find_number_of_datastates(self._file_dict["base_folder"], flag_verbose = self.flag_verbose)
# #         
# #         n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
# # 
# #         self.ds = spds[:,1]
# #         for i in range(n_ds):
# #             if self.ds[i] == "-1":
# #                 self.ds[i] = 1
# #             else:
# #                 self.ds[i] = 0
# #         self.ds = numpy.array(self.ds, dtype = "int")
# #         
# #         self.b_intf_n = [n_t1_bins, n_ds, n_sp]
# #         
# #         self.b_intf = numpy.empty(self.b_intf_n)
# #         self.b_intf_axes = [t1_bins, spds[:,1], self.ds]
# #         self.b_intf_units = ["T1 (bins)", "Datastates", "Spectra"]
# 
# 
# 
# class FT2DIR(MH.MosquitoHelperMethods):
#     """
# 
#     """
# 
#     def __init__(self, objectname, flag_verbose = False):
#         self.verbose("New pe_col class", flag_verbose)
#         DCC.dataclass.__init__(self, objectname = objectname, flag_verbose = flag_verbose)
# 
# 
# 
#     def make_plots(self, x_range = [0,0], invert_colors = False, flip_spectrum = False):
#  
#         inv = 1
#         if invert_colors:
#             inv = -1
# 
# 
#         for sp in range(self.s_n[3]):
#             for sm in range(self.s_n[4]):
#                 for de in range(self.s_n[5]): 
#                     for du in range(self.s_n[6]):
#                         for sc in range(self.s_n[7]):
# 
#                             fig = plt.figure()
#                             ax = fig.add_subplot(111)        
#                             ax.set_aspect("equal")
#                         
#                             title = "%s %s fs" % (self._basename, self.s_axes[5][de])
#                             if flip_spectrum:
#                                 PL.contourplot(inv * self.s[:, :, 0, sp, sm, de, du, sc], self.s_axes[1], self.s_axes[0], y_range = [0,0], x_range = x_range, y_label = "w3 (cm-1)", x_label = "w1 (cm-1)", title = title, ax = ax)                        
#                             else:
#                                 PL.contourplot(inv * self.s[:, :, 0, sp, sm, de, du, sc].T, self.s_axes[0], self.s_axes[1], x_range = x_range, y_range = [0,-1], x_label = "w3 (cm-1)", y_label = "w1 (cm-1)", title = title, ax = ax)
# 
#         plt.show()
#    



if __name__ == "__main__": 


    pass