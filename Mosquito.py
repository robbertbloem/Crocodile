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
import Crocodile.Resources.Plotting as PL
import Crocodile.Resources.Equations as EQ
import Crocodile.Resources.Mathematics as M

import imp

imp.reload(MH)
imp.reload(IOM)
imp.reload(PL)

class VCD(MH.MosquitoHelperMethods):

    def __init__(self, objectname, flag_verbose = 0):
        self.verbose("New VCD class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = DCC.MeasurementMethod["vcd"], flag_verbose = flag_verbose)
        


class show_shots(MH.MosquitoHelperMethods):
    """
    show_shots
    
    r: the individual shots, separated by their datastate
    f: average shots for a datastate. I.e. in the VCD experiment, some datastates are measured twice and others only once. 
    s: the signals. 
    
    NORMAL CHOPPED MEASUREMENT  
    100 shots
    datastates = [0,1]
    r: 4 datastates, 50 shots for each
    f: as r
    s: 50 signals
    
    VCD MEASUREMENT
    120 shots
    datastates = [0,1,0,2,3,2]
    r: 4 datastates, 40 shots for each
    f: 4 datastates with either 40 or 80 shots
    s: 40 signals
    
    
    
    """

    def __init__(self, objectname, flag_verbose = 0):
        self.verbose("New show_shots class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = DCC.MeasurementMethod["show_shots"], flag_verbose = flag_verbose)
        
        self.ss_colors = ["k", "r", "green", "blue", "yellow", "orange", "gold", "purple", "brown", "pink", "darkgreen", "lightblue", "grey"]

    def import_data(self, reload_data = False):
        self.import_data_show_shots(reload_data = reload_data)


    def plot_shots(self, ax = False, pixels = 15, shots = -1):  
        """
        
        INPUT:
        ax (matplotlib axes object, default: False): If False, a new figure will be made.
        pixels (int, list or ndarray, default: 15): the pixel or pixels you want to view.
        shots (int, list or ndarray, default: -1): the range of shots to be plotted. If shots < 0, all will be plotted. If type(shots) = int, shots 0 to int are plotted. If shots is a list [a,b], shots a to b will be plotted. 
        
        """

        if type(shots) == int:
            if shots < 0:
                s = 0
                e = self.r_n[1]
            else:
                s = 0
                e = shots
        elif type(shots) == numpy.ndarray or type(shots) == list:
            s = shots[0]
            e = shots[1]
        else:
            s = 0
            e = self.r_n[1]
        
        shot_axes = numpy.arange(self.r_n[1])[s:e]
        
        if type(pixels) == "int":
            pixel_range = numpy.array([pixels])
        elif type(pixels) == numpy.ndarray:
            pixel_range = pixels
        elif type(pixels) == "list":
            pixel_range = numpy.array(pixels)
        else:
            pixel_range = numpy.array([15])
        
        ax, new_axis = self.check_axis(ax)
            
        if ax == False:
            fig = plt.figure()
            ax = fig.add_subplot(111)    
        
        for ds in range(self.r_n[2]):
            color = self.ss_colors[ds]
            for pi in pixel_range:

                data = self.r[pi,s:e,ds,0,0,0,0,0]
                mask = numpy.isfinite(data)   

                ax.plot(shot_axes[mask], data[mask], marker = ".", linestyle = "none", color = color)
                
        ax.set_xlabel("Shots")
        ax.set_ylabel("V")
            
        plt.show()


    def plot_signals(self, ax = False, pixels = 15, range = -1):  
        """
        
        INPUT:
        ax (matplotlib axes object, default: False): If False, a new figure will be made.
        pixels (int, list or ndarray, default: 15): the pixel or pixels you want to view.
        signals (int, list or ndarray, default: -1): the range of signals to be plotted. If signals < 0, all will be plotted. If type(signals) = int, signals 0 to int are plotted. If signals is a list [a,b], signals a to b will be plotted. 
        
        """

        if type(range) == int:
            if range < 0:
                s = 0
                e = self.s_n[1]
            else:
                s = 0
                e = range
        elif type(range) == numpy.ndarray or type(range) == list:
            s = range[0]
            e = range[1]
        else:
            s = 0
            e = self.s_n[1]
        
        shot_axes = numpy.arange(self.s_n[1])[s:e]
        
        if type(pixels) == "int":
            pixel_range = numpy.array([pixels])
        elif type(pixels) == numpy.ndarray:
            pixel_range = pixels
        elif type(pixels) == "list":
            pixel_range = numpy.array(pixels)
        else:
            pixel_range = numpy.array([15])
        
        ax, new_axis = self.check_axis(ax)
            
        if ax == False:
            fig = plt.figure()
            ax = fig.add_subplot(111)    

        for pi in pixel_range:

            data = self.s[pi,s:e,0,0,0,0,0,0]
            mask = numpy.isfinite(data)   

            ax.plot(shot_axes[mask], data[mask], marker = ".", linestyle = "none", color = "b")

        ax.set_xlabel("Signals")
        ax.set_ylabel("mOD")
                    
        plt.show()


#     def print_stats(self):
#         mean_s = numpy.mean(self.s)
#         std_s = numpy.std(self.s)
#         
#         print(mean_s, std_s)

            

class pump_probe(MH.MosquitoHelperMethods):
    """
    pump probe
    """

    def __init__(self, objectname, flag_verbose = 0):
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


class scan_spectrum(MH.MosquitoHelperMethods):
    """
    
    """
    
    def __init__(self, objectname, flag_verbose = 0):
        self.verbose("New scan_spectrum class", flag_verbose)
        DCC.dataclass.__init__(self, objectname = objectname, measurement_method = DCC.MeasurementMethod["scan_spectrum"], flag_verbose = flag_verbose)

    def import_data(self):
        self.import_data_scan_spectrum()
        
    def make_plot(self, ax = False, normalize = False, fit = False):
        
        if ax == False:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        
        if normalize:
            for ds in range(self.r_n[2]):
                self.r[:,0,ds,0,0,0,0,0] -= numpy.nanmin(self.r[:,0,ds,0,0,0,0,0])
                self.r[:,0,ds,0,0,0,0,0] /= numpy.nanmax(self.r[:,0,ds,0,0,0,0,0])

        ax.plot(self.r_axes[0], self.r[:,0,0,0,0,0,0,0], color = "g")
        ax.plot(self.r_axes[0], self.r[:,0,1,0,0,0,0,0], color = "r")
             
             
        
        if fit:
            colors = ["lightgreen", "orange"]
            labels = ["probe", "reference"]
            sigma = (self.r_axes[0][0] - self.r_axes[0][-1]) / 4
            A = [sigma,numpy.mean(self.r_axes[0]), 0, 1] # initial guess
            
            
            print("           mu       sigma   offset    scale")
            for ds in range(self.r_n[2]):
                A_final = M.fit(self.r_axes[0], self.r[:,0,ds,0,0,0,0,0], EQ.rb_gaussian, A)
                ax.plot(self.r_axes[0], EQ.rb_gaussian(A_final, self.r_axes[0]), color = colors[ds])
                
                print("{label:10} {mu:.5}   {sigma:.3}   {offset:.3}   {scale:.3}".format(label = labels[ds], mu = A_final[1], sigma = A_final[0], offset = A_final[2], scale = A_final[3]))
                
#                 print("mu: {mu}, sigma: {sigma}, offset: {os}, scale: {sc}".format())

        


        


class FT2DIR(MH.MosquitoHelperMethods):
    """

    """

    def __init__(self, objectname, flag_verbose = 0):
        self.verbose("New FT2DIR class", flag_verbose)
        DCC.dataclass.__init__(self, objectname = objectname, measurement_method = DCC.MeasurementMethod["ft_2d_ir"], flag_verbose = flag_verbose)

    def import_data(self, t1_offset = 0, import_temp_scans = False):
        self.import_data_2dir(import_temp_scans = import_temp_scans, t1_offset = t1_offset)

    def make_plots(self, aspect = "equal", single_plot = False, **kwargs):
            
        if single_plot:
            pass
            
        for sp in range(self.s_n[3]):
            for sm in range(self.s_n[4]):
                for de in range(self.s_n[5]): 
                    for du in range(self.s_n[6]):
                        for sc in range(self.s_n[7]):
                            
#                             if single_plot == False:
                            fig = plt.figure()
                            ax = fig.add_subplot(111)  
                            if aspect:
                                ax.set_aspect(aspect)
                        
                            title = "%s %s fs" % (self._basename, self.s_axes[5][de])
                            if "flip_spectrum" in kwargs and kwargs["flip_spectrum"]:
                                PL.contourplot(self.s[:, :, 0, sp, sm, de, du, sc], self.s_axes[1], self.s_axes[0], x_label = "w1 (cm-1)", y_label = "w3 (cm-1)", ax = ax, **kwargs)
                     
                            else:
                                PL.contourplot(self.s[:, :, 0, sp, sm, de, du, sc].T, self.s_axes[0], self.s_axes[1], x_label = "w3 (cm-1)", y_label = "w1 (cm-1)", ax = ax, **kwargs)


        plt.show()
   





if __name__ == "__main__": 


    pass