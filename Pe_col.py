"""
Handling of 2D-IR measurements

"""


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

class pe_col(DCC.dataclass):
    """
    pe_col: implementation class for colinear 2D-IR setup. This includes all importing functions.
    """

    def __init__(self, objectname, flag_verbose = False):
        self.verbose("New pe_col class", flag_verbose)
        DCC.dataclass.__init__(self, objectname = objectname, flag_verbose = flag_verbose)

# 
#     def set_file_info(self, data_folder, date, basename, timestamp):
#         self.set_file_dict(data_folder, date, basename, timestamp)    


    def import_data(self, import_temp_scans = False):
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
        
        # import data
        self.import_supporting_data(import_temp_scans)
        self.import_measurement_data()
        return True  


    def import_supporting_data(self, import_temp_scans = False):
        
        w3_axis, n_w3 = IOM.import_wavenumbers(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index = IOM.import_bins(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        self.bin_sign = bin_sign
        self.t1_zero_index = t1_zero_index
        
        n_ds = IOM.import_ndatastates(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        n_sp = IOM.import_nspectra(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)

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
        
        sm, sm_names, n_sm = IOM.import_slow_modulation(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        de, n_de = IOM.import_delays(self._file_dict, self.file_format, flag_verbose = self.flag_verbose)
        
        du = numpy.array([0])
        n_du = 1
        
        if import_temp_scans:
            n_sc = IOM.find_number_of_scans(self._file_dict["base_folder"], self._file_dict["base_filename"], self._file_dict["extension"], flag_verbose = self.flag_verbose)
        else:
            n_sc = 1
        sc = numpy.arange(n_sc)

        if self.measurement_type == "signal":
            # data already divided by count
            # 1 data state: the signal
            _n_ds = 1
        else:
            # data not yet divided by count
            # 2*n_datastates
            _n_ds = 2 * n_ds
        
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

#         if False:
#             x = (180 - 47) * numpy.pi / 180 #####
#             f = f * numpy.exp(1j*x)#####
# 
#         if False:
#             x = (90 - 47) * numpy.pi / 180 #####
#             f = f * numpy.exp(1j*x)#####
# 
#         if False:
#             x = (270 - 47) * numpy.pi / 180 #####
#             f = f * numpy.exp(1j*x)#####
#             
#         if False:
#             x = (360 - 47) * numpy.pi / 180 #####
#             f = f * numpy.exp(1j*x)#####

        angle = numpy.angle(f)
        check_angle = numpy.angle(f * numpy.exp(1j*numpy.pi/2))

        if numpy.all(angle[idx] > -3*numpy.pi/4) and numpy.all(angle[idx] < 3*numpy.pi/4):
            self.phase_rad = numpy.mean(angle[idx])
        else:
            self.phase_rad = numpy.mean(check_angle[idx]) - numpy.pi/2
        
#         print(angle[idx] * 180 / numpy.pi)
#         print(check_angle[idx] * 180 / numpy.pi)

        print("Phase in degrees: %.1f" % (self.phase_rad * 180 / numpy.pi))

#         plt.plot(w_axis, angle)
#         plt.plot(w_axis[idx], angle[idx], lw = 2)
#         plt.show()



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
                     


    def make_plots(self, x_range = [0,0], invert_colors = False, flip_spectrum = False):
 
        inv = 1
        if invert_colors:
            inv = -1


        for sp in range(self.s_n[3]):
            for sm in range(self.s_n[4]):
                for de in range(self.s_n[5]): 
                    for du in range(self.s_n[6]):
                        for sc in range(self.s_n[7]):

                            fig = plt.figure()
                            ax = fig.add_subplot(111)        
                            ax.set_aspect("equal")
                        
                            title = "%s %s fs" % (self._basename, self.s_axes[5][de])
                            if flip_spectrum:
                                PL.contourplot(inv * self.s[:, :, 0, sp, sm, de, du, sc], self.s_axes[1], self.s_axes[0], y_range = [0,0], x_range = x_range, y_label = "w3 (cm-1)", x_label = "w1 (cm-1)", title = title, ax = ax)                        
                            else:
                                PL.contourplot(inv * self.s[:, :, 0, sp, sm, de, du, sc].T, self.s_axes[0], self.s_axes[1], x_range = x_range, y_range = [0,-1], x_label = "w3 (cm-1)", y_label = "w1 (cm-1)", title = title, ax = ax)

        plt.show()












