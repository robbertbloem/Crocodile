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
    """
    Yet-to-be-implemented
    
    """
    

    def __init__(self, objectname, flag_verbose = 0):
        self.verbose("New VCD class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = "VCD", flag_verbose = flag_verbose)
        


class show_shots(MH.MosquitoHelperMethods):
    """
    Show shots
    
    When data is saved in Show Shots, all shots for all pixels are saved. 
    
    r: Although the data is written to different files (scans), they are imported as if it is one long scan. The shots dimension is thus (shots/scan) * (number of scans). The data is separated by datastates. Probe and reference are stored separately.  
    
    f: I call this half-signal. The probe is divided by the reference. The length is n_signals: (total number of shots) / (number of shots per signal). (number of shots per signal) is usually the same as the number of datastates, but not for the VCD. In that case the shots with the same datastate are averaged. 
    
    s: The signal. 
    
    r, f and s are saved as a numpy binary file. This saves time importing all the files all the time. When the binary file is used, it will be memory mapped -- this means it is not actually read until it is needed. This saves yet more time. 
    
    
    """

    def __init__(self, objectname, flag_verbose = 0):
        """
        Initialize show shots. 
    
        INPUT:
        - objectname (str): a name
        - flag_verbose (int, default 0): how verbose the program should be. The larger the value, the more verbose. 

        CHANGELOG:
        20160414-RB: started function
    
        """
        self.verbose("New show_shots class", flag_verbose)
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = "Show Shots", flag_verbose = flag_verbose)
        
        self.ss_colors = ["k", "r", "green", "blue", "yellow", "orange", "gold", "purple", "brown", "pink", "darkgreen", "lightblue", "grey"]

    def import_data(self, reload_data = False):
        """
        Import the show shots data
    
        INPUT:
        - reload_data (Bool, False): if False, try to import the numpy binary file first. If that is not present, or if True, it will import the original csv files. 
    
        CHANGELOG:
        201604-RB: started function
    
        """
        self.import_data_show_shots(reload_data = reload_data)


    def plot_shots(self, ax = False, pixels = 15, shots = -1):  
        """
        Make a plot of the shots.

        INPUT:
        - ax (matplotlib axes object, default: False): If False, a new figure will be made.
        - pixels (int, list or ndarray, default: 15): the pixel or pixels you want to view.
        - shots (int, list or ndarray, default: -1): the range of shots to be plotted. If shots < 0, all will be plotted. If type(shots) = int, shots 0 to int are plotted. If shots is a list [a,b], shots a to b will be plotted. 

        CHANGELOG:
        201604-RB: started function
     
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
        Plot the signals.
        
        INPUT:
        ax (matplotlib axes object, default: False): If False, a new figure will be made.
        pixels (int, list or ndarray, default: 15): the pixel or pixels you want to view.
        signals (int, list or ndarray, default: -1): the range of signals to be plotted. If signals < 0, all will be plotted. If type(signals) = int, signals 0 to int are plotted. If signals is a list [a,b], signals a to b will be plotted. 

        CHANGELOG:
        201604-RB: started function
   
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
        MH.MosquitoHelperMethods.__init__(self, objectname = objectname, measurement_method = "Pump Probe", flag_verbose = flag_verbose)


class scan_spectrum(MH.MosquitoHelperMethods):
    """
    Scan spectrum. 
    
    r: the pixel axis is now the for the wavelengths that have been scanned. Datastates is 2: probe and reference. 
    
    
    """
    
    def __init__(self, objectname, flag_verbose = 0):
        self.verbose("New scan_spectrum class", flag_verbose)
        DCC.dataclass.__init__(self, objectname = objectname, measurement_method = "Scan Spectrum", flag_verbose = flag_verbose)

    def import_data(self):
        """
        Import Scan Spectrum data.

        CHANGELOG:
        201604-RB: started function
    
        """    

        self.import_data_scan_spectrum()
        
    def make_plot(self, ax = False, normalize = False, fit = False):
        """
        Make a plot of scan spectrum data. 
    
        INPUT:
        - ax (plt axis instance, or False): If False, a new figure and axis instance will be made. 
        - normalize (bool, False): If True, the minimum is subtract from the data, then it is divided by the maximum. 
        - fit (Bool, False): If True, a fit will be made and will also be plotted. The fitting parameters are written to the terminal. 
    
        CHANGELOG:
        201604-RB: started function
    
        """
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

                      
            print("           mu       sigma   offset    scale")
            for ds in range(self.r_n[2]):
                A = [sigma, self.r_axes[0][numpy.argmax(self.r[:,0,ds,0,0,0,0,0])], 0, 1] # initial guess
            
                A_final = M.fit(self.r_axes[0], self.r[:,0,ds,0,0,0,0,0], EQ.rb_gaussian, A)
                ax.plot(self.r_axes[0], EQ.rb_gaussian(A_final, self.r_axes[0]), color = colors[ds])
                
                print("{label:10} {mu:.5}   {sigma:.3}   {offset:.3}   {scale:.3}".format(label = labels[ds], mu = A_final[1], sigma = A_final[0], offset = A_final[2], scale = A_final[3]))
                


        


        


class FT2DIR(MH.MosquitoHelperMethods):
    """
    Process 2D-IR data. 
    
    b: only when data is saved as 'intensity'. The length is n_bins. Probe and reference are not yet divided. 
    
    r: Response. The length of 'bins' is the number of bins where t1 >= 0. 'r' is for response: delta absorption has already been calculated. 
    
    f: Fourier transformed data. The length of 'bins' depends on the amount of zeropadding. These are complex numbers. 
    
    s: Spectrum. 
    
    


    """

    def __init__(self, objectname, flag_verbose = 0):
        self.verbose("New FT2DIR class", flag_verbose)
        DCC.dataclass.__init__(self, objectname = objectname, measurement_method = "FT_2DIR", flag_verbose = flag_verbose)

    def import_data(self, **kwargs): #t1_offset = 0, import_temp_scans = False):
        """
        2D-IR import data. Wrapper for the Mosuqito Helper import method. 

        INPUT:
        - t1_offset (int, 0): if t1=0 is not where it should be, adjust it. The value is in integer and is the number of bins.
        - import_temp_scans (Bool, False): If True, import individual scans instead of the averaged one.

        CHANGELOG:
        201604-RB: started function
    
        """    

        self.import_data_2dir(**kwargs) #import_temp_scans = import_temp_scans, t1_offset = t1_offset)

    def make_plots(self, **kwargs):
        """
        Plot all 2D-IR spectra as separate plots. 
    
        INPUT:
        - kwargs:
            - sp, sm, de, du, sc (lists, ndarray): lists that limit the range of spectra to be plotted. By default all spectra will be plotted. 
            - aspect (str, 'equal'): the aspect ratio of the axes. 
            - flip_spectrum (bool, False): by default (False) the x axis is w_3. If True, it will be w_1. Make sure the axes ranges are changed as well. 
    
        DESCRIPTION:
        - - 
    
        CHANGELOG:
        201604-RB: started function
    
        """
        
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)
        
        flag_make_title = False
        if "title" not in kwargs:
            flag_make_title = True
        
        for _sp in sp:
            for _sm in sm:
                for _de in de:
                    for _du in du:
                        for _sc in sc:

                            fig = plt.figure()
                            ax = fig.add_subplot(111)  
                            if "aspect" in kwargs:
                                if kwargs["aspect"] != False:
                                    ax.set_aspect(kwargs["aspect"])
                            else:
                                ax.set_aspect("equal")
                            
                            if flag_make_title:
                                kwargs["title"] = "%s %s fs" % (self._basename, self.s_axes[5][_de])

                            if "flip_spectrum" in kwargs and kwargs["flip_spectrum"]:
                                PL.contourplot(self.s[:, :, 0, _sp, _sm, _de, _du, _sc], self.s_axes[1], self.s_axes[0], x_label = "w1 (cm-1)", y_label = "w3 (cm-1)", ax = ax, **kwargs)
                     
                            else:
                                PL.contourplot(self.s[:, :, 0, _sp, _sm, _de, _du, _sc].T, self.s_axes[0], self.s_axes[1], x_label = "w3 (cm-1)", y_label = "w1 (cm-1)", ax = ax, **kwargs)


        plt.show()
   

    def make_overview_plot(self, **kwargs):
        """
        Make a single figure with all plots. 
        
        INPUT:
        - kwargs:
            - sp, sm, de, du, sc (lists, ndarray): lists that limit the range of spectra to be plotted. By default all spectra will be plotted.
            - aspect (str, 'equal'): the aspect ratio of the axes. 
            - flip_spectrum (bool, False): by default (False) the x axis is w_3. If True, it will be w_1. Make sure the axes ranges are changed as well. 
        
        OUTPUT:
        - -
        
        DESCRIPTION:
        - - 
        
        CHANGELOG:
        2016014/RB: started function
        
        """
    
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)
    
        n_plots = len(sp) * len(sm) * len(de) * len(du) * len(sc)

        x, y = FU.find_subplots(n_plots, flag_verbose = self.flag_verbose)
        print(x,y)
        fig = plt.figure()
        ax = [0] * n_plots
        for ax_i in range(n_plots):
            ax[ax_i] = fig.add_subplot(y, x, ax_i + 1)  
            if "aspect" in kwargs:
                if kwargs["aspect"] != False:
                    ax[ax_i].set_aspect(kwargs["aspect"])
            else:
                ax[ax_i].set_aspect("equal")

        ax_i = 0
        for _sp in sp:
            for _sm in sm:
                for _de in de:
                    for _du in du:
                        for _sc in sc:

                            title = "%s %s fs" % (self._basename, self.s_axes[5][_de])
                            if "flip_spectrum" in kwargs and kwargs["flip_spectrum"]:
                                PL.contourplot(self.s[:, :, 0, _sp, _sm, _de, _du, _sc], self.s_axes[1], self.s_axes[0], x_label = "w1 (cm-1)", y_label = "w3 (cm-1)", ax = ax[ax_i], title = title, **kwargs)
                     
                            else:
                                PL.contourplot(self.s[:, :, 0, _sp, _sm, _de, _du, _sc].T, self.s_axes[0], self.s_axes[1], x_label = "w3 (cm-1)", y_label = "w1 (cm-1)", ax = ax[ax_i], title = title, **kwargs)

                            ax_i += 1
                            

    def make_Z_table(self, x_range = [0,0], y_range = [0,-1], **kwargs):
        """
        
        INPUT:
        - -
        
        OUTPUT:
        - -
        
        DESCRIPTION:
        - - 
        
        CHANGELOG:
        2016014/RB: started function
        
        """
        
        print("sp sm de sc   min    max    delta")

        # determine the range to be plotted
        x_min, x_max, y_min, y_max = FU.find_axes(self.s_axes[0], self.s_axes[1], x_range = x_range, y_range = y_range, flag_verbose = self.flag_verbose)
        
#         print(x_min, x_max, y_min, y_max)
    
        # find the area to be plotted
        x_min_i, x_max_i = FU.find_axes_indices(self.s_axes[0], x_min, x_max)
        y_min_i, y_max_i = FU.find_axes_indices(self.s_axes[1], y_min, y_max)
    
#         print(x_min_i, x_max_i, y_min_i, y_max_i)
        
        for _sp in range(self.s_n[2]):
            for _sm in range(self.s_n[4]):
                for _de in range(self.s_n[5]):
                    for _sc in range(self.s_n[7]):
                        data, x_axis, y_axis = FU.truncate_data(self.s[:,:,_sp,0,_sm,_de,0,_sc].T, self.s_axes[0], self.s_axes[1], x_min_i, x_max_i, y_min_i, y_max_i)
                        
#                         print(data)
                    
                        a = numpy.amin(data)
                        b = numpy.amax(data)
                        print(" {sp}  {sm}  {de}  {sc}  {min:5.1f} {max:5.1f} {delta:5.1f}".format(sp = _sp, sm = _sm, de = _de, sc = _sc, min = a, max = b, delta = (b - a)))
        
        
    def export_2d_data_for_gnuplot(self, s = True, **kwargs):
        
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)
    
        self.check_results_folder()
    
        paf = self.file_dict["result_folder"] + self._file_dict["basename"] + "_" + self._file_dict["timestamp"]
    
        ax_i = 0
        for _sp in sp:
            for _sm in sm:
                for _de in de:
                    for _du in du:
                        for _sc in sc:
                            
                            if s: 
                                
                                x_axis = self.s_axes[0]
                                y_axis = self.s_axes[1]
                                data = self.s[:,:, 0, _sp, _sm, _de, _du, _sc]
                                path_and_filename =  "{paf}_s_gnu_sp{sp}_sm{sm}_de{de}_du{du}_sc{sc}{ext}".format(paf = paf , sp = _sp, sm = _sm, de = _de, du = _du, sc = _sc, ext = self.file_dict["extension"])  
                            
                                self.export_2d_data_for_gnuplot_helper(x_axis, y_axis, data, path_and_filename)        


if __name__ == "__main__": 


    pass