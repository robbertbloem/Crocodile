from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.DataClass as DC
import Crocodile.Resources.Plotting as PL
import Crocodile.Resources.Mathematics as MATH
import Crocodile.Resources.IOMethods as IOM

# import imp
# imp.reload(DC)

class pe(DC.dataclass):
    
    def __init__(self, objectname, measurements = 2, flag_verbose = False):
        
        self.verbose("New Pe.pe class", flag_verbose)
        
        DC.dataclass.__init__(self, objectname, measurements = measurements, dimensions = 3, flag_verbose = flag_verbose)



    def absorptive_helper(self, array, axes, window_function = "none", window_length = 0):
        """
        helper function for absorptive

        This is a function to Fourier Transform experimental 2DIR spectra, ie, spectra with a time and frequency axis. It basically repeats the Fourier function for all pixels.

        INPUT:
        - array (numpy.ndarray): a 2d array of time * pixels
        - window_function (name, "none"): "none", "guassian" etc
        - window_length (int, 0): length of the window. 0 means the whole range
        - flag_plot (BOOL, False): will plot the windowed time domain

        CHANGELOG:
        20101204/RB: started as fourier_helper
        20110909/RB: continued
        20130131/RB: re-implemented as absorptive-helper. The FFT is now done in this function instead of in the Mathematics module. That function had a load of overhead which was not needed.
        
        """
        
        x, y = numpy.shape(array)

        if axes == 2 or axes == [1,1]: 
            array[:,0] /= 2
            array[0,:] /= 2

            # prepare output
            # it needs to have the correct size and it should be complex
            if self.zeropad_to != None:
                x = self.zeropad_to
            ft_array = numpy.reshape(numpy.zeros(x*y, dtype=numpy.cfloat),(x,y))          
            
            # do fft
            if window_function != "none":
                self.printWarning("No 2D window functions have been implemented", inspect.stack())

            ft_array[:,:] = numpy.fft.fft2(array, n = self.zeropad_to)
        
        else:

            if axes == 1 or axes == [0,1]:                  
                array = array.T
    
            # prepare output
            # it needs to have the correct size and it should be complex
            if self.zeropad_to != None:
                x = self.zeropad_to
            ft_array = numpy.reshape(numpy.zeros(x*y, dtype=numpy.cfloat),(x,y))          
            
            # do fft
            for i in range(y):   
                array[:,i] /= 2     # correct first element
                if window_function != "none":
                    array[:,i] = MATH.window_functions(array[:,i], window_function, window_length)
                ft_array[:,i] = numpy.fft.fft(array[:,i], n = self.zeropad_to)
        
            if axes == 1 or axes == [0,1]:                  
                ft_array = ft_array.T                
        
        return ft_array







    def super_absorptive(self, axes, window_function = "none", window_length = 0, flag_verbose = False):
        """
        Calculate the absorptive spectrum.
        
        This function does the Fourier transform. It phases the spectrum. It makes the axes.
        
        INPUT:
        - axes (number or list): 0 or [1,0] for first axis, 1 or [0,1] for second axis, 2 or [1,1] for both axes
        - window_function (name, "none"): "none", "guassian" etc
        - window_length (int, 0): length of the window. 0 means the whole range
        - flag_plot (BOOL, False): will plot the windowed time domain
        
        CHANGELOG:
        20101204/RB: started
        20110909/RB: continued 
        20130131/RB: copied from croc to Crocodile. The function should now work for time-freq, freq-time and time-time domain. Only the FFT will give an error, the other problems give a warning.
        20130208/RB: now checks used variables before the calculation instead of try/except during the calculation
        
        WARNING:
        Realistically, only the time-freq FFT is performed. The others may have problems.
        
        """
        
        self.verbose("Super absorptive", flag_verbose)


        # VARIABLE CHECKING
        
        # r should be a list with numpy.ndarrays
        # if not the case, return False
        if type(self.r) != list:
            self.printError("r should be a list.", inspect.stack())  
            return False
        elif type(self.r[0]) != numpy.ndarray or type(self.r[1]) != numpy.ndarray:
            self.printError("r[0] or r[1] should be a numpy.ndarray. Did you assign r?", inspect.stack())  
            return False
        
        # r_axis should be a list
        # index 0 and 2 should be a numpy.ndarray
        # if not the case, skip making the axes
        flag_axis = True
        if type(self.r_axis) != list:
            flag_axis = False
            self.printWarning("r_axis should be a list. s_axis will not be calculated.", inspect.stack())
        elif type(self.r_axis[0]) != numpy.ndarray or type(self.r_axis[2]) != numpy.ndarray:
            flag_axis = False
            self.printWarning("r_axis[0] or r_axis[2] should be a numpy.ndarray. s_axis will not be calculated.", inspect.stack()) 
        
        # undersampling has to be an integer
        # if it is a BOOL, it can be converted to 0 or 1
        if type(self.undersampling) == bool: 
            if self.undersampling:
                self.undersampling = 1
            else:
                self.undersampling = 0
            self.printWarning("undersampling is a BOOL, will use " + str(self.undersampling), inspect.stack())
        if type(self.undersampling) != int:
            self.printError("undersampling is NaN, not a valid value. ", inspect.stack())
            return False
        
        # the phase should not be a numpy.nan
        if numpy.isnan(self.phase_degrees):
            self.printError("phase_degrees is NaN, not a valid value. ", inspect.stack())    
            return False            

        # ACTUAL FUNCTION

        # do the FFT
        try: 
            for i in range(len(self.r)):
                self.f[i] = self.absorptive_helper(array = numpy.copy(self.r[i]), axes = axes, window_function = window_function, window_length = window_length)
        except ValueError:
            self.printError("Problem with the Fourier Transforms. Are r[0] and r[1] assigned?", inspect.stack())
            return False          
        
        self.verbose("  done with FFT", flag_verbose)
        
        # phase the spectrum
        self.s = numpy.real(numpy.exp(1j * self.phase_rad) * self.f[0] + numpy.exp(-1j * self.phase_rad) * self.f[1])
        
        # select part of the data
        # calculate the axes
        if axes == 0 or axes == [1,0] or axes == 2 or axes == [1,1]:
            if self.undersampling % 2 == 0:
                self.f[0] = self.f[0][:(len(self.f[0])/2)][:]
                self.f[1] = self.f[1][:(len(self.f[1])/2)][:]
                self.s = self.s[:(len(self.s)/2)][:]
            else:
                self.f[0] = self.f[0][(len(self.f[0])/2):][:]
                self.f[1] = self.f[1][(len(self.f[1])/2):][:]
                self.s = self.s[(len(self.s)/2):][:]  
                
            if flag_axis:
                self.s_axis[0] = MATH.make_ft_axis(length = 2*numpy.shape(self.s)[0], dt = self.r_axis[0][1]-self.r_axis[0][0], undersampling = self.undersampling)
                self.s_axis[0] = self.s_axis[0][0:len(self.s_axis[0])/2]
                         
        if axes == 1 or axes == [0,1] or axes == 2 or axes == [1,1]:
            if self.undersampling % 2 == 0:
                self.f[0] = self.f[0][:][:(len(self.f[0])/2)]
                self.f[1] = self.f[1][:][:(len(self.f[1])/2)]
                self.s = self.s[:][:(len(self.s)/2)]
            else:
                self.f[0] = self.f[0][:][(len(self.f[0])/2):]
                self.f[1] = self.f[1][:][(len(self.f[1])/2):]
                self.s = self.s[:][(len(self.s)/2):]
            
            if flag_axis:
                self.s_axis[2] = MATH.make_ft_axis(length = 2*numpy.shape(self.s)[0], dt = self.r_axis[2][1]-self.r_axis[2][0], undersampling = self.undersampling)
                self.s_axis[2] = self.s_axis[2][0:len(self.s_axis[0])/2]
        
        if flag_axis:
            if axes == 0 or axes == [1,0]:
                self.s_axis[2] = self.r_axis[2] + self.r_correction[2]   
            
            if axes == 1 or axes == [0,1]:
                self.s_axis[0] = self.r_axis[0] + self.r_correction[0]   
        
            
        # add some stuff to self
        self.s_units = ["cm-1", "fs", "cm-1"]
        
        if flag_axis:
            self.s_resolution = [(self.s_axis[0][1] - self.s_axis[0][0]), 0, (self.s_axis[2][1] - self.s_axis[2][0])]
            
        return True     
                  
    
            

    
    # PLOTTING 
    def plot(self, 
        plot_type = "S", 
        ax = False, 
        x_range = [0,0],
        y_range = [0,-1],
        pixel = -1,
        x_label = "", 
        y_label = "", 
        title = "", 
        flipxy = False,
        flag_verbose = False,
        **kwargs):

        """
        Plot the data.
        
        INPUT:
        - plot_type ('S' (default), 'R', 'NR', 'T'): plot the spectrum, the rephasing, non-rephasing spectra or time domain. Note that R and NR only work if the Fourier transform has been done and saved.
        - ax (False (default) or matplotlib axes instance): if False, it will make a new figure, otherwise it will use the axes instance, allowing subplots etc.
        - x_label, y_label, title (string, default=''): the labels for the axes. If no label is set, it will use the default. Use 'no_label' or 'no_title' to show no label.
        - pixel (int, -1): if pixel is an element of w_3, it will plot only this pixel, otherwise it will plot a 2D-plot
        - flipxy (BOOL, False): will flip the data and the axes. Non-default labels will not be flipped. 
        - x_range, y_range (array with 2 elements, [0,0], [0,-1]): the range to be plotted. 
        Possible cases:
        - [min, max]: plot range min to max
        - [0, 0]: plot the whole range
        - [0, -1]: use the range from the other axis. If both have this, it will plot both axes complete. (ie. it is identical to both having [0,0])
        
        **kwargs (stuff for the plotting function)
        - zlimit (number or list, -1): the z-range that will be used
            Possible cases:
            zlimit = 0, show all, not don't care about centering around zero
            zlimit = -1, show all, centered around zero
            zlimit = all else, use that, centered around zero
            zlimit = [a,b], plot from a to b
        - contours (number): number of contours to be used
        - invert_colors (BOOL, False): data = -data

        CHANGELOG:
        20110910/RB: started as croc-function
        20130131/RB: rewrote function for Crocodile. Integrated it with plot_T
        20130213/RB: some arguments are now kwargs: they are not explicitly listed but are passed on to the plotting function.        
        """

        self.verbose("Plot", flag_verbose)

        if plot_type == "S":
            data = self.s
        elif plot_type == "R":
            data = numpy.real(numpy.exp(-1j * self.phase_rad) * self.f[0])
        elif plot_type == "NR":
            data = numpy.real(numpy.exp(1j * self.phase_rad) * self.f[1])
        elif plot_type == "T":
            data = numpy.concatenate((numpy.flipud(self.r[1][1:,:]), self.r[0])).T
        
        if pixel < 0 or pixel > len(self.s_axis[2]):
            
            self.verbose("  -> contourplot", flag_verbose)
            
            # S, R and NR are freq-freq, T is time-freq
            # flipxy will change the xy-axis
            # the labels and range may be changed, if they are not at the defaults
            
            if plot_type == "S" or plot_type == "R" or plot_type == "NR":
                if flipxy:
                    # not the normal way
                    data = data.T
                    x_axis = self.s_axis[0]
                    y_axis = self.s_axis[2]
                    if x_label == "":
                        x_label = r"$\omega_1 (cm^{-1})$"
                    if y_label == "":
                        y_label = r"$\omega_3 (cm^{-1})$"   
                else:
                    x_axis = self.s_axis[2]
                    y_axis = self.s_axis[0]
                    if x_label == "":
                        x_label = r"$\omega_3 (cm^{-1})$"
                    if y_label == "":
                        y_label = r"$\omega_1 (cm^{-1})$"         
            else: # time domain
                if flipxy:
                    # not the normal way
                    data = data.T    
                    x_axis = self.r_axis[2]
                    y_axis = numpy.concatenate((-numpy.flipud(self.r_axis[0][1:]), self.r_axis[0]))
                    if x_label == "":
                        x_label = r"$\omega_3 (cm^{-1})$" 
                    if y_label == "":
                        y_label = r"$t_1 (fs)$" 
                    if y_range == [0,-1]: 
                        y_range = [0,0]
                else:     
                    x_axis = numpy.concatenate((-numpy.flipud(self.r_axis[0][1:]), self.r_axis[0]))
                    y_axis = self.r_axis[2]
                    if x_label == "":
                        x_label = r"$t_1 (fs)$" 
                    if y_label == "":
                        y_label = r"$\omega_3 (cm^{-1})$" 
                    if y_range == [0,-1]:
                        y_range = [0,0]            
            
            PL.contourplot(data, x_axis, y_axis, ax = ax, x_range = x_range, y_range = y_range, x_label = x_label, y_label = y_label, title = title, flag_verbose = flag_verbose, **kwargs)

        else:
            
            self.verbose("  -> linear plot", flag_verbose)
            
            if plot_type == "S" or plot_type == "R" or plot_type == "NR":
                x_axis = self.s_axis[0]
                x_label = r"$\omega_1 (cm^{-1})$"
            else:
                x_axis = self.r_axis[0]
                x_label = r"$\t_1 (fs)$"
                
            PL.linear(data[pixel], x_axis, x_range = x_range, y_range = y_range, ax = ax, x_label = x_label, y_label = y_label, title = title, flag_verbose = flag_verbose, **kwargs)

        if not ax:
            plt.show()




    def plot_R(self, **kwargs):
        """
        Wrapper for plot_type = R. See plot for more info.
        """
        self.plot(plot_type = "R", **kwargs)

    def plot_NR(self, **kwargs):
        """
        Wrapper for plot_type = R. See plot for more info.
        """
        self.plot(plot_type = "NR", **kwargs)
        
    def plot_T(self, **kwargs):
        """
        Wrapper for plot_type = R. See plot for more info.
        """
        self.plot(plot_type = "T", **kwargs)



        





    def save_data(self, export_path, r = False, s = True, flag_verbose = False):
        
        if s and numpy.shape(self.s) != (1,):
            IOM.save_data_PE(export_path, self.base_filename, s = self.s, s_axis = self.s_axis)
        if r and numpy.shape(self.r) != (3,):
            IOM.save_data_PE(export_path, self.base_filename, r = self.r, r_axis = self.r_axis)














### IMPLEMENTATION CLASSES ###
# If you really want to do something with them, move them to a separate file.

class pe_wt(pe):

    def __init__(self, objectname, flag_verbose = False):

        self.verbose("New Pe.pe.pe_wt class", flag_verbose)

        pe.__init__(self, objectname)


    def absorptive(self, window_function = "none", window_length = 0):

        self.super_absorptive(axes = [0,1], window_function = "none", window_length = 0)




class pe_tt(pe):

    def __init__(self, objectname, flag_verbose = False):

        self.verbose("New Pe.pe.pe_tt class", flag_verbose)

        pe.__init__(self, objectname)


    def absorptive(self, window_function = "none", window_length = 0):

        self.super_absorptive(axes = [1,1], window_function = "none", window_length = 0)



class pe_ww(pe):

    def __init__(self, objectname, flag_verbose = False):

        self.verbose("New Pe.pe.pe_wt class", flag_verbose)

        pe.__init__(self, objectname)


    def absorptive(self, window_function = "none", window_length = 0):

        self.super_absorptive(axes = [0,1], window_function = "none", window_length = 0)




if __name__ == "__main__": 
  
    pass
  
    # flag_verbose = True
    # 
    # a = pe("Fiets", flag_verbose)
    # 
    # print(a)