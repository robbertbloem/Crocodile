from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.ClassTools as CT

import Crocodile.Resources.Mathematics as MATH
import Crocodile.Resources.Equations as EQ




def find_indices(axis, range):
    """
    Find the minimum and maximum indices of the values given by range in axis.

    CHANGELOG:
    20120622/RB: generalization of function in line shape method

    USAGE:
    - axis (ndarray): an array with values. The values should be successive.
    - range (list of 2 int): the minimum and maximum of the range

    You have an axis AXIS from 2000 to 2200 cm-1 in 32 elements. If range RANGE = [2050, 2150] it will find the indices of the points between these two values, including those values. It will return the first and last index. 

    """
    try:
        min_i = numpy.where(axis < range[0])[0][-1]
    except: 
        min_i = 0

    try:
        max_i = numpy.where(axis > range[1])[0][0] + 1
    except: 
        max_i = len(axis)

    return min_i, max_i


class Lineshape(CT.ClassTools):
    """
    Some tools to calculate the line shape of spectra. The methods work for the spectra we are looking at right now, but I'm not sure how it will work out in the future. 

    CHANGELOG:
    20120622/RB: started the method, based on functions that were in an earlier separate script.
    20130212/RB: copied to Crocodile

    USAGE:
    When you initialize the class, you should give it a 'mess': a spectrum. This is a leftover from when it was an independent script. Then you should set some variables. The ones with an 'i' at the end indicate a minimum and maximum index. These can be found using 'croc.Resources.Function.find_indices(axis, range). You should also set fitting parameters. 


    """


    def __init__(self, mess, flag_verbose = False):
        """

        """

        self.verbose("Init lineshape class", flag_verbose = flag_verbose)

        self.mess = mess
        self.objectname = self.mess.objectname

        # double Lorentzian fit
        self.dl_x_i = [0,0]   # indices
        self.dl_y_i = [0,0]

        self.dl_A_in = [0,0,0, 0,0,0] # input fitting parameters
        self.dl_A = []      # fitting parameters
        self.dl_ble = []    # freq of minima 
        self.dl_esa = []    # freq of maxima

        # linear fit on min/max lorentzian fit
        self.l_i = [4, 8]       # linear fit, range, in indices OF THE INDICES OF THE DOUBLE LORENTZIAN

        self.l_A_in = [0, 1]    # fitting parameter
        self.l_A_ble = []
        self.l_A_esa = []

        self.l_angle_ble = 0
        self.l_angle_esa = 0

        self.l_slope_ble = 0
        self.l_slope_esa = 0

        ### center frequency using contour method ###
        # requires ellipse to run
        self.cf_ble = [0,0]
        self.cf_esa = [0,0]


        ### find peak heights ###
        self.ph_ble_x_i = [0,0]
        self.ph_ble_y_i = [0,0]

        self.ph_ble_max = 0
        self.ph_ble_min = 0

        self.ph_esa_x_i = [0,0]
        self.ph_esa_y_i = [0,0]

        self.ph_esa_max = 0
        self.ph_esa_min = 0

        ### find peak along w1 ###
        self.w1_peaks_x_i = [0,0]
        self.w1_peaks_y_i = [0,0]       

        self.w1_peaks_A_in = [0,0,0,0]

        self.w1_peaks = [0]

        ### plot related stuff ###
        self.color_array = ["b", "g", "r", "c", "m", "y", "k"]
        self.color_overlap_array = ["b", "r"]

        self.plot_x_i = [0,0]
        self.plot_y_i = [0,0]


    def fit_double_lorentzian(self, flag_plot = False, flag_verbose = False):
        """
        For a selection of points on the w1-axis, take a cut (giving w3 vs z (intensity) plot) and fit it with a double Lorentzian. 
        self.dl_x_i[0] etc are the min/max indices to be fitted
        """

        if flag_verbose:
            self.verbose("Fit double Lorentzian for " + self.objectname, flag_verbose = flag_verbose)
            self.verbose("  x_min: " + str(self.dl_x_i[0]) + " " + str(self.mess.s_axis[2][self.dl_x_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  x_max: " + str(self.dl_x_i[1]) + " " + str(self.mess.s_axis[2][self.dl_x_i[1]]), flag_verbose = flag_verbose)
            self.verbose("  y_min: " + str(self.dl_y_i[0]) + " " + str(self.mess.s_axis[0][self.dl_y_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  y_max: " + str(self.dl_y_i[1]) + " " + str(self.mess.s_axis[0][self.dl_y_i[1]]), flag_verbose = flag_verbose)

        # select the part of the data to be fitted
        data = self.mess.s[self.dl_y_i[0]:self.dl_y_i[1], self.dl_x_i[0]:self.dl_x_i[1]]
        x_axis = self.mess.s_axis[2][self.dl_x_i[0]:self.dl_x_i[1]]
        y_axis = self.mess.s_axis[0][self.dl_y_i[0]:self.dl_y_i[1]]
        
        # arrays for the results
        n_y, n_x = numpy.shape(data)
        y_max = numpy.zeros(n_y)                # index of the maximum
        y_min = numpy.zeros(n_y)                # index of the minimum
        y_out_array = numpy.zeros((n_y, 8))     # fitting parameters

        if flag_plot:
            plt.figure()
            color_array = ["b", "g", "r", "c", "m", "y", "k"]

        # calculate the fit for the cut of w1
        for i in range(n_y):

            y = data[i,:]

            A_out = MATH.fit(x_axis, y, EQ.rb_two_lorentzians, self.dl_A_in)        

            y_out_array[i,:] = A_out

            x_fit = numpy.arange(x_axis[0], x_axis[-1], 0.1)
            y_fit = EQ.rb_two_lorentzians(A_out, x_fit)

            if flag_plot:  
                plt.plot(x_fit, y_fit, c = color_array[i%len(color_array)])
                plt.plot(x_axis, y, ":", c = color_array[i%len(color_array)])

            y_max[i] = x_fit[numpy.argmax(y_fit)]
            y_min[i] = x_fit[numpy.argmin(y_fit)]

        self.dl_ble = y_min
        self.dl_esa = y_max
        self.dl_A = y_out_array

        if flag_plot:
            plt.show()


    def fit_tilt(self, flag_verbose = False):

        if flag_verbose:
            self.verbose("Fit tilt for " + self.objectname, flag_verbose = flag_verbose)
            self.verbose("  x_min: " + str(self.dl_x_i[0]) + " " +  str(self.mess.s_axis[2][self.dl_x_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  x_max: " + str(self.dl_x_i[1]) + " " +  str(self.mess.s_axis[2][self.dl_x_i[1]]), flag_verbose = flag_verbose)
            self.verbose("  y_min: " + str(self.dl_y_i[0] + self.l_i[0]) + " " +  str(self.mess.s_axis[0][self.dl_y_i[0] + self.l_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  y_max: " + str(self.dl_y_i[0] + self.l_i[1]) + " " +  str(self.mess.s_axis[0][self.dl_y_i[0] + self.l_i[1]]), flag_verbose = flag_verbose)

        y = self.mess.s_axis[0][self.dl_y_i[0] + self.l_i[0]:self.dl_y_i[0] + self.l_i[1]]     

        x = self.dl_ble[self.l_i[0]:self.l_i[1]]
        self.l_A_ble = MATH.fit(x, y, EQ.linear, self.l_A_in)

        x = self.dl_esa[self.l_i[0]:self.l_i[1]]
        self.l_A_esa = MATH.fit(x, y, EQ.linear, self.l_A_in)

        self.l_angle_ble = 90 - numpy.arctan(self.l_A_ble[1]) * 180 / numpy.pi
        self.l_angle_esa = 90 - numpy.arctan(self.l_A_esa[1]) * 180 / numpy.pi

        self.l_slope_ble = 1 / self.l_A_ble[1]
        self.l_slope_esa = 1 / self.l_A_esa[1]


    def find_peak_heights(self, flag_verbose = False):

        if flag_verbose:
            self.verbose("Find peak heights for bleach of " + self.objectname, flag_verbose = flag_verbose)
            self.verbose("  x_min: " + str(self.ph_ble_x_i[0]) + " " + str(self.mess.s_axis[2][self.ph_ble_x_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  x_max: " + str(self.ph_ble_x_i[1]) + " " + str(self.mess.s_axis[2][self.ph_ble_x_i[1]]), flag_verbose = flag_verbose)
            self.verbose("  y_min: " + str(self.ph_ble_y_i[0]) + " " + str(self.mess.s_axis[0][self.ph_ble_y_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  y_max: " + str(self.ph_ble_y_i[1]) + " " + str(self.mess.s_axis[0][self.ph_ble_y_i[1]]), flag_verbose = flag_verbose)

        data = self.mess.s[self.ph_ble_y_i[0]:self.ph_ble_y_i[1], self.ph_ble_x_i[0]:self.ph_ble_x_i[1]]

        self.ph_ble_min = numpy.min(data)
        self.ph_ble_max = numpy.max(data)

        if flag_verbose:
            self.verbose("Find peak heights for ESA of " + self.objectname, flag_verbose = flag_verbose)
            self.verbose("  x_min: " + str(self.ph_esa_x_i[0]) + " " + str(self.mess.s_axis[2][self.ph_esa_x_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  x_max: " + str(self.ph_esa_x_i[1]) + " " + str(self.mess.s_axis[2][self.ph_esa_x_i[1]]), flag_verbose = flag_verbose)
            self.verbose("  y_min: " + str(self.ph_esa_y_i[0]) + " " + str(self.mess.s_axis[0][self.ph_esa_y_i[0]]), flag_verbose = flag_verbose)
            self.verbose("  y_max: " + str(self.ph_esa_y_i[1]) + " " + str(self.mess.s_axis[0][self.ph_esa_y_i[1]]), flag_verbose = flag_verbose)

        data = self.mess.s[self.ph_esa_y_i[0]:self.ph_esa_y_i[1], self.ph_esa_x_i[0]:self.ph_esa_x_i[1]]

        self.ph_esa_min = numpy.min(data)
        self.ph_esa_max = numpy.max(data)

    def find_w1_peaks(self, flag_verbose = False):

        data = self.mess.s[self.w1_peaks_y_i[0]:self.w1_peaks_y_i[1], self.w1_peaks_x_i[0]:self.w1_peaks_x_i[1]]
        y_axis = self.mess.s_axis[0][self.w1_peaks_y_i[0]:self.w1_peaks_y_i[1]]        

        y = numpy.sum(data,1)

        A_out = MATH.fit(y_axis, y, EQ.rb_lorentzian, self.w1_peaks_A_in)

        self.w1_peaks[0] = A_out[1]





