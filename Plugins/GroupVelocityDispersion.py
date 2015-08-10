"""
This script calculates the index of refraction (n), group velocity (VG) and group velocity dispersion (GVD). 

The calculation of the index of refraction is done using the Sellmeier Equation, the parameters for this are stored in GVD_Materials. The VG is the first derivative of this, the GVD the second derivative. This is calculated using a fairly simple (y2-y1)/(x2-x1) calculation. To accomodate this, these are all calculated for a range of values. 

All wavelengths are in micron!

SOURCES (for calculations):
http://www.rp-photonics.com/chromatic_dispersion.html

SOURCES (for data):
RI: http://refractiveindex.info -> WARNING: C coefficients need to be squared SOMETIMES!
MG: http://www.cvimellesgriot.com


"""


from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Plugins.GVD_Materials as GVDM
import Crocodile.Resources.Constants as C
import Crocodile.Resources.Equations as E
import Crocodile.Resources.Mathematics as M
import Crocodile.Resources.Functions as F


def calculate_gvd(range_L, n_steps, SC):
    """
    Function to calculate the Group Velocity (VG) and Group Velocity Dispersion (GVD) using the Sellmeier Equation.
    
    The index of refraction as a function of wavelength is described by the Sellmeier equation. The VG is the first derivative of this, the GVD the second derivative. 
    
    This is calculated using a fairly simple (y2-y1)/(x2-x1) calculation. To accomodate this, these are all calculated for a range of values.
    
    CHANGELOG:
    20110909/RB: started
    20150805/RB: copied to Crocodile. Added documentation. 
    
    INPUT:
    range_L (array, 2 elements): minimum and maximum in micron
    n_steps (int): number of steps for the calculation. It is best to use a large number of steps. 
    SC (array): Sellmeier Coefficients
    
    OUPUT:
    n_x: x-axis for Sellmeier (wavelengths in micron) 
    n_y: index of refraction
    vg_x: x-axis (wavelengths in nm) for group velocity
    vg_y: group velocity
    gvd_x: x-axis (wavelengths in nm) 
    gvd_y: group velocity dispersion
    
    """

    # calculate the index of refraction over the defined range
    n_x = numpy.linspace(range_L[0], range_L[1], num = n_steps)
    n_y = E.Sellmeier(SC, n_x)
    
    # calculate the derivatives
    [d1x, d1y] = numpy.array(M.derivative(n_x, n_y))
    [d2x, d2y] = numpy.array(M.derivative(d1x, d1y))
    
    # define the x axes
    vg_x = d1x
    gvd_x = d2x
    
    vg_y = (C.c_ms / (vg_x)) / (1 - (vg_x / n_y[1:-1]) * d1y)
    gvd_y = (gvd_x**3/(2*numpy.pi*C.c_ms*C.c_ms)) * d2y
    
    gvd_y *= 1e21
    
    return n_x, n_y, vg_x, vg_y, gvd_x, gvd_y


def calculate_effect_GVD(t0, gvd):
    """
    Calculates the effect the group velocity dispersion has on an unchirped Gaussian pulse. 
    
    Source:
    http://www.rp-photonics.com/chromatic_dispersion.html
    numpy.log() is natural log
    
    CHANGELOG:
    20120221/RB: started the function
    20150807/RB: copied to Crocodile. Added some documentation.
    
    INPUT:
    t0 (array or int): the pulse length(s) of the pulse before the gvd
    gvd (numpy.float or ndarray): the GVD in fs^2 (!!!). You have to multiply the GVD with the thickness of the material. 
    
    OUTPUT:
    t1 (ndarray): the new pulse length as a function of GVD.
    
    """
    
    if type(t0) == int:
        t0 = [t0]
    
    for i in range(len(t0)):
        t1 = t0[i] * numpy.sqrt(1 + (4 * numpy.log(2) * gvd / t0[i]**2)**2 )
    
    return t1


def check_wavelength_range(plot_range_um, print_for_um):
    """
    Sanity check for the plot_range_um and print_for_um. 
    
    INPUT:
    plot_range_um: the range to be plotted, for example from 3 to 8 micron
    print_for_um: a list for which the results are printed. 
    
    OUTPUT:
    plot_range_um: the plot range, with the lower boundary the lowest wavelength
    range_um: a range to be calculated. A selection of the results can be printed or plotted. 
    
    """
    
    
    if plot_range_um[0] > plot_range_um[1]:
        plot_range_um = plot_range_um[::-1]
    
    range_um = plot_range_um[:]
    
    for w in print_for_um:
        if w < range_um[0]:
            range_um[0] = w
        if w > range_um[1]:
            range_um[1] = w
            
        range_um[0] *= 0.9
        range_um[1] *= 1.1
        
    return plot_range_um, range_um


def gvd_function_wavelength(material, plot_range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0, n_steps = 1000):

    """
    Calculate the GVD as a function of wavelength for a particular material. Prints the results for particular wavelengths and plots the results. 
    
    INPUT:
    material (string): a string corresponding to the list in GVD_Materials
    plot_range_um (list, len=2): list with lower and upper wavelength (in micron)  
    print_for_um (list): list with wavelengths to be printed (in micron)
    flag_plot (bool): whether to make a plot
    ax (0 or axis instance): if zero, a new plot will be generated. If a pyplot axis instance, it will be added to that particular axis.  
    n_steps (int): the number of steps to be calculated. A larger number will give a higher resolution, especially when calculating VG and GVD. 
    
    """

    [n, SC] = GVDM.MaterialProperties(material)
    
    plot_range_um, range_um = check_wavelength_range(plot_range_um, print_for_um)
    
    if SC != []:
    
        n_x, n_y, vg_x, vg_y, gvd_x, gvd_y = calculate_gvd(range_um, n_steps, SC)
    
        if print_for_um != []:        
            print("Material: %s" % (material))
            print(" um     gvd (fs^2/mm)")
            for i in range(len(print_for_um)):
                index = numpy.where(gvd_x >= print_for_um[i])[0]

                print("%3.3f %12.3f" % (gvd_x[index[0]], gvd_y[index[0]]))
        
        if flag_plot:
        
            if ax == 0:
                own_plot = True
                fig = plt.figure()
                ax = fig.add_subplot(111)
            else:
                own_plot = False
        
            label = "Group velocity dispersion of %s" % (material)
            
            ax.plot(gvd_x, gvd_y, label = label)
            ax.set_xlabel("Wavelength (micron)")
            ax.set_ylabel("Group velocity dispersion (fs^2/mm)")
            ax.legend(loc = 4)
            
            if own_plot:
                ax.set_xlim(plot_range_um)
                ax.set_title(label)
                plt.show()

              
def pulse_length_function_wavelength(material, material_mm = 1.0, pulse_length_fs = 100, plot_range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0, n_steps = 1000):

    """
    Calculate how the pulse length is affected by a material, as a function of wavelength for a particular material. Prints the results for particular wavelengths and plots the results. 
    
    INPUT:
    material (string): a string corresponding to the list in GVD_Materials
    material_mm (float): thickness of the material
    pulse_length_fs (float): length of the pulse before entering the material.
    plot_range_um (list, len=2): list with lower and upper wavelength (in micron)  
    print_for_um (list): list with wavelengths to be printed (in micron)
    flag_plot (bool): whether to make a plot
    ax (0 or axis instance): if zero, a new plot will be generated. If a pyplot axis instance, it will be added to that particular axis.  
    n_steps (int): the number of steps to be calculated. A larger number will give a higher resolution, especially when calculating VG and GVD. 
    
    """


    [n, SC] = GVDM.MaterialProperties(material)
    
    plot_range_um, range_um = check_wavelength_range(plot_range_um, print_for_um)
    
    if SC != []:
    
        n_x, n_y, vg_x, vg_y, gvd_x, gvd_y = calculate_gvd(range_um, n_steps, SC)
        
        p_y = numpy.zeros(len(gvd_y))
        for i in range(len(gvd_y)):
            p_y[i] = calculate_effect_GVD(pulse_length_fs, gvd_y[i] * material_mm)
    
        if print_for_um != []:        
            print("Material: %s" % (material))
            print("Initial pulse length: %.1f fs" % (pulse_length_fs))
            print(" um     gvd (fs^2/mm)  pulse length (fs)")
            for i in range(len(print_for_um)):
                index = numpy.where(gvd_x >= print_for_um[i])[0]

                print("%3.3f %12.3f %12.1f" % (gvd_x[index[0]], gvd_y[index[0]], p_y[index[0]]))
        
        if flag_plot:
        
            if ax == 0:
                own_plot = True
                fig = plt.figure()
                ax = fig.add_subplot(111)
            else:
                own_plot = False
  
            label = "%.0f fs pulse after %.1f mm %s" % (pulse_length_fs, material_mm, material)
            
            ax.plot(gvd_x, p_y, label = label)
            ax.set_xlabel("Wavelength (micron)")
            ax.set_ylabel("Pulse length (fs)")
            ax.legend(loc = 4)
            
            if own_plot:                  
                ax.set_xlim(plot_range_um)
                ax.set_ylim(0, ax.get_ylim()[1])
                ax.set_title(label)
                plt.show()


def n_function_wavelength(material, plot_range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0, n_steps = 1000):

    """
    Calculate the index of refraction as a function of wavelength for a particular material. Prints the results for particular wavelengths and plots the results. 
    
    INPUT:
    material (string): a string corresponding to the list in GVD_Materials
    plot_range_um (list, len=2): list with lower and upper wavelength (in micron)  
    print_for_um (list): list with wavelengths to be printed (in micron)
    flag_plot (bool): whether to make a plot
    ax (0 or axis instance): if zero, a new plot will be generated. If a pyplot axis instance, it will be added to that particular axis.  
    n_steps (int): the number of steps to be calculated. A larger number will give a higher resolution, especially when calculating VG and GVD. 
    
    """
    
    [n, SC] = GVDM.MaterialProperties(material)
    
    plot_range_um, range_um = check_wavelength_range(plot_range_um, print_for_um)
    
    if SC != []:
        
        if print_for_um != []:        
            n = numpy.sqrt(E.Sellmeier(SC, numpy.array(print_for_um)))
            print("Material: %s" % (material))
            print(" um        n")
            for i in range(len(print_for_um)):
                print("%3.3f %12.3f" % (print_for_um[i], n[i]))
        
        if flag_plot:

            if ax == 0:
                own_plot = True
                fig = plt.figure()
                ax = fig.add_subplot(111)
            else:
                own_plot = False

            x = numpy.linspace(plot_range_um[0], plot_range_um[1], n_steps)
            n = numpy.sqrt(E.Sellmeier(SC, x))
  
            label = "Index of refraction of %s" % (material)
            
            ax.plot(x, n, label = label)
            ax.set_xlabel("Wavelength (micron)")
            ax.set_ylabel("Index of refraction")
            ax.legend(loc = 4)
            
            if own_plot:
                ax.set_title(label)
                plt.show()


def gvd_function_wavelength_for_materials(materials = [], plot_range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0):
    
    """
    Calculate the GVD as a function of wavelength, for a range of materials. 
    
    materials (list): a list with materials, selected from GVD_Materials
    plot_range_um (list, len=2): list with lower and upper wavelength (in micron)  
    print_for_um (list): list with wavelengths to be printed (in micron)
    flag_plot (bool): whether to make a plot
    ax (0 or axis instance): if zero, a new plot will be generated. If a pyplot axis instance, it will be added to that particular axis.  

    """

    if ax == 0 and flag_plot == True:
        own_plot = True
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        own_plot = False

    for material in materials:
        
        gvd_function_wavelength(material, plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = flag_plot, ax = ax, n_steps = 1000)

    if own_plot:
       plt.show()


def pulse_length_function_wavelength_for(materials = [], material_mms = [], pulse_length_fss = 100, plot_range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0):
    """
    Calculate the increase in pulse length, for given initial pulse lengths, materials, and thicknesses. 

    materials (list with strings): a list with materials, selected from GVD_Materials
    material_mms (list with floats): thicknesses to be calculated.
    pulse_length_fss (list with floats): initial pulse lengths. 
    plot_range_um (list, len=2): list with lower and upper wavelength (in micron)  
    print_for_um (list): list with wavelengths to be printed (in micron)
    flag_plot (bool): whether to make a plot
    ax (0 or axis instance): if zero, a new plot will be generated. If a pyplot axis instance, it will be added to that particular axis.  

    """   
    if ax == 0 and flag_plot == True:
        own_plot = True
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        own_plot = False


    for pulse_length_fs in pulse_length_fss:
        for material in materials:
            for material_mm in material_mms:
        
                pulse_length_function_wavelength(material, material_mm = material_mm, pulse_length_fs = pulse_length_fs, plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = flag_plot, ax = ax, n_steps = 1000)
        
    ax.set_title("Pulse length for different situations")
    ax.set_ylim(0, ax.get_ylim()[1])

    if own_plot:
       plt.show()







if __name__ == "__main__": 
    
    plt.close("all")
    
    # get material properties
    n, SC = GVDM.MaterialProperties("caf2")

    # set some values
    material = "caf2"
    plot_range_um = [3,8]
    print_for_um = [4,5,6]
    pulse_length_fs = 50
    material_mm = 2.0
        
#     plot_range_um, range_um = check_wavelength_range(plot_range_um = plot_range_um, print_for_um = print_for_um)

#     n_function_wavelength(material, plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = True, ax = 0)

#     n_function_wavelength(material = "baf2", plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)

#     gvd_function_wavelength(material, plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = True, ax = 0)

#    gvd_function_wavelength(material = "baf2", plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)
    
#     pulse_length_function_wavelength(material, pulse_length_fs = pulse_length_fs, material_mm = material_mm, plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = True, ax = 0)

#    pulse_length_function_wavelength(material = material, pulse_length_fs = 100, material_mm = material_mm, plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)

#     gvd_function_wavelength_for_materials(materials = ["baf2", "caf2", "znse"], plot_range_um = plot_range_um, print_for_um = print_for_um, flag_plot = True, ax = 0)

#    pulse_length_function_wavelength_for(materials = ["caf2", "znse"], material_mms = [1.0, 2.0], pulse_length_fss = [100], plot_range_um = [3,7], print_for_um = [6], flag_plot = True, ax = 0)

    plt.show()
    