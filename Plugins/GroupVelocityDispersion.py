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
    
    CHANGELOG:
    20110909/RB: started
    20150805/RB: copied to Crocodile
    
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
    n_y = numpy.sqrt(E.Sellmeier(SC, n_x))
    
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
    
    INPUT:
    t0 (array or int): the pulse length(s) of the pulse before the gvd
    gvd (numpy.float or ndarray): the GVD in fs^2 (!!!)
    
    OUTPUT:
    t1 (ndarray): the new pulse length as a function of GVD.
    
    """
    
    if type(t0) == int:
        t0 = [t0]
    
    for i in range(len(t0)):
        t1 = t0[i] * numpy.sqrt(1 + (4 * numpy.log(2) * gvd / t0[i]**2)**2 )
    
    return t1

def gvd_function_wavelength(material, range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0, n_steps = 1000):

    [n, SC] = GVDM.MaterialProperties(material)
    
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
                ax.set_ylim(0, numpy.max(p_y)*1.05)
                ax.set_title(label)
                plt.show()
                
def pulse_length_function_wavelength(material, material_mm = 1.0, pulse_length_fs = 100, range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0, n_steps = 1000):

    [n, SC] = GVDM.MaterialProperties(material)
    
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
                ax.set_ylim(0, numpy.max(p_y)*1.05)
                ax.set_title(label)
                plt.show()
                


def n_function_wavelength(material, range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0, n_steps = 1000):

    [n, SC] = GVDM.MaterialProperties(material)
    
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

            x = numpy.linspace(range_um[0], range_um[1], n_steps)
            n = numpy.sqrt(E.Sellmeier(SC, x))
  
            label = "Index of refraction of %s" % (material)
            
            ax.plot(x, n, label = label)
            ax.set_xlabel("Wavelength (micron)")
            ax.set_ylabel("Index of refraction")
            ax.legend(loc = 4)
            
            if own_plot:
                ax.set_title(label)
                plt.show()

def gvd_function_wavelength_for_materials(materials = [], range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for material in materials:
        
        gvd_function_wavelength(material, range_um = range_um, print_for_um = print_for_um, flag_plot = flag_plot, ax = ax, n_steps = 1000)

def pulse_length_function_wavelength_for(materials = [], material_mms = [], pulse_length_fss = 100, range_um = [3,7], print_for_um = [], flag_plot = True, ax = 0):
    
    if ax == 0:
        own_plot = True
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        own_plot = False


    for pulse_length_fs in pulse_length_fss:
        for material in materials:
            for material_mm in material_mms:
        
                pulse_length_function_wavelength(material, material_mm = material_mm, pulse_length_fs = pulse_length_fs, range_um = range_um, print_for_um = print_for_um, flag_plot = flag_plot, ax = ax, n_steps = 1000)
        
    ax.set_title("Pulse length for different materials")
    ax.set_ylim(0, ax.get_ylim()[1])

    if own_plot:
       plt.plot()







if __name__ == "__main__": 
    
    plt.close("all")
    
    
    n, SC = GVDM.MaterialProperties("caf2")
#   print(n, SC)
    
    
#    range_L = [3,7]
#    n_steps = 10
#    print(calculate_gvd(range_L, n_steps, SC))
    
    material = "caf2"
    range_um = [3,8]
    print_for_um = [4,5,6]
    pulse_length_fs = 50
    material_mm = 2.0
    
    
    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax = 0
    
#    n_function_wavelength(material, range_um = range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)
#
#    n_function_wavelength(material = "baf2", range_um = range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)

    
#    gvd_function_wavelength(material, range_um = range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)
#
#    gvd_function_wavelength(material = "baf2", range_um = range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)
    
#    pulse_length_function_wavelength(material, pulse_length_fs = pulse_length_fs, material_mm = material_mm, range_um = range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)

#    pulse_length_function_wavelength(material = material, pulse_length_fs = 100, material_mm = material_mm, range_um = range_um, print_for_um = print_for_um, flag_plot = True, ax = ax)

#    ax.set_ylim(0, ax.get_ylim()[1])

    gvd_function_wavelength_for_materials(materials = ["baf2", "caf2", "znse"], range_um = range_um, print_for_um = print_for_um, flag_plot = True, ax = 0)

#    pulse_length_function_wavelength_for(materials = ["caf2", "znse"], material_mms = [1.0, 2.0], pulse_length_fss = [100], range_um = [3,7], print_for_um = [6], flag_plot = True, ax = 0)

#    GroupVelocityDispersion(material = "caf2", print_for_um = [4,5,6], range_um = [3,7], material_path_mm = 2, pulse_length_fs = 100, n_steps = 100, flag_plot_n = False, flag_plot_gvd = False, y_range_gvd = [0,0], flag_plot_pulse_length = True)

    plt.show()
    