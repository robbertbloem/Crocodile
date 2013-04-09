from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.Mathematics as MA
import Crocodile.Resources.Equations as EQ
import Crocodile.Resources.Constants as CONST

def envelope(A,t):
    """
    Wrapper for rb_gaussian
    
    INPUT:
    t: numpy.array
    A[0]: sigma (sigma^2 = variance)
    A[1]: mu (mean)
    A[2]: offset 
    A[3]: scale, before offset
    
    CHANGELOG:
    20130408/RB: started function
    """
    return EQ.rb_gaussian(A,t)

def pulse(A,t):
    """
    
    INPUT:
    t: numpy.array
    Gaussian:
    A[0]: sigma (sigma^2 = variance)
    A[1]: mu (mean)
    A[2]: offset 
    A[3]: scale, before offset
    Cosine:
    A[4]: frequency
    A[5]: phase
    
    Adapted from EQ.rb_gaussian() and EQ.rb_cos().
    
    CHANGELOG:
    20130408/RB: started function
    """
    return (A[3]/(A[0]*numpy.sqrt(2*numpy.pi))) * numpy.exp(-(t-A[1])**2/(2*A[0]**2)) * (numpy.cos(2 * numpy.pi * A[4] * t + numpy.pi*A[5]))


def fit_convolve(A,t):
    """
    Fit the convolution of pulse()

    INPUT:
    See pulse()

    CHANGELOG:
    20130408/RB: started function
    """
    y_s = pulse(A, t)
    return numpy.convolve(y_s, y_s, mode='same')   
    

def plot_res(t, y_data, y_conv, y_single, y_env):
    """
    Plot the results of the fitting procedure.
    
    INPUT:
    t: array for time axis
    y_data: original data
    y_conv: the fitted convoluted single pulse
    y_single: the single pulse
    y_env: the envelope
    All inputs are arrays and should have the same length
    
    CHANGELOG:
    20130408/RB: started function
    """
    fig = plt.figure()
    ax = []
    for i in range(2):
        ax.append(fig.add_subplot(2,1,i+1))
    
    ax[0].plot(t, y_data)
    ax[0].plot(t, y_conv)
    
    ax[1].plot(t, y_single)
    ax[1].plot(t, y_env)

    for i in range(2):
        ax[i].set_xlim(t[0], t[-1])
        ax[i].set_xlabel("time (fs)")
    
    ax[0].set_title("Original data (blue) and convoluted pulse (green)")
    ax[1].set_title("Single pulse (blue) and envelope (green)")

    plt.show()

def print_res(A, FWHM, dt):
    """
    Print the results of the fitting.
    
    INPUT:
    A: the fitting parameters
    FWHM: full-width-half-maximum, in fs
    dt: time step
    
    CHANGELOG:
    20130408/RB: started function
    """
    base = 20
    print("Fitting results")
    print("%10s = %.2f fs" % ("Sigma", A[0]))
    print("%10s = %.2f" % ("Mu", A[1]))
    print("%10s = %.2f" % ("Offset", A[2]))
    print("%10s = %.2f" % ("Scale", A[3]))
    print("%10s = %.2f 1/fs" % ("Frequency", A[4]))
    print("%10s = %.2f fs" % ("Period", 1/A[4]))
    print("%10s = %.2f" % ("Phase", A[5]))
    print()
    print("FWHM")
    print(u"%10s = %.1f fs (\u00B1 %.1f)" % ("FWHM", FWHM, dt))


def default_A():
    sigma = 50 # fs
    mu = 0      
    offset = 0
    scale = 50
    frequency = 0.055 # 1/fs
    phase = 0
    return [sigma, mu, offset, scale, frequency, phase]


def calculate_pulse_length(data, A = [], dt = CONST.hene_fringe_fs, plot_results = True, print_results = True, fwhm_limit = 200, mu_shift = 10, flag_recursive_limit = False):

    """
    Calculate the pulse length
    
    INPUT:
    data (array): the result of a measure phase measurement for one pulse pair. Usually measure phase is saved as two measurements of two pulse pairs. 
    A (list, default = []): [sigma, mu, offset, scale, frequency, phase] if A is empty, some reasonable default values are used.
    dt (float, HeNe-fringe in fs): to convert from indices to time
    plot (BOOL, True): plot the results
    print_res (BOOL, True): print the results
    fwhm_limit (float, 200): if FWHM is above this limit, it will try again by shifting the mean a bit (by mu_shift)
    mu_shift (float, 10): shift the mean a bit
    flag_recursive_limit (BOOL, False): function will try again if set to False. To prevent endless recursion.
    
    OUTPUT:
    - the output can be plotted to confirm the fit
    - the fitting results can be printed
    - the FWHM is returned
    
    CHANGELOG:
    20130408/RB: started function
    """

    # input for a single IR pulse    
    if A == []:
        A = default_A()
    
    data /= numpy.amax(data)
    
    # make the time axis
    t = numpy.arange(len(data)) - len(data)/2
    t *= dt
    
    # fit the single pulse by convoluting it with itself and fitting it to the data
    A_out = MA.fit(t, data, fit_convolve, A)
    
    # calculate the fitted pulse and convoluted pulse
    y_single = pulse(A_out, t)
    y_conv = fit_convolve(A_out, t)

    # calculate the envelope
    sigma = A_out[0] 
    mu = A_out[1]
    offset = A_out[2] 
    scale = A_out[3] 
    A = [sigma, mu, offset, scale]
    y_env = envelope(A,t)

    # for the FWHM we need a higher time resolution
    t2 = numpy.arange(len(data)*10) - 5*len(data)
    t2 *= dt / 10
    y_fwhm = envelope(A,t2)
     # use the envelope to calculate the FWHM
    # make it positive
    if y_fwhm[int(len(t2)/2)] < 0:
        y_fwhm = -y_fwhm
    max_fwhm = numpy.amax(y_fwhm)
    # select the part of the list
    l = numpy.where(y_fwhm > max_fwhm/2)[0]
    # the length of the list in indices, multiply with time distance
    FWHM = len(l) * dt / 10
 
    if FWHM > fwhm_limit and flag_recursive_limit == False:
        # FWHM is above the limit, try again by shifting the mean a bit
        print("Trying again...")
        A = default_A()
        A[1] = mu_shift
        FWHM = calculate_pulse_length(data, A = A, dt = dt, plot_results = plot_results, print_results = print_results, flag_recursive_limit = True)
    
    else:
        # yay! a reasonable result. plot and print the results.
        if plot_results:
            plot_res(t, data, y_conv, y_single, y_env)
        if print_results:
            print_res(A_out, FWHM, dt)
    
    return FWHM


