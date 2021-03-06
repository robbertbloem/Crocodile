"""
croc.Resources.Mathematics



"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect

import numpy

import numpy
import matplotlib 
import matplotlib.pyplot as plt

try:
    from scipy.optimize.minpack import leastsq
    scipy_import = True
except ImportError:
    scipy_import = False

import itertools

import Crocodile.Resources.Constants as CONST
import PythonTools.Debug as DEBUG


def square_formula(a, b, c):
    """
    Calculates 
    ax^2 + bx + c = 0

    """
    x0 = (-b + numpy.sqrt(b**2 - 4*a*c))/(2*a)
    x1 = (-b - numpy.sqrt(b**2 - 4*a*c))/(2*a)
    return x0, x1


### CODE FOR FITTING PROCEDURE ###

def minimize(A, t, y0, function):
    """
    Needed for fit
    """
    return y0 - function(A, t)

def fit(x_array, y_array, function, A_start, return_all = False):
    """
    Fit data
    
    20101209/RB: started
    20130131/RB: imported in Crocodile, added example to doc-string

    INPUT:
    x_array: the array with time or something
    y-array: the array with the values that have to be fitted
    function: one of the functions, in the format as in the file "Equations"
    A_start: a starting point for the fitting
    return_all: the function used to return only the final result. The leastsq method does however return more data, which may be useful for debugging. When the this flag is True, it will return these extras as well. For legacy purposes the default is False. See reference of leastsq method for the extra output: http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html
    
    OUTPUT:
    A_final: the final parameters of the fitting
    When return_all == True:
    - cov_x (ndarray): Uses the fjac and ipvt optional outputs to construct an estimate of the jacobian around the solution. None if a singular matrix encountered (indicates very flat curvature in some direction). This matrix must be multiplied by the residual variance to get the covariance of the parameter estimates - see curve_fit.
    - infodict (dict): a dictionary of optional outputs with the key s:
        - "nfev" : the number of function calls
        - "fvec" : the function evaluated at the output
        - "fjac" : A permutation of the R matrix of a QR
                 factorization of the final approximate
                 Jacobian matrix, stored column wise.
                 Together with ipvt, the covariance of the
                 estimate can be approximated.
        - "ipvt" : an integer array of length N which defines
                 a permutation matrix, p, such that
                 fjac*p = q*r, where r is upper triangular
                 with diagonal elements of nonincreasing
                 magnitude. Column j of p is column ipvt(j)
                 of the identity matrix.
        - "qtf"  : the vector (transpose(q) * fvec).
    - mesg (str): A string message giving information about the cause of failure.
    - ier (int): An integer flag. If it is equal to 1, 2, 3 or 4, the solution was found. Otherwise, the solution was not found. In either case, the optional output variable "mesg" gives more information.


    EXAMPLE:
    Fit some data to this function from Crocodile.Resources.Equations:
    def linear(A, t):
        return A[0] + A[1] * t  
    
    ### 
    x = x-axis
    y = some data
    A = [0,1] # initial guess
    A_final = fit(x, y, Crocodile.Resources.Equations.linear, A)
    ###
    
    WARNING:
    Always check the result, it might sometimes be sensitive to a good starting point.

    """
    if scipy_import:
        param = (x_array, y_array, function)
    
        A_final, cov_x, infodict, mesg, ier = leastsq(minimize, A_start, args=param, full_output=True)

        if return_all:
            return A_final, cov_x, infodict, mesg, ier
        else:
            return A_final
    else:
        DEBUG.printError("Scipy.leastsq is not loaded. Fit is not done", inspect.stack())
        return False



### FOURIER TRANSFORM STUFF ###

# the general Fourier transform method
def fourier(array, 
    zero_in_middle = False, 
    first_correction = False, 
    zeropad_to = None, 
    window_function = "none", 
    window_length = 0, 
    flag_plot = False, 
    flag_verbose = False):
    """
    A Fourier transform for any dimension.

    INPUT:
    - array (x-dimensions ndarray): to be FFT'd
    - zero_in_middle (BOOL): for FFT the zero-time should be the first element of the array. If the zero is in the middle, it will be shifted first
    - first_correction (BOOL): if the first element of the array has to be halved, check this as True
    - zeropad_to (number): Length of the transformed axis of the output. If n is smaller than the length of the input, the input is cropped. If it is larger, the input is padded with zeros. If n is None, the length of the input (along the axis specified by axis) is used.

    OUTPUT:
    array (x-dimensions ndarray): Fourier transformed array

    CHANGELOG:
    20101204 RB: started
    20110909 RB: added zeropadding

    """
    
    DEBUG.verbose("FFT", flag_verbose)
    
    # shift time = 0 to first element
    if zero_in_middle == True:
        array = numpy.fft.ifftshift(array)

    # half the first element
    if first_correction == True: 
        dim = len(numpy.shape(array))
        if dim == 1:
            array[0] /= 2
        elif dim == 2:
            array[0,:] /= 2
            array[:,0] /= 2
        elif dim > 2:
            DEBUG.printError("Correction of the first element is not done!", inspect.stack())

    # window function
    if window_function != "none": 
        array = window_functions(array, window_function, window_length, flag_plot = flag_plot)

    # the fft
    array = numpy.fft.fft(array, n = zeropad_to)

    # move the array back if it was shifted
    if zero_in_middle == True:
        array = numpy.fft.fftshift(array)

    return array 




def window_functions(array, 
    window_function, 
    window_length = 0, 
    flag_plot = False,
    flag_verbose = False):
    """
    croc.Absorptive.window_functions

    Different window functions.

    INPUT:
    - array (ndarray): the array where the window-functions will be applied to
    - window_length (int): the length of the window. If the length is 0 or equal or larger than the length of array, this will be set to the length of the array.
    - window_function: the function
        - none: will apply a rectangular window with 1's for all elements
        - ones: will make a rectangular window with a certain length and pads with zeros.
        - triangular: will make a triangular window with a certain length and pads with zeros
        - gaussian: will make a gaussian window for the full length of array, but will go to zero at around window_length.
        
    CHANGELOG:
    20110909/RB: added zeropadding

    """
    DEBUG.verbose("window function: " + window_function, flag_verbose)

    dim = len(numpy.shape(array))

    # for single dimensions
    if dim == 1:
        # the window function should end up with the same length as the array
        array_length = numpy.shape(array)[0]

        # if it is smaller than the length, make it that length
        if window_length > 0 and window_length < array_length:
            n_max = window_length
            zeros = numpy.zeros(array_length - window_length) 
        else:
            n_max = array_length
            zeros = []

        # the windows
        if window_function == "none":
            window = numpy.ones(array_length)

        elif window_function == "ones":
            window = numpy.concatenate((numpy.ones(n_max).T, zeros)) 

        elif window_function == "triangle":
            window = numpy.concatenate((numpy.linspace(1, 0, n_max).T, zeros))   

        elif window_function == "gaussian":
            window = numpy.exp(-(2.2*numpy.arange(0, array_length)/(n_max))**2)
            
        elif window_function == "experimental": 
            window = numpy.exp(-(2.2*numpy.arange(0, array_length)/(n_max))**2)

        else:
            DEBUG.printError("Unknown window function.", inspect.stack())
            window = numpy.ones(array_length)

        if flag_plot:
            m = numpy.max(array)

            plt.figure()
            plt.plot(array)
            plt.plot(window * m)
            plt.plot(array*window)
            plt.title("window function is scaled")
            plt.show()

        return array * window

    # for higher dimensions
    else:
        DEBUG.printError("Not implemented yet for multiple dimensions.", inspect.stack())
        return 0      



# make FT axis
def make_ft_axis(length, dt, undersampling = 0, normalized_to_period = 0, zero_in_middle = False, flag_verbose = False):
    """
    fourier transform the time axis
    20101204/RB: started
    20130131/RB: now uses speed of light value from Crocodile.Resources.Constants

    INPUT:
    length: amount of samples
    dt: time between samples
    undersampling: 
    normalized_to_period: will normalize the spectrum to 1 for the specified period. If the number is 0, the spectrum will not be normalized.
    zero_in_middle: determines the frequency axes.

    OUTPUT:
    A frequency axis.     

    """
    DEBUG.verbose("make FT axis", flag_verbose)
    if normalized_to_period == 0:   
        resolution = 1 / ( CONST.wavenumberToInvFs * length * dt)
    else:
        resolution = normalized_to_period / (length * dt)

    array = numpy.arange((undersampling)*length/2, (undersampling+1)*length/2)*resolution

    if zero_in_middle == False:
        return numpy.concatenate((array,-numpy.flipud(array)))
    else:
        return numpy.concatenate((-numpy.flipud(array), array))





### CORRELATION STUFF ###


def correlation(array, maxtau = 200, step_type = "tau", flag_normalize = True, flag_verbose = False):
    """
    Calculation of the correlation using the method Jan used.
    The method is slow and probably wrong. Use correlation_fft instead.

    For every iteration, the copied array will be 'rolled' or 'rotated' left by 1 for maxtau times. The copied array will be multiplied with the original array, but only the elements with a certain step between them will be used. The larger the step size, the quicker the method but also the more noisy the result will become.

    INPUT:
    array (ndarray): 1-d array with the data
    maxtau (int): the maximum shift, also the maximum to where the correlation is calculated. This directly affects the speed of the calculation. (with N^2?)
    step_type ("1", "tau"): The step size. With "1" the step size is 1, this will result in a longer calculation but less noise. With "tau" the step size is the current "tau". The calculation will be faster but the result will be noisier.
    flag_normalize (BOOL, True): see note below.

    CHANGELOG:
    20120215: Introduced step_type
    20130131/RB: introduced flag_normalize
    20130204/RB: tested if 'array2 = numpy.roll(array2, -1)' is better nested in the itertools call, but it makes no real change on the speed of the function. No changes made.
    """
    DEBUG.verbose("Correlation Jan-style", flag_verbose)
    
    array = array - numpy.mean(array)

    array2 = numpy.copy(array)

    c = numpy.zeros(maxtau)

    for i in range(0, maxtau):

        array2 = numpy.roll(array2, -1)

        if step_type == "tau":
            step = i+1
        elif step_type == "1":
            step = 1
        else:   
            DEBUG.printWarning("step_type is not recognized, will use 'tau'", inspect.stack())
            step = i+1

        a = list(itertools.islice(array * array2, None, len(array)-i-1, step))

        c[i] = numpy.sum(a) / len(a)

    if flag_normalize:
        return c/c[0]
    else:
        return c



# def correlation_fft(array, flag_normalize = True, flag_verbose = False):
#     """
#     Calculate the autocorrelation using fft.
#     
#     This method was verified using a naive implementation in C. 
#     
#     INPUT:
#     - array (ndarray): the data
#     - flag_normalizae (Bool, True): if True, the starting value of the autocorrelation is 1. If not, it is an absolute value.
#     - flag_verbose (Bool, False): if True, print some debugging stuff
#     
#     OUTPUT:
#     - autocorrelation of array, normalized to the length or to 1, the real values   
#     
#     201202xx/RB: started function
#     20130205/RB: the function now uses an actual Fourier transform
#     20130207/RB: take the first part of the array, not the last part reversed. This was done to agree with Jan's correlation method, but now it seems that one is wrong.
#     20130515/RB: the result is now always divided by the length of the array and it gives the actual absolute value. Added some documentation
#     
#     """
#     
#     DEBUG.verbose("correlation_fft", flag_verbose)
# 
#     # by subtracting the mean, the autocorrelation decays to zero
#     array -= numpy.mean(array)
# 
#     # zeropad to closest 2^n to prevent aliasing
#     l = 2 ** int(1+numpy.log2(len(array) * 2))
# 
#     # calculate autocorrelation
#     s = numpy.fft.fft(array, n=l)
#     r = numpy.fft.ifft(s * numpy.conjugate(s))
# 
#     # normalize to length
#     r = r[:len(array)] / len(array)
# 
#     # return value
#     if flag_normalize:
#         return numpy.real(r/r[0])
#     else:
#         return numpy.real(r)

def correlation_fft(a, v = -1, flag_normalize = True, flag_verbose = False):
    """
    Calculate the autocorrelation using fft.
    
    This method was verified using a naive implementation in C. 
    
    INPUT:
    - a, v (ndarray): the data, 1D array. If v == 1, the autocorrelation of a with a will be calculated. 
    - flag_normalizae (Bool, True): if True, the starting value of the autocorrelation is 1. If not, it is an absolute value.
    - flag_verbose (Bool, False): if True, print some debugging stuff
    
    OUTPUT:
    - autocorrelation of array, normalized to the length or to 1, the real values   
    
    201202xx/RB: started function
    20130205/RB: the function now uses an actual Fourier transform
    20130207/RB: take the first part of the array, not the last part reversed. This was done to agree with Jan's correlation method, but now it seems that one is wrong.
    20130515/RB: the result is now always divided by the length of the array and it gives the actual absolute value. Added some documentation
    20160317/RB: added v, to calculate the correlation between two arrays. 
    
    """
    
#     DEBUG.verbose("correlation_fft", flag_verbose)

    if v == -1:
        v = a[:]

    # by subtracting the mean, the autocorrelation decays to zero
    a -= numpy.mean(a)
    v -= numpy.mean(v)

    l_a = len(a)
    l_v = len(v)
    
    if l_a > l_v:
        l_pad = 2 ** int(1+numpy.log2(l_a * 2))
        l_min = l_v
    else:
        l_pad = 2 ** int(1+numpy.log2(l_v * 2))
        l_min = l_a

    # zeropad to closest 2^n to prevent aliasing
    a = numpy.pad(a, (0, l_pad-l_a), "constant", constant_values = 0)
    v = numpy.pad(v, (0, l_pad-l_v), "constant", constant_values = 0)

    # calculate autocorrelation
    s_a = numpy.fft.fft(a)
    s_v = numpy.conjugate(numpy.fft.fft(v))
    r = numpy.fft.ifft(s_a * s_v)

    # normalize to length
    r = r[:l_min] / l_min

    # return value
    if flag_normalize:
        return numpy.real(r/r[0])
    else:
        return numpy.real(r)

### DERIVATIVE ###



def derivative(x, y):
    """
    20110909/RB: rudimentary method to calculate the derivative
    """

    dx = x[1] - x[0]

    l = len(y)

    x_temp = numpy.zeros(l-2)
    y_temp = numpy.zeros(l-2)

    for i in range(1,l-1):
        x_temp[i-1] = x[i]

        dy = (y[i] - y[i-1] + y[i+1] - y[i]) / 2

        y_temp[i-1] = dy / dx

    return x_temp, y_temp    



























