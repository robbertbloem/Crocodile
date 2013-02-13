from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.Debug as DEBUG

cdict = {'red':   [(0.0,  0.0, 0.0),(0.475,  1.0, 1.0),(0.525,  1.0, 1.0),
            (1.0,  1.0, 1.0)],
         'green': [(0.0,  0.0, 0.0),(0.475,  1.0, 1.0),(0.525,  1.0, 1.0),
            (1.0,  0.0, 0.0)],
         'blue':  [(0.0,  1.0, 1.0),(0.475,  1.0, 1.0),(0.525,  1.0, 1.0),
            (1.0,  0.0, 0.0)]
        }
rwb_cmap = matplotlib.colors.LinearSegmentedColormap('rwb_colormap', cdict, 256)


def find_axes(x_axis, y_axis, x_range, y_range, flag_verbose = False):
    """
    croc.Plotting.find_axis

    INPUT:
    - x_axis, y_axis (ndarray): the axes
    - x_range, y_range (array): the ranges to be plotted
        Possible cases:
        - [min, max]: plot range min to max
        - [0, 0]: plot the whole range
        - [0, -1]: use the range from the other axis. If both have this, it will plot both axes complete. (ie. it is identical to both having [0,0])

    OUTPUT:
    - xmin, xmax, ymin, ymax: a value

    CHANGELOG
    - 20110910 RB: init

    """
    DEBUG.verbose("find axes", flag_verbose)

    # x_min and x_max are determined by y
    if x_range[0] == 0 and x_range[1] == -1:

        # they are both dependent on each other, plot everything
        if y_range[0] == 0 and y_range[1] == -1:
            x_min = x_axis[0]
            x_max = x_axis[-1]
            y_min = y_axis[0]
            y_max = y_axis[-1]     

        # take the whole y-axis, the same for the x-axis
        # it may extend beyond the range of the data 
        elif y_range[0] == 0 and y_range[1] == 0:
            x_min = y_axis[0]
            x_max = y_axis[-1]
            y_min = y_axis[0]
            y_max = y_axis[-1]         

        # take part of the y-axis, the same for the x-axis
        # it may extend beyond the range of the data
        else:
            x_min = y_range[0]
            x_max = y_range[1]
            y_min = y_range[0]
            y_max = y_range[1]  

    else:
        # plot the whole x-axis
        if x_range[0] == 0 and x_range[1] == 0:
            x_min = x_axis[0]
            x_max = x_axis[-1]
        # plot only part of the x-axis
        else:
            x_min = x_range[0]
            x_max = x_range[1]

        # plot the whole y-axis
        if y_range[0] == 0 and y_range[1] == 0:
            y_min = y_axis[0]
            y_max = y_axis[-1]
        # use the range of the x-axis
        elif y_range[0] == 0 and y_range[1] == -1:
            y_min = x_min
            y_max = x_max
        # use the specified range
        else:
            y_min = y_range[0]
            y_max = y_range[1]        
    if flag_verbose:
        DEBUG.verbose("  axes found: x_min: " + str(x_min) + ", x_max: " + str(x_max) + ", y_min: " + str(y_min) + ", y_max: " + str(y_max), flag_verbose)

    return x_min, x_max, y_min, y_max    


def find_axes_indices(axis, val_min, val_max, flag_verbose = False):
    """
    Find the indices of values val_min and val_max in axis, such that val_min and val_max are included. 
    Because the intended use is to slice an array data[val_min_i:val_max_i], 1 is added to val_max_i.
    If the val_min or val_max exceed the axis, then the index is 0 or -1, respectively.
    
    BEHAVIOR:
    axis = [3,4,5,6]
    if val_min == 4.5: val_min_i = 1
    if val_min == 5: val_min_i = 1
    if val_min == 1: val_min_i = 0
    
    if val_max = 4.5: val_max_i = 3: axis[0:3] = [3,4,5]
    if val_max = 4: val_max_i = 2: axis[0:2] = [3,4]
    if val_max = 10: val_max_i = -1: axis[0:-1] = [3,4,5]
    
    CHANGELOG:
    201108xx/RB: originated in contourplot
    20130213/RB: moved to separate function. Better handling of edge cases.
    
    """
    if val_min > val_max:
        DEBUG.printWarning("val_min > val_max, will give strange result.", inspect.stack())
    
    temp = numpy.where(axis < val_min)
    if len(temp[0]) == 0:
        val_min_i = 0
    else:
        val_min_i = temp[0][-1]

    temp = numpy.where(axis > val_max)
    if len(temp[0]) == 0:
        val_max_i = -1
    else:
        val_max_i = temp[0][0] + 1
        if val_max_i == len(axis):
            val_max_i = -1

    return val_min_i, val_max_i


def make_contours_2d(data, zlimit = 0, contours = 21, flag_verbose = False):
    """
    zlimit = 0, show all, not don't care about centering around zero
    zlimit = -1, show all, centered around zero
    zlimit = all else, use that, centered around zero
    zlimit = [a,b], plot from a to b
    
    CHANGELOG:
    201108xx/RB: started in Plotting module
    20130213/RB: moved to Functions module
    """

    DEBUG.verbose("make contours 2D", flag_verbose)

    if zlimit == 0:
        ma = numpy.amax(data)
        mi = numpy.amin(data)
        if flag_verbose:
            DEBUG.verbose("  zlimit = 0, min: " + str(mi) + ", max: " + str(ma), flag_verbose)
        return numpy.linspace(mi, ma, num=contours)

    elif zlimit == -1:
        ma = numpy.amax(data)
        mi = numpy.amin(data)
        if flag_verbose:
            DEBUG.verbose("  zlimit = 0, min: " + str(mi) + ", max: " + str(ma), flag_verbose)

        if abs(mi) > abs(ma):
            ma = abs(mi)
        else:
            ma = abs(ma)
        return numpy.linspace(-ma, ma, num=contours) 

    elif type(zlimit) == list:
        DEBUG.verbose("  zlimit == list", flag_verbose)
        return numpy.linspace(zlimit[0], zlimit[1], num=contours)   

    else:
        DEBUG.verbose("  zlimit is linspace", flag_verbose)
        return numpy.linspace(-abs(zlimit), abs(zlimit), num=contours) 


def truncate_data(data, x_axis, y_axis, x_min_i, x_max_i, y_min_i, y_max_i, flag_verbose = False):
    """
    Truncate data in a non-default way. 
    x_axis = [1,2,3,4,5]
    x_min_i = 0
    if x_max_i == 2: x_axis = [1,2]
    if x_max_i == 4: x_axis = [1,2,3,4]
    if x_max_i == -1: x_axis = [1,2,3,4,5]
    
    CHANGELOG:
    20130213: was in Plotting.contourplot function, moved to Functions module, now handles x_max_i == -1
    """
    
    DEBUG.verbose("Truncate data", flag_verbose)

    if x_max_i == -1:
        DEBUG.verbose("  x_max == -1", flag_verbose)
        data = data[:,x_min_i:]
        x_axis = x_axis[x_min_i:]
    else:
        DEBUG.verbose("  x_max != -1", flag_verbose)
        data = data[:,x_min_i:x_max_i]
        x_axis = x_axis[x_min_i:x_max_i]        

    if y_max_i == -1:
        DEBUG.verbose("  y_max == -1", flag_verbose)
        data = data[y_min_i:,:]
        y_axis = y_axis[y_min_i:]   
    else:
        DEBUG.verbose("  y_max != -1", flag_verbose)
        data = data[y_min_i:y_max_i,:]
        y_axis = y_axis[y_min_i:y_max_i]         
        
    return data, x_axis, y_axis

