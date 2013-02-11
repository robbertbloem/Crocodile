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




def make_contours_2d(data, zlimit = 0, contours = 21, flag_verbose = False):
    """
    zlimit = 0, show all, not don't care about centering around zero
    zlimit = -1, show all, centered around zero
    zlimit = all else, use that, centered around zero
    zlimit = [a,b], plot from a to b
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




def contourplot(data, x_axis, y_axis,
    ax = False,
    x_range = [0,0],
    y_range = [0,-1],
    zlimit = -1,
    contours = 12,
    filled = True,
    black_contour = True, 
    x_label = "", 
    y_label = "", 
    title = "", 
    diagonal_line = True, 
    invert_colors = False, 
    linewidth = 1,
    flag_verbose = False):
    
    DEBUG.verbose("contour plot", flag_verbose)


    y, x = numpy.shape(data)
    if len(x_axis) != x and len(y_axis) != y:
        DEBUG.printError("The data should have the same shape as the axes, wrong for both axes", inspect.stack())
        return False
    elif len(x_axis) != x:
        DEBUG.printError("The data should have the same shape as the axes, wrong for the x-axis", inspect.stack())
        return False
    elif len(y_axis) != y:
        DEBUG.printError("The data should have the same shape as the axes, wrong for the y-axis", inspect.stack())  
        return False          

    if invert_colors:
        data = -data

    # determine the range to be plotted
    x_min, x_max, y_min, y_max = find_axes(x_axis, y_axis, x_range, y_range, flag_verbose)
    
    # make the contours
    # first find the area to be plotted
    # not the most elegant way I guess
    try:
        y_min_i = numpy.where(y_axis < y_min)[0][-1]
    except: 
        y_min_i = 0
    
    try:
        y_max_i = numpy.where(y_axis > y_max)[0][0] + 1
    except: 
        y_max_i = len(y_axis)
    
    try:
        x_min_i = numpy.where(x_axis < x_min)[0][-1]
    except: 
        x_min_i = 0
    
    try:
        x_max_i = numpy.where(x_axis > x_max)[0][0] + 1
    except: 
        x_max_i = len(x_axis)
        
    # truncate the data, this speeds up the plotting
    data = data[y_min_i:y_max_i,x_min_i:x_max_i]
    x_axis = x_axis[x_min_i:x_max_i]
    y_axis = y_axis[y_min_i:y_max_i]

    # now make the actual contours   
    V = make_contours_2d(data, zlimit, contours, flag_verbose)        

    # make sure there is an axis-object
    if ax == False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        flag_show = True
    else:
        flag_show = False

    # actually plot the thing
    if filled:
        ax.contourf(x_axis, y_axis, data, V, cmap = rwb_cmap)
    if black_contour:
        if filled:
            ax.contour(x_axis, y_axis, data, V, linewidths = linewidth, linestyles = "solid", colors = "k")
        else:
            ax.contour(x_axis, y_axis, data, V, linewidths = linewidth, colors = "k")
    
    # the diagonal line
    if diagonal_line:
        ax.plot([x_axis[0]-100,x_axis[-1]+100], [x_axis[0]-100,x_axis[-1]+100], "k", linewidth = linewidth)

    # we only want to see a certain part of the spectrum   
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    # add some text
    if x_label != "" and x_label != "no_label":
        ax.set_xlabel(x_label)
    
    if y_label != "" and y_label != "no_label":
        ax.set_ylabel(y_label)
    
    if title != "":
        ax.set_title(title)    

    if flag_show:
        plt.show()

    return True


def linear(data, x_axis, 
    x_range = [0, 0], 
    y_range = [0, 0], 
    ax = False, 
    x_label = "", 
    y_label = "", 
    title = "", 
    legend = "", 
    plot_real = True,
    flag_verbose = False):

    DEBUG.verbose("linear plot", flag_verbose)

    if plot_real:
        data = numpy.real(data)

    # make the x-axis
    if x_range == [0, 0]:
        x_min = x_axis[0]
        x_max = x_axis[-1]
    else:
        x_min = x_range[0]
        x_max = x_range[1]

    # select the appropriate data range
    x_min_i = numpy.where(x_axis > x_min)[0][0]
    x_max_i = numpy.where(x_axis < x_max)[0][-1] 

    x_axis = x_axis[x_min_i:x_max_i]
    data = data[x_min_i:x_max_i]     

    # make sure there is an axis-object
    if ax == False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        flag_show = True
    else:
        flag_show = False

    # the actual plot
    ax.plot(x_axis, data, label = legend)

    ax.set_xlim(x_min, x_max)
    if y_range != [0,-1]:
        ax.set_ylim(y_range[0], y_range[1])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.set_title(title)

    if legend != "":
        plt.legend()

    if flag_show:
        plt.show()
        
    return True







