from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.Functions as FU
import PythonTools.Debug as DEBUG






def contourplot(data, x_axis, y_axis,
#     ax = False,
#     x_range = [0,0],
#     y_range = [0,-1],
#     zlimit = -1,
#     contours = 12,
#     filled = True,
#     black_contour = True, 
#     x_label = "", 
#     y_label = "", 
#     title = "", 
#     diagonal_line = True, 
#     invert_colors = False, 
#     linewidth = 1,
#     flag_verbose = False,
    **kwargs):
    
    """
    - data, x_axis, y_axis: data and axes
    - ax (bool (False) or matplotlib axes instance): if False, it will make a new figure, otherwise it will use the axes instance, allowing subplots etc.
    - x_label, y_label, title (string, default=''): the labels for the axes. If no label is set, it will use the default. Use 'no_label' or 'no_title' to show no label.
    - x_range, y_range (array with 2 elements, [0,0], [0,-1]): the range to be plotted. Possible cases:
        - [min, max]: plot range min to max
        - [0, 0]: plot the whole range
        - [0, -1]: use the range from the other axis. If both have this, it will plot both axes complete. (ie. it is identical to both having [0,0])
    - zlimit (number or list, -1): the z-range that will be used
        Possible cases:
        zlimit = 0, show all, not don't care about centering around zero
        zlimit = -1, show all, centered around zero
        zlimit = all else, use that, centered around zero
        zlimit = [a,b], plot from a to b
    - contours (int, 16): number of contours to be used
    - invert_colors (BOOL, False): data = -data   
    
    CHANGELOG:
    201108xx/RB: started function
    20130213/RB: moved some things out as separate functions
    20160418/RB: replaced arguments with kwargs
    
    """
    if "flag_verbose" in kwargs:
        flag_verbose = kwargs["flag_verbose"]
    else:
        flag_verbose = False
        
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

    if "invert_colors" in kwargs and kwargs["invert_colors"]:
        data = -data
    
    if "x_range" in kwargs:
        x_range = kwargs["x_range"]
    else:
        x_range = [0,0]

    if "y_range" in kwargs:
        y_range = kwargs["y_range"]
    else:
        y_range = [0,-1]

    # determine the range to be plotted
    x_min, x_max, y_min, y_max = FU.find_axes(x_axis, y_axis, x_range, y_range, flag_verbose)
    
    # find the area to be plotted
    x_min_i, x_max_i = FU.find_axes_indices(x_axis, x_min, x_max)
    y_min_i, y_max_i = FU.find_axes_indices(y_axis, y_min, y_max)
    
    # truncate the data, this speeds up the plotting
    data, x_axis, y_axis = FU.truncate_data(data, x_axis, y_axis, x_min_i, x_max_i, y_min_i, y_max_i)

    if "zlimit" in kwargs:
        zlimit = kwargs["zlimit"]
    else:
        zlimit = -1
        
    if "contours" in kwargs:
        contours = kwargs["contours"]
    else:
        contours = 16

    # now make the actual contours   
    V = FU.make_contours_2d(data, zlimit, contours, flag_verbose)        

    # make sure there is an axis-object
    if "ax" in kwargs:
        ax = kwargs["ax"]
        flag_show = False
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        flag_show = True
#     else:
#         flag_show = False

    if "linewidth" in kwargs:
        linewidth = kwargs["linewidth"]
    else:
        linewidth = 1

    # actually plot the thing
    if "filled" in kwargs and kwargs["filled"] == False:
        pass
    else:
        ax.contourf(x_axis, y_axis, data, V, cmap = FU.rwb_cmap)
        
    if "black_contour" in kwargs and kwargs["black_contour"] == False:
        pass
    else:
        if "filled" in kwargs and kwargs["filled"] == False:
            ax.contour(x_axis, y_axis, data, V, linewidths = linewidth, colors = "k")
        else:
            ax.contour(x_axis, y_axis, data, V, linewidths = linewidth, linestyles = "solid", colors = "k")
#         else:
            
    
    # the diagonal line
    if "diagonal_line" in kwargs and kwargs["diagonal_line"]:
        ax.plot([x_axis[0]-100,x_axis[-1]+100], [x_axis[0]-100,x_axis[-1]+100], "k", linewidth = linewidth)

    # we only want to see a certain part of the spectrum   
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    
    # add some text
    if "x_label" in kwargs and kwargs["x_label"] != "" and kwargs["x_label"]  != "no_label":
        ax.set_xlabel(kwargs["x_label"] )
    
    if "y_label" in kwargs and kwargs["y_label"] != "" and kwargs["y_label"] != "no_label":
        ax.set_ylabel(kwargs["y_label"])
    
    if "title" in kwargs and kwargs["title"] != "":
        ax.set_title(kwargs["title"])    

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







