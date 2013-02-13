from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.Debug as DEBUG
# import Crocodile.Resources.Plotting as PL
import Crocodile.Resources.Functions as FU


def plot_overlap(
    ma, 
    ax,
    la = [],
    x_range = [0,0],
    y_range = [0,-1],
    zlimit = -1,
    contours = 12,
    colors = ["b", "r"],
    x_label = "", 
    y_label = "", 
    title = "", 
    diagonal_line = True, 
    invert_colors = False, 
    ma_linewidth = 1,
    la_linewidth = 2,
    flag_verbose = False):
    
    DEBUG.verbose("contour plot", flag_verbose)

    for i  in range(2):
        # check for correct lengths
        y, x = numpy.shape(ma[i].s)
        if len(ma[i].s_axis[2]) != x and len(ma[i].s_axis[0]) != y:
            DEBUG.printError("The data should have the same shape as the axes, wrong for both axes", inspect.stack())
            return False
        elif len(ma[i].s_axis[2]) != x:
            DEBUG.printError("The data should have the same shape as the axes, wrong for the x-axis", inspect.stack())
            return False
        elif len(ma[i].s_axis[0]) != y:
            DEBUG.printError("The data should have the same shape as the axes, wrong for the y-axis", inspect.stack())  
            return False          
        
        # invert colors
        if invert_colors:
            ma[i].s = ma[i].s
    
    x_axis = ma[0].s_axis[2]
    y_axis = ma[0].s_axis[0]
    
    # determine the range to be plotted
    x_min, x_max, y_min, y_max = FU.find_axes(x_axis, y_axis, x_range, y_range, flag_verbose)
    
    # find the area to be plotted
    x_min_i, x_max_i= FU.find_axes_indices(x_axis, x_min, x_max)
    y_min_i, y_max_i= FU.find_axes_indices(y_axis, y_min, y_max)

    x_axis = x_axis[x_min_i:x_max_i]
    y_axis = y_axis[y_min_i:y_max_i]   
    
    for i in range(2):
     
        # truncate the data, this speeds up the plotting
        data = ma[i].s[y_min_i:y_max_i,x_min_i:x_max_i]

        # now make the actual contours   
        V = FU.make_contours_2d(data, zlimit, contours, flag_verbose)        
        
        # actually plot the thing
        ax.contour(x_axis, y_axis, data, V, linewidths = ma_linewidth, colors = colors[i])
    
    if len(la) == 4:
        ax.plot(la[0],la[1], c = colors[0], lw = la_linewidth)
        ax.plot(la[2],la[3], c = colors[1], lw = la_linewidth)
    
    # the diagonal line
    if diagonal_line:
        ax.plot([x_axis[0]-100,x_axis[-1]+100], [x_axis[0]-100,x_axis[-1]+100], "k", linewidth = ma_linewidth)
    
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
    
    return True