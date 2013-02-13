# for compatibility with Python 3
# in python 2, you can write 
# >> print "string"
# in python 3, and with the __future__ thingy, it becomes:
# >> print("string")
from __future__ import print_function
# in Python 2: 
# 2.0/3.0 = 0, because it is modulo or something
# 2.0//3.0 = 1.5, // is a normal division
# in Python 3, and with the __future__ thingy, it is the other way around. 
from __future__ import division

from __future__ import absolute_import

# unicode_literals break all python 2.x code and with python 3.3 this __future__ is not needed anymore
# from __future__ import unicode_literals

# needed 'reload' for python 3
import imp

import numpy
import matplotlib 
import matplotlib.pyplot as plt

# for arrays with objects: adding, saving, importing
import PythonTools.ObjectArray 
# usually imported as OA
# import PythonTools.ObjectArray as OA

# photon echo stuff
import Crocodile.Pe
# import Crocodile.Pe as PE

# photon echo time-frequency stuff
import Crocodile.Pe_tw
# import Crocodile.Pe as PETW

# photon echo merge stuff
import Crocodile.Pe_merge
# import Crocodile.Pe as PEME

# if True, you can follow what the script is doing
flag_verbose = False

#######################
# IMPORT OBJECT ARRAY #
#######################

# create an object array
name = ""   # a name
oa = PythonTools.ObjectArray.objectarray(name, flag_verbose = flag_verbose)

# use import_db() to import with unknown order
oa.import_db(pickle_path_and_filename, flag_verbose = flag_verbose)

# use load_objectarray() to import a pickle, with the order as in obj_id_array
pickle_path_and_filename = ""           # path of the pickle
obj_id_array = ["obj_id1", "obj_id2"]   # array with the unique object ids
oa.load_objectarray(pickle_path_and_filename, obj_id_array, flag_verbose = flag_verbose)

# note that from mess_array we only need the first column
mess_array = [
    ["obj_id1", "sample",   1200, 300],
    ["obj_id2", "solvent",  1210, 300],
    ["obj_id3", "sample",   1220, 300],
]
obj_id_array = numpy.array(mess_array)[:,0]

#######################
# ADD TO OBJECT ARRAY #
#######################

# add a single object: obj
# the obj should have an obj_id, this limits it to Crocodile objects
obj = Crocodile.Pe("pe_object", flag_verbose = flag_verbose)
oa.add_object(obj, flag_verbose = flag_verbose)

# use add_array_with_objects() to add an array with objects, each having an obj_id
obj_array = [obj, obj]
oa.add_array_with_objects(obj_array, flag_verbose = flag_verbose)


#####################
# SAVE OBJECT ARRAY #
#####################

# use save_objectarray() to save 
flag_overwrite = False      # False: add data, True: overwrite data
oa.save_objectarray(pickle_path_and_filename, flag_overwrite = flag_overwrite, flag_verbose = flag_verbose)

######################
# PRINT OBJECT ARRAY #
######################

# print the object in an ordered way
# when index is a valid index in oa, it will print only that object, otherwise it will print all
oa.print_objects(index = -1, flag_verbose = flag_verbose)

# print only the index, obj_id and sub_type
oa.print_object_ids(flag_verbose = flag_verbose)

#############################
# OBJECTS FROM OBJECT ARRAY #
#############################

# return a list with the indices of objects with sub_type. It has to be string or an array with strings
sub_type = "subtype"
sub_type = ["sub_type_1", "sub_type_2"]
list_with_indices = oa.list_objects_with_sub_type(sub_type, flag_verbose = flag_verbose)

# return a list with objects
# for those who paid attention, this list can indeed be used for add_array_with_objects()
list_with_objects = oa.objects_with_sub_type(sub_type, flag_verbose = flag_verbose)



################################################
# INTERACTING WITH OBJECTS IN THE OBJECT ARRAY #
################################################

# let's say we added some Pe_LV objects to the objectarray oa 

# the phase_degrees for object at index 2 can be found like this:
print(oa.obj_array[2].phase_degrees)

# if we want to set zeropad_by and then calculate absorptive for all objects
for i in range(len(oa.obj_array)):
    oa.obj_array[i].zeropad_by = 2
    oa.obj_array[i].absorptive(flag_verbose = flag_verbose)

# or:
for obj in oa.obj_array:
    obj.zeropad_by = 2
    obj.absorptive(flag_verbose = flag_verbose)




######################################
# FIND INDEX FOR A CERTAIN FREQUENCY #
######################################

# see Crocodile.Resources.Functions.find_axes_indices() to find the indices for a range from val_min to val_max
# see Crocodile.Resources.Functions.truncate_data() to truncate data to more managable portions

# use numpy.where(condition) to find the indices where condition is True. For example
# numpy.where(axis > value)


############
# PLOTTING #
############

# again, there are some Pe_LV objects in objectarray oa 
# plotting has been made a little bit more complex and a lot more versatile

# for a single plot, you can still use
oa.obj_array[0].plot(flag_verbose = flag_verbose)

# to plot an array with measurements, use the following (for 5 plots)
n = 5               # number of plots we want to make
x = 3               # number of columns
y = 2               # number of rows
fig = plt.figure()  # make plot object
ax = [0] * n        # initialize the axes-array
for i in range(n):
    # add a new subplot
    ax[i] = fig.add_subplot(y,x,i+1)
    # plot the subplot
    oa.obj_array[i].plot(ax = ax[i], flag_verbose = flag_verbose)

# oops, forgot to set the x_lim <- advantage is we can change this now
for i in range(n):
    ax[i].set_xlim(2000,2100)

# show the whole thing
plt.show()

#################################
# MORE COOL STUFF WITH PLOTTING #
#################################
# or: the reason to uses axes

### MORE CONTROL OF SUB PLOTS ###

# instead of 
ax[i] = fig.add_subplot(y,x,i+1)
# use add_axes()
# now you can precisely determine the positions of plots
# add_axes() takes a tupple with (left_edge, bottom_edge, width, height)
# these values are given as fraction of the total width and height

axes_inch_per_unit = 0.3        # inch per unit

axes_x_units = 7                # horizontal nr of units
axes_y_units = 8                # vertical nr of units

axes_l = 1.7 / axes_x_units     # left edge of subplot, in units

axes_w = 4 / axes_x_units       # width of subplot, in units

axes_t = 5.5 / axes_y_units     # lower edge of top row, in units
axes_m = 3.5 / axes_y_units     # lower edge of middle row, in units
axes_b = 1.5 / axes_y_units     # lower edge of bottom row, in units

axes_h = 2 / axes_y_units       # height of subplot, in units

axes_coords = [                 # array with coords for all subplots
    (axes_l, axes_t, axes_w, axes_h),
    (axes_l, axes_m, axes_w, axes_h),
    (axes_l, axes_b, axes_w, axes_h),
]

fig_width = axes_inch_per_unit * axes_x_units   # width of figure, in units
fig_height = axes_inch_per_unit * axes_y_units  # height of figure, in units

fig = plt.figure(figsize = (fig_width, fig_height), dpi = 300) # init figure

ax = [0] * 3                    # init axes array

for i in range(3):
    ax[i] = fig.add_axes(axes_coords[i])    # init axes


### ADD A SECOND Y-AXIS ###

ax[0] = fig.add_subplot(1,1,1)
ax[1] = ax[0].twinx()


### LOOK OF THE PLOT ###

x = [1,2,3]     # some data
y = [1,2,3]
color = "g"     # color of line "b", "g", "r", "k" (black)
linewidth = 1   # set linewidth
linestyle = "-" # set linestyle ":", "-."
ax[i].plot(x, y, c = color, linestyle = linestyle, linewidth = linewidth)


### LOOK OF THE TICKS AND LABELS ###

fontsize = 8    # fontsize
# set ticks at positions (in data-units)
ax[i].set_yticks([0, 5, 10, 15, 20])
# set tick labels for the ticks
ax[i].set_yticklabels(["0", "", "10", "", "20"], fontsize = fontsize)
# or
ax[i].set_yticklabels(["bla", "", "di", "", "bla"], fontsize = fontsize)

# change the size of the ticks "x" or "y" axis, which = "both", "major" or "minor"
ax[i].tick_params("x", which = "both", length = 2)

# change the width of (a selection of) the edges. Use 0 to show no edges. 
# iteritems() fails in Python 3
for loc, spine in ax[i].spines.iteritems():
    if loc in ['left', 'right', 'bottom', 'top']:
        spine.set_linewidth(0.5)

### LATEX TEXT ###

# Matplotlib can use latex to make fancy labels. The r"" indicates a raw string: Python will not parse it ("\n" will make a newline, r"\n" wil print \n. Use $$ as you would in LaTex. 
text = r"$cm^{-1}$" 



### PRINT TEXT IN COLUMNS ###

# see for more details: http://docs.python.org/2/library/string.html#format-specification-mini-language
a = 5
b = "string"
c = "another string"
print("{0:3d} {1:10s} {2:10s}".format(a, b, c))














