"""
Script to process data
This script uses Crocodile and PythonTools

COMMAND LINE INPUT:
-a <string>: do activity <string>. See __main__ for details.
-v: more verbose output
-r: reload modules. Mostly for debugging purposes
-i: plot non-default measurement
-p: print all object ids

INPUT (in this script):
- activity: do this activity. See __main__ for details.
- mess_date: the date of the measurements. This is used to construct the file name.
- mess_array: obj_id (unique), base_filename, time stamp, population time. 
Note that the order is different from before. The ordering now is consistent with the naming of the file.
- sub_type_plus/min: 'base_filename' is used to indicate the sub_type. This is the simplest way to select files for adding/subtracting, given that you keep the naming consistent. It accepts a string or a list with strings.
- some variables for processing and plotting. Note that plot_x_range and plot_y_range now use the same notation as pe.plot()
- flag_verbose (BOOL): get a load of information 
- pickle_base_name: best to keep it to mess_date, but you can add something here
- pickle_path:
- data_path:

ACTIVITIES:
- import_only (1,2): import data
- r_to_s (1,2): calculate absorptive
- remove_r (1,2): remove time domain
- mess_s_to_sub_s/merge (2): data addition/subtraction
- plot_m (2): plot merged 
- plot_all: plot all spectra
- plot: plot a single spectrum, given by plot_index, use -i to set non-default
- print: print all object ids (also -p)

use 'import' to do (1) together
use 'all' to do (2) together

CHANGELOG
20130212/RB: first version
20130214/RB: lots of debugging, new activities 'plot_all', 'print'. 'r_to_s' now does not save s in _r.pickle, this prevents phase corrections to be applied over and over.
20130219/RB: introduced window functions, corrected multiple plots (n_x >= n_y).


"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import os

try:
    import argparse
    flag_parse = True
except:
    flag_parse = False

import numpy
import matplotlib.pyplot as plt

import PythonTools.ObjectArray as OA
import Crocodile.Pe_tw as PETW
import Crocodile.Pe_merge as PEME

# global variables, no need to pass them around all the time
global invert_colors
global contours
global plot_x_range
global plot_y_range

##########################################
### NO NEED TO MAKE CHANGES ABOVE HERE ###
##########################################

### INPUT ###

activity = "plot_all" # "plot_m" # "merge" # "import" #  

mess_date = 20010101

mess_array = [
    ["N3_5_1",  "N3_5", 1712,   300], # c1
    ["N3_5_2",  "N3_5", 1713,   300], # c2
    ["N3_5_3",  "N3_5", 1714,   300], # c4
    ["h2o_1",   "h2o",  1636,   300], # c3
    ["h2o_2",   "h2o",  1638,   300], # c3
    ["h2o_3",   "h2o",  1640,   300], # c4
]

sub_type_plus = ["N3_5"]
sub_type_min = ["h2o"]

### PROCESSING VARIABLES ###
window_function = "none" # "gaussian"
zeropad_by = 4
w3_correction = -14
# use int to apply for all, use list to apply to specific measurements
phase_correction = 0

zlimit_array = [-1]
invert_colors = True
contours = 14
plot_x_range = [1990,2100]
plot_y_range = [0,-1]

### LESS USED VARIABLES ###

plot_index = 0

flag_verbose = False

pickle_base_name = str(mess_date)

location = os.uname()[1]
if location == "rbmbp.local":
    path_to_PML3 = "/Volumes/public_hamm/PML3/"
elif location == "pcipdc.uzh.ch":
    path_to_PML3 = "/home/hamm_storage/Hamm_Groupshare/PML3/"
else: # maybe klemens?
    path_to_PML3 = "/home/hamm_storage/Hamm_Groupshare/PML3/"

data_path = path_to_PML3 + "data/" + str(mess_date) + "/"
pickle_path = path_to_PML3 + "analysis/" + str(mess_date) + "/"

##########################################
### NO NEED TO MAKE CHANGES BELOW HERE ###
##########################################

# parse arguments
if flag_parse:
    parser = argparse.ArgumentParser(description='Input from command line')
    parser.add_argument("-a", type = str, default = activity, help = False)
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
    parser.add_argument("-p", "--print_object_ids", action="store_true", help="Print object ids")
    parser.add_argument("-i", type = int, default = plot_index, help="Plot single measurement")
    args = parser.parse_args()

    activity = args.a
    plot_index = args.i
    flag_verbose = args.verbose
    flag_reload = args.reload

# reload modules
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)

# some functions
def import_pickle(oa_object, pickle_path_and_filename, mess_array = [], flag_verbose = False):
    if len(mess_array) > 0:
        obj_id_array = numpy.array(mess_array)[:,0]
        oa_object.load_objectarray(pickle_path_and_filename, obj_id_array, flag_verbose = flag_verbose)
    else:
        oa_object.import_db(pickle_path_and_filename, flag_verbose = flag_verbose)
    
def find_subplot_size(i, flag_verbose = False):
    s = numpy.sqrt(i)
    x = numpy.ceil(s)
    y = numpy.floor(s)
    if x*y < i:
        y = y + 1
    return x, y

def plot_all(oa_object, zlimit_array = [-1], flag_verbose = False):  

    n = len(oa_object.obj_array)

    if len(zlimit_array) < n:
        zlimit_array = numpy.concatenate((zlimit_array, numpy.zeros(n - len(zlimit_array))-1))

    ax = [0] * n
    x, y = find_subplot_size(n) 
    fig = plt.figure()
    for i in range(n):
        ax[i] = fig.add_subplot(y,x,i+1)
        oa_object.obj_array[i].plot(ax = ax[i], 
            zlimit = zlimit_array[i],
            invert_colors = invert_colors,
            contours = contours,
            x_range = plot_x_range,
            y_range = plot_y_range,
            flag_verbose = flag_verbose,
            title = oa_object.obj_array[i].obj_id)
    plt.show()



if __name__ == "__main__": 

    # time domain - r and f
    # spectrum - s
    # merged - m
    
    pickle_paf_r = pickle_path + pickle_base_name + "_r.pickle"
    pickle_paf_s = pickle_path + pickle_base_name + "_s.pickle"
    pickle_paf_m = pickle_path + pickle_base_name + "_m.pickle"
    
    oar = False 
    oas = False 
    oam = False 
    
    ###############
    # IMPORT DATA #
    ###############
    if activity == "import" or activity == "import_only" or activity == "all":
        # if oar already exists, don't bother to re-import it
        if not oar:
            oar = OA.objectarray(str(mess_date) + "_r", flag_verbose = flag_verbose)
            import_pickle(oar, pickle_paf_r, mess_array, flag_verbose = flag_verbose)

        # import the data, then save
        for m in mess_array:
            if m[0] not in oar.obj_id_array:
                mess = PETW.pe_LV(m[0], m[1], m[2], m[3], flag_verbose = flag_verbose)
                mess.path = data_path + mess.base_filename + "/"
                mess.import_data(flag_verbose = flag_verbose)
                mess.r_correction[2] = w3_correction
                oar.add_object(mess, flag_verbose = flag_verbose) 

        oar.save_objectarray(pickle_paf_r, flag_verbose = flag_verbose) 
      
    ########################
    # CALCULATE ABSORPTIVE #
    ########################
    if activity == "import" or activity == "r_to_s" or activity == "all":
        if not oar:
            oar = OA.objectarray(str(mess_date) + "_r", flag_verbose = flag_verbose)
            import_pickle(oar, pickle_paf_r, mess_array, flag_verbose = flag_verbose)

        # calculate absorptive and save
        for i in range(len(oar.obj)):
            if type(phase_correction) == int:
                oar.obj[i].phase_degrees += phase_correction
            elif i < len(phase_correction): 
                oar.obj[i].phase_degrees += phase_correction[i]
            oar.obj[i].zeropad_by = zeropad_by
            oar.obj[i].absorptive(window_function = window_function, flag_verbose = flag_verbose)
        
        for m in oar.obj:
            m.r = [0,0]
            m.f = [0,0]
        oar.save_objectarray(pickle_paf_s, flag_overwrite = True, flag_verbose = flag_verbose)
        oar = False
        
    #########
    # MERGE #
    #########
    if activity == "mess_s_to_sub_s" or activity == "merge" or activity == "all":
        if len(sub_type_plus) > 0 and len(sub_type_min) > 0:
        
            if not oas:
                oas = OA.objectarray(str(mess_date) + "_s", flag_verbose = flag_verbose)
                import_pickle(oas, pickle_paf_s, mess_array, flag_verbose = flag_verbose)        
            
            # select object with sub_types
            class_plus = oas.objects_with_sub_type(sub_type_plus, flag_verbose = flag_verbose)
            class_min = oas.objects_with_sub_type(sub_type_min, flag_verbose = flag_verbose)        
            mess = PEME.pe_merge(str(mess_date), class_plus, class_min, flag_verbose = flag_verbose)
            
            # make a merged object
            oam = OA.objectarray(str(mess_date) + "_m", flag_verbose = flag_verbose)
            oam.add_object(mess, flag_verbose = flag_verbose)
            oam.save_objectarray(pickle_paf_m, flag_overwrite = True, flag_verbose = flag_verbose)

    #######################
    # PLOT MERGED SPECTRA #
    #######################
    if activity == "plot_m" or activity == "all":
        if not oam:
            oam = OA.objectarray(str(mess_date) + "_m", flag_verbose = flag_verbose)
            import_pickle(oam, pickle_paf_m, flag_verbose = flag_verbose)   
        
        # plot all - usually only one
        plot_all(oam, zlimit_array = zlimit_array, flag_verbose = flag_verbose)
    
    ####################
    # PLOT ALL SPECTRA #
    ####################
    if activity == "plot_all":
        if not oas:
            oas = OA.objectarray(str(mess_date) + "_s", flag_verbose = flag_verbose)
            import_pickle(oas, pickle_paf_s, mess_array, flag_verbose = flag_verbose) 

        # plot all     
        plot_all(oas, zlimit_array = zlimit_array, flag_verbose = flag_verbose)

    ##########################
    # PLOT SPECIFIC SPECTRUM #
    ##########################
    # use -i n to plot index n, default is 0
    if activity == "plot":
        if not oas:
            oas = OA.objectarray(str(mess_date) + "_s", flag_verbose = flag_verbose)
            import_pickle(oas, pickle_paf_s, mess_array, flag_verbose = flag_verbose) 

        oas.obj_array[plot_index].plot(zlimit = zlimit_array[0],
            invert_colors = invert_colors,
            contours = contours,
            x_range = plot_x_range,
            y_range = plot_y_range,
            flag_verbose = flag_verbose,
            title = oas.obj_array[plot_index].obj_id)

                    
    #####################
    # PRINT ALL OBJECTS #
    #####################
    if activity == "print" or args.print_object_ids:
        if not oar:
            oar = OA.objectarray(str(mess_date) + "_r", flag_verbose = flag_verbose)
            import_pickle(oar, pickle_paf_r, mess_array, flag_verbose = flag_verbose) 
        oar.print_object_ids(flag_verbose = flag_verbose)

        if not oas:
            oas = OA.objectarray(str(mess_date) + "_s", flag_verbose = flag_verbose)
            import_pickle(oas, pickle_paf_s, mess_array, flag_verbose = flag_verbose) 
        oas.print_object_ids(flag_verbose = flag_verbose)

        if not oam:
            oam = OA.objectarray(str(mess_date) + "_m", flag_verbose = flag_verbose)
            import_pickle(oas, pickle_paf_m, flag_verbose = flag_verbose) 
        oam.print_object_ids(flag_verbose = flag_verbose)



        
        



