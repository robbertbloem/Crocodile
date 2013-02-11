"""
Script to process data
This script uses Crocodile and PythonTools

COMMAND LINE INPUT:
-a <string>: do activity <string>. See __main__ for details.
-v: more verbose output
-r: reload modules. Mostly for debugging purposes

INPUT (in this script):
- activity: do this activity. See __main__ for details.
- mess_date: the date of the measurements. This is used to construct the file name.
- mess_array: obj_id (unique), base_filename, time stamp, population time. 
Note that the order is different from before. The ordering now is consistent with the naming of the file.
- sub_type_plus/min: 'base_filename' is used to indicate the sub_type. This is the simplest way to select files for adding/subtracting, given that you keep the naming consistent.
- some variables for processing and plotting. Note that plot_x_range and plot_y_range now use the same notation as pe.plot()
- flag_verbose (BOOL): get a load of information 
- pickle_base_name: best to keep it to mess_date, but you can add something here
- pickle_path:
- data_path:
"""



from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os

try:
    import argparse
    flag_parse = True
except:
    flag_parse = False

import numpy
import matplotlib 
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

activity = "plot_all" # "plot_m" # "mess_s_to_sub_s" #"import" #  

mess_date = 20130208

mess_array = [
    ["Aha1",    "aha",  1505,   300],
    ["Buf1",    "buf",  1509,   300],
    ["Aha2",    "aha",  1512,   300],
    ["Buf2",    "buf",  1515,   300],
]

sub_type_plus = "aha"
sub_type_min = "buf"

### PROCESSING VARIABLES ###
zeropad_by = 4
w3_correction = -14

zlimit_array = [-1]
invert_colors = True
contours = 14
plot_x_range = [0,0]
plot_y_range = [0,-1]

### LESS USED VARIABLES ###

flag_verbose = False

pickle_base_name = str(mess_date) 

location = os.uname()[1]
if location == "rbmbp.local":
    path_to_PML3 = "/Volumes/public_hamm/PML3/"
elif location == "pcipdc.uzh.ch":
    path_to_PML3 = "/home/hamm_storage/Hamm_Groupshare/PML3/"
else: # maybe klemens?
    path_to_PML3 = "/home/hamm_storage/Hamm_Groupshare/PML3/"

data_path = path_to_PML3 + "/data/" + str(mess_date) + "/"
# pickle_path = path_to_PML3 + "/analysis/" + str(mess_date) + "/"
pickle_path = "/Volumes/public_hamm/Robbert/python_test/"

##########################################
### NO NEED TO MAKE CHANGES BELOW HERE ###
##########################################

# parse arguments
if flag_parse:
    parser = argparse.ArgumentParser(description='Input from command line')
    parser.add_argument("-a", type = str, default = activity, help = False)
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
    args = parser.parse_args()

    activity = args.a
    flag_verbose = args.verbose
    flag_reload = args.reload

# reload modules
if flag_reload:
    reload(OA)
    reload(PETW)
    reload(PEME)

# some functions
def import_pickle(oa_object, pickle_path_and_filename, flag_verbose = False):
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
        ax[i] = fig.add_subplot(x,y,i+1)
        oa_object.obj_array[i].plot(ax = ax[i], 
            zlimit = zlimit_array[i],
            invert_colors = invert_colors,
            contours = contours,
            x_range = plot_x_range,
            y_range = plot_y_range,
            flag_verbose = flag_verbose)
    plt.show()



if __name__ == "__main__": 

    # time domain - r and f
    # spectrum - s
    # merged - m
    
    pickle_paf_r = pickle_path + "test_r.pickle"
    pickle_paf_s = pickle_path + "test_s.pickle"
    pickle_paf_m = pickle_path + "test_m.pickle"
    
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
            import_pickle(oar, pickle_paf_r, flag_verbose = flag_verbose)

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
            import_pickle(oar, pickle_paf_r, flag_verbose = flag_verbose)

        # calculate absorptive and save
        for m in oar.obj:
            m.zeropad_by = zeropad_by
            m.absorptive(flag_verbose = flag_verbose)
        oar.save_objectarray(pickle_paf_r, flag_verbose = flag_verbose) 
    
    ##################
    # REMOVE R AND F #
    ##################
    if activity == "import" or activity == "remove_r" or activity == "all":
        if not oar:
            oar = OA.objectarray(str(mess_date) + "_r", flag_verbose = flag_verbose)
            import_pickle(oar, pickle_paf_r, flag_verbose = flag_verbose)
        
        # make new oas objectarray, add objects, remove time domain, save
        oas = OA.objectarray(str(mess_date) + "_s", flag_verbose = flag_verbose)
        oas.add_array_with_objects(oar.obj, flag_verbose = flag_verbose)
        for m in oas.obj:
            m.r = [0,0]
            m.f = [0,0]
        oas.save_objectarray(pickle_paf_s, flag_overwrite = True, flag_verbose = flag_verbose)
        
    #########
    # MERGE #
    #########
    if activity == "mess_s_to_sub_s" or activity == "merge" or activity == "all":
        if not oas:
            oas = OA.objectarray(str(mess_date) + "_s", flag_verbose = flag_verbose)
            import_pickle(oas, pickle_paf_s, flag_verbose = flag_verbose)        
        
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
            import_pickle(oas, pickle_paf_s, flag_verbose = flag_verbose) 
        
        # plot all     
        plot_all(oas, zlimit_array = zlimit_array, flag_verbose = flag_verbose)
               



    

























    
