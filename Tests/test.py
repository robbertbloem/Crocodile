from __future__ import print_function
from __future__ import division

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

reload(OA)
reload(PETW)
reload(PEME)


def import_pickle(oa_object, pickle_path_and_filename):
    oa_object.import_db(pickle_path_and_filename)

def import_data(oa_object, m_array, data_path):
    for m in m_array:
        if m[0] not in oa_object.obj_id_array:
            mess = PETW.pe_LV(m[0], m[1], m[3], m[2])
            mess.path = data_path + mess.base_filename + "/"
            mess.import_data()
            oa_object.add_object(mess)    

def r_to_s(oa_object):
    for m in oa_object.obj:
        m.absorptive()

def remove_r(oa_object):
    for m in oa_object.obj:
        m.r = [0,0]
        m.f = [0,0]

def find_subplot_size(i):
    s = numpy.sqrt(i)
    x = numpy.ceil(s)
    y = numpy.floor(s)
    if x*y < i:
        y = y + 1
    return x, y

def plot_all(oa_object):  
    n = len(oa_object.obj_array)
    ax = [0] * n
    x, y = find_subplot_size(n) 
    fig = plt.figure()
    for i in range(n):
        ax[i] = fig.add_subplot(x,y,i+1)
        oa_object.obj_array[i].plot(ax = ax[i])
    plt.show()


flag_verbose = False

activity = "plot_all" # "plot_m" # "mess_s_to_sub_s" #"import" #  

mess_date = 20130208

mess_array = [
    ["Aha1", "aha", 300, 1505, 1],
    ["Buf1", "buf", 300, 1509, 1],
    ["Aha2", "aha", 300, 1512, 1],
    ["Buf2", "buf", 300, 1515, 1],
]

pickle_base_name = str(mess_date) 

pickle_path = "/Users/robbert/Developer/Crocodile/temp/"
data_path = "/Volumes/public_hamm/PML3/data/" + str(mess_date) + "/"


if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description='Input from command line')
    parser.add_argument("-a", type = str, default = activity, help = False)
    args = parser.parse_args()
    activity = args.a
      
    pickle_paf_r = pickle_path + "test_r.pickle"
    pickle_paf_s = pickle_path + "test_s.pickle"
    pickle_paf_m = pickle_path + "test_m.pickle"
    
    oar = False # time domain - r and f
    oas = False # spectrum - s
    oam = False # merged - m
    
    if activity == "import" or activity == "import_only" or activity == "all":
        # if oar already exists, don't bother to re-import it
        if not oar:
            oar = OA.objectarray(str(mess_date) + "_r")
            import_pickle(oar, pickle_paf_r)

        # import the data, then save
        import_data(oar, mess_array, data_path)  
        oar.save_objectarray(pickle_paf_r) 
      
    if activity == "import" or activity == "r_to_s" or activity == "all":
        if not oar:
            oar = OA.objectarray(str(mess_date) + "_r")
            import_pickle(oar, pickle_paf_r)

        # calculate absorptive and save
        r_to_s(oar)
        oar.save_objectarray(pickle_paf_r) 
    
    if activity == "import" or activity == "remove_r" or activity == "all":
        if not oar:
            oar = OA.objectarray(str(mess_date) + "_r")
            import_pickle(oar, pickle_paf_r)
        
        # make new oas objectarray, add objects, remove time domain, save
        oas = OA.objectarray(str(mess_date) + "_s")
        oas.add_array_with_objects(oar.obj)
        remove_r(oas)
        oas.save_objectarray(pickle_paf_s, flag_overwrite = True)
        
    if activity == "mess_s_to_sub_s" or activity == "to_ss" or activity == "all":
        if not oas:
            oas = OA.objectarray(str(mess_date) + "_s")
            import_pickle(oas, pickle_paf_s)        
        
        # select object with sub_types
        class_plus = oas.objects_with_sub_type("aha")
        class_min = oas.objects_with_sub_type("buf")        
        mess = PEME.pe_merge(str(mess_date), class_plus, class_min)
        
        # make a merged object
        oam = OA.objectarray(str(mess_date) + "_m")
        oam.add_object(mess)
        oam.save_objectarray(pickle_paf_m, flag_overwrite = True)

    if activity == "plot_m" or activity == "all":
        if not oam:
            oam = OA.objectarray(str(mess_date) + "_m")
            import_pickle(oam, pickle_paf_m)   
        
        # there can be more objects in the merged file... better solution?
        oam.obj[0].plot()    
    
    # not included in "all"
    if activity == "plot_all":
        if not oas:
            oas = OA.objectarray(str(mess_date) + "_s")
            import_pickle(oas, pickle_paf_s)      
        
        # plot all
        plot_all(oas)
               


    

























    
