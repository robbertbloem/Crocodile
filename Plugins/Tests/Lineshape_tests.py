from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import argparse
import unittest
import os

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe_tw as PETW
import Crocodile.Pe_merge as PEME
import Crocodile.Plugins.Lineshape as LS
import PythonTools.ObjectArray as OA
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

# add arguments
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
parser.add_argument("-s1", "--skip1", action="store_true", help="Skip testing suite 1: ")
# parser.add_argument("-s2", "--skip2", action="store_true", help="Skip testing suite 2: check_value_set_key")
# parser.add_argument("-s3", "--skip3", action="store_true", help="Skip testing suite 3: print objects")
# parser.add_argument("-s4", "--skip4", action="store_true", help="Skip testing suite 4: add array with objects")

# process
args = parser.parse_args()

# reload
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    reload(Crocodile.Resources.ReloadCrocodile)
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)
 


def generate_pickle(pickle_path_and_filename, flag_verbose = False):
    
    DEBUG.verbose("Generating pickle", flag_verbose)
    
    aha = PETW.pe_LV("aha", "PI99aha", 1626, 300)
    aha.path = "/Users/robbert/Developer/Crocodile/Tests/Test_resources/PI99aha_1626_T300/"
    aha.import_data()
    aha.zeropad_by = 4
    aha.absorptive()
    aha.r = [0,0]
    
    buf = PETW.pe_LV("buf", "PI99buf", 1617, 300)
    buf.path = "/Users/robbert/Developer/Crocodile/Tests/Test_resources/PI99buf_1617_T300/"
    buf.import_data()
    buf.zeropad_by = 4
    buf.absorptive()  
    buf.r = [0,0]
    
    m = PEME.pe_merge("merge", class_plus = [aha], class_min = [buf])  
    
    oa = OA.objectarray("merge")
    oa.add_object(m)
    oa.save_objectarray(pickle_path_and_filename, flag_overwrite = True)
    
    

   
class Test_Generate_Pickle(unittest.TestCase):
    """

    CHANGELOG:
    20130212/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):

        self.flag_verbose = args.verbose
        self.pickle_path_and_filename = "/Users/robbert/Developer/Crocodile/temp/lineshape_test.pickle"
        if not os.path.exists(self.pickle_path_and_filename):
            generate_pickle(self.pickle_path_and_filename, flag_verbose = self.flag_verbose)
        self.oa = OA.objectarray("object") 
        self.oa.import_db(self.pickle_path_and_filename, flag_verbose = self.flag_verbose)  



    def test_do_stuff(self):
    
        # make lineshape object
        ls = LS.Lineshape(self.oa.obj_array[0], flag_verbose = self.flag_verbose)
        
        # set some values
        # these are the values used for the AHA paper (with some exceptions)
        
        # Double Lorentzian
        ls.dl_x_i = LS.find_indices(ls.mess.s_axis[2], [2040, 2150])
        ls.dl_y_i = LS.find_indices(ls.mess.s_axis[0], [2090, 2132])
        ls.dl_A_in = [10, 2070, 0, 1, 10, 2120, 0, 1] 
        
        # linear fit on points of Double Lorentzian
        ls.l_i = [4, 8] 
        
        # peak heigth
        ls.ph_ble_x_i = LS.find_indices(ls.mess.s_axis[2], [2100, 2125])
        ls.ph_ble_y_i = LS.find_indices(ls.mess.s_axis[0], [2100, 2125])
        
        ls.ph_esa_x_i = LS.find_indices(ls.mess.s_axis[2], [2060, 2095])
        ls.ph_esa_y_i = LS.find_indices(ls.mess.s_axis[0], [2100, 2125])
        
        # peaks along w1 
        ls.w1_peaks_x_i = LS.find_indices(ls.mess.s_axis[2], [2110, 2130])
        ls.w1_peaks_y_i = LS.find_indices(ls.mess.s_axis[0], [2070, 2140])  
        ls.w1_peaks_A_in = [10, 2110, 0, -1]
        
        # plot range
        ls.plot_x_i = LS.find_indices(ls.mess.s_axis[2], [2040, 2150])
        ls.plot_y_i = LS.find_indices(ls.mess.s_axis[0], [2050, 2160]) 
          
        
        # calculate stuff
        ls.fit_double_lorentzian(flag_verbose = self.flag_verbose, flag_plot = False)
        ls.fit_tilt(flag_verbose = self.flag_verbose)
        
        ls.find_peak_heights(flag_verbose = self.flag_verbose)
        
        ls.find_w1_peaks(flag_verbose = self.flag_verbose)
        
        
        

        
        
        
        
        
            
    
if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_Generate_Pickle)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skip generating pickle", True)

    # if args.skip2 == False:
    #     suite = unittest.TestLoader().loadTestsFromTestCase(Test_check_value_set_key)
    #     unittest.TextTestRunner(verbosity=1).run(suite)    
    # else:
    #     DEBUG.verbose("Skipping suite 2: check_value_set_key", True)