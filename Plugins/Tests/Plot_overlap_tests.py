from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
# from __future__ import unicode_literals

import argparse
import unittest
import os

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe_tw as PETW
import Crocodile.Pe_merge as PEME
import Crocodile.Plugins.Plot_overlap as PO
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




class Test_Plot_overlap(unittest.TestCase):
    """

    CHANGELOG:
    20130212/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):

        self.flag_verbose = args.verbose

        self.flag_verbose = args.verbose
        self.pickle_path_and_filename = "/Users/robbert/Developer/Crocodile/temp/lineshape_test.pickle"
        if not os.path.exists(self.pickle_path_and_filename):
            generate_pickle(self.pickle_path_and_filename, flag_verbose = self.flag_verbose)
        
        self.oa1 = OA.objectarray("object1") 
        self.oa1.import_db(self.pickle_path_and_filename, flag_verbose = self.flag_verbose)  
        self.oa2 = OA.objectarray("object2") 
        self.oa2.import_db(self.pickle_path_and_filename, flag_verbose = self.flag_verbose)  

        self.la = [[2075,2125],[2100,2125],[2075,2125],[2075,2100]]


    def test_do_stuff(self):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ma = [self.oa1.obj_array[0], self.oa2.obj_array[0]]
        # make a small change so we can see the difference
        ma[1].s = numpy.roll(ma[1].s, 1, 0)
        
        
        
        PO.plot_overlap(ma, ax, self.la, ma_linewidth = 2, la_linewidth = 4, contours = 16, x_range = [2050,2150])
        
        plt.show()













if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_Plot_overlap)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skip generating pickle", True)

    # if args.skip2 == False:
    #     suite = unittest.TestLoader().loadTestsFromTestCase(Test_check_value_set_key)
    #     unittest.TextTestRunner(verbosity=1).run(suite)    
    # else:
    #     DEBUG.verbose("Skipping suite 2: check_value_set_key", True)