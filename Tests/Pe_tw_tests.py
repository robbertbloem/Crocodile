from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe as PE
import Crocodile.Pe_tw as PETW
import Crocodile.Resources.IOMethods as IOM
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

# add arguments
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
parser.add_argument("-s1", "--skip1", action="store_true", help="Skip testing suite 1: importing LabView")
parser.add_argument("-s2", "--skip2", action="store_true", help="Skip testing suite 2: init Pe_tw")
# parser.add_argument("-s3", "--skip3", action="store_true", help="Skip testing suite 3: print objects")
# parser.add_argument("-s4", "--skip4", action="store_true", help="Skip testing suite 4: add array with objects")

# process
args = parser.parse_args()

# reload
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)





class Test_Pe_LV_import(unittest.TestCase):
    """

    CHANGELOG:
    20130211/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose
        self.mess = PETW.pe_LV("Test", "PI99aha", 1626, 300, flag_verbose = self.flag_verbose)

    def test_initialized_values(self):
        self.assertEqual(self.mess.r_axis[1], 300)
        self.assertEqual(self.mess.time_stamp, "1626")
        self.assertEqual(self.mess.base_filename, "PI99aha_1626_T300")
        self.assertEqual(self.mess.path, "PI99aha_1626_T300/")
        self.assertEqual(self.mess.sub_type, "PI99aha")
        self.assertEqual(self.mess.undersampling, 0)
        self.assertEqual(self.mess.n_pixels, 32)
        self.assertEqual(self.mess.data_type_version, 0)
        self.assertEqual(self.mess.n_shots, 0)
        self.assertEqual(self.mess.n_steps, 0)
        self.assertEqual(self.mess.n_scans, 0)
        self.assertTrue(numpy.all(self.mess.reference))
        self.assertEqual(numpy.shape(self.mess.reference), (32,))
       
    def test_import_data_correct(self):
        """
        Data is imported. Returns True. Check if phase is indeed imported.
        """
        self.mess.path = "Test_resources/PI99aha_1626_T300/"
        res = self.mess.import_data(flag_verbose = self.flag_verbose)
        self.assertTrue(res)
        self.assertAlmostEqual(self.mess.phase_degrees, 319.705649)

    def test_import_data_incorrect_path(self):
        """
        Path is incorrect, no data imported. Returns False
        """
        self.mess.path = "Test_resources/FolderDoesNotExist/"
        DEBUG.verbose("\nError that directory can't found is intentional", True)
        res = self.mess.import_data(flag_verbose = self.flag_verbose)
        self.assertFalse(res)
        
    def test_import_data_wrong_file_format(self):
        """
        LV_file_format.666 should not be recognized.
        """
        self.mess.path = "Test_resources/petw_test_folder_1/"
        DEBUG.verbose("\nError about unknown file format is intentional", True)
        res = self.mess.import_data(flag_verbose = self.flag_verbose)
        self.assertFalse(res)

    def test_import_data_no_data(self):
        """
        folder exists, LV_file_format.1 exists and is correct, but contains no other data
        """
        self.mess.path = "Test_resources/petw_test_folder_2/"
        DEBUG.verbose("\nError about importing is intentional", True)
        res = self.mess.import_data(flag_verbose = self.flag_verbose)
        self.assertFalse(res)

    def test_plot(self):
        """
        Data is imported. Returns True. Check if phase is indeed imported.
        """
        self.mess.path = "Test_resources/PI99aha_1626_T300/"
        res = self.mess.import_data(flag_verbose = self.flag_verbose)
        self.mess.absorptive()
        self.mess.plot()
        # self.mess.plot_R()
        # self.mess.plot_NR()
        # self.mess.plot_T()
        self.assertTrue(res)
        # self.assertAlmostEqual(self.mess.phase_degrees, 319.705649)


class Test_Pe_tw_init(unittest.TestCase):
    """

    CHANGELOG:
    20130211/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose
        self.mess = PETW.pe_tw("Test", "Test", 1626, 300, undersampling = 4, flag_verbose = self.flag_verbose)

    def test_initialized_values(self):
        self.assertEqual(self.mess.r_axis[1], 300)
        self.assertEqual(self.mess.time_stamp, "1626")
        self.assertEqual(self.mess.base_filename, "Test_1626_T300")
        self.assertEqual(self.mess.path, "Test_1626_T300/")
        self.assertEqual(self.mess.sub_type, "Test")
        self.assertEqual(self.mess.undersampling, 4)
        self.assertEqual(self.mess.n_pixels, 32)
        self.assertEqual(self.mess.data_type_version, 0)
        self.assertEqual(self.mess.n_shots, 0)
        self.assertEqual(self.mess.n_steps, 0)
        self.assertEqual(self.mess.n_scans, 0)
        self.assertTrue(numpy.all(self.mess.reference))
        self.assertEqual(numpy.shape(self.mess.reference), (32,))
        



















if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_Pe_LV_import)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 1: importing LabView", True)
        
        
    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_Pe_tw_init)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 2: init Pe_tw", True)











