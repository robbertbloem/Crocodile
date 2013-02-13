from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import argparse
import unittest

import numpy

import Crocodile.Pe as PE
import PythonTools.Debug as DEBUG



# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

# add arguments
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
parser.add_argument("-s1", "--skip1", action="store_true", help="Skip testing suite 1: absorptive")
parser.add_argument("-s2", "--skip2", action="store_true", help="Skip testing suite 2: save data")


# process
args = parser.parse_args()

# reload
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)
    



class Test_Pe_absorptive(unittest.TestCase):
    """

    CHANGELOG:
    20130208/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose
        self.mess = PE.pe("Test", flag_verbose = self.flag_verbose)
        
        self.mess.r = [numpy.zeros((10,20)), numpy.zeros((10,20))]
        self.mess.r_axis = [numpy.arange(10), 10, numpy.arange(20)]
        
        self.mess.phase_degrees = 0
        
        self.mess.undersampling = 0

    
    def test_pe_r_1(self):
        """
        Correct input
        """
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(numpy.shape(self.mess.s), (5,20))

    def test_pe_r_2(self):
        """
        r = [0,0], ie the data is not loaded yet
        """
        self.mess.r = [0,0]
        DEBUG.verbose("\nError is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)
        self.assertFalse(res)
        self.assertEqual(self.mess.s, [0])
        
    def test_pe_r_3(self):
        """
        r = int, ie some wrong assignment has been made
        """
        self.mess.r = 0
        DEBUG.verbose("\nError is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)
        self.assertFalse(res)
        self.assertEqual(self.mess.s, [0])

    def test_pe_r_axis_1(self):
        """
        r = [0,okay], ie the data is not loaded yet
        """
        self.mess.r_axis[0] = 0
        DEBUG.verbose("\nWarning is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(self.mess.s_axis, [0,0,0])

    def test_pe_r_axis_2(self):
        """
        r = [okay,0], ie the data is not loaded yet
        """
        self.mess.r_axis[2] = 0
        DEBUG.verbose("\nWarning is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(self.mess.s_axis, [0,0,0])
        
    def test_pe_r_axis_3(self):
        """
        r_axis = int, ie some wrong assignment has been made
        """
        self.mess.r_axis = 0
        DEBUG.verbose("\nWarning is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(self.mess.s_axis, [0,0,0])

    def test_pe_undersampling_false(self):
        """
        undersampling is False, should be treated as 0
        """
        self.mess.undersampling = False
        DEBUG.verbose("\nWarning is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(self.mess.undersampling, 0)

    def test_pe_undersampling_true(self):
        """
        undersampling is True, should be treated as 1
        """
        self.mess.undersampling = True
        DEBUG.verbose("\nWarning is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(self.mess.undersampling, 1)
         
    def test_pe_undersampling_nan(self):
        """
        undersampling is numpy.nan, should return False
        """
        self.mess.undersampling = numpy.nan
        DEBUG.verbose("\nError is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)    
        self.assertFalse(res)

    def test_pe_phase_false(self):
        """
        _phase_degrees is False (override the setter method), should be treated as 0
        """
        self.mess._phase_degrees = False
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(self.mess.phase_degrees, 0)
    
    def test_pe_phase_true(self):
        """
        _phase_degrees is True (override the setter method), should be treated as 1
        """
        self.mess._phase_degrees = True
        res = self.mess.super_absorptive(axes = 0)
        self.assertTrue(res)
        self.assertEqual(self.mess.phase_degrees, 1)
         
    def test_pe_phase_nan(self):
        """
        _phase_degrees is numpy.nan (override the setter method), should return False
        """
        self.mess._phase_degrees = numpy.nan
        DEBUG.verbose("\nError is intentional", True) 
        res = self.mess.super_absorptive(axes = 0)    
        self.assertFalse(res)      




class Test_Pe_save_data(unittest.TestCase):


    def setUp(self):
        self.flag_verbose = args.verbose
        self.mess = PE.pe("Test", flag_verbose = self.flag_verbose)
        
        self.mess.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        self.mess.r_axis = [numpy.arange(10), 10, numpy.arange(20)]
        
        self.mess.s = numpy.ones((10,20))
        self.mess.s_axis = [numpy.arange(10), 10, numpy.arange(20)]
        
        self.mess.base_filename = "test"
        
        self.export_path = "/Users/robbert/Developer/Crocodile/temp/"
        
    def test_save(self):
        print(self.mess)
        self.mess.save_data(self.export_path, s = True, r = True, flag_verbose = self.flag_verbose)






if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_Pe_absorptive)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 1: absorptive", True)

    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_Pe_save_data)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 2: save data", True)


