from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.IOMethods as IOM
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

imp.reload(IOM)

global suite_list
suite_list = [
    "Suite 1: Zeropad setters/getters",
    "Suite 2: LV file format",
    "Suite 3: find number of scans",
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])
parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])
parser.add_argument("-s3", "--skip3", action = "store_true", help = suite_list[2])

# process
args = parser.parse_args()

# reload
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)


def execute(args):

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_check_and_make_list)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    elif args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_LV_fileformat)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    elif args.skip3 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_number_of_scans)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping :" + suite_list[0], True)



class Test_find_number_of_scans(unittest.TestCase):
    """

    CHANGELOG:
    20151103/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose    

    def test_individual_strings(self):
        base_folder = ""
        test_input = [
            # simple case
            {"input": ["scan_0.csv","scan_1.csv","scan_2.csv","scan_3.csv", "fiets_0.csv"], "result": 4, "base_filename":"scan", "extension":".csv", "test": "equal"},
            # unequal case
            {"input": ["scan_0.csv","scan_1.csv","scan_2.csv","scan_3.csv", "fiets_0.csv"], "result": 6, "base_filename":"scan", "extension":".csv", "test": "not equal"},
            # high number
            {"input": ["scan_0.csv","scan_1.csv","scan_2.csv","scan_200.csv", "fiets_0.csv"], "result": 201, "base_filename":"scan", "extension":".csv", "test": "equal"},
            # more numbers
            {"input": [     "/Users/robbert/Dropbox/Amsterdam/20151029/find_t0_crystal_f1_141458/find_t0_crystal_f1_141458_wavenumbers_0.csv",            "/Users/robbert/Dropbox/Amsterdam/20151029/find_t0_crystal_f1_141458/find_t0_crystal_f1_141458_wavenumbers_1.csv",
"/Users/robbert/Dropbox/Amsterdam/20151029/find_t0_crystal_f1_141458/find_t0_crystal_f1_141458_wavenumbers_2.csv", 
"/Users/robbert/Dropbox/Amsterdam/20151029/find_t0_crystal_f1_141458/find_t0_crystal_f1_141458_wavenumbers_3.csv",
"/Users/robbert/Dropbox/Amsterdam/20151029/find_t0_crystal_f1_141458/find_t0_crystal_f1_141458_delays_3.csv"], 
"result": 4, 
"base_filename":"/Users/robbert/Dropbox/Amsterdam/20151029/find_t0_crystal_f1_141458/find_t0_crystal_f1_141458_wavenumbers", 
"extension":".csv", 
"test": "equal"},
        ]
        
        for i in range(len(test_input)):
            res = IOM.find_number_of_scans(base_filename = test_input[i]["base_filename"], extension = test_input[i]["extension"], verbose = False, test_input = test_input[i]["input"])
            if test_input[i]["test"] == "equal":
                self.assertEqual(res, test_input[i]["result"])
            elif test_input[i]["test"] == "not equal":
                self.assertNotEqual(res, test_input[i]["result"])
            else:
                print("Invalid test")



class Test_find_LV_fileformat(unittest.TestCase):
    """

    CHANGELOG:
    20151103/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose    

    def test_individual_strings(self):
        base_folder = ""
        test_input = [
            {"input": "fiets_LV_fileformat.3", "result": 34, "test": "not equal"},
            {"input": "fiets_LV_fileformat.3", "result": 3, "test": "equal"},
            {"input": "fiets_LV_fileformat.345", "result": 345, "test": "equal"},
            {"input": "fiets_123_LV_fileformat.345", "result": 345, "test": "equal"},
            {"input": "fiets_123LV_fileformat.345", "result": 345, "test": "equal"},
            {"input": "/Users/robbert/Dropbox/Amsterdam/20151029/find_t0_crystal_f1_141458/find_t0_crystal_f1_141458_LV_fileformat.3", "result": 3, "test": "equal"},
        ]
        
        for i in range(len(test_input)):
            res = IOM.find_LV_fileformat(base_folder, verbose = False, test_input = test_input[i]["input"])
            if test_input[i]["test"] == "equal":
                self.assertEqual(res, test_input[i]["result"])
            elif test_input[i]["test"] == "not equal":
                self.assertNotEqual(res, test_input[i]["result"])
            else:
                print("Invalid test")


    def test_list_of_strings(self):
        base_folder = ""
        test_input = [
            {"input": ["fiets_123456.csv", "fiets_123456_LV_fileformat.345"], "result": 345, "test": "equal"},
        ]
        
        for i in range(len(test_input)):
            res = IOM.find_LV_fileformat(base_folder, verbose = False, test_input = test_input[i]["input"])
            if test_input[i]["test"] == "equal":
                self.assertEqual(res, test_input[i]["result"])
            elif test_input[i]["test"] == "not equal":
                self.assertNotEqual(res, test_input[i]["result"])
            else:
                print("Invalid test")

class Test_check_and_make_list(unittest.TestCase):
    """

    CHANGELOG:
    20151028/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose


    def test_ndarray_2d(self):      
        """
        In is same as out.
        """  
        n_var_in = (4,8)
        var_in = numpy.zeros(n_var_in)
        var_out, n_var_out = IOM.check_and_make_list(var_in, verbose = self.flag_verbose)        
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)

    def test_ndarray_1d(self):
        """
        In is same as out.
        """  
        n_var_in = (4)
        var_in = numpy.zeros(n_var_in)
        var_out, n_var_out = IOM.check_and_make_list(var_in, verbose = self.flag_verbose)      
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)

    def test_ndarray_1d_2d(self):
        """
        In is same as out.
        """  
        n_var_in = (4,1)
        var_in = numpy.zeros(n_var_in)
        var_out, n_var_out = IOM.check_and_make_list(var_in, verbose = self.flag_verbose)    
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)

    def test_ndarray_1d_2d_inv(self):
        """
        In is same as out.
        """  
        n_var_in = (1,4)
        var_in = numpy.zeros(n_var_in)
        var_out, n_var_out = IOM.check_and_make_list(var_in, verbose = self.flag_verbose)      
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)


if __name__ == '__main__':

    execute(args)
