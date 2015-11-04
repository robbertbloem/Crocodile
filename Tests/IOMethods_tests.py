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
    "Suite 4: find number of datastates",
    "Suite 5: import supporting files", 
    "Suite 6: slow modulation", 
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])
parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])
parser.add_argument("-s3", "--skip3", action = "store_true", help = suite_list[2])
parser.add_argument("-s4", "--skip4", action = "store_true", help = suite_list[3])
parser.add_argument("-s5", "--skip5", action = "store_true", help = suite_list[4])
parser.add_argument("-s6", "--skip6", action = "store_true", help = suite_list[5])

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
    else:
        DEBUG.verbose("Skipping: " + suite_list[0], True)
        
    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_LV_fileformat)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping: " + suite_list[1], True)
        
    if args.skip3 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_number_of_scans)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping: " + suite_list[2], True)
        
    if args.skip4 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_number_of_datastates)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping: " + suite_list[3], True)
        
    if args.skip5 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_import_supporting_files_LV3)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping: " + suite_list[4], True)

    if args.skip6 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_import_slowModulation)
        unittest.TextTestRunner(verbosity=1).run(suite)  
    else:
        DEBUG.verbose("Skipping: " + suite_list[5], True)


class Test_import_slowModulation(unittest.TestCase):
    """

    CHANGELOG:
    20151103/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        
        self.flag_verbose = args.verbose 
        
        self.file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "A",
            "timestamp": "123456",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/A_123456/",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/A_123456/A_123456", 
        }
        
        self.file_dict_B = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "B",
            "timestamp": "234506",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/B_234506/",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/B_234506/B_234506", 
        }
        
        self.fileformat = 3
        

    def test_import_slowModulation_A(self):
        
        sm, sm_names, n_sm = IOM.import_slow_modulation(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)

        print(sm, sm_names, n_sm)
        
        self.assertEqual(sm, numpy.array([]))
        self.assertEqual(sm_names, numpy.array([]))
        self.assertEqual(n_sm, 1)

        
        
        
        
#     def test_lines(self):
#     
#         tests = [
#             {"input": ["1", "Mod A: [NaN]", "Mod B: [NaN]", "Mod C: [NaN]", "Spectrometer: [NaN]"], "sm":"", "sm_names":"", "n_sm":""}
#             {"input": ["2", "Mod A: [NaN,NaN]", "Mod B: [0,1]", "Mod C: [NaN,NaN]", "Spectrometer: [NaN,NaN]"], "sm":"", "sm_names":"", "n_sm":""}
#         ]
#         
#         for test in tests:
#             sm, sm_names, n_sm = IOM.extract_slow_modulation(test["input"], flag_verbose = self.flag_verbose)
#             self.assertEqual(sm, test["sm"])
#             self.assertEqual(sm_names, test["sm_names"])
#             self.assertEqual(n_sm, test["n_sm"])


class Test_import_supporting_files_LV3(unittest.TestCase):
    """

    CHANGELOG:
    20151103/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        
        self.flag_verbose = args.verbose  

        self.file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "A",
            "timestamp": "123456",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/A_123456/",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/A_123456/A_123456", 
        }
        
        self.file_dict_B = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "B",
            "timestamp": "234506",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/B_234506/",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/B_234506/B_234506", 
        }
        
        self.fileformat = 3
        
        
    def test_import_bins(self):
        """
        Bins and times have the same direction
        """
        
        t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index = IOM.import_bins(file_dict = self.file_dict_A, fileformat = self.fileformat, flag_verbose = self.flag_verbose)
        
        self.assertFalse(bin_sign)
        
        self.assertEqual(t1_bins[0], 0)
        self.assertEqual(t1_bins[-1], 854)
        self.assertEqual(numpy.shape(t1_bins)[0], 855)
        
        self.assertEqual(t1_fs[0], -301.843484)
        self.assertEqual(t1_fs[-1], 1500.774246)
        self.assertEqual(numpy.shape(t1_fs)[0], 855)
        
        self.assertEqual(n_t1_bins, 855)
        self.assertEqual(n_t1_fs, 712)
        self.assertEqual(t1_zero_index, 143)


    def test_import_bins_inverted(self):
        """
        Bins and times have the opposite direction
        """       
        t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index = IOM.import_bins(file_dict = self.file_dict_B, fileformat = self.fileformat, flag_verbose = self.flag_verbose)
        
        self.assertTrue(bin_sign)
        
        self.assertEqual(t1_bins[0], 0)
        self.assertEqual(t1_bins[-1], 854)
        self.assertEqual(numpy.shape(t1_bins)[0], 855)
        
        self.assertEqual(t1_fs[0], -301.843484)
        self.assertEqual(t1_fs[-1], 1500.774246)
        self.assertEqual(numpy.shape(t1_fs)[0], 855)
        
        self.assertEqual(n_t1_bins, 855)
        self.assertEqual(n_t1_fs, 712)
        self.assertEqual(t1_zero_index, 143)        


    def test_import_nspectra(self):
        
        n_sp = IOM.import_nspectra(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_sp, 1)


    def test_import_ndatastates_1(self):
        """
        1 datastate
        """
        n_ds = IOM.import_ndatastates(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_ds, 1)
        

    def test_import_ndatastates_2(self):
        """
        2 datastates
        """
        n_ds = IOM.import_ndatastates(self.file_dict_B, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_ds, 2)
    
    
    def test_import_wavenumbers(self):
        w3_axis_wn, n_w3 = IOM.import_wavenumbers(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_w3, 32)
        self.assertEqual(w3_axis_wn[0], 1969.110381)
        self.assertEqual(w3_axis_wn[-1], 2144.561837)


    def test_import_delays_A(self):    
        """
        Multiple delays
        """
        de, n_de = IOM.import_delays(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_de, 5)
        self.assertEqual(de[0], 0)
        self.assertEqual(de[-1], 400)

    def test_import_delays_B(self):  
        """
        Single delay
        """
        de, n_de = IOM.import_delays(self.file_dict_B, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_de, 1)
        self.assertEqual(de[0], 300)
        self.assertEqual(de[-1], 300)


    def test_import_wavelengths(self):
        wl, n_wl = IOM.import_wavelengths(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)    
        self.assertEqual(n_wl, 32)
        self.assertEqual(wl[0], 5078.435469)
        self.assertEqual(wl[-1], 4662.957173)

    def test_import_nshots(self):
        n_sh = IOM.import_nshots(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_sh, 100)

    def test_import_spectraAndDatastates_A(self):
        spds, n_sp, n_ds = IOM.import_spectraAndDatastates(self.file_dict_A, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_ds, 1)
        self.assertEqual(n_sp, 1)
        self.assertTrue(numpy.all(spds == numpy.array([[0,1]])))
        
    def test_import_spectraAndDatastates_B(self):
        spds, n_sp, n_ds = IOM.import_spectraAndDatastates(self.file_dict_B, self.fileformat, flag_verbose = self.flag_verbose)
        self.assertEqual(n_ds, 2)
        self.assertEqual(n_sp, 1)
        self.assertTrue(numpy.all(spds == numpy.array([[0,1],[0,-1]])))



class Test_find_number_of_datastates(unittest.TestCase):
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
            # simplest case
            {"input": [            
"azide_intensity_175815_count_ds0_sp0_sm0_de0_du0_0.csv",
], 
"result": 1, 
"base_filename":"scan", 
"extension":".csv", 
"test": "equal"},
            # fail test
            {"input": [            
"azide_intensity_175815_count_ds0_sp0_sm0_de0_du0_0.csv",
], 
"result": 2, 
"base_filename":"scan", 
"extension":".csv", 
"test": "not equal"},
            # simple case
            {"input": [            
"azide_intensity_175815_count_ds0_sp0_sm0_de0_du0_0.csv",
"azide_intensity_175815_count_ds0_sp0_sm0_de0_du0_1.csv",
"azide_intensity_175815_count_ds1_sp0_sm0_de0_du0_0.csv",
"azide_intensity_175815_count_ds1_sp0_sm0_de0_du0_1.csv",
"azide_intensity_175815_intensity_ds0_sp0_sm0_de0_du0_0.csv",
"azide_intensity_175815_intensity_ds0_sp0_sm0_de0_du0_1.csv",
"azide_intensity_175815_intensity_ds1_sp0_sm0_de0_du0_0.csv",
"azide_intensity_175815_intensity_ds1_sp0_sm0_de0_du0_1.csv",
], 
"result": 2, 
"base_filename":"scan", 
"extension":".csv", 
"test": "equal"},
            # long filenames
            {"input": [            
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds0_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds0_sp0_sm0_de0_du0_1.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds1_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds1_sp0_sm0_de0_du0_1.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds0_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds0_sp0_sm0_de0_du0_1.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds1_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds1_sp0_sm0_de0_du0_1.csv",
], 
"result": 2, 
"base_filename":"scan", 
"extension":".csv", 
"test": "equal"},
            # long filenames, higher numbers
            {"input": [            
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds10_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds10_sp0_sm0_de0_du0_1.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds11_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_count_ds11_sp0_sm0_de0_du0_1.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds10_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds10_sp0_sm0_de0_du0_1.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds11_sp0_sm0_de0_du0_0.csv",
"/Users/robbert/Desktop/20151006/azide_intensity_175815/azide_intensity_175815_intensity_ds11_sp0_sm0_de0_du0_1.csv",
], 
"result": 12, 
"base_filename":"scan", 
"extension":".csv", 
"test": "equal"},
        ]
        
        for i in range(len(test_input)):
            res = IOM.find_number_of_datastates(base_folder = "", flag_verbose = self.flag_verbose, test_input = test_input[i]["input"])
            if test_input[i]["test"] == "equal":
                self.assertEqual(res, test_input[i]["result"])
            elif test_input[i]["test"] == "not equal":
                self.assertNotEqual(res, test_input[i]["result"])
            else:
                print("Invalid test")



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
            res = IOM.find_number_of_scans(base_folder = "", base_filename = test_input[i]["base_filename"], extension = test_input[i]["extension"], flag_verbose = self.flag_verbose, test_input = test_input[i]["input"])
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
            res = IOM.find_LV_fileformat(base_folder, flag_verbose = self.flag_verbose, test_input = test_input[i]["input"])
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
            res = IOM.find_LV_fileformat(base_folder, flag_verbose = self.flag_verbose, test_input = test_input[i]["input"])
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
        var_out, n_var_out = IOM.check_and_make_list(var_in, flag_verbose = self.flag_verbose)        
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)

    def test_ndarray_1d(self):
        """
        In is same as out.
        """  
        n_var_in = (4)
        var_in = numpy.zeros(n_var_in)
        var_out, n_var_out = IOM.check_and_make_list(var_in, flag_verbose = self.flag_verbose)      
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)

    def test_ndarray_1d_2d(self):
        """
        In is same as out.
        """  
        n_var_in = (4,1)
        var_in = numpy.zeros(n_var_in)
        var_out, n_var_out = IOM.check_and_make_list(var_in, flag_verbose = self.flag_verbose)    
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)

    def test_ndarray_1d_2d_inv(self):
        """
        In is same as out.
        """  
        n_var_in = (1,4)
        var_in = numpy.zeros(n_var_in)
        var_out, n_var_out = IOM.check_and_make_list(var_in, flag_verbose = self.flag_verbose)      
        self.assertEqual(n_var_in, n_var_out)
        res = numpy.all(var_in == var_out)
        self.assertTrue(res)


if __name__ == '__main__':

    execute(args)
