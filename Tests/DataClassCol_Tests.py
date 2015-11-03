from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.DataClassCol as DCC
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

imp.reload(DCC)

global suite_list
suite_list = [
    "Suite 1: Set data folder",
    "Suite 2: LV file format",
    "Suite 3: find number of scans",
    "Suite 4: find number of datastates",
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])
# parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])
# parser.add_argument("-s3", "--skip3", action = "store_true", help = suite_list[2])
# parser.add_argument("-s4", "--skip4", action = "store_true", help = suite_list[3])

# process
args = parser.parse_args()

# reload
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)


def execute(args):

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_setting_file_stuff)
        unittest.TextTestRunner(verbosity=1).run(suite) 
    else:
        DEBUG.verbose("Skipping: " + suite_list[0], True)
        
#     if args.skip2 == False:
#         suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_LV_fileformat)
#         unittest.TextTestRunner(verbosity=1).run(suite)  
#     else:
#         DEBUG.verbose("Skipping: " + suite_list[1], True)
#         
#     if args.skip3 == False:
#         suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_number_of_scans)
#         unittest.TextTestRunner(verbosity=1).run(suite)  
#     else:
#         DEBUG.verbose("Skipping: " + suite_list[2], True)
#         
#     if args.skip4 == False:
#         suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_number_of_datastates)
#         unittest.TextTestRunner(verbosity=1).run(suite)  
#     else:
#         DEBUG.verbose("Skipping: " + suite_list[3], True)
        


class Test_setting_file_stuff(unittest.TestCase):
    """

    CHANGELOG:
    20151103/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose    
        

    def test_construct_file_paths(self):
    
        test_input = [
            {"data_folder": "data/",
            "date": "20150102",
            "basename": "azide",
            "timestamp": "123456",
            "extension": ".csv",
            "res_base_folder": "data/20150102/azide_123456/",
            "res_base_filename": "data/20150102/azide_123456/azide_123456", 
            "test": "equal"}, 
        ]
        
        for i in range(len(test_input)):
            x = DCC.dataclass("test class")
            
            x.data_folder = test_input[i]["data_folder"]
            x.date = test_input[i]["date"]
            x.basename = test_input[i]["basename"]
            x.timestamp = test_input[i]["timestamp"]

            res = x._file_dict["base_folder"]
            print(res)
            if test_input[i]["test"] == "equal":
                self.assertEqual(res, test_input[i]["res_base_folder"])
            elif test_input[i]["test"] == "not equal":
                self.assertNotEqual(res, test_input[i]["res_base_folder"])
            else:
                print("Invalid test")        

            res = x._file_dict["base_filename"]
            print(res)
            if test_input[i]["test"] == "equal":
                self.assertEqual(res, test_input[i]["res_base_filename"])
            elif test_input[i]["test"] == "not equal":
                self.assertNotEqual(res, test_input[i]["res_base_filename"])
            else:
                print("Invalid test")  



    def test_data_folder_mac(self):
            
        test_input = [
            {"input": "fiets", "result": "fiets/", "test": "equal"},
            {"input": "fiets/", "result": "fiets/", "test": "equal"},
            {"input": "c:/data/", "result": "c:/data/", "test": "equal"},
            {"input": "/Users/robbert/Desktop/", "result": "/Users/robbert/Desktop/", "test": "equal"},
        ]
        
        for i in range(len(test_input)):
            x = DCC.dataclass("test class")
            x.data_folder = test_input[i]["input"]
            res = x.data_folder
            if test_input[i]["test"] == "equal":
                self.assertEqual(res, test_input[i]["result"])
            elif test_input[i]["test"] == "not equal":
                self.assertNotEqual(res, test_input[i]["result"])
            else:
                print("Invalid test")




        
 
#     def test_data_folder_windows(self):
#         x = DCC.dataclass("test class")
#         x.data_folder = "fiets"
#         print(x._file_dict)
 
if __name__ == '__main__':

    execute(args)
