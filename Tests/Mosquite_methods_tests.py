from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Mosquito_methods as MM
import PythonTools.Debug as DEBUG

plt.close("all")

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

imp.reload(MM)

global suite_list
suite_list = [
    "Suite 1: Show shots",
    "Suite 2: Show spectrum",
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])
parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])
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
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_show_shots)
        unittest.TextTestRunner(verbosity=1).run(suite) 
    else:
        DEBUG.verbose("Skipping: " + suite_list[0], True)
        
    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_show_spectrum)
        unittest.TextTestRunner(verbosity=1).run(suite) 
    else:
        DEBUG.verbose("Skipping: " + suite_list[1], True)

class Test_show_spectrum(unittest.TestCase):
    """

    CHANGELOG:
    20151104/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose  

        
    def test_f1(self):

        file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "show_spectrum_f1",
            "timestamp": "141014",
            "extension": ".csv",
        }  
        objectname = "test A"   
        A = MM.show_spectrum(objectname, flag_verbose = self.flag_verbose)
        A.set_file_dict(file_dict_A["data_folder"], file_dict_A["date"], file_dict_A["basename"], file_dict_A["timestamp"])
        A.import_data()
#         A.make_plot()
        print(A.file_format)

    def test_pp(self):

        file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "show_spectrum_pp",
            "timestamp": "143628",
            "extension": ".csv",
        }  
        objectname = "test A"   
        A = MM.show_spectrum(objectname, flag_verbose = self.flag_verbose)
        A.set_file_dict(file_dict_A["data_folder"], file_dict_A["date"], file_dict_A["basename"], file_dict_A["timestamp"])
        A.import_data()
#         A.make_plot()
        print(A.file_format)

    def test_pp_sepds(self):

        file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "show_spectrum_pp_sepds",
            "timestamp": "145025",
            "extension": ".csv",
        }  
        objectname = "test A"   
        A = MM.show_spectrum(objectname, flag_verbose = self.flag_verbose)
        A.set_file_dict(file_dict_A["data_folder"], file_dict_A["date"], file_dict_A["basename"], file_dict_A["timestamp"])
        A.import_data()
#         A.make_plot()
        print(A.file_format)    

class Test_show_shots(unittest.TestCase):
    """

    CHANGELOG:
    20151104/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose  

        
    def test_f1(self):

        file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "show_shots_f1",
            "timestamp": "140921",
            "extension": ".csv",
        }  
        objectname = "test A"   
        A = MM.show_shots(objectname, flag_verbose = self.flag_verbose)
        A.set_file_dict(file_dict_A["data_folder"], file_dict_A["date"], file_dict_A["basename"], file_dict_A["timestamp"])
        A.import_data()
        A.make_plot()
        print(A.file_format)

        
    def test_pp(self):

        file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "show_shots_pp",
            "timestamp": "143606",
            "extension": ".csv",
        }  
        objectname = "test A"   
        A = MM.show_shots(objectname, flag_verbose = self.flag_verbose)
        A.set_file_dict(file_dict_A["data_folder"], file_dict_A["date"], file_dict_A["basename"], file_dict_A["timestamp"])
        A.import_data()
        A.make_plot()
        print(A.file_format)
    
 
if __name__ == '__main__':

    execute(args)
   