from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import imp

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe_col as PC
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

imp.reload(PC)

global suite_list
suite_list = [
    "Suite 1: Test of init",
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
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_of_init)
        unittest.TextTestRunner(verbosity=1).run(suite) 
    else:
        DEBUG.verbose("Skipping: " + suite_list[0], True)
        



class Test_of_init(unittest.TestCase):
    """

    CHANGELOG:
    20151104/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose    
        

    def test_1(self):
        
        objectname = "test_class"
        data_folder = "/User/Robbert/Bla"
        date = "20010101"
        basename = "foo"
        timestamp = "123456"
        MacOSX = True
        
        x = PC.pe_col(objectname, flag_verbose = self.flag_verbose)

        x.set_file_info(data_folder, date, basename, timestamp)
        
        print(x)
        
        
 
if __name__ == '__main__':

    execute(args)
