from __future__ import print_function
from __future__ import division

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.DataClass as DC
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

# add arguments
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
parser.add_argument("-s1", "--skip1", action="store_true", help="Skip testing suite 1: pe_merge")
# parser.add_argument("-s2", "--skip2", action="store_true", help="Skip testing suite 2: check_value_set_key")
# parser.add_argument("-s3", "--skip3", action="store_true", help="Skip testing suite 3: print objects")
# parser.add_argument("-s4", "--skip4", action="store_true", help="Skip testing suite 4: add array with objects")

# process
args = parser.parse_args()

# reload
if args.reload:
    reload(DC)



class Test_dataclass_setters(unittest.TestCase):
    """

    CHANGELOG:
    20130208/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):

        self.flag_verbose = args.verbose

        self.dc = DC.dataclass("test", 3, 2, flag_verbose = self.flag_verbose)



    def test_zeropad_to_correct(self):
        """
        This is correct
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        self.dc.zeropad_to = 40
        self.assertEqual(self.dc._zeropad_to, 40)
        self.assertEqual(self.dc._zeropad_by, 4.0)

    def test_zeropad_to_r_uninit(self):
        """
        without r set, it should give an error
        """
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_to = 40
        self.assertEqual(self.dc._zeropad_to, None)
        self.assertEqual(self.dc._zeropad_by, 1.0)

    def test_zeropad_to_nan(self):
        """
        set zeropad_to to numpy.nan
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_to = numpy.nan
        self.assertEqual(self.dc._zeropad_to, None)
        self.assertEqual(self.dc._zeropad_by, 1.0)


    def test_zeropad_by_correct(self):
        """
        This is correct
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        self.dc.zeropad_by = 4
        self.assertEqual(self.dc._zeropad_to, 40)
        self.assertEqual(self.dc._zeropad_by, 4.0)
    
    def test_zeropad_by_r_uninit(self):
        """
        without r set, it should give an error
        """
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_by = 4
        self.assertEqual(self.dc._zeropad_to, None)
        self.assertEqual(self.dc._zeropad_by, 1.0)
    
    def test_zeropad_by_nan(self):
        """
        set zeropad_by to numpy.nan
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_by = numpy.nan
        self.assertEqual(self.dc._zeropad_to, None)
        self.assertEqual(self.dc._zeropad_by, 1.0)

    def test_zeropad_by_float(self):
        """
        set zeropad_by to float
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        self.dc.zeropad_by = 3.05
        self.assertEqual(self.dc._zeropad_to, 30)
        self.assertEqual(self.dc._zeropad_by, 3.05)

    def test_zeropad_by_read(self):
        """
        The use of zeropad_to is prefered
        reading zeropad_by should give a warning
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        DEBUG.verbose("\nWarning is intentional", True)
        print(self.dc.zeropad_by)







if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_dataclass_setters)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 1: pe init", True)

    # if args.skip2 == False:
    #     suite = unittest.TestLoader().loadTestsFromTestCase(Test_check_value_set_key)
    #     unittest.TextTestRunner(verbosity=1).run(suite)    
    # else:
        # DEBUG.verbose("Skipping suite 2: check_value_set_key", True)