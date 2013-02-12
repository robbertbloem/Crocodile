from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.DataClass as DC
import PythonTools.Debug as DEBUG
# import PythonTools.ReloadAll as RA




# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

# add arguments
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
parser.add_argument("-s1", "--skip1", action="store_true", help="Skip testing suite 1: zeropad setters/getters")
parser.add_argument("-s2", "--skip2", action="store_true", help="Skip testing suite 2: phase setters/getters")

# process
args = parser.parse_args()

# reload
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)



class Test_dataclass_zeropad(unittest.TestCase):
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

    ##################
    ### ZEROPAD_TO ###
    ##################    
    def test_zeropad_to_correct(self):
        """
        This is correct
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        self.dc.zeropad_to = 40
        self.assertEqual(self.dc.zeropad_to, 40)
        self.assertEqual(self.dc.zeropad_by, 4.0)

    def test_zeropad_to_r_uninit(self):
        """
        without r set, it should give an error
        """
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_to = 40
        self.assertEqual(self.dc.zeropad_to, None)
        self.assertEqual(self.dc.zeropad_by, 1.0)

    def test_zeropad_to_nan(self):
        """
        set zeropad_to to numpy.nan
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_to = numpy.nan
        self.assertEqual(self.dc.zeropad_to, None)
        self.assertEqual(self.dc.zeropad_by, 1.0)

    ##################
    ### ZEROPAD_BY ###
    ##################  
    def test_zeropad_by_correct(self):
        """
        This is correct
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        self.dc.zeropad_by = 4
        self.assertEqual(self.dc.zeropad_to, 40)
        self.assertEqual(self.dc.zeropad_by, 4.0)
    
    def test_zeropad_by_r_uninit(self):
        """
        without r set, it should give an error
        """
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_by = 4
        self.assertEqual(self.dc.zeropad_to, None)
        self.assertEqual(self.dc.zeropad_by, 1.0)
    
    def test_zeropad_by_nan(self):
        """
        set zeropad_by to numpy.nan
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        DEBUG.verbose("\nError is intentional", True)
        self.dc.zeropad_by = numpy.nan
        self.assertEqual(self.dc.zeropad_to, None)
        self.assertEqual(self.dc.zeropad_by, 1.0)

    def test_zeropad_by_float(self):
        """
        set zeropad_by to float
        """
        self.dc.r = [numpy.ones((10,20)), numpy.ones((10,20))]
        self.dc.zeropad_by = 3.05
        self.assertEqual(self.dc.zeropad_to, 30)
        self.assertEqual(self.dc.zeropad_by, 3.0)



class Test_dataclass_phase(unittest.TestCase):
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

    #####################
    ### PHASE_DEGREES ###
    #####################
    def test_phase_degrees_float(self):
        """
        This is correct
        """
        self.dc.phase_degrees = 31.1
        self.assertEqual(self.dc.phase_degrees, 31.1)
        self.assertAlmostEqual(self.dc.phase_rad, 31.1 * numpy.pi / 180)

    def test_phase_degrees_int(self):
        """
        This is correct
        """
        self.dc.phase_degrees = 42
        self.assertEqual(self.dc.phase_degrees, 42)
        self.assertAlmostEqual(self.dc.phase_rad, 42 * numpy.pi / 180)

    def test_phase_degrees_false(self):
        """
        False should be interpreted as 0
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.phase_degrees = False
        self.assertEqual(self.dc.phase_degrees, 0)
        self.assertEqual(self.dc.phase_rad, 0)

    def test_phase_degrees_true(self):
        """
        True should be interpreted as 1
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.phase_degrees = True
        self.assertEqual(self.dc.phase_degrees, 1)
        self.assertAlmostEqual(self.dc.phase_rad, 1 * numpy.pi / 180)

    def test_phase_degrees_nan(self):
        """
        numpy.nan should give error, no value is set
        """
        DEBUG.verbose("\nError is intentional", True)
        self.dc.phase_degrees = numpy.nan
        self.assertEqual(self.dc.phase_degrees, None)
        self.assertEqual(self.dc.phase_rad, None)

    def test_phase_degrees_none(self):
        """
        None should will reset phase_degrees to init value
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.phase_degrees = None
        self.assertEqual(self.dc.phase_degrees, None)
        self.assertEqual(self.dc.phase_rad, None)

    #################
    ### PHASE_RAD ###
    #################  
    def test_phase_rad_float(self):
        """
        This is correct
        """
        self.dc.phase_rad = 1.1
        self.assertAlmostEqual(self.dc.phase_degrees, 1.1 * 180 / numpy.pi)
        self.assertEqual(self.dc.phase_rad, 1.1)
    
    def test_phase_rad_int(self):
        """
        This is correct
        """
        self.dc.phase_rad = 2
        self.assertAlmostEqual(self.dc.phase_degrees, 2 * 180 / numpy.pi)
        self.assertEqual(self.dc.phase_rad, 2)
    
    def test_phase_rad_false(self):
        """
        False should be interpreted as 0
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.phase_rad = False
        self.assertEqual(self.dc.phase_degrees, 0)
        self.assertEqual(self.dc.phase_rad, 0)
    
    def test_phase_rad_true(self):
        """
        True should be interpreted as 1
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.phase_rad = True
        self.assertAlmostEqual(self.dc.phase_degrees, 1 * 180 / numpy.pi)
        self.assertEqual(self.dc.phase_rad, 1)
    
    def test_phase_rad_nan(self):
        """
        numpy.nan should give error, no value is set
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.phase_rad = numpy.nan
        self.assertEqual(self.dc.phase_degrees, None)
        self.assertEqual(self.dc.phase_rad, None)

    def test_phase_rad_nan(self):
        """
        None sets phase_degrees to init value
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.phase_rad = None
        self.assertEqual(self.dc.phase_degrees, None)
        self.assertEqual(self.dc.phase_rad, None)

    ##################
    ### TIME_STAMP ###
    ##################  
    
    # time_stamp should be have a 'hhmm' format. If only 'hmm' is given (measurement before 12:00), prepend a zero. Give a warning for non-sensical time stamps.
    
    def test_time_stamp_4_int(self):
        """
        This is correct
        """
        self.dc.time_stamp = 1234
        self.assertEqual(self.dc.time_stamp, "1234")

    def test_time_stamp_3_int(self):
        """
        If only 3 digits are given, prepend a zero.
        """
        self.dc.time_stamp = 234
        self.assertEqual(self.dc.time_stamp, "0234")

    def test_time_stamp_4_str(self):
        """
        This is correct
        """
        self.dc.time_stamp = "1234"
        self.assertEqual(self.dc.time_stamp, "1234")

    def test_time_stamp_5_int(self):
        """
        This is correct
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.time_stamp = 12345
        self.assertEqual(self.dc.time_stamp, "12345")

    def test_time_stamp_False(self):
        """
        This is correct
        """
        DEBUG.verbose("\nWarning is intentional", True)
        self.dc.time_stamp = False
        self.assertEqual(self.dc.time_stamp, "False")


if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_dataclass_zeropad)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 1: zeropad getters/setters", True)

    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_dataclass_phase)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 2: phase getters/setters", True)