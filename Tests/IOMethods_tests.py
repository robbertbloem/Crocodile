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
    # "Suite 2: Phase_degrees setters/getters",
]

# add arguments
parser.add_argument("-v", "--verbose", action = "store_true", help = "Increase output verbosity")
parser.add_argument("-r", "--reload", action = "store_true", help = "Reload modules")
parser.add_argument("-s1", "--skip1", action = "store_true", help = suite_list[0])
# parser.add_argument("-s2", "--skip2", action = "store_true", help = suite_list[1])

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
        DEBUG.verbose("Skipping :" + suite_list[0], True)




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
