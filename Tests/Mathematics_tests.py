from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from imp import reload

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.Mathematics as MATH
import Crocodile.Resources.Equations as EQ
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

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
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_zeropad_set_get)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping :" + suite_list[0], True)




class Test_zeropad_set_get(unittest.TestCase):
    """

    CHANGELOG:
    20130208/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose

    def test_fit(self):
        x = numpy.arange(10)
        y = numpy.sin(x)
        A = [0,1,2,3]
        A = MATH.fit(x, y, EQ.rb_cos, A)
        print(A)



     
        
if __name__ == '__main__':

    execute(args)
