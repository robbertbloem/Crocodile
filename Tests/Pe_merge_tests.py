from __future__ import print_function
from __future__ import division

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe as PE
import Crocodile.Pe_merge as PEME
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

# add arguments
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
parser.add_argument("-s1", "--skip1", action="store_true", help="Skip testing suite 1: pe_merge")
parser.add_argument("-s2", "--skip2", action="store_true", help="Skip testing suite 2: check_value_set_key")
# parser.add_argument("-s3", "--skip3", action="store_true", help="Skip testing suite 3: print objects")
# parser.add_argument("-s4", "--skip4", action="store_true", help="Skip testing suite 4: add array with objects")

# process
args = parser.parse_args()

# reload
if args.reload:
    reload(PE)
    reload(PEME)



class Test_Pe_merge(unittest.TestCase):
    """

    CHANGELOG:
    20130208/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        
        self.flag_verbose = args.verbose
        
        self.mess = [0] * 3

        for i in range(len(self.mess)):
            
            fake_data = (i+1) * numpy.ones((10,20))
            
            self.mess[i] = PE.pe("Test", flag_verbose = self.flag_verbose)
            
            self.mess[i].r = [fake_data, fake_data]
            self.mess[i].r_axis = [numpy.arange(10), 10, numpy.arange(20)]   

            self.mess[i].undersampling = 0
            # assign after r to enabled calculation of zeropad_by
            # assign before s to prevent a warning that it changed
            self.mess[i].zeropad_to = 0 
            
            self.mess[i].s = fake_data
            self.mess[i].s_axis = [numpy.arange(10), 10, numpy.arange(20)]    
            
        self.mess[0].obj_id = "A"
        self.mess[1].obj_id = "B"
        self.mess[2].obj_id = "C"



    def test_rs_1(self):
        """
        This should work correct
        r and s exist
        """
        DEBUG.verbose("\nZeropad warning is intentional", True)
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0]], class_min = [self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertTrue(numpy.all(mer.s))

      
    def test_s_1(self):
        """
        for 1 object s is still as initialized: [0]
        no zeropad warning because the s is not set and there is no need
        """
        self.mess[0].s = [0]
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0]], class_min = [self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertEqual(mer.s, [0])        

    def test_s_2(self):
        """
        for 1 object, s has a different size
        no zeropad warning because the s is not set and there is no need
        """
        self.mess[0].s = numpy.ones((20,20))
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0]], class_min = [self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertEqual(mer.s, [0])      



    def test_r_1(self):
        """
        for 1 object, only r exists: s is still as initialized: [0]
        no zeropad warning because it is not set
        """
        self.mess[0].r = [0]
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0]], class_min = [self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertEqual(mer.r, [0,0])           
        
    def test_r_2(self):
        """
        for 1 object, r has a different size
        no zeropad warning because it is not set
        """
        self.mess[0].r[0] = numpy.ones((20,20))
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0]], class_min = [self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertEqual(mer.r, [0,0])     







class Test_check_value_set_key(unittest.TestCase):
    """

    CHANGELOG:
    20130208/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):

        self.flag_verbose = args.verbose
        
        self.A = example_class()
        self.B = example_class()
        
        self.B.var_int = 1
        self.B.var_list = [0,0]
        self.B.zeropad_to = 2

       
    def test_int_uninit(self):      
        PEME.check_value_set_key(self.A, self.B, "var_int")    
        self.assertEqual(self.A.var_int, 1)

    def test_int_equal(self):    
        self.A.var_int = 1 
        PEME.check_value_set_key(self.A, self.B, "var_int")    
        self.assertEqual(self.A.var_int, 1)

    def test_int_not_equal(self):    
        self.A.var_int = 2  
        PEME.check_value_set_key(self.A, self.B, "var_int")    
        self.assertTrue(numpy.isnan(self.A.var_int))


    def test_list_uninit(self):      
        PEME.check_value_set_key(self.A, self.B, "var_list")    
        self.assertEqual(self.A.var_list, [0,0])
    
    def test_list_equal(self):    
        self.A.var_list = [0,0] 
        PEME.check_value_set_key(self.A, self.B, "var_list")    
        self.assertEqual(self.A.var_list, [0,0])

    def test_list_not_equal(self):    
        self.A.var_list = [1,1] 
        PEME.check_value_set_key(self.A, self.B, "var_list")  
        self.assertTrue(numpy.isnan(self.A.var_list))


    def test_zeropad_uninit(self):      
        PEME.check_value_set_key(self.A, self.B, "zeropad_to")    
        self.assertEqual(self.A.zeropad_to, 2)
    
    def test_zeropad_equal(self):    
        self.A.zeropad_to = 2
        PEME.check_value_set_key(self.A, self.B, "zeropad_to")    
        self.assertEqual(self.A.zeropad_to, 2)
    
    def test_zeropad_not_equal(self):    
        self.A.zeropad_to = 1
        DEBUG.verbose("\nZeropad warning is intentional", True)
        PEME.check_value_set_key(self.A, self.B, "zeropad_to")  
        self.assertFalse(self.A.zeropad_to)



class example_class():
    
    def __init__(self):
        
        self.var_int = False
        self.var_list = False
        self.zeropad_to = False















if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_Pe_merge)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 1: pe init", True)

    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_check_value_set_key)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 2: check_value_set_key", True)