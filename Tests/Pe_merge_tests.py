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
            self.mess[i].phase_degrees = 42
            self.mess[i].n_scans = 1
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
        DEBUG.verbose("\nTwo intentional zeropad warnings", True)
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
        for 1 object, only s exists: r is still as initialized: [0]
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

    def test_2p(self):
        """
        mess[0].s = 1 + [1].s: 2
        (1 + 2) / 2 = 1.5
        """
        DEBUG.verbose("\nWarning about zeropad intentional", True)
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0], self.mess[1]], class_min = [], flag_verbose = self.flag_verbose)
        self.assertTrue(numpy.all(mer.s == 1.5))

    def test_2m(self):
        """
        - mess[0].s = 1 - [1].s: 2
        -(1 + 2) / 2 = -1.5
        """
        DEBUG.verbose("\nWarning about zeropad intentional", True)
        mer = PEME.pe_merge("Test", class_plus = [], class_min = [self.mess[0], self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertTrue(numpy.all(mer.s == -1.5))

    def test_2p_1m(self):
        """
        mess[0].s = 1 + [1].s: 2 - [2].s: 3
        (1 + 2) / 2 - 3 = -1.5
        """
        DEBUG.verbose("\nWarning about zeropad intentional", True)
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0], self.mess[1]], class_min = [self.mess[2]], flag_verbose = self.flag_verbose)
        self.assertTrue(numpy.all(mer.s == -1.5))


    def test_phase_degrees_uninit_1(self):
        """
        mess[0].phase_degrees = None
        s=[0] to suppress zeropad warning
        """
        self.mess[0].s = [0]
        DEBUG.verbose("\nIntentional phase warning", True)
        self.mess[0].phase_degrees = None
        DEBUG.verbose("\nIntentional phase warning", True)
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0]], class_min = [self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertEqual(mer.phase_degrees, None)  

    def test_phase_degrees_uninit_2(self):
        """
        mess[1].phase_degrees = None
        s=[0] to suppress zeropad warning
        Tests if the order makes a difference. It shouldn't
        """
        self.mess[0].s = [0]
        DEBUG.verbose("\nIntentional phase warning", True)
        self.mess[1].phase_degrees = None
        DEBUG.verbose("\nError is intentional", True)
        mer = PEME.pe_merge("Test", class_plus = [self.mess[0]], class_min = [self.mess[1]], flag_verbose = self.flag_verbose)
        self.assertEqual(mer.phase_degrees, None)  



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

       
    def test_int_uninit_0(self):      
        index = 0
        PEME.check_value_set_key(self.A, self.B, "var_int", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_int, 1)

    def test_int_uninit_1(self):      
        index = 1
        PEME.check_value_set_key(self.A, self.B, "var_int", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_int, False)

    def test_int_equal_0(self): 
        index = 0   
        self.A.var_int = 1 
        PEME.check_value_set_key(self.A, self.B, "var_int", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_int, 1)

    def test_int_equal_1(self): 
        index = 1  
        self.A.var_int = 1 
        PEME.check_value_set_key(self.A, self.B, "var_int", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_int, 1)

    def test_int_not_equal_0(self):
        index = 0    
        self.A.var_int = 2  
        PEME.check_value_set_key(self.A, self.B, "var_int", index, flag_verbose = self.flag_verbose)    
        self.assertTrue(numpy.isnan(self.A.var_int))

    def test_int_not_equal_1(self):
        index = 1  
        self.A.var_int = 2  
        PEME.check_value_set_key(self.A, self.B, "var_int", index, flag_verbose = self.flag_verbose)    
        self.assertTrue(numpy.isnan(self.A.var_int))

    def test_list_uninit_0(self): 
        index = 0     
        PEME.check_value_set_key(self.A, self.B, "var_list", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_list, [0,0])
        
    def test_list_uninit_1(self): 
        index = 1
        PEME.check_value_set_key(self.A, self.B, "var_list", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_list, False)
    
    def test_list_equal_0(self): 
        index = 0   
        self.A.var_list = [0,0] 
        PEME.check_value_set_key(self.A, self.B, "var_list", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_list, [0,0])
        
    def test_list_equal_1(self): 
        index = 1
        self.A.var_list = [0,0] 
        PEME.check_value_set_key(self.A, self.B, "var_list", index, flag_verbose = self.flag_verbose)    
        self.assertEqual(self.A.var_list, [0,0])

    def test_list_not_equal_0(self):  
        index = 0
        self.A.var_list = [1,1] 
        PEME.check_value_set_key(self.A, self.B, "var_list", index, flag_verbose = self.flag_verbose)  
        self.assertTrue(numpy.isnan(self.A.var_list))

    def test_list_not_equal_1(self):  
        index = 1
        self.A.var_list = [1,1] 
        PEME.check_value_set_key(self.A, self.B, "var_list", index, flag_verbose = self.flag_verbose)  
        self.assertTrue(numpy.isnan(self.A.var_list))
    


class example_class():
    
    def __init__(self):
        
        self.var_int = False
        self.var_list = False
















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