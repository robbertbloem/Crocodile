from __future__ import print_function
from __future__ import division

import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.ObjectArray as OA
import Crocodile.Resources.DataClass as DC



reload(OA)
reload(DC)



class Test_ObjectArray(unittest.TestCase):
    
    def setUp(self):
        
        self.path_and_filename = "/Users/robbert/Developer/Crocodile/temp/test.pickle"
        self.flag_verbose = False

        oa = OA.objectarray("test")
        
        a = OA.testobject("Auto", "a", "power", flag_verbose = self.flag_verbose)
        b = OA.testobject("Boot", "b", "power", flag_verbose = self.flag_verbose)
        c = OA.testobject("Fiets", "c", "human", flag_verbose = self.flag_verbose)
        
        oa.add_object(a, flag_verbose = self.flag_verbose)
        oa.add_object(b, flag_verbose = self.flag_verbose)
        oa.add_object(c, flag_verbose = self.flag_verbose)
        
        oa.save_objectarray(self.path_and_filename, flag_verbose = self.flag_verbose)
        
        self.obj_id_array = oa.obj_id_array

    def test_load_objectarray_1(self):
        """
        Test if object are correctly ordered.
        obj_id_array_in is same as pickle
        """
        oa = OA.objectarray("test_new")
        oa.load_objectarray(self.path_and_filename, obj_id_array_in = self.obj_id_array, flag_verbose = self.flag_verbose) 
        self.assertTrue(["a", "b", "c"] ==  oa.obj_id_array)      

    def test_load_objectarray_2(self):  
        """
        Test if object are correctly ordered.
        obj_id_array_in misses one obj_id compared to pickle. The object array should miss that object.
        """
        obj_id_array = ["a", "b"]
        oa = OA.objectarray("test_new")
        oa.load_objectarray(self.path_and_filename, obj_id_array_in = obj_id_array, flag_verbose = self.flag_verbose) 
        self.assertTrue(["a", "b"] ==  oa.obj_id_array)      

    def test_load_objectarray_3(self):  
        """
        Test if object are correctly ordered.
        obj_id_array_in has one element too much compared to pickle. The object array will only contain the elements of the pickle
        """
        obj_id_array = ["a", "b", "c", "d"]
        oa = OA.objectarray("test_new")
        oa.load_objectarray(self.path_and_filename, obj_id_array_in = obj_id_array, flag_verbose = self.flag_verbose) 
        self.assertTrue(["a", "b", "c"] ==  oa.obj_id_array)   

    def test_load_objectarray_4(self):  
        """
        Test if object are correctly ordered.
        obj_id_array_in has one element too much compared to pickle. The object array will only contain the elements of the pickle
        """
        obj_id_array = ["a", "b", "d"]
        oa = OA.objectarray("test_new")
        oa.load_objectarray(self.path_and_filename, obj_id_array_in = obj_id_array, flag_verbose = self.flag_verbose) 
        self.assertTrue(["a", "b"] ==  oa.obj_id_array)

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(Test_ObjectArray)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    