from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import argparse
import unittest

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.Functions as FU
import PythonTools.Debug as DEBUG

# init argument parser
parser = argparse.ArgumentParser(description='Command line arguments')

# add arguments
parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
parser.add_argument("-r", "--reload", action="store_true", help="Reload modules")
parser.add_argument("-s1", "--skip1", action="store_true", help="Skip testing suite 1: pe_merge")
parser.add_argument("-s2", "--skip2", action="store_true", help="Skip testing suite 2: check_value_set_key")
parser.add_argument("-s3", "--skip3", action="store_true", help="Skip testing suite 3: print objects")
# parser.add_argument("-s4", "--skip4", action="store_true", help="Skip testing suite 4: add array with objects")

# process
args = parser.parse_args()

# reload
if args.reload:
    import Crocodile.Resources.ReloadCrocodile
    Crocodile.Resources.ReloadCrocodile.reload_crocodile(flag_verbose = args.verbose)
    



# find_axes(x_axis, y_axis, x_range, y_range, flag_verbose = False)


class Test_find_axes_x(unittest.TestCase):
    """

    CHANGELOG:
    20130213/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose

        self.axis_1 = numpy.arange(1000, 2000, 51)
        self.axis_2 = numpy.arange(0, 3000, 100)
        
        # all the tests and the results
        # x_range, y_range, results for x_axis = axis_1, x_axis = axis_2
        
        self.param = [
# both complete
[[0,0],         [0,0],          (1000,1969,0,2900),     (0,2900,1000,1969)],  
# one like the other       
[[0,-1],        [0,0],          (0,2900,0,2900),        (1000,1969,1000,1969)],           
[[0,0],         [0,-1],         (1000,1969,1000,1969),  (0,2900,0,2900)],     
[[0,-1],        [0,-1],         (1000,1969,0,2900),     (0,2900,1000,1969)],  
# for a defined set     
[[1300,1500],   [0,0],          (1300,1500,0,2900),     (1300,1500,1000,1969)],   
[[0,0],         [1600,1800],    (1000,1969,1600,1800),  (0,2900,1600,1800)],
[[1300,1500],   [1600,1800],    (1300,1500,1600,1800),  (1300,1500,1600,1800)],  
# weird input
[[0,-2],        [0,0],          (0,-2,0,2900),          (0,-2,1000,1969)], 
        ]  

    def function(self, x_axis, y_axis, param):
        return FU.find_axes(x_axis, y_axis, param[0], param[1], flag_verbose = self.flag_verbose)


    def test_params_x(self):

        for i in range(len(self.param)):
            res = self.function(self.axis_1, self.axis_2, self.param[i])
            self.assertEqual(res, self.param[i][2])

    def test_params_y(self):

        for i in range(len(self.param)):
            res = self.function(self.axis_2, self.axis_1, self.param[i])
            self.assertEqual(res, self.param[i][3])            



# make_contours_2d


class Test_make_contours_2d(unittest.TestCase):
    """

    CHANGELOG:
    20130213/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose
        self.data = numpy.array([[-3,-2],[1,2]])

    # contours = odd
    def test_c3_zm1(self):
        """
        contours = 3, zlimit = -1
        res == [-3,0,3]
        """        
        res = FU.make_contours_2d(self.data, contours = 3, zlimit = -1)
        self.assertTrue(numpy.all(res == [-3,0,3]))

    def test_c3_z0(self):
        """
        contours = 3, zlimit = 0
        res == [-3,-0.5,2]
        """        
        res = FU.make_contours_2d(self.data, contours = 3, zlimit = 0)
        self.assertTrue(numpy.all(res == [-3,-0.5,2]))
        
    def test_c3_z2(self):
        """
        contours = 3, zlimit = 2
        res == [-2,0,2]
        """        
        res = FU.make_contours_2d(self.data, contours = 3, zlimit = 2)
        self.assertTrue(numpy.all(res == [-2,0,2]))

    def test_c3_zlist(self):
        """
        contours = 3, zlimit = 2
        res == [-2,0,2]
        """        
        res = FU.make_contours_2d(self.data, contours = 3, zlimit = [-3,1])
        self.assertTrue(numpy.all(res == [-3,-1,1]))

    # contours = even
    def test_c4_zm1(self):
        """
        contours = 4, zlimit = -1
        res == [-3,-1,1,3]
        """        
        res = FU.make_contours_2d(self.data, contours = 4, zlimit = -1)
        self.assertTrue(numpy.all(res == [-3,-1,1,3]))

    def test_c4_z0(self):
        """
        contours = 4, zlimit = 0
        res == [-3,-1.33333,0.33333,2]
        using numpy.allclose() for rounding error
        """        
        res = FU.make_contours_2d(self.data, contours = 4, zlimit = 0)
        self.assertTrue(numpy.allclose(res, [-3,-1.33333,0.33333,2]))

    def test_c4_z2(self):
        """
        contours = 4, zlimit = 2.4
        res == [-2.4,-0.8,0.8,2.4]
        using numpy.allclose() for rounding error in third element
        """        
        res = FU.make_contours_2d(self.data, contours = 4, zlimit = 2.4)
        self.assertTrue(numpy.allclose(res, [-2.4, -0.8 , 0.8 , 2.4]))

    def test_c4_zlist(self):
        """
        contours = 4, zlimit = [-4,2]
        res == [-4,-2,0,2]
        """        
        res = FU.make_contours_2d(self.data, contours = 4, zlimit = [-4,2])
        self.assertTrue(numpy.all(res == [-4,-2,0,2]))

    # weird input, not checked for both contours
    def test_c4_zm2(self):
        """
        zlimit is negative, 
        zlimit = -2 should give same result as zlimit = 2
        contours = 4, zlimit = -2
        res == [-3,-1,1,3]
        """        
        res = FU.make_contours_2d(self.data, contours = 4, zlimit = -2)
        self.assertTrue(numpy.allclose(res, [-2, -0.66666667, 0.66666667,2]))

    def test_c4_zlist_reverse(self):
        """
        list from high to low
        contours = 4, zlimit = [2,-4]
        res == [2,0,-2,-4]
        """        
        res = FU.make_contours_2d(self.data, contours = 4, zlimit = [2,-4])
        self.assertTrue(numpy.all(res == [2,0,-2,-4]))



class Test_find_axes_indices(unittest.TestCase):
    """

    CHANGELOG:
    20130213/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose
        self.axis = numpy.linspace(20,30,12)

# self.axis
# 0 20.0
# 1 20.9090909091
# 2 21.8181818182
# 3 22.7272727273
# 4 23.6363636364
# 5 24.5454545455
# 6 25.4545454545
# 7 26.3636363636
# 8 27.2727272727
# 9 28.1818181818
# 10 29.0909090909
# 11 30.0
     
    def test_upper_edge_1(self):
        val_min = 22.0
        val_max = 27.2
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (2,8))

    def test_upper_edge_2(self):
        self.axis = numpy.linspace(20,30,11)
        val_min = 22.5
        val_max = 27.0
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (2,8))

    def test_upper_edge_3(self):
        val_min = 22.0
        val_max = 27.3
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (2,9))

    def test_upper_edge_over(self):
        val_min = 22.0
        val_max = 37.3
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (2,11))
        
    def test_lower_edge_1(self):
        val_min = 22.8
        val_max = 27.0
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (3,8))
    
    def test_lower_edge_2(self):
        val_min = 22.7
        val_max = 27.0
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (2,8))

    def test_lower_edge_3(self):
        self.axis = numpy.linspace(20,30,11)
        val_min = 22.0
        val_max = 27.5
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (1,8))
     
    def test_lower_edge_over(self):
        val_min = 1.0
        val_max = 27.2
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (0,8))

    def test_reverse_edges(self):
        """
        The function does no sanity-checking. 
        """
        val_min = 27.0
        val_max = 22.0
        DEBUG.verbose("\nWarning is intentional", True)
        res = FU.find_axes_indices(self.axis, val_min, val_max)
        self.assertEqual(res, (7,3))
   



if __name__ == '__main__':

    if args.skip1 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_axes_x)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 1: importing LabView", True)

    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_make_contours_2d)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 2: importing LabView", True)

    if args.skip3 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_find_axes_indices)
        unittest.TextTestRunner(verbosity=1).run(suite)    
    else:
        DEBUG.verbose("Skipping suite 3: init Pe_tw", True)