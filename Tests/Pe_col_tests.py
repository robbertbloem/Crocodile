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
    "Suite 2: Test of importing",
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
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_of_init)
        unittest.TextTestRunner(verbosity=1).run(suite) 
    else:
        DEBUG.verbose("Skipping: " + suite_list[0], True)
        
    if args.skip2 == False:
        suite = unittest.TestLoader().loadTestsFromTestCase(Test_of_importing)
        unittest.TextTestRunner(verbosity=1).run(suite) 
    else:
        DEBUG.verbose("Skipping: " + suite_list[1], True)



class Test_of_importing(unittest.TestCase):
    """

    CHANGELOG:
    20151104/RB: started the suite

    """
    #############
    ### SETUP ###
    #############
    def setUp(self):
        self.flag_verbose = args.verbose  
        
        self.file_dict_A = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "A",
            "timestamp": "123456",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/A_123456/",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/A_123456/A_123456", 
        }  
        objectname = "test A"   
        self.A = PC.pe_col(objectname, flag_verbose = self.flag_verbose)
        self.A.set_file_info(self.file_dict_A["data_folder"], self.file_dict_A["date"], self.file_dict_A["basename"], self.file_dict_A["timestamp"])

        self.file_dict_B = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "B",
            "timestamp": "234506",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/B_234506/",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/B_234506/B_234506", 
        }       
        objectname = "test B"      
        self.B = PC.pe_col(objectname, flag_verbose = self.flag_verbose)
        self.B.set_file_info(self.file_dict_B["data_folder"], self.file_dict_B["date"], self.file_dict_B["basename"], self.file_dict_B["timestamp"])

        self.file_dict_C = {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "azide",
            "timestamp": "174608",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/azide_174608",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/azide_174608/azide_174608", 
        }       
        objectname = "test C"      
        self.C = PC.pe_col(objectname, flag_verbose = self.flag_verbose)
        self.C.set_file_info(self.file_dict_C["data_folder"], self.file_dict_C["date"], self.file_dict_C["basename"], self.file_dict_C["timestamp"])

        self.file_dict_D= {
            "data_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources",
            "date": "20010101",
            "basename": "azide_intensity",
            "timestamp": "175321",
            "extension": ".csv",
            "base_folder": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/azide_intensity_175321",
            "base_filename": "/Users/robbert/Developer/Crocodile/Tests/Test_resources/20010101/azide_intensity_175321/azide_intensity_175321", 
        }       
        objectname = "test D"      
        self.D = PC.pe_col(objectname, flag_verbose = self.flag_verbose)
        self.D.set_file_info(self.file_dict_D["data_folder"], self.file_dict_D["date"], self.file_dict_D["basename"], self.file_dict_D["timestamp"])

#         print(self.A)


    def test_file_dict_not_set(self):
        """
        When the file_dict is not set, importing should return False
        """
        self.A._file_dict = {"data_folder": "", "date": "", "basename": "", "timestamp": "", "base_folder": "", "base_filename": "", "extension": ""}
        res = self.A.import_data()
        self.assertFalse(res)


    def test_file_dict_with_incorrect_folder(self):
        """
        When the file_dict is incorrect (it points to a non-existing folder) it should return False.
        """       
        self.A.set_file_info("Not a folder", self.file_dict_A["date"], self.file_dict_A["basename"], self.file_dict_A["timestamp"])
        res = self.A.import_data()
        self.assertFalse(res)

       
    def test_importing_but_file_format_file_does_not_exist(self):
        """
        When the file_dict is not found, importing should return True
        """        
        res = self.B.import_data()    
        self.assertFalse(res)


#     def test_importing_A(self):
#         """
#         When the file_dict is set, importing should return True
#         """        
#         res = self.A.import_data()         
#         self.assertTrue(res)

    def test_importing_C(self):
        """
        When the file_dict is set, importing should return True
        """        
        res = self.C.import_data(import_temp_scans = False) 
        self.C.b_to_r()
        self.C.calculate_phase()
          
#         print(self.C)    
        print(numpy.shape(self.C.r))
        self.assertTrue(res)


    def test_importing_D(self):
        """
        When the file_dict is set, importing should return True
        """        
        res = self.D.import_data(import_temp_scans = False) 
        self.D.b_to_r()
        self.D.calculate_phase()
          
#         print(self.D)    
        print(numpy.shape(self.D.r))
        self.assertTrue(res)

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
        
        self.assertEqual(x.file_dict["base_filename"], "/User/Robbert/Bla/20010101/foo_123456/foo_123456")




        
        
 
if __name__ == '__main__':

    execute(args)
