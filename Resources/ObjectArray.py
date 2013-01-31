from __future__ import print_function
from __future__ import division

import inspect
import shelve
import os

import Crocodile
import Crocodile.Resources.ClassTools as CT

reload(CT)


class objectarray(CT.ClassTools):
    
    def __init__(self, name = "name", obj_id = "objectarray", flag_verbose = False):
        
        self.verbose("Created object array", flag_verbose)
        
        self.obj_id = obj_id    # identifier
        
        self.name = name        # name - optional
        
        self.obj_array = []     # array with objects
        self.obj_id_array = []  # array with object-identifiers, to prevent duplicates. The array is not necesarily in the correct order.



    def add_object(self, obj, flag_verbose = False):
        """
        Add an object to the array
        
        """
    
        self.verbose("add_object", flag_verbose)
        
        # test if object-id is unique
        if obj.obj_id in self.obj_id_array:
            self.printError("obj_id already exists, will not add new object.", inspect.stack())
            return False
        else:
            self.verbose("  New object is appended.", flag_verbose)
            self.obj_array.append(obj)
            self.obj_id_array.append(obj.obj_id)
            return True
        

    def save_objectarray(self, path_and_filename, flag_overwrite = False, flag_verbose = False):
        """
        Save the array
        """

        self.verbose("Save object array", flag_verbose)

        # check filename
        if path_and_filename[-7:] != ".pickle":
            path_and_filename += ".pickle"
            self.verbose("  Added .pickle to path_and_filename", flag_verbose)
        
        # check if you want to overwrite it
        if flag_overwrite:
            self.printWarning("make_db: overwrite flag is True", inspect.stack())
            flag_overwrite = "n"
        else:
            flag_overwrite = "c"  

        # save it
        db = shelve.open(path_and_filename, flag = flag_overwrite)
        for object in self.obj_array:
            db[object.obj_id] = object
        db.close()
        os.system("chmod 777 " + path_and_filename)
        
        return True

        
    def import_db(self, path_and_filename, flag_verbose = False):
        """
        Imports a database. The function checks for the existence of the database. It returns "False" if the file doesn't exist. Otherwise, it will return an array with class instances.
        """
        
        if path_and_filename[-7:] != ".pickle":
            path_and_filename += ".pickle"
            self.verbose("  Added .pickle to path_and_filename", flag_verbose)
     
        if os.path.isfile(path_and_filename) == True:
            db=shelve.open(path_and_filename)
            obj_array = []
            for key in db:
                if flag_verbose:
                    self.verbose(key, flag_verbose)
                self.obj_array.append(db[key])
                self.obj_id_array.append(db[key].obj_id)
            db.close() 
            return True
        else:     
            self.printError("The file doesn't exist!", inspect.stack())
            return False
      

    def print_objects(self, flag_verbose = False):
        """
        Print the objects of the array. This will print all the details of the objects.
        """
        self.verbose("Print objects", flag_verbose)
        for i in self.obj_array:
            print(i)

    def print_object_ids(self, flag_verbose = False):
        """
        Print the objects of the array. This will only print the object ids.
        """
        self.verbose("Print object ids", flag_verbose)
        for i in self.obj_id_array:
            print(i)


class testobject(CT.ClassTools):
    
    def __init__(self, name, obj_id, flag_verbose = False):
        
        self.verbose("Create test object", flag_verbose)     
        self.name = name
        self.obj_id = obj_id



if __name__ == "__main__": 

    pass
  
    # flag_verbose = True
    # 
    # a = objectarray(flag_verbose)
    # 
    # b = testobject("Auto", "a")
    # c = testobject("Fiets", "b")
    # 
    # a.add_object(b, flag_verbose)
    # a.add_object(c, flag_verbose)
    # # a.add_object(c, flag_verbose = True)
    # 
    # path = "/Users/robbert/Developer/Crocodile/temp/test.pickle"
    # a.save_objectarray(path)
    # 
    # a = 0
    # 
    # a = objectarray(flag_verbose)
    # 
    # a.import_db(path, flag_verbose)
    # 
    # a.print_objects(flag_verbose)
    # a.print_object_ids(flag_verbose)
    # 
    # print(a)
    
    
    
    
    
    
    # print(a)
    # print(b)
    # print(c)