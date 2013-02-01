from __future__ import print_function
from __future__ import division

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Resources.ObjectArray as OA

reload(OA)

flag_verbose = True

a = OA.objectarray(flag_verbose)

b = OA.testobject("Auto", "a")
c = OA.testobject("Fiets", "b")
d = OA.testobject("Boot", "c")

a.add_object(b, flag_verbose)
a.add_object(c, flag_verbose)
a.add_object(d, flag_verbose)


a.print_object_ids(flag_verbose)


path = "/Users/robbert/Developer/Crocodile/temp/test.pickle"
a.save_objectarray(path)

a = 0

a = OA.objectarray(flag_verbose)

a.import_db(path, flag_verbose)

a.print_objects(flag_verbose)
a.print_object_ids(flag_verbose)

print(a)


