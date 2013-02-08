from __future__ import print_function
from __future__ import division



import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.ObjectArray as OA
import Crocodile.Pe_tw as PETW

reload(OA)
reload(PETW)

flag_verbose = False

oa = OA.objectarray("20130206")
        
# a = OA.testobject("Auto", "a", "power", flag_verbose = self.flag_verbose)
# self.oa.add_object(a, flag_verbose = flag_verbose)

mess = PETW.pe_LV("aha_1", "aha", 1748, 300, flag_verbose = flag_verbose)
oa.add_object(mess, flag_verbose = flag_verbose)
oa.obj[0].path = "/Volumes/public_hamm/PML3/data/20130206/aha_1748_T300/"
oa.obj[0].import_data(flag_verbose)
oa.obj[0].zeropad_by = 4
oa.obj[0].absorptive()
# oa.obj[0].plot(flag_verbose = flag_verbose)

mess = PETW.pe_LV("aha_2", "aha", 1750, 300, flag_verbose = flag_verbose)
oa.add_object(mess, flag_verbose = flag_verbose)
oa.obj[1].path = "/Volumes/public_hamm/PML3/data/20130206/aha_1750_T300/"
oa.obj[1].import_data(flag_verbose)
oa.obj[1].zeropad_by = 4
oa.obj[1].absorptive()
# oa.obj[1].plot(flag_verbose = flag_verbose)

# oa.print_objects()
oa.print_objects(index=-1)
oa.print_objects(index=2)

oa.print_object_ids(flag_verbose = flag_verbose)
