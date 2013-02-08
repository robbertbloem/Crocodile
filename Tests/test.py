from __future__ import print_function
from __future__ import division



import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.ObjectArray as OA
import Crocodile.Pe_tw as PETW
import Crocodile.Pe_merge as PEME

reload(OA)
reload(PETW)
reload(PEME)

flag_verbose = False

oa = OA.objectarray("20130206")
        
mess = PETW.pe_LV("aha_1", "aha", 1748, 300, flag_verbose = flag_verbose)
oa.add_object(mess, flag_verbose = flag_verbose)
oa.obj[0].path = "/Volumes/public_hamm/PML3/data/20130206/aha_1748_T300/"
oa.obj[0].import_data(flag_verbose)
oa.obj[0].zeropad_by = 2
oa.obj[0].absorptive()
oa.obj[0].comment = "aha_1 comment aha_1"
# oa.obj[0].plot(flag_verbose = flag_verbose)

mess = PETW.pe_LV("aha_2", "aha", 1750, 300, flag_verbose = flag_verbose)
oa.add_object(mess, flag_verbose = flag_verbose)
oa.obj[1].path = "/Volumes/public_hamm/PML3/data/20130206/aha_1750_T300/"
oa.obj[1].import_data(flag_verbose)
oa.obj[1].zeropad_by = 2
oa.obj[1].absorptive()
# oa.obj[1].r = [0]
# oa.obj[1].phase_degrees = 50
oa.obj[1].comment = "aha_2 comment aha_2"
# oa.obj[1].plot(flag_verbose = flag_verbose)

mess = PETW.pe_LV("buf_1a", "buf", 1720, 300, flag_verbose = flag_verbose)
oa.add_object(mess, flag_verbose = flag_verbose)
oa.obj[2].path = "/Volumes/public_hamm/PML3/data/20130206/buf_1720_T300/"
oa.obj[2].import_data(flag_verbose)
oa.obj[2].zeropad_by = 2
oa.obj[2].absorptive()
oa.obj[2].comment = "buf_1 comment buf_1"
# oa.obj[2].plot(flag_verbose = flag_verbose)

# print(type(oa.obj[0].phase_degrees))

# oa.print_objects(flag_verbose = flag_verbose)
# 
# print(oa.object_with_sub_type(sub_type = "aha"))
# 
# flag_verbose = True
# 
mess = PEME.pe_merge("fiets", class_plus = [oa.obj[0], oa.obj[1]], class_min = [oa.obj[2]], flag_verbose = flag_verbose)
print(mess)

# mess.plot()