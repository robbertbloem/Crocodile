from __future__ import print_function
from __future__ import division

import inspect

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe_tw as PETW

reload(PETW)

flag_verbose = True

mess = PETW.pe_LV("test", "VC_1", 1538, 300, flag_verbose = flag_verbose)

mess.path = "/Volumes/public_hamm/PML3/data/20121210/VC_1_1538_T300/"

mess.import_data(flag_verbose)

# mess.plot_T(pixel = -1, flag_verbose = flag_verbose)


mess.zeropad_by = 4

mess.absorptive()

mess.zeropad_by = 2

mess.plot(flag_verbose = flag_verbose)

# print(mess)

# plt.figure()
# plt.contour(mess.s)
# plt.show()