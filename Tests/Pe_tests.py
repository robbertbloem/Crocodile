from __future__ import print_function
from __future__ import division

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile
import Crocodile.Pe as PE

reload(PE)

flag_verbose = True

a = PE.pe("Fiets", flag_verbose)

x = numpy.arange(10)
y = numpy.arange(10)

X, Y = numpy.meshgrid(x,y)

Z = numpy.cos(X+Y)

a.s = Z
a.s_axis = [y, 0, x]

a.plot(pixel = 2, flag_verbose = flag_verbose)

fig = plt.figure()
ax = [0,0]
ax[0] = fig.add_subplot(211)
ax[1] = fig.add_subplot(212)

a.plot(ax=ax[0], flag_verbose = flag_verbose)
a.plot(ax=ax[1], flipxy=True, flag_verbose = flag_verbose)

print(a)