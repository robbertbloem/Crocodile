# python 3 compatibility
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals
from imp import reload 

import Crocodile.Resources.DataClass

import Crocodile.Pe
import Crocodile.Pe_tw
import Crocodile.Pe_merge

import Crocodile.Resources.Functions
import Crocodile.Resources.Constants
import Crocodile.Resources.Equations
import Crocodile.Resources.IOMethods
import Crocodile.Resources.Mathematics
import Crocodile.Resources.Pe_tw_vb6
import Crocodile.Resources.Plotting

import Crocodile.Plugins.Lineshape
import Crocodile.Plugins.Plot_overlap

import PythonTools.Debug as DEBUG

def reload_crocodile(flag_verbose = False):
    
    DEBUG.verbose("Reloading Crocodile", flag_verbose)
    
    reload(Crocodile.Resources.DataClass)
    
    reload(Crocodile.Pe)
    reload(Crocodile.Pe_tw)
    reload(Crocodile.Pe_merge)
    
    reload(Crocodile.Resources.Functions)
    reload(Crocodile.Resources.Constants)
    reload(Crocodile.Resources.Equations)
    reload(Crocodile.Resources.IOMethods)
    reload(Crocodile.Resources.Mathematics)
    reload(Crocodile.Resources.Pe_tw_vb6)
    reload(Crocodile.Resources.Plotting)
    
    reload(Crocodile.Plugins.Lineshape)
    reload(Crocodile.Plugins.Plot_overlap)