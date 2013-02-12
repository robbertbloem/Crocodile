import Crocodile.Resources.DataClass

import Crocodile.Pe
import Crocodile.Pe_tw
import Crocodile.Pe_merge

import Crocodile.Resources.Constants
import Crocodile.Resources.Equations
import Crocodile.Resources.IOMethods
import Crocodile.Resources.Mathematics
import Crocodile.Resources.Pe_tw_vb6
import Crocodile.Resources.Plotting

import PythonTools.Debug as DEBUG

def reload_crocodile(flag_verbose = False):
    
    DEBUG.verbose("Reloading Crocodile", flag_verbose)
    
    reload(Crocodile.Resources.DataClass)
    
    reload(Crocodile.Pe)
    reload(Crocodile.Pe_tw)
    reload(Crocodile.Pe_merge)
    
    reload(Crocodile.Resources.Constants)
    reload(Crocodile.Resources.Equations)
    reload(Crocodile.Resources.IOMethods)
    reload(Crocodile.Resources.Mathematics)
    reload(Crocodile.Resources.Pe_tw_vb6)
    reload(Crocodile.Resources.Plotting)