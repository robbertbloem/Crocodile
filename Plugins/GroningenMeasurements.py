from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Mosquito as M
import Crocodile.Resources.Functions as FU
import Crocodile.Resources.Plotting as PL
import Crocodile.Resources.Equations as EQ
import Crocodile.Resources.Mathematics as MATH

import imp

imp.reload(M)

class FT2DIR(M.FT2DIR): 

    def make_overview_plot_Groningen(self, saving = False, **kwargs):
        """
        Make a single figure with all plots. 
        
        INPUT:
        - saving: for saving the size is a lot bigger.
        - kwargs:
            - sp, sm, de, du, sc (lists, ndarray): lists that limit the range of spectra to be plotted. By default all spectra will be plotted.
            - aspect (str, 'equal'): the aspect ratio of the axes. 
            - flip_spectrum (bool, False): by default (False) the x axis is w_3. If True, it will be w_1. Make sure the axes ranges are changed as well. 
        
        OUTPUT:
        - -
        
        DESCRIPTION:
        - - 
        
        CHANGELOG:
        2016014/RB: started function
        20160530/RB: copied it to the Groningen class for adjustments
        
        """
    
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)
    
        n_plots = len(sp) * len(sm) * len(de) * len(du) * len(sc)

        x, y = FU.find_subplots(n_plots, flag_verbose = self.flag_verbose)

        if saving:
            fig = plt.figure(figsize = (30,20), dpi = 300)
        else:
            fig = plt.figure()
        ax = [0] * n_plots
        for ax_i in range(n_plots):
            ax[ax_i] = fig.add_subplot(y, x, ax_i + 1)  
            if "aspect" in kwargs:
                if kwargs["aspect"] != False:
                    ax[ax_i].set_aspect(kwargs["aspect"])
            else:
                ax[ax_i].set_aspect("equal")

        ax_i = 0
        for _sp in sp:
            for _sm in sm:
                for _de in de:
                    for _du in du:
                        for _sc in sc:
                        
                            if self.s_axes[4][:,_sm] == [1]:
                            
                                title = "{name}\nparallel, {dex} fs".format(name = self.objectname, spx = self.s_axes[3][_sp], smx = self.s_axes[4][:,_sm], dex = self.s_axes[5][_de])
                            else:
                                title = "{name}\nperpendicular, {dex} fs".format(name = self._basename, spx = self.s_axes[3][_sp], smx = self.s_axes[4][:,_sm], dex = self.s_axes[5][_de])
                                                        
                            if "flip_spectrum" in kwargs and kwargs["flip_spectrum"]:
                                PL.contourplot(self.s[:, :, 0, _sp, _sm, _de, _du, _sc], self.s_axes[1], self.s_axes[0], x_label = "w1 (cm-1)", y_label = "w3 (cm-1)", ax = ax[ax_i], title = title, **kwargs)
                     
                            else:
                                PL.contourplot(self.s[:, :, 0, _sp, _sm, _de, _du, _sc].T, self.s_axes[0], self.s_axes[1], x_label = "w3 (cm-1)", y_label = "w1 (cm-1)", ax = ax[ax_i], title = title, **kwargs)

                            ax_i += 1
                            
        return fig


if __name__ == "__main__": 
    pass