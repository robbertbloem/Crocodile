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


        
    def export_2d_data_for_Maxim(self, paf, s = True, **kwargs):
        
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)
    
#         self.check_results_folder()
    
#         paf = self.file_dict["result_folder"] + self._file_dict["basename"] + "_" + self._file_dict["timestamp"]
    
        ax_i = 0
        for _sp in sp:
            for _sm in sm:
                for _de in de:
                    for _du in du:
                        for _sc in sc:
                            
                            if s: 
                                
                                x_axis = self.s_axes[0]
                                y_axis = self.s_axes[1]
                                data = self.s[:,:, 0, _sp, _sm, _de, _du, _sc]
                                
                                if _sm == 0:
                                    sm_name = "perp"
                                else:
                                    sm_name = "para"
                                
                                path_and_filename =  "{paf}_s_{sm_name}_{de}fs{ext}".format(paf = paf , sm_name = sm_name, de = int(self.s_axes[5][_de]), ext = self.file_dict["extension"])  
                            
                                self.export_2d_data_for_Maxim_helper(x_axis, y_axis, data, path_and_filename, **kwargs)     


    def export_2d_data_for_Maxim_helper(self, x_axis, y_axis, data, path_and_filename, **kwargs):
        """
        
    
        INPUT:
        - data:
        - x_range, y_range (optional)
    
        OUTPUT:

    
        CHANGELOG:
        20160502-RB: started function

        """
        
        if "x_range" in kwargs:
            x_range = kwargs["x_range"]
        else:
            x_range = [0,0]

        if "y_range" in kwargs:
            y_range = kwargs["y_range"]
        else:
            y_range = [0,-1]

        # determine the range to be plotted
        x_min, x_max, y_min, y_max = FU.find_axes(x_axis, y_axis, x_range, y_range, self.flag_verbose)
    
        # find the area to be plotted
        x_min_i, x_max_i = FU.find_axes_indices(x_axis, x_min, x_max)
        y_min_i, y_max_i = FU.find_axes_indices(y_axis, y_min, y_max)
          
        # truncate the data
        data, x_axis, y_axis = FU.truncate_data(data, y_axis, x_axis, y_min_i, y_max_i, x_min_i, x_max_i) 
        
        n_x = len(x_axis)
        n_y = len(y_axis)
            
        with open(path_and_filename, 'w') as f:
        
            for x_i in range(n_x):
                for y_i in range(n_y):
                    f.write("{y:4.3f}, {x:4.3f}, {z:.8f}\n".format(y = y_axis[y_i], x = x_axis[x_i], z = data[y_i, x_i]))
            
#                 f.write("\n")
            
        f.close()


class show_spectrum(M.show_spectrum): 

    def export_1d_data_for_Maxim(self, paf, s = True, **kwargs):
        
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)
    
        ax_i = 0
        for _ds in ds:

            x_axis = self.s_axes[0]
            probe = self.r[:,0,0,0, 0,0,0,0]
            ref = self.r[:,0,1,0, 0,0,0,0]
            
            path_and_filename =  "{paf}_r_probe{ext}".format(paf = paf, ext = self.file_dict["extension"])  
            self.export_1d_data_for_Maxim_helper(x_axis, probe, path_and_filename, **kwargs)   
            
            path_and_filename =  "{paf}_r_ref{ext}".format(paf = paf, ext = self.file_dict["extension"]) 
            self.export_1d_data_for_Maxim_helper(x_axis, ref, path_and_filename, **kwargs)  


    def export_1d_data_for_Maxim_helper(self, x_axis, data, path_and_filename, **kwargs):
        """
        
    
        INPUT:
        - data:
        - x_range, y_range (optional)
    
        OUTPUT:

    
        CHANGELOG:
        20160502-RB: started function

        """
    
        with open(path_and_filename, 'w') as f:
        
            for x_i in range(len(x_axis)):
                f.write("{x:4.3f}, {y:1.8f}\n".format(x = x_axis[x_i], y = data[x_i]))
            
        f.close()


class pump_probe(M.pump_probe):

    def make_contour_plot(self, **kwargs):
    
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)

            
        if "y_range" not in kwargs:
            temp = numpy.where(self.s_axes[5] > 0)[0]
            kwargs["y_range"] = [self.s_axes[5][temp[0]], self.s_axes[5][temp[-1]]]
        
        for _sp in sp:
            for _sm in sm:
                for _du in du:
                    for _sc in sc:
                        
                        if "ax" not in kwargs:
                            fig = plt.figure()
                            ax = fig.add_subplot(111)  
                        else:
                            ax = kwargs["ax"]
                            del(kwargs["ax"])
                        

                        if _sm == 0:
                            smx = "perpendicular"
                        else:
                            smx = "parallel"
                        
                        
                        kwargs["title"] = "{name}\n{smx}".format(name = self.objectname, smx = smx)

                        PL.contourplot(self.s[:, 0, 0, _sp, _sm, :, _du, _sc].T, self.s_axes[0], self.s_axes[5], x_label = "w3 (cm-1)", y_label = self.s_units[5], ax = ax, **kwargs)
                        
                        if "ylog" in kwargs and kwargs["ylog"]:                    
                            ax.set_yscale('log')


    def linear_plot(self, **kwargs):

        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)

        n_de = self.s_n[5]
        colors = [0] * n_de
        
        for _de in range(n_de): 
            temp = _de / n_de
#             self.s_axes[5][_de] / self.s_axes[5][-2]
            colors[_de] = plt.get_cmap("jet")(temp)
# nipy_spectral
#         http://matplotlib.org/examples/color/named_colors.html
#         colors = [
#             "firebrick",
#             "r",
#             "salmon",
#             "navajowhite", #"tisque",
#             "tan",
#             "gold",
#             "darkkhaki",
#             "olivedrab",
#             "chartreuse",
#             "green",
#             "cadetblue",
#             "skyblue",
#             "dodgerblue",
#             "slategrey",
#         ]
    
#         print(colors)
            
        if "y_range" not in kwargs:
            temp = numpy.where(self.s_axes[5] > 0)[0]
            kwargs["y_range"] = [self.s_axes[5][temp[0]], self.s_axes[5][temp[-1]]]
        
        for _sp in sp:
            for _sm in sm:
                for _du in du:
                    for _sc in sc:
                    
                        if "ax" not in kwargs:
                            fig = plt.figure()
                            ax = fig.add_subplot(111)  
                        else:
                            ax = kwargs["ax"]
                            del(kwargs["ax"])
                    

                        if _sm == 0:
                            smx = "perpendicular"
                        else:
                            smx = "parallel"
                    
                    
                        kwargs["title"] = "{name}\n{smx}".format(name = self.objectname, smx = smx)
                        for _de in de:
                            ax.plot(self.s_axes[0], self.s[:,0,0,_sp, _sm,_de,_du,_sc], color = colors[_de], label = self.s_axes[5][_de])
                            
                        ax.set_xlim(self.s_axes[0][0], self.s_axes[0][-1])
                        ax.set_title(kwargs["title"])
                        
                        ax.set_ylabel("mOD")
                        ax.set_xlabel("w (cm-1)")

        ax.legend(fontsize=10, loc="upper left")

    def export_pp_data_for_Maxim(self, paf, s = True, **kwargs):
        
        pi, bish, sp, ds, sm, de, du, sc = self.multiplot_ranges(**kwargs)
    
        ax_i = 0
        for _sp in sp:
            for _sm in sm:
#                 for _de in de:
                for _du in du:
                    for _sc in sc:
                        
                        if s: 
                            
                            x_axis = self.s_axes[0]
                            y_axis = self.s_axes[5]
                            data = self.s[:,0, 0, _sp, _sm, :, _du, _sc]
                            
                            if _sm == 0:
                                sm_name = "perp"
                            else:
                                sm_name = "para"
                            
                            path_and_filename =  "{paf}_s_{sm_name}{ext}".format(paf = paf , sm_name = sm_name, ext = self.file_dict["extension"])  
                        
                            self.export_pp_data_for_Maxim_helper(x_axis, y_axis, data, path_and_filename, **kwargs)     


    def export_pp_data_for_Maxim_helper(self, x_axis, y_axis, data, path_and_filename, **kwargs):
        """
        
    
        INPUT:
        - data:
        - x_range, y_range (optional)
    
        OUTPUT:

    
        CHANGELOG:
        20160502-RB: started function

        """
        
        n_x = len(x_axis)
        n_y = len(y_axis)
            
        with open(path_and_filename, 'w') as f:
        
            for x_i in range(n_x):
                for y_i in range(n_y):
                    f.write("{y}, {x:4.3f}, {z:1.8f}\n".format(y = y_axis[y_i], x = x_axis[x_i], z = data[x_i, y_i]))
            
#                 f.write("\n")
            
        f.close()   

    
        
    


def merge_linear_spectra(m1, m2):
    """
    
    """
    
    
    # copy lengths
    r_n = m1.r_n[:]
    # sm to 2
    r_n[4] = r_n[4] * 2
    # new r
    r = numpy.zeros(r_n)
    r_noise = numpy.zeros(r_n)
    
    # copy axes
    r_axes = m1.r_axes[:]
    # set sm to [0,1]
    r_axes[4] = numpy.arange(2)
    
    r[:,:,:,:, 0,:,:,:] = m1.r[:,:,:,:, 0,:,:,:]
    r[:,:,:,:, 1,:,:,:] = m2.r[:,:,:,:, 0,:,:,:]

    r_noise[:,:,:,:, 0,:,:,:] = m1.r_noise[:,:,:,:, 0,:,:,:]
    r_noise[:,:,:,:, 1,:,:,:] = m2.r_noise[:,:,:,:, 0,:,:,:]

    m1.r_n = r_n
    m1.r = r
    m1.r_noise = r_noise
    m1.r_axes = r_axes

    # copy lengths
    s_n = m1.s_n[:]
    # sm to 2
    s_n[4] = 2
    # new r
    s = numpy.zeros(s_n)
    s_noise = numpy.zeros(s_n)
    
    # copy axes
    s_axes = m1.s_axes[:]
    # set sm to [0,1]
    s_axes[4] = numpy.arange(2)
    
    s[:,:,:,:, 0,:,:,:] = m1.s[:,:,:,:, 0,:,:,:]
    s[:,:,:,:, 1,:,:,:] = m2.s[:,:,:,:, 0,:,:,:]

    s_noise[:,:,:,:, 0,:,:,:] = m1.s_noise[:,:,:,:, 0,:,:,:]
    s_noise[:,:,:,:, 1,:,:,:] = m2.s_noise[:,:,:,:, 0,:,:,:]

    m1.s_n = s_n
    m1.s = s
    m1.s_noise = s_noise
    m1.s_axes = s_axes
    
    return m1


if __name__ == "__main__": 
    pass