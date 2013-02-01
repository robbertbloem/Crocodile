from __future__ import print_function
from __future__ import division

import inspect
import re

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import Crocodile.Pe_tw as PETW
import Crocodile.Resources.IOMethods as IOM
import Crocodile.Resources.Constants as CONST





class pe_VB6(PETW.pe_tw):
    """


    VB6 data comes in three types:
    old_style: before 2009
    new_style/1.1: before fast scanning, use import_data
    1.2: raw data and fringes
    1.3: raw data, fringes and the counter values
    1.4: binned data
    1.2-1.4: use add_data

    """

    def __init__(self, objectname, root_filename, time_stamp, population_time, flag_verbose = False):

        self.verbose("New Pe.pe.pe_VB6 class", flag_verbose)

        PETW.pe_tw.__init__(self, objectname = objectname, root_filename = root_filename, time_stamp = time_stamp, population_time = population_time, undersampling = 0, flag_verbose = flag_verbose)  

        self.n_fringes = 0 
        self.extra_fringes = 20  
        self.n_channels = 37

        self.imported_scans = []





    def import_data(self, scans = [0], noise = True, meta = True):
        """
        croc.pe.import_data

        Imports data from experimentally measured spectra.

        INPUT
        - scans: 
            - [0]: import file where everything is already averaged
            - [n]: import scan n
            - [a, b, c]: import scans a, b, c and average them
        - noise: will also load the noise files (only for new-style)
        - meta: will also load the meta-file (only for new-style)

        CHANGELOG:
        20110910 RB: started
        20130201/RB: function has been transfered from croc, but is not tested yet.

        """

        # determine old-style or new-style
        if self.time_stamp == 0000:
            self.printError("Old style not implemented yet", inspect.stack())

        # new-style
        else:  
            # construct the file name  
            filebase = self.path + self.base_filename + "_" + str(self.time_stamp) + "_T" + str(self.r_axis[1])

             # import a single scan
            if len(scans) == 1:

                # import the averaged data
                if scans[0] == 0:
                    file_R = filebase + ".dat"
                    file_NR = filebase + "_NR.dat"

                    if noise:
                        file_R_noise = filebase + "_R_noise.dat"
                        file_NR_noise = filebase + "_NR_noise.dat"                    

                # import a single scan
                else:
                    file_R = filebase + "_R_" + str(scans[0]) + ".dat"
                    file_NR = filebase + "_NR_" + str(scans[0]) + ".dat"

                    if noise:
                        file_R_noise = filebase + "_R_" + str(scans[0]) + "_noise.dat"
                        file_NR_noise = filebase + "_NR_" + str(scans[0]) + "noise_.dat"                     

                # load the actual data
                try:
                    temp = numpy.loadtxt(file_R)
                except IOError:
                    D.printError(("Unable to load file: " + file_R), inspect.stack())
                    raise
                    return 0   
                self.r_axis[0] = temp[1:,0]
                self.r_axis[2] = temp[0,1:self.n_pixels+1]
                self.r[0] = temp[1:,1:self.n_pixels+1]

                try:
                    temp = numpy.loadtxt(file_NR)
                except IOError:
                    D.printError(("Unable to load file: " + file_NR), inspect.stack())
                    raise
                    return 0
                self.r[1] = temp[1:,1:self.n_pixels+1]               

                if noise:
                    try:
                        temp = numpy.loadtxt(file_R_noise)
                    except IOError:
                        self.printError(("Unable to load file: " + file_R_noise), inspect.stack())
                        raise
                        return 0
                    self.r_noise[0] = temp[1:,1:self.n_pixels+1]                   

                    try:
                        temp = numpy.loadtxt(file_NR_noise)
                    except IOError:
                        self.printError(("Unable to load file: " + file_NR_noise), inspect.stack())
                        raise
                        return 0
                    self.r_noise[1] = temp[1:,1:self.n_pixels+1]     

                # fill in some details 
                self.r_units = ["fs", "fs", "cm-1"]
                self.r_resolution = [(self.r_axis[0][1] - self.r_axis[0][0]), 0, (self.r_axis[2][1] - self.r_axis[2][0])]   
                self.n_steps = len(self.r_axis[0])         


            # import multiple scans and average them
            else:
                self.printError("The ability to import specific scans is not yet implemented.", inspect.stack())
                return 0

            # import the meta data
            if meta: 
                path_and_filename = self.path + self.base_filename + "_meta.txt"
                out = self.import_meta_VB6(path_and_filename)
                self.date = out[0]
                self.data_type_version = out[1]
                self.n_shots = out[2] 
                self.n_fringes = out[3]
                self.phase_degrees = out[4]
                self.comment = out[5]






    def add_data(self, scan, flag_construct_r = False, flag_verbose = False):

        # import the reference and some other stuff
        if len(self.imported_scans) == 0:

            self.verbose("  importing reference and meta", flag_verbose)

            path_and_filename = self.path + self.base_filename + "_ref.dat"
            out = self.import_reference_VB6(path_and_filename, self.n_pixels)
            self.r_axis[2] = out[0] 
            self.reference = out[1]

            path_and_filename = self.path + self.base_filename + "_meta.txt"
            out = self.import_meta_VB6(path_and_filename)
            self.date = out[0]
            self.data_type_version = out[1]
            self.n_shots = out[2] 
            self.n_fringes = out[3]
            self.phase_degrees = out[4]
            self.comment = out[5]
            # self.temp_scans = out[6]

        # see if we already imported this file
        if scan in self.imported_scans:
            self.printWarning("Scan is already imported, but will be imported anyway because flag_import_override is set to True.", inspect.stack()) 
            return False    
        else:
            self.verbose("  scan has not been imported before ", flag_verbose)

        filename = self.make_filenames(scan)        

        # for the 4 files
        for k in range(4):

            flag_skip_binning = False

            # to distuinguish between the two diagrams
            if k < 2:
                diagram = 0
            else:
                diagram = 1

            # import the data
            try:
                self.verbose("  importing data, data_type_version: " + self.data_type_version, flag_verbose)

                if self.data_type_version == "1.2":
                    [m, fringes] = IOM.import_raw_data(filename[k], n_shots = self.n_shots, n_channels = self.n_channels, flag_counter = False)
                elif self.data_type_version == "1.3":
                    [m, c, fringes] = IOM.import_raw_data(filename[k], n_shots = self.n_shots, n_channels = self.n_channels, flag_counter = True)
                elif self.data_type_version == "1.4":
                    result = self.import_binned_data(filename[k], diagram)
                    flag_skip_binning = True

                    if result == False:
                        try:
                            [m, c, fringes] = self.import_raw_data(filename[k], flag_counter = True)
                            flag_skip_binning = False
                        except IOError:
                            return False

            except IOError:
                return False

            # if the number of fringes can not be set using the meta file, find it here 
            if self.n_fringes == 0:
                self.n_fringes = int(numpy.abs(fringes[1] - fringes[0]))

            # counter
            if self.data_type_version == "1.2":
                # need to reconstruct it
                m_axis, counter, correct_count = self.reconstruct_counter(m, self.x_channel, self.y_channel, fringes[0], fringes[1], flag_plot = False)

            elif self.data_type_version == "1.3":   
                # only need to check the value
                m_axis = c + fringes[0]

                if m_axis[-1] == fringes[1]:
                    correct_count = True
                else:
                    correct_count = False

            elif self.data_type_version == "1.4": 
                # the case where the binning failed and VB6 reverts to data version 1.3
                if flag_skip_binning == False:
                    m_axis = c + fringes[0]

                    if m_axis[-1] == fringes[1]:
                        correct_count = True
                    else:
                        correct_count = False   
                else:
                    correct_count = True

            else:
                self.printError("Unknown data type", inspect.stack())
                correct_count = False

            # check for consistency
            if correct_count == False:
                print("Scan: " + str(scan) + ", File: " + str(k) + ": Miscount!")

                self.incorrect_count[k] += 1

            # if it is consistent, continue to bin the data
            elif flag_skip_binning == False:
                print("Scan: " + str(scan) + ", File: " + str(k) + ": Count is correct!")

                # make b the correct size, if it isn't already
                if numpy.shape(self.b_axis)[-1] == 2:
                    self.make_arrays()

                # bin the data
                self.bin_data(m, m_axis, diagram)

            else:
                print("Scan: " + str(scan) + ", File: " + str(k) + ": Scan imported")


                # all the data is now written into self.b* 

        # construct the actual measurement
        # this should catch the situation that there are 4 miscounts in the first scan
        if type(self.b_count) == numpy.ndarray:
            if flag_construct_r:
                self.construct_r()
            self.imported_scans.append(scan)
            self.n_scans = len(self.imported_scans)            
            return True
        else:
            return False




    def bin_data(self, m, m_axis, diagram):

        i_range = range(self.n_shots)

        for i in i_range:            
            # find the fringe
            j = (-1)**diagram * int(m_axis[i]) + self.extra_fringes - (-1)**diagram * 4000

            # add it to the bin, depending on pem-state and diagram
            # and add 1 to counter 
            if m[self.chopper_channel, i] < 2.5:         
                self.b[2 * diagram][j, :] += m[:,i] 
                self.b_count[2 * diagram, j] += 1
            else:
                self.b[2 * diagram + 1][j, :] += m[:,i] 
                self.b_count[2 * diagram + 1, j] += 1






    def make_filenames(self, scan):
        """
        croc.Pe.pefs.make_filenames()

        Makes the filenames. Apart from the path, all the variables should have been entered in the init. 

        INPUT:
        - scan (int): the number of the scan

        IMPLICIT REQUIREMENTS:
        In the data structure, the following variables should be set:
        - self.path
        - self.base_filename
        - self.time_stamp
        - self.r_axis[1]: (population time)

        """
        filename = [0] * 4  # the filenames

        filename[0] = self.path + self.base_filename + "_R1" + "_" + str(scan) + ".bin"
        filename[1] = self.path + self.base_filename + "_R2" + "_" + str(scan) + ".bin"
        filename[2] = self.path + self.base_filename + "_NR1" + "_" + str(scan) + ".bin"
        filename[3] = self.path + self.base_filename + "_NR2" + "_" + str(scan) + ".bin" 

        return filename



    def make_arrays(self):
        """
        croc.Pe.pefs.make_arrays()

        Creates correctly sized arrays for self.b, self.b_count, self.b_axis and self.r.

        INPUT:
        - none

        IMPLICIT REQUIREMENTS:
        In the data structure, the following variables should be set:
        - self.n_fringes
        - self.extra_fringes
        - self.n_channels

        """

        self.verbose("make arrays", True)

        self.b = [numpy.zeros((self.n_fringes + 2 * self.extra_fringes, self.n_channels), dtype = "cfloat"),numpy.zeros((self.n_fringes + 2 * self.extra_fringes, self.n_channels), dtype = "cfloat"),numpy.zeros((self.n_fringes + 2 * self.extra_fringes, self.n_channels), dtype = "cfloat"),numpy.zeros((self.n_fringes + 2 * self.extra_fringes, self.n_channels), dtype = "cfloat")] 
        # from -20 to 920
        self.b_axis[0] = numpy.arange(- self.extra_fringes, self.n_fringes + self.extra_fringes)
        # from -20 to 920
        self.b_axis[1] = numpy.arange(- self.extra_fringes, self.n_fringes + self.extra_fringes)
        self.b_count = numpy.zeros((4, self.n_fringes + 2 * self.extra_fringes))
        self.r = [numpy.zeros((self.n_fringes, 32), dtype = "cfloat"),numpy.zeros((self.n_fringes, 32), dtype = "cfloat")] 





    def import_binned_data(self, path_and_filename, diagram):

        b, b_count, b_axis = self.import_binned_data_2(path_and_filename, self.n_pixels, diagram)

        # if b etc are integers, the importing went wrong
        if type(b) == int:
            return False

        # if numpy.shape(self.b_axis)[-1] == 2:
        if self.b == [0,0]:
            self.make_arrays()     

        for i in range(len(b_axis)):
            # rephasing
            if diagram == 0:
                index_self = int(b_axis[i] - 4000 + self.extra_fringes)
            elif diagram == 1:
                index_self = int(-(b_axis[i] - 4000) + self.extra_fringes)

            for j in range(4):
                self.b[j][index_self, :self.n_pixels] += b[j][i,:]

                if diagram == 0:
                    if j == 0 or j == 1:
                        self.b_count[j][index_self] += b_count[j][i]
                else:
                    if j == 2 or j == 3:
                        self.b_count[j][index_self] += b_count[j-2][i]

        return True


    def construct_r(self, flag_no_logarithm = False):
        """
        croc.Pe.pefs.construct_r()

        This function will construct the rephasing and non-rephasing diagrams. 

        It will first average the data. Then it will select the part where the fringes are not negative. Then it will take the difference between the two PEM-states and calculate the optical density. Data points that are not-a-number will be converted to zero.


        IMPLICIT REQUIREMENTS:
        In the data structure, the following variables should be set:
        - self.n_fringes
        - self.extra_fringes
        - self.b
        - self.b_count
        - self.b_axis
        - self.n_pixels
        - self.r: this will be overwritten


        OUTPUT:
        - self.r_axis[0]: the times
        - self.r
        - self.r_units

        """

        b_fringes = self.n_fringes + 2 * self.extra_fringes

        temp = numpy.zeros((4, b_fringes, self.n_channels), dtype = "cfloat")

        # average the data for the two diagrams
        i_range = range(b_fringes)
        for j in range(4):
            for i in i_range:
                if self.b_count[j][i] != 0:
                    temp[j,i,:] = self.b[j][i,:] / self.b_count[j,i]    
                else:    
                    temp[j,i,:] = numpy.zeros(self.n_channels)                

        # select n_fringes, ie. discard the extra fringes
        temp = temp[:,self.extra_fringes:(self.n_fringes+self.extra_fringes),:self.n_pixels]

        # make the r_axis
        self.r_axis[0] = self.b_axis[0][self.extra_fringes:(self.n_fringes+self.extra_fringes)] * CONST.hene_fringe_fs

        if flag_no_logarithm:
            # for testing purposes
            for j in range(2):
                self.r[j][:,:self.n_pixels] = temp[2*j,:,:self.n_pixels] - temp[2*j+1,:,:self.n_pixels]

        else:
            # convert it to mOD
            for j in range(2):
                self.r[j][:,:self.n_pixels] = -numpy.log10(1+ 2*(temp[2*j,:,:self.n_pixels] - temp[2*j+1,:,:self.n_pixels])/self.reference[:self.n_pixels]) 

        self.r_units = ["fs", "fs", "cm-1"]


    def reconstruct_counter(self, data, x_channel, y_channel, start_counter = 0, end_counter = 0):
        """
        croc.Resources.PEFunctions.reconstruct_counter()
    
        This function will use the feedback from the HeNe's and reconstruct the fringes. It will check whether y > 0 and whether x changes from x[i-1] < 0 to x[i+1] > 0 and whether x[i-1] < x[i] < x[i+1] (or the other way around for a count back). 
        After a count in a clockwise direction, it can only count again in the clockwise direction after y < 0. 
    
        INPUT:
        - data (2darray, channels x samples): data. 
        - start_counter (int, 0): where the count starts. 
        - flag_plot (BOOL, False): plot the x and y axis and the counts. Should be used for debugging purposes.
        - x_channel, y_channel (int): the channel with the x and y data
    
        OUTPUT:
        - m_axis (1d-ndarray, length of samples): the exact fringe for that sample
        - counter (int): the last value of the fringes. It is the same as m_axis[-1]. For legacy's sake.
    
    
        CHANGELOG:
        20110920 RB: started the function
        20111003 RB: change the way it counts. It will now not only check if the x goes through zero, it will also make sure that the point in between is actually in between. This reduced the miscounts from 80/400 t0 30/400.
        20120227 RB: moved the function to croc.Resources.PEFunctions
        20130201/RB: Crocodile, moved it to legacy functions.
    
    
        """
    
    
        # put the required data in some better readable arrays
        x = data[x_channel,:]
        y = data[y_channel,:]
    
        # determine the median values
        med_x = numpy.min(x) + (numpy.max(x) - numpy.min(x))/2
        med_y = numpy.min(y) + (numpy.max(y) - numpy.min(y))/2
    
        # some stuff
        length = len(x)
        counter = start_counter
    
        # the fringe count will be written in this array
        m_axis = numpy.zeros(length)
    
        # this is the (counter-) clockwise lock
        c_lock = False
        cc_lock = False
    
        count_c = 0
        count_cc = 0
    
        count_limit = 20
        length_limit = 1500
    
        # where did the counter start
        m_axis[0] = counter
    
        # do the loop
        i_range = range(1, length - 1)
        for i in i_range:
            # count can only change when y > 0
            if y[i] > med_y:
                if c_lock == False:
                    if x[i-1] < med_x and x[i] > med_x: 
                        # add 1 to the counter
                        counter += 1        
    
                        # lock clockwise count       
                        c_lock = True
    
                        count_c += 1
    
                        # if we are the beginning or end, unlock  counterclockwise 
                        if count_c < count_limit or i > length - length_limit:
                            cc_lock = False
                        else:
                            cc_lock = True
    
                if cc_lock == False:
                    if x[i-1] > med_x and x[i] < med_x:
                        counter -= 1
    
                        cc_lock = True  
    
                        count_cc += 1
    
                        if count_cc < count_limit or i > length - length_limit:
                            c_lock = False
                        else:
                            c_lock = True
    
            else:
                # we are at y<0
                # only unlock clockwise if we have more than 100 counts up
                if count_c > count_limit:              
                    c_lock = False
                if count_cc > count_limit:
                    cc_lock = False
                # or if we are the beginning or end
                if i < length_limit or i > length - length_limit:
                    c_lock = False
                    cc_lock = False
    
            m_axis[i] = counter
    
        m_axis[-1] = counter
    
        if counter == end_counter:
            correct_count = True
        else:
            correct_count = False
    
        return m_axis, counter, correct_count  





    def import_reference_VB6(self, path_and_filename, n_pixels):
    
        try:       
            data = numpy.loadtxt(path_and_filename)
        except IOError:
            self.printError("Unable to load reference", inspect.stack())
            return 0,0
    
        r_axis_2 = data[0,1:(n_pixels+1)]
        reference = data[1,1:]
    
        return r_axis_2, reference
    
    
    
    def import_meta_VB6(self, path_and_filename):
        """
        croc.pe.import_meta
    
        Import data from the meta-data files.
    
        INPUT: 
        - a class
    
        COMMENTS
        - it will scan the meta-file for n_scans, n_shots, phase and comments
    
        CHANGELOG
        RB 20110909: started function        
    
    
        """
    
        temp_scans = 0
        date = False
        data_type_version = False
        n_shots = False
        n_fringes = False
        phase_degrees = numpy.nan
        comment = ""
    
        # found at: http://docs.python.org/library/re.html#matching-vs-searching , is equivalent to scanf(..%f..)
        regex = re.compile("[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?")
    
        # scan the file
        try:
            F = open(path_and_filename)
            for line in F:
    
                if re.match("Start_time", line):
    
                    date = int(str(line[11:15]) + str(line[16:18]) + str(line[19:21]))
    
                if re.match("mess2Dheterodyne_meta_format", line): 
                    data_type_version = str(line[-5:-2])               
    
                if re.match("Shots", line): 
                    n_shots = int((re.search(regex, line)).group())
    
                if re.match("Fringes", line): 
                    n_fringes = int((re.search(regex, line)).group())            
    
                try:
                    if re.match("Phase", line):
                        phase_degrees = float((re.search(regex, line)).group())
                except AttributeError:
                    DEBUG.printWarning("No or invalid phase", inspect.stack())
    
                if re.match("Comments", line):
                    comment = line[9:]
    
                if re.match("Scan", line):
                    temp_scans = int((re.search(regex, line[4:7])).group())
    
            F.close()
        except IOError:
            self.printError(("Unable to load file: " + path_and_filename), inspect.stack())
    
            return 0,0,0,0,0,0,0         
    
    
    
        # number of scans is (number of scans started) - 1, because the last one wasn't finished
        # if self.mess_type != "FastScan":
        #    if temp_scans:
        #        self.n_scans = temp_scans - 1
        #    else:
        #        self.n_scans = 1
    
        return date, data_type_version, n_shots, n_fringes, phase_degrees, comment, temp_scans
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def import_data_FS(self, path_and_filename, n_shots = 30000, n_channels = 37, flag_counter = False):
        """
        This method is a derivative of croc.Pe.PeFS.import_raw_data, but without the reliance on the class structure.
    
        INPUT:
        - path_and_filename (string): where the file can be found
        - n_shots (int): number of shots
        - n_channels (int): number of channels
    
        OUTPUT:
        - m (2d ndarray): array with (channels x shots)
        - fringes (array with two elements): the begin and end fringe, as appended to the data.
    
        CHANGELOG:
        - 201110xx RB: wrote function
        - 201202xx RB: rewrote the function for wider purpose, moved it to IOMethods
    
        """    
        try:
            data = numpy.fromfile(path_and_filename)
    
            # remove the two fringes at the end
            fringes = [data[-2], data[-1]]
            data = data[:-2]
    
            # construct m
            m = numpy.zeros((n_channels, n_shots), dtype = "cfloat")
    
            if flag_counter:
                c = data[-n_shots:]
                if c[-1] == 0:
                    # this is to repair a (now fixed) bug from VB6. 
                    c = data[-n_shots-1:-1]
    
                data = data[:-n_shots]
    
            # order the data in a 2d array
            i_range = range(n_shots)
            j_range = range(n_channels)
            for i in i_range:
                for j in j_range:
                    m[j, i] = data[j + i * n_channels] 
    
            if flag_counter:
                return m, c, fringes
            else:
                return m, fringes
    
        except IOError:
            self.printError("Unable to import raw data from file " + path_and_filename, inspect.stack())
            if flag_counter:
                return 0,0,0
            else:
                return 0,0
    
    
    
    def import_binned_data_2(self, path_and_filename, n_pixels, diagram):
    
        try:
            data = numpy.fromfile(path_and_filename)
    
            fringes = [data[-2], data[-1]]
            b_axis = numpy.arange(fringes[0], fringes[1] + 1)
            n_fringes = len(b_axis)
    
            # without the two fringes at the end, the number of elements should be equal to (2*n_pixels + 1), for each chopper state the pixels and the count
    
            if ((len(data)-2)/n_fringes) == (2 * n_pixels + 2):
    
                b_count = [0]*2
                b_count[0] = data[-2 - 2*n_fringes:-2 - n_fringes]
                b_count[1] = data[-2 - n_fringes:-2]
    
                data = data[:2*n_pixels*n_fringes]         
    
                # rearrange the data
                b = numpy.zeros((4, n_fringes, n_pixels), dtype = "cfloat")
    
                for i in range(2): # the two chopper states
                    for j in range(n_fringes):
                        for p in range(n_pixels):
                            b[i+diagram*2,j,p] = data[i * n_fringes * n_pixels + j * n_pixels + p]
    
                return b, b_count, b_axis
    
            else:
                self.printWarning("The file does not contain binned data: " + path_and_filename, inspect.stack())
                return 1, 1, 1
    
        except IOError:
            self.printError("Unable to import binned data from file " + path_and_filename, inspect.stack())
            raise
            return 0, 0, 0   

