from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import inspect
import re
import os

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import PythonTools.Debug as DEBUG
import Crocodile.Resources.Constants as CONST

def import_data_LV_A(path, base_filename):
    """
    Imports data for 'LV_file_format.1'
    
    """

    try:
        # data
        path_and_filename = path + base_filename + ".csv"
        data = numpy.loadtxt(path_and_filename, dtype = "float", delimiter = ",")
        data = data.T

        # time axis in fringes
        path_and_filename = path + base_filename + "_t1.csv"
        t1_axis = numpy.loadtxt(path_and_filename, dtype = "float", delimiter = ",") 

        # frequency axis
        path_and_filename = path + base_filename + "_w3.csv"
        w3_axis = numpy.loadtxt(path_and_filename, dtype = "float", delimiter = ",") 

        # phase
        path_and_filename = path + base_filename + "_phase.txt"
        f = open(path_and_filename)
        for line in f:
            temp = line
        f.close()
        if temp == "NaN":
            DEBUG.printWarning("Phase is Not-a-Number, will be set to 0 ", inspect.stack())
            phase = 0
        else:
            phase = float(temp)

        # last pump
        path_and_filename = path + base_filename + "_lastpump.txt"
        f = open(path_and_filename)
        for line in f:
            lastpump = line
        f.close()

        # determine number of fringes
        n_fringes = int((len(t1_axis)+1)/2)
        n_pixels = len(w3_axis)

        # convert NaN to zeros
        data = numpy.nan_to_num(data)

        # labview measures 4000-N to 4000+N, we want the data split into 4000-N to 4000 (non-rephasing) and 4000 to 4000+N (rephasing)
        R = data[n_fringes-1:, :]
        NR = numpy.flipud(data[:n_fringes, :])

        # for the FFT, we don't want 4000 to be zero. The axes run from 0 to N
        # also: calculate the axis in fs        
        t1fr_axis = numpy.arange(n_fringes)
        t1fs_axis = numpy.arange(n_fringes) * CONST.hene_fringe_fs

        # return everything
        return R, NR, t1fs_axis, t1fr_axis, w3_axis, phase, lastpump, n_fringes, n_pixels

    except IOError:
        DEBUG.printError("Unable to import LabView data from file " + path + base_filename, inspect.stack())
        # raise
        return False
        

def save_data_PE(path, base_filename, s = False, s_axis = False, r = False, r_axis = False):
    
    if type(s) != bool:
        path_and_filename = path + base_filename + "_s.csv"
        numpy.savetxt(path_and_filename, s, delimiter = ",")
    
    if type(s_axis) != bool:
        path_and_filename = path + base_filename + "_s_w1.csv"
        numpy.savetxt(path_and_filename, s_axis[0], delimiter = ",")        
        path_and_filename = path + base_filename + "_s_w3.csv"
        numpy.savetxt(path_and_filename, s_axis[2], delimiter = ",")  

    if type(r) != bool:
        path_and_filename = path + base_filename + "_R.csv"
        numpy.savetxt(path_and_filename, r[0], delimiter = ",")
        path_and_filename = path + base_filename + "_NR.csv"
        numpy.savetxt(path_and_filename, r[1], delimiter = ",")
    
    if type(r_axis) != bool:
        path_and_filename = path + base_filename + "_r_t1.csv"
        numpy.savetxt(path_and_filename, r_axis[0], delimiter = ",")        
        path_and_filename = path + base_filename + "_r_w3.csv"
        numpy.savetxt(path_and_filename, r_axis[2], delimiter = ",")      

def check_and_make_list(var, flag_verbose = False):
    """
    numpy.loadtxt returns a 1 or more dimensional array. Here I make some sense from it. It also calculates the lengths of arrays. 
    
    if var is a python-list, return the list and the length
    if var is an ndarray
    
    """

    if flag_verbose:
        print("IOMethods.check_and_make_list:")
        print(var)

    if type(var) == "list":
        n_var = len(var)
    elif type(var) == numpy.ndarray:
        n_var = numpy.shape(var)
        if len(n_var) == 1:
            n_var = n_var[0]
    else:
        var = numpy.array([var])
        n_var = 1

    return var, n_var 
 
def find_LV_fileformat(base_folder, flag_verbose = False, test_input = False):
    """
    Scans the directory for the LV_fileformat file. 
    
    INPUT:
    - base_folder: folder with the measurement
    
    OUTPUT:
    - fileformat_version, an integer
    
    
    CHANGELOG:
    20151028-RB: function uses regex
    
    """

    fileformat_version = -1
    
    p = re.compile(r"LV_fileformat\.[0-9]*\Z")

    if test_input == False:
        files = os.listdir(base_folder)
    else:
        if type(test_input) == str:
            files = [test_input]
        elif type(test_input) == list:
            files = test_input
            
    for file in files:
        m = p.search(file)
        if m != None:
            fileformat_version = int(m.group()[14:])

    if flag_verbose:
        print("IOMethods.find_LV_fileformat: %i" % fileformat_version)     

    return fileformat_version

def find_number_of_scans(base_folder, base_filename, extension, flag_verbose = False, test_input = False):
    """
    Looks for the last number in the file string. It will return the largest number + 1 (because we start at zero). 
    """

    if extension[0] == ".":
        extension = "\\" + extension
    else:
        extension = r"\\\." + extension

    n_scans = -1
    
#     s = base_filename + 
    s = r"_[0-9]*" + extension + "\Z"
    p1 = re.compile(s)
    s = r"[0-9]*" + extension + "\Z"
    p2 = re.compile(s)

    if test_input == False:
        files = os.listdir(base_folder)
    else:
        if type(test_input) == str:
            files = [test_input]
        elif type(test_input) == list:
            files = test_input

    for file in files:
        m1 = p1.search(file)
        if m1 != None:
            m2 = p2.search(m1.group(0))
            if m2 != None:
                temp = int(m2.group(0)[:-(len(extension)-1)])
                if temp > n_scans:
                    n_scans = temp

    n_scans += 1

    if flag_verbose:    
        print("IOMethods.find_number_of_scans: number of scans found: %i" % n_scans)

    return n_scans

def find_number_of_datastates(base_folder, flag_verbose = False, test_input = False):

    n_ds = -1
    
    s = r"_ds[0-9]*_"
    p1 = re.compile(s)

    if test_input == False:
        files = os.listdir(base_folder)
    else:
        if type(test_input) == str:
            files = [test_input]
        elif type(test_input) == list:
            files = test_input

    for file in files:
        m1 = p1.search(file)
        if m1 != None:
            temp = int(m1.group(0)[3:-1])
            if temp > n_ds:
                n_ds = temp
    n_ds += 1 
    return n_ds    

def check_basename_extension_suffix(file_dict, extension, suffix, flag_verbose = False):
    """
    Normalizes the basename, extension and suffix. 
    """
    
    if extension[0] != ".":
        extension = "." + extension

    if basename[-1] != "_":
        basename += "_"

    if suffix[0] == "_":
        suffix = suffix[1:]
    
    if suffix[-1] == "_":
        suffix = suffix[:-1]
        
    return basename, extension, suffix


def import_file(file_dict, suffix, flag_verbose = False):
    """
    This does the actual importing. It is very general and does no checking.
    """

#     basename, extension, suffix = check_basename_extension_suffix(basename, extension, suffix, verbose)

    filename = file_dict["base_filename"] + "_" + suffix +  file_dict["extension"]
    
    if flag_verbose == True or flag_verbose > 0:
        print("IOMethods.import_file: %s" % filename)
    
    data = numpy.loadtxt(filename, delimiter = ",", ndmin = 1)

    if flag_verbose == True or flag_verbose > 0:
        print("IOMethods.import_file: done")

    return data  

def import_bins(file_dict, fileformat, t1_offset = 0, flag_verbose = 0):
    """
    Import the bins. Column 0 are the bins, column 1 are the associated times. 
    
    The bins are actually boring: they start at 0 and continue to bin N-1. 
    
    INPUT:
    t1_offset (int, 0): if t0 is not where it should be, adjust it. The value is in integer and is the number of bins. 
    
    OUTPUT:
    t1_bins (1D ndarray): list with bins
    t1_fs (1D ndarray): list with times
    bin_sign (Boolean): The times can run in two directions: from low to high (bin_sign = False) or from high to low (bin_sign = True). In these scripts it is assumed that the time runs from low to high. When importing fast scans bin_sign should be checked and, when needed, the time axis should be reversed. 
    n_t1_bins (int): number of bins 
    n_t1_fs (int): number of steps where time >= 0
    t1_zero_index: index of the bin where t1 = 0
    
    """
    if flag_verbose:
        print("IOMethods.import_bins")

    suffix = "bins"  
    data = import_file(file_dict, suffix, flag_verbose)
    data, shape = check_and_make_list(data, flag_verbose)

    # extract the bins and the times
    t1_bins = data[:,0]
    t1_fs = data[:,1] + t1_offset * CONST.hene_fringe_fs
    
    # if the bins and times are opposite, reverse the time axis. 
    # set bin_sign to True, reverse other time axes as well.
    if t1_fs[0] > t1_fs[1]:
        bin_sign = True
        t1_fs = t1_fs[::-1]
    else:
        bin_sign = False

    # bin where t1 = 0
    t1_zero_index = numpy.where(t1_fs >= 0)[0][0]
    
    # length of bins and length 
    n_t1_bins = len(t1_bins)
    n_t1_fs = n_t1_bins - t1_zero_index
    
    temp = int(n_t1_fs / 2)
    n_t1_fs = 2 * temp

    if flag_verbose:
        print("t1_bins: {a}, {b}, ..., {c}".format(a=t1_bins[0], b=t1_bins[1], c=t1_bins[-1]))
        print("t1_fs: {a}, {b}, ..., {c}".format(a=t1_fs[0], b=t1_fs[1], c=t1_fs[-1]))
        print("bin_sign: {bin_sign}".format(bin_sign = bin_sign))
        print("n_t1_bins: {n_t1_bins}".format(n_t1_bins = n_t1_bins))
        print("n_t1_fs: {n_t1_fs}".format(n_t1_fs = n_t1_fs))
        print("t1_zero_index: {t1_zero_index}".format(t1_zero_index = t1_zero_index))
        print("DONE: IOMethods.import_bins\n")
    
    return t1_bins, t1_fs, bin_sign, n_t1_bins, n_t1_fs, t1_zero_index

def import_nspectra(file_dict, fileformat, flag_verbose = 0):
    """
    Import the number of spectra
    """
    if flag_verbose:
        print("IOMethods.import_nspectra")
    suffix = "Nspectra"
    data = import_file(file_dict, suffix, flag_verbose)
    n_spectra, shape = check_and_make_list(data, flag_verbose)  
    if flag_verbose:
        print("n_sp: {n_sp}".format(n_sp = int(n_spectra)))
        print("DONE: IOMethods.import_nspectra\n")
    return int(n_spectra)


def import_ndatastates(file_dict, fileformat, flag_verbose = 0):
    """
    Import the number of datastates
    """
    if flag_verbose:
        print("IOMethods.import_ndatastates")
    suffix = "Ndatastates"
    data = import_file(file_dict, suffix, flag_verbose)
    n_ds, shape = check_and_make_list(data, flag_verbose)   
    if flag_verbose:
        print("n_ds: {a}".format(a = int(n_ds)))
        print("DONE: IOMethods.import_ndatastates\n")
    return int(n_ds)


def import_wavenumbers(file_dict, fileformat, flag_verbose = 0):
    """
    Import the spectrometer axis
    """
    if flag_verbose:
        print("IOMethods.import_wavenumbers")
    suffix = "wavenumbers"
    data = import_file(file_dict, suffix, flag_verbose)
    w3_axis_wn, n_w3 = check_and_make_list(data, flag_verbose)  
    if flag_verbose:
        print("n_w3: {a}".format(a = n_w3))
        print("w3: {a}".format(a = w3_axis_wn))
        print("DONE: IOMethods.import_wavenumbers\n")
    return w3_axis_wn, n_w3

def import_delays(file_dict, fileformat, flag_verbose = 0):
    """
    Import the delays of the IR delay. 
    """
    if flag_verbose:
        print("IOMethods.import_delays")
    suffix = "delays"
    data = import_file(file_dict, suffix, flag_verbose)
    delays, n_delays = check_and_make_list(data, flag_verbose) 
    if flag_verbose:
        print("n_de: {a}".format(a = n_delays))
        print("de: {a}".format(a = delays))
        print("DONE: IOMethods.import_delays\n")
    return delays, n_delays

def import_nshots(file_dict, fileformat, flag_verbose = 0):
    """
    Import the number of shots
    """
    if flag_verbose:
        print("IOMethods.import_nshots")
    suffix = "Nshots"
    data = import_file(file_dict, suffix, flag_verbose)
    n_shots, temp = check_and_make_list(data, flag_verbose)     
    if flag_verbose:
        print("    n_sh: {a}".format(a = int(n_shots)))
        print("DONE: IOMethods.import_nshots\n")
    return int(n_shots)

def import_wavelengths(file_dict, fileformat, flag_verbose = 0):
    """
    The center wavelengths used in Scan Spectrum.
    """
    if flag_verbose:
        print("IOMethods.import_wavelengths (for Scan Spectrum)")
    suffix = "wavelengths"
    data = import_file(file_dict, suffix, flag_verbose)
    wl, n_wl = check_and_make_list(data, flag_verbose)     
    if flag_verbose:
        print("n_wl: {a}".format(a = n_wl))
        print("wl: {a}".format(a = wl))
        print("DONE: IOMethods.import_wavelengths (for Scan Spectrum)\n")
    return wl, n_wl

def import_spectraAndDatastates(file_dict, fileformat, flag_verbose = 0):

    if flag_verbose:
        print("IOMethods.import_spectraAndDatastates")

    suffix = "spectraAndDatastates"
    spds = import_file(file_dict, suffix, flag_verbose)

    if fileformat < 4:
    
        # we want it to be a 2D array
        temp = numpy.shape(spds)
        if len(temp) == 1:
            spds = numpy.array([spds])
            
        n_sig = len(numpy.unique(spds[:,1]))
        add_ds = True

    else:
        n_sig = int(spds[-1,1])
        if spds[-1,0]:
            add_ds = True
        else:
            add_ds = False
        spds = spds[:-1,:]
        

    n_sp = len(numpy.unique(spds[:,0]))
    n_ds = len(numpy.unique(spds[:,1]))

    if flag_verbose:
        print("spds: {a}".format(a = spds))
        print("n_sp: {a}".format(a = n_sp))
        print("n_ds: {a}".format(a = n_ds))
        print("n_sig: {a}".format(a = n_sig))
        print("add_ds: {a}".format(a = add_ds))
        print("DONE: IOMethods.import_spectraAndDatastates\n")

    return spds, n_sp, n_ds, n_sig, add_ds


def import_slow_modulation(file_dict, fileformat, flag_verbose = 0):
    """
    Import the slow modulation file.
    
    OUTPUT:
    sm (2D ndarray): values 
    sm_names (list): the names
    n_sm (int): number of slow modulation states
    """
    if flag_verbose:
        print("IOMethods.import_slow_modulation")
    filename = file_dict["base_filename"] + "_slowModulation" +  file_dict["extension"]
    # read all lines and strip the newline character         
    lines = [line.rstrip('\n') for line in open(filename)]
    sm, sm_names, n_sm = extract_slow_modulation(lines, fileformat, flag_verbose = flag_verbose)
    
    if flag_verbose:
        print("n_sm: {a}".format(a = n_sm))
        print("SM name: states")
        for _sm in range(len(sm_names)):
            print("{a}: {b}".format(a = sm_names[_sm], b = sm[_sm]))
        print("DONE: IOMethods.import_slow_modulation\n")    

    return sm, sm_names, n_sm



def extract_slow_modulation(lines, fileformat, flag_verbose = 0):     

    n_sm = int(lines[0])
    n_lines = len(lines)-1
    sm = numpy.zeros((n_lines, n_sm))       
    sm_names = [""] * n_lines # if we make this an ndarray, it will become a string

    p = re.compile(r'\[|\]')
    r = re.compile(r',')
    q = re.compile(r':')
    for _l in range(n_lines):
        m = p.split(lines[_l+1])            
        n = r.split(m[1])            
        sm_names[_l] = q.split(m[0])[0]    
        for _sm in range(n_sm):
            if "NaN" in n[_sm]:
                sm[_l,_sm] = numpy.nan
            else:
                sm[_l,_sm] = int(float(n[_sm]))
    
    sm_names = numpy.array(sm_names)

    idx = []
    for _l in range(n_lines):
        if not numpy.any(numpy.isnan(sm[_l])):
            idx.append(_l)

    if idx == []:
        sm = numpy.array([[]]) 
        sm_names = numpy.array([[]])
    else:
        sm = sm[idx,:] 
        sm_names = sm_names[idx]       
    
    return sm, sm_names, n_sm



if __name__ == '__main__':

    pass
#     find_LV_fileformat("")







 