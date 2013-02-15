CROCODILE README
A collection of code to process measurements


0 - DISCLAIMER
This software is experimental and may come with bugs. Please report those.


1 - GOAL ETC

Crocodile is an improved version of the croc scripts. Since croc was written a lot has changed and Crocodile is designed from the start to work with those changes. 
Most notable is the inclusion of a separate class to keep track of the array of measurements. In croc this was scattered around. These tools are part of PythonTools.
The classes have also been reorganized. 
A lot of the deeper code (the actual calculations) have not changed (at least, not too much).
The changes mean that the scripts used until now are not compatible.

2 - INSTALLATION/SETUP:

- make sure you have PythonTools installed
- make sure you have Numpy, Scipy, Matplotlib and ipython installed. This is done most conveniently using the "Enthought Python Distribution" (Google it). There is a full version with academic license. It also works with the free version. There are versions for Windows, Mac and Linux. 
- copy everything to a convenient place
- add that place to the PYTHONPATH in for example .bash_profile: 
export PYTHONPATH=$HOME/Developer:$HOME/path_to_crocodile/Crocodile:$PYTHONPATH 


3 - ORGANIZATION


3.1 - DATACLASS, PE

The highest class is DataClass. This has been scaled down from croc, a lot of the experimental details have been moved to the classes where experiments are described.

PE is the photon echo class. It contains a method to plot 2D contour plots (plot) and to calculate the absorptive spectrum (super_absorptive). 

PE has 5 subclasses: pe_tw, pe_wt, pe_tt, pe_ww and pe_merge, where 'w' stands for frequency and 't' for time. pe_tw has t_1 and omega_3, for example. These subclasses determine how the absorptive spectrum is calculated. (1D-FFT over one or the other axis, or 2D)
pe_tw is for the measurements
pe_ww can be used when combining already-calculated spectra
pe_tt is for simulations
pe_wt is to complete the set, although I don't see a purpose for it

pe_tw has subclasses for the different types of measurements, for example pe_LV for LabView data, pe_VB6 for Visual Basic 6. Because the VB6 code is not used anymore and there is quite a lot of it, it has been collected and put in a separate file.


3.2 - RESOURCES, PLUGINS

The folder "Resources" contains some code that supports DataClass and Pe, think about plotting standard 2DIR-contourplots. Other modules depend on the code here and changes should be made with care.
The folder "Plugins" are modules that depend on the DataClass, Pe and Resources modules, but are more hacked together for a special purpose. Think about a special routine that makes a fancy plot with two overlapping 2DIR-measurements. If the method is widely used, it can be moved to (a module in) the "Resources" folder. 


3.3 - TESTS

Most modules come with tests to ensure:
- that a function works correctly with the correct input
- that a functions elegantly handles edge-cases and incorrect input
- to ensure that the behavior of functions does not change when functions are modified. 
Test should be made for the smallest possible functions, from the start. The tests describes the behavior of the function: what to do when the input is this or that. 
In general, tests should only be added. When test are removed or changed the behavior may change, leading to errors elsewhere.
To put it differently: changing tests is usually not the correct way to fix a test that gives an error. 
Read more about test driven development: https://en.wikipedia.org/wiki/Test-driven_development


3.4 - CODE SNIPPETS

The module code_snippets.py contains some useful functions and short descriptions. 


4 - DEVELOPMENT NOTES

4.1 - TESTS

Use tests, as described in 3.3.


4.2 - PYTHON 3 COMPATIBILITY

Summary: the scripts were developed in Python 2.7 and were tested with 3.3. The scripts don't work in Python 3.0-3.2.

The scripts are written in Python 2.7. Using the following __future__ imports the scripts have maximum compatibility for Python 3.
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

A fourth __future__ import
from __future__ import unicode literals
is not used. It is a bag of hurt to make sure something is a string or unicode. From Python 3.3 this is changed back to the situation of 2.7. 

Note that 
from imp import reload 
is needed to reload modules in Python 3.

Note that I wasn't able to install Scipy in Python3. This is only used in the Crocodile.Resources.Mathematics.fit() function. A flag indicates if Scipy is imported. If not, the function will return an error.



















