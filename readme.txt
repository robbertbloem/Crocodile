CROCODILE README
A collection of code to process measurements


0 - DISCLAIMER
This software is experimental and may come with bugs. Please report those.


1 - GOAL ETC

Crocodile is an improved version of the croc scripts. Since croc was written a lot has changed and Crocodile is designed from the start to work with those changes. 
Most notable is the inclusion of a separate class to keep track of the array of measurements. In croc this was scattered around. These tools are part of PythonTools.
The classes have also been reorganized. 
A lot of the deeper code (the actual calculations) has not been changed (at least, not too much).
The changes mean that the scripts used until now are not compatible.

2 - INSTALLATION/SETUP:

- make sure you have Numpy, Scipy, Matplotlib and ipython installed. This is done most conveniently using the "Enthought Python Distribution" (Google it). There is an academic license and versions for Windows, Mac and Linux. 
- copy everything to a convenient place
- add that place to the PYTHONPATH in for example .bash_profile: 
export PYTHONPATH=$HOME/Developer:$HOME/Developer/Crocodile:$PYTHONPATH 


3 - ORGANIZATION

3.1 - OBJECTARRAY


3.2 - MEASUREMENTS

The highest class is DataClass. This has been scaled down from croc, a lot of the experimental details have been moved to the classes where experiments are described.

PE is the photon echo class. It contains a method to plot 2D contour plots (plot) and to calculate the absorptive spectrum (super_absorptive). 

PE has 4 subclasses: pe_tw, pe_wt, pe_tt and pe_ww, where 'w' stands for frequency and 't' for time. pe_tw has t_1 and omega_3, for example. These subclasses determine how the absorptive spectrum is calculated. (1D-FFT over one or the other axis, or 2D)
pe_tw is for the measurements
pe_ww can be used when combining already-calculated spectra
pe_tt is for simulations
pe_wt is to complete the set, although I don't see a purpose for it

pe_tw has subclasses for the different types of measurements. 


3.3 - DEBUG

The debug class contains the methods:
- verbose: will print stuff (in blue) when the program runs. The function is always called, flag_verbose will determine if it will print the message. For complicated messages, you might want to do this before you call the function, but otherwise it cleans up a lot of if-statements in the code.
- printWarning, printError: will print stuff (purple or red) as warnings or errors. Both indicate an exception, after warning means that the script can continue running, an error means that it will stop. Part of the argument is location = inspect.stack(), which will print the file, function and line number where the error occurred. 














