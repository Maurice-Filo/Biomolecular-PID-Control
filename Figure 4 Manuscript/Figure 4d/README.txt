#### Generate manuscript Figure 4d ####

# Software Requirements
 - MATLAB R2016a (academic use) or above (Tested on R2021a)

# System Requirements
 - OS: Windows 7 or above (Tested on Windows 7 and 10)
       macOS 10.12 or above (Tested on macOS 10.15.7)

# Installation Requirements
 - No separate installation required, please run the scripts in MATLAB environment

# Instructions
 - Open "Generate_APIF1_Simulations.m" in MATLAB environment and set the variable "Load" to 0
 - Run "Generate_APIF1_Simulations.m" to generate the simulated data saved as "APIF1_DeterministicSimulations.mat"
 	> Expected run time with N = 100: ~15-20 seconds
 - Set the variable "Load" to 1 and rerun "Generate_APIF1_Simulations.m" to generate the figures.

# NOTE 
The smoothness of the settling time and overshoot plots is dictated by the grid resolution N of the corresponding parameter. Higher values for N yield smoother plots at the expense of higher computation time. In the paper, N was set to 200. 

   