#### Generate manuscript Figure 7b ####

# Software Requirements
 - MATLAB R2016a (academic use) or above (Tested on R2021a)

# System Requirements
 - OS: Windows 7 or above (Tested on Windows 7 and 10)
       macOS 10.12 or above (Tested on macOS 10.15.7)

# Installation Requirements
 - No separate installation required, please run the scripts in MATLAB environment

# Instructions
 - Open "APIDF2_Deg_Star_PerformanceTracking.m" in MATLAB environment and set the variable "Load" to 0
 - Run "APIDF2_Deg_Star_PerformanceTracking.m" to generate the simulated data saved as "APIDF2_Deg_Star_PerformanceTracking_Sweep.mat"
 	> Expected run time with N = 30: ~15-25 seconds
 - Set the variable "Load" to 1 and rerun "PIDF2_Deg_Star_PerformanceTracking.m" to generate the figures.

# NOTE 
The smoothness of the performance intensity plots is dictated by the 2D grid resolution N of the corresponding parameters. Higher values for N yield smoother plots at the expense of higher computation time. In the paper, N was set to 300. 

   