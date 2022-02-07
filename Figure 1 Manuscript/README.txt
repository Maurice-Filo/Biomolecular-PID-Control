#### Generate manuscript Figure 1b ####

# Software Requirements
 - MATLAB R2016a (academic use) or above (Tested on R2021a)

# System Requirements
 - OS: Windows 7 or above (Tested on Windows 7 and 10)
       macOS 10.12 or above (Tested on macOS 10.15.7)

# Installation Requirements
 - No separate installation required, please run the scripts in MATLAB environment

# Instructions
 - Open "Trajectories_Performance.m" in MATLAB environment and set the variable "Load" to 0
 - Run "Trajectories_Performance.m" to generate the simulated data saved as "Simulations_Concept.mat"
 	> Expected run time with N_SamplePaths = 100: ~5-10 seconds
 - Set the variable "Load" to 1 and rerun "Trajectories_Performance.m" to generate the figures.

# NOTE 
The smoothness of the mean and variance plots is dictated by the number of stochastic trajectories stored in the variable N_SamplePaths. Higher values for N_SamplePaths yield a smoother plot at the expense of higher computation time. In the paper, N_SamplePaths was set to 1000. 

   