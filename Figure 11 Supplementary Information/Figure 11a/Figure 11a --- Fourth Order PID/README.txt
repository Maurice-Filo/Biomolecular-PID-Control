#### Generate Supplementary Figure 11a (Fourth Order PID) ####

# Software Requirements
 - MATLAB R2016a (academic use) or above (Tested on R2021a)

# System Requirements
 - OS: Windows 7 or above (Tested on Windows 7 and 10)
       macOS 10.12 or above (Tested on macOS 10.15.7)

# Installation Requirements
 - No separate installation required, please run the scripts in MATLAB environment

# Instructions
 - Open "APIDF4_Deg_GeneExp_Variance.m" in MATLAB environment and set the variable "Load" to 0
 - Run "APIDF4_Deg_GeneExp_Variance.m" to generate the simulated data saved as "APIDF4_Deg_GeneExp_Variance_Sweep.mat"
 	> Expected run time: ~250-300 seconds
 - Set the variable "Load" to 1 and rerun "APIDF4_Deg_GeneExp_Variance.m" to generate the figures.

# NOTE 
The accuracy and smoothness of the mean and variance plots is dictated by the grid resolution N of the derivative gain K_D, the number of (parallel) stochastic trajectories N_Trajectories and Iterations and the maximum allowed number of events N_MaxRxnEvents_Long. Higher values yield more accurate and smoother plots at the expense of higher computation time. In the attached Matlab file, the numbers are chosen to demonstrate the simulation results quickly with low accuracy. In the paper the following values were used:
 	N_Trajectories = 64
	Iterations = 512
	N_MaxRxnEvents_Long = 1e8
	N = 50


   