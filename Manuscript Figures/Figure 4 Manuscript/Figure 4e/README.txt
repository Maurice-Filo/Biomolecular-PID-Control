#### Generate manuscript Figure 4e ####

# Software Requirements
 - MATLAB R2016a (academic use) or above (Tested on R2021a)

# System Requirements
 - OS: Windows 7 or above (Tested on Windows 7 and 10)
       macOS 10.12 or above (Tested on macOS 10.15.7)

# Installation Requirements
 - No separate installation required, please run the scripts in MATLAB environment

# Instructions
 - Open "Generate_APIF1_StochasticSimulations" in MATLAB environment and set the variable "Load" to 0
 - Run "Generate_APIF1_StochasticSimulations" to generate the simulated data saved as "Results.mat"
 	> Expected run time: ~80-100 seconds
 - Set the variable "Load" to 1 and rerun "Generate_APIF1_StochasticSimulations" to generate the figures.

# NOTE 
The accuracy and smoothness of the mean and variance plots is dictated by the grid resolution N_Simulations of the corresponding parameter, the number of stochastic trajectories N_Trajectories used to compute the dynamics of the variance, the final time tf_Long used to compute the stationary variance for each parameter value on the grid and the maximum allowed number of events N_MaxRxnEvents_Long. Higher values yield more accurate and smoother plots at the expense of higher computation time. In the attached Matlab file, the numbers are chosen to demonstrate the simulation results quickly with low accuracy. In the paper the following values were used:
     	N_Simulations = 100
 	N_Trajectories = 10000	
	tf_Long = 2e4
	n_MaxRxnEvents_Long = 4e7


   