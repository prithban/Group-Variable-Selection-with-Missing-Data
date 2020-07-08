# Group Variable Selection with Missing Data
This is the instruction file to run the simulations for 'A Comparative Study of Variable Selection in presence of Missing Data via Different Grouped Penalties'.

## Overview of the R codes
* ***degree_freedom.R*** : function for calculating degrees of freedom for each penalty under weighted least squares
* ***degree_freedom_noweight.R*** : function for calculating degrees of freedom for each penalty without any weights
* ***parallel_compute.R*** : code for running the simulation with stacked multiple imputed data in parallel
* ***parallel_compute_fulldata.R*** : code for running the simulation with full data in parallel
* ***sim_aoas_data_generate.R*** : function for generating the data for simulation
* ***sim_imput_combined.R*** : function for obtaining results of the simulated multiple imputed data for all four data weighing scheme
* ***sim_imput_combined.R*** : function for obtaining results of the simulated full data for all four data weighing scheme

## Instructions for Simulation Run
* Install the packages required for running the codes. The packeages are ***parallel***, ***MASS***, ***mice***, ***mpath***, and ***grpreg***.
* Run the ***parallel_compute_fulldata.R*** code to obtain the results for full data in Tables 1 and 2
* Run the ***parallel_compute.R*** code to obtain the results for stacked multiple imputed data in Tables 3 to 12
