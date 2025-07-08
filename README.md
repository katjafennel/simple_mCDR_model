# README

Code repository accompanying (and extending) scripts used in the article "The Verification Challenge of Marine Carbon Dioxide Removal" by Katja Fennel in Annual Review of Marine Science, Vol. 18, 2026 (https://doi.org/10.1146/annurev-marine-032123-025717). There are separate folders for case 1 (`case_01`) and case 2 (`case_02`) described in the article. The code depends on basic routines for calculating carbonate system parameters in seawater from this repository ENTER LINK.

## Cases 1.1 and 1.2

The model for case 1, described in section 4.1, represents a homogenous, fully mixed surface ocean that is in exchange with a homogenous, fully mixed atmosphere and subjected to an OAE intervention at time t = 0. Case 1.1 starts from a surface ocean that is in equilibrium with the atmosphere before being perturbed. Case 1.2 starts from a disequilibrium (surface ocean is over- or undersaturated in CO2). 

Both cases can be run from within directory `case_01` by calling the script: `>> case_01`

The script initializes a scructure of input parameters for the base case and for a number of cases that explore the sensitivity to input parameters and parameterizations. The model time stepping is carried out by the functions `f_case_01p?.m` which are called from within the script for each case. The ouput from each of these function calls is returned as a structure, plotted by calling script `plot_case_01p?.m` (if the logical variable `plot_output` is set to 1) and saved in the directory `output`. All parameters can be adjusted by editing the script `case_01.m`.  

## Case 2

The model for case 2, described in section 4.2, represents a horizontally resolved, inhomogeneous surface ocean with explicit inclusion of horizontal dispersion and vertical mixing. Because horizontal dispersion is assumed to be isotropic, the two horizontal dimensions can be collapsed into one spatial dimension which measures the distance from the injection site. 

### Case 2.1
Case 2.1 examins a spatially resolved surface ocean where seawater is again assumed to be in equilibrium with the atmosphere before being perturbed, and vertical exchange is neglected. Different horizontal dispersion coefficients are explored. This case is run by calling the script: `>> case_02p1`

The script initializes a structure of input parameters for the base case and for a few sensitivity runs exploring the effect of changes in horizontal dispersion. The time stepping ocurrs by calling the function `f_case_02p1.m`. This function solves the tracer dispersion via 3 one-way coupled spatio-temporal domains. The first has a horizontal resolution of a few meters, only represents dispersion, and is run for only a few minutes. It serves to disperse the sharp delta-function-shaped injection into a bell curve. The second domain has a horizontal resolution of 100 m and a timestep of about 1 minute. In this domain dispersion and gas exchange are simulated for a few days, before transitioning to the third domain with a horizontal resolution of 1 km and a time step of about 1 hour. In this domain dispersion and gas exchange are simulated for several years.

The script `case_02p1.m` calls the function `pCO2_dictionary.m` to create a look-up table for pCO2 given the predefined temperature and salinity of this case and given inputs of DIC and alkalinity. The time-stepping function `f_case_02p1.m` calls the function `pCO2_fun.m` to retrieve pCO2 values form the look-up table. More specifically, pCO2 is calculated through combination of the lookup table (passed to the function as a dictionary) and bilinear interpolation. This is done because individual calls of `f_csys_alk_DIC.m` are relativley slow. The pCO2 has to be calculated for every gridcell at every timestep, and thus has to be fast for the code to run efficiently. This is accomplished by `pCO2_fun.m` function because dictionary calls and bilinear interpolation are fast. 

### Case 2.2
Case 2.2 examines a spatially resolved surface ocean with a deep ocean reservoir. DIC and alkalinity of the deep reservoir are prescribed and the vertical exchange between surface and deep ocean depends on the vertical gradient of DIC and the vertical exchange coefficient. This case is run by calling the script: `>> case_02p2`

The time stepping ocurrs by calling the function `f_case_02p2.m` in a manner similar to case 2.1. 


## Plotting script
A plotting script `plotting_exmaple.m` is included and can be modified as desired.


## History
The content of the README file is up to date as of July 7, 2025. -KF














