# CO2CHAIN
Matlab code implementing the CO2CHAIN method for channel head extraction in high resolution DEMs, developed in Lurin et al, 2023.

This repository contains the main CO2CHAIN function (CO2CHAIN.m) as well as four functions that are called by CO2CHAIN.m

CO2CHAIN requires the MatlabTopoToolBox by Wolfgang Schwanghart (Schwanghart, W., & Kuhn, N. J. (2010). TopoToolbox: A set of Matlab functions for topographic analysis. Environmental Modelling & Software, 25(6), 770-781., https://github.com/wschwanghart/topotoolbox), and thus your Matlab version must fulfill all of the TopoToolBox requirements. 

As for now, CO2CHAIN works on 3 m and 10 m grid resolution DEM, therefore the input of the CO2CHAIN should be in one of these resolution. If you have 1 m data,we recommend dowsampling them first using the 'resample' function. 

THe manual for the 5 scripts are included in the scripts. 
