ETAS and Inversion for Beginners.
The codes only contain basics for the ETAS simulation and MLE parameter inversion.

Download fminsearchbnd at the following address and put the 'fminsearchbnd' in the same folder.
https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon

The ETAS simulator is a simplified version of AFTSimulator.m:  

RUN_ETAS_Simulation.m -> generate catalog

RUN_MLE_Inversion -> Invert parameters from the generated catalog

For ETAS simulation, refer to (AFTSimulator.m)
Felzer, K. R., T. W. Becker, R. E. Abercrombie, G. Ekstrom, and J. R.
Rice, Triggering of the 1999 Mw 7.1 Hector Mine earthquake by aftershocks
of the 1992 Mw 7.3 Landers earthquake, J. Geophys. Res., 107, 2190,
doi:10.1029/2001JB000911, 2002.

The inversion code is created for:
Im, K. and Avouac, J.-P., Cascading Foreshocks, Aftershocks, and Earthquake 
Swarms in a Discrete Fault Network, GJI, 202?
