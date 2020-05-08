## Wind turbulence generation using text-based input files

For a more robust and time-efficient Matlab implementation, see https://se.mathworks.com/matlabcentral/fileexchange/68632-wind-field-simulation-the-fast-version.

[![View Wind field simulation (the user-friendly version) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/50041-wind-field-simulation-the-user-friendly-version)
[![DOI](https://zenodo.org/badge/260773326.svg)](https://zenodo.org/badge/latestdoi/260773326)

## Summary

A method to simulate spatially correlated turbulent wind histories is implemented following [1,2]. 
Two possible vertical wind profiles and two possible wind spectra are implemented. The user is free to implement new ones. The wind co-coherence is a simple exponential decay as done by Davenport [3]. If the wind field is simulated in a grid, the function windSim.m should be used (cf. Examples 1 and 2). For a more complex geometry, such as a radial grid, the function windSim.m has an optional parameter to include two inputs (cf. Example3.mlx): The first one contains the wind properties, and the second one contains the coordinates of the nodes where wind histories are simulated (cf. Example 3).

## Content

The folder windSim.zip contains:
-	1 input file INPUT.txt for Example1.m
-	1 input file INPUT_MAST.txt for Example2.m
-	2 input files windData.txt and circle.txt for Example3.m
-	The function windSim.m
-	3 examples files Example1.m, Example2.m, Example3.m
-	The function coherence.m that computes the co-coherence.
Notes:
-	Simulating the wind field in a high number of points with a high sampling frequency may take a lot of time. 
-	This code aims to be highly customizable 
-	A faster version of the present submission has been used to simulate the turbulent wind load on a floating suspension bridge [4]. 

## References

[1] Shinozuka, M., Monte Carlo solution of structural dynamics, Computers and Structures, Vol. 2, 1972, pp. 855 – 874 

[2] Deodatis, G., Simulation of ergodic multivariate stochastic processes, Journal of Engineering Mechanics, ASCE, Vol. 122 No. 8, 1996, pp. 778 – 787. 

[3] Davenport, A. G. (1961), The spectrum of horizontal gustiness near the ground in high winds. Q.J.R. Meteorol. Soc., 87: 194–211 

[4] Wang, J., Cheynet, E., Snæbjörnsson, J. Þ., & Jakobsen, J. B. (2018). Coupled aerodynamic and hydrodynamic response of a long span bridge suspended from floating towers. Journal of Wind Engineering and Industrial Aerodynamics, 177, 19-31.
