# Quasi-periodic-orbit-generator

*This repository contains the code relative to the article Optimal Escape from Sun-Earth and Earth-Moon L2 with Electric Propulsion (doi: https://doi.org/10.3390/aerospace9040186) and to my PhD thesis "Low-Thrust Optimal Escape Trajectories from Lagrangian Points and Quasi-Periodic Orbits in a High-Fidelity Model" (https://iris.polito.it/handle/11583/2976595). Please cite as 

"Mascolo, Luigi, and Lorenzo Casalino. 2022. "Optimal Escape from Sun-Earth and Earth-Moon L2 with Electric Propulsion" Aerospace 9, no. 4: 186. https://doi.org/10.3390/aerospace9040186" 
or as
"Low-Thrust Optimal Escape Trajectories from Lagrangian Points and Quasi-Periodic Orbits in a High-Fidelity Model / Mascolo, Luigi. - (2023 Feb 03), pp. 1-217." 
when referring to this repository.*

This code comes from a regroup of my efforts to automate the search of Quasi-Periodic Orbits in higher fidelity models during my PhD. 
It is based on different codes and inherited some expertise from previous efforts at the Politecnico di Torino (in FORTRAN).
It mixes a great part of Matlab computations and a single executable file compiled in FORTRAN (to protect its intellectual property). 

Due to size-limit issues, JPL's DE430 ephemerides (implemented for the high-fidelity model and used for the state of celestial bodies throughout the analysis) have to be downloaded separately at the following link:
https://www.dropbox.com/sh/15vc4ddqwddgrxp/AAAx_WBNzVPxvTYjTen83wh8a?dl=0
These two files have to be in the same folder containing all the Matlab and FORTRAN files. 

Ciao, have fun!

I'm very open to suggestions: please reach out to me at luigi.mascolo@polito.it :) 
