# Physics

******************
./Elldot/
******************
Code applying Kohn-Sham equations to a circular quantum dot. 
Elliptic and rectangular perturbation can be added to the system.
The code also features a simple predictor-corrector real-time solver.

******************
./Fortran-Localization/
******************
Application of the Foster-Boys localization method on various one-dimensional well systems. 
I have used Numerical Recipes' Jacobi rotation routine to minimize the Foster-Boys objective function. 
After applying the localization routine, the system excitations are calculated using linear response, more exactly the Casida Equations.
Output is written as the normalized transition density matrices and particle-hole maps of the given system, which can then be plotted in gnuplot to visualize the electron excitation density differences.

******************
./Gaussian/
******************
Compilation files for Gaussian09 computations of organic photovoltaic electron excitations and charge-transfer processes.

******************
./Hubbard/
******************
Primarily Fortran 90 files applying the Schroedinger equation on a simple multi-site lattice system with non-collinear magnetism and solving via direct diagonalization of the interacting (exact) and noninteracting electron systems.
The noninteracting system is then iteratively optimized to calculate the Kohn-Sham fields required to reproduce the exact density of the system. 
I utilize the Complex Gradient method outlined in Numerical recipes to optimize the noninteracting solution, which is then solved to yield the exact exchange-correlation magnetic fields and energies of the lattice system.

******************
./NWChem-Codes/
******************
Particle-hole map and localized particle-hole map python generator files. The files are separated by the program used to calculate the excitation energies of the system. The differences between Gaussian and NWChem output files have been accounted for and the PHM generators produce particle-hole maps for each output file case. 
Localized maps have not been completed for Gaussian systems as the rotation matrices generated from the Foster-Boys mimization have not been accurately applied and computed yet.
Research to accurately apply the rotation matrices on the canonical systems is ongoing.
