This package provides a framework that calculates the trajectory of a test charge in an electric field that is predetermined
by user input. There are three possible inputs: "null", "single", and "double"; the inputs must be given in the Build file, on
line 31. Also on this line, the dimensions of the axes can be changed to increase/decrease resolution. The three inputs correspond
to no charge distribution, a single Gaussian peak, and a pair of Gaussian peaks.

Contributions:
GitHub author tracking should be representative of the work done, but to be explicit, MD wrote the code for the Gauss-Seidel module,
trajectory plotter Python script, and the build file, and GS wrote the code for the verlet module, and netCDF writer, and we both
contributed to the main file.
