# 2S_ESS Code and Documentation

The 2S-ESS code was originally developed by Robert Spurr (RT Solutions, Inc.) and optimized by Vijay Natraj (Jet Propulsion 
Laboratory, California Institute of Technology).

The code is in subdirectory sourcecode. Documentation of the 2S and ESS codes are in subdirectory docs. Refer also to the paper 
to be published in JQSRT (citation to be added after paper is accepted and has a doi).

Driver programs for a solar forcing scenario (no thermal emission) are in subdirectory TEST. There is one program each for the 
regular (original) code and for the optimzed version. The results are saved in results_saved.dat and results_opt_saved.dat, 
respectively.

The test case has the following main parameters:

114 layers, all with the same optical properties (i.e. identical layers)
total optical depth = 1.14 (0.01 per layer)
single scattering albedo = 0.5
asymmetry parameter = 0.9 (Henyey-Greenstein phase function assumed)
surface albedo = 0.3
SZA = 30 degrees
VZA = 0 degrees (nadir)
AZM = 0 degrees (principal plane)
layer height = 1 km (TOA altitude = 114 km)
"enhanced" spherical option used

For other use cases, the user can contact the author at vijay.natraj@jpl.nasa.gov.
 
