# simu_pkg
simu_pkg is program that i write in my spare time mainly focus on simulating some simple things such as XRD and molecular dynamics and so on if possible in the future.

how to run
----------
1. to run complete package, execute mainProgram.py. All execution setting parameters should be provided in advance in an unified format, for instance, please see input.inp provided.
also, atomic coordinate file is required. Format shoule be standard xyz file.
2. note that actually all subroutine can run seperately and can be flexible if you want to employ one certain function in your own program or, just some data pre/post-processing, it is okay.

mainProgram.py
--------------
integrated entry of simu_pkg, usage can be found in input file -> input.inp. Note that the filename of input file can be arbitrary. Just type-in filename when this program tells you to do.
presently functions supported are:
1. supercell duplication, atomic coordinates centering... that all basic preprocessing functions needed for X-ray diffraction simulation
2. X-ray diffraction simulation
3. explicit use of supercell duplication

centerCoords.py
---------------
center atomic coordinates, input format is 2 dimensional list, REMEMBER to check datatype of coordinates. FLOAT type is required.

for short:
input requirement: 2 dimensional list, coordinates must be float type
output description: 2 dimensional list

coordsRead.py
-------------
read standard formatted xyz files and, (optional, default as false) convert coordinates from string to float.

for short:
input requirement: standard xyz file, and filename
output description: 2 dimensional list

coordsWrite.py
--------------
write standard formatted xyz files, a 2-dimensional list that including atomic coordinates is compulsory. It also supports input element information from external standard input (i.e., from keyboard) if an n*3 list is detected in input parameter list.

for short:
input requirement: 2 dimensional list, coordinates must be float type, also a filename.
output description: no print-out or return, just save a file.

diffraction.py
--------------
shot a beam from one source to one acceptor, phase information is preserved.
WARNING: in final XRD pattern visualization both 1d and 2d, I forget to take absolute value of datapoints, please, take care, however preservation of phase information may be useful and I will add another option that controls preservation of phase information.

for short:
input requirement: wavelength, beam propagating distance, (not implemented yet) angular resoluted amplitude parameters
output description: intensity

distance.py
-----------
calculate point to set distance, atomic coordinates list and a point xyz is required.

for short:
input requirement: atomic coordinates 2 dimensional list, point coordinate 1 dimensional list
output description: distance list 1 dimensional

pxrd.py
-------
post-processing of 2d-xrd pattern, add intensities of all points with same radius.

for short:
input requirement: 2d xrd pattern dataset 2 dimensional list
output description: 1d pxrd pattern 1 dimensional list

readinput.py
------------
read inputfile and save keywords and parameters in dictionary.

for short:
input requirement: unified input file, filename
output description: a dictionary containing all possible keywords and values

smooth.py
---------
smooth datasets, 1d and 2d datasets are supported. For gaussian smooth, a cutoff method for accelerating smooth that avoids O(N^(2n)) time complexity is supported.

for short:
input requirement: 1d or 2d datasets, sigma value and other optional parameters.
output description: same-sized datasets as input

supercell.py
------------
duplicate atomic coordinates in certain given times, at most six crystal parameters are required. symmetry check routine will be added in the future.

for short:
input requirement: atomic coordinates list 2 dimensional list, at most six crystal parameters, number of replicas in 3 dimensions
output description: atomic coordinates list 2 dimensional list

xrd.py
------
main xrd subroutine. directly use is not recommended for negligible external application possibility.

xrd_madescr.py
--------------
initialize screen that x ray projects to. directly use is not recommended for negligible external application possibility.

other files
-----------
input.inp: parameter file that contains keywords supported by mainProgram.py presently.
sketch.xyz: example atomic coordinates file, in standard xyz format.

