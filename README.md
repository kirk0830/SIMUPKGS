---------------------
Welcome to simu-pkgs!
---------------------

simu-pkgs is a integrated module that contains some basic functions such as:
1) supercell duplication
2) free energy calculation of isolated molecule (based on frequencies result yiedled by other software, but I will add frequency calculation into my package in the future)
3) finite cell x-ray diffraction simulation
4) q6-order parameter analysis


To run this package, you only need to open main.py and run it. Note: numpy, matplotlib package are needed, please install them in advance:\r\n
1) pip install numpy
2) pip install matplotlib

There is also another way to use this package, if you have understood how this package works, you can run all subroutine of it, for instance, pp_smooth.py provides easy-to-use smmoth function to smooth any 1 or 2 dimensional graphs, so you can directy import this module in you OWN program.

As this package is just I write in my spare time (yes this is only one of my hobbies), update may be somehow slow, so if you are not satisfied with it you can also have contribution in this package, it is welcome, certainly!

Footnote: calculation of distance matrix has just been implemented, so the next function that will come soon is dynamics simulation based on classical mechanics. I also have plan to write static calculation based on quantum mechanics, but maybe large amount of time is required.
