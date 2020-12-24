---------------------
Welcome to SIMUPKGS
---------------------

UPDATE at 2020/12/24:
A basic framework of Hartree-Fock method is added into directory. Modules support function-integration, more specificly, one-electron integrator, two-, three and four will be coded and upload. For I am busy recently so this will be a long process.
Also there seem some module that I just code them in several month ago but failed to upload them.

UPDATE at 2020/6:
SIMUPKGS is a integrated module that contains some basic functions such as:
1) supercell duplication
2) free energy calculation of isolated molecule (based on frequencies result yiedled by other software, but I will add frequency calculation into my package in the future)
3) finite cell x-ray diffraction simulation
4) q4-order parameter analysis

To run this package, you only need to open main.py and run it. Note: numpy, matplotlib package are needed, please install them in advance:\r\n
1) pip install numpy (for most modules)
2) pip install matplotlib (for XRD related modules)
3) pip install pandas (for absToolkit)

There is also another way to use this package, if you have understood how this package works, you can run all subroutine of it, for instance, pp_smooth.py provides easy-to-use smmoth function to smooth any 1 or 2 dimensional graphs, so you can directy import this module in you OWN program.

As this package is just I write in my spare time (yes this is only one of my hobbies), update may be somehow slow, so if you are not satisfied with it you can also have contribution in this package, it is welcome, certainly!

Footnote: calculation of distance matrix has just been implemented, so the next function that will come soon is dynamics simulation based on classical mechanics. I also have plan to write static calculation based on quantum mechanics, but maybe large amount of time is required.
