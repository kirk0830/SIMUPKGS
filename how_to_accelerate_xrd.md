yesterday i was finding ways to accelerate XRD simulation. Because for now, if denote length (also width, for there is no reason to consider a non-square shaped screen) projected screen where diffracted X ray is accepted, as N, number of atoms as M, therefore distance calculation has O(M*N^2) time complexity, diffraction has O(M*N^2) time complexity too. Note that N will be significantly larger than M, commonly, in actually use, there will be at least millions of steps, which seems horrible and, performance is quite inferior to that implemented in some matural simulation packages such as GROMACS, etc.

Actually, crystal diffraction is just inference of plenty of single-atom x ray scattering. And if distance from atom to screen in perpendicular direction is unchanged, 
ONLY ONE LINE OF DATA IS NEEDED TO CALCULATED!
and all intensity information can be obtained simply by adding some terms in "the line of data".

this is, because:
1. for one atom, single-atom scattering patter can be obtained by rotating the "line" for 360 degree, i.e., single-atom scattering is isotropic.
2. for atoms with the same distance from itself to screen in perpendicular direction, anisotropic type pattern can be obtained by copy, move and summation isotropic pattern obtained in last step.
3. for atom with different distance to screen, some scattering may need to re-calculate because distance to certain pixel may be not strictly equal to any one of other atom.
therefore for a single element crystal, one may only need to calculate one atom for each atomic layer.

However, I also meet problems when I want to realize this method. One main problem is the inconsistency between resolution of atomic coordinates and that of screen pixels. For atoms, their distance is in angstrom scale, 1E-10 meters, for pixels, 1E-5 is already VERY expensive for simulation, however, can not distinguish atom at all, yet.

therefore there is a THING will happen that although it is physically equivalent that two-atom diffraction is just the result of overlapping of two single-atom scattering, on screen, it can not be recognized and datas will be saved in the same pixel.
Now you may want to ask why this problem did not emerge in present O(M*N^2) algorithm, answer is, bacause wavelength is in angstrom scale, intensity is quite sensitive to the change of distance. However if we store this information in low resolution (the screen), nm-resoluted information will be lost.

so I will find other ways or, I will decrease the perpendicular distance from crystal to screen, because in time complexity view, this method seems quite attracting.
