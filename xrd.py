import numpy as np

from xrd_madescr import xrd_madescr
from distance import dist
from diffraction import diffraction

def xrd(
    wavelength,
    maxangle,
    dist2scr,
    scrresol,
    if_angle_resol,
    atomcoords
):
    natom = np.shape(atomcoords)[0]
    print('XRD| main subroutine activated, totally '+str(natom)+' in simulation box. An accelerated version will come soon.')
    print('XRD| periodic boundary condition method is not implemented, for accuracy, please use SUPERCELL instead.')
    wavelength_ = wavelength * 1E-9
    print('XRD| units conversion: wavelength value input: '+str(wavelength)+' -> '+str(wavelength_)+' (m)')

    scr2d = xrd_madescr(angle_max=maxangle, dist=dist2scr, resol=scrresol, give2d=True)
    scrcoords = xrd_madescr(angle_max=maxangle, dist=dist2scr, resol=scrresol)

    scr_width = np.shape(scr2d)[0]
    print('DISTANCE| i will be iteratively called...')
    print('DIFFRACTION| i will be iteratively called too...')
    print('DIFFRACTION| time complexity estimation: totally '+str(scr_width**2 * scr_width**2 * natom)+' steps.')

    #exit()
    for ix in range(scr_width):
        for iy in range(scr_width):
            pixelcoord = scrcoords[ix+iy*scr_width][:]
            dist2atoms = dist(atomlist=atomcoords, pointxyz=pixelcoord)
            for idist in dist2atoms:
                scr2d[ix][iy] += diffraction(wavelength=wavelength_, 
                                             dist=idist,
                                             anglemode=if_angle_resol,
                                             )
    
    return scr2d
