import numpy as np
#from time import sleep
# distance from point to set, input coordinates of set and point, return a list of distances in the same order of input.
import periodic_boundary as pbc

def dist(
        atomlist, 
        pointxyz, 
        point_in_coord = True,
        pbcFlag = False, 
        pbc_a = 0, 
        pbc_b = 0, 
        pbc_c = 0,
        pbc_alpha = 0,
        pbc_beta = 0,
        pbc_gamma = 0):

    """
    distance calculation function, pbc correction is optional, which is designed to be compatible with either isolated or periodic
    systems.\n
    \n
    # Input description\n
    atomlist: atoms collection that yielded from function "readcoords", a correct format is required.\n
    pointxyz: a one-dimensional list containing element symbol (optional), x, y, z.\n
    point_in_coord and other parameters: see options below.\n
    # Functions\n
    1. calculate distance from one atom to atoms for isolated system, set point_in_coord as True and pbcFlag as False\n
    2. calculate distance from atoms to a macroscopic point, useful when simulate XRD, lights need to be projected on screen pixel.\n
    3. calculate distance from one atom to atoms for periodic system, set point_in_coord as True and pbcFlag as True, 
    remember to input six crystal parameters: a, b, c and alpha, beta, gamma. lengths are in Angstrom and angles are in degree.\n
    WARNING: PBC is NOT compatible with point_in_coord as False!!! (but I will modify this incompatibility in the future)\n
    \n
    # Output\n
    a distance list in the same order of atoms list.\n
    """

    shapeREF = [0, 0, 0]
    #print('DISTANCE| distance calculation is processing, O(N) time complexity algorithm.')
    if np.shape(pointxyz) != np.shape(shapeREF):
        print('DISTANCE| ***error*** [x, y, z] format is expected for point input. QUIT.')
        exit()

    if np.ndim(atomlist) != 2:
        print('DISTANCE| ***warning*** [x, y, z] format is expected for atom-set input.\n'
             +'          this warning may emerge in some special conditions such as calculating a diagonal distance matrix.\n'
             +'          if so, it is safe to ignore it, but check your input otherwise.')
        #sleep(5)
        #exit()

    natom = np.shape(atomlist)[0]
    # ncolumn = np.shape(atomlist)[1], 
    # I directly use: 
    # x = atomlist[iatom][-3]
    # y = atomlist[iatom][-2]
    # z = atomlist[iatom][-1]

    if point_in_coord:
        rescale = 1
    else:
        rescale = 1E-10

    dist_set2point = []

    if pbcFlag:
        cellOper = pbc.cell_operator(
                                     a=pbc_a, 
                                     b=pbc_b, 
                                     c=pbc_c, 
                                     alpha=pbc_alpha, 
                                     beta=pbc_beta, 
                                     gamma=pbc_gamma)

    for iatom in range(natom):

        if pbcFlag:
            disp_r = [
                atomlist[iatom][-3]*rescale-pointxyz[-3],
                atomlist[iatom][-2]*rescale-pointxyz[-2],
                atomlist[iatom][-1]*rescale-pointxyz[-1]
            ]
            disp_r = pbc.vectorPbcCorr(cellOper, disp_r)
            dist_r = np.linalg.norm(disp_r)
        else:
            dist_x = np.sqrt((atomlist[iatom][-3]*rescale-pointxyz[-3])**2)
            dist_y = np.sqrt((atomlist[iatom][-2]*rescale-pointxyz[-2])**2)
            dist_z = np.sqrt((atomlist[iatom][-1]*rescale-pointxyz[-1])**2)
            dist_r = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)

        dist_set2point.append(dist_r)

    return dist_set2point

    # return description:
    # [
    # dist1, dist2, dist3, ...
    # ]
    # if point_in_coord set to True, output will in Angstrom unit, if set to false, output will in SI unit
    # also SI unit is required for pointxyz. Note that this is only needed in XRD projection.