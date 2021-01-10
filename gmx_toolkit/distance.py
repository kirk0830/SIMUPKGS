import numpy as np
#from time import sleep
# distance from point to set, input coordinates of set and point, return a list of distances in the same order of input.
import periodic_boundary as pbc

def dist(
        atomlist, 
        pointxyz, 
        point_in_unit = 'Angstrom',
        set_in_unit = 'Angstrom',
        dist_in_unit = 'Angstrom',
        pbcFlag = False, 
        pbc_a = 0, 
        pbc_b = 0, 
        pbc_c = 0,
        pbc_alpha = 0,
        pbc_beta = 0,
        pbc_gamma = 0,
        append_ele = False,
        ):

    """
    distance calculation function, pbc correction is optional, which is designed to be compatible with either isolated or periodic
    systems.\n
    \n
    # Input description\n
    atomlist: atoms collection that yielded from function "readcoords", a correct format is required.\n
    pointxyz: a one-dimensional list containing element symbol (optional), x, y, z, if not specify .\n
    point_in_unit: unit of referring point\n
    set_in_unit: units of coordinates\n
    dist_in_unit: units of output distances\n
    PBC parameters: see options below.\n
    append_ele: append pair elemental information at the end of array. Input this parameter by directly inputing element symbol\n
    i.e.: append_ele = 'C'\n
    # Functions\n
    1. calculate distance from one atom to atoms for isolated system, set pbcFlag as False\n
    2. calculate distance from atoms to a macroscopic point, useful when simulate XRD, lights need to be projected on screen pixel.\n
    3. calculate distance from one atom to atoms for periodic system, set pbcFlag as True, 
    remember to input six crystal parameters: a, b, c and alpha, beta, gamma. lengths are in Angstrom and angles are in degree.\n
    WARNING: PBC is NOT compatible with set_in_unit as False!!! (but I will modify this incompatibility in the future)\n
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

    if set_in_unit == 'Angstrom':
        set_rescale = 1
    elif set_in_unit == 'nm':
        set_rescale = 10
    elif set_in_unit == 'SI':
        set_rescale = 1E10
    
    if point_in_unit == 'Angstrom':
        point_rescale = 1
    elif point_in_unit == 'nm':
        point_rescale = 10
    elif point_in_unit == 'SI':
        point_rescale = 1E10
    
    if dist_in_unit == 'Angstrom':
        dist_rescale = 1
    elif dist_in_unit == 'nm':
        dist_rescale = 1E-1
    elif dist_in_unit == 'SI':
        dist_rescale = 1E-10

    dist_set2point = []
    elem_pair = []
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
                atomlist[iatom][-3]*set_rescale-pointxyz[-3]*point_rescale,
                atomlist[iatom][-2]*set_rescale-pointxyz[-2]*point_rescale,
                atomlist[iatom][-1]*set_rescale-pointxyz[-1]*point_rescale
            ]
            disp_r = pbc.vectorPbcCorr(cellOper, disp_r)
            dist_r = np.linalg.norm(disp_r)*dist_rescale
        else:
            dist_x = np.sqrt((atomlist[iatom][-3]*set_rescale-pointxyz[-3]*point_rescale)**2)
            dist_y = np.sqrt((atomlist[iatom][-2]*set_rescale-pointxyz[-2]*point_rescale)**2)
            dist_z = np.sqrt((atomlist[iatom][-1]*set_rescale-pointxyz[-1]*point_rescale)**2)
            dist_r = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)*dist_rescale

        if append_ele:
            pair = [str(append_ele), str(atomlist[iatom][0])]
            elem_pair.append(pair)
        
        dist_set2point.append(dist_r)

    if append_ele:
        return [dist_set2point, elem_pair]
    else:
        return dist_set2point

    # return description:
    # [
    # dist1, dist2, dist3, ...
    # ]
    # if set_in_unit set to True, output will in Angstrom unit, if set to false, output will in SI unit
    # also SI unit is required for pointxyz. Note that this is only needed in XRD projection.