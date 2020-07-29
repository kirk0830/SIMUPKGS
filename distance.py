import numpy as np

# distance from point to set, input coordinates of set and point, return a list of distances in the same order of input.

def dist(atomlist, pointxyz):

    shapeREF = [0, 0, 0]
    #print('DISTANCE| distance calculation is processing, O(N) time complexity algorithm.')
    if np.shape(pointxyz) != np.shape(shapeREF):
        print('DISTANCE| ***error*** [x, y, z] format is expected for point input. QUIT.')
        exit()

    if np.ndim(atomlist) != 2:
        print('DISTANCE| ***error*** [x, y, z] format is expected for atom-set input. QUIT.')
        exit()

    natom = np.shape(atomlist)[0]
    # ncolumn = np.shape(atomlist)[1], 
    # I directly use: 
    # x = atomlist[iatom][-3]
    # y = atomlist[iatom][-2]
    # z = atomlist[iatom][-1]

    dist_set2point = [ 
        np.sqrt(
            (atomlist[iatom][-3]*1E-10-pointxyz[0])**2
            +(atomlist[iatom][-2]*1E-10-pointxyz[1])**2
            +(atomlist[iatom][-1]*1E-10-pointxyz[2])**2
            ) for iatom in range(natom) 
        ]
    return dist_set2point

    # return description:
    # [
    # dist1, dist2, dist3, ...
    # ]