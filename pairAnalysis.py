from distance import dist
import numpy as np

# this function is to generate distance diagnoal matrix
def distTriangleMat(
                    coordslist, 
                    full_matrix = True,
                    pbc = False,
                    cell_a = 0,
                    cell_b = 0,
                    cell_c = 0,
                    cell_alpha = 90,
                    cell_beta = 90,
                    cell_gamma = 90
                    ):
    
    """
    function to calculate distance between ALL atoms, return a OD matrix, 2 dimensional.\n
    # input requirement\n
    coordlist: atomic coordinates, standardized and widely used in this package\n
    full_matrix: [bool] output a fully diagonal matrix or not. If set as false, half of elements in output matrix will be zero cuz have not been calculated\n
    pbc: [bool] if use periodic boundary condition correction. Useful in periodic system, set to False if called in an isolated system calculation\n
    other six parameters: MUST be set reasonable values if set pbc as True, these are cell parameters\n
    """

    natom = np.shape(coordslist)[0]
    # dont know how many columns the coordlist has, either 3 or 4, use -1, -2 and -3 to visit
    distmat = np.zeros((natom, natom))
    print('DISTANCE| distance matrix generated, size: '+str(natom)+'x'+str(natom)+'.')
    for i in range(natom):

        if pbc:
            distlist = dist(
                            atomlist=coordslist[i+1:][:], 
                            pointxyz=coordslist[i][-3:],
                            point_in_coord=True,
                            pbcFlag=True,
                            pbc_a=cell_a,
                            pbc_b=cell_b,
                            pbc_c=cell_c,
                            pbc_alpha=cell_alpha,
                            pbc_beta=cell_beta,
                            pbc_gamma=cell_gamma
                            )
        
        else:
            distlist = dist(
                            atomlist=coordslist[i+1:][:], 
                            pointxyz=coordslist[i][-3:]
                            )

        print('DISTANCE| calculate distance list from atom '+str(i)+' to other atoms...')

        for j in range(natom-i-1):

            distmat[i][i+j+1] = distlist[j]
            if full_matrix:
                distmat[i+j+1][i] = distmat[i][i+j+1]

    return distmat

# for a given cutoff value, return the dict that contains nearest neighboring particles' indices
def nearestNeighbors(distmat, cutoff):

    print('PAIR| count nearest neighboring particles... current cutoff value in Angstrom: '+str(cutoff))
    natom = np.shape(distmat)[0]

    # initialization of neighborDict:
    neighbordict = {}
    for i in range(natom):

        neighbordict[i] = []

    for i in range(natom):
        for j in range(natom-i-1):

            if distmat[i][i+j+1] < cutoff:

                # structured data
                neighbordict[i].append((i+j+1, distmat[i][i+j+1]))
                neighbordict[i+j+1].append((i, distmat[i][i+j+1]))

    print('PAIR| nearestNeighbors information has been saved in dict. Note that data in parenthesis has\n'
         +'      no difference with that in bracket. It is just for convienence of later pair-calculation.')
    return neighbordict

def triangleExist(pointA, pointB, pointC, collection):

    aTriangle = [pointA, pointB, pointC]
    nComb = np.shape(collection)[0]

    for iComb in range(nComb):

        if np.array_equal(aTriangle, collection[iComb]):

            return True
            break

def q6_order_parameter(coords, cutoff):

    """
    q6 order parameter calculation\n
    # Definition\n
    q6 of atom i is denoted as q6_i,\n
    q6_i = 1 - 3/8 sum((cos(theta)+1/3)^2, theta), theta)\n
    theta is the angle of j-i-k, j and k are atoms neighboring to atom i, the summation above goes over all possible atoms j and k.
    # Input requirement\n
    coords: atomic coordinates, 2 dimensional list, as used over all functions in this package\n
    cutoff: cutoff value given manually to calculate only neigboring atoms within a certain distance, unit is Angstrom\n
    # Output description\n
    a dictionary that contains all possible atoms' q6
    """
    print('PAIR| q6 order parameter analysis is activated, a series of calculations will be carried out.')
    distMatrix = distTriangleMat(coordslist=coords)
    neighors = nearestNeighbors(distmat=distMatrix, cutoff=cutoff)

    description = [('index',int),('distance',float)]
    natom = len(neighors)

    q6_dict = {}
    for iatom in range(natom):

        connectAtoms = neighors[iatom]
        q6_thisAtom = 0

        # format: [(1, dist1), (2, dist2), ...]
        if len(connectAtoms) >= 4:
            print('PAIR| q6-calculation is performed on atom '+str(iatom))
            data_to_sort = np.array(connectAtoms, dtype=description)
            data_sorted = np.sort(data_to_sort, order='distance')
            points_to_calc = [data_sorted[0][0], data_sorted[1][0], data_sorted[2][0], data_sorted[3][0]]
            print('PAIR| q6: atom '+str(iatom)+' -> atoms '+str(points_to_calc))
            # start
            for atom_i in range(3):
                for disp_i in range(1,4-atom_i):
                    # vector_i
                    vector_i = [
                        coords[points_to_calc[atom_i]][-3]-coords[iatom][-3],
                        coords[points_to_calc[atom_i]][-2]-coords[iatom][-2],
                        coords[points_to_calc[atom_i]][-1]-coords[iatom][-1]
                        ]
                    vector_j = [
                        coords[points_to_calc[atom_i+disp_i]][-3]-coords[iatom][-3],
                        coords[points_to_calc[atom_i+disp_i]][-2]-coords[iatom][-2],
                        coords[points_to_calc[atom_i+disp_i]][-1]-coords[iatom][-1]  
                    ]
                    cos_ij = np.dot(vector_i, vector_j)/np.linalg.norm(vector_i)/np.linalg.norm(vector_j)
                    q6_thisAtom -= 3/8*(cos_ij+1/3)**2

            q6_thisAtom += 1
            q6_dict[iatom] = q6_thisAtom

        else:
            print('PAIR| for atom '+str(iatom)+' has neighboring atoms less than 4, q6-calculation skips.')

    return q6_dict
