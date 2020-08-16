import numpy as np

def centercoords(coordsinp, masspower = False, totmass = 1, masslist = []):

    """
    This function is to center all atoms. If you want to center them based on their mass, set masspower as True,
    and total mass in relative atomic mass is required, also, mass list is required.\n
    coordsinp: 2 dimensional list, format like standard xyz file\n
    masspower: bool\n
    totmass: float\n
    masslist: 1 dimensional float list, order should be in accord with coordsinp
    """

    print('COORD| i hope you have convert coordinates to float type, if not, no one knows what will happen...')

    natom = np.shape(coordsinp)[0]
    if masspower:
        if np.shape(masslist)[0] != natom:
            print('COORD| ***error*** number of atoms in atomlist is not consistent with masslist. quit.')
            exit()
        else:
            print('COORD| mass-powered atom centering is activated, useful when calculate rotation interia.')

    x_center = 0
    y_center = 0
    z_center = 0

    coordsout = coordsinp

    for iatom in range(natom):

        mass_power = 1

        if masspower:
            mass_power = masslist[iatom]

        x_center += coordsinp[iatom][-3]/natom * (mass_power/totmass)
        y_center += coordsinp[iatom][-2]/natom * (mass_power/totmass)
        z_center += coordsinp[iatom][-1]/natom * (mass_power/totmass)
    
    for iatom in range(natom):

        coordsout[iatom][-3] -= x_center
        coordsout[iatom][-2] -= y_center
        coordsout[iatom][-1] -= z_center

    return coordsout

# one-atom operations

def rotate_x(coord, angle):

    # rotate coordinates around x-axis, which means, rotate coordinates in yz plane.
    # derivation:
    # original vector: rcosA, rsinA
    # rotated vector: rcos(A+B), rsin(A+B) = rcosAcosB - rsinAsinB, rsinAcosB + rcosAsinB
    # re-write in matrix form:
    # 1   0     0  | x       x                       x
    # 0 cosB -sinB | rcosA = rcosAcosB - rsinAsinB = rcos(A+B)
    # 0 sinB  cosB | rsinA   rsinAcosB + rcosAsinB   rsin(A+B)

    angle = angle/180*np.pi

    operator = [
        [1, 0, 0],
        [0, np.cos(angle), -np.sin(angle)],
        [0, np.sin(angle), np.cos(angle)]
    ]

    return np.matmul(coord, operator)

def rotate_y(coord, angle):

    # rotate coordinates around y-axis, which means, rotate coordinates in xz plane.

    # re-write in matrix form:
    # cosB   0   -sinB |
    #  0     1      0  |
    # sinB   0    cosB |

    angle = angle/180*np.pi

    operator = [
        [np.cos(angle), 0, -np.sin(angle)],
        [0, 1, 0],
        [np.sin(angle), 0, np.cos(angle)]
    ]

    return np.matmul(coord, operator)

def inversion(coord, mode = 'xyz'):

    if mode == 'xyz':
        coord[-3] = -coord[-3] # x
        coord[-2] = -coord[-2] # y
        coord[-1] = -coord[-1] # z
    elif mode == 'x':
        # same as reflection towards yz-plane.
        coord[-3] = -coord[-3]
    elif mode == 'y':
        # same as reflection towards xz-plane.
        coord[-2] = -coord[-2]
    elif mode == 'z':
        # same as reflection towards xy-plane.
        coord[-1] = -coord[-1]
    elif mode == 'xy':
        coord[-3] = -coord[-3] # x
        coord[-2] = -coord[-2] # y
    elif mode == 'yz':
        coord[-2] = -coord[-2] # y
        coord[-1] = -coord[-1] # z
    elif mode == 'xz':
        coord[-3] = -coord[-3] # x
        coord[-1] = -coord[-1] # z
    else:
        exit()
    return coord

def mirror_xy(coord):

    return inversion(coord, mode='z')

def mirror_xz(coord):

    return inversion(coord, mode='y')

def mirror_yz(coord):

    return inversion(coord, mode='x')

# poly-atom operations
# visit the 1st atom in list1, find atom_i in list2, residue: list2[0:i-1], list2[i+1:], add together name as list2
# visit the 2nd atom in list1, find atom_j in list2...
# ...

def findinlist(atom, atomlist):

    # do not call this subroutine from external programs except you know what you are doing.

    natom = np.shape(atomlist)[0]
    nTrue = 0
    index = -1

    for iatom in range(natom):
        if np.array_equal(atom, atomlist[iatom][:]):
            nTrue += 1
            index = iatom

    if nTrue > 1:
        print('COORD| ***error*** atoms overlap! quit.')
        exit()
    else:
        return index

def isCoordSame(coords1, coords2, acc = 6):

    natom = np.shape(coords1)[0]

    if np.shape(coords2)[0] != natom:
        print('COORD| ***error*** atom number is not in consistency.')
        return False
    else:
        print('COORD| atom number check passed. :) Good luck.')
    # pre-processing of coordinates
    for iatom in range(natom):

        coords1[iatom][-1] = round(coords1[iatom][-1], acc)
        coords1[iatom][-2] = round(coords1[iatom][-2], acc)
        coords1[iatom][-3] = round(coords1[iatom][-3], acc)
        coords2[iatom][-1] = round(coords2[iatom][-1], acc)
        coords2[iatom][-2] = round(coords2[iatom][-2], acc)
        coords2[iatom][-3] = round(coords2[iatom][-3], acc)
    
    print('COORD| convert two atomlists to the same precision level: '+str(10**-acc))
    list2find = coords2
    index = -2

    identicalFlag = False
    for iatom in range(natom):

        if index == -2:
            index = findinlist(coords1[iatom][:], list2find)
                
        if iatom > 0:
            if index >= 0:
                list2find = list2find[0:index]+list2find[index+1:]
                index = findinlist(coords1[iatom][:], list2find)
                if iatom == natom-1 and index >= 0:
                    identicalFlag = True
            else:
                break

    return identicalFlag

def cell_duplication(coords, a, b, c, alpha = 90, beta = 90, gamma = 90, nx=1, ny=1, nz=1):
    print('SUPERCELL| supercell function is activated. I hope you have type-in parameters\n'
         +'           that all required. alpha, beta, gamma will use 90 deg as default if \n'
         +'           not explicitly defined.'
    )
    if np.ndim(coords) != 2:
        print('SUPERCELL| ***ERROR*** may be you forget to give coordinates list. Program is terminated.')
        exit()

    if a<=0 or b<=0 or c<=0 or alpha<=0 or beta<=0 or gamma<=0 or alpha>=180 or beta>=180 or gamma>=180:
        print('SUPERCELL| ***ERROR*** invalid parameter input! Program is terminated.')
        exit()
    
    a = float(a)
    b = float(b)
    c = float(c)
    alpha_r = float(alpha/180) * np.pi
    beta_r = float(beta/180) * np.pi
    gamma_r = float(gamma/180) * np.pi
    a_vec = [a, 0, 0]
    b_vec = [b*np.cos(alpha_r), b*np.sin(alpha_r), 0]
    c1 = c*np.cos(gamma_r)
    c2 = c*(np.cos(beta_r)-np.cos(alpha_r)*np.cos(gamma_r))/np.sin(alpha_r)
    c3 = np.sqrt(c**2 - c1**2 - c2**2)
    c_vec = [c1, c2, c3]

    print('SUPERCELL| cell parameters information:\n'
         +'           a = '+str(a)+' Angstrom | alpha = '+str(alpha)+' deg\n'
         +'           b = '+str(b)+' Angstrom | beta = '+str(beta)+' deg\n'
         +'           c = '+str(c)+' Angstrom | gamma = '+str(gamma)+' deg\n'
         +'           cell(1) = '+str(a_vec)+'\n'
         +'           cell(2) = '+str(b_vec)+'\n'
         +'           cell(3) = '+str(c_vec)+'\n'
    )
    
    natom = np.shape(coords)[0]
    print('SUPERCELL| original cell information: '+str(natom)+' atoms in total.')
    ncol = np.shape(coords)[1]

    if nx != int(nx) or ny != int(ny) or nz != int(nz) or nx <= 0 or ny <= 0 or nz <= 0:
        print('SUPERCELL| invalid duplicate number, they should have been integars. quit.')
        exit()
    
    coords_readin = coords
    coords_add = []
    if ncol == 3:
        print('SUPERCELL| there seems no element information in coordinates input, single-\n'
             +'           element treatment will be employed...'
        )
        for iatom in range(natom):
            for inx in range(nx):
                for iny in range(ny):
                    for inz in range(nz):
                        atom_dup = coords[iatom][:]
                        atom_dup[0] += (inx*a_vec[0] + iny*b_vec[0] + inz*c_vec[0])
                        atom_dup[1] += (inx*a_vec[1] + iny*b_vec[1] + inz*c_vec[1])
                        atom_dup[2] += (inx*a_vec[2] + iny*b_vec[2] + inz*c_vec[2])
                        if inx or iny or inz:
                            coords_add.append(atom_dup)
    elif ncol == 4:
        for iatom in range(natom):
            for inx in range(nx):
                for iny in range(ny):
                    for inz in range(nz):
                        # atom_dup = coords[iatom] # do not execute this line!
                        atom_dup = coords[iatom][:]
                        atom_dup[1] += (inx*a_vec[0] + iny*b_vec[0] + inz*c_vec[0])
                        atom_dup[2] += (inx*a_vec[1] + iny*b_vec[1] + inz*c_vec[1])
                        atom_dup[3] += (inx*a_vec[2] + iny*b_vec[2] + inz*c_vec[2])
                        if inx or iny or inz:
                            coords_add.append(atom_dup)
# developer notes here: it is interesting that if use:
# atom_dup = coords[iatom], atom_dup is not a new object, instead, it is, the constant pointer,
# or say, the re-named coords[iatom]. So latter steps will directly change data in coords[][],
# which is not what we expect.
# instead, if we write atom_dup = coords[iatom][:], a new object will be created, all following
# steps will proceed as usual.
    else:
        print('SUPERCELL| ***ERROR*** wrong coordinates format. Program is terminated.')
        exit()

    return coords_readin+coords_add
