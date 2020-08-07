import numpy as np

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
