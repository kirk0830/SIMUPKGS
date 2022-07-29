# this program is compatible with NPT ensemble.
# during simulation, the fractional coordinates are invariant with rescaling of cell dimension
# therefore program will directly make the [0,1]^3 space to grid
# but cell dimension is better to be provided to rationally make grid size
from lmp_readin import lmp2list
import numpy as np

def ndigitsgen(frame, nx, ny, nz):
    
    '''
    Generate a list in 3 dimension, from fractional xyz data, respect to certain level of griding\n
    For example dividing [0,1] along x-axis in 1001 slices, say, 0, 0.001, 0.002, 0.003, ..., 0.999, 1.000\n
    frame (2d list): in format of [['element1', x1, y1, z1], ['element2', x2, y2, z2], ...]\n
    nx, ny, nz (int): grid size
    '''
    dx = 1./nx # e.g., nx = 1000, dx = 0.001
    dy = 1./ny
    dz = 1./nz
    natom = len(frame)
    # xyz format: standard xyz format: element, x, y, z, without number of atoms and the second comment line
    n = np.zeros((nx+1, ny+1, nz+1), dtype=int)
    print('n-digit list generation information:\nGrid size: {}x{}x{}\nList information: {}'.format(
        nx, ny, nz, type(n)
    ))
    for iatom in frame:
        n[int(iatom[1]/dx)][int(iatom[2]/dy)][int(iatom[3]/dz)] += 1
    
    if np.max(n) > 1:
        
        print('WARNING: present precision of griding cannot generate required n-digits list, re-try with higher grid size!')
        return False
    # return n[nx+1, ny+1, nz+1]
    return n
# from a frame of xyz return a frame of digits, n[nx+1, ny+1, nz+1]
def ndigits_t(trj, nx, ny, nz):
    
    nframe = len(trj) # where trj is a Python dict
    ndigits_t_list = []
    for idx_frame in range(nframe):
        ndigit_frame = ndigitsgen(trj[idx_frame], nx, ny, nz)
        if ndigit_frame:
            ndigits_t_list.append(ndigit_frame)
        else:
            print('WARNING: increase nx, ny and nz to eliminate number larger than 1 in n-digits list!')
            raise TypeError
    return ndigits_t_list
# from a trajectory of xyz return a trajectory of digits, nlist[nt, nx+1, ny+1, nz+1]
def slicegen(list_txyz, axis, index):
    
    '''
    Slice generation\n
    list_txyz (4-d list): list_txyz[t, x, y, z]\n
    axis (str or int): 'x' or 0 will be interpreted as x-axis, 'y' or 1 corresponds to y-axis and otherwise z-axis by default\n
    index (int): the index along axis to slice 4d list
    '''
    if (axis=='x') or (axis==0):
        # (t, y, z)
        #return np.transpose(list_txyz, (1,0,2,3))[index][:][:][:]
        return list_txyz[:,index,:,:]
    elif axis=='y' or (axis==1):
        # (t, x, z)
        #return np.transpose(list_txyz, (2,0,1,3))[index][:][:][:]
        return list_txyz[:,:,index,:]
    else:
        # (t, x, y)
        #return np.transpose(list_txyz, (3,0,1,2))[index][:][:][:]
        return list_txyz[:,:,:,index]
# slice one list of (t, x, y, z) to (t, x, y), (t, x, z) or (t, y, z)
# but it may be better only perform transposition for once instead of performing it every time

def qzt(list_t2d, idx_t_disp):
    
    # developer notes: here use indices to count frames instead of real time interval
    len_t, len_dim1, len_dim2 = np.shape(list_t2d)
    qzt_val = 0.
    for idx_dim1 in range(len_dim1):
        for idx_dim2 in range(len_dim2):
            A = 0.
            B = 0.
            for idx_t in range(len_t - idx_t_disp + 1):
                # kernal loop
                A += list_t2d[idx_t][idx_dim1][idx_dim2]*list_t2d[idx_t+idx_t_disp][idx_dim1][idx_dim2]
                B += list_t2d[idx_t][idx_dim1][idx_dim2]
            qz_val += A/B
    return qzt_val
# calculate one point of q(z, t)
def qz(list_t2d, idx_t_max):
    
    qz_val = []
    for idx_t in range(idx_t_max):
        
        qz_val.append(qzt(list_t2d, idx_t))
    return qz_val
# calculate a series of t of q(z,t)
def q(lmp_trjfile, elements, axis, zlist, idx_t_max, nx = 1000, ny = 1000, nz = 1000):
    
    lmp_txyz = ndigits_t(
        trj = lmp2list(
            filename = lmp_trjfile,
            elements = elements,
            isort = True,
            vc = True
            ), 
        nx = nx, 
        ny = ny, 
        nz = nz
        )
    q_val = []
    for z in zlist:
        # zlist is also given by fractional z
        
        q_val.append(
                     qz(
                         list_t2d = slicegen(
                             list_txyz = lmp_txyz,
                             axis = axis, 
                             index = int(z*1000)
                             ),
                         idx_t_max = idx_t_max
                     )
                     )
    return q_val
# calculate q, the spatial-time correlation function
