# read coords from lammps
def transpose(M, mode = 'hard'):
    
    if mode == 'hard':
        nrow = len(M)
        ncol = len(M[0])
        Mt = []
        #print('nrow, ncol = {}, {}'.format(nrow, ncol))
        for icol in range(ncol):
            Mt_line = []
            for irow in range(nrow):
                Mt_line.append(M[irow][icol])
            Mt.append(Mt_line)
        return Mt
    else:
        print('Soft mode of transposing matrix, just interchange subscripts ij to ji')
        return M

def delete_sortcol(M):
    
    return transpose(transpose(M)[1::])

def lmp2dict(filename, elements=[], isort = True):
    
    frames = {}
    idx_frame = -1
    xyz_frame = []
    cell_parameters = [100., 100., 100., 90., 90., 90.]
    # plan to output all frames in one dict, whose keys are indices of frames, 1, 2, 3, ...
    
    natom = 0
    with open(filename, mode = 'r', encoding='utf-8') as LMPf:
        
        line = 'start'
        while line:
            line = LMPf.readline()[0:-1] # first line of frame read: ITEM: TIMESTEP
            if line.startswith('ITEM:'):
                # save and re-initialization
                if idx_frame >= 0:
                    if isort:
                        # sort all lines according to the first column
                        xyz_frame = sorted(xyz_frame, key = lambda x: x[0])
                    xyz_frame = delete_sortcol(xyz_frame)
                    frames[idx_frame] = xyz_frame
                idx_frame += 1
                xyz_frame = []
                cell_parameters = [100., 100., 100., 90., 90., 90.]
                # but actually non-orthodosal box is not implemented
                # read a new header
                header = line[6::]+' '
                line = LMPf.readline()[0:-1] # read the second line: 0
                header += line+' | '
                line = LMPf.readline()[0:-1] # read the third line: ITEM: NUMBER OF ATOMS
                header += line[6::]+' '
                line = LMPf.readline()[0:-1] # read the fourth line: 1728
                natom = int(line)
                header += line+' | '
                line = LMPf.readline()[0:-1] # read the fifth line: ITEM: BOX BOUNDS xy xz yz pp pp pp
                line = LMPf.readline()[0:-1] # 0.0000000000000000e+00 3.9103107999999999e+01 0.0000000000000000e+00

                cell_parameters[0] = float(line.split(' ')[1])
                line = LMPf.readline()[0:-1] # 0.0000000000000000e+00 2.7933903000000001e+01 0.0000000000000000e+00
                cell_parameters[1] = float(line.split(' ')[1])
                line = LMPf.readline()[0:-1] # 0.0000000000000000e+00 2.7933903000000001e+01 0.0000000000000000e+00
                cell_parameters[2] = float(line.split(' ')[1])
                # read 6-8th line
               
                line = LMPf.readline()[0:-1] # read the nineth line: ITEM: ATOMS id type xs ys zs
            # complete reading header of frame
            else:
                # read regular coords info
                words = line.split(' ')
                words = [word for word in words if word != '']
                if len(words) == 0:
                    continue
                xyz_line = [
                    int(words[0]),
                    elements[int(words[1])-1], 
                    float(words[2])*cell_parameters[0], 
                    float(words[3])*cell_parameters[1], 
                    float(words[4])*cell_parameters[2]
                    ]
                xyz_frame.append(xyz_line)
                    
    return frames
