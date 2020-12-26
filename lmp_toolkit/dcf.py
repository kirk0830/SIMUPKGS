# displacement correlation function
# ---------------------------------
# this program will use output file from command:
# dump dcf all atom dumpfreq dump.trajectory.trj
# -----------------------------------------------

# STEP1: input file information required...
filename = 'dump.trajectory.trj'
dump_scale_mode = True
# STEP2: MD run parameters required...
dumpfreq = 50000
dt = 0.5
# STEP3: parameters required for this program...
index_i = 1
extractList = [index_i, 6, 23]
extractMode = True
listFromFile = True
listFilename = 'id.txt'
window_width = 20
# STEP4: file name of output file...
dcfFilename = 'dcf-out.dat'

def is_inlist(atomindex, list):

    for item in list:
        if atomindex == item:
            return item
    return False

step = dt * dumpfreq
if listFromFile:
    idfiletag = open(listFilename, 'r', encoding = 'utf-8')
    words = idfiletag.readlines()
    words = [int(word.split('\n')[0]) for word in words]
    words.insert(0, index_i)
    idfiletag.close()
    extractList = words
# atom type
Pd_type_num = 2
# to recognize
head_of_str = str(index_i)+' '+str(Pd_type_num)+' '
# so that one can use 'startswith' function to extract its coordinates.
frame_firstline = 'ITEM: TIMESTEP'
frame_thirdline = 'ITEM: NUMBER OF ATOMS'
frame_fifthline = 'ITEM: BOX BOUNDS xy xz yz pp pp pp'

line = 'go'
filetag = open(filename, 'r', encoding = 'utf-8')
natom = -1
frame_read = 0

coords_hist = {}
dcf_dict = {}
dcf_result_dict = {}
if extractMode:
    # initialize coords_hist
    for iatom in extractList:
        coords_hist[iatom]=[]
        if iatom != index_i:
            dcf_dict[iatom]=[]
            dcf_result_dict[iatom]=0

    print('SIMUPKGS-LAMMPS| complete intialization of coords_hist and dcf_dict')

while line:

    if line.startswith(frame_firstline):
        frame_read += 1
        # start with a new frame
        print('SIMUPKGS-LAMMPS| parse a new frame: No. '+str(frame_read))
        line = filetag.readline() # content: exact number of TIMESTEP
        line = filetag.readline() # content: ITEM: NUMBER OF ATOMS
        line = filetag.readline() # content: exact number of atoms
        if frame_read == 1:
            words = line.split('\n')
            natom = int(words[0])
            print('SIMUPKGS-LAMMPS| There are '+str(natom)+' atoms in total.\n'
                 +'                 number of atoms calculated for correlation: '+str(len(extractList)))
        line = filetag.readline() # content: ITEM: BOX BOUNDS xy xz yz pp pp pp
        line = filetag.readline() # content: xlo, xhi, tilt
        words = line.split(' ')
        xlo = float(words[0])
        xhi = float(words[1])
        X = xhi - xlo
        line = filetag.readline() # content: ylo, yhi, tilt
        words = line.split(' ')
        ylo = float(words[0])
        yhi = float(words[1])
        Y = yhi - ylo
        line = filetag.readline() # content: zlo, zhi, tilt
        words = line.split(' ')
        zlo = float(words[0])
        zhi = float(words[1])
        Z = zhi - zlo
        line = filetag.readline() # content: ITEM: ATOMS id type xs ys zs
        if natom == -1:
            print('SIMUPKGS-LAMMPS| ***I/O error, negative number of atoms!')
            exit()
        for iatom in range(natom):

            line = filetag.readline().split('\n')[0]
            words = line.split(' ')
            index = int(words[0])
            atomtype = int(words[1])
            xs_atom = float(words[2])
            ys_atom = float(words[3])
            zs_atom = float(words[4])

            if dump_scale_mode:
                x_atom = X*xs_atom
                y_atom = Y*ys_atom
                z_atom = Z*zs_atom
            else:
                x_atom = xs_atom
                y_atom = ys_atom
                z_atom = zs_atom

            icoord = [x_atom, y_atom, z_atom]
            #print('atom number: '+str(index)+' present coordinates parsed: '+str(icoord))
            if extractMode:
                iatom_flag = is_inlist(index, extractList)
                if iatom_flag:
                    coords_hist[iatom_flag].append(icoord)
            else:
                print('SIMUPKGS-LAMMPS| ***error, not implemented yet.')
                exit()
        line = filetag.readline()
    elif frame_read == 0:
        line = filetag.readline()
    else:
        print('SIMUPKGS-LAMMPS| ***line read error! present line:')
        print(line)
        exit()

print('-------------trajectory file parsing complete---------------')
filetag.close()

tot_frame = frame_read
num_window = tot_frame - window_width + 1
# upper is file I/O
# now we get dataset.
jatoms = extractList[1:]

def getDisplacement(r, r0):

    x = r[0] - r0[0]
    y = r[1] - r0[1]
    z = r[2] - r0[2]
    return [x, y, z]

from math import sqrt
def correlation_core(r1, r2):

    scalar_prod = r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2]
    norm_r1 = sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2])
    norm_r2 = sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2])
    if scalar_prod == 0 and (norm_r1 == 0 or norm_r2 == 0):
        res = 1
        print('SIMUPKGS-LAMMPS| ***warning: 0/0 exception emerges! use 0/0 = 1 to skip...')
    else:
        res = scalar_prod/norm_r1/norm_r2
    return res

def getDistance(r1, r2):

    r = getDisplacement(r1, r2)
    dist = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
    return dist

def averageOverWindow(twoDimList, nrow = num_window, ncol= window_width - 1):
    # across line
#    print('SIMUPKGS-LAMMPS| averaging across time windows function is called...\n'
#         +'                 present number of windows: '+str(nrow))
    result = []
    for icol in range(ncol):
        row_val = 0
        for irow in range(nrow):
            row_val += twoDimList[irow][icol]/nrow
        result.append(row_val)
    return result

init_dist_list = []
# time correlation function definition recall:
# <B(0)C(t)> = 1/T*integration{0,T, B(x)C(x+t)dx}
for iwindow in range(num_window):
    # ith window
    # t0
    r0_i = coords_hist[index_i][iwindow][:]
    r_i = coords_hist[index_i][iwindow+1][:]
    d_i = getDisplacement(r=r_i, r0=r0_i)
    for index_j in jatoms:
        dcf = []
        r0_j = coords_hist[index_j][iwindow][:]
        if iwindow == 0:
            distij = getDistance(r0_i, r0_j)
            init_dist_list.append(distij)
            print('SIMUPKGS-LAMMPS| initial distance between particle '
                 +str(index_i)
                 +' and '
                 +str(index_j)
                 +': '
                 +str(distij)+' Angstrom(s)')
        # select one atom from extractlist
        for iframe in range(window_width):
            r_j = coords_hist[index_j][iwindow+iframe][:]
            #print('                 frame '+str(iframe)+', rj: '+str(r_j))
            if iframe != 0:
                d_j = getDisplacement(r=r_j, r0=r0_j)
                dcf_value = correlation_core(d_i, d_j)
                dcf.append(dcf_value)
        dcf_dict[index_j].append(dcf)
        #print('SIMUPKGS-LAMMPS| check data size: length of dcf of one single time window: '+str(len(dcf)))

from os.path import isfile
from os import remove


if isfile(dcfFilename):
    remove(dcfFilename)

dcffileTag = open(dcfFilename, 'a+', encoding = 'utf-8')
# dcfFile format: atom id, initial distance, dcf values

idist = 0
from re import split as rsplit
for j_index in jatoms:

    dcf_result_dict[j_index] = averageOverWindow(twoDimList = dcf_dict[j_index])
    line = str(j_index) +' '+str(init_dist_list[idist]) + ' ' + str(dcf_result_dict[j_index])+'\n'
    words = rsplit(',|\[|\]', line)
    for word in words:
        dcffileTag.writelines(' '+word)
    idist += 1

dcffileTag.close()
print('-'*100)
print('SIMUPKGS-LAMMPS| addtional information: if you want to plot correlation function, remember timestep:\n'
     +'                 dt = '+str(step)+' fs, total length of window = '+str(step*window_width)+' fs.')