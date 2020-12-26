# velocity correlation function
# ------------------------------
# illustration for using vcf.py:
# ------------------------------
# you should explicitly use command:
# dump vcf all custom dumpfreq vcf_dumpfile.dat id type x y z vx vy vz

# STEP1: input file information required...
filename = 'dump.vcf.dat'
# STEP2: MD run parameters required...
dumpfreq = 2000
dt = 0.5
# STEP3: parameters required for this program...
index_i = 1260
extractMode = True
extractList = [index_i, 6, 23]
listFromFile = True
listFilename = 'id.txt'
window_width = 10
velocity_gen = True
inteval_anal = True
window_start = 87
window_end = 99
# STEP4: file name of output file...
vcfFilename = 'vcf-out.dat'

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
vel_hist = {}
vcf_dict = {}
vcf_result_dict = {}
if extractMode:
    # initialize coords_hist
    for iatom in extractList:
        coords_hist[iatom]=[]
        vel_hist[iatom]=[]
        vcf_dict[iatom]=[]
        vcf_result_dict[iatom]=0

    print('SIMUPKGS-LAMMPS| complete intialization of coords_hist and vcf_dict')

while line:

    if line.startswith(frame_firstline):
        frame_read += 1
        # start with a new frame
        #print('SIMUPKGS-LAMMPS| parse a new frame: No. '+str(frame_read))
        line = filetag.readline() # content: exact number of TIMESTEP
        line = filetag.readline() # content: ITEM: NUMBER OF ATOMS
        line = filetag.readline() # content: exact number of atoms
        if frame_read == 1:
            words = line.split('\n')
            natom = int(words[0])
            print('SIMUPKGS-LAMMPS| There are '+str(natom)+' atoms in total.\n'
                 +'                 number of atom-pairs calculated for correlation: '+str(len(extractList)))
        line = filetag.readline() # content: ITEM: BOX BOUNDS xy xz yz pp pp pp
        line = filetag.readline() # content: xlo, xhi, tilt
        line = filetag.readline() # content: ylo, yhi, tilt
        line = filetag.readline() # content: zlo, zhi, tilt
        line = filetag.readline() # content: ITEM: ATOMS id type x y z vx vy vz
        if natom == -1:
            print('SIMUPKGS-LAMMPS| ***I/O error, negative number of atoms!')
            exit()
        for iatom in range(natom):

            line = filetag.readline().split('\n')[0]
            words = line.split(' ')
            index = int(words[0])
            atomtype = int(words[1])
            x_atom = float(words[2])
            y_atom = float(words[3])
            z_atom = float(words[4])
            
            vx_atom = float(words[5])
            vy_atom = float(words[6])
            vz_atom = float(words[7])

            icoord = [x_atom, y_atom, z_atom]
            ivel = [vx_atom, vy_atom, vz_atom]
            #print('atom number: '+str(index)+' present coordinates parsed: '+str(icoord))
            if extractMode:
                iatom_flag = is_inlist(index, extractList)
                if iatom_flag:
                    # create an atom-specific history for every atom that want to calculate correlation for.
                    coords_hist[iatom_flag].append(icoord)
                    vel_hist[iatom_flag].append(ivel)
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
if inteval_anal:
    num_window = window_end - window_start -window_width + 1
    if window_end == -1:
        print('SIMUPKGS-LAMMPS| ***error: haven\'t give correct value to parameter window_end.')
        exit()
else:
    num_window = tot_frame - window_width + 1
    window_start = 0
    window_end = num_window
if num_window <= 0:
    print('SIMUPKGS_LAMMPS| ***error: wrong window definition, check revelant input parameters!')
print('SIMUPKGS-LAMMPS| parse configuration: window_start:      '+str(window_start)+'\n'
     +'                                      window_end:        '+str(window_end)+'\n'
     +'                                      window_width:      '+str(window_width)+'\n'
     +'                                      number of windows: '+str(num_window))
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
        print('SIMUPKGS-LAMMPS| ***warning: 0/0 exception emerges! use 0/0 = 1 to skip...\n'
             +'                          -> this will happen when measuring correlations between the FIRST frame and whatever other\n'
             +'                             frames especially when you didnt generate velocities for particles at first timestep of run,\n'
             +'                             if so, it will be safe to ignore it and will not affect accuracy of result too much if you\n'
             +'                             choose a medium or large number of windows.\n'
             +'                             ***IF NOT, PLEASE RE-CHECK YOUR INPUT FILE!!!')
    else:
        res = scalar_prod/norm_r1/norm_r2
    return res

def getDistance(r1, r2):

    r = getDisplacement(r1, r2)
    dist = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
    return dist

def averageOverWindow(twoDimList, nrow = num_window, ncol= window_width):
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
for iwindow in range(window_start, window_end-window_width+1):
    # ith window
    # t0
    v_i = vel_hist[index_i][iwindow][:]

    for index_j in jatoms:
        vcf = []
        if iwindow == window_start:
            r0_i = coords_hist[index_i][0][:]
            r0_j = coords_hist[index_j][0][:]
            distij = getDistance(r0_i, r0_j)
            init_dist_list.append(distij)
#            print('SIMUPKGS-LAMMPS| initial distance between particle '
#                 +str(index_i)
#                 +' and '
#                 +str(index_j)
#                 +': '
#                 +str(distij)+' Angstrom(s)')
        # select one atom from extractlist
        for iframe in range(window_width):
            try:
                v_j = vel_hist[index_j][iwindow+iframe][:]
            except IndexError:
                print('ERROR INDEX INFORMATION: iwindow = '+str(iwindow)+', iframe = '+str(iframe))
                exit()
            vcf_value = correlation_core(v_i, v_j)
            vcf.append(vcf_value)
        vcf_dict[index_j].append(vcf)
        #print('SIMUPKGS-LAMMPS| check data size: length of vcf of one single time window: '+str(len(vcf)))

from os.path import isfile
from os import remove

if isfile(vcfFilename):
    remove(vcfFilename)

vcffileTag = open(vcfFilename, 'a+', encoding = 'utf-8')
# vcfFile format: atom id, initial distance, vcf values

idist = 0
from re import split as rsplit
for j_index in jatoms:

    vcf_result_dict[j_index] = averageOverWindow(twoDimList = vcf_dict[j_index])
    line = str(j_index) +' '+str(init_dist_list[idist]) + ' ' + str(vcf_result_dict[j_index])+'\n'
    words = rsplit(',|\[|\]', line)
    for word in words:
        vcffileTag.writelines(' '+word)
    idist += 1

vcffileTag.close()
print('-'*100)
print('SIMUPKGS-LAMMPS| addtional information: if you want to plot correlation function, remember timestep:\n'
     +'                 dt = '+str(step)+' fs, total length of window = '+str(step*window_width)+' fs.')