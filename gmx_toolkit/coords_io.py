import numpy as np
import re

def readcoords(filename = None, f_type = 'xyz', convert = True):

    """
    read standard xyz file and convert all coordinates into float data type (optional).\n
    # Input requirement\n
    filename: string\n
    f_type: file type of coordinates, options: xyz. lmp.dump and gro\n
    (for gro file, coodinate inforamtion, residue information and cell parameters will be retured)\n
    (also note that unit of coorinate in gro file is nm, not traditionally Angstrom, output will convert
    data from nm to Angstrom automatically when set 'convert' to True)\n
    convert: bool, select if convert data type from str to float, default value is True\n
    # Output description\n
    2-dimensional list, may contain element symbol in the first column, but must contain 3-dimensional coordinates in last 3 column.
    """
    coords_inp = []

    gro_scale = 1

    if f_type == 'xyz':

        print('COORD| standard xyz file read-in function is activated. Parsing file: '+filename)
        if filename == None:
            print('COORD| ***error*** no coordinate file is read, program is terminated.')
            exit()
        try:
            fileTag = open(filename, 'r', encoding='utf-8')
            print('COORD| the first two lines in *.xyz file will be omitted. Although it may be a\n'
                +'       good way to check completeness of file.')
        except FileNotFoundError:
            print('COORD| ***error*** no coordinate file is found, program is terminated.')
            exit()

        lineskip = 2
        iline = fileTag.readline()
        lineind = 1
        while iline:
            if lineind>lineskip:
                fragments = re.split(' |\n', iline)
                words = [word for word in fragments if word != '']
                coords_inp.append(words)
            
            iline = fileTag.readline()
            lineind += 1

        fileTag.close()

    elif f_type == 'lmp.dump':

        try:
            infile = open(filename, 'r', encoding='utf-8')
            print('COORD| LAMMPS coordinate file read-in: '+str(filename))
        except FileNotFoundError:
            print('COORD| ***error*** lammps coordinate file error: file is non-exist. quit.')
            exit()

        # wash
        # read till find mass information:
        line = infile.readline()
        words = re.split(' |\n', line)
        
        while words[0] != 'Masses':

            line = infile.readline()
            words = re.split(' |\n', line)
        
        print('COORD| #interact# for lammps coordinate files, you need to specify atom types manually. please input their\n'
            +'       element symbols, sperate them with space, so that can be directly used to write standart xyz file.')

        # print till meet atom coordinates:
        while words[0] != 'Atoms':

            line = infile.readline()
            words = re.split(' |\n', line)
            if len(words) >= 3:
                print('       -> element: '+line[:-1])

        elements = input('COORD| waiting for element symbols input: ').split(' ')
        elements = [element for element in elements if element != '']
        # proceed coordinates:
        coords_inp = []
        iloop = 0
        while line and iloop < 50:

            line = infile.readline()
            words = re.split(' |\n', line)
            if len(words) >= 3:

                words = [word for word in words if word != '']
                line2coord = [elements[int(words[1])-1], float(words[-3]), float(words[-2]), float(words[-1])]
                coords_inp.append(line2coord)

        infile.close()

    elif f_type == 'gro':

        gro_scale = 10

        print('COORD| GROMACS format file read-in function is activated. Parsing file: '+filename)
        if filename == None:
            print('COORD| ***error*** no coordinate file is read, program is terminated.')
            exit()
        try:
            fileTag = open(filename, 'r', encoding='utf-8')
            print('COORD| the first two lines in *.gro file will be omitted. Although it may be a\n'
                +'       good way to check completeness of file.')
        except FileNotFoundError:
            print('COORD| ***error*** no coordinate file is found, program is terminated.')
            exit()

        lineskip = 2
        iline = fileTag.readline()
        lineind = 1
        while iline:
            if lineind>lineskip:
                fragments = re.split(' |\n', iline)
                words = [word for word in fragments if word != '']
                if len(words) <= 6:
                    words = [words[1], words[-3], words[-2], words[-1]]
                    coords_inp.append(words)
                elif len(words) == 9:
                    cell_param = words
            
            iline = fileTag.readline()
            lineind += 1

        fileTag.close()

# final process

    if len(coords_inp) == 0:
        print('COORD| ***error*** there is no valid atomic information in *.xyz file, please check your file. QUIT.')
        exit()
        
    if convert:
        print('COORD| read-in complete, converting coordinates from str to float. Here atoms will be counted.')
        coordShape = np.shape(coords_inp)
        print('COORD| there are '+str(coordShape[0])+' atoms in total.')

        for iatom in range(coordShape[0]):
            for iaxis in range(1,4):
                coords_inp[iatom][iaxis] = float(coords_inp[iatom][iaxis])*gro_scale
            
        print('COORD| conversion complete.')
        
    else:
        print('COORD| read-in complete, remember that convert data from str to float-type for future use.')

    return coords_inp

def writecoords(
                coordlist, 
                projectname = 'simu_pkg',
                f_type = 'XYZ',
                res_parti = {},
                cell_param = []
                ):

    """
    output coordinates to external standard formatted xyz file.\n
    # input requirement\n
    coordlist: 2-dimensional list in the format of output of function readcoords\n
    projectname: any string that will be prefix of xyz file\n
    full filename output: projectname-pos-out.xyz\n
    f_type: supported format: standard xyz file -> xyz\n
                              GROMACS gro file -> gro\n
    resi_parti: if file_format is set to gro, this parameter must be specified! Residue information of every atom\n
                data format: {'Residue_A': [1, 3, 142, 153, ...], 'Residue_B': [2, 5, 124, ...], ...}\n
                a dictionary must be specified! If not, residue information of all atoms will be the same...\n
    cell_param: if file_format is set to gro, this parameter must be specified! Cell parameters will be printed\n
                at the end of file, just like what regular gro file looks like.
    """
    natom = np.shape(coordlist)[0]

    if np.shape(coordlist)[1] == 3:

        print('COORD| ***warning*** 3 column-atomic coordinates list detected! Please make sure that it'
             +'      is what you want. A raw will be printed out for check:')
        print(coordlist[0])
        ifquit = input('COORD| ***warning*** you can decide now if quit program. (y/n)')

        if ifquit == 'y':

            exit()
            print('COORD| ***error*** program quit from CoordWrite subroutine.')

        else:

            print('COORD| ***warning*** program will proceed continuely. Coordinates will be treated as'
                 +'      single-element, element symbol you need to input...')
            element = input('COORD| ***warning*** waiting for an element symbol input, all elements are supported.')

            if element:

                coordoutput = []
                for i in range(natom):

                    line = coordlist[i][:]
                    line.insert(0, element)
                    coordoutput.append(line)
    
    elif np.shape(coordlist)[1] == 4:

        coordoutput = coordlist

    else:

        print('COORD| ***error*** bad input format. program stoppted at CoordWrite routine. quit.')
        exit()

    print('COORD| write coordinates to file. name after '+str(projectname)+' , standard xyz file(s) output.')

    filename = projectname+'-pos-out.xyz'
    filetag = open(filename, 'a+', encoding='utf-8')

    filetag.writelines(str(natom)+'\n')
    filetag.writelines('Generated by SIMUPKGS coords_io.coordWrite routine.\n')
    for iline in range(natom):
        if f_type == 'xyz' or f_type == 'standard' or f_type == '':
            filetag.writelines(
                "%3.3s"%coordoutput[iline][0]+' '
                +"%20.10f"%coordoutput[iline][1]
                +"%20.10f"%coordoutput[iline][2]
                +"%20.10f"%coordoutput[iline][3]
                +'\n'
                )
        elif f_type == 'gro':
            # GROMACS gro file output is activated! Add the last line as cell parameters!
            pass
    
    filetag.close()
