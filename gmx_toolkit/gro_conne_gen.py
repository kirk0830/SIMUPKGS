import coords_io
from distance import dist as dist
import re
import numpy as np

def parse_conf_file(filename, verbosity = 'silent'):

    """
    GROMACS *.gro to *.top, structure to topology file conversion\n
    # input requirement:\n
    filename: filename of configuration file that contains information about
    elements, bond, angle, dihedral forcefield parameters.
    """
    # create dictionaries to store data read from configuration file
    bond_dict = {}
    angle_dict = {}
    dihe_dict = {}

    filetag = open(file = filename, encoding = 'utf-8')
    line = 'line0'
    while line:
        line = filetag.readline()
        if line.startswith('$'):
            # read a new section
            # parse keywords of this line
            words = re.split(pattern = '\$|\n| ', string = line)
            words = [word for word in words if word != '']
            if verbosity == 'debug':
                print(line)
                print(words)
            keyword = words[0]
            if keyword == 'ELEMENT':
                # start to parse element section
                line = filetag.readline()
                words = re.split(pattern = ' |\n', string = line)
                elements = [word for word in words if word != '']
            elif keyword == 'BOND':
                bond_parse = True
            elif keyword == 'ANGLE':
                angle_parse = True
            elif keyword == 'DIHEDRAL':
                dihe_parse = True
            elif keyword == 'END':
                bond_parse = False
                angle_parse = False
                dihe_parse = False
        elif line.startswith('%'):
            continue
        else:
            # for lines dont start with section starting signal '$'
            # line is already read-in
            if bond_parse:
                words = re.split(pattern = ' |\n', string = line)
                words = [word for word in words if word != '']
                _key = words[0]+' '+words[1]
                _type = words[2]
                if _type == '1':
                    # basic harmonic type:
                    # format to be printed later in topology file:
                    #     i     j  func  bond_len_nm      k_J/nm2
                    #     3     4     1 1.000000e-01 3.744680e+05
                    _value = [float(words[3]), float(words[4])]
                elif _type == '2':
                    pass
                else:
                    _value = [float(words[3]), float(words[4])]
                bond_dict[_key] = _value
            elif angle_parse:
                words = re.split(pattern = ' |\n', string = line)
                words = [word for word in words if word != '']
                _key = words[0]+' '+words[1]+' '+words[2]
                _type = words[3]
                if _type == '1':
                    # basic harmonic type:
                    # format to be printed later in topology file:
                    #     i     j     k  func   theta0_deg   k_J/theta2
                    #     1     3     4     1 1.200000e+02 2.928800e+02
                    _value = [float(words[4]), float(words[5])]
                elif _type == '2':
                    pass
                else:
                    # default as func == 1
                    _value = [float(words[4]), float(words[5])]
                angle_dict[_key] = _value
            elif dihe_parse:
                words = re.split(pattern = ' |\n', string = line)
                words = [word for word in words if word != '']
                _key = words[0]+' '+words[1]+' '+words[2]+' '+words[3]
                _type = words[4]
                if _type == '1':
                    # dihedral type 1, regular dihedral:
                    # format to be printed later in topology file:
                    #     i     j     k     J  func      phi_deg   k_J/theta2            n
                    #     2     1     3     4     1 1.800000e+02 3.347200e+01 2.000000e+00
                    _value = [float(words[5]), float(words[6]), float(words[7])]
                elif _type == '2':
                    # dihedral type 2, improper dihedral:
                    # format to be printed later in topology file:
                    #     i     j     k     J  func      phi_deg   k_J/theta2
                    #     3     4     5     1     2 0.000000e+00 1.673600e+02
                    _value = [float(words[5]), float(words[6])]

                dihe_dict[_key] = _value

    return [elements, bond_dict, angle_dict, dihe_dict]

def get_data_from_dict(handle, dataFrom, verbosity = 'silent'):
    """
    data search function\n
    # input requiement:\n
    handle: 0. handle should be provided by an array-like data: ['C', 'H'] or something\n
            1. bond type information handle: 'C H': get forcefield parameters of C-H bond\n
            2. angle type information handle: 'H C H': get ... of H-C-H bond\n
            3. dihedral ...: 'H C C H': get ... of H-C-C-H dihedral\n
    dataFrom: Python specific data type dictionary, a dictionary stores forcefield infomation
    """
    key = ''
    rev_key = ''
    handleLen = len(handle)
    for iatom in range(handleLen):
        key += str(handle[iatom])+' '
        rev_key += str(handle[-iatom-1])+' '
    key = key[0:-1]
    rev_key = rev_key[0:-1]

    try:
        data = dataFrom[str(key)]
    except KeyError:
        try:
            data = dataFrom[str(rev_key)]
        except KeyError:
            if verbosity == 'debug':
                print('SIMUPKGS| ***warning*** present measured atomic pair is not record in configuration file.\n'
                    +'          pair: '+str(key))
            return False
    if verbosity == 'everything':
        print('SIMUPKGS_DEBUG_MODE| data extracted from dict:')
        print(data)
        
    return data

def topol_gen(
    gro_file, 
    param_file, 
    d_tol = 0.1, 
    verbosity = 'debug'
    ):
    """
    main function of\n
    conversion from gro and configuration file to topology file\n
    # input requirement:\n
    gro_file: GROMACS structure file\n
    param_file: configuration file that contains setting on how to make topology file\n
    d_tol: tolerance parameter that controls the precision of judging connectivity between
    atoms, where basic criteria are read from configuration file
    """
    # unlike structure file, this may be the only place that will need parse configuration file
    # that contains settings about generating topol information from gro.
    [_, bond_dict, angle_dict, dihe_dict] = parse_conf_file(filename = param_file)
    if verbosity == 'debug':
        print('SIMUPKGS_DEBUG_MODE| bond:')
        print(bond_dict)
        print('SIMUPKGS_DEBUG_MODE| angle:')
        print(angle_dict)
        print('SIMUPKGS_DEBUG_MODE| dihedral:')
        print(dihe_dict)
    # read-in structure and parse file, convert into array
    data_import = coords_io.readcoords(
                                        filename='template.gro', 
                                        f_type='gro', 
                                        convert=True
                                        )
    natom = coords_io.np.shape(data_import)[0]
    # in GROMACS, coordinates are in unit nm, not Angstrom!

    dist_matrix = []
    pair_matrix = []
    for iatom in range(natom):

        look_at_this = data_import[iatom][1::]
        [dist_array_iatom, pair_array_iatom] = dist(
                                                    atomlist = data_import, 
                                                    pointxyz = look_at_this, 
                                                    point_in_unit = 'Angstrom',
                                                    set_in_unit = 'Angstrom',
                                                    dist_in_unit = 'Angstrom',
                                                    pbcFlag = False,
                                                    append_ele = data_import[iatom][0]
                                                    )
        dist_matrix.append(dist_array_iatom)
        pair_matrix.append(pair_array_iatom)

    # dist_matrix: (Natom, Natom) shaped
    # pair_matrix: (Natom, Natom) shaped
    if verbosity == 'debug':
        print('SIMUPKGS_DEBUG_MODE| dist_matrix_output: (off)')
        #print(dist_matrix)
        print('SIMUPKGS_DEBUG_MODE| pair_matrix_output: (off)')
        #print(pair_matrix)

        print('SIMUPKGS| bond_bool_matrix initialization...')
    
    # generate connectivity matrix
    bond_bool_matrix = np.zeros(shape = (natom, natom))

    for iatom in range(natom):
        for jatom in range(iatom+1, natom):
            # skip diagonal element
            # extract distance info.
            pairLength = dist_matrix[iatom][jatom]
            pairName = pair_matrix[iatom][jatom] # ['C', 'O']
            if verbosity == 'everything':
                print('SIMUPKGS_DEBUG_MODE| check parameters input in following get_data_from_dict:\n'
                     +'                     handle = '+str(pairName)+'\n'
                     +'                     dataFrom = bond_dict')
            crit_bond = get_data_from_dict(handle = pairName, dataFrom = bond_dict)
            if crit_bond:
                if (pairLength >= crit_bond[1]-d_tol) and (pairLength <= crit_bond[1]+d_tol):
                    # suppose there is a bond!
                    if verbosity == 'debug':
                        print('SIMUPKGS_DEBUG_MODE| one bond is determined by distance with tolerence value '+str(d_tol))
                        print('                     '+str(pairName)+': critical bond length = {} Angstrom'.format(str(crit_bond[1])))
                        print('                     atomic index pair: {}-{}'.format(iatom, jatom))
                    bond_bool_matrix[iatom][jatom] = 1
                    bond_bool_matrix[jatom][iatom] = 1
    # bond_bool_matrix is generated

    bond_print = []
    angle_print = []
    dihe_print = []

    for iatom in range(natom):
        for jatom in range(natom):

            if jatom == iatom:
                continue

            if bond_bool_matrix[iatom][jatom]:
                # find i-j bond
                # double count is avoided by check wether j is larger than i
                if iatom > jatom:
                    # bond is already counted
                    pass
                else:
                    # this bond is not counted yet!
                    # >>>
                    bond2append = [iatom, jatom]
                    get_this_bond = [data_import[iatom][0], data_import[jatom][0]]
                    write_this_bond = get_data_from_dict(handle = get_this_bond, dataFrom = bond_dict)
                    if write_this_bond:
                        bond2append += write_this_bond
                        bond_print.append(bond2append)

                for katom in range(natom):

                    if (katom == iatom) or (katom == jatom):
                        continue

                    if bond_bool_matrix[jatom][katom]:
                        # find j-k bond and i-j-k angle

                        # for angle counting
                        if katom < iatom:
                            # this angle is already counted
                            pass
                        else:
                            # this angle is not counted yet!
                            # >>>
                            angle2append = [iatom, jatom, katom]
                            get_this_angle = [data_import[iatom][0], data_import[jatom][0], data_import[katom][0]]
                            write_this_angle = get_data_from_dict(handle = get_this_angle, dataFrom = angle_dict)
                            if write_this_angle:
                                angle2append += write_this_angle
                                angle_print.append(angle2append)

                        for latom in range(natom):

                            if (latom == iatom) or (latom == jatom) or (latom == katom):
                                continue

                            if bond_bool_matrix[katom][latom]:
                                # find k-l bond and j-k-l angle and i-j-k-l dihedral
                                if latom < iatom:
                                    # this dihedral is already counted
                                    pass
                                else:
                                    # this dihedral is not counted yet!
                                    # >>>
                                    dihe2append = [iatom, jatom, katom, latom]
                                    get_this_dihe = [
                                        data_import[iatom][0], 
                                        data_import[jatom][0], 
                                        data_import[katom][0],
                                        data_import[latom][0]
                                        ]
                                    write_this_angle = get_data_from_dict(handle = get_this_dihe, dataFrom = dihe_dict)
                                    if write_this_angle:
                                        dihe2append += write_this_angle
                                        dihe_print.append(dihe2append)

    return [bond_print, angle_print, dihe_print]

def topol_print(
                generdict = {},
                includelist = [],
                atomlist = [], 
                bondlist = [], 
                anglelist = [], 
                dihelist = [],
                pos_restr_list = [],
                comments = [],
                filename = ''
                ):

    """
    GROMACS type topology file output file in required format\n
    # input requirement:\n
    atomlist: 2-d array contains information of atom name, atom mass, atom charge, atom c6 and c12 coefficients\n
    bondlist: 2-d array contains all bonds and parameters to define a bond type data line\n
    anglelist, dihelist: 2-d ..., contain all angle/dihedral parameters...\n
    filename: name of the file output\n
    # output description in detail\n
    all data will be identical as they can to form topology file. It is known that standard format
    of topology file is strict in GROMACS, as a summary, format of top file and also gro file are recorded
     here:\n
    1. top file:\n
    6-element space for atomic index, function type\n
    8-element space for type name of atoms, residue number, name of residue, name of atom, charge number
    and charge, where charges are 3-decimal number\n
    2. gro file:\n
    8-element space for atomic coordinates with 3-decimal precision\n

    """
    # I would like to check standard format of topology file firstly,
    # also I remember that GROMACS has strict rule of gro-type file.
    # to see standard format of topology file, plz see urea.itp

    pass