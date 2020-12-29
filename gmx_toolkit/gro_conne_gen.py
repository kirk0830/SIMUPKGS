import coords_io
from distance import dist as dist
import re
import numpy as np

def parse_conf_file(filename, verbosity = 'silent'):

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
        else:
            # for lines dont start with section starting signal '$'
            # line is already read-in
            if bond_parse:
                words = re.split(pattern = ' |\n', string = line)
                words = [word for word in words if word != '']
                _key = words[0]+' '+words[1]
                _value = [float(words[2]), float(words[3])]
                bond_dict[_key] = _value
            elif angle_parse:
                words = re.split(pattern = ' |\n', string = line)
                words = [word for word in words if word != '']
                _key = words[0]+' '+words[1]+' '+words[2]
                _value = [float(words[3]), float(words[4])]
                angle_dict[_key] = _value
            elif dihe_parse:
                words = re.split(pattern = ' |\n', string = line)
                words = [word for word in words if word != '']
                _key = words[0]+' '+words[1]+' '+words[2]+' '+words[3]
                _value = [float(words[4]), float(words[5])]
                dihe_dict[_key] = _value

    return [elements, bond_dict, angle_dict, dihe_dict]

def get_data_from_dict(handle, dataFrom, verbosity = 'silent'):

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
    
    return data

def topol_gen(
    gro_file, 
    param_file, 
    d_tol = 0.1, 
    verbosity = 'debug'
    ):

    # unlike structure file, this may be the only place that will need parse configuration file
    # that contains settings about generating topol information from gro.
    [_, bond_dict, angle_dict, dihe_dict] = parse_conf_file(filename = param_file)

    # read-in structure and parse file, convert into array
    data_import = coords_io.readcoords(filename='template.gro', f_type='gro', convert=True)
    natom = coords_io.np.shape(data_import)[0]

    dist_matrix = []
    pair_matrix = []
    for iatom in range(natom):

        look_at_this = data_import[iatom][1::]
        [dist_array_iatom, pair_array_iatom] = dist(
                                                    atomlist = data_import, 
                                                    pointxyz = look_at_this, 
                                                    append_ele = data_import[iatom][0]
                                                    )
        dist_matrix.append(dist_array_iatom)
        pair_matrix.append(pair_array_iatom)

    # dist_matrix: (Natom, Natom) shaped
    # pair_matrix: (Natom, Natom, 2) shaped

    # generate connectivity matrix
    bond_bool_matrix = np.zeros(shape = (natom, natom))

    for iatom in range(natom):
        for jatom in range(iatom+1, natom):
            # skip diagonal element
            # extract distance info.
            pairLength = dist_matrix[iatom][jatom]
            pairName = pair_matrix[iatom][jatom][:] # ['C', 'O']
            crit_bond = get_data_from_dict(handle = pairName, dataFrom = bond_dict)
            if crit_bond:
                if (pairLength >= crit_bond[1]-d_tol) and (pairLength <= crit_bond[1]+d_tol):
                    # suppose there is a bond!
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
