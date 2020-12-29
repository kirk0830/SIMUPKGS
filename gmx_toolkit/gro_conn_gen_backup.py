'''
    # generate bond list:
    bond_data_dict = {}
    # data structure:
    # bond_data_dict['1 2'] = [1, 2, kb, lb]
    # store key in value is for later formatted output
    for iatom in range(natom):
        for jatom in range(natom):
            # warning: in GROMACS, atomic number starts from 1, rather than 0
            pair = str(iatom+1)+' '+str(jatom+1)
            rev_pair = str(jatom+1)+' '+str(iatom+1)
            try:
                bond_data_dict[pair]
                # if call value of present key successfully, this data is already stored, skip
                continue
            except KeyError:
                try:
                    bond_data_dict[rev_pair]
                    # if call value of present key successfully, this data is already stored, skip
                    continue
                except KeyError:
                    # a new line will stored in this dict that will be printed to external file with
                    # required format...
                    bond_data_dict[pair] = [str(iatom+1), str(jatom+1)]
            # create a new bond type in dict, add data by using append!

            if bond_bool_matrix[iatom][jatom] == 1:
                pairName = pair_matrix[iatom][jatom][:] # ['C', 'O']
                pair2check = pairName[0]+' '+pairName[1]
                rev_pair2check = pairName[0]+' '+pairName[1]
                try:
                    data2append = bond_dict[pair2check]
                except KeyError:
                    try:
                        data2append = bond_dict[rev_pair2check]
                    except KeyError:
                        print('SIMUPKGS| ***error*** conflict found here!'
                             +'          1. these two atoms bond with each other\n'
                             +'          2. there is no data for them two in configure file\n'
                             +'          QUIT.')
                        exit()
                bond_data_dict[pair] = bond_data_dict[pair] + data2append
'''