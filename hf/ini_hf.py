
import numpy as np
import gau_intg_lib as gint

def overlap_matrix(basis_dict, coords):

    """
    overlap matrix integral\n
    # input requirement\n
    basis_dict: basis information that has been saved in a dictionary\n
    coords: atomic coordinates that have been saved in a 2-d array\n
    # output description\n
    Sij: overlap matrix, shape in (natom, natom)

    """
    # first get the number of atoms as the number of rows and columns of overlap matrix
    # and for each atom, there will be a independent number of Gaussian functions needed
    # to fit atomic wavefunction
    # it may be good to store basis information in dictionary,
    # and get basis information by calling one certain key
    # therefore the number of atoms: natom = len(basis_dict)
    # number of Gaussian functions for each atom: nbasis_i = len(basis_dict[i])

    # basis information is stored in order: 
    # alpha, d, l, m, n, alpha, d, l, m, n, ...
    natom = len(basis_dict)
    Sij = np.zeros(shape = (natom, natom))
    for irow in range(natom):
        for icol in range(natom):
            for igau in range(0, len(basis_dict[irow]), 5):
                for jgau in range(0, len(basis_dict[icol]), 5):
                    alpha_A = basis_dict[irow][igau]
                    alpha_B = basis_dict[icol][jgau]
                    # make sure that Sij_Gau has been generalized!!!

                    Sij[irow][icol] += gint.Sij_Gau(
                                                    a = alpha_A,
                                                    b = alpha_B,
                                                    atom_A = coords[irow][:],
                                                    atom_B = coords[icol][:],
                                                    aglr_a = basis_dict[irow][igau+2:igau+4],
                                                    aglr_b = basis_dict[icol][jgau+2:jgau+4]
                                                    )*basis_dict[irow][igau+1]*basis_dict[icol][jgau+1]

    return Sij

def kinetic_matrix(basis_dict, coords):

    """
    kinetic matrix integral\n
    # input requirement\n
    basis_dict: basis information that has been saved in a dictionary\n
    coords: atomic coordinates that have been saved in a 2-d array\n
    # output description\n
    Tij: kinetic matrix, shape in (natom, natom)

    """
    # first get the number of atoms as the number of rows and columns of kinetic matrix
    # and for each atom, there will be a independent number of Gaussian functions needed
    # to fit atomic wavefunction
    # it may be good to store basis information in dictionary,
    # and get basis information by calling one certain key
    # therefore the number of atoms: natom = len(basis_dict)
    # number of Gaussian functions for each atom: nbasis_i = len(basis_dict[i])

    # basis information is stored in order: alpha, d, alpha, d, ...
    natom = len(basis_dict)
    Tij = np.zeros(shape = (natom, natom))
    for irow in range(natom):
        for icol in range(natom):
            for igau in range(len(basis_dict[irow])):
                for jgau in range(len(basis_dict[icol])):
                    alpha_A = basis_dict[irow][igau]
                    alpha_B = basis_dict[icol][jgau]

                    Tij[irow][icol] += gint.Tij_Gau(
                                                    a = alpha_A,
                                                    b = alpha_B,
                                                    atom_A = coords[irow][:],
                                                    atom_B = coords[icol][:],
                                                    aglr_a = basis_dict[irow][igau+2:igau+4],
                                                    aglr_b = basis_dict[icol][jgau+2:jgau+4]
                                                    )*basis_dict[irow][igau+1]*basis_dict[icol][jgau+1]

    return Tij

def Vne_matrix(basis_dict, coords, chargelist, mode = 'all'):

    """
    nuclear attraction matrix integral\n
    # input requirement\n
    basis_dict: basis information that has been saved in a dictionary\n
    coords: atomic coordinates that have been saved in a 2-d array\n
    mode: specify if only attraction matrix of one atom is calculated. If set to 
    'all', attraction matrices will be summed up altogether, if set to one specific
    number, then only attraction matrix of that atom will be calculated and returned.\n
    # output description\n
    Vneij: nuclear attraction matrix, shape in (natom, natom)
    """
    # first get the number of atoms as the number of rows and columns of kinetic matrix
    # and for each atom, there will be a independent number of Gaussian functions needed
    # to fit atomic wavefunction
    # it may be good to store basis information in dictionary,
    # and get basis information by calling one certain key
    # therefore the number of atoms: natom = len(basis_dict)
    # number of Gaussian functions for each atom: nbasis_i = len(basis_dict[i])

    # basis information is stored in order: alpha, d, alpha, d, ...
    natom = len(basis_dict)
    Vneij = np.zeros(shape = (natom, natom))

    if mode == 'all':
        Vne_sum_lo = 0
        vne_sum_hi = natom
    else:
        Vne_sum_lo = int(mode)
        vne_sum_hi = Vne_sum_lo + 1
    
    for icore in range(Vne_sum_lo, vne_sum_hi):
        for irow in range(natom):
            for icol in range(natom):
                for igau in range(len(basis_dict[irow])):
                    for jgau in range(len(basis_dict[icol])):
                        alpha_A = basis_dict[irow][igau]
                        alpha_B = basis_dict[icol][jgau]
                        
# ----------------------still under construction------------------------
                        # make sure that Vneij_Gau has been generalized!!!
                        Vneij[irow][icol] += gint.Vneij_Gau(
                                                        A = alpha_A,
                                                        B = alpha_B,
                                                        rAB2 = 0,
                                                        rCP2 = 0,
                                                        ZC = chargelist[icore])
# ----------------------still under construction------------------------
    return Vneij

def Vee_matrix(basis_dict, coords):

    """
    electron-electron interaction matrix integral\n
    # input requirement\n
    basis_dict: basis information that has been saved in a dictionary\n
    coords: atomic coordinates that have been saved in a 2-d array\n
    # output description\n
    Veeijkl: electron-electron interaction matrix, shape in (natom, natom, natom, natom)

    """
    # first get the number of atoms as the number of rows and columns of electron-electron
    # interaction matrix and for each atom, there will be a independent number of Gaussian 
    # functions needed to fit atomic wavefunction.
    # it may be good to store basis information in dictionary,
    # and get basis information by calling one certain key
    # therefore the number of atoms: natom = len(basis_dict)
    # number of Gaussian functions for each atom: nbasis_i = len(basis_dict[i])

    # basis information is stored in order: alpha, d, alpha, d, ...
    natom = len(basis_dict)
    Veeijkl = np.zeros(shape = (natom, natom, natom, natom))
    for id1 in range(natom):
        for id2 in range(natom):
            for id3 in range(natom):
                for id4 in range(natom):
                    for igau in range(len(basis_dict[id1])):
                        for jgau in range(len(basis_dict[id2])):
                            for kgau in range(len(basis_dict[id3])):
                                for lgau in range(len(basis_dict[id4])):
                                    alpha_A = basis_dict[id1][igau]
                                    alpha_B = basis_dict[id2][jgau]
                                    alpha_C = basis_dict[id3][kgau]
                                    alpha_D = basis_dict[id4][lgau]
# ----------------------still under construction------------------------
                                    # symmetry should be used here to accelerate calculation
                                    
                                    Veeijkl[id1][id2][id3][id4] += gint.Veeijkl_Gau(
                                        A = alpha_A,
                                        B = alpha_B,
                                        C = alpha_C,
                                        D = alpha_D,
                                        rAB2 = 0,
                                        rCD2 = 0,
                                        rPQ2 = 0
                                    )
# ----------------------still under construction------------------------
    return Veeijkl