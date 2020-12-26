# It is meaningless for users to build wheels
# But learn how to build wheels is meaningful to me
# -ykhuang 2020/12/26

import numpy as np

def npDiag(mat, sortFlag = False, verbosity = 'low'):
    
    if verbosity == 'debug':
        print('SIMUPKGS_DEBUG_MODE | matrix_to_diagonalize: '+str(mat))
        print('SIMUPKGS_DEBUG_MODE | matrix shape information: '+str(np.shape(mat)))
        print('SIMUPKGS_DEBUG_MODE | matrix data-type information: '+str(type(mat)))

    mat = np.array(mat, dtype = float)
    [eigenvalsList, U] = np.linalg.eig(mat)

    if sortFlag:

        print('DIAG | diagnolization with sorting eigenvalues increasingly enabled.')
        positionRecord = range(len(eigenvalsList))
        label = [('record', int), ('data', float)]
        eigenvalLabelled = [(positionRecord[i], eigenvalsList[i]) for i in positionRecord ]
        eigenval2sort = np.array(object = eigenvalLabelled, dtype = label)
        eigenvalsSorted = np.sort(eigenval2sort, order = 'data')
        order2rearrange = [ eigenvalsSorted[i][0] for i in positionRecord ]
        eigenvalsList = [ eigenvalsSorted[i][1] for i in positionRecord ]
        U_trans = np.transpose(U)
        U2trans = []
        for iline in order2rearrange:

            U2trans.append(U_trans[iline])
        
        U = np.transpose(U2trans)
        
    eigenMat = np.diag(eigenvalsList)

    return [eigenMat, U]
