from copy import deepcopy
import _mat_lib as mlib
from hermite_diag_lib import hermite_diag as hdiag

from qr_decomp import qr
from gram_schmidt import gs_orth as gs

import numpy as np

from lu_decomp import ludecomp as lu
from triangle_solv import triang_solv as sbssolv

def rayleigh_ritz_diag(
    mat_in, 
    num_eigen, 
    preconditioner = 'full', 
    diag_mode = 'householder', 
    dj_solver = 'lu',
    sort_eigval = 'lowest', 
    batch_size = -1, 
    conv_thr = 1E-6, 
    conv_calc_mode = 'norm1',
    max_iter = 50, 
    verbosity = 'silent'
    ):

    """
    Subspace diagonalization (Rayleigh-Ritz subspace method)\n
    # input requirement\n
    mat_in: must be SQUARED matrix and in FLOAT data type\n
    num_eigen: number of eigen vectors and values want to find\n
    preconditioner: preconditioner of residual vector, avaliable options: 'full', 'single', 'dj' or 'none'.\n
    >For 'full' mode (recommended, most stable), (D_A - theta_i*I)|t_i> = |r_i>\n
    >, where DA is diagonal matrix that only has non-zero element on its diagonal, D_A[i][i] = A[i][i]\n
    >For 'single' mode (simple but always cannot converge), (A[i][i] - theta_i)|t_i> = |r_i>\n
    >For 'dj' (Davidson-Jacobi) mode (accurate but singluarity-unstable):\n
    > (I-|y_i><y_i|)(D_A - theta_i*I)(I-|y_i><y_i|)|t_i> = |r_i>\n
    > |t_i> will be solved by LU-decomposition and forward/back substitution method, relatively time-costly\n
    >For 'none' mode, preconditioner won't be used, i.e.: |t_i> = |r_i>\n
    diag_mode: 'jacobi', 'householder' or 'np' (numpy integrated). Basic algorithm for diagonalize matrix in subspace\n
    dj_solver: 'lu' or 'np', the most two fast algorithm for solving linear equation\n
    >For 'lu', use LU-decompsition and forward/backsubstitution\n
    >For 'np', use numpy.linalg.solve function\n
    batch_size: total number of dimensions of subspace, only will be read-in if mode is set to
    'batch'\n
    sort_eigval: 'lowest', 'highest' or 'None'\n
    >>For 'lowest', sort eigenvalues in an increasing order\n
    >>For 'highest', sort eigenvalues in an decreasing order\n
    >>For 'None', will not sort eigenvalues\n
    conv_thr: convergence threshold of eigen values for 'batch' mode, has no effect in other modes\n
    conv_calc_mode: 'abs', 'sqr', 'sum', 'norm1' or 'norm2', for measuring lists of eigenvalues of adjacent two iteration steps\n
    >>For 'abs', use absolute value of difference between element in old list and corresponding one in the new list\n
    >>For 'sqr', use squared value of ...\n
    >>For 'sum', just sum up all differences\n
    >>For 'norm1', use norm of difference between two lists that treated as vectors\n
    >>For 'norm2', only measure difference between norms of vectorized old and new list\n
    max_iter: maximum number of iterations for 'batch' mode, has no effect in other modes\n
    # output description\n
    [eigval_list, eigvec_list]\n
    eigval_list: list of eigenvalues\n
    eigvec_list: list of eigenvectors, arranged according to eigval_list\n
    """
    dim = len(mat_in)
    mat_op_on = deepcopy(mat_in)
    I = np.eye(dim)

    buffer_state = 'YES'
    if batch_size < num_eigen:

        batch_size = num_eigen
        buffer_state = 'NO'

    U = mlib.eye(n = dim, m = batch_size)

    eigval_list0 = mlib.zeros(n = 1, m = num_eigen)[0][:]
    eigval_list = mlib.zeros(n = 1, m = num_eigen)[0][:]

    eigval_conv = 1

    if (verbosity == 'high') or (verbosity == 'debug'):

        print(
             'Diagonalization in subspace Initialization Information\n'
            +'-'*50+'\n'
            +'Preconditioner:                               {}\n'.format(preconditioner)
            +'Number of eigenvalue-vector pairs to find:    {}\n'.format(num_eigen)
            +'If buffer vectors used for batch mode:        {}\n'.format(buffer_state)
            )
        
    istep = 0
    while (istep < max_iter) and (eigval_conv > conv_thr):

        # --------------------subspace generation--------------------
        submat = mlib.unitary_transform(U = U, mat = mat_op_on)
        # -----------------------------------------------------------

        # -----------------subspace diagonalization------------------
        if diag_mode == 'np':

            [eigval_list, eigvec_set] = np.linalg.eig(submat)
        else:

            [eigval_list, eigvec_set] = hdiag(
                hmat_in = submat,
                mode = diag_mode,
                eigval_format = 'list', 
                conv_level = 8,
                max_iter = 50,
                verbosity = verbosity
                )
        # -----------------------------------------------------------

        # ---------subspace eigenvalues and vectors sorting---------
        # sort eigenvalues
        if sort_eigval == 'None':
            # will not sort eigenvalues...
            pass
        elif (sort_eigval == 'lowest') or (sort_eigval == 'highest'):
            # sort eigenvalues anyway...
            sort_idx = np.argsort(
                a = eigval_list,
                axis = -1
            )
            # np.argsort will return a list of indices of elements in list to sort but
            # in a sorted order, ascendent as default
            if sort_eigval == 'highest':

                sort_idx = sort_idx[::-1]
                # reverse the indices list
        # -----------------------------------------------------------

        # --------------------eigenvalues storage--------------------
        templist = [eigval_list[idx] for idx in sort_idx]
        eigval_list = templist

        eigval_conv = mlib.list_diff(
            list_old = eigval_list0,
            list_new = eigval_list[0:num_eigen],
            mode = conv_calc_mode
        )
        eigval_list0 = eigval_list[0:num_eigen]
        # -----------------------------------------------------------


        # -----------------------preprocessing-----------------------
        # rearrange of eigenvectors in subspace
        for idim in range(batch_size):

            templist = [eigvec_set[idim][idx] for idx in sort_idx]
            eigvec_set[idim][:] = templist

        U_new = []
        for ivec in range(batch_size):
                
            s = eigvec_set[:][ivec]                                 # no. ivec subspace eigenvector,       bra
            y = mlib.dot(U, mlib.bra2ket(s), mode = 'matket')       # no. ivec original space Ritz vector, ket
            Resi_mat = mlib.minus(
                mat_op_on,
                mlib.eye(n = dim, m = dim, amplify = eigval_list[ivec])
            )
            r = mlib.dot(
                Resi_mat,
                y,
                mode = 'matket'
            )                                                       # no. ivec residual vector,            ket
            
            if preconditioner == 'full':

                t = []
                for icompo in range(len(r)):

                    ti = r[icompo][0]/(mat_op_on[icompo][icompo]-eigval_list[ivec])
                    t.append(ti)                                    # new vector to append to U,           bra
            elif preconditioner == 'single':

                orig_idx = sort_idx[ivec]
                t = [-ri[0]/(mat_op_on[orig_idx][orig_idx]-eigval_list[ivec]) for ri in r]#                         bra
            elif preconditioner == 'dj':
                r = [[-r[idim][0]] for idim in range(dim)]
                # (I-|y_i><y_i|)(D_A - theta_i*I)(I-|y_i><y_i|)|t_i> = |r_i>
                perp_op = mlib.minus(
                    I, 
                    mlib.ketbra(
                        ket = y,
                        bra = y,
                        mode = '2ket',
                        amplify = 1.0
                    )
                    )
                # preconditioner of full mode
                full_prcdtnr = mlib.zeros(dim, dim)
                for idim in range(dim):
                    full_prcdtnr[idim][idim] = mat_op_on[idim][idim] - eigval_list[ivec]
                # final assembly
                dj_op = mlib.unitary_transform(U = perp_op, mat = full_prcdtnr)
                # solve D|t> = |r>
                if dj_solver == 'lu':
                    
                    if verbosity == 'high':
                        print('RAYLEIGH-RITZ| Start LU-decomposition...\nRAYLEIGH-RITZ| LU-decomposition carried out on matrix:')
                    [_, L_dj, U_dj] = lu(mat_in = dj_op)
                    if verbosity == 'high':
                        print('RAYLEIGH-RITZ| Start forwardsubstitution...')
                    y_dj = sbssolv(
                        triang_mat = L_dj,
                        b = mlib.ket2bra(r),
                        mode = 'lower'
                    )
                    if verbosity == 'high':
                        print('RAYLEIGH-RITZ| Start backsubstitution...')
                    t = sbssolv(
                        triang_mat = U_dj,
                        b = y_dj,
                        mode = 'upper'
                    )                                                 # new vector to append to U,           bra
                elif dj_solver == 'np':

                    t = np.linalg.solve(dj_op, mlib.ket2bra(r))
                if verbosity == 'high':

                    print('RAYLEIGH-RITZ| New |t> generated!')
            elif preconditioner == 'none':

                t = mlib.ket2bra(r)
            else:

                print('RAYLEIGH-RITZ| ***WARNING***: preconditioner required is not recognized, use \'none\' instead.')
                t = mlib.ket2bra(r)

            t = mlib.normalize(t)
            U_new.append(t)
            
        U_new = mlib.transpose(U_new)
        #[_, U, _] = qr(mat_in = U_new, mode = 'householder', verbosity = verbosity)
        #[U, _] = np.linalg.qr(U_new)
        [_, U] = gs(U_new)
        istep += 1
        if verbosity != 'silent':

            print('RAYLEIGH-RITZ| Step {}: conv = {}, conv_thr = {}'.format(istep, eigval_conv, conv_thr))

    return [eigval_list[0:num_eigen], U[:][0:num_eigen]]
