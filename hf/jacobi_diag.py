import _mat_lib as mlib
from copy import deepcopy
from diag_so2 import diag as orthso2

def op_gen(size, U, irow, icol):

    # generate one Jacobi operator by reading one SO(2) operator
    # and one (irow, icol)-pair, where irow >= icol
    P = mlib.eye(n = size)
    P[icol][icol] = U[0][0]
    P[icol][irow] = U[0][1]
    P[irow][icol] = U[1][0]
    P[irow][irow] = U[1][1]

    return P

def ijacobi(sym_mat_in, tri_diag = False, skipthr = 1E-10, verbosity = 'silent'):

    nline = len(sym_mat_in)
    mat_op_on = deepcopy(sym_mat_in)

    op_accum = mlib.eye(n = nline)

    if tri_diag:
        row_skip = 1
    else:
        row_skip = 0

    for icol in range(nline-1):
        for irow in range(icol+1+row_skip, nline):

            if abs(mat_op_on[irow][icol]/nline**2) < skipthr:
                continue
            # guarantee irow >= icol, i.e., only concentrate on lower triangonal part
            sub_mat = [
                [mat_op_on[icol][icol], mat_op_on[icol][irow]],
                [mat_op_on[irow][icol], mat_op_on[irow][irow]]
                ]
            
            [diag_mat, Uso2] = orthso2(sub_mat)
            if verbosity == 'debug':
                print('\nJacobi (sequential 2-dimensional) diagonalization report\n'+'-'*50)
                print('Present matrix:')
                mlib.matrix_print(mat_op_on, decimal = 4)
                print('Matrix to diagonalize:')
                mlib.matrix_print(sub_mat, decimal = 4)
                print('Diagonalized matrix:')
                mlib.matrix_print(diag_mat, decimal = 4)
                print('Unitary operator at present step:')
                mlib.matrix_print(Uso2, decimal = 4)
            
            
            #mat_op_on[icol][icol] = diag_mat[0][0]
            #mat_op_on[icol][irow] = diag_mat[0][1]
            #mat_op_on[irow][icol] = diag_mat[1][0]
            #mat_op_on[irow][irow] = diag_mat[1][1]

            op = op_gen(size = nline, U = Uso2, irow = irow, icol = icol)

            mat_op_on = mlib.unitary_transform(U = op, mat = mat_op_on)
            
            op_accum = mlib.dot(op_accum, op)

            if verbosity == 'debug':
                print('Unitary operator that operates on whole matrix:')
                mlib.matrix_print(op, decimal = 4)
                print('Accumulative unitary operator:')
                mlib.matrix_print(op_accum, decimal = 4)
    
    if verbosity == 'debug':
        print('-'*30)
        print('JACOBI-DIAG| Final report:\nDiagonalized matrix:')
        mlib.matrix_print(mat_op_on, decimal = 4)
        print('Accumulative unitary operator:')
        mlib.matrix_print(op_accum, decimal = 4)
    
    return [sym_mat_in, mat_op_on, op_accum]

def _off_diag_abs_sum(mat):

    # in-build function, sum up all off-diagonal elements
    nrow = len(mat)
    ncol = len(mat[-1][:])
    if nrow != ncol:

        raise TypeError
    S = 0
    for irow in range(nrow):
        for icol in range(ncol):
            if irow != icol:
                S += abs(mat[irow][icol])
    return S

def jacobi_diag(mat_in, tri_diag = False, conv_level = 3, max_iter = 10, verbosity = 'default'):

    """
    Jacobi diagonalization (sequential SO2-diagonalization, cyclic, iterative)\n
    # input requirement\n
    mat_in: real, symmetric, squared matrix\n
    tri_diag: if set to True, will return a tri-diagonalized matrix, rather than fully
    diagonalized one\n
    conv_level: convergence level, an positive integar is wanted. The higher this value
    is set, the clearer the off-diagonal parts of produced matrix will be\n
    max_iter: although this algorithm is promised to converge, maximum iteration step is
    still used to prohibit a high time-consuming but in vain calculation\n
    verbosity: if set to anything except 'silent', iteration information will be printed
    on screen. if set to 'debug', huge amount of information will be printed, use with
    caution!\n
    # output description\n
    [orig_mat, diag, U]\n
    orig_mat: original matrix\n
    diag: diagonalized matrix\n
    U: unitary operator\n
    # formula\n
    mat_in = U·diag·U', where U' denotes transpose of U
    """
    if not mlib.symm_check(mat = mat_in):

        raise TypeError

    if tri_diag:

        tri_diag_str = '-'*25+'Tri-diagonalization mode'+'-'*25
    else:

        tri_diag_str = ''

    conv_thr = 10.0**(-conv_level)
    sod_log = [] # s.o.d: sum of off-diagonal elements

    isweep = 0
    [_, rawDiag, op] = ijacobi(
        sym_mat_in = mat_in,
        tri_diag = tri_diag,
        skipthr = conv_thr,
        verbosity = verbosity
    )

    op_accum = op

    sod = _off_diag_abs_sum(mat = rawDiag)
    conv = sod
    sod_log.append(sod)

    if verbosity != 'silent':
        print('JACOBI-DIAG| Iteration {}:\n'.format(isweep)
             +'             Sum of all off-diagonal elements S = {}\n'.format(sod)
             +'             Convergence:                         {}\n'.format(conv)
             +'             Convergence threshold:               {}\n'.format(conv_thr)
             +tri_diag_str
             )

    while (conv > conv_thr) and (isweep < max_iter):

        isweep += 1

        [_, rawDiag, op] = ijacobi(
            sym_mat_in = rawDiag,
            tri_diag = tri_diag,
            verbosity = verbosity
        )

        sod = _off_diag_abs_sum(mat = rawDiag)
        conv = abs(sod - sod_log[-1])
        sod_log.append(sod)
        op_accum = mlib.dot(op_accum, op)

        if verbosity != 'silent':
            print('JACOBI-DIAG| Iteration {}:\n'.format(isweep)
                 +'             Sum of all off-diagonal elements S = {}\n'.format(sod)
                 +'             Convergence:                         {}\n'.format(conv)
                 +'             Convergence threshold:               {}\n'.format(conv_thr)
                 +tri_diag_str
                 )
    
    if (conv > conv_thr) and (isweep == max_iter):

        print('JACOBI-DIAG| ***WARNING*** Diagonalization non-converged!')
    
    return [rawDiag, op_accum]
