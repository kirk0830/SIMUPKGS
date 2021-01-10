from givens import tri_diag as givens_tdiag
from householder import householder as hh
import _mat_lib as mlib
from copy import deepcopy

def trid(mat_in, mode = 'givens'):

    """
    Tri-diagonalization of symmetrical real matrix\n
    # input requirement\n
    mat_in: symmetrical real matrix, use FLOAT data type, or will emerge huge roundoff error\n
    mode: 'givens' or 'householder', 'householder' works faster than 'givens' for 'givens'
    will seperately work on every off-tri-diagonal element in an iterative manner, not recommended
    for use other than diagonalization, although there are more efficient diagonalization method\n
    # output description\n
    [tridiag, U]\n
    tridiag: tri-diagonalized matrix\n
    U: unitary operator\n
    # formula\n
    mat_in = U·tridiag·U', where U' denotes transpose of U
    """

    if not mlib.symm_check(mat = mat_in):
        
        print('***error*** symmetric matrix is demanded.')
        exit()

    if mode == 'givens':
        
        [tridiag, U] = givens_tdiag(
            mat_in = mat_in,
            verbosity = 'silent'
        )
    elif mode == 'householder':

        nline = len(mat_in)

        mat_op_on = deepcopy(mat_in)
        mat_op_on_T = mlib.transpose(mat_op_on)
        U = mlib.eye(nline)

        for iline in range(nline-1):

            vec_op_on = mat_op_on_T[iline][iline+1::]
            reflect_op = hh(vec_in = vec_op_on, verbosity = 'silent')
            identi_op = mlib.eye(iline+1)
            op = mlib.combine_block(identi_op, reflect_op)

            #mat_op_on = mlib.dot(op, mat_op_on)
            mat_op_on = mlib.unitary_transform(op, mat_op_on)
            mat_op_on_T = mlib.transpose(mat_op_on)
            U = mlib.dot(op, U)

        U = mlib.transpose(U)
        tridiag = mat_op_on
    
    else:

        exit()

    return [tridiag, U]
    
