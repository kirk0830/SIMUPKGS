from householder import householder as hh
from gram_schmidt import gs_orth as gs
import _mat_lib as mlib
from copy import deepcopy

def qr(mat_in, mode = 'householder', verbosity = 'silent'):

    """
    QR decomposition\n
    A = QR, where R is upper-triangonal matrix\n
    # input requirement\n
    mat_in: matrix to perform QR decomposition, only squared matrix is supported\n
    mode: 'householder' (recommended) or 'schmidt' (not recommended, bug exists, unsolved)\n
    # output description\n
    [original matrix, Q matrix, R matrix]
    """
    nline = len(mat_in)
    mat_op_on = deepcopy(mat_in)

    if mode == 'householder':

        mat_op_on_T = mlib.transpose(mat_op_on)
        op_log = []
        for iline in range(nline):

            vec_op_on = mat_op_on_T[iline][iline::]
            reflect_op = hh(vec_in = vec_op_on, verbosity = verbosity)
            identi_op = mlib.eye(iline)
            op = mlib.combine_block(identi_op, reflect_op)
            if verbosity == 'debug':
                print('QR| (householder) matrix before operation:\n{}'.format(mat_op_on))
            mat_op_on = mlib.dot(op, mat_op_on)
            if verbosity == 'debug':
                print('QR| (householder) matrix after operation:\n{}\n OPERATION:\n{}'.format(mat_op_on, op))
            mat_op_on_T = mlib.transpose(mat_op_on)
            op_log.append(op)

        Q = mlib.eye(nline)
        for iop in op_log:

            Q = mlib.dot(iop, Q)
        
        Q = mlib.transpose(Q)
        return [mat_in, Q, mat_op_on]

    elif mode == 'schmidt':

        #print('QR| ***warning*** There seems one bug that has not been discovered yet, although Gram-Schmidt works well.')
        [_, Q] = gs(mat_in, mode = 'column', verbosity = verbosity)

        Q_t = mlib.transpose(Q)
        R = mlib.dot(Q_t, mat_op_on)

        return [mat_in, Q, R]
