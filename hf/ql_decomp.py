# QL-decomposition, L is left/lower matrix

from householder import householder as hh
import _mat_lib as mlib
from copy import deepcopy

def ql(mat_in, verbosity = 'silent'):

    # use householder by default
    mat_op_on = deepcopy(mat_in)
    
    nline = len(mat_op_on)
    Qt = mlib.eye(nline)

    for iline in range(nline):

        mat_op_on_T = mlib.transpose(mat_op_on)
        vec_in = mat_op_on_T[-iline-1][:(nline-iline)]
        norm_vec = mlib.mod_of_vec(vec_in)
        vec_desti = mlib.zeros(n = 1, m = nline-iline)[0][:]
        vec_desti[-1] = norm_vec

        reflect_op = hh(vec_in, vec_desti, mode = 'L', verbosity = verbosity)
        identi_op = mlib.eye(iline)
        op = mlib.combine_block(reflect_op, identi_op)

        if verbosity == 'debug':
            print('QL| vector read-in: {}\nQL| vector reflected to: {}\n'.format(vec_in, vec_desti))
            print('QL| Reflection operator:')
            mlib.matrix_print(reflect_op, decimal = 4)
            print('QL| Integrated Householder operator:')
            mlib.matrix_print(op, decimal = 4)
            print('QL| Present matrix before operation:')
            mlib.matrix_print(mat_op_on, decimal = 4)

        Qt = mlib.dot(op, Qt)

        mat_op_on = mlib.dot(op, mat_op_on)

        if verbosity == 'debug':
            print('Present matrix after operation:')
            mlib.matrix_print(mat_op_on, decimal = 4)

    Q = mlib.transpose(Qt)

    return [Q, mat_op_on]
