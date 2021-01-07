import _mat_lib as mlib
from copy import deepcopy
from math import sqrt

def cholesky(mat_in, print_transpose = False, iexpr = False):

    """
    Cholesky decomposition\n
    # input requirement\n
    mat_in: input matrix, MUST BE: 1. POSITIVE DEFINITE, 2. SYMMETRY\n
    print_transpose: for Cholesky decomposition is to decompose positive-definite symmetry matrix into
    one lower-triangonal matrix and its transpose, i. e., A = LL', if set this flag to True, L' will 
    also be returned\n
    iexpr: if expressions of elements will be output by stdout\n
    # output description\n
    if print_transpose: [original matrix, L, L']\n
    if !print_transpose: [original matrix, L]
    """
    if mlib.det(mat_in = mat_in) == 0:

        print('***error*** Singular matrix! Cholesky decomposition doesn\'t support present matrix!')
        exit()

    mat_in_bak = deepcopy(mat_in)
    nline = len(mat_in)

    L = mlib.zeros(nline)
    U = mlib.zeros(nline)

    for i in range(nline):

        # diagonal element first to calculate
        expr_str_diag = 'L({}, {}) = sqrt[A({}, {})'.format(i, i, i, i)
        term_Aii = 0
        for j in range(i):

            term_Aii += L[i][j]**2
            word = ' - L({}, {})^2'.format(i, j)
            expr_str_diag += word

        if iexpr:
            print(expr_str_diag+']')
        mat_in[i][i] -= term_Aii
        L[i][i] = sqrt(mat_in[i][i])
        U[i][i] = L[i][i]

        for j in range(i+1, nline):

            expr_str_off_diag = 'L({}, {}) = [A({}, {})'.format(j, i, i, j)
            term_Aij = 0
            for k in range(i):

                word = ' - L({}, {})*L({}, {})'.format(k, i, k, j)
                expr_str_off_diag += word
                term_Aij += L[k][i]*L[k][j]
            
            if iexpr:
                print(expr_str_off_diag+']/L({}, {})'.format(i, i))

            mat_in[i][j] -= term_Aij
            L[j][i] = mat_in[i][j]/L[i][i]
            U[i][j] = L[j][i]

    if print_transpose:

        return [mat_in_bak, L, U]
    else:

        return [mat_in_bak, L]
