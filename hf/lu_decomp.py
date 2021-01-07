import _mat_lib as mlib
from copy import deepcopy

def ludecomp(mat_in, pivot = False, iexpr = False):
    """
    Lower triangonal-Upper triangonal matrix decomposition, Crout Algorithm\n
    # input requirement\n
    mat_in: squared matrix\n
    pivot (not implemented): for actual use, there may be zero diagonal element that will cause numerical failure, 
    use pivot to swap rows of original matrix, therefore it is possible when zero diagonal element emerges, 
    program will return a result of LU decomposition of row-swapped original matrix, not that of original matrix. 
    If so, program will pop a warning on this unexpected condition\n
    iexpr: only available if mode == 'recursive', print expression of every element of lower and upper 
    triangonal matrices\n
    # output description\n
    [original matrix, Lower triangonal matrix, Upper triangonal matrix]
    """
    # in principle, should check if singular in advance.

    if mlib.det(mat_in = mat_in) == 0:

        print('***error*** Singular matrix! LU decomposition doesn\'t support present matrix!')
        exit()

    # deepcopy = save as
    mat_in_bak = deepcopy(mat_in)
    nline = len(mat_in)

    L = mlib.eye(nline)
    U = mlib.zeros(nline)

    if iexpr:
        print('-'*80+'\n'+'>> Lower-Upper triangonal matrix decomposition (recursive mode): EXPRESSIONS <<\n'+'-'*80)
    for i in range(nline):
        for j in range(i, nline):

            term_Aij = 0
            expr_str_U = 'U({}, {}) = A({}, {})'.format(i, j, i, j)
            for k in range(i):
                term_Aij += L[i][k]*U[k][j]
                word = ' - L({}, {})*U({}, {})'.format(i, k, k, j)
                expr_str_U += word
            mat_in[i][j] -= term_Aij
            U[i][j] = mat_in[i][j]
            if iexpr:
                print(expr_str_U)

        for j in range(i, nline):
            if j != i:

                term_Aji = 0
                expr_str_L_off_diag = 'L({}, {}) = [A({}, {})'.format(j, i, j, i)
                for k in range(i):
                    term_Aji += L[j][k]*U[k][i]
                    word = ' - L({}, {})*U({}, {})'.format(j, k, k, i)
                    expr_str_L_off_diag += word
                mat_in[j][i] -= term_Aji
                L[j][i] = mat_in[j][i]/U[i][i]
                expr_str_L_off_diag = expr_str_L_off_diag+']/U({}, {})'.format(i, i)
                if iexpr:
                    print(expr_str_L_off_diag)
    if iexpr:
        print('Numerical results can be collected by line: [A0, L, U] = ludecomp(A, mode = \'recursive\', iexpr = True)')

    return [mat_in_bak, L, U]
