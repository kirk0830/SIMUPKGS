# Numerical recipe, Chapter 2-1
from _mat_lib import det, dot, eye
from copy import deepcopy

def gauss_jordan(mat_in, mode = 'inversion', b = [], verbosity = 'silent'):

    """
    Gauss-Jordan inversion\n
    # input requirement\n
    mat_in: squared matrix in 2-dimension\n
    mode: 'inversion' or 'backsubstitution', for the former, will return an inversed mat_in
    while for the latter, will return semi-inversed mat_in, transformed b-vector\n
    b: if mode == 'backsubstitution', this keyword must be specified explicitly. However, if mode == 'inversion'
    and b is given correctly, equation Ax = b will be solved by x = A-1b, x will also be returned.\n
    # output description\n
    Has been introduced in section "input description"
    """

    mat_op_on = deepcopy(mat_in)

    if det(mat_in = mat_op_on):
        # main text of this function
        if verbosity == 'debug':
            print('GAU-JOR INV| non-singluar matrix, safe to calculate inverse...')
        nline = len(mat_op_on)
        record = eye(n = nline)

        zero_diag = 0
        for iline in range(nline):
            if mat_op_on[iline][iline] == 0:
                zero_diag += 1
        if zero_diag > 0:
            print('GAU-JOR INV| zero diagonal element(s) detected: {}\nPartial pivot method will be used...'.format(zero_diag))
            # pivot
        for irow in range(nline):
            # OPERATION 1: PIVOT SWAP
            if mat_op_on[irow][irow] == 0:
                print('GAU-JOR INV| zero diagonal element encounted at row {}, try to swap with other row...'.format(irow))
                pivot = 0
                row2swap = irow
                while pivot == 0 and row2swap <= nline:
                    if mat_op_on[row2swap][irow] != 0:
                        pivot = mat_op_on[row2swap][irow]
                        if verbosity == 'debug':
                            print('GAU-JOR INV| [OPERATION 1] exchange row {} with row {}'.format(row2swap, irow))
                        temp = mat_op_on[row2swap][:]
                        mat_op_on[row2swap][:] = mat_op_on[irow][:]
                        mat_op_on[irow][:] = temp
                        # same exchange operates on element matrix...
                        temp = record[row2swap][:]
                        record[row2swap][:] = record[irow][:]
                        record[irow][:] = temp
                        if mode == 'backsubstitution':
                            temp = b[row2swap]
                            b[row2swap] = b[irow]
                            b[irow] = temp
                    else:
                        row2swap += 1
            else:
                # OPERATION 2: NORMALIZE DIAGONAL ELEMENT
                diag = mat_op_on[irow][irow]
                for icol in range(nline):
                    if verbosity == 'debug':
                        print('GAU-JOR INV| [OPERATION 2] STATUS: normalize row {} with factor {}.'.format(irow, diag))
                        print('GAU-JOR INV| [OPERATION 2] row {}: {}'.format(irow, mat_op_on[irow][:]))
                    mat_op_on[irow][icol] /= diag
                    if verbosity == 'debug':
                        print('GAU-JOR INV| [OPERATION 2] RESULTANT row: {}'.format(mat_op_on[irow][:]))
                    record[irow][icol] /= diag
                if mode == 'backsubstitution':
                    b[irow] /= diag

                # OPERATION 3: ELIMINATE ALL OFF-DIAGONAL ELEMENT TO 0
                if mode == 'backsubstitution':
                    line2sub_start = irow
                else:
                    line2sub_start = 0

                for irow2sub in range(line2sub_start, nline):
                    if irow2sub == irow:
                        continue
                    else:

                        factor = mat_op_on[irow2sub][irow]
                        if verbosity == 'debug':
                            print('GAU-JOR INV| [OPERATION 3] STATUS: subtracting row {} from row {}'.format(mat_op_on[irow][:], mat_op_on[irow2sub][:]))
                            print('GAU-JOR INV| [OPERATION 3] row number: {}, {}'.format(irow, irow2sub))
                            print('GAU-JOR INV| [OPERATION 3] FACTOR = {}'.format(factor))
                        for icol in range(nline):

                            mat_op_on[irow2sub][icol] -= mat_op_on[irow][icol]*factor
                            record[irow2sub][icol] -= record[irow][icol]*factor
                        if mode == 'backsubstitution':
                            b[irow2sub] -= b[irow]*factor
                        if verbosity == 'debug':
                            print('GAU-JOR INV| [OPERATION 3] RESULTANT ROW: {}'.format(mat_op_on[irow2sub][:]))
                            print('GAU-JOR INV| [OPERATION 3] RESULTANT MATRIX: {}'.format(mat_op_on))
        if mode == 'inversion':
            if len(b) != nline:
                return record
            else:
                b_in_2d = [[b[i]] for i in range(nline)]
                x = dot(record, b_in_2d)
                x_in_1d = [x[i][0] for i in range(nline)]
                return [record, x_in_1d]
        elif mode == 'backsubstitution':

            return [mat_op_on, b]
    else:
        print('Singular matrix! Not suitable for Gauss-Jordan method to find its inverse, quit.')
        exit()
