# Numerical recipe, Chapter 2-1

def zeros(n, m = -1):

    line = []
    zero_mat = []
    if m <= 0:
        for _ in range(n):
            line.append(0)
    else:
        for _ in range(m):
            line.append(0)
    for _ in range(n):
        zero_mat.append(line)
    return zero_mat

def ones(n, m = -1):

    line = []
    ones_mat = []
    if m <= 0:
        for _ in range(n):
            line.append(1)
    else:
        for _ in range(m):
            line.append(1)
    for _ in range(n):
        ones_mat.append(line)
    return ones_mat

def eye(n):

    eye_mat = []

    for irow in range(n):
        line = []
        for icol in range(n):
            if irow == icol:
                line.append(1)
            else:
                line.append(0)
        eye_mat.append(line)
    return eye_mat

def _cut_mat(mat_in, i, j):
    
    nline = len(mat_in)
    mat_out = []
    for irow in range(nline):
        new_row = []
        for icol in range(nline):
            if (irow != i) and (icol != j):
                new_row.append(mat_in[irow][icol])
        if len(new_row):
            mat_out.append(new_row)

    return mat_out

def det(mat_in):

    size = len(mat_in)
    if size == 2:

        val = mat_in[0][0]*mat_in[1][1]-mat_in[0][1]*mat_in[1][0]
        return val
    else:

        expr = 0
        for iline in range(size):

            expr += (-1)**(iline+2)*mat_in[iline][0] * det(_cut_mat(mat_in = mat_in, i = iline, j = 0))
        return expr

def dot(mat_a, mat_b):

    nrow_a = len(mat_a)
    ncol_a = len(mat_a[0][:])
    nrow_b = len(mat_b)
    ncol_b = len(mat_b[0][:])

    if ncol_a == nrow_b:
        mat_out = []
        for irow in range(nrow_a):
            line = []
            for icol in range(ncol_b):
                element = 0
                for i in range(ncol_a):
                    element += mat_a[irow][i]*mat_b[i][icol]
                line.append(element)
            mat_out.append(line)
        
        return mat_out
    else:
        print('ERROR! Number of columns of matrix A does not equal to number of rows of matrix B!')
        raise TypeError

def gaussian_jordan(mat_in, mode = 'inversion', b = [], verbosity = 'silent'):

    """
    Gaussian-Jordan inversion\n
    # input requirement\n
    mat_in: squared matrix in 2-dimension\n
    mode: 'inversion' or 'backsubstitution', for the former, will return an inversed mat_in
    while for the latter, will return semi-inversed mat_in, transformed b-vector\n
    b: if mode == 'backsubstitution', this keyword must be specified explicitly. However, if mode == 'inversion'
    and b is given correctly, equation Ax = b will be solved by x = A-1b, x will also be returned.\n
    # output description\n
    Has been introduced in section "input description"
    """
    if det(mat_in = mat_in):
        # main text of this function
        print('GAU-JOR INV| non-singluar matrix, safe to calculate inverse...')
        nline = len(mat_in)
        record = eye(n = nline)

        zero_diag = 0
        for iline in range(nline):
            if mat_in[iline][iline] == 0:
                zero_diag += 1
        if zero_diag > 0:
            print('GAU-JOR INV| zero diagonal element(s) detected: {}\nPartial pivot method will be used...'.format(zero_diag))
            # pivot
        for irow in range(nline):
            # OPERATION 1: PIVOT SWAP
            if mat_in[irow][irow] == 0:
                print('GAU-JOR INV| zero diagonal element encounted at row {}, try to swap with other row...'.format(irow))
                pivot = 0
                row2swap = irow
                while pivot == 0 and row2swap <= nline:
                    if mat_in[row2swap][irow] != 0:
                        pivot = mat_in[row2swap][irow]
                        if verbosity == 'debug':
                            print('GAU-JOR INV| [OPERATION 1] exchange row {} with row {}'.format(row2swap, irow))
                        temp = mat_in[row2swap][:]
                        mat_in[row2swap][:] = mat_in[irow][:]
                        mat_in[irow][:] = temp
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
                diag = mat_in[irow][irow]
                for icol in range(nline):
                    if verbosity == 'debug':
                        print('GAU-JOR INV| [OPERATION 2] STATUS: normalize row {} with factor {}.'.format(irow, diag))
                        print('GAU-JOR INV| [OPERATION 2] row {}: {}'.format(irow, mat_in[irow][:]))
                    mat_in[irow][icol] /= diag
                    if verbosity == 'debug':
                        print('GAU-JOR INV| [OPERATION 2] RESULTANT row: {}'.format(mat_in[irow][:]))
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

                        factor = mat_in[irow2sub][irow]
                        if verbosity == 'debug':
                            print('GAU-JOR INV| [OPERATION 3] STATUS: subtracting row {} from row {}'.format(mat_in[irow][:], mat_in[irow2sub][:]))
                            print('GAU-JOR INV| [OPERATION 3] row number: {}, {}'.format(irow, irow2sub))
                            print('GAU-JOR INV| [OPERATION 3] FACTOR = {}'.format(factor))
                        for icol in range(nline):

                            mat_in[irow2sub][icol] -= mat_in[irow][icol]*factor
                            record[irow2sub][icol] -= record[irow][icol]*factor
                        if mode == 'backsubstitution':
                            b[irow2sub] -= b[irow]*factor
                        if verbosity == 'debug':
                            print('GAU-JOR INV| [OPERATION 3] RESULTANT ROW: {}'.format(mat_in[irow2sub][:]))
                            print('GAU-JOR INV| [OPERATION 3] RESULTANT MATRIX: {}'.format(mat_in))
        if mode == 'inversion':
            if len(b) != nline:
                return record
            else:
                b_in_2d = [[b[i]] for i in range(nline)]
                x = dot(record, b_in_2d)
                x_in_1d = [x[i][0] for i in range(nline)]
                return [record, x_in_1d]
        elif mode == 'backsubstitution':

            return [mat_in, b]
    else:
        print('Singular matrix! Not suitable for Gaussian-Jordan method to find its inverse, quit.')
        exit()
