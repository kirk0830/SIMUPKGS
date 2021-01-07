# homemade
from copy import deepcopy

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
        zero_mat.append(deepcopy(line))
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
        ones_mat.append(deepcopy(line))
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
        eye_mat.append(deepcopy(line))
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
