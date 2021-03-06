# homemade
from copy import deepcopy

def zeros(n, m = -1):

    line = []
    zero_mat = []
    if m <= 0:
        for _ in range(n):
            line.append(0.0)
    else:
        for _ in range(m):
            line.append(0.0)
    for _ in range(n):
        zero_mat.append(deepcopy(line))
    return zero_mat

def ones(n, m = -1):

    line = []
    ones_mat = []
    if m <= 0:
        for _ in range(n):
            line.append(1.0)
    else:
        for _ in range(m):
            line.append(1.0)
    for _ in range(n):
        ones_mat.append(deepcopy(line))
    return ones_mat

def eye(n, m = -1, amplify = 1.0):

    eye_mat = []

    if (n == m) or (m <= 0):
        m = n
    
    for irow in range(n):
        line = []
        for icol in range(m):
            if irow == icol:
                line.append(1.0 * amplify)
            else:
                line.append(0.0)
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

def dot(mat_a, mat_b, mode = 'matmat'):

    """
    # Usage\n
    mode: 'matmat', 'matket' or 'bramat'\n
    >'matmat': multiply one matrix with the other\n
    >'matket': rotate one ket (mat_b) |> by matrix (mat_a)\n
    >'bramat': rotate one bra (mat_a) <| by matrix (mat_b)
    """
    if mode == 'matmat':

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
            print('MLIB| (DOT) ERROR! Number of columns of matrix A does not equal to number of rows of matrix B!')
            print('Matrix A: ({}, {})'.format(nrow_a, ncol_a))
            print('Matrix B: ({}, {})'.format(nrow_b, ncol_b))
            exit()

    elif mode == 'matket':

        nrow_a = len(mat_a)
        ncol_a = len(mat_a[0][:])
        nrow_b = len(mat_b)

        if ncol_a == nrow_b:

            ket_out = []
            for irow in range(nrow_a):
                component = 0
                for icol in range(nrow_b):
                    component += mat_a[irow][icol]*mat_b[icol][0]
                ket_out.append([component])

            return ket_out
        else:
            print('MLIB| (DOT) ERROR! Number of columns of matrix A does not equal to number of rows of matrix B!')
            print('Matrix A: ({}, {})'.format(nrow_a, ncol_a))
            print('Matrix B: ({}, {})'.format(nrow_b, '1'))
            exit()

    elif mode == 'bramat':

        ncol_a = len(mat_a)
        nrow_b = len(mat_b)
        ncol_b = len(mat_b[0][:])

        if ncol_a == nrow_b:

            bra_out = []
            for icol in range(ncol_a):
                component = 0
                for irow in range(ncol_a):
                    component += mat_a[irow]*mat_b[irow][icol]
                bra_out.append(component)

            return bra_out
        else:
            print('MLIB| (DOT) ERROR! Number of columns of matrix A does not equal to number of rows of matrix B!')
            print('Matrix A: ({}, {})'.format('1', ncol_a))
            print('Matrix B: ({}, {})'.format(nrow_b, ncol_b))
            exit()

def pivot(mat_in, row_num, mode = 'partial'):

    nline = len(mat_in)

    if mode == 'partial':

        temp_row = mat_in[row_num][:]
        irow = row_num
        while mat_in[irow][row_num] == 0:
            irow += 1
            if irow == nline:
                print('***error*** partial pivoting failed due to non-zero element in present column not founded.')
                exit()
        mat_in[row_num][:] = mat_in[irow][:]
        mat_in[irow][:] = temp_row

        return mat_in
    elif mode == 'full':

        pass

from math import sqrt

def mod_of_vec(vec):

    mod = 0
    for icompo in vec:

        mod += icompo**2
    return sqrt(mod)

def braket(bra, ket, mode = '2bra', decimal = False):

    """
    scalar (inner) product of two vectors.
    <|>\n
    bra: <|, row vector, left vector\n
    ket: |>, column vector, right vector\n
    This function is for calculating scalar product between two vectors,
    for calculation (span) of two vectors, use function 'ketbra'.
    """

    prod = 0
    n = len(bra)
    m = len(ket)

    if n != m:
        print('***error*** shape of |> is not consistent with <|, quit.')
        exit()

    if mode == 'braket':

        pass
    elif mode == '2ket':

        pass
    elif mode == '2bra':
        # two 1d list
        for i in range(n):

            prod += bra[i]*ket[i]

    if decimal:
        return round(prod, decimal)
    else:
        return prod
        
def transpose(mat_in):

    nrow = len(mat_in)
    ncol = len(mat_in[-1][:])

    mat_out = zeros(n = ncol, m = nrow)
    for irow in range(nrow):
        for icol in range(ncol):

            mat_out[icol][irow] = mat_in[irow][icol]
    return mat_out

# things will be easier if import numpy but here I wont import it
from _in_matrix_op import matrix_minus as minus
from _in_matrix_op import matrix_plus as plus

def ketbra(ket, bra, mode = '2bra', amplify = 1.0):

    """
    tensor product |><|\n
    mode: 'ketbra', '2ket' or '2bra'
    """
    # standard format...
    # bra: row vector: <|: [b1, b2, b3, ...]
    # ket: col vector: |>: [[k1], [k2], [k3], ...]

    ncol = len(bra)
    nrow = len(ket)

    tensor = zeros(n = len(ket), m = len(bra))

    if mode == 'ketbra':

        for irow in range(nrow):
            for icol in range(ncol):

                tensor[irow][icol] = bra[icol]*ket[irow][-1]*amplify
    elif mode == '2ket':

        for irow in range(nrow):
            for icol in range(ncol):

                tensor[irow][icol] = bra[icol][-1]*ket[irow][-1]*amplify
    elif mode == '2bra':

        for irow in range(nrow):
            for icol in range(ncol):

                tensor[irow][icol] = bra[icol]*ket[irow]*amplify
    
    return tensor

def combine_block(mat1, mat2, mode = 'diag'):

    nline1 = len(mat1)
    nline2 = len(mat2)
    nline = nline1 + nline2

    mat_out = zeros(nline)
    if mode == 'diag':

        for irow in range(nline):
            for icol in range(nline):
                if (irow < nline1) and (icol < nline1):
                    mat_out[irow][icol] = mat1[irow][icol]
                elif (irow >= nline1) and (icol >= nline1):
                    mat_out[irow][icol] = mat2[irow-nline1][icol-nline1]

    return mat_out

def trun_square_mat(mat_in, mode = 'leftup'):

    mat_out = deepcopy(mat_in)

    if mode == 'leftup':
        line = mat_out[0][:]
        mat_out.remove(line)
        mat_out = transpose(mat_out)
        line = mat_out[0][:]
        mat_out.remove(line)
        mat_out = transpose(mat_out)
    
    return mat_out
    
def roundoff(mat_in, decimal = 10):

    mat_out = deepcopy(mat_in)
    nrow = len(mat_out)
    ncol = len(mat_out[-1][:])
    for irow in range(nrow):
        for icol in range(ncol):
            mat_out[irow][icol] = round(mat_out[irow][icol], decimal)
    return mat_out

def matrix_print(mat_in, decimal = False, comma = False):

    nline = len(mat_in)
    aprxed_mat = deepcopy(mat_in)
    if decimal:

        aprxed_mat = roundoff(mat_in = mat_in, decimal = decimal)

    for iline in range(nline):

        if iline == 0:
            print('MATRIX PRINT>>> \n[')
        if comma:
            if iline != nline-1:
                ifcomma = ', '
            else:
                ifcomma = ''
            print(aprxed_mat[iline][:]+ifcomma)
        else:
            print(aprxed_mat[iline][:])
        if iline == nline-1:
            print(']')

def unitary_transform(U, mat):

    mat_out = deepcopy(mat)
    Ut = transpose(U)
    mat_out = dot(Ut, mat_out)
    mat_out = dot(mat_out, U)

    return mat_out

def symm_check(mat):

    nrow = len(mat)
    ncol = len(mat[-1][:])
    if nrow != ncol:

        return False
    for irow in range(nrow):
        for icol in range(irow, ncol):
            if mat[irow][icol] != mat[icol][irow]:
                return False
    
    return True

def diag_to_list(diag_mat_in):

    nrow = len(diag_mat_in)
    eigval_list = []

    for irow in range(nrow):

        eigval_list.append(diag_mat_in[irow][irow])
    
    return eigval_list

def list_diff(list_old, list_new, mode = 'abs'):

    if mode == 'norm2':

        norm_old = mod_of_vec(list_old)
        norm_new = mod_of_vec(list_new)
        return abs(norm_new-norm_old)

    list_diff = minus(list_new, list_old)
    result = 0.0
    for idiff in list_diff:

        if mode == 'abs':
            result += abs(idiff)
        elif mode == 'sqr':
            result += idiff**2
        elif mode == 'sum':
            result += idiff
        elif mode == 'norm1':
            result += idiff**2
        else:
            print('MLIB| (List-diff) present mode \'{}\' selected is not supported yet'.format(mode))
    if mode == 'norm1':
        return sqrt(result)
    else:
        return result

def bra2ket(bra):

    return [[ibra] for ibra in bra]

def ket2bra(ket):

    return [iket[0] for iket in ket]

def normalize(vec, mode = 'bra'):

    if mode == 'bra':

        vec_op_on = deepcopy(vec)
    elif mode == 'ket':

        vec_op_on = vec[:][0]

    norm = mod_of_vec(vec_op_on)

    bralist = []
    for ivec in vec_op_on:

        if (ivec == 0) and (norm == 0):
            bralist.append(1.0)
        else:
            bralist.append(ivec/norm)

    if mode == 'bra':

        return bralist
    elif mode == 'ket':
        
        return(bra2ket(bralist))
