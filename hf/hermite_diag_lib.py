# only for diagonalizing Hermite matrix

from householder import householder as hh
from jacobi_diag import jacobi_diag as ja
import _mat_lib as mlib

from qr_decomp import qr

def hermite_diag(hmat_in, mode = 'householder', eigval_format = 'mat', conv_level = 3, max_iter = 50, verbosity = 'silent'):

    """
    Hermite matrix diagonalization\n
    # input requirement\n
    hmat_in: Hermite matrix, must be REAL, SYMMETRY\n
    mode: 'householder' (recommended) or 'jacobi'\n
    eigval_format: 'mat' or 'list', for 'mat', will return a diagonalized matrix, for 'list', will return list of diagonal\n
    conv_level: (only valid if mode set to 'jacobi') convergence threshold, 10**(-conv_level)\n
    max_iter: (only valid if mode set to 'jacobi') maximum number of iteration\n
    # output description\n
    [diagH, U]\n
    diagH: diagonalized Hermite matrix\n
    U: unitary matrix, also the set of column eigenvectors: U = [|1>|2>|3>...]
    """
    diagH = 0
    U = 0

    if mode == 'householder':

        [_, Q, R] = qr(
            mat_in = hmat_in, 
            mode = 'householder',
            verbosity = verbosity
            )
        # A = Q·R
        # Q'·A = R
        # Q'·A·Q = R·Q = diag(A)
        diag_hmat = mlib.dot(R, Q)
        
        U = Q
        diagH = diag_hmat

    elif mode == 'jacobi':

        [diag_hmat, U] = ja(
            mat_in = hmat_in,
            tri_diag = False,
            conv_level = conv_level,
            max_iter = max_iter,
            verbosity = verbosity
        )

        diagH = diag_hmat
    else:

        exit()
    
    if eigval_format == 'mat':

        return [diagH, U]
    elif eigval_format == 'list':

        nline = len(hmat_in)
        eigval_list = [diagH[iline][iline] for iline in range(nline)]
        return [eigval_list, U]
