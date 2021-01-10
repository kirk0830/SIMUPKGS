import _mat_lib as mlib
from jacobi_diag import op_gen
from copy import deepcopy
from math import sqrt


def tri_diag(mat_in, verbosity='silent'):
    """
    Tri-diagonalization (Givens method)\n
    # input requirement\n
    mat_in: symmetrical, real matrix, FLOAT data type\n
    # output description\n
    [tdiag, U]\n
    tdiag: tri-diagonalized matrix\n
    U: unitary operator (accumulative Givens operator, in present context)\n
    # formulation\n
    see LATEX formatted text at the end of present function and\n
    mat_in = U·tdiag·U', where U' denotes transpose of unitary operator
    """
    nline = len(mat_in)
    mat_op_on = deepcopy(mat_in)

    op_accum = mlib.eye(nline)
    for irow in range(1, nline-1):
        for icol in range(irow + 1, nline):

            s = mat_op_on[irow-1][icol]/sqrt(
                mat_op_on[irow-1][irow]**2 + mat_op_on[irow-1][icol]**2
            )
            if s == 0:
                c = 1.0
            else:
                c = -(mat_op_on[irow-1][irow]/mat_op_on[irow-1][icol])*s

            P = [[c, s], [-s, c]]
            op = op_gen(size=nline, U=P, irow=icol, icol=irow)

            op_accum = mlib.dot(op_accum, op)
            mat_op_on = mlib.unitary_transform(U=op, mat=mat_op_on)

            if verbosity == 'debug':

                print('\nGivens tri-diagonalization comprehensive report\n'
                      + '-'*50+'\nGivens operator (2x2) P({}, {}) print:'.format(irow, icol))
                mlib.matrix_print(P, decimal=4)
                print('Embedded Givens operator print:')
                mlib.matrix_print(op, decimal=4)
                print('Matrix print:')
                mlib.matrix_print(mat_op_on, decimal=4)

    return [mat_op_on, op_accum]


# LATEX FORMULATION starts from this line
# \left( \begin{matrix}
# 	1&		0&		0&		0\\
# 	0&		c&		0&		-s\\
# 	0&		0&		1&		0\\
# 	0&		s&		0&		c\\
# \end{matrix} \right) \left( \begin{matrix}
# 	a_{11}&		a_{12}&		a_{13}&		a_{14}\\
# 	a_{21}&		a_{22}&		a_{23}&		a_{24}\\
# 	a_{31}&		a_{32}&		a_{33}&		a_{34}\\
# 	a_{41}&		a_{42}&		a_{43}&		a_{44}\\
# \end{matrix} \right) \left( \begin{matrix}
# 	1&		0&		0&		0\\
# 	0&		c&		0&		s\\
# 	0&		0&		1&		0\\
# 	0&		-s&		0&		c\\
# \end{matrix} \right)

# =\left( \begin{matrix}
# 	a_{11}&		a_{12}&		a_{13}&		a_{14}\\
# 	ca_{21}-sa_{41}&		ca_{22}-sa_{42}&		ca_{23}-sa_{43}&		ca_{24}-sa_{44}\\
# 	a_{31}&		a_{32}&		a_{33}&		a_{34}\\
# 	sa_{21}+ca_{41}&		sa_{22}+ca_{42}&		sa_{23}+ca_{43}&		sa_{24}+ca_{44}\\
# \end{matrix} \right) \left( \begin{matrix}
# 	1&		0&		0&		0\\
# 	0&		c&		0&		s\\
# 	0&		0&		1&		0\\
# 	0&		-s&		0&		c\\
# \end{matrix} \right)

# =\left( \begin{matrix}
# 	a_{11}&		ca_{12}-sa_{14}&		a_{13}&		sa_{12}+ca_{14}\\
# 	ca_{21}-sa_{41}&		c^2a_{22}-sca_{42}-sca_{24}+s^2a_{44}&		ca_{23}-sa_{43}&		sca_{22}-s^2a_{42}+c^2a_{24}-sca_{44}\\
# 	a_{31}&		ca_{32}-sa_{34}&		a_{33}&		sa_{32}+ca_{34}\\
# 	sa_{21}+ca_{41}&		sca_{22}+c^2a_{42}-s^2a_{24}-sca_{44}&		sa_{23}+ca_{43}&		s^2a_{22}+sca_{42}+sca_{24}+c^2a_{44}\\
# \end{matrix} \right)

# irow:\ 4\ icol:\ 2,\ let\ \left( 4,1 \right) \ and\ \left( 1,4 \right) \ zero

# sa_{21}+ca_{41}=0
# sa_{12}+ca_{14}=0

# c=-\frac{a_{21}}{a_{41}}s
# c^2+s^2=1

# \frac{a_{21}^{2}+a_{41}^{2}}{a_{41}^{2}}s^2=1

# s=\frac{a_{41}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}
# c=-\frac{a_{21}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}

# \left( \begin{matrix}
# 	1&		0&		0&		0\\
# 	0&		-\frac{a_{21}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}&		0&		-\frac{a_{41}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}\\
# 	0&		0&		1&		0\\
# 	0&		\frac{a_{41}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}&		0&		-\frac{a_{21}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}\\
# \end{matrix} \right) \left( \begin{matrix}
# 	a_{11}&		a_{12}&		a_{13}&		a_{14}\\
# 	a_{21}&		a_{22}&		a_{23}&		a_{24}\\
# 	a_{31}&		a_{32}&		a_{33}&		a_{34}\\
# 	a_{41}&		a_{42}&		a_{43}&		a_{44}\\
# \end{matrix} \right) \left( \begin{matrix}
# 	1&		0&		0&		0\\
# 	0&		-\frac{a_{21}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}&		0&		\frac{a_{41}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}\\
# 	0&		0&		1&		0\\
# 	0&		-\frac{a_{41}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}&		0&		-\frac{a_{21}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}\\
# \end{matrix} \right)

# \mathbf{P}_{24}:

# s=\frac{a_{41}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}
# c=-\frac{a_{21}}{\sqrt{a_{21}^{2}+a_{41}^{2}}}

# \mathbf{P}_{ij}:

# s=\frac{a_{j,i-1}}{\sqrt{a_{i,i-1}^{2}+a_{j,i-1}^{2}}}
# c=-\frac{a_{i,i-1}}{\sqrt{a_{i,i-1}^{2}+a_{j,i-1}^{2}}}
