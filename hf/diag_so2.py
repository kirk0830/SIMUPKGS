import math

def diag(mat):

    """
    Diagonalization function in 2-dimension\n
    # input requirement\n
    mat: input matrix in 2-dimension, shape = (2, 2)\n
    ***WARNING*** check if mat is in float data type!\n
    # output description\n
    [diag, U]\n
    diag: diagonalized matrix\n
    U: unitary operator\n
    # Formula\n
    mat = U·diag·U', where U' means transpose of U
    """

    if abs(mat[0][0]-mat[1][1]) > 1E-5:
        theta = 0.5*math.atan(2*mat[0][1]/(mat[0][0]-mat[1][1]))
    else:
        theta = math.pi/4

    eigenVal1 = mat[0][0]*(math.cos(theta)**2) + mat[1][1]*(math.sin(theta)**2) + mat[0][1]*math.sin(2*theta)
    eigenVal2 = mat[0][0]*(math.sin(theta)**2) + mat[1][1]*(math.cos(theta)**2) - mat[0][1]*math.sin(2*theta)
    U = [[math.cos(theta), math.sin(theta)], [math.sin(theta), -math.cos(theta)]]
    diagMat = [[eigenVal1, 0.0], [0.0, eigenVal2]]

    return [diagMat, U]
