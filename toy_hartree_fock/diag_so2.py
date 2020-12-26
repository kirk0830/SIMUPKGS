import math

def diag(mat):

    if abs(mat[0][0]-mat[1][1]) > 1E-5:
        theta = 0.5*math.atan(2*mat[0][1]/(mat[0][0]-mat[1][1]))
    else:
        theta = math.pi/4

    eigenVal1 = mat[0][0]*(math.cos(theta)**2) + mat[1][1]*(math.sin(theta)**2) + mat[0][1]*math.sin(2*theta)
    eigenVal2 = mat[0][0]*(math.sin(theta)**2) + mat[1][1]*(math.cos(theta)**2) - mat[0][1]*math.sin(2*theta)
    U = [[math.cos(theta), math.sin(theta)], [math.sin(theta), -math.cos(theta)]]
    diagMat = [[eigenVal1, 0.0], [0.0, eigenVal2]]

    return [diagMat, U]
