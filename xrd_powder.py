import math
import numpy as np

def pxrd(xrd2d):

    width = np.shape(xrd2d)[0]
    pxrdout=np.zeros(math.ceil(width*math.sqrt(2)))
    for ix in range(width):
        for iy in range(width):
            pxrdout[math.floor(math.sqrt(ix**2+iy**2))] += xrd2d[ix][iy]

    return pxrdout
