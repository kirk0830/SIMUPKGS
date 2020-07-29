import numpy as np

def gau1d(dataInp, sigma, fast = False, fast_cutoff = 0):

    shapeInp = np.shape(dataInp)
    dataOut = np.zeros_like(dataInp)

    print('SMOOTH| 1d-Gaussian function type smooth is activated.')
    if (fast_cutoff - np.floor(fast_cutoff)):
        print('SMOOTH! ***WARNING*** non-integar cutoff value detected, an non-integar cutoff\n'
             +'                     is not meaningful. I will transfer it using floor method.')
        fast_cutoff = np.floor(fast_cutoff)

    if fast and (int(fast_cutoff)>0):
        print('SMOOTH| fast smooth method is activated, cutoff = '+str(fast_cutoff))

        for ix in range(shapeInp[0]):
            for ix_ix in range(max(0, ix-fast_cutoff), min(shapeInp[0], ix+fast_cutoff)):
                dataOut[ix] += dataInp[ix_ix] * np.e ** (-(ix-ix_ix)**2/2/sigma**2)
        
    else:

        for ix in range(shapeInp[0]):
            for ix_ix in range(shapeInp[0]):
                dataOut[ix] += dataInp[ix_ix] * np.e ** (-(ix-ix_ix)**2/2/sigma**2)

    return dataOut

def gau2d(dataInp, sigma, fast = False, fast_cutoff = 0):

    shapeInp = np.shape(dataInp)
    dataOut = np.zeros_like(dataInp)

    print('SMOOTH| 2d-Gaussian function type smooth is activated.')
    if (fast_cutoff - np.floor(fast_cutoff)):
        print('SMOOTH! ***WARNING*** non-integar cutoff value detected, an non-integar cutoff\n'
             +'                     is not meaningful. I will transfer it using floor function.')
        fast_cutoff = np.floor(fast_cutoff)

    if fast and (int(fast_cutoff)>0):
        print('SMOOTH| fast smooth method is activated, cutoff = '+str(fast_cutoff))

        for ix in range(shapeInp[0]):
            for iy in range(shapeInp[1]):
                for ix_ix in range(max(0, ix-fast_cutoff), min(shapeInp[0], ix+fast_cutoff)):
                    for iy_iy in range(max(0, iy-fast_cutoff), min(shapeInp[1], iy+fast_cutoff)):
                        dataOut[ix][iy] += dataInp[ix_ix][iy_iy] * np.e ** (-((ix-ix_ix)**2+(iy-iy_iy)**2)/2/sigma**2)
        
    else:

        for ix in range(shapeInp[0]):
            for iy in range(shapeInp[1]):
                for ix_ix in range(shapeInp[0]):
                    for iy_iy in range(shapeInp[1]):
                        dataOut[ix][iy] += dataInp[ix_ix][iy_iy] * np.e ** (-((ix-ix_ix)**2+(iy-iy_iy)**2)/2/sigma**2)

    return dataOut

# I will add more smooth method soon.