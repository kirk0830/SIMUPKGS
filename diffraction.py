import numpy as np

def diffraction(wavelength, dist, anglemode = False, angleParalist = [], angle = 0):
    # atom type should be recognized before enabling this function.
    # maybe an element-recogizing function is needed.
    amplitude = 1
    if anglemode:
        if angleParalist == []:
            print('DIFFRACTION| ***error*** diffraction paramters are required, however not found. QUIT.')
            exit()
        else:
            pass
            # here amplitude will be recalculated and will not be 1 simply.
            # reference I have forgotten. :(
    else:
        amplitude = 1
        #print('DIFFRACTION| ideal diffraction ~~> not accurate but may be qualitivelly reliable.')
    
    return amplitude*np.cos(2*np.pi*dist/wavelength)
