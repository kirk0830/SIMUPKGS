import numpy as np

def centercoords(coordsinp):

    print('COORD| i hope you have convert coordinates to float type, if not, no one knows what will happen...')

    natom = np.shape(coordsinp)[0]

    x_center = 0
    y_center = 0
    z_center = 0

    coordsout = coordsinp

    for iatom in range(natom):

        x_center += coordsinp[iatom][-3]/natom
        y_center += coordsinp[iatom][-2]/natom
        z_center += coordsinp[iatom][-1]/natom
    
    for iatom in range(natom):

        coordsout[iatom][-3] -= x_center
        coordsout[iatom][-2] -= y_center
        coordsout[iatom][-1] -= z_center

    return coordsout
