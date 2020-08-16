import numpy as np

def centercoords(coordsinp, masspower = False, totmass = 1, masslist = []):

    """
    This function is to center all atoms. If you want to center them based on their mass, set masspower as True,
     and total mass in relative atomic mass is required, also, mass list is required.\n
    coordsinp: 2 dimensional list, format like standard xyz file\n
    masspower: bool\n
    totmass: float\n
    masslist: 1 dimensional float list, order should be in accord with coordsinp
    """

    print('COORD| i hope you have convert coordinates to float type, if not, no one knows what will happen...')

    natom = np.shape(coordsinp)[0]
    if masspower:
        if np.shape(masslist)[0] != natom:
            print('COORD| ***error*** number of atoms in atomlist is not consistent with masslist. quit.')
            exit()
        else:
            print('COORD| mass-powered atom centering is activated, useful when calculate rotation interia.')

    x_center = 0
    y_center = 0
    z_center = 0

    coordsout = coordsinp

    for iatom in range(natom):

        mass_power = 1

        if masspower:
            mass_power = masslist[iatom]

        x_center += coordsinp[iatom][-3]/natom * (mass_power/totmass)
        y_center += coordsinp[iatom][-2]/natom * (mass_power/totmass)
        z_center += coordsinp[iatom][-1]/natom * (mass_power/totmass)
    
    for iatom in range(natom):

        coordsout[iatom][-3] -= x_center
        coordsout[iatom][-2] -= y_center
        coordsout[iatom][-1] -= z_center

    return coordsout
