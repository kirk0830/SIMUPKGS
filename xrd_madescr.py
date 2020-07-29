import numpy as np
import math

def xrd_madescr(angle_max, dist, resol, give2d = False):

    scr_edge = dist*math.tan(angle_max*np.pi/180)
    scr_size = math.ceil(scr_edge) * math.floor(10**resol)
    if give2d:
        print('XRD| A '+str(scr_size)+'x'+str(scr_size)+' screen matrix will be created!')

    scr = np.zeros((scr_size, scr_size))
    scr_dl = scr_edge/scr_size
    if give2d:
        print('XRD| pixels coordinates calculated as: \n'
            +'      ix_scr = (index_x-1)*scr_dl, iy_scr = (index_y-1)*scr_dl\n'
            +'      scr_dl = '+str(scr_dl)+' (m)'
        )

    scrcoords = np.zeros((scr_size**2, 3))
    for ipix_x in range(scr_size):
        for ipix_y in range(scr_size):
            scrcoords[ipix_x+scr_size*ipix_y][0] = ipix_x*scr_dl
            scrcoords[ipix_x+scr_size*ipix_y][1] = ipix_y*scr_dl
            scrcoords[ipix_x+scr_size*ipix_y][2] = dist

    if give2d:
        return scr
    else:
        return scrcoords
