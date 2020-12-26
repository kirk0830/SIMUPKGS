import numpy as np
import _in_matrix_op as mop

def _cache_refresh(varlog, var):
    varlog.append(var)
    return varlog[-3:]


def direct_search(
    xfuncy, 
    x_ini, 
    dx_ini, 
    dx_rescale, 
    scaleAlpha = False,
    regionSearch = False, 
    conservatism = True, 
    maxloop = 500, 
    verbosity = 'low'
    ):

    x = x_ini
    dx = dx_ini

    y = xfuncy(x)

    xlog = [x, x, x]
    ylog = [y, y, y]

    x_cache = xlog
    y_cache = ylog
   
    if regionSearch:

        # tranverse all dimensions
        ndim = len(x)
        x_idim = x
        region = []
        for idim in range(ndim):

            # pick up one dimension
            idx = np.zeros_like(dx)
            idx[idim] = dx[idim]
            if scaleAlpha:
                 x_idim = direct_search(
                    xfuncy = xfuncy,
                    x_ini = x_idim,
                    dx_ini = idx,
                    dx_rescale = dx_rescale,
                    scaleAlpha = True,
                    regionSearch = False,
                    conservatism = conservatism,
                    maxloop = maxloop,
                    verbosity = verbosity
                )

            else:
                region_idim = direct_search(
                    xfuncy = xfuncy,
                    x_ini = x,
                    dx_ini = idx,
                    dx_rescale = dx_rescale,
                    scaleAlpha = True,
                    regionSearch = True,
                    conservatism = conservatism,
                    maxloop = maxloop,
                    verbosity = verbosity
                )
                region.append(region_idim)

            if scaleAlpha:
                return x_idim
            else:
                return region

    else:
        iloop = 0
        loopFlag = True
        while loopFlag and iloop <= maxloop:
            # just use line along gradient, not search for a square-shaped region
            x = mop.matrix_plus(x, dx)
            y = xfuncy(x)

            if y > y_cache[-1]:

                if iloop == 0:
                    # if go larger, back and restart
                    dx = np.dot(-1, dx)
                    x = x_cache[-1][:]
                    y = xfuncy(x)
                else:
                    # should output, temporarily omitted here
                    loopFlag = False
                    x_cache = _cache_refresh(xlog, x)
                    y_cache = _cache_refresh(ylog, y)
                    xout = [x_cache[0], x_cache[2]]
                    if scaleAlpha:
                        # recommended for multi-dimensional, also can be used in 1-dimensional problem
                        return np.dot(0.5, mop.matrix_plus(xout[0], xout[1]))
                    else:
                        # doesnt work in multi-dimension, for comparing between n-dimensional is ill-defined
                        # instead, will use norm
                        norm_dx_1 = np.linalg.norm(mop.matrix_minus(x, xout[0]))
                        norm_dx_2 = np.linalg.norm(mop.matrix_minus(x, xout[1]))
                        if norm_dx_1 > norm_dx_2:
                            return [xout[1], xout[0]]
                        else:
                            return [xout[0], xout[1]]
            else:
                # the situation where y_new < y_old
                dx = np.dot(dx_rescale, dx)
                # save...
                x_cache = _cache_refresh(xlog, x)
                y_cache = _cache_refresh(ylog, y)
            iloop += 1
            if conservatism == False:
                iloop = 1
