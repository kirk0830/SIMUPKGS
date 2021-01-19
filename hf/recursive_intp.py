import _mat_lib as mlib
import random
from intp_findNearest import find_nearest as fn

def calc_C(
           order, 
           idx, 
           xdata, 
           ydata, 
           x, 
           C_mat, 
           D_mat, 
           fit_mode,
           ini_key
           ):

    if order == 0:

        return ydata[idx]
    elif C_mat[order][idx] != ini_key:

        return C_mat[order][idx]
    else:

        iC = calc_C(
                    order = order - 1, 
                    idx = idx + 1, 
                    xdata = xdata, 
                    ydata = ydata, 
                    x = x,
                    C_mat = C_mat,
                    D_mat = D_mat,
                    fit_mode = fit_mode,
                    ini_key = ini_key
                    )
        iD = calc_D(
                    order = order - 1, 
                    idx = idx, 
                    xdata = xdata, 
                    ydata = ydata, 
                    x = x,
                    C_mat = C_mat,
                    D_mat = D_mat,
                    fit_mode = fit_mode,
                    ini_key = ini_key
                    )

        if fit_mode == 'poly':

            xfrac = (xdata[idx] - x)/(xdata[idx] - xdata[idx + order])
            C_result = xfrac*(iC - iD)
        elif fit_mode == 'rational':

            xfrac = (x - xdata[idx])/(x - xdata[idx + order])

            try:
                C_result = xfrac*iD*(iC - iD)/(xfrac*iD - iC)
            except ZeroDivisionError:
                print('INTEP| ***WARNING*** ZeroDivisonError raised! Status:')
                print('iC = {}\niD = {}\nxfrac = {}\nAccording to Hopital rule, set result as 0'.format(iC, iD, xfrac))
                C_result = 0

        C_mat[order][idx] = C_result

        return C_result

def calc_D(
           order, 
           idx, 
           xdata, 
           ydata, 
           x,
           C_mat,
           D_mat,
           fit_mode,
           ini_key
           ):

    if order == 0:

        return ydata[idx]
    elif D_mat[order][idx] != ini_key:

        return D_mat[order][idx]
    else:

        iC = calc_C(
                    order = order - 1, 
                    idx = idx + 1, 
                    xdata = xdata, 
                    ydata = ydata, 
                    x = x,
                    C_mat = C_mat,
                    D_mat = D_mat,
                    fit_mode = fit_mode,
                    ini_key = ini_key
                    )
        iD = calc_D(
                    order = order - 1, 
                    idx = idx, 
                    xdata = xdata, 
                    ydata = ydata, 
                    x = x,
                    C_mat = C_mat,
                    D_mat = D_mat,
                    fit_mode = fit_mode,
                    ini_key = ini_key
                    )

        if fit_mode == 'poly':

            xfrac = (xdata[idx + order] - x)/(xdata[idx] - xdata[idx + order])
            D_result = xfrac*(iC - iD)
        elif fit_mode == 'rational':

            xfrac = (x - xdata[idx])/(x - xdata[idx + order])
            try:
                D_result = iC*(iC - iD)/(xfrac*iD - iC)
            except ZeroDivisionError:
                print('INTEP| ***WARNING*** ZeroDivisonError raised! Status:')
                print('iC = {}\niD = {}\nxfrac = {}\nAccording to Hopital rule, set result as 0'.format(iC, iD, xfrac))
                D_result = 0

        D_mat[order][idx] = D_result

        return D_result

def recur_intp(
               xdata, 
               ydata, 
               x, 
               intp_order = -1,
               fit_mode = 'poly', 
               calc_mode = 'point',
               start_from = 'nearest'
               ):

    """
    Intepoltation, recursive method\n
    # input requirement\n
    xdata: 1d array\n
    ydata: 1d array, should has the same length with that of xdata\n
    x: the point whose y is demanded\n
    intp_order: order of intepoltation, the higher order is, the more accurate result will get, but it is not
    always needed to use FULL order. If order of polynomial is known, just type-in, this will save plenty of
    time\n
    fit_mode: two kinds of functions are supported, 'poly' or 'rational'. For 'poly' mode, it is named Neville's
    Algorithm, for 'rational', it is named Bulirsch-Stoer Algorithm\n
    calc_mode: 'point' or 'line'. If 'point' mode is set, not all correction terms will be calculated, but only
    calculate terms will be of use. For 'line', all correction terms will be calculated and returned with a single
    y-term demanded\n
    start_from: 'nearest' or something else. For 'nearest', will start from the nearest x point that provided in
    xdata, otherwise will start from the first element of xdata\n
    # output description\n
    y: scalar, for calc_mode == 'point'\n
    [C, D, y]: matrix, matrix and scalar, for calc_mode == 'line'
    """

    lx = len(xdata)
    ly = len(ydata)
    rndkey = random.random()

    C = mlib.ones(n = lx, m = ly, amplify = rndkey)
    D = mlib.ones(n = lx, m = ly, amplify = rndkey)

    if lx != ly:

        raise TypeError

    if intp_order < 0:

        intp_order = lx

    if start_from == 'nearest':

        idx = fn(xdata = xdata, x = x, method = 'bisection')
    else:

        idx = 0

    if x == xdata[idx]:

        return ydata[idx]

    y = 0
    if calc_mode == 'point':

        for iorder in range(intp_order):

            if idx > 0:
                y += calc_D(
                    order = iorder,
                    idx = idx,
                    xdata = xdata,
                    ydata = ydata,
                    x = x,
                    C_mat = C,
                    D_mat = D,
                    fit_mode = fit_mode,
                    ini_key = rndkey
                )
                idx -= 1

            else:
                y += calc_C(
                    order = iorder,
                    idx = idx,
                    xdata = xdata,
                    ydata = ydata,
                    x = x,
                    C_mat = C,
                    D_mat = D,
                    fit_mode = fit_mode,
                    ini_key = rndkey
                )

        return y

    elif calc_mode == 'line':

        pass
