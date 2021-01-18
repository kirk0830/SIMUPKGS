import _mat_lib as mlib


def diamond_intp(x_inp, xdata, ydata, mode = 'polynomial', speedup = False):

    """
    Recursive intepolation\n
    # input requirement\n
    x_inp: x coordinate\n
    xdata: present descrete data of x, in 1d-array\n
    ydata: present descrete data of y, in 1d-array\n
    mode: 'polynomial' or 'rational'\n
    >For 'polynomial', Neville's algorithm is used\n
    >For 'rational', Bulirsch-Stoer algorithm is used\n
    # output description\n
    one-point data of y corresponding to x
    """
    nx = len(xdata)
    ny = len(ydata)

    if nx != ny:
        print('Length of xdata is not consistent with that of ydata!\nLength of xdata: {}, length of ydata: {}'.format(nx, ny))
        raise TypeError
    
    C = mlib.zeros(nx, nx)
    D = mlib.zeros(nx, nx)

    # initialize C and D, left and right corrections

    for ix in range(nx):
        C[0][ix] = ydata[ix]
        D[0][ix] = ydata[ix]
    
    dist = abs(x_inp - xdata[0])
    xs = 0

    # find xs that resides near x:
    for i in range(nx):
        idist = abs(x_inp - xdata[i])
        if idist < dist:
            dist = idist
            xs = i
    
    y = ydata[xs]

    if speedup:

        pass
    else:

        # will calculate the whole C and D matrix...

        for m in range(1, nx):
            for i in range(0, nx-m):

                if mode == 'polynomial':

                    delta = (C[m-1][i+1] - D[m-1][i+1])/(xdata[i] - xdata[i+m])
                    C[m][i] = delta * (xdata[i] - x_inp)
                    D[m][i] = delta * (xdata[i+m] - x_inp)
                elif mode == 'rational':

                    pass

    # add correction
    im = 1

    for icorr in range(nx-1):

        try:
            y += D[im+1][xs-1]
            xs -= 1
        except IndexError:
            y += C[im+1][xs]

        im += 1

    return y
