def find_nearest(xdata, x, method = 'tranverse'):

    dist = abs(xdata[0] - x)
    nearest = 0

    lx = len(xdata)

    if method == 'tranverse':

        for i in range(lx):

            temp_dist = abs(xdata[i] - x)
            if temp_dist < dist:

                nearest = i
                dist = temp_dist
    elif method == 'bisection':

        bd1 = 0
        bd2 = lx - 1
        if x > xdata[bd2]:

            return bd2
        elif x < xdata[bd1]:

            return bd1
        else:

            while abs(bd1 - bd2) > 1:

                temp = bd1
                bd1 = round((bd1 + bd2)/2)
                if (xdata[bd1] - x)*(xdata[bd2] - x) > 0:

                    bd2 = temp
            nearest = min(bd1, bd2)
    return nearest