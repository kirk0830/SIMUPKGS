from intp_findNearest import find_nearest as fn
import numpy as np

def spline(xdata, ydata, x, calc_mode = 'line', left_bound_D2y = 0, right_bound_D2y = 0):

    lx = len(xdata)
    ly = len(ydata)

    if lx != ly:
        raise TypeError

    if calc_mode == 'line':

        D2y_coeff_mat = np.zeros(shape = (lx, lx))
        D2y_b = np.zeros(shape = (1, lx))

        for i in range(1, lx - 1):

            D2y_coeff_mat[i][i-1] = (xdata[i] - xdata[i-1])/6
            D2y_coeff_mat[i][i] = (xdata[i+1] - xdata[i-1])/3
            D2y_coeff_mat[i][i+1] = (xdata[i+1] - xdata[i])/6
            D2y_b[i] = (ydata[i+1] - ydata[i])/(xdata[i+1] - xdata[i]) - (ydata[i] - ydata[i-1])/(xdata[i] - xdata[i-1])
        
        D2y_coeff_mat[0][0] = 1
        D2y_coeff_mat[-1][-1] = 1
        D2y_b[0] = left_bound_D2y
        D2y_b[-1] = right_bound_D2y

        D2y = np.linalg.solve(D2y_coeff_mat, D2y_b)

        idx = fn(xdata = xdata, x = x, method = 'bisection')

        A = (xdata[idx+1] - x)/(xdata[idx+1] - xdata[idx])
        B = 1 - A
        C = 1/6*(A**3 - A)*(xdata[idx+1]-xdata[idx])**2
        D = 1/6*(B**3 - B)*(xdata[idx+1]-xdata[idx])**2

        y = A*ydata[idx] + B*ydata[idx+1] + C*D2y[idx] + D*D2y[idx+1]

        return [D2y, y]
    
    elif calc_mode == 'point':

        pass
    
