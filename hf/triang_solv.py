# including both back and forward substitution method to solve linear equations

from copy import deepcopy

def triang_solv(triang_mat, b, mode = 'lower'):

    b_op_on = deepcopy(b)
    back_recog = [
        'upper', 'U', 'back'
    ]
    forward_recog = [
        'lower', 'L', 'forward'
    ]
    for imode in back_recog:
        if imode == mode:
            mode = 'b'
    for imode in forward_recog:
        if imode == mode:
            mode = 'f'
    
    if (mode != 'b') and (mode != 'f'):

        print('***error*** invalid mode selected for triangonal matrix linear equations solving!')
        exit()
    
    x_log = []
    nline = len(triang_mat)
    for irow in range(nline):

        nx = len(x_log)

        if mode == 'b':

            irow *= -1
            irow -= 1
            ia = -nx-1
        else:

            ia = nx
        
        for icol in range(nx):

            ix = icol

            if mode == 'b':

                icol *= -1
                icol -= 1

            b_op_on[irow] -= x_log[ix]*triang_mat[irow][icol]
    
        x = b_op_on[irow]/triang_mat[ia][ia]
        x_log.append(x)
    
    if mode == 'b':

        return x_log[::-1]

    else:
        return x_log
    
