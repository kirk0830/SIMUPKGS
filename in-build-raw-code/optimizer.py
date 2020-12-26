
import numpy as np
import _in_matrix_op as mop

def memoryUpdate(varlog, var):
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

    """
    # simple minima searching algorithm\n
    This function is used to find local minima by simple accurate line searching method\n
    # input requirement\n
    xfuncy: a function that should be defined external and, just type-in function name to use it in this function\n
    x_ini: one point for initial guess of variable\n
    dx_ini: initial guess of step. a small value is recommanded for very first optimization\n
    dx_rescale: if set to number larger than 1, an aggressive search will be activated and each step is twice of the latter one. If set as 1, a normal searching. If set to a number smaller than 1 but still positive, this search will finally reach a fix point.\n
    conservatism: True is default and recommanded, in case of unknown non-convergence\n
    maxloop: only needed when 'conservatism = True', the max step to find an inteval that may contains minima\n
    verbosity: 'low', 'medium' and 'high' are available options\n
    # output description\n
    a list that contains two element, the first is lower boundary, the second is upper boundary
    """

    # the first parameter must be a function that can yield function value from input x_ini
    x = x_ini
    dx = dx_ini
    try:
        y = xfuncy(x)
    except TypeError:
        print('LinearSearch| ***error*** xfuncy parameter has not obtained correct value. quit.')
        exit()
    # variables storage
    xlog = [x, x, x]
    ylog = [y, y, y]
    xmem = xlog
    ymem = ylog
    # loop
    iloop = 0
    loopFlag = True

    if verbosity == 'low':
        print('LinearSearch| verbosity has been set to \'low\', no information will be printed out on screen.')
    while loopFlag and iloop <= maxloop:
        
        x = mop.matrix_plus(x, dx)
        y = xfuncy(x)
        if verbosity == 'high':
            print('LinearSearch| searching point: x = '+str(x)+', y = '+str(y))
        if y > ymem[-1]:
            if iloop == 0:
                # if go larger, back and restart
                if verbosity == 'high':
                    print('LinearSearch| search direction is wrong, redirected.')
                dx = np.dot(-1, dx)
                x = xmem[-1][:]
                y = xfuncy(x)
            else:
                # should output, temporarily omitted here
                loopFlag = False
                xmem = memoryUpdate(xlog, x)
                ymem = memoryUpdate(ylog, y)
                xout = [xmem[0], xmem[2]]
                if verbosity == 'medium' or verbosity == 'high':
                    print('LinearSearch| interval has been found! x_min = '+str(min(xout))+', x_max = '+str(max(xout))+'.')
                return [min(xout), max(xout)]
                
        else:
            # the situation where y_new < y_old
            if verbosity == 'high':
                print('LinearSearch| step '+str(iloop)+', larger step is used to find minimium.')
            dx = np.dot(dx_rescale, dx)
            # save...
            xmem = memoryUpdate(xlog, x)
            ymem = memoryUpdate(ylog, y)
        iloop += 1
        if conservatism == False:
            iloop = 1

def golden_search(
    xfuncy, 
    x_min, 
    x_max, 
    epsilon, 
    conservatism = True, 
    maxloop = 500, 
    otherRatio = False, 
    ratio = 0.618, 
    verbosity = 'low'
    ):

    """
    # Golden search algorithm\n
    a simple method that finds one local minima.\n
    # Input requirement\n
    xfuncy: external function of y(x)\n
    x_min: lower boundary of inteval where to find minima\n
    x_max: upper boundary of inteval where to find minima\n
    epsilon: accurancy control parameter, the larger you set, the more precise you will get\n
    conservatism: see document of function "direct_search" function\n
    maxloop: see document of function "direct_search" function\n
    otherRatio: [bool] although this function is named Golden_search, take in mind that Golden just means 0.618, so other proportion is also supported\n
    ratio: [float] if set otherRatio = True, you can use an arbitrary value to shrink your inteval\n
    # Output description\n
    a list contains both modified lower and upper boundaries
    """
    x_left = x_min
    x_right = x_max
    try:
        y_left = xfuncy(x_left)
        y_right = xfuncy(x_right)
    except TypeError:
        print('LinearSearch| ***error*** xfuncy parameter has not obtained correct value. quit.')
        exit()
    
    if verbosity != 'low':
        print('LinearSearch| golden search method is activated, search domain: ['+str(x_min)+', '+str(x_max)+'].')
    
    x_left_mem = []
    x_right_mem = []
    y_left_mem = []
    y_right_mem = []

    if otherRatio:
        t = ratio
        if verbosity != 'low':
            print('LinearSearch| golden echo: a modified ratio to segment inteval is used. ratio input: '+str(t))
    else:
        t = 0.618
    
    # we use epsilon as the magnitude of inteval, i.e., if epsilon = 1, convergence requirement: x_max - x_min <= 1E-1*(x_max_0 - x_min_0)
    iloop = 0
    while (iloop <= maxloop) and (
        np.linalg.norm(mop.matrix_minus(x_right, x_left))
        >= 
        10**(-epsilon)*np.linalg.norm(mop.matrix_minus(x_max, x_min))
        ):
        iloop += 1

        if y_right >= y_left:
            x_right_mem.append(x_right)
            x_right = mop.matrix_plus(
                x_left, 
                np.dot(max(t, 1-t), mop.matrix_minus(x_right, x_left))
                )
            y_right = xfuncy(x_right)
            y_right_mem.append(y_right)
            if verbosity == 'high':
                print('LinearSearch| shrink inteval upper boundary: '+str(x_right_mem[-1])+' -> '+str(x_right))
        else:
            x_left_mem.append(x_left)
            x_left = mop.matrix_plus(
                x_left, 
                np.dot(min(t, 1-t), mop.matrix_minus(x_right, x_left))
                )
            y_left = xfuncy(x_left)
            y_left_mem.append(y_left)
            if verbosity == 'high':
                print('LinearSearch| shrink inteval lower boundary: '+str(x_left_mem[-1])+' -> '+str(x_left))
        
        if verbosity != 'low':
            print('LinearSearch| present inteval: ['+str(x_left)+', '+str(x_right)+'].')
        if conservatism == False:
            iloop = 1

    return [x_left, x_right]

def quadratic_search(
    xfuncy, 
    x_ini, 
    dx_ini, 
    conservatism = True, 
    maxloop = 500, 
    verbosity = 'low'
    ):
    """
    # Usage Warning\n
    This function is based on consequential 2-spline fitting. So it may always happen if local curvature is negative that this function will find a maximum instead of
    minimum, so this function is not recommended and will only coded in the future. Instead, I will most probably write a ploynomial_search that supports 2- and 3-spline
    both.
    """
    pass

def armijo_fuzzy_search_1d(
    xfuncy, 
    xgrady, 
    x_ini, 
    dx_ini, 
    dxInterface = False,
    sigma = 0.2, 
    beta = 0.5, 
    converLevel = 20, 
    conservatism = True, 
    maxloop = 500, 
    verbosity = 'low'
    ):
    """
    # Original Armijo fuzzy search algorithm\n
    This is original version of Armijo algorithm, where step in it is quite arbitrary and may be unreasonable under some circumstances.\n
    # Formulation (Armijo criterion)\n
    Given parameters beta in (0, 1) and sigma in (0, 0.5), for m belonging to non-negative integars:\n
    if inequality f(x_k + beta**m * dx) <= f(x_k) + sigma * beta**m * Df(x_k) * dx can be satisfied,\n
    use the minimal m as output and update x_k <- x_k + beta**m * dx is admitted.\n
    This criterion is equivalent with:\n
    [f(x_k + beta**m * dx) - f(x_k)]/(beta**m * Df(x_k) * dx) ~=~ Df(x_k) <= sigma*Df(x_k),\n
    that is, ensuring Df(x_k) is negative.\n
    # Input requirement\n
    xfuncy: external function y(x)\n
    xgrady: derivative of y(x)\n
    x_ini: starting point of x\n
    dx_ini: for internal usage, i.e., dxInterface = False, updating step magnitude, dx = 1, stepsize = 1E-1*x_ini;
    for external usage that dx provided by other algorithm, set dxInterface = True and dx is directly given the proper value.\n
    sigma: derivative shrinking coefficent\n
    beta: step shrinking coeffient\n
    converLevel: maximum of order of m to try.
    """

    # Type chcek:
    try:
        y_ini = xfuncy(x_ini)
        yy_ini = xgrady(x_ini)
        if verbosity == 'high':
            print('LinearSearch| initial point information: x = '+str(x_ini)+', y = '+str(y_ini)+', dy/dx = '+str(yy_ini))
    except TypeError:
        print('LinearSearch| ***error*** xfuncy parameter has not obtained correct value. quit.')
        exit()
    
    if dxInterface:
        dx = dx_ini
    else:
        dx = x_ini * 10**(-dx_ini) # will be deprecated in advanced oprimization methods
    
    x = x_ini
    x_mem = [x]
    m_mem = [] 

    iloop = 0
    if_reversed_dx = False
    while iloop <= maxloop:

        m=0

        while xfuncy(x+beta**m*dx) >= (xfuncy(x) + sigma*(beta**m*dx)*xgrady(x)) and m<=converLevel:
            print('LinearSearch| Armijo order m = '+str(m)
                 +', terms of inequality: '+str(xfuncy(x+beta**m*dx))+', '+str(xfuncy(x) + sigma*beta**m*dx*xgrady(x)))
            m += 1

        if m<(converLevel-1):
            m_mem.append(m)
            x += beta**m*dx
            
            if verbosity != 'low':
                print('LinearSearch| x updated: '+str(x_mem[-1])+' -> '+str(x))
            x_mem.append(x)
        elif iloop == 0 and if_reversed_dx == False:
            print('LinearSearch| ***warning*** max Armijo order has been reached at the initial opt-step, searching direction is set\n'
                 +'                  reversed to find if possible to get lower value...')
            dx *= -1
            if_reversed_dx = True
            continue
        else:
            print('LinearSearch| ***warning*** max Armijo order has been reached, Armijo exit.')
            break
        iloop += 1
        if conservatism == False:
            iloop = 1
    
    return x

def armijo_step_revision(
    xfuncy, 
    xgrady, 
    x, 
    dx, 
    sigma = 0.2, 
    beta = 0.5, 
    m_thre = 20,
    verbosity = 'low'):
    # internal function, not be expected to be called directly from external function
    m = 0

    while (m <= m_thre) and (
        xfuncy(
         mop.matrix_plus(x, np.dot(beta**m, dx))
         ) 
         >= 
         (xfuncy(x) + np.dot(sigma*beta**m, np.dot(dx, xgrady(x))))
         ):
        
        if verbosity == 'high' or verbosity == 'debug':
            if verbosity == 'debug':
                print('OPT| Armijo: sigma = '+str(sigma)+', beta = '+str(beta))
            print('OPT| Armijo-asisted convergence method, step shrinking degree m = '+str(m))
        m += 1
    if m <= m_thre:
        return beta**m*dx
    else:
        print('OPT| ***warning*** Armijo method failed to find suitable stepsize, this may happen when\n'
             +'                   optimization task is already converged but inappropriate values are given\n'
             +'                   such as convergence threshold, Armijo beta, Armijo sigma.\n'
             +'     operation --> return stepsize = 0, to expire number of steps of the whole optimization.')
        return 0

def newton_1d(
    xfuncy, 
    x_ini, 
    dx_ini, 
    epsilon,
    conservatism = True,
    maxloop = 50,
    verbosity = 'low'):
    # note: epsilon here measures convergence threshold of variable value
    dx = dx_ini
    x_1 = x_ini
    x_2 = x_1 + dx
    x_3 = x_2 + dx
    # Newton assumes that it is already near minima, so curvature must be positive
    # make fitting quadratic function a*x**2 + b*x + c = f

    A = [
        [x_1**2, x_1, 1],
        [x_2**2, x_2, 1],
        [x_3**2, x_3, 1]
    ]
    B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]
    if verbosity == 'debug':
        print('OPT| on-the-fly info: [type: initial info]\n'
             +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
             +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
             +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
             +'                      forward step: {}'.format(str(dx)))
    # solve for a, b, c
    para = np.linalg.solve(A, B)
    dx = -para[1]/2/para[0] - x_1
    iloop = 0
    while dx > float(epsilon) and iloop < maxloop:

        x_1 = -para[1]/2/para[0]
        x_2 = x_1 + dx
        x_3 = x_2 + dx
        A = [
            [x_1**2, x_1, 1],
            [x_2**2, x_2, 1],
            [x_3**2, x_3, 1]
        ]
        B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]
        if verbosity == 'debug':
            print('OPT| on-the-fly info: [type: step info, step {}]\n'.format(str(iloop))
                +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
                +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
                +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
                +'                      forward step: {}'.format(str(dx)))
        para = np.linalg.solve(A, B)
        dx = -para[1]/2/para[0] - x_1
        iloop += 1

    return x_1

def steep_1d(
    xgrady, 
    x_ini, 
    acc_level, 
    shrink = True, 
    alpha = 0.1, 
    conservatism = True,
    maxloop = 50, 
    verbosity = 'low'):
    # original steep method
    # analytical expressions are needed
    x_mem = [x_ini, x_ini, x_ini]
    if shrink == False:
        alpha = 1
    dx = -alpha * xgrady(x_ini)

    iloop = 0
    while dx > 10**(-acc_level) and iloop < maxloop:
        x = x_mem[-1] + dx
        if verbosity == 'high' or verbosity == 'debug':
            print('OPT| steep method updates variable: '+str(x_mem[-1])+' -> '+str(x))
        x_mem = memoryUpdate(x_mem, x)
        dx = -alpha * xgrady(x)
        iloop += 1
        if conservatism == False:
            iloop = 1
    return x

# if function is multi-dimensional, its output is still 1-dimensional, but its derivative has the same 
# size as varaible, multi-dimensional. So for optimization method that not uses derivative, codes may 
# not need to be changed largely...

# from here, n-dimensional functions are supported

def armijo_steep(
    xfuncy,
    xgrady, 
    x_ini, 
    acc_level, 
    Armijo_sigma = 0.2, 
    Armijo_beta = 0.5, 
    max_Armijo = 20, 
    conservatism = True,
    maxloop = 50, 
    verbosity = 'low'
    ):
    x_mem = [x_ini, x_ini, x_ini]

    dx_0 = -xgrady(x_ini)
    dx = armijo_step_revision(
        xfuncy = xfuncy, 
        xgrady = xgrady,
        x = x_ini,
        dx = dx_0,
        sigma = Armijo_sigma,
        beta = Armijo_beta,
        m_thre = max_Armijo,
        verbosity = verbosity
        )

    iloop = 0
    while np.linalg.norm(dx) > 10**(-acc_level) and iloop < maxloop:
        x = mop.matrix_plus(x_mem[-1], dx)
        if verbosity == 'high' or verbosity == 'debug':
            print('OPT| steep method updates variable: '+str(x_mem[-1])+' -> '+str(x))
        x_mem = memoryUpdate(x_mem, x)
        dx_0 = -xgrady(x)
        dx = armijo_step_revision(
            xfuncy = xfuncy, 
            xgrady = xgrady,
            x = x_ini,
            dx = dx_0,
            sigma = Armijo_sigma,
            beta = Armijo_beta,
            m_thre = max_Armijo,
            verbosity = verbosity
            )
        iloop += 1
        if conservatism == False:
            iloop = 1
    return x

def armijo_steep_newton(
    xfuncy,
    xgrady, 
    x_ini, 
    dx_ini, 
    acc_level, 
    Armijo_sigma = 0.2, 
    Armijo_beta = 0.5, 
    max_Armijo = 20, 
    conservatism = True,
    maxloop = 50, 
    verbosity = 'low'
    ):
    epsilon = 10**(-acc_level)

    ndim = len(x_ini)
    dx = dx_ini
    x_1 = x_ini
    x_2 = mop.matrix_plus(x_1, dx)
    x_3 = mop.matrix_plus(x_2, dx)
    x_1_sqr = mop.matrix_dot_power(x_1, 2)
    x_2_sqr = mop.matrix_dot_power(x_2, 2)
    x_3_sqr = mop.matrix_dot_power(x_3, 2)

    B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]

    if verbosity == 'debug':
        print('OPT| on-the-fly info: [type: initial info]\n'
             +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
             +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
             +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
             +'                      forward step: {}'.format(str(dx)))
    
    para = np.zeros((ndim, 3))

    for i in range(ndim):
        A = [
            [x_1_sqr[i], x_1[i], 1],
            [x_2_sqr[i], x_2[i], 1],
            [x_3_sqr[i], x_3[i], 1]
        ]
        para_this = np.linalg.solve(A, B)
        para[i][:] = para_this
        dx[i] = -para[i][1]/2/para[i][0] - x_1[i]

    iloop = 0

    ifSteep = False
    steeplist = []
    while np.linalg.norm(dx) > float(epsilon) and iloop < maxloop:

        if ifSteep:
            x_1_prev = x_1
            x_1 = np.dot(-0.5, mop.matrix_dot_division(para[:][1], para[:][0]))
            for index in steeplist:
                x_1[index] = x_1_prev[index] + dx[index]
            ifSteep = False
            steeplist = []
        else:
            x_1 = np.dot(-0.5, mop.matrix_dot_division(para[:][1], para[:][0]))

        x_2 = mop.matrix_plus(x_1, dx)
        x_3 = mop.matrix_plus(x_2, dx)
        x_1_sqr = mop.matrix_dot_multiply(x_1, x_1)
        x_2_sqr = mop.matrix_dot_multiply(x_2, x_2)
        x_3_sqr = mop.matrix_dot_multiply(x_3, x_3)

        B = [xfuncy(x_1), xfuncy(x_2), xfuncy(x_3)]
        for i in range(ndim):
            A = [
                [x_1_sqr[i], x_1[i], 1],
                [x_2_sqr[i], x_2[i], 1],
                [x_3_sqr[i], x_3[i], 1]
            ]
            para_this = np.linalg.solve(A, B)
            para[i][:] = para_this
            # save parameters, anyway.
            if para[i][0] >= 0:
                dx[i] = -para[i][1]/2/para[i][0] - x_1[i]
            else:
                print('OPT| Armijo-steep-newton: negative curvature detected, optimization switches to Armijo-steep.')
                dx_0 = -xgrady(x_1[i])
                dx[i] = armijo_step_revision(
                        xfuncy = xfuncy, 
                        xgrady = xgrady,
                        x = x_ini,
                        dx = dx_0,
                        sigma = Armijo_sigma,
                        beta = Armijo_beta,
                        m_thre = max_Armijo,
                        verbosity = verbosity
                        )
                steeplist.append(i)
        if len(steeplist) > 0:
            ifSteep = True

        if verbosity == 'debug':
            print('OPT| on-the-fly info: [type: step info, step {}]\n'.format(str(iloop))
                +'                      point1: {}, value1: {}\n'.format(str(x_1), str(B[0]))
                +'                      point2: {}, value2: {}\n'.format(str(x_2), str(B[1]))
                +'                      point3: {}, value3: {}\n'.format(str(x_3), str(B[2]))
                +'                      forward step: {}'.format(str(dx)))

        iloop += 1

    return x_1

def cg(
    xfuncy, 
    xgrady, 
    x0, 
    acc, 
    betamode = 'FR',
    conservatism = True, 
    maxloop = 50,
    verbosity = 'low', 
    memoryout = False):

# Developer note:
# ---------------
# conjugate gradient method is widely used and accepted as the most rubust optimization algorithm due to it
# must optimize variables and will never cause any raises in other dimension -> it is because between any 
# two steps, optimization goes in orthogonal directions. Also note that CG is based on gradient, so its speed
# may be slower than Newton. In case of oscillation, an accurate linear searching method is recommeded and
# also it is what I implement in CG.

# core of CG includes two aspects,
# 1. orthogonal displacement on manifold, or, like others say that, A-conjugated
# 2. orthogonal gradient in bottom space

# especially for multi-dimensional quadratic function, f(x) = 1/2*x'Ax - b'x, where A is positive defined 
# symmetric matrix, b is an arbitrary vector and has the same length as x, varaible vector in n-dimensional
# space.
# Therefore A can be decomposed into product of one matrix U and its transposition U':
# A = U'U
# also b can be treated as a vector rotated by inversion of matrix U, and original vector is denoted as b0
# f(x) = 0.5*x'Ax - b'x
#      = 0.5*(x'U'Ux) - b'U^{-1}(Ux)
# so we have a new variable Ux now. It is, vector on A mainfold. Denote Ux as X, f(x) = F(X)
# F(x) = 0.5*X'X - b0'X. we can even absorb 0.5 into A but it is trivial. So up to now we obtain a function
# , whose contour is on manifold f(x) and if it is projected directly on xy-plane, that contour will be that
# of f(x). In short: we draw a contour, which is ellipse on xy-plane but perfect circle on manifold. Transition
# matrix is just U.

# So CG keeps orthogonal proceeds on manifold, but keeps gradients orthogonal in real space.
# If all in the same space, that will be a revision of Steep Descent, I think.

# core of iteration of CG:
# position update: 
# x_{k+1} = x_{k} + alpha_{k}*d_{k}             ...(1)
# stepsize update: 
# d_{k+1} = -g_{k} + beta_{k}*d_{k}             ...(2)
# orthogonal relation:
# d_{i}'Ad_{j} = Kronecker_{ij}                 ...(3)
# g_{i}'g_{j} = Kronecker_{ij}                  ...(4)
# precise linear searching:
# alpha_{k} = -(g_{k}'d_{k})/(d_{k}'Ad_{k})     ...(5)
# derivated orthodox:
# g_{i}'d_{j} = Kronecker_{ij}                  ...(6)

    '''
    # conjugate gradient algorithm\n
    1. input requirement:\n
    xfuncy: external function that has form like f(x), this code is based on formulation of quadratic type function\n
    xgrady: gradient of external function, also should be input as an object\n
    x0: initial guess of variable that to be optimized\n
    epsilon: convergence threshold measuring gradient\n
    betamode: for nonlinear optimization task, available options: FR, Dixon, DY, CW, HS, PRP. FR is default for linear
    optimization, others are for nonlinear optimization\n
    conservatism: if set to false, program will not stop if no convergence is reached, use with care\n
    maxloop: if conservatism = True, for optimization step exceeds maxloop, optimization will stop and output 
    variable immediately\n
    verbosity: information level that program prints on screen, available options: 'low', 'medium', 'high', 'debug'\n
    memoryout: during optimization, all parameters will be stored temporarily in dictionary, if set to true, output
    will contain such information instead of outputing optimized variable solely\n
    2. output description\n
    see memoryout parameter in input requirement section
    '''

    epsilon = 10**(-acc)
    alpha_mem = []
    x_mem = []
    grad_mem = []
    beta_mem = []
    disp_mem = []
    zero_x = np.zeros_like(x0)
    b = np.dot(-1, xgrady(zero_x))

    def cg_alpha():
        return np.dot(-1, np.dot(g, d)/(xfuncy(d)+np.dot(b, d)))
    def cg_beta(mode):

        if mode == 'FR':
            return np.dot(g, g)/np.dot(grad_mem[-1], grad_mem[-1])
        elif mode == 'Dixon':
            return -np.dot(g, g)/np.dot(d,grad_mem[-1])
        elif mode == 'DY':
            return np.dot(g, g)/np.dot(d, mop.matrix_minus(g, grad_mem[-1]))
        elif mode == 'CW' or mode == 'HS':
            return np.dot(g, mop.matrix_minus(g, grad_mem[-1]))/np.dot(d, mop.matrix_minus(g, grad_mem[-1]))
        elif mode == 'PRP':
            return np.dot(g, mop.matrix_minus(g, grad_mem[-1]))/np.dot(grad_mem[-1], grad_mem[-1])

    x = x0
    g = xgrady(x)
    d = np.dot(-1, g)
    alpha = 0
    beta = 0

    iloop = 0
    while np.linalg.norm(d) > epsilon and iloop < maxloop:

        alpha_mem.append(alpha)
        alpha = cg_alpha()
        x_mem.append(x)
        if verbosity == 'high' or verbosity == 'debug':
            print('\nCG| optimization step '+str(iloop)+'\n'
                 +'    variable = '+str(x)+'\n'
                 +'    present value = '+str(xfuncy(x)))
        x = mop.matrix_plus(x, np.dot(alpha, d))
        grad_mem.append(g)
        g = xgrady(x)
        beta_mem.append(beta)
        beta = cg_beta(mode = betamode)
        disp_mem.append(d)
        d = mop.matrix_minus(np.dot(beta, d), g)
        iloop += 1
        if conservatism == False:
            iloop = 1

    if memoryout:
        outdict = {}
        outdict['result'] = x
        outdict['cache_alpha'] = alpha_mem
        outdict['cache_beta'] = beta_mem
        outdict['cache_x'] = x_mem
        outdict['cache_grad'] = grad_mem
        outdict['cache_disp'] = disp_mem
        return outdict
    else:
        return x

def newton(
    x0,
    xfuncy,
    xgrady,
    xhessy = 0,
    quasi = False,
    rank = 1,
    wbroyden = 1,
    linearSearch = False,
    armijo_beta = 0.5, 
    armijo_sigma = 0.2, 
    acc = 3,
    conservatism = True,
    maxloop = 50,
    alpha_sd = 1, 
    verbosity = 'low'):

# Developer note:
# --------------
# the other main optimization algorithm style is Newton, which will estimate Hessian matrix directly or
# indirectly. For classical Newton method, Hessian is precisely calculated. For quasi-Newton, Hessian
# will never be calculated but use error vector to correct quasi-Hessian at very step. There are also
# different method to generate quasi-Hessian correction matrix. One names rank1 method that correction
# matrix is generated from only on error vector: alpha*|u_{k}><u_{k}|, which is just rank one because in
# space described by this correction matrix, there is only one vector! so it calls rank1.
# also rank2 will use 2 error vector.

# both classical and quasi version Newton methods are based on second order Talor expansion, of multi-d
# function f(x) at x_{k}:
# f(x) ~= f(x_{k}) + g_k'(x-x_{k}) + 0.5*(x-x_{k})'G_k(x-x_{k})
# gradient of approximated f(x):
# g(x) = g_k + G_k(x-x_{k})

# 1. classical Newton:
# there must be one point nearby where g(x) = 0 and will be reached at the next step:
# -g_{k} = G_{k}(x_{k+1}-x{k})
# denote displacement of x_{k} as d_{k}, and
# -g_{k} = G_{k}d_{k}, solve it and find d_{k} then x_{k+1} = x_{k} + alpha_{k}*d_{k}
# alpha_{k} can either be estimated by precise linear searching or, fuzzy Armijo

# 2. quasi-Newton (rank1):
# g(x) = g_k + G_k(x-x_{k}), use expansion where x_{k+1}, x:= x_{k}, then:
# g_{k} - g_{k+1} = -y_{k} = -G_{k+1}s_{k}

# for G_{k+1}, G_{k+1} ~= B_{k} + E_{k}
# where E_{k} = alpha*|u_{k}><u_{k}|
# rewrite: (B_{k} + E_{k})|s_{k}> = |y_{k}>
# expand: B_{k}|s_{k}> + alpha|u_{k}><u_{k}|s_{k}> = |y_{k}>
#         alpha<u_{k}|s_{k}>|u_{k}> = (|y_{k}> - B_{k}|s_{k}>)
# therefore, |u_{k}> = beta*(|y_{k}> - B_{k}|s_{k}>)
#            E_{k} = alpha*beta**2*(|y_{k}> - B_{k}|s_{k}>)(|y_{k}> - B_{k}|s_{k}>)'
# back to (B_{k} + E_{k})|s_{k}> = |y_{k}>, <- E_{k}
#         E_{k}|s_{k}> = |y_{k}> - B_{k}|s_{k}>
#         alpha*beta**2*(|y_{k}> - B_{k}|s_{k}>)[(|y_{k}> - B_{k}|s_{k}>)'|s_{k}>] = |y_{k}> - B_{k}|s_{k}>
# therefore, alpha*beta**2 = 1/[(|y_{k}> - B_{k}|s_{k}>)'|s_{k}>]
#            E_{k} = (|y_{k}> - B_{k}|s_{k}>)(|y_{k}> - B_{k}|s_{k}>)'/(|y_{k}> - B_{k}|s_{k}>)'|s_{k}>
# use I (unit matrix) as initial guess of Hessian

# 3. Broyden-Fletcher-Goldfarb-Shanno (BFGS)
# they use E_{k} = alpha*|u_{k}><u_{k}| + beta*|u_{k}><u_{k}|
# all others are the same, iterative expression is directly given as below:
#

# 4. Davidson-Fletcher-Powell (DFP)
# define inversion of B_{k} as H_{k}, iterative expression is directly given as below:
# 

# 5. Broyden
# linear combination of BFGS and DFP


    '''
    # newton type algorithm\n
    1. input requirement\n
    x0: initial guess of variable that to be optimized\n
    xfuncy: external function that has form like f(x), this code is based on formulation of quadratic type function\n
    xgrady: gradient function of xfuncy, also should be input as an object\n
    xhessy: hessian matrix function of xfuncy, also should be input as an object\n
    [quasi, rank, wbroyden]: different Newton type optimization are accessible by setting these three parameters\n
    classical Newton: quasi = False\n
    quasi Newton: quasi = True, rank = 1\n
    BFGS (Broyden-Fletcher-Goldfarb-Shanno): quasi = True, rank = 2, wbroyden = 1\n
    DFP (Davidson-Fletcher-Powell): quasi = True, rank = 2, wbroyden = 0\n
    broyden: quasi = True, rank = 2, 0 < wbroyden < 1\n
    linearSearch: linear search method for determining scaling factor of step, activate if set to exact name of method,
    linear search method supported: 'simple', 'golden', 'quadratic', 'armijo', 'accurate' or fractional number directly\n
    armijo_beta: special for fuzzy linear search, Armijo method\n
    armijo_sigma: special for fuzzy linear search, Armijo method\n
    acc: epsilon = 10**(-acc)\n
    conservatism: if set to false, program will not stop if no convergence is reached, use with care\n
    maxloop: if conservatism = True, for optimization step exceeds maxloop, optimization will stop and output 
    variable immediately\n
    verbosity: information level that program prints on screen, available options: 'low', 'medium', 'high', 'debug'\n
    '''
    epsilon = 10**(-acc)

    if quasi == False:

        if verbosity == 'high' or verbosity == 'debug':
            print('NEWTON| classical Newton method is activated...')
        if xhessy == 0:
            if verbosity != 'low':
                print('NEWTON| ***error*** for classical Newton method, analytical Hessian matrix must be given! quit.')
            exit()
        
        iloop = 0
        x = x0
        while np.linalg.norm(xgrady(x)) > epsilon and iloop < maxloop:

            g = xgrady(x)
            hessian = xhessy(x)
            d = np.linalg.solve(hessian, np.dot(-1, g))

            if linearSearch == False:
                alpha = 1
            elif linearSearch == 'simple':
                pass
            elif linearSearch == 'golden':
                pass
            elif linearSearch == 'quadratic':
                pass
            elif linearSearch == 'armijo':
                s = armijo_step_revision(
                        xfuncy = xfuncy,
                        xgrady = xgrady,
                        x = x,
                        dx = d,
                        sigma = armijo_sigma,
                        beta = armijo_beta,
                        verbosity = verbosity
                    )
            elif linearSearch == 'accurate':
                alpha = -np.dot(g, d)/np.dot(np.matmul(hessian, d), d)
                s = np.dot(alpha, d)
            elif type(linearSearch) == type(0.1):
                alpha = linearSearch
                s = np.dot(alpha, d)
            else:
                print('NEWTON| invalid linear search method input, quit.')
                exit()
            
            x = mop.matrix_plus(x, s)
            iloop += 1
            if conservatism == False:
                iloop = 1

    else:
        # quasi-Newton method...
        if rank == 1:
            # classical quasi-newton
            x = x0
            hessian_inv_apx = np.eye(len(x))
            hessian_inv_apx_0 = hessian_inv_apx
            iloop = 0
            g = xgrady(x)

            while np.linalg.norm(g) > epsilon and iloop < maxloop:
                
                if verbosity != 'low':
                    print('\nNEWTON| optimization step '+str(iloop+1))
                d = np.dot(-1, np.matmul(hessian_inv_apx, g))
                if verbosity == 'debug':
                    print('NEWTON| unscaled displacement updated: '+str(d))
                
                if linearSearch == False:
                    alpha = 1
                elif linearSearch == 'simple':
                    pass
                elif linearSearch == 'golden':
                    pass
                elif linearSearch == 'quadratic':
                    pass
                elif linearSearch == 'armijo':
                    pass
                elif linearSearch == 'accurate':
                    hessian = np.linalg.inv(hessian_inv_apx)
                    alpha = -np.dot(g, d)/np.dot(np.matmul(hessian, d), d)
                elif type(linearSearch) == type(0.1):
                    alpha = linearSearch
                else:
                    print('NEWTON| invalid linear search method input, quit.')
                    exit()
                s = np.dot(alpha, d)
                if verbosity == 'debug':
                    print('NEWTON| (linear search based) scaled displacement updated: '+str(s))
                x = mop.matrix_plus(x, s)
                if verbosity == 'debug':
                    print('NEWTON| variable updated: '+str(x))
                _g = g
                g = xgrady(x)
                if verbosity == 'debug':
                    print('NEWTON| gradient updated: '+str(g))
                y = mop.matrix_minus(g, _g)
                if verbosity == 'debug':
                    print('NEWTON| gradient displacement updated: '+str(y))

                u = mop.matrix_minus(s, np.matmul(hessian_inv_apx, y))
                if np.linalg.norm(u) == 0:
                    if verbosity == 'high' or verbosity == 'debug':
                        print('NEWTON| Hessian matrix is accurately obtained')
                    break
                if verbosity == 'debug':
                    print('NEWTON| error basis vector updated: '+str(u))
                hessian_inv_update = mop.ketbra(u, u)/np.dot(u, y)
                hessian_inv_apx_0 = hessian_inv_apx

                if np.dot(y, y) <= 0:
                    hessian_inv_apx = hessian_inv_apx_0
                else:
                    hessian_inv_apx = mop.matrix_plus(hessian_inv_apx_0, hessian_inv_update)
                if verbosity == 'debug':
                    print('NEWTON| approximated inversion of Hessian: '+str(hessian_inv_apx))
                if verbosity != 'low':
                    print('NEWTON| optimization value: '+str(xfuncy(x)))
                
                iloop += 1

                if conservatism == False:
                    iloop = 1

        elif rank == 2:
            if wbroyden == 0:
                algostr = 'Broyden-Fletcher-Goldfarb-Shanno (BFGS)'
            elif wbroyden == 1:
                algostr = 'Davidson-Fletcher-Powell (DFP)'
            else:
                algostr = 'Broyden class'
            x = x0
            hessian_inv_apx = np.eye(len(x0))
            hessian_inv_apx_0 = hessian_inv_apx
            g = xgrady(x)

            iloop = 0
            while np.linalg.norm(g) > epsilon and iloop < maxloop:

                if conservatism == False:
                    iloop = 0
                if verbosity != 'low':
                    print('\nNEWTON| '+algostr+' step '+str(iloop + 1)+' value = '+str(xfuncy(x)))

                d = np.dot(-1, np.matmul(hessian_inv_apx, g))
                if verbosity == 'high' or verbosity == 'debug':   
                    print('        unscaled displacement d = '+str(d))             
                if linearSearch == False:
                    alpha = 1
                elif linearSearch == 'simple':
                    pass
                elif linearSearch == 'golden':
                    pass
                elif linearSearch == 'quadratic':
                    pass
                elif linearSearch == 'armijo':
                    pass
                elif linearSearch == 'accurate':
                    hessian = np.linalg.inv(hessian_inv_apx)
                    alpha = -np.dot(g, d)/np.dot(np.matmul(hessian, d), d)
                elif type(linearSearch) == type(0.1):
                    alpha = linearSearch
                else:
                    print('NEWTON| invalid linear search method input, quit.')
                    exit()

                s = np.dot(alpha, d)
                if verbosity == 'high' or verbosity == 'debug':   
                    print('        linear search scaled displacement s = '+str(s))
                x = mop.matrix_plus(x, s)
                if verbosity == 'high' or verbosity == 'debug':   
                    print('        variable updated x = '+str(x))
                _g = g
                g = xgrady(x)
                if verbosity == 'high' or verbosity == 'debug':   
                    print('        gradient updated g = '+str(g))
                y = mop.matrix_minus(g, _g)

                hessian_inv_apx_0 = hessian_inv_apx

                hessian_inv_update_2_up = mop.ketbra(s, s)
                hessian_inv_update_2_down = np.dot(s, y)
                if hessian_inv_update_2_down <= 0:
                    hessian_inv_apx = hessian_inv_apx_0
                    iloop += 1
                    continue
                else:
                    hessian_inv_update_1_up = mop.ketbra(np.matmul(hessian_inv_apx_0, y), np.matmul(hessian_inv_apx_0, y))
                    hessian_inv_update_1_down = -np.dot(y, np.matmul(hessian_inv_apx_0, y))
                    hessian_inv_update_1 = np.dot(1/hessian_inv_update_1_down, hessian_inv_update_1_up)
                    
                    hessian_inv_apx = mop.matrix_plus(hessian_inv_apx, hessian_inv_update_1)

                    hessian_inv_update_2 = np.dot(1/hessian_inv_update_2_down, hessian_inv_update_2_up)

                    hessian_inv_apx = mop.matrix_plus(hessian_inv_apx, hessian_inv_update_2)
                    if wbroyden > 0:
                        # DFP and Broyden
                        v_factor = np.sqrt(np.dot(y, np.matmul(hessian_inv_apx_0, y)))
                        v_vector_1 = np.dot(s, 1/np.dot(y, s))
                        v_vector_2 = np.dot(np.matmul(hessian_inv_apx_0, y), hessian_inv_update_1_down)
                        v = np.dot(v_factor, mop.matrix_plus(v_vector_1, v_vector_2))

                        hessian_inv_update_3 = np.dot(wbroyden, mop.ketbra(v, v))

                        hessian_inv_apx = mop.matrix_plus(hessian_inv_apx, hessian_inv_update_3)
                    else:
                        # BFGS
                        pass
                    iloop += 1

    return x
