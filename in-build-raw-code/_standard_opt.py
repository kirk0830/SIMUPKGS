import numpy as np
import _in_matrix_op as mop

def cg(
    xfuncy, 
    xgrady, 
    x0, 
    epsilon, 
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
                pass
            elif linearSearch == 'accurate':
                alpha = -np.dot(g, d)/np.dot(np.matmul(hessian, d), d)
            elif type(linearSearch) == type(0.1):
                alpha = linearSearch
            else:
                print('NEWTON| invalid linear search method input, quit.')
                exit()
            
            x = mop.matrix_plus(x, np.dot(alpha, d))
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