# householder reflection
import _mat_lib as mlib

def householder(vec_in, vec_desti = [], mode = 'reduce', verbosity = 'silent'):

    """
    Householder algorithm\n
    Householder build a reflection operator that can transform one vector to one wanted
    direction\n
    # input requirement\n
    vec_in: one vector that want to reflect, SHOULD BE INPUT AS <|, 1d list\n
    vec_desti: one vector whose direction will vec_in be reflected onto, if not given
    explicitly, program will quit unless KEYWORD mode is set to 'reduce'\n
    mode: 'reduce' or anything else. 'reduce' mode means to reflect vec_in to direction
    along x-axis, i.e., [1, 0, 0, ...]\n
    # output description\n
    Householder operator P
    """
    len_vec = len(vec_in)
    mod_vec_in = mlib.mod_of_vec(vec_in)

    if mode == 'reduce':

        vec_desti = mlib.zeros(n = 1, m = len_vec)[0][:]
        vec_desti[0] = mod_vec_in
        
    u = mlib.minus(vec_in, vec_desti)

    mod_u = mlib.mod_of_vec(u)
    v = mlib.zeros(n = 1, m = len_vec)[0][:]
    for i in range(len(u)):

        if (u[i] == 0) and (mod_u == 0):
            v[i] = 1.0 # normalize manually
        else:
            v[i] = u[i]/mod_u
    # v = [iterm/mod_u for iterm in u]

    I = mlib.eye(len_vec)
    vvT = mlib.ketbra(v, v, mode = '2bra', amplify = 2)
    P = mlib.minus(I, vvT)
    if verbosity == 'debug':
        print('HOUSEHOLDER| Comprehensive report\n'
             +'input check:\n'
             +'vector = {}\n'.format(vec_in)
             +'destination = {}\n'.format(vec_desti)
             +'norm vector (original) = {}\n'.format(u)
             +'normalized norm vector = {}'.format(v))
        print('2|v><v| tensor product  =')
        mlib.matrix_print(vvT)
        print('identity operator generated:\n{}\n'.format(I)
             +'Householder operator:\n{}'.format(P))
    return P

def hshldr_op(mat, hhop):

    pass