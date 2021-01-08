# householder reflection
import _mat_lib as mlib

def householder(vec_in, vec_desti = [], mode = 'reduce', verbosity = 'silent'):

    len_vec = len(vec_in)
    if mode == 'reduce':

        vec_desti = mlib.zeros(n = 1, m = len_vec)[0][:]
        mod_vec_in = mlib.mod_of_vec(vec_in)
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
                +'normalized norm vector = {}\n'.format(v)
                +'identity operator generated:\n{}\n'.format(I)
                +'Householder operator:\n{}'.format(P))
        return P

