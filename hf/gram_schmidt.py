from copy import deepcopy
import _mat_lib as mlib

def gs_orth(vec_set, mode = 'column', start = 0, normalize = True, verbosity = 'silent'):

    """
    Gram-Schmidt orthogonalization\n
    WARNING: input matrix (I mean vector set) must be float data type. An integar data type
    will cause unexpected large roundoff error, even totally wrong result will be returned!\n
    # input requirement\n
    vec_set: set of vectors that need to be orthogonalized, 2d array (list), FLOAT\n
    mode: 'column' or 'row', specify how vectors are placed in matrix the first parameter you
    entered\n
    start: which vector is selected as a reference to orthogonalize all other vectors, note 
    that this number corresponds to index of vector, so the first vector is '0' rather than '1
    '\n
    normalize: True or False, choose whether all vectors output are already normalized\n
    # output description\n
    [original vector set, resultant matrix containing orthogonalized vectors]

    """
    # Gram-Schmidt method
    len_vec = len(vec_set)
    nvec = len(vec_set[-1][:])
    if mode == 'row':
        temp = nvec
        nvec = len_vec
        len_vec = temp
    
    norm_orth_set = mlib.zeros(nvec, m = len_vec)
    orth_set = mlib.zeros(nvec, m = len_vec)

    # original vector collection, no matter how they arrange, <| or |>, save them as <|
    vec_colle = []
    vec_set_bak = deepcopy(vec_set)

    if mode == 'column':

        for icol in range(nvec):

            _vec = []
            _vec = [vec_set[i][icol] for i in range(len_vec)]
            vec_colle.append(_vec)
        if verbosity == 'debug':
            print('GRAM-SCHMIDT| |i> => <i|, column vectors have been saved as <| for easy calculation:{}'.format(vec_colle))
    elif mode == 'row':

        for irow in range(nvec):

            _vec = []
            _vec = [vec_set[irow][i] for i in range(len_vec)]
            vec_colle.append(_vec)
        if verbosity == 'debug':
            print('GRAM-SCHMIDT| row vectors have been saved:{}'.format(vec_colle))

    else:
        print('***error*** invalid mode required in reading input vectors set.')
        exit()
    

    # will save mod of every vector
    mod_vec = []
    for ivec in range(nvec):

        mod = mlib.mod_of_vec(vec = vec_colle[ivec][:])
        mod_vec.append(mod)
    
    # select no. start as the first basis and directly normalize it
    orth_set[start][:] = vec_colle[start][:]
    norm_orth_set[start][:] = [vec_colle[start][i]/mod_vec[start] for i in range(len_vec)]
    if verbosity == 'debug':
        print('GRAM-SCHMIDT| The first basis has been fixed as:\n{}\nGRAM-SCHMIDT| Normalized:\n{}'.format(orth_set[start][:], norm_orth_set[start][:]))
    orthlog = [start]

    for i in range(nvec):

        if i == start:
            continue

        vec2orth = vec_colle[i][:]
        cut = mlib.zeros(n = 1, m = len_vec)[0][:]

        for index in orthlog:
            
            innerprod = mlib.braket(vec2orth, norm_orth_set[index][:])
            compo2cut = [innerprod*item for item in norm_orth_set[index][:]]
            cut = mlib.plus(compo2cut, cut)
            if verbosity == 'debug':
                print('GRAM-SCHMIDT| Present vector to perform Gram-Schmidt:\n{}\nGRAM-SCHMIDT| Basis:\n{}'.format(vec2orth, norm_orth_set[index][:]))
                print('GRAM-SCHMIDT| scalar product between present vector and basis: {}'.format(innerprod))
                print('GRAM-SCHMIDT| vector needed to be subtract from original vector is:\n{}'.format(compo2cut))
        orth_set[i][:] = mlib.minus(vec2orth, cut)
        mod_vec[i] = mlib.mod_of_vec(orth_set[i][:]) #refresh mod info
        norm_orth_set[i][:] = [item/mod_vec[i] for item in orth_set[i][:]]
        print('GRAM-SCHMIDT| Yield new unnormalized and normalized basis:\n{}\n{}'.format(orth_set[i][:], norm_orth_set[i][:]))
        orthlog.append(i)

    if mode == 'column':

        norm_orth_set = mlib.transpose(norm_orth_set)
        orth_set = mlib.transpose(orth_set)

    if normalize:

        return [vec_set_bak, norm_orth_set]
    else:

        return [vec_set_bak, orth_set]