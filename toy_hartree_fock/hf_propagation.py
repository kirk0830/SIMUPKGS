#from diag_numpy import npDiag as diag
from diag_so2 import diag as diag
import numpy as np
import _in_matrix_op as mop
def Gij(denMat, Vee, mode = 'restricted'):

    """
    electron-electron interaction part of Fock operator\n
    Multiple summation comes from expansion of wavefunction-dependent operator\n
    g_ij := <i|1/r|j>, where i and j are needed to be expanded with atomic wavefunction basis.\n
    During which, density matrix is induced (from mapping operation from canonical HF equation \n
    solution spanned space to atomic wavefunction spanned space).\n
    Therefore for every subsript pair (i, j) of G, there are N*N (k, l)-pairs summation needed\n
    to be carried out.\n
    Density matrix element denoted as P_kl,\n
    multiplies <k| and |l>, and then <i| and |j>\n
    denMat: density matrix\n
    Vee: electron-electron four-subscript interaction matrix\n
    mode: (not fully implemented yet) restricted or unrestricted Hartree-Fock method. \n
    For unrestricted HF method, two density matrices and Vee should be provided.
    """
    if mode == 'restricted':
        NBasis = np.shape(denMat)[0]
        # NBasis is number of atomic basis wavefunction, not that of Gaussian functions!

        G = np.zeros(shape = (NBasis, NBasis))
        for i in range(NBasis):
            for j in range(NBasis):
                # definition of every element in G matrix,
                # needs summation over all (k, l)-pairs
                for k in range(NBasis):
                    for l in range(NBasis):
                        # Vee: exchange j <-> l to convert Columb to Exchange
                        G[i][j] += denMat[k][l]*(Vee[i][j][k][l]-0.5*Vee[i][l][k][j])
                    
    elif mode == 'unrestricted':
        pass
    else:
        print('ERROR | Mode selected not implemented yet.')

    return G

def hf_fp(Tij, natom, VneList, Gij, printH = True):

    """
    Hartree-Fock forward propagation function\n
    output description:\n
    printH == True: [Hcore, F]\n
    printH == False: F\n
    input requirement:\n
    Tij: kinetic integral matrix\n
    natom: number of atoms, also number of Vne integral matrix\n
    VneList: list of Vne, nuclear attraction integral matrix\n
    Gij: electron-electron interaction integral matrix
    """
    Hij = np.zeros_like(Tij)
    Fij = np.zeros_like(Tij)
    Hij = mop.matrix_plus(Hij, Tij)
    
    if natom > 1:
        for i in range(natom):
            Hij = mop.matrix_plus(Hij, VneList[i])
    elif natom == 1:
        Hij = mop.matrix_plus(Hij, VneList)
    
    Fij = mop.matrix_plus(Hij, Gij)

    if printH:
        return [Hij, Fij]
    else:
        return Fij

# Given FC = eSC, to find C
# S always has non-zero off-diagonal terms, cuz atomic wavefunction is uncessarily orthonormal
# e is diagonal matrix cuz it is canonical HF(R) equation
# find U'SU = s (diagonalization)
# find X^{-1}SX = I (diagonalization + normalization)
# thus we use canonical orthogonalization, form X = U/sqrt(s)
# i.e. s^{-1/2} U' S U s^{-1/2} = I
# U <- classical diagonalization of S
# s^{1/2}

# X^{-1}F        C = X^{-1}eS        C
# X^{-1}F        C = eX^{-1}S        C
# X^{-1}FX X^{-1}C = eX^{-1}SX X^{-1}C
#       F'       C'= e               C' <- classical eigenvalue problem, diagonalize F' to find C'
# C' = X^{-1}C -> C = XC'


# this function gives new coefficient matrix and density matrix

def hf_bp(
    Fij, 
    Sij, 
    printP = True, 
    restricted = True, 
    singular_exclu = False, 
    eigValThre = 1E-5, 
    verbosity = 'debug'
    ):

    """
    Hartree-Fock backward propagation function\n
    output description:\n
    printP == True: [C, P]\n
    printP == False: C\n
    input requirement:\n
    Fij: Fock matrix\n
    Sij: overlap matrix\n
    printP: if calculate density matrix\n
    restricted: if present calculation uses restricted scheme\n
    singluar_exclu: excludes extremely small eigenvalues of Sij by truncating, not implemented yet\n
    eigValThre: threshold for truncating diagonalized Sij, not implemented yet\n
    """
    [s, U] = diag(mat = Sij)
    sqrtInvDiagS = np.zeros_like(s)
    for i in range(np.shape(s)[0]):

        sqrtInvDiagS[i][i] = 1/np.sqrt(s[i][i])

    X = np.matmul(U, sqrtInvDiagS)

    Fij_prime = np.matmul(np.transpose(X), Fij)
    Fij_prime = np.matmul(Fij_prime, X)

    [e_prime, C_prime] = diag(mat = Fij_prime)
    C = np.matmul(X, C_prime)
    if verbosity == 'debug':
        print('SIMUPKGS_DEBUG_MODE | Uij = '+str(U))
        print('SIMUPKGS_DEBUG_MODE | s^{-1/2}ij = '+str(sqrtInvDiagS))
        print('SIMUPKGS_DEBUG_MODE | Xij = '+str(X))
        print('SIMUPKGS_DEBUG_MODE | Fij\' = '+str(Fij_prime))
        print('SIMUPKGS_DEBUG_MODE | Cij\' = '+str(C_prime))
        print('SIMUPKGS_DEBUG_MODE | Cij = '+str(C))

    if printP:
        if restricted:
            # only half of C will be used to calculated P
            transC = np.transpose(C)
            restrictedC = np.zeros_like(C)
            for i in range(0, np.shape(C)[0], 2):

                restrictedC[i][:] = transC[i][:]
                restrictedC[i+1][:] = transC[i][:]

            restrictedC = np.transpose(restrictedC)
            P = np.matmul(restrictedC, np.transpose(restrictedC))
        else:
            P = np.matmul(C, np.transpose(C))
        return [C, P]

    else:

        return C

def hf_elecEner(Fij, Hij, Pij):

    return 0.5*np.sum(mop.matrix_dot_multiply(Pij, mop.matrix_plus(Fij, Hij)))
