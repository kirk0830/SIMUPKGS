#import basisReadIn

#print(basisReadIn.getBasis('STO-3G', verbosity= 'debug'))
import numpy as np
A = [[1, 0.45077041164778764], [0.45077041164778764, 1]]
print(type(A))

'''
A = np.zeros(shape = (2, 2))
B = np.ones_like(A)
C = [A, B]
print(C[1])
print(C[1][:][:])

'''
from A_Szabo_RHF_HeH import scf as scf

iop = 2
N = 3
R = 1.4632
zeta1 = 2.0925
zeta2 = 1.24
ZA = 2.0
ZB = 1.0

result = scf(
    N_GauFun = N,
    R_bond = R,
    zeta_A = zeta1,
    zeta_B = zeta2,
    charge_A = ZA,
    charge_B = ZB,
    conv_thre = 1E-6,
    maxscf = 50,
    verbosity = 'debug'
)
print(result)
