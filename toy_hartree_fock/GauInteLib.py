import numpy as np

# only s-orbital is supported now :( --2020/12/26

def Sij_Gau(A, B, rAB2, angularA = [0, 0, 0], angularB = [0, 0, 0]):

    """
    Overlap integration function\n
    A and B are Gaussian function specific parameters\n
    rAB2: squared distance of center of two integrated Gaussian functions\n
    angularA and -B: angular index (powers of x, y and z at head of Gaussian function) list. For s-orbital, use [0, 0, 0]. For p-orbital, use [1, 0, 0] or [0, 1, 0] or [0, 0, 1], as px, py and pz, respectively.\n
    >>--Up to 2020/12/26, integrations for high angular momentum orbitals are not implemented yet--
    """
    return ((np.pi/(A+B))**1.5)*np.exp(-A*B/(A+B)*rAB2)

def Tij_Gau(A, B, rAB2, angularA = [0, 0, 0], angularB = [0, 0, 0]):

    """
    Kinetic integration function\n
    A and B are Gaussian function specific parameters\n
    rAB2: squared distance of center of two integrated Gaussian functions.\n
    angularA and -B: angular index (powers of x, y and z at head of Gaussian function) list. For s-orbital, use [0, 0, 0]. For p-orbital, use [1, 0, 0] or [0, 1, 0] or [0, 0, 1], as px, py and pz, respectively.\n
    >>--Up to 2020/12/26, integrations for high angular momentum orbitals are not implemented yet--
    """
    return (A*B)/(A+B)*(3-2*A*B*rAB2/(A+B))*((np.pi/(A+B))**1.5)*np.exp(-A*B*rAB2/(A+B))

def numeric_erf(x):

    """
    in built approximated erf function, integrated Gaussian function\n
    denote t = 1/(1+px)\n
    erf(x) = 1 - sum{a_{i}*t^i}*exp(-x^2), i goes from 1 to 5.\n
    denote f = sum{a_{i}*t^i}\n
    p and a_{i} values are given directly...
    """
    A = [0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429]
    p = 0.3275911
    t = 1/(1+p*x)
    tN = t
    f = A[0]*tN
    for i in range(1, 5):
        tN *= t
        f += A[i]*tN
    return 1 - f*np.exp(-x*x)

def F0(t):

    """
    Special function F0(t)\n
    prviding solution for integrating 1/r terms\n
    F0(t) := t^{-1/2} * integral{0}{t^{1/2}}{exp{-y^2} dy}
    """
    if t < 1.0E-6:
        return 1-t/3
    else:
        return np.sqrt(np.pi/t)*numeric_erf(np.sqrt(t))/2

def Vneij_Gau(A, B, rAB2, rCP2, ZC, angularA = [0, 0, 0], angularB = [0, 0, 0]):
    """
    Nuclear attraction integration for s-orbital\n
    Parameter list:\n
    A and B: parameters that represent property of Gaussian function,\n
    more explicitly, A and B are oribtal shrinking coefficients.\n
    
    In the following distance input list, A and B are positions of electrons,\n
    C is the center of two electrons, P is position of core.\n
    rAB2: squared distance between two electrons.\n
    rCP2: squared distance between middle point of two electrons and one core.\n
    ZC: nuclear charge, unit in e\n
    angularA and -B: angular index (powers of x, y and z at head of Gaussian function) list. For s-orbital, use [0, 0, 0]. For p-orbital, use [1, 0, 0] or [0, 1, 0] or [0, 0, 1], as px, py and pz, respectively.\n
    >>--Up to 2020/12/26, integrations for high angular momentum orbitals are not implemented yet--
    """
    return -ZC*(2*np.pi/(A+B))*F0((A+B)*rCP2)*np.exp(-A*B*rAB2/(A+B))

def Veeijkl_Gau(A, B, C, D, rAB2, rCD2, rPQ2, angularA = [0, 0, 0], angularB = [0, 0, 0]):

    """
    Electron-electron repulsion integration, for s-orbital\n
    A, B, C, and D are orbitial specfic parameters.\n
    In following explanation, e-e repulsion will be denoted as <12|34>.\n
    rAB2: squared distance between 1 and 2 electron.\n
    rCD2: squared distance between 3 and 4 electron.\n
    rPQ2: squared distance between centers of 1 and 2 and 3 and 4.\n
    angularA and -B: angular index (powers of x, y and z at head of Gaussian function) list. For s-orbital, use [0, 0, 0]. For p-orbital, use [1, 0, 0] or [0, 1, 0] or [0, 0, 1], as px, py and pz, respectively.\n
    >>--Up to 2020/12/26, integrations for high angular momentum orbitals are not implemented yet--
    """
    Vee = 2*np.pi**2.5/((A+B)*(C+D)*np.sqrt(A+B+C+D))
    Vee *= F0((A+B)*(C+D)/(A+B+C+D)*rPQ2)
    Vee *= np.exp(-A*B*rAB2/(A+B))*np.exp(-C*D*rAB2/(C+D))

    return Vee

