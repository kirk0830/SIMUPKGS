# integrals in this file are directly referred from textbook
# Quantum Chemistry -- basic concept and ab initio method (2nd edn), Xu Guangxian, et al.,
# 2009, Beijing: Science Press.

# Chapter 10: integrals in quantum chemistry: Gaussian functions

import numpy as np

def df(n):
    '''
    full name: double factorial function: n!!\n
    # input requirement\n
    n: integar\n
    # output description\n
    integar
    '''
    if n <= 0:
        return 1
    else:
        return n * df(n-2)

def C(N, n):

    return np.math.factorial(N)/np.math.factorial(n)/np.math.factorial(N-n)

def ic_gamma_pade(m, w):

    """
    incompleted gamma integral, pade approximated version, polynomial\n
    # input requirement\n
    m: order\n
    w: variable\n
    # output description\n
    float number
    """
    a0 = 1/(2*m + 1)**(2/(2*m + 1))

    if m == 0:
        a = [
            0.213271302431420E0,
            0.629344460255614E-1,
            0.769838037756759E-2,
            0.758433197127160E-3,
            0.564691197633667E-4
        ]
        b = [
            0.879937801660182E0,
            0.338450368470103E0,
            0.738522953299624E-1,
            0.101431553402629E-1,
            0.955528842975585E-3,
            0.720266520392572E-4
        ]
        wmax = 16.3578
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 1:
        a = [
            0.295195994716045E-1,
            0.128790985465415E-1,
            0.998165499553218E-3,
            0.970927983276419E-4,
            0.493839847029699E-5
        ]
        b = [
            0.461403194679124E0,
            0.108494164372449E0,
            0.171462934845042E-1,
            0.196918657845508E-2,
            0.160138863265254E-3,
            0.857708713007233E-5
        ]
        wmax = 17.4646
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 2:
        a = [
            -0.575763488635418E-2,
            0.731474973333076E-2,
            0.251276149443393E-3,
            0.264336244559094E-4
        ]
        b = [
            0.274754154712841E0,
            0.425364830353043E-1,
            0.493902790955943E-2,
            0.437251500927601E-3,
            0.288914662393981E-4
        ]
        wmax = 15.2368
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 3:
        a = [
            -0.290110430424666E0,
            0.561884370781462E-2,
            0.301628267382713E-3,
            0.110671035361856E-4
        ]
        b = [
            0.171637608242892E0,
            0.187571417256877E-1,
            0.178536829675118E-2,
            0.137360778130936E-3,
            0.791915206883054E-5
        ]
        wmax = 16.0419
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 4:
        a = [
            -0.452693111179624E-1,
            0.490070062899003E-2,
            -0.561879719979307E-4,
            0.550814626951998E-5,
        ]
        b = [
            0.108051989937231E0,
            0.855924943430755E-2,
            0.724968571389473E-3,
            0.502338223156067E-4,
            0.249107837399141E-5
        ]
        wmax = 16.8955
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 5:
        a = [
            -0.566143259316101E-1,
            0.455916894577203E-2,
            -0.894152721395639E-4,
            0.328096732308082E-5
        ]
        b = [
            0.662932958471386E-1,
            0.383724443872493E-2,
            0.327167659811839E-3,
            0.210490437682548E-4,
            0.883562935089333E-6
        ]
        wmax = 17.7822
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 6:
        a = [
            -0.503259167534352E-1,
            0.273135625430953E-2,
            -0.310733624819100E-4
        ]
        b = [
            0.586609328033371E-1,
            0.194044691497128E-2,
            0.109442742602192E-3,
            0.613406236401726E-5
        ]
        wmax = 15.8077
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 7:
        a = [
            -0.548201062615785E-1,
            0.253099908233175E-2,
            -0.333589469427863E-4
        ]
        b = [
            0.389873128779298E-1,
            0.569890860832729E-3,
            0.422187129333708E-4,
            0.286010059144633E-5
        ]
        wmax = 16.5903
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 8:
        a = [
            -0.581618006078160E-1,
            0.238525529084601E-2,
            -0.329989020317093E-4
        ]
        b = [
            0.240929282666615E-1,
            -0.202677647499956E-3,
            0.119820675974460E-4,
            0.145792086904409E-5
        ]
        wmax = 17.3336
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 9:
        a = [
            -0.334843993901400E-1,
            0.846637494147059E-3,
        ]
        b = [
            0.495875606944471E-1,
            0.946642302340943E-3,
            0.108367772249790E-4
        ]
        wmax = 15.6602
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 10:
        a = [
            -0.335292171805959E-1,
            0.749168957990503E-3
        ]
        b = [
            0.421492932021683E-1,
            0.582840985360327E-3,
            0.237676790577455E-5
        ]
        wmax = 16.5258
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 11:
        a = [
            -0.332669773770348E-1,
            0.668720489602687E-3
        ]
        b = [
            0.363057685289467E-1,
            0.345646100984643E-3,
            -0.19087233373450E-5
        ]
        wmax = 17.5395
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 12:
        a = [
            -0.326241966410798E-1,
            0.598705175467956E-3
        ]
        b = [
            0.318680048277695E-1,
            0.202419662347765E-3,
            -0.362095173837973E-5
        ]
        wmax = 18.5783
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 13:
        a = [
            -0.317754368014894E-1,
            0.537678595933584E-3
        ]
        b = [
            0.284036027081815E-1,
            0.113673420662576E-3,
            -0.416076810552774E-5
        ]
        wmax = 19.6511
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 14:
        a = [
            -0.308755854748829E-1,
            0.485046451960769E-3
        ]
        b = [
            0.255694625434059E-1,
            0.542010192055080E-4,
            -0.424759498527876E-5
        ]
        wmax = 20.7839
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 15:
        a = [
            -0.300143806671997E-1,
            0.439983032427912E-3
        ]
        b = [
            0.231478878674366E-1,
            0.105546581596674E-4,
            -0.418932957034726E-5
        ]
        wmax = 21.9998
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    elif m == 16:
        a = [
            -0.288346417991609E-1,
            0.397161796318408E-3
        ]
        b = [
            0.215021933120724E-1,
            -0.128592457453950E-5,
            -0.362120651688135E-5
        ]
        wmax = 20.9225
        if w > wmax:
            return df(2*m-1)/(2*w)**(m+0.5)*(np.pi)**0.5
        else:
            term1 = a0
            term2 = 1
            for ia in range(len(a)):
                term1 += w**(ia+1) * a[ia]
            for ib in range(len(b)):
                term2 += w**(ib+1) * a[ib]
            return (term1/term2)**(m+0.5)
    else:
        print('SIMUPKGS| ***error*** m parameter of incompleted gamma function Fm(w) is not reachable now.')
        exit()

def Sij_Gau(a, b, atom_A, atom_B, aglr_a = [0, 0, 0], aglr_b = [0, 0, 0]):
    
    """
    generalized overlap integral\n
    # input requirement\n
    a: alpha coefficient of Gaussian function a\n
    b: alpha coefficient of Gaussian function b\n
    atom_A: coordinate xyz of atom A\n
    atom_B: coordinate xyz of atom B\n
    aglr_a: angular momentum list along three direction, x, y and z of atom A\n
    aglr_b: same as aglr_a, that of atom B\n
    # output description\n
    overlap integral, element data
    """

    # geometry...
    rABx = abs(atom_A[-3] - atom_B[-3])
    rABy = abs(atom_A[-2] - atom_B[-2])
    rABz = abs(atom_A[-1] - atom_B[-1])

    rAB2 = rABx**2 + rABy**2 + rABz**2
    rAB = np.sqrt(rAB2)

    rPA = (b/(a+b))*rAB
    rPAx = (b/(a+b))*rABx
    rPAy = (b/(a+b))*rABy
    rPAz = (b/(a+b))*rABz

    rPB = (a/(a+b))*rAB
    rPBx = (a/(a+b))*rABx
    rPBy = (a/(a+b))*rABy
    rPBz = (a/(a+b))*rABz

    [la, ma, na] = aglr_a
    [lb, mb, nb] = aglr_b

    Ix = 0
    Iy = 0
    Iz = 0

    sum_lam_mu_x = la + lb
    sum_lam_mu_y = ma + mb
    sum_lam_mu_z = na + nb

    K = np.exp(-a*b/(a+b)*rAB2)

    for i in range(1+np.math.floor(sum_lam_mu_x/2)):
        f2i = 0
        for lam in range(1+la):
            mu = 2*i - lam
            if mu <= lb and mu >= 0:
                f2i += C(la, lam)*C(lb, mu)*rPAx**(la-lam)*rPBx**(lb-mu)
        Ix += f2i*np.math.gamma(i+0.5)/(a+b)**(i+0.5)

    for i in range(1+np.math.floor(sum_lam_mu_y/2)):
        f2i = 0
        for lam in range(1+ma):
            mu = 2*i - lam
            if mu <= mb and mu >= 0:
                f2i += C(ma, lam)*C(mb, mu)*rPAy**(ma-lam)*rPBy**(mb-mu)
        Iy += f2i*np.math.gamma(i+0.5)/(a+b)**(i+0.5)

    for i in range(1+np.math.floor(sum_lam_mu_z/2)):
        f2i = 0
        for lam in range(1+na):
            mu = 2*i - lam
            if mu <= nb and mu >= 0:
                f2i += C(na, lam)*C(nb, mu)*rPAz**(na-lam)*rPBz**(nb-mu)
        Iz += f2i*np.math.gamma(i+0.5)/(a+b)**(i+0.5)

    return K*Ix*Iy*Iz

def Tij_Gau(a, b, atom_A, atom_B, aglr_a = [0, 0, 0], aglr_b = [0, 0, 0]):

    """
    generalized kinetic integral\n
    # input requirement\n
    a: alpha coefficient of Gaussian function a\n
    b: alpha coefficient of Gaussian function b\n
    atom_A: coordinate xyz of atom A\n
    atom_B: coordinate xyz of atom B\n
    aglr_a: angular momentum list along three direction, x, y and z of atom A\n
    aglr_b: same as aglr_a, that of atom B\n
    # output description\n
    overlap integral, element data
    """
    [la, ma, na] = aglr_a
    [lb, mb, nb] = aglr_b

    if lb > 1:
        Ix_term1 = -0.5*lb*(lb-1)*Sij_Gau(a = a,
                                          b = b,
                                          atom_A = atom_A,
                                          atom_B = atom_B,
                                          aglr_a = [la, ma, na],
                                          aglr_b = [lb - 2, mb, nb]
                                          )
    else:
        Ix_term1 = 0
    Ix_term2 = b*(2*lb+1)*Sij_Gau(a = a,
                                  b = b,
                                  atom_A = atom_A,
                                  atom_B = atom_B,
                                  aglr_a = [la, ma, na],
                                  aglr_b = [lb, mb, nb]
                                  )
    Ix_term3 = -2*b**2*Sij_Gau(a = a,
                               b = b,
                               atom_A = atom_A,
                               atom_B = atom_B,
                               aglr_a = [la, ma, na],
                               aglr_b = [lb + 2, mb, nb]
                               )
    Ix = Ix_term1 + Ix_term2 + Ix_term3

    if mb > 1:
        Iy_term1 = -0.5*mb*(mb-1)*Sij_Gau(a = a,
                                          b = b,
                                          atom_A = atom_A,
                                          atom_B = atom_B,
                                          aglr_a = [la, ma, na],
                                          aglr_b = [lb, mb - 2, nb]
                                          )
    else:
        Iy_term1 = 0
    Iy_term2 = b*(2*mb+1)*Sij_Gau(a = a,
                                  b = b,
                                  atom_A = atom_A,
                                  atom_B = atom_B,
                                  aglr_a = [la, ma, na],
                                  aglr_b = [lb, mb, nb]
                                  )
    Iy_term3 = -2*b**2*Sij_Gau(a = a,
                               b = b,
                               atom_A = atom_A,
                               atom_B = atom_B,
                               aglr_a = [la, ma, na],
                               aglr_b = [lb, mb + 2, nb]
                               )
    Iy = Iy_term1 + Iy_term2 + Iy_term3

    if nb > 1:
        Iz_term1 = -0.5*nb*(nb-1)*Sij_Gau(a = a,
                                          b = b,
                                          atom_A = atom_A,
                                          atom_B = atom_B,
                                          aglr_a = [la, ma, na],
                                          aglr_b = [lb, mb, nb - 2]
                                          )
    else:
        Iz_term1 = 0
    Iz_term2 = b*(2*nb+1)*Sij_Gau(a = a,
                                  b = b,
                                  atom_A = atom_A,
                                  atom_B = atom_B,
                                  aglr_a = [la, ma, na],
                                  aglr_b = [lb, mb, nb]
                                  )
    Iz_term3 = -2*b**2*Sij_Gau(a = a,
                               b = b,
                               atom_A = atom_A,
                               atom_B = atom_B,
                               aglr_a = [la, ma, na],
                               aglr_b = [lb, mb, nb + 2]
                               )
    Iz = Iz_term1 + Iz_term2 + Iz_term3

    return Ix + Iy + Iz

# denote center (atom) coordinate as Ax, Ay and Az
# Gaussian function as G(a, A, l, m, n)
# denote d/dAx as Dx, similarily, Dy and Dz
# DxG(l, m, n) = -lG(l-1, m, n) + 2aG(l+1, m, n)
# DyG(l, m, n) = -mG(l, m-1, n) + 2aG(l, m+1, n)
# DzG(l, m, n) = -nG(l, m, n-1) + 2aG(l, m, n+1)

# DxG(000) = 2aG(100) >>> G2px
# DyG(000) = 2aG(010) >>> G2py
# DzG(000) = 2aG(001) >>> G2pz

# DyG(100) = 2aG(110) >>> G3dxy
# DzG(100) = 2aG(101) >>> G3dxz

# DyG(001) = 2aG(011) >>> G3dyz

# DxG(100) = -G(000) + 2aG(200) >>> G3dx2
# DyG(010) = -G(000) + 2aG(020) >>> G3dy2
# DzG(001) = -G(000) + 2aG(002) >>> G3dz2

