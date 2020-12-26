import numpy as np
import GauInteLib as gint
from basisReadIn import getBasis as Pople
from diag_numpy import npDiag as diag
import hf_propagation as hfprop

def iniIntg(N_GauFun, R_bond, zeta_A, zeta_B, charge_A, charge_B):

    """
    Geometry-dependent electrons integrals initialization\n
    output information: [T, Vne_A, Vne_B, Vee, S]\n
    >---ALL FOLLOWING OUTPUT MATRICES ARE IN ATOMIC BASIS SPANNED SPACE---\n
    T: kinetic integral, under present HeH circumstance, 2*2 shaped\n
    Vne_A and Vne_B: nuclear attraction integrals, 2*2 shaped, A: He, B: H\n
    Vee: electron-electron interaction integral, 2*2*2*2 shaped\n
    S: overlap matrix, 2*2 shaped\n
    Input requirement:\n
    N_GauFun: number of Gaussian functions used for expanding every atomic\n
    wavefunction, will be integrated into "basisRReadIn" subroutine and read\n
    number of Gaussian functions automatically\n
    R_bond: system-limited, dependent geometry descriptive parameter. In \n
    present system, defines all geometry information\n
    charge_A and charge_B: atomic charges of core A (He) and core B (H)\n
    -\n
    [!!!WARNING!!!]: benefits from Hartree-Fock-Roothaan method, this function is\n
    only needed to be called for once at very beginning of calculation, if\n
    performing SINGLE-POINT calculation. However, will be called for many times\n
    if during geometry changing calculation.
    """
    # use 1 atomic wavefunction per atom to expand solution of canonical HF equation
    # create descriptive parameter lists for those two atomic wavefunctions
    # this can be generailized to polyatomic system
    # -----------------------GEOMETRY RELATED INFO SEC.-----------------------------
    R2 = R_bond * R_bond
    # -------------------------ATOMIC BASIS INFO SEC.-------------------------------
    d1 = np.zeros(shape = (3,1))
    a1 = np.zeros_like(d1)
    d2 = np.zeros_like(d1)
    a2 = np.zeros_like(d1)
    # ----------------------NUMERICAL BASIS INFO SEC.-------------------------------
    [unscaled_d, unscaled_a] = Pople(basisName = 'STO-3G', angular = 0)

    for i in range(N_GauFun):

        a1[i] = unscaled_a[i]*(zeta_A**2)
        d1[i] = unscaled_d[i]*((2*a1[i]/np.pi)**0.75)
        a2[i] = unscaled_a[i]*(zeta_B**2)
        d2[i] = unscaled_d[i]*((2*a2[i]/np.pi)**0.75)

    # ------------------------------------------------------------------------------
    # MATRIX ELEMENT INITIALIZATION, 1 and 2 are atomic basis function to expand
    # single-electron spin-orbital wavefunction of HeH molecule.
    # That function, is solution of canonical Hartree-Fock equation.
    # FC = eSC
    # To wrap lines, use if True, no other meanings. Contents: Matrix initialization
    if True:
        # Normal overlap integration...
        S11 = 1
        S12 = 0
        S22 = 1
        # Kinetic integration...
        T11 = 0
        T12 = 0
        T22 = 0
        # Nuclear attraction integration...
        # Atom A
        Vne11_A = 0
        Vne12_A = 0
        Vne22_A = 0
        # Atom B
        Vne11_B = 0
        Vne12_B = 0
        Vne22_B = 0
        # electron-electron integration...
        # Remember that operator G_ij of F_ij posses to canonical HF solution inside, so the 
        # medium two come from G_ij operator directly
        V1111 = 0
        V2111 = 0
        V2121 = 0
        V2211 = 0
        V2221 = 0
        V2222 = 0
        # V1111 = V1(11)1, V2111 = V2(11)1, ...
        # elements above are independent components of V matrix
    # To wrap lines. Contents: two-subscript matrix element calculation
    for i in range(N_GauFun):
        for j in range(N_GauFun):

            rAP = a2[j]*R_bond/(a1[i]+a2[j])
            # use R_{AP} represents a somewhat portion of R_bond
            rAP2 = rAP * rAP
            rBP2 = (R_bond-rAP)**2
            # and R_{BP} is the other portion of bond length

            S12 += gint.Sij_Gau(
                a1[i], 
                a2[j], 
                R2
                )*d1[i]*d2[j]
            T11 += gint.Tij_Gau(
                a1[i], 
                a1[j], 
                0
                )*d1[i]*d1[j]
            T12 += gint.Tij_Gau(
                a1[i], 
                a2[j], 
                R2
                )*d1[i]*d2[j]
            T22 += gint.Tij_Gau(
                a1[i], 
                a2[j], 
                0
                )*d1[i]*d2[j]

            Vne11_A += gint.Vneij_Gau(
                a1[i], 
                a1[j], 
                0, 
                0, 
                charge_A
                )*d1[i]*d1[j]
            Vne12_A += gint.Vneij_Gau(
                a1[i], 
                a2[j], 
                R2, 
                rAP2, 
                charge_A
                )*d1[i]*d2[j]
            Vne22_A += gint.Vneij_Gau(
                a2[i], 
                a2[j], 
                0, 
                R2, 
                charge_A
                )*d2[i]*d2[j]
            Vne11_B += gint.Vneij_Gau(
                a1[i], 
                a1[j], 
                0, 
                R2, 
                charge_B
                )*d1[i]*d1[j]
            Vne12_B += gint.Vneij_Gau(
                a1[i], 
                a2[j], 
                R2, 
                rBP2, 
                charge_B
                )*d1[i]*d2[j]
            Vne22_B += gint.Vneij_Gau(
                a2[i], 
                a2[j], 
                0, 
                0, 
                charge_B
                )*d2[i]*d2[j]
    # To wrap lines. Contents: four-subscript matrix element calculation
    for i in range(N_GauFun):
        for j in range(N_GauFun):
            for k in range(N_GauFun):
                for l in range(N_GauFun):

                    if i == 0 and j == 0 and k == 0 and l == 0:
                        print(
                            '           Four-center electron-electron integration activated. For some physical interpretation,\n'
                           +'           please see source code, where atomic basis functions used for expanding solution of canonical\n' 
                           +'           Hartree-Fock function are employed to understand matrix element of 4-dimensional Vee matrix, \n'
                           +'           that is, atomic basis function can be a viewpoint of electron localization.\n'
                            )
                    rAP = R_bond*a2[i]/(a2[i]+a1[j])
                    rBP = R_bond - rAP
                    rAQ = R_bond*a2[k]/(a2[k]+a1[l])
                    rBQ = R_bond - rAQ
                    rPQ = rAP-rAQ
                    rAP2 = rAP * rAP
                    rBP2 = rBP * rBP
                    rAQ2 = rAQ * rAQ
                    rBQ2 = rBQ * rBQ
                    rPQ2 = rPQ * rPQ

                    V1111 = V1111 + gint.Veeijkl_Gau(
                        a1[i], 
                        a1[j], 
                        a1[k], 
                        a1[l], 
                        0, 
                        0, 
                        0
                        )*d1[i]*d1[j]*d1[k]*d1[l]
                    V2111 = V2111 + gint.Veeijkl_Gau(
                        a2[i], 
                        a1[j], 
                        a1[k], 
                        a1[l], 
                        R2, 
                        0, 
                        rAP2
                        )*d2[i]*d1[j]*d1[k]*d1[l]
                    V2121 = V2121 + gint.Veeijkl_Gau(
                        a2[i], 
                        a1[j], 
                        a2[k], 
                        a1[l], 
                        R2, 
                        R2, 
                        rPQ2
                        )*d2[i]*d1[j]*d2[k]*d1[l]
                    V2211 = V2211 + gint.Veeijkl_Gau(
                        a2[i], 
                        a2[j], 
                        a1[k], 
                        a1[l], 
                        0, 
                        0, 
                        R2
                        )*d2[i]*d2[j]*d1[k]*d1[l]
                    V2221 = V2221 + gint.Veeijkl_Gau(
                        a2[i], 
                        a2[j], 
                        a2[k], 
                        a1[l], 
                        0, 
                        R2, 
                        rBQ2
                        )*d2[i]*d2[j]*d2[k]*d1[l]
                    # rBQ indicates integration (exchange) between k and l, stored as ijkl, (ij|kl)
                    V2222 = V2222 + gint.Veeijkl_Gau(
                        a2[i], 
                        a2[j], 
                        a2[k], 
                        a2[l], 
                        0, 
                        0, 
                        0
                        )*d2[i]*d2[j]*d2[k]*d2[l]
    # To wrap lines, use if True, no other meanings. Contents: form FULL Vee
    if True:
        Vee = np.zeros(shape = (2, 2, 2, 2))

        Vee[0][0][0][0] = V1111 # Columb

        # exchange format: (ii, jj), also means (electron1 motion, electron2 motion)
        Vee[1][0][0][0] = V2111 # electron 1 exchanges from spatial orbital 2 to 1
        Vee[0][1][0][0] = V2111 # Exchange, bi-directional
        Vee[0][0][1][0] = V2111 # electron 2 exchanges from spatial orbital 2 to 1
        Vee[0][0][0][1] = V2111 # Exchange, bi-directional

        Vee[1][0][1][0] = V2121 # electron 1 and 2 exchange from spatial orbital from 2 to 1
        Vee[0][1][0][1] = V2121 # Exchange, bi-directional
        Vee[0][1][1][0] = V2121 # electron 1 and 2 exchange from spatial orbital from 1 to 2 and 2 to 1 respectively
        Vee[1][0][0][1] = V2121 # Exchange, bi-directional

        Vee[1][1][0][0] = V2211 # Columb, electron 1 stays at spatial orbital 2, electron 2 stays at spatial orbital 1
        Vee[0][0][1][1] = V2211 # Columb, permute 2211

        Vee[1][1][1][0] = V2221 # electron 2 exchanges from spatial orbital 2 to 1
        Vee[1][1][0][1] = V2221 # Exchange, bi-directional
        Vee[1][0][1][1] = V2221 # electron 1 exchanges from spatial orbital 2 to 1
        Vee[0][1][1][1] = V2221 # Exchange, bi-directional

        Vee[1][1][1][1] = V2222 # Columb

        S = [[S11, S12[0]], [S12[0], S22]]
        # developer notes here: I dont know why T11 is (1,)-shaped while T is (2, 2, 1)-shaped
        # I dont know!
        T = [[T11[0], T12[0]], [T12[0], T22[0]]]
        Vne_A = [[Vne11_A[0], Vne12_A[0]], [Vne12_A[0], Vne22_A[0]]]
        Vne_B = [[Vne11_B[0], Vne12_B[0]], [Vne12_B[0], Vne22_B[0]]]
    # ------------------------------------------------------------------------------
    return [T, Vne_A, Vne_B, Vee, S]

def scf(N_GauFun, 
        R_bond, 
        zeta_A, 
        zeta_B, 
        charge_A, 
        charge_B, 
        conv_thre = 1E-6, 
        maxscf = 100,
        verbosity = 'debug'
        ):
    if verbosity == 'debug':
        print('SIMUPKGS | input parameters collection:\n'
             +'           Number of Gaussian function used for expanding atomic wfc:'+str(N_GauFun)+'\n'
             +'           He-H bond length in Angstrom: '+str(R_bond)+'\n'
             +'           Orbital shrink parameters A: '+str(zeta_A)+'\n'
             +'           Orbital shrink parameters B: '+str(zeta_B)+'\n'
             +'           Core charge A: '+str(charge_A)+'\n'
             +'           Core charge B: '+str(charge_B)+'\n'
             +'           convergence threshold: '+str(conv_thre)+'\n'
             +'           maximun number of SCF iterations: '+str(maxscf)
        )
    [Tij, Vne_A_ij, Vne_B_ij, Veeijlk, Sij] = iniIntg(
                                                      N_GauFun = N_GauFun,
                                                      R_bond = R_bond,
                                                      zeta_A = zeta_A,
                                                      zeta_B = zeta_B,
                                                      charge_A = charge_A,
                                                      charge_B = charge_B
                                                      )
    if verbosity == 'debug':
        print('SIMUPKGS | electronic integrals have been initialized successfully.\n'
             +'           shape information of matrices:\n'
             +'           Tij:      '+str(np.shape(Tij))+'\n'
             +'           Vne_A_ij: '+str(np.shape(Vne_A_ij))+'\n'
             +'           Vne_B_ij: '+str(np.shape(Vne_B_ij))+'\n'
             +'           Veeijlk:  '+str(np.shape(Veeijlk))+'\n'
             +'           Sij:      '+str(np.shape(Sij))+'\n'
        )
        print('SIMUPKGS_DEBUG_MODE | Tij = \n'+str(Tij))
        print('SIMUPKGS_DEBUG_MODE | Vne_A_ij = \n'+str(Vne_A_ij))
        print('SIMUPKGS_DEBUG_MODE | Vne_B_ij = \n'+str(Vne_B_ij))
        print('SIMUPKGS_DEBUG_MODE | Veeijlk = \n'+str(Veeijlk))
        print('SIMUPKGS_DEBUG_MODE | Sij = \n'+str(Sij))
    # initial guess, SCF_GUESS = ATOMIC
    Pij = np.zeros(shape = (2, 2))
    Cij = np.zeros_like(Pij)
    if verbosity == 'debug':
        print('SIMUPKGS | density matrix has been initialized successfully.')
    # forward propagation...
    Gij = hfprop.Gij(denMat = Pij, Vee = Veeijlk, mode = 'restricted')
    [Fij, Hij] = hfprop.hf_fp(
                              Tij = Tij,
                              natom = 2,
                              VneList = [Vne_A_ij, Vne_B_ij],
                              Gij = Gij,
                              printH = True,
                              )
    if verbosity == 'debug':
        print('SIMUPKGS_DEBUG_MODE | Gij = \n'+str(Gij))
        print('SIMUPKGS_DEBUG_MODE | Fij = \n'+str(Fij))
        print('SIMUPKGS_DEBUG_MODE | Hij = \n'+str(Hij))

    E = hfprop.hf_elecEner(Fij = Fij, Hij = Hij, Pij = Pij)
    iscf = 0
    E0 = 10086
    if verbosity == 'debug':
        print('SIMUPKGS | initial guess of energy E = '+str(E)+' Hartree')
    while (abs(E-E0) > conv_thre) and (iscf < maxscf):
        # solve FC=eSC
        [Cij, Pij] = hfprop.hf_bp(
                                  Fij = Fij, 
                                  Sij = Sij,
                                  printP = True,
                                  restricted = True
                                  )
        E0 = E
        E = hfprop.hf_elecEner(Fij = Fij, Hij = Hij, Pij = Pij)
        # forward propagation...
        Gij = hfprop.Gij(Pij, Vee = Veeijlk, mode = 'restricted')
        [Fij, Hij] = hfprop.hf_fp(
                                Tij = Tij,
                                natom = 2,
                                VneList = [Vne_A_ij, Vne_B_ij],
                                Gij = Gij,
                                printH = True,
                                )

        if verbosity != 'silent':
            print('SIMUPKGS | SCF iteration on-the-fly print:\n'
                 +'           Step '+str(iscf)+'\n'
                 +'           Energy = '+str(E)+' Hartree\n'
                 +'           Delta E = '+str(E-E0)+' Hartree\n'
                 +'           Convergence threshold = '+str(conv_thre)+'\n'
                 )
        if verbosity == 'debug':
            print('SIMUPKGS_DEBUG_MODE | Gij = \n'+str(Gij))
            print('SIMUPKGS_DEBUG_MODE | Cij = \n'+str(Cij))
            print('SIMUPKGS_DEBUG_MODE | Pij = \n'+str(Pij))
        iscf += 1

    
    output_dict = {}
    output_dict['coefficents'] = Cij
    output_dict['density_matrix'] = Pij
    output_dict['electronic_ener'] = E

    return output_dict
