import math
import numpy as np
from coordsCenter import centercoords
from distance import dist
import element_info as eleinfo

# to use this code, coordinate file must be specified in input file.
# also all parameters required: temperature, pressure, frequencies
# bond length is not necessary.

def freeEnergy(atomlist,
               enerinfo = '0 J',
               temp = 0,
               pressure = 101325,
               freqlist = [],
               unituse = '2'
              ):
    #--------------------constants definition-----------------------
    R = 8.314
    N_A = 6.02E23
    kb = R/N_A
    h = 6.626E-34
    #--------------------input preprocessing------------------------
    enerinfo = [ word for word in enerinfo.split(' ') if word != '' ]
    ener0 = float(enerinfo[0])
    enerUnit = enerinfo[1]
    if enerUnit == 'j' or enerUnit == 'J':
        ener0 = ener0
    elif enerUnit == 'kj' or enerUnit == 'kJ' or enerUnit == 'KJ' or enerUnit == 'Kj':
        ener0 *= 1E3
    elif enerUnit == 'eV' or enerUnit == 'EV' or enerUnit == 'ev':
        ener0 *= 1.6E-19
    elif enerUnit == 'a.u.' or enerUnit == 'au' or enerUnit == 'A.U.' or enerUnit == 'AU':
        ener0 *= 27.2113838565563*1.6E-19
    else:
        print('FreeEner| energy unit can not be recognized, please input unit such as J/kJ/eV/a.u., quit.')
        exit()
    natom = np.shape(atomlist)[0]
    masslist = []

    symbol_ref = 'Ar'
    for iatom in range(natom):

        if type(atomlist[iatom][0]) != type(symbol_ref):
            print('FreeEner| ***error*** invalid format of atomlist. quit.')
            exit()
        masslist.append(eleinfo.getElement(atomlist[iatom][0]))
    
    mass = np.sum(masslist)
    masslist_kg = [ mass/N_A/1000 for mass in masslist ]
    mass_kg = mass/N_A/1000
    #----------------rotation interia calculation---------------------------------
    coords = centercoords(atomlist, masspower=True, totmass=mass, masslist=masslist)
    # mass center corrected to 0, 0, 0.
    center = [0, 0, 0]
    r = dist(atomlist=coords, pointxyz=center) # unit in m

    rot_interia = 0
    for iatom in range(natom):

        rot_interia += masslist_kg[iatom]*r[iatom]**2
    
    #----------------------------------vibrations-----------------------------------
    E_viblist = [ freq*100*3E8*h for freq in freqlist ]
    ZPE = 0
    E_vib = 0 # averaged vibrational energy, not obtained from 2*ZPE
    q_vib = 1
    S_vib = 0
    for e_vib in E_viblist:

        ZPE += 0.5*e_vib
        if temp == 0:
            continue
        e_vib_norm = e_vib/kb/temp
        q_vib *= 1/(1-math.exp(-e_vib_norm))
        E_vib += e_vib*math.exp(-e_vib_norm)
        S_vib += e_vib_norm*math.exp(-e_vib_norm)/(1-math.exp(-e_vib_norm))-math.log(1-math.exp(-e_vib_norm))

    # exception of 0K input
    if temp == 0:

        print('FreeEner| 0K-temperature input detected! ZPE correction mode is activated, mass, bond length, pressure'
            +'          frequencies information will be discarded.')
        U = ener0 + ZPE # unit in J
        if unituse == '1':
            ener0 *= N_A
            U *= N_A
            unituse = ' J/mol\n'
        elif unituse == '2':
            ener0 /= 1.6E-19
            U /= 1.6E-19
            unituse = ' eV\n'
        elif unituse == '3':
            ener0 /= 1.6E-19*27.2113838565563
            U /= 1.6E-19*27.2113838565563
            unituse = ' a.u.\n'
        else:
            print('FreeEner| wrong output unit input, quit.')
            exit()
        
        print('FreeEner| ZPE correction mode output information:'
            +'          electronic energy    = '+str(ener0)+unituse
            +'          zero-point corrected = '+str(U)+unituse
            )
        exit()

    E_vib /= q_vib

    # partition functions
    q_trans = ((2*math.pi*mass_kg*kb*temp)/(h**2))**1.5
    q_rot = 8*math.pi**2*kb*temp*rot_interia/h
    if masslist[0]==masslist[-1]:
        q_rot /= 2

    print('FreeEner| partition functions information:\n'
        +'          translation q_trans (volume omitted) = '+str(q_trans)+'\n'
        +'          rotation q_rot                       = '+str(q_rot)+'\n'
        +'          vibration q_vib                      = '+str(q_vib)+'\n'
        +'          overall Q = q_trans*q_rot*q_vib      = '+str(q_trans*q_rot*q_vib))

    print('FreeEner| energy calculation starts. Ideal gas assumption is used:\n'
        +'          E_trans + E_rot = 5/2NkT')
    E_rigid = 5/2*kb*temp # E_rigid = E_trans + E_rot
    U = ener0 + ZPE + E_rigid + E_vib
    H = U + kb*temp

    # entropy calculation, ideal stastical mechanics is used.
    S_trans = R*(
        5/2+math.log(
            R*temp/pressure*q_trans
            )
        )
    S_rot = R*(
        3/2+math.log(
            math.sqrt(math.pi)/2*(8*math.pi**2*kb*temp/h**2)**1.5*math.sqrt(q_rot)
            )
        )
    S_vib = S_vib
    S = S_trans + S_rot + S_vib

    # J·mol-1
    ener0 *= N_A
    ZPE *= N_A
    E_rigid *= N_A
    E_vib *= N_A
    U *= N_A
    H *= N_A

    # g_elec = 1
    # q_elec = g_elec*math.exp(-ener0/kb/temp) too large to calculate
    A = -N_A*kb*temp*math.log(q_trans*q_rot*q_vib)+ener0

    G = H-temp*S

    if unituse == '1':
        unituse = ' J/mol\n'
    elif unituse == '2':
        # eV
        ener0 /= 1.6E-19*6.02E23
        ZPE /= 1.6E-19*6.02E23
        E_rigid /= 1.6E-19*6.02E23
        E_vib /= 1.6E-19*6.02E23
        U /= 1.6E-19*6.02E23
        H /= 1.6E-19*6.02E23
        G /= 1.6E-19*6.02E23
        A /= 1.6E-19*6.02E23
        unituse = ' eV\n'
    elif unituse == '3':
        # a.u.
        ener0 /= 1.6E-19*6.02E23*27.2113838565563
        ZPE /= 1.6E-19*6.02E23*27.2113838565563
        E_rigid /= 1.6E-19*6.02E23*27.2113838565563
        E_vib /= 1.6E-19*6.02E23*27.2113838565563
        U /= 1.6E-19*6.02E23*27.2113838565563
        H /= 1.6E-19*6.02E23*27.2113838565563
        G /= 1.6E-19*6.02E23*27.2113838565563
        A /= 1.6E-19*6.02E23*27.2113838565563
        unituse = ' a.u.\n'

    print('\nFreeEner| energies information:\n'
        +'          electronic energy              = '+str(ener0)+str(unituse)
        +'          zero-point energy              = '+str(ZPE)+str(unituse)
        +'          translation + rotation         = '+str(E_rigid)+str(unituse)
        +'          vibrational energy             = '+str(E_vib)+str(unituse)
        +'          -------------------------------------------------------------\n'
        +'          total energy (internal energy) = '+str(U)+str(unituse)
        +'          total enthalpy                 = '+str(H)+str(unituse)
        )

    print('FreeEner| entropies information:\n'
        +'          translation                    = '+str(S_trans)+' J/mol/K\n'
        +'          rotation                       = '+str(S_rot)+' J/mol/K\n'
        +'          vibration                      = '+str(S_vib)+' J/mol/K\n'
        +'          -------------------------------------------------------------\n'
        +'          total entropy                  = '+str(S)+' J/mol/K\n'
        )


    print('FreeEner| free energy information:\n'
        +'          Gibbs free energy     = '+str(G)+str(unituse)
        +'          Helmholtz free energy = '+str(A)+str(unituse)
        +'          temperature           = '+str(temp)+' K\n'
        +'          pressure              = '+str(pressure)+' Pa'
        )