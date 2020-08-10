import math

print('=------------------------------------------------------------=\n'
     +'|                                                            |\n'
     +'|      Independent subroutine of Free energy calculation     |\n'
     +'|                                                            |\n'
     +'=------------------------------------------------------------=')
print('FreeEner| vibarational information is needed, please remember.\n'
     +'FreeEner| molecules supported: linear structure, 2 or 3 atoms.\n'
     +'FreeEner| formulation used:\n'
     +'          G = H - TS\n'
     +'          H = U + pV, ideal gas, pV = nRT = NkT\n'
     +'          U = E_dft + ZPE + E_trans + E_vib + E_rot\n'
     +'          S = S_trans + S_rot + S_vib')

R = 8.314
N_A = 6.02E23
kb = R/N_A
h = 6.626E-34
#----------------------------------data input--------------------------------------
enerinfo = input('FreeEner| energy calculated by electronic method and unit of it are required,\n'
                +'          remember to seperate them with space (e.g. 20 kJ): ')
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
natom = int(input('FreeEner| number of atoms input: (2/3) '))
masslist = input('FreeEner| atomic mass is required, input them and seperate with space: ')
bondlist = input('FreeEner| bond length(s) is required, input them and seperate with space (in Angstrom): ')
freqlist = input('FreeEner| vibrational frequencies are required, in cm-1: ')
temp = float(input('FreeEner| temperature in Kelvin is required: '))
pressure = float(input('FreeEner| pressure in Pa is required: '))
outunit = input('FreeEner| output unit is required, [1] J·mol-1; [2] eV; [3] a.u. ')
#----------------------------------------------------------------------------------
masslist = [ float(mass) for mass in masslist.split(' ') if mass != '' ]
bondlist = [ float(bond) for bond in bondlist.split(' ') if bond != '' ]
freqlist = [ float(freq) for freq in freqlist.split(' ') if freq != '' ]

if natom != len(masslist) or natom != len(bondlist)+1:
     print('FreeEner| ***error*** data input is not consistent with itself. quit.'
          +'          number of atoms: '+str(natom)+'\n'
          +'          number of mass:  '+str(len(masslist))+'\n'
          +'          number of bonds: '+str(len(bondlist)))
     exit()

# create atomic coordinates in 1 dimension
coordlist = []
distlist = bondlist
distlist.insert(0, float(0))
coord = 0
massAll = 0
massCenter = 0

for iatom in range(natom):

     coord += distlist[iatom]
     coordlist.append(coord)
     # mass center
     massAll += masslist[iatom]
     massCenter += masslist[iatom]*coord

massCenter /= massAll

print('FreeEner| atomic coordinates in 1-d are: '+str(coordlist)+'\n'
     +'          mass center coordinate: '+str(massCenter))

mass_kg = [ mass/N_A/1000 for mass in masslist ]
massAll_kg = massAll/N_A/1000

# calculate distance from atoms to mass center
dist2center = []
rot_interia = 0
for iatom in range(natom):

     dist = abs(coordlist[iatom]-massCenter)*1E-10
     rot_interia += mass_kg[iatom]*dist**2

# vibrations
E_viblist = [ freq*100*3E8*h for freq in freqlist ]
ZPE = 0
E_vib = 0 # averaged vibrational energy, not obtained from 2*ZPE
q_vib = 1
S_vib = 0
for e_vib in E_viblist:

     ZPE += 0.5*e_vib
     e_vib_norm = e_vib/kb/temp
     q_vib *= 1/(1-math.exp(-e_vib_norm))
     E_vib += e_vib*math.exp(-e_vib_norm)
     S_vib += e_vib_norm*math.exp(-e_vib_norm)/(1-math.exp(-e_vib_norm))-math.log(1-math.exp(-e_vib_norm), math.pi)

E_vib /= q_vib

# partition functions
q_trans = ((2*math.pi*massAll_kg*kb*temp)/(h**2))**1.5
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
          R*temp/pressure*q_trans, math.e
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

g_elec = 1
# q_elec = g_elec*math.exp(-ener0/kb/temp) too large to calculate
A = -N_A*kb*temp*math.log(q_trans*q_rot*q_vib, math.e)+ener0

G = H-temp*S

if outunit == '1':
     outunit = ' J/mol\n'
elif outunit == '2':
     # eV
     ener0 /= 1.6E-19*6.02E23
     ZPE /= 1.6E-19*6.02E23
     E_rigid /= 1.6E-19*6.02E23
     E_vib /= 1.6E-19*6.02E23
     U /= 1.6E-19*6.02E23
     H /= 1.6E-19*6.02E23
     G /= 1.6E-19*6.02E23
     A /= 1.6E-19*6.02E23
     outunit = ' eV\n'
elif outunit == '3':
     # a.u.
     ener0 /= 1.6E-19*6.02E23*27.2113838565563
     ZPE /= 1.6E-19*6.02E23*27.2113838565563
     E_rigid /= 1.6E-19*6.02E23*27.2113838565563
     E_vib /= 1.6E-19*6.02E23*27.2113838565563
     U /= 1.6E-19*6.02E23*27.2113838565563
     H /= 1.6E-19*6.02E23*27.2113838565563
     G /= 1.6E-19*6.02E23*27.2113838565563
     A /= 1.6E-19*6.02E23*27.2113838565563
     outunit = ' a.u.\n'

print('\nFreeEner| energies information:\n'
     +'          electronic energy              = '+str(ener0)+str(outunit)
     +'          zero-point energy              = '+str(ZPE)+str(outunit)
     +'          translation + rotation         = '+str(E_rigid)+str(outunit)
     +'          vibrational energy             = '+str(E_vib)+str(outunit)
     +'          -------------------------------------------------------------\n'
     +'          total energy (internal energy) = '+str(U)+str(outunit)
     +'          total enthalpy                 = '+str(H)+str(outunit)
     )

print('FreeEner| entropies information:\n'
     +'          translation                    = '+str(S_trans)+' J/mol/K\n'
     +'          rotation                       = '+str(S_rot)+' J/mol/K\n'
     +'          vibration                      = '+str(S_vib)+' J/mol/K\n'
     +'          -------------------------------------------------------------\n'
     +'          total entropy                  = '+str(S)+' J/mol/K\n'
     )


print('FreeEner| free energy information:\n'
     +'          Gibbs free energy     = '+str(G)+str(outunit)
     +'          Helmholtz free energy = '+str(A)+str(outunit)
     +'          temperature           = '+str(temp)+' K\n'
     +'          pressure              = '+str(pressure)+' Pa'
     )
