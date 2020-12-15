import math

print('FreeEner| absorption / desorption free energy calculated under low pressure approximation is activated.\n'
      '          reference: Statistical mechanics——for theoretical chemistry (Chinese version), Chen Minbo, Science Press.')
print('FreeEner| formulation:\n'
     +'          K = q_prod/q_reac * exp(-Delta(E+ZPE)/kb/T)\n'
     +'          G = -kbTlnK = Delta(E+ZPE)-kbTln(q_prod/q_reac)')
R = 8.314
N_A = 6.02E23
kb = R/N_A
h = 6.626E-34

unitInp = input('FreeEner| unit of energies is required: [1] J; [2] J/mol; [3] eV; [4] a.u.: ')
# note that input type is str.
enerini = float(input('FreeEner| energy of isolated slab: '))
enermole = float(input('FreeEner| energy of isolated molecule: '))
enerfin = float(input('FreeEner| energy of composited slab (slab-gas molecule): '))

if unitInp == '1':
    enerini *= 1
    enermole *= 1
    enerfin *= 1
elif unitInp == '2':
    enerini /= N_A
    enermole /= N_A
    enerfin /= N_A
elif unitInp == '3':
    enerini *= 1.6E-19
    enermole *= 1.6E-19
    enerfin *= 1.6E-19
elif unitInp == '4':
    enerini *= 1.6E-19*27.2113838565563
    enermole *= 1.6E-19*27.2113838565563
    enerfin *= 1.6E-19*27.2113838565563
else:
    print('FreeEner| wrong requirement of unit. quit.')
    exit()

natom = int(input('FreeEner| number of atoms input: (2/3) '))
masslist = input('FreeEner| atomic mass is required, input them and seperate with space: ')
bondlist = input('FreeEner| bond length(s) is required, input them and seperate with space (in Angstrom): ')
freqlist = input('FreeEner| vibrational frequencies are required, in cm-1: ')
temp = float(input('FreeEner| temperature in Kelvin is required: '))
outunit = input('FreeEner| output unit is required, [1] J; [2] J/mol; [3] eV; [4] a.u.: ')

masslist = [ float(mass) for mass in masslist.split(' ') if mass != '' ]
bondlist = [ float(bond) for bond in bondlist.split(' ') if bond != '' ]
freqlist = [ float(freq) for freq in freqlist.split(' ') if freq != '' ]

if natom != len(masslist) or natom != len(bondlist)+1:
     print('FreeEner| ***error*** data input is not consistent with itself. quit.\n'
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
     +'          mass center coordinate: '+str(massCenter)+'\n')

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
Q = q_trans*q_rot*q_vib
print('FreeEner| partition functions information:\n'
     +'          translation q_trans (volume omitted) = '+str(q_trans)+'\n'
     +'          rotation q_rot                       = '+str(q_rot)+'\n'
     +'          vibration q_vib                      = '+str(q_vib)+'\n'
     +'          overall Q = q_trans*q_rot*q_vib      = '+str(Q)+'\n')

deltaE = enerfin - enerini - enermole
deltaU = deltaE + ZPE
K = math.exp(-deltaU/kb/temp)*1/Q
G = -kb*temp*math.log(K)

if outunit == '2':
    G *= N_A
    deltaE *= N_A
    deltaU *= N_A
    outunit = ' J/mol'
elif outunit == '3':
    G /= 1.6E-19
    deltaE /= 1.6E-19
    deltaU /= 1.6E-19
    outunit = ' eV'
elif outunit == '4':
    G /= 1.6E-19*27.2113838565563
    deltaE /= 1.6E-19*27.2113838565563
    deltaU /= 1.6E-19*27.2113838565563
    outunit = ' a.u.'
else:
    print('FreeEner| ***warning*** wrong output requirement on unit is detected. will use J as unit.')


print('FreeEner| free energy information:\n'
     +'          absorption energy change (E)         = '+str(deltaE)+str(outunit)+'\n'
     +'          absorption energy change (E+ZPE)     = '+str(deltaU)+str(outunit)+'\n'
     +'          absorption free energy change        = '+str(G)+str(outunit)+'\n'
     +'          desorption free energy change        = '+str(-G)+str(outunit))
