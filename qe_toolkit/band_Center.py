# center calculation for Quantumespresso

prefix = 'VLAB'
atom_range = [1, 2, 3, 4]
element = 'Au'
band = 'd'
band_wfc = 4
fermi = 3.9

print(
    '--------------------------Information collection--------------------------\n'
   +'prefix of QuantumEspresso job:                           '+prefix+'\n'
   +'atomic index of element that interested in:              '+str(atom_range)+'\n'
   +'element selected:                                        '+element+'\n'
   +'speific angular momentum atomic orbital selected:        '+band+'\n'
   +'atomic orbital corresponds to the index of wavefunction: '+str(band_wfc)+'\n'
   +'Fermi energy calculated from QuantumEspresso:            '+str(fermi)+' eV\n'
   +'-'*74
)

pdos_up = []
pdos_dw = []
ener = []
iloop = 0

print('BAND-CENTER| on-the-fly info. output...')
for iatom in atom_range:

    iloop += 1
    filename = prefix+'.'+'pdos_atm#'+str(iatom)+'('+element+')_wfc#'+str(band_wfc)+'('+band+')'
    print('           | info.: read new file: '+filename)

    fflag = open(filename)
    line = fflag.readline()
    
    if iloop == 1:
        while line:

            line = fflag.readline()
            words = line.split(' ')
            words = [word for word in words if word != '']
            if len(words) < 3:
                break
            ener.append(float(words[0]))
            pdos_up.append(float(words[1]))
            pdos_dw.append(float(words[2]))
            #print('ener read: {}, pdos_up_read: {}, pdos_dw_read: {}'.format(float(words[0]), float(words[1]), float(words[2])))
    else:
        dataidx = 0
        while line:

            line = fflag.readline()
            words = line.split(' ')
            words = [word for word in words if word != '']
            if len(words) < 3:
                break
            
            pdos_up[dataidx] += float(words[1])
            pdos_dw[dataidx] += float(words[2])
            dataidx += 1
            #print('ener read: {}, pdos_up_read: {}, pdos_dw_read: {}'.format(float(words[0]), float(words[1]), float(words[2])))

    fflag.close()

# calculate band center...
term1 = 0
term2 = 0


for idx in range(len(ener)):

    if ener[idx] <= fermi:
        term1 += (pdos_up[idx] + pdos_dw[idx])*ener[idx]
        term2 += (pdos_up[idx] + pdos_dw[idx])

center = term1/term2
print(
    'BAND-CENTER| results:\n'
   +'             Band center: {} eV\n'.format(center)
   +'             Band center respect to Fermi level: {} eV\n'.format(fermi-center)
    )

'''
import matplotlib.pyplot as plt

fig = plt.plot(ener, pdos_up)
plt.show()
'''
