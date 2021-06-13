# alternative version of cp2k_bs2csv.py

import numpy as np
import matplotlib.pyplot as plt

filename = 'graphene.bs'
plotBandsMode = 'balance'
# balance: nLUMO = nHOMO, all: all
emphasizeOcc = True

f1tag = open(filename, 'r', encoding = 'utf-8')
line = f1tag.readline()
print('First line read: {}'.format(line))

words = [word for word in line.split(' ') if word != '']
#print(words)
nSpecialPoints = int(words[3])
nPoints = int(words[6])
nBands = int(words[-2])
nOccBands = 0

#exit()
bands = np.zeros(shape = (nBands, nPoints))

for i in range(nSpecialPoints):

    print('Special points information: No. {}'.format(i+1))
    line = f1tag.readline()
    print(line)

readPoint = 1
iband = 0
ibandOcc = 0

while line:

    # enter the main read loop...
    line = f1tag.readline()
    if line.startswith('#  Point'):

        print('Read No. {} k-point, information:'.format(readPoint))
        words = [word for word in line.split(' ') if word != '']
        print('band spin:          {}'.format(words[4][0]))
        print('k-space coordinate: ({}, {}, {})'.format(words[-4], words[-3], words[-2]))
        print('k-point weight:     {}'.format(words[-1][0:-2]))
        readPoint += 1
    elif line.startswith('#   Band'):

        print('Read energy information of present k-point...\n')
        iband = 0
        ibandOcc = 0
    elif len(line) == 0:

        nOccBands = ibandOcc
        print('Reach EOF')
    else:

        line = line[0:-2]
        numbers = [float(number) for number in line.split(' ') if number != '']
        bands[iband][readPoint-2] = numbers[1]
        iband += 1
        if numbers[-1] > 0:

            ibandOcc += 1

f1tag.close()
print('*.bs file {} read complete! Draw bands...'.format(filename))

if plotBandsMode == 'balance':

    nBandsPlot = 2*nOccBands
else:

    nBandsPlot = nBands

if emphasizeOcc:

    for idx_band in range(nBandsPlot):

        if idx_band >= nOccBands:
            plt.plot(bands[idx_band][:], 'b-')
        else:
            plt.plot(bands[idx_band][:], 'r-')
else:

    for idx_band in range(nBandsPlot):
        plt.plot(bands[idx_band][:], 'b-')
    
plt.show()
