# alternative version of cp2k_bs2csv.py

import numpy as np
import matplotlib.pyplot as plt

filename = 'graphene.bs'

f1tag = open(filename, 'r', encoding = 'utf-8')
line = f1tag.readline()
print('First line read: {}'.format(line))

words = [word for word in line.split(' ') if word != '']
print(words)
nSpecialPoints = int(words[3])
nPoints = int(words[6])
nBands = int(words[-2])

bands = np.zeros(shape = (nBands, nPoints))

for i in range(nSpecialPoints):

    print('Special points information: No. {}'.format(i+1))
    line = f1tag.readline()
    print(line)

readPoint = 1
iband = 0

while line:

    # enter the main read loop...
    line = f1tag.readline()
    if line.startswith('#  Point'):

        print('Read No. {} k-point, information:'.format(readPoint))
        print(line[0:-2])
        words = [word for word in line.split(' ') if word != '']
        print('band spin:          {}'.format(words[4][0]))
        print('k-space coordinate: {}, {}, {}'.format(words[-4], words[-3], words[-2]))
        print('k-point weight:     {}'.format(words[-1][0:-2]))
        readPoint += 1
    elif line.startswith('#   Band'):

        print('Read energy information...')
        iband = 0
    elif len(line) == 0:

        pass
    else:

        line = line[0:-2]
        numbers = [float(number) for number in line.split(' ') if number != '']
        bands[iband][readPoint-2] = numbers[1]
        iband += 1

f1tag.close()
print('*.bs file {} read complete! Draw bands...'.format(filename))

for idx_band in range(nBands):

    plt.plot(bands[idx_band][:], 'b-')
plt.show()
