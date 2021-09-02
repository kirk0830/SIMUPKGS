# alternative version of cp2k_bs2csv.py

from sys import stdout
import numpy as np
import matplotlib.pyplot as plt
import sys
STDOUT_BKUP = sys.stdout

# SETTINGs
filename        = 'Layer-up-up-x-1-1.bs'
plotBandsMode   = 'balance' # plotBandsMode: how many bands to print
                            # [options]  balance: nLUMO = nHOMO
                            #            all: all
emphasizeOcc    = True      # emphasizeOcc: whether draw occupied bands in different color
verbosity       = 'high'  # verbosity: amount of information print
                            # [options]  silent: no information, only plot bands
                            #            low: silent + possible warning(s)
                            #            medium: low + I/O information
                            #            high: medium + formatted data output
                            #            debug: DO NOT USE IT
redirect        = 'no'      # redirect: whether print information to external file
                            # [options]  no: only print on screen
                            #            file: only print in file

filelog = filename[:-3] + '.log'
f1tag = open(filename, 'r', encoding = 'utf-8')
f2tag = open(filelog, 'a+', encoding = 'utf-8')
if redirect == 'file':
    sys.stdout = f2tag

line = f1tag.readline()
if verbosity != 'silent' and verbosity != 'low':
    print('First line read: {}'.format(line))

words = [word for word in line.split(' ') if word != '']
if verbosity == 'debug':
    print(words)
nSpecialPoints = int(words[3])
nPoints = int(words[6])
nBands = int(words[-2])
nOccBands = 0

nWarning = 0

bands = np.zeros(shape = (nBands, nPoints))

for i in range(nSpecialPoints):

    line = f1tag.readline()
    if verbosity == 'high' or verbosity == 'debug':
        print('Special points information: No. {}'.format(i+1))
        print(line)

readPoint = 1
iband = 0
ibandOcc = 0

while line:

    # enter the main read loop...
    line = f1tag.readline()
    if line.startswith('#  Point'):

        words = [word for word in line.split(' ') if word != '']
        if verbosity == 'high' or verbosity == 'debug':
            print('Read No. {} k-point, information:'.format(readPoint))
            print('band spin:          {}'.format(words[4][0]))
            print('k-space coordinate: ({}, {}, {})'.format(words[-4], words[-3], words[-2]))
            print('k-point weight:     {}'.format(words[-1][0:-2]))
        readPoint += 1
    elif line.startswith('#   Band'):

        if verbosity != 'silent' and verbosity != 'low':
            print('Read energy information of present k-point...\n')
        iband = 0
        ibandOcc = 0
    elif len(line) == 0:

        nOccBands = ibandOcc
        if verbosity != 'silent' and verbosity != 'low':
            print('Reach EOF')
    else:

        line = line[0:-2]
        numbers = [float(number) for number in line.split(' ') if number != '']
        bands[iband][readPoint-2] = numbers[1]
        iband += 1
        if numbers[-1] > 0:

            ibandOcc += 1

f1tag.close()
if verbosity != 'silent' and verbosity != 'low':
    print('*.bs file {} read complete! Draw bands...'.format(filename))

if plotBandsMode == 'balance':

    nBandsPlot = 2*nOccBands
    if nBandsPlot > nBands:
        if verbosity != 'silent' and verbosity != 'low':
            print(
                '---'*50
            +'\n***Warning*** Demand: \"Plotting bands in the way nHOMO = nLUMO\" can not be satisfied due to'
            +' number of LUMO supplied in *.bs file is insufficient.\n'
            +'---'*50
            )
            print(
                '---'*50
            +'\n***MOTION*** Adjust to plot ALL bands (change parameter \'plotBandsMode\' to \'all\')...\n'
            +'---'*50
            )
        nBandsPlot = nBands
        nWarning += 1
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

if verbosity != 'silent':
    print('Program summary: {} Warning(s)'.format(nWarning))
if redirect == 'file':
    sys.stdout = STDOUT_BKUP
f2tag.close()
