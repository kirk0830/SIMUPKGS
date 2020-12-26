# subroutine for read-in basis
# and convert to normalized form
import re
basisInfoFileTag = open(file = 'basis.dat', encoding = 'utf-8')

def getBasis(basisName = 'STO-1G', angular = 0, verbosity = 'low'):

    """
    Gaussian-type basis functions pop-up function\n
    basisName: name of basis, e.g., STO-1G, etc.\n
    angular: integar required, for s orbital, l = 0, d, l = 2, etc.
    """
    coeff = []
    expoCoeff = []
    line = 'Starting Line 0'
    while line:

        line = basisInfoFileTag.readline()
        if line.startswith('#'):
            # comment line, skip
            continue
        if line.startswith('$'):
            # basis starting line read-in
            words = re.split(pattern = '#| |\$|\n', string = line)
            basisReadIn = words[1]
            if basisReadIn == basisName:
                # if find the correct basis name
                print('SIMUPKGS | Basis information found, name of basis required: '+basisName)
                # start to read basis information...
                line = basisInfoFileTag.readline() # standard format: 1 0.27095 0 0 0
                if verbosity == 'debug':
                    print('BASIS_DEBUG_MODE | line read: '+str(line))
                while not line.startswith('$'):
                # stop reading in case of start to read a new section of basis...

                    words = re.split(pattern = ' |\n', string = line)
                    words = [float(word) for word in words if word != ''] # standard format: [1, 0.27095, 0, 0, 0]
                    if verbosity == 'debug':
                        print('BASIS_DEBUG_MODE | line split: '+str(words))
                    # check angular momentum
                    if (words[-1]+words[-2]+words[-3]) != angular:
                        # read the next line, or, if in the future samely-named basis are stored in different
                        # section, I will change this to "continue" directly to jump out present section and
                        # find another section that starts with the basis name demanded...
                        pass
                    else:
                        # read-in!
                        coeff.append(words[0])
                        expoCoeff.append(words[1])
                        # and read new lines directly...

                    line = basisInfoFileTag.readline()
    
    return [coeff, expoCoeff]
