import re

def readInputFile(inputfilename):

    """
    read standard input file from external, all keywords intepretation will come soon as an isolated manual.
    Any comments can be added into input files as long as lines start with symbol '&'.\n
    # Input requirement\n
    filename, in string type\n
    # output description\n
    a dictionary type dataset containing all keywords and their values.
    """

    keywordslist = {}
    try:
        inputfile = open(inputfilename, 'r+', encoding='utf-8')
        line = inputfile.readline()
        while line:
            if line.startswith('&'):
                line = inputfile.readline()
                continue
            oneline = []
            words = re.split(' |=|\n', line)
            for word in words:
                if word != '':
                    oneline.append(word)
            line = inputfile.readline()
            if len(oneline)>0:
                keyword = oneline[0]
                value = oneline[1:]
                if len(value)==1:
                    value = oneline[1]
                keywordslist[keyword] = value
        
        return keywordslist

    except FileNotFoundError:

        print('READ-IN| no such input file. file parser breaks down. QUIT.')
        exit()

import coords_io