import re

def readin(inputfilename):

    keywordslist = {}
    try:
        inputfile = open(inputfilename, 'r+', encoding='utf-8')
        line = inputfile.readline()
        while line:
            oneline = []
            words = re.split(' |=|\n', line)
            for word in words:
                if word != '':
                    oneline.append(word)
            line = inputfile.readline()
            if len(oneline)>0:
                keyword = oneline[0]
                value = oneline[1]
                keywordslist[keyword] = value
        
        return keywordslist

    except FileNotFoundError:

        print('READ-IN| no such input file. file parser breaks down. QUIT.')
        exit()
