from pandas import ExcelFile as pd_ExcelFile
from re import split as re_split
from numpy import shape as np_shape

def absToolkit(searchWidth = 5, debugMode = False, indep_mode = False):
    filetag = pd_ExcelFile('lib.xlsx')
    sheetlist = filetag.sheet_names

    print('absToolkit| Welcome, absorption lib is activated. Currently only molecule absorbs on metal is supported,\n'
         +'            but I will update and add more available dataset in the future.\n'
         +'            Advantage of this program: you can just add any data into lib.xlsx without hestitation such as\n'
         +'            messing up order of data.\n'
         )

    print('absToolkit| usage: just do what program tells you to do, you will obtain what you want to get.\n'
         +'absToolkit| absToolkit is a module belonging to open source SIMUPKGS package, it provides custormize\n'
         +'            usage. Any data you trust and suppose to be useful in your future work can be added in\n'
         +'            database saved in the same directory of this program -> lib.xlsx\n'
         +'            to guarantee program will read your data properly, please read quotes in that file.\n'
         +'            or you can just send email to developer of this tool, lib will be updated at Github\n'
         +'            --->\n'
         +'            For more tools (source codes), visit https://github.com/kirk0830/SIMUPKGS\n\n'
         +'absToolkit| Developer notes:\n\n'
         +'           *Dipole-dipole coupling, see -> Chin. J. Catal., 2017, 38: 1473–1480, for distinguishing isolated\n'
         +'            atoms and cluster, one important feature is whether CO frequencies are loading-dependent.\n'
         +'            dipole-dipole coupling (actually dynamically correlation) will cause higher vibrational \n'
         +'            frequencies at higher loading. And this correlation will decrease intensity of overall absorbance.\n\n'
         +'           *Stark effect, which is Zeeman effect in electrostatic version, see -> J. Phys.: Conf. Ser. 117 012003.\n'
         +'            Stark effect is a field-dipole interaction that breaks degeneracy of M-CO vibration frequencies.\n\n'
         +'           *Failure of DFT prediction of CO\n'
         +'            as extraordinary delocalization yielded by DFT, a compensating DFT + molecular U on C 2p and O 2p\n'
         +'            should be used. Also under several circumstances, traditional functionals such as PBE, revPBE, rPBE\n'
         +'            will overestimate absorption energy, vdw-correction will even make problem worse. So do not blindly\n'
         +'            trust so-called first-principle.\n\n'
         +'           *Chemical environment distinguished by FWHM\n'
         +'            full width of half max, can be used as an qualitative indicator of homogenity of PGM_iso-CO sites\n'
         +'            such as Ir(CO)2 on various supports, FWHM ~5cm-1, similar with organmetallic complexes ~4cm-1, \n'
         +'            indicating an nearly isolated site, see J. Phys. Chem. Lett., 2016, 7, 3854–3860.\n'
         +'----------------------------------------------------------------------------------------------------\n'
         +'absToolkit| step 1: choose absorbate:'
        )
    if debugMode:
        exit()
    
    for isheet in range(len(sheetlist)):
        print('            ['+str(isheet)+'] '+str(sheetlist[isheet]))
    mole = input('absToolkit| (to directly quit program, enter 999) >> ')
    if mole == '999' or mole == '':
        exit()
    else:
        mole = int(mole)
    
    sheetObj = filetag.parse(sheet_name=mole)
    elementlist = set(sheetObj.iloc[:,0].values)
    print('absToolkit| for molecule you choose, there are '+str(sheetObj.shape[0])+' lines of data available.\n'
         +'----------------------------------------------------------------------------------------------------\n'
         +'absToolkit| step 2: type-in element symbol of your metal, currently supported metals: '+str([ielement for ielement in elementlist]))
    
    element = input('absToolkit| (to directly quit program, enter 999) >> ')
    if element == '999' or element == '':
        exit()

    print('----------------------------------------------------------------------------------------------------\n'
         +'absToolkit| step 3: type-in frequence in cm-1 you want to look up.')

    try:
        freqIn = float(input('absToolkit| >> '))
    except ValueError:
        print('absToolkit| ***error*** for no valid input, first five rows of dataset will be printed, quit.')
        # sheetOut = filetag.parse(sheet_name=mole).iloc[[0, 2, 3]] # to print out several specified lines
        print('----------------------------------------SEARCH RESULT-----------------------------------------------')
        print(sheetObj.iloc[:5,:6])
        exit()

    print('----------------------------------------------------------------------------------------------------\n'
         +'absToolkit| step 4: (press Enter directly or choose \'normal mode\' if you dont know what this step means)\n'
         +'            [1] normal mode (0)\n'
         +'            [2] very tight (0.0001)\n'
         +'            [3] tight (0.001)\n'
         +'            [4] default (0.01)\n'
         +'            [5] loose (0.1)\n'
         +'            [6] very loose (0.3)\n'
         +'            [7] awfully loose (0.5)\n'
         +'            [8] Do you feel lucky? (1.0)\n'
         +'            ----->>> I just salute to Gaussian software :)')
    rescale = input('absToolkit| >> ')
    if rescale == '' or rescale == '1':
        rescale = 0
    elif rescale == '2':
        rescale = 0.0001
    elif rescale == '3':
        rescale = 0.001
    elif rescale == '4':
        rescale = 0.01
    elif rescale == '5':
        rescale = 0.1
    elif rescale == '6':
        rescale = 0.3
    elif rescale == '7':
        rescale = 0.5
    elif rescale == '8':
        rescale = 1
    else:
        rescale = 0
    
    rescale_min = 1-rescale
    rescale_max = 1+rescale

    sheetdata = sheetObj.iloc[:].values # 2 dimensional list obtained.
    nline = np_shape(sheetdata)[0]
    
    outlist = []
    for iline in range(nline):

        if sheetdata[iline][0] == element:

            freqCheck = re_split(', ', str(sheetdata[iline][2]))
            freqCheck = [ float(freq) for freq in freqCheck ]
            if len(freqCheck) == 1:
                freqmin = freqCheck[0] - searchWidth
                freqmax = freqCheck[0] + searchWidth
                if freqIn >= freqmin and freqIn <= freqmax:
                    outlist.append(iline)
            elif len(freqCheck) == 2:
                freqmin = min(freqCheck[0], freqCheck[1])*rescale_min
                freqmax = max(freqCheck[0], freqCheck[1])*rescale_max
                if freqIn >= freqmin and freqIn <= freqmax:
                    outlist.append(iline)
    print('----------------------------------------SEARCH RESULT-----------------------------------------------')
    print(sheetObj.iloc[outlist, :6])
    #strout = sheetObj.iloc[0].values
    
    # use
    # .loc for label based indexing or
    # .iloc for positional indexing
    # .values for get data from it

    # developer notes:
    # this may be the first time that I use pandas pacakge... So I take note here for future use. Parameters without any explanations
    # are evaluated as parameters that wont be used by me. But you can still find detailed description at website:
    # https://pandas.pydata.org/pandas-docs/stable/reference

    # parameters of ExcelFile.parse(
    # sheet_name=0, : specify sheet you want to read, list is also supported. if not specified, all sheeets will be read in.
    # header=0, : determine which line pandas will start with, header = 0 for a full read
    # names=None, : specify names of column as you like, will overwrite results read-in.
    # index_col=None, : select one column as the first
    # usecols=None, : number of columns you want pandas to parse, there are also other usage, I omit here.
    # squeeze=False, 
    # converters=None, 
    # true_values=None, 
    # false_values=None, 
    # skiprows=None, 
    # nrows=None, : number of rows to skip at the beginning
    # na_values=None, 
    # parse_dates=False, 
    # date_parser=None, 
    # thousands=None, 
    # comment=None, 
    # skipfooter=0, 
    # convert_float=True, 
    # mangle_dupe_cols=True, 
    # **kwds)
    # 
    print('----------------------------------------------------------------------------------------------------')
    print('absToolkit| thank you, absorption toolkit will quit. Results above only the first two reference is listed.\n'
         +'            If you need more reference to cite, search the number emerges at head of line, it is the line number\n'
         +'            of datafile lib.xlsx.\n\n'
         +'            If you want to have contribution on this program, \n'
         +'            -> send your customized lib.xlsx to me ykhuang@dicp.ac.cn\n'
         +'            -> make comments and report issues at https://github.com/kirk0830/SIMUPKGS\n')
    if indep_mode:
        quitTag = input('absToolkit| press ENTER to quit. >> ')
    
absToolkit(indep_mode=True)