header_r110 = ' [Molden Format]\n [Cell]\n 11.8274000 0.00000000 0.00000000\n 0.00000000 13.0889000 0.00000000\n 0.00000000 0.00000000 24.1003000\n'
header_r110 += ' [Nval]\n H 1\n C 4\n O 6\n Ti 12\n Pd 18\n [Atoms] AU\n'

header_a101 = ' [Molden Format]\n [Cell]\n 10.3034000 0.00000000 0.00000000\n 0.00000000 11.3902000 0.00000000\n 0.00000000 0.00000000 24.4100000\n'
header_a101 += ' [Nval]\n H 1\n C 4\n O 6\n Ti 12\n Pd 18\n [Atoms] AU\n'

def cp2k_molden_multiwfn(filein, header, fileout = ''):

    if len(fileout) < 1:

        fileout = filein[:30] + '-MTWFN.molden'
    filetag1 = open(filein, mode = 'r', encoding = 'utf-8')
    filetag2 = open(fileout, mode = 'a+', encoding = 'utf-8')

    filetag2.writelines(header)

    for iline in range(2):

        print('Skip line {} in file {} and wait to write new header:'.format(iline, filein))
        _prtline = filetag1.readline()
        print(_prtline)
    filetag2.writelines(filetag1.readlines())

    filetag1.close()
    filetag2.close()

a101FileList = [
    'a101.Pd-1-ace_pbe0_cdft_h_c_9.954-MOS-1_0.molden',
    'a101.Pd-1-ace-2H_pbe0_cdft_h_c_9.848-MOS-1_0.molden',
    'a101.Pd-1-H2_cdft_h_c_1.89-MOS-1_0.molden',
    'a101.Pd-2-ace_pbe0_cdft_h_c_9.874-MOS-1_0.molden',
    'a101.Pd-2-ace-2H_pbe0_cdft_h_c_9.802-MOS-1_0.molden',
    'a101.Pd-2-H2_cdft_h_c_1.88-MOS-1_0.molden'
]
r110FileList = [
    'r110.Pd-ace_pbe0_cdft_h_c_9.846-MOS-1_0.molden',
    'r110.Pd-ace-2H_pbe0_cdft_h_c_9.764-MOS-1_0.molden',
    'r110.Pd-H2_pbe0_cdft_h_c_1.836-MOS-1_0.molden'
]
for ifile in a101FileList:

    cp2k_molden_multiwfn(
        filein = ifile,
        header = header_a101
    )
for jfile in r110FileList:

    cp2k_molden_multiwfn(
        filein = jfile,
        header = header_r110
    )
