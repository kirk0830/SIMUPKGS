import main_io as io
import coords_transform as ct
#-------------------------------------------------------------------------
prefacefile = open('preface.txt','r',encoding='utf-8')
line = prefacefile.readline().split('\n')[0]
while line:
    print(line)
    line = prefacefile.readline().split('\n')[0]
prefacefile.close()
#-------------------------------------------------------------------------

filename = input('\n\nSIMU_PKGS| waiting for type-in input filename...')

print('SIMU_PKGS| up to now, x-ray diffraction simulation, supercell duplication is implemented.')
if len(filename) < 1:
    print('SIMU_PKGS| ***error*** no input file information, QUIT.')
    exit()

keywordsdict = io.readInputFile(inputfilename=filename)
simulation = keywordsdict['run_type']
coordsfile = keywordsdict['coordinate_file']
ifsupercell = keywordsdict['supercell']
projectname = keywordsdict['project']

print(
    'SIMU_PKGS| simulation package initialization information:\n'
   +'           simulation job:             '+str(simulation)+'\n'
   +'           atomic coordinate filename: '+str(coordsfile)+'\n'
   +'           supercell duplication:      ['+str(ifsupercell)+']\n'
    )


print(
    'SIMU_PKGS| reading atomic coordinates...'
)
coords = io.coords_io.readcoords(filename=coordsfile, convert=True)

if ifsupercell == 'true' or ifsupercell == 'True' or ifsupercell == 'TRUE':
    dupx = int(keywordsdict['nx'])
    dupy = int(keywordsdict['ny'])
    dupz = int(keywordsdict['nz'])
    print(
        'SUPERCELL| duplication information:      '+str(dupx)+'x'+str(dupy)+'x'+str(dupz)+'\n'
        )
    cell_a = float(keywordsdict['cell_a'])
    cell_b = float(keywordsdict['cell_b'])
    cell_c = float(keywordsdict['cell_c'])
    cell_alpha = float(keywordsdict['cell_alpha'])
    cell_beta = float(keywordsdict['cell_beta'])
    cell_gamma = float(keywordsdict['cell_gamma'])

    print('SUPERCELL| for now, symmetry check is not implemented yet. <- 2020/07/28')

    supercoords = ct.cell_duplication(
                                        coords=coords,
                                        a=cell_a,
                                        b=cell_b,
                                        c=cell_c, 
                                        alpha=cell_alpha, 
                                        beta=cell_beta, 
                                        gamma=cell_gamma,
                                        nx=dupx,
                                        ny=dupy,
                                        nz=dupz)
    finalcoordsread = supercoords
else:
    finalcoordsread = coords

if simulation == 'xrd':

    print('SIMU_PKGS| X-RAY DIFFRACTION SIMULATION activated. XRD-specific information read-in.')

    centeredcoords = ct.centercoords(finalcoordsread)

    xrd_wavelength = float(keywordsdict['xrd_lambda'])
    xrd_maxangle = float(keywordsdict['xrd_maxangle'])
    xrd_scrdist = float(keywordsdict['xrd_scrz'])
    xrd_scrresol = float(keywordsdict['xrd_resolution'])

    xrd_angle_amp_flag = keywordsdict['angle_amp']
    if xrd_angle_amp_flag == 'true' or xrd_angle_amp_flag == 'True' or xrd_angle_amp_flag == 'TRUE':
        xrd_angle_amp_flag = True
        pass
    else:
        xrd_angle_amp_flag = False
        # angle parameters required here and will read-in internal database in future.
    
    import xrd_main
    xrd2dpattern = xrd_main.xrd(
                                wavelength=xrd_wavelength,
                                maxangle=xrd_maxangle,
                                dist2scr=xrd_scrdist,
                                scrresol=xrd_scrresol,
                                if_angle_resol=xrd_angle_amp_flag,
                                atomcoords=centeredcoords)

    # expected output is 2d intensity matrix.

    xrd_smoothflag = keywordsdict['smooth_2dpattern']
    if xrd_smoothflag == 'true' or xrd_smoothflag == 'True' or xrd_smoothflag == 'TRUE':
        xrd_smoothflag = True
        xrd_smooth_sigma = float(keywordsdict['smooth_sigma'])

        xrd_smooth_fast_flag = keywordsdict['fast_smooth']
        if xrd_smooth_fast_flag == 'true' or xrd_smooth_fast_flag == 'True' or xrd_smooth_fast_flag == 'TRUE':
            xrd_smooth_fast_flag = True
            xrd_smooth_fast_cutoff = int(keywordsdict['fast_smooth_cutoff'])
        else:
            xrd_smooth_fast_flag = False
        
        from pp_smooth import gau2d
        smth_xrd2d = gau2d(
                           dataInp=xrd2dpattern, 
                           sigma=xrd_smooth_sigma,
                           fast=xrd_smooth_fast_flag,
                           fast_cutoff=xrd_smooth_fast_cutoff
                           )
        xrd2dpattern = smth_xrd2d

    xrd_denoiseflag = keywordsdict['denoise_2dpattern']
    if xrd_denoiseflag == 'true' or xrd_denoiseflag == 'True' or xrd_denoiseflag == 'TRUE':
        pass
        # i will add-in denoise module later, this is not hard.
    
    ifplotxrd2d = keywordsdict['plotxrd2d']
    if ifplotxrd2d == 'true' or ifplotxrd2d == 'True' or ifplotxrd2d == 'TRUE':
        ifplotxrd2d = True

    ifpxrd = keywordsdict['pxrd']
    if ifpxrd == 'true' or ifpxrd == 'True' or ifpxrd == 'TRUE':
        ifpxrd = True
        from xrd_powder import pxrd
        pxrddata = pxrd(xrd2d=xrd2dpattern)

        ifplotpxrd = keywordsdict['plotpxrd']
        if ifplotpxrd == 'true' or ifplotpxrd == 'True' or ifplotpxrd == 'TRUE':
            ifplotpxrd = True

    if ifpxrd:
        if ifplotpxrd and ifplotxrd2d:
            import matplotlib.pyplot as plt
            plt.subplot(121)
            plt.imshow(xrd2dpattern) # additional: , cmap=plt.cm.gray
            plt.subplot(122)
            plt.plot(pxrddata)
            plt.show()
        elif ifplotpxrd and not (ifplotxrd2d):
            import matplotlib.pyplot as plt
            plt.plot(pxrddata)
            plt.show()
    elif not (ifpxrd) and ifplotxrd2d:
        import matplotlib.pyplot as plt
        plt.imshow(xrd2dpattern)
        plt.show()

elif simulation == 'mechanics':

    # i plan to add all kinds of molecular dynamics simulations here, different integration schemes like
    # frog leap, velocity-verlet, verlet and simple scheme will be supported.

    # also, some optimization like ionic optimization and cell optimization will be treated as 
    # NVT (T=0) and NPT (T=0)

    # optimizer is the THING that i am working on now.

    # the way energy calculated can be categorized into classical mechanics and quantum mechanics,
    # for what i major in, i will write classical and quantum in periodic system first.
    pass

elif simulation == 'supercell':
    
    io.coords_io.writecoords(finalcoordsread, projectname)

elif simulation == 'free_energy_mole':

    from fe import freeEner_isolated
    ener_elec = keywordsdict['ener0']
    ener_elec = ener_elec[0]+' '+ener_elec[1]
    temp = float(keywordsdict['temperature'])
    press = float(keywordsdict['pressure'])
    freqs = keywordsdict['frequencies']
    #freqs = [ float(freq) for freq in freqs.split(' ') if freq != '' ]
    freqs = [ float(freq) for freq in freqs ]
    unit = keywordsdict['outunit']
    freeEner_isolated(atomlist=finalcoordsread,
                      enerinfo=ener_elec,
                      temp=temp,
                      pressure=press,
                      freqlist=freqs,
                      unituse=unit
                      )

elif simulation == 'post-pros-q6':
    import pair
    q6_cutoff = float(keywordsdict['q6_cutoff'])
    print(pair.q6_order_parameter(coords=finalcoordsread, cutoff=q6_cutoff))
    
print('SIMU_PKG| simulation complete and terminates successfully.\n'
     +'SIMU_PKG| any comments will be appreciated, contact me at https://github.com/kirk0830/simu_pkg\n'
     +'SIMU_PKG| HAVE A NICE DAY, GOODBYE. :)')
