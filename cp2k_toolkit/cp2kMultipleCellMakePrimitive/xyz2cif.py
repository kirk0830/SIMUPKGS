# version 2022/06/15
# auther: kirk0830
import numpy as np

def xyz2cif(
    XYZfname, 
    cell_a = 10, 
    cell_b = 10, 
    cell_c = 10, 
    alpha = 90, 
    beta = 90, 
    gamma = 90,
    multiple_cell = 1,
    tx = 0.,
    ty = 0.,
    tz = 0.
    ):
    '''
    # XYZ -> CIF converter\n
    - version 2022/6/15\n
    XYZfname: full name of *.xyz file\n
    cell_a, cell_b and cell_c: a, b and c of cell in Angstrom\n
    alpha, beta, and gamma: alpha, beta and gamma of cell in deg\n
    multiple_cell: default as 1, enable common conversion from *.xyz to *.cif. If use CP2K MULTIPLE_CELL, and want to obtain
    primitive cell again after cell optimization, use multiple_cell = [A, B, C], where A, B and C are times multiplied by 
    primitive cell.
    '''
    if type(multiple_cell) == int:
        scale = [1., 1., 1.]
    else:
        scale = multiple_cell
    cell_a /= scale[0]
    cell_b /= scale[1]
    cell_c /= scale[2]
    a1 = cell_a
    a2 = 0.
    a3 = 0.
    b1 = cell_b*np.cos(gamma/180*np.pi)
    b2 = cell_b*np.sin(gamma/180*np.pi)
    b3 = 0.
    c1 = cell_c*np.cos(beta/180*np.pi)
    c2 = 1/b2*(cell_b*cell_c*np.cos(alpha/180*np.pi)-b1*c1)
    c3 = np.sqrt(cell_c**2 - c1**2 - c2**2)
    O = [
        [a1, b1, c1],
        [a2, b2, c2],
        [a3, b3, c3]
    ]
    O = np.array(O)
    O_inv = np.linalg.inv(O)
    CIFfname = ''
    if XYZfname.endswith('.xyz'):
        CIFfname = XYZfname[0:-3] + 'cif'
    else:
        CIFfname = XYZfname+'.cif'
        XYZfname += '.xyz'
    with open(file = CIFfname, mode = 'w', encoding='utf-8') as CIFf:
        CIFf.writelines(
            'data_image0\n'
           +'_cell_length_a       '+str(round(cell_a, ndigits=4))+'\n'
           +'_cell_length_b       '+str(round(cell_b, ndigits=4))+'\n'
           +'_cell_length_c       '+str(round(cell_c, ndigits=4))+'\n'
           +'_cell_angle_alpha    '+str(alpha)+'\n'
           +'_cell_angle_beta    '+str(beta)+'\n'
           +'_cell_angle_gamma    '+str(gamma)+'\n'
           +'\n'
           +'_symmetry_space_group_name_H-M    \"P 1\"\n'
           +'_symmetry_int_tables_number       1\n'
           +'\n'
           +'loop_\n'
           +'  _symmetry_equiv_pos_as_xyz\n'
           +'  \'x, y, z\'\n'
           +'\n'
           +'loop_\n'
           +'  _atom_site_label\n'
           +'  _atom_site_occupancy\n'
           +'  _atom_site_fract_x\n'
           +'  _atom_site_fract_y\n'
           +'  _atom_site_fract_z\n'
           +'  _atom_site_thermal_displace_type\n'
           +'  _atom_site_B_iso_or_equiv\n'
           +'  _atom_site_type_symbol\n'
           )
        # currently only orthogonal cell is supported.
        with open(file = XYZfname, mode = 'r', encoding='utf-8') as XYZf:

            line = 'start'
            idx_line = -1
            idx_atom = 0
            while line:

                line = XYZf.readline()[0:-1]
                idx_line += 1
                if idx_line == 0:
                    print('Number of atoms in present XYZ file:'+line)
                if idx_line > 1 and len(line) > 1:
                    # read XYZ data
                    words = line.split(' ')
                    words = [word for word in words if word != '']
                    r = [
                        [float(words[1]) + tx],
                        [float(words[2]) + ty],
                        [float(words[3]) + tz]
                        ]
                    r = np.array(r)
                    frac_r = np.matmul(O_inv, r)
                    idx_atom += 1
                    line2write = '  {:<{}}'.format(words[0]+str(idx_atom),9)
                    line2write += '1.0000 '
                    line2write += '%.5f' %(frac_r[0][0])
                    line2write += '  '
                    line2write += '%.5f' %(frac_r[1][0])
                    line2write += '  '
                    line2write += '%.5f' %(frac_r[2][0])
                    line2write += '  Biso   1.000  '
                    line2write += words[0]
                    line2write += '\n'
                    if type(multiple_cell) == int:
                        CIFf.writelines(line2write)
                    else:
                        if (abs(frac_r[0][0]) >= 1) or (abs(frac_r[1][0]) >= 1) or (abs(frac_r[2][0]) >= 1):
                            pass
                        else:
                            CIFf.writelines(line2write)
'''
from sys import argv

if len(argv) == 2:
    xyz2cif(
        XYZfname = argv[1]
    )
elif len(argv) == 5:
    xyz2cif(
        XYZfname = argv[1],
        cell_a = argv[2],
        cell_a = argv[3],
        cell_a = argv[4]
    )
elif len(argv) == 6:
    xyz2cif(
        XYZfname = argv[1],
        cell_a = argv[2],
        cell_a = argv[3],
        cell_a = argv[4],
        multiple_cell = argv[5]
    )
elif len(argv) == 8:
    xyz2cif(
        XYZfname = argv[1],
        cell_a = argv[2],
        cell_a = argv[3],
        cell_a = argv[4],
        alpha = argv[5],
        alpha = argv[6],
        alpha = argv[7]
    )
'''
# USAGE

a = 10.8681
b = 10.8681
c = 30.3527
alpha = 90
beta = 90
gamma = 120
multiple_cell = 1

xyz2cif('Rh1.YSZ-3M-site1-pos-fin.xyz', a, b, c, alpha, beta, gamma)

# Contents of file graphite-pos-fin.xyz
'''
     200
 i =       10, E =     -1139.5619820685
  C         0.0000144063       -0.0000079179        0.0000000219
  C        -0.0000333096        1.4249121074        0.0000000522
  C        -0.0000234552        0.0000213071        3.3093801073
  C        -1.2339670668        0.7124387438        3.3093806926
  C         2.4679964087       -0.0000072362        0.0000000935
  C         2.4679454766        1.4249104098       -0.0000003459
  C         2.4679600313        0.0000223568        3.3093805528
  C         1.2340179221        0.7124398925        3.3093804749
  C         4.9359747717       -0.0000085479       -0.0000004462
  C         4.9359251233        1.4249100202        0.0000000278
  C         4.9359420699        0.0000238348        3.3093805670
  C         3.7020003096        0.7124413451        3.3093804716
  C         7.4039572658       -0.0000095093        0.0000001583
  C         7.4039073310        1.4249100632        0.0000002041
  C         7.4039239462        0.0000231210        3.3093804771
  C         6.1699825952        0.7124416195        3.3093804185
  C         9.8719398075       -0.0000090615        0.0000002307
  C         9.8718911103        1.4249107076        0.0000003195
  C         9.8719022574        0.0000209770        3.3093808042
  C         8.6379633329        0.7124403758        3.3093805605
  C        -1.2339764955        2.1373282825       -0.0000001390
  C        -1.2340242185        3.5622487179       -0.0000001514
  C        -1.2340151562        2.1373578285        3.3093808686
  C        -2.4679580079        2.8497707723        3.3093808177
  C         1.2340070159        2.1373289648       -0.0000001886
  C         1.2339547293        3.5622472507        0.0000001268
  C         1.2339688122        2.1373589239        3.3093808584
  C         0.0000264276        2.8497722932        3.3093803419
  C         3.7019845291        2.1373269396        0.0000003598
  C         3.7019353071        3.5622456410        0.0000003058
  C         3.7019515973        2.1373594016        3.3093802531
  C         2.4680101752        2.8497729811        3.3093802640
  C         6.1699655250        2.1373268010        0.0000004236
  C         6.1699164418        3.5622462256        0.0000000772
  C         6.1699324599        2.1373593328        3.3093805572
  C         4.9359916978        2.8497728598        3.3093806996
  C         8.6379490903        2.1373274514       -0.0000003732
  C         8.6379000108        3.5622476815       -0.0000000312
  C         8.6379110539        2.1373572667        3.3093806378
  C         7.4039722389        2.8497720192        3.3093807414
  C        -2.4679674962        4.2746651072        0.0000003248
  C        -2.4680159237        5.6995851543        0.0000000789
  C        -2.4680062706        4.2746890593        3.3093805592
  C        -3.7019496095        4.9871041364        3.3093806164
  C         0.0000158728        4.2746656460       -0.0000000685
  C        -0.0000358683        5.6995832657       -0.0000000102
  C        -0.0000225488        4.2746905364        3.3093804261
  C        -1.2339649751        4.9871057810        3.3093806274
  C         2.4679943343        4.2746628444       -0.0000000772
  C         2.4679441065        5.6995814026       -0.0000003499
  C         2.4679603374        4.2746907114        3.3093807047
  C         1.2340183449        4.9871070054        3.3093808922
  C         4.9359752620        4.2746627834       -0.0000001557
  C         4.9359264913        5.6995817273       -0.0000002487
  C         4.9359421487        4.2746894926        3.3093807930
  C         3.7020008910        4.9871059362        3.3093803658
  C         7.4039576346        4.2746638054        0.0000000954
  C         7.4039092788        5.6995833040        0.0000000383
  C         7.4039198112        4.2746881374        3.3093808314
  C         6.1699810250        4.9871051889        3.3093807494
  C        -3.7019593314        6.4120010575        0.0000000205
  C        -3.7020073037        7.8369193827       -0.0000001937
  C        -3.7019975242        6.4120233759        3.3093806062
  C        -4.9359399095        7.1244395129        3.3093806371
  C        -1.2339753223        6.4120006531       -0.0000002215
  C        -1.2340281112        7.8369174035       -0.0000000732
  C        -1.2340138963        6.4120249010        3.3093808207
  C        -2.4679561511        7.1244409488        3.3093808003
  C         1.2340031336        6.4119986980        0.0000000441
  C         1.2339530142        7.8369158471       -0.0000000022
  C         1.2339690386        6.4120255484        3.3093805977
  C         0.0000274581        7.1244420582        3.3093802286
  C         3.7019844018        6.4119976951        0.0000002126
  C         3.7019349346        7.8369158578        0.0000002004
  C         3.7019507910        6.4120242440        3.3093806496
  C         2.4680094755        7.1244420274        3.3093804891
  C         6.1699673676        6.4119988545        0.0000002939
  C         6.1699188951        7.8369169233       -0.0000001048
  C         6.1699298521        6.4120222681        3.3093803837
  C         4.9359908118        7.1244403703        3.3093808556
  C        -4.9359497350        8.5493318195       -0.0000001442
  C        -4.9359975744        9.9742499207        0.0000000736
  C        -4.9359873604        8.5493592582        3.3093808603
  C        -6.1699299769        9.2617759241        3.3093805633
  C        -2.4679671505        8.5493328406        0.0000001201
  C        -2.4680192614        9.9742494615        0.0000002096
  C        -2.4680051304        8.5493610394        3.3093801194
  C        -3.7019466787        9.2617776371        3.3093802503
  C         0.0000119464        8.5493302086        0.0000002103
  C        -0.0000382206        9.9742481565        0.0000000397
  C        -0.0000216953        8.5493617112        3.3093804593
  C        -1.2339635344        9.2617788831        3.3093804066
  C         2.4679933162        8.5493295685       -0.0000002503
  C         2.4679442482        9.9742475522       -0.0000000544
  C         2.4679599344        8.5493610182        3.3093806575
  C         1.2340190857        9.2617789353        3.3093808773
  C         4.9359768170        8.5493301524       -0.0000002218
  C         4.9359276060        9.9742489726       -0.0000003613
  C         4.9359387792        8.5493588728        3.3093806446
  C         3.7019995891        9.2617777203        3.3093803736
  C         0.0000130252       -0.0000077080        6.6187609822
  C        -0.0000351528        1.4249129020        6.6187609448
  C        -0.0000234517        0.0000210873        9.9281413046
  C        -1.2339664138        0.7124382910        9.9281411496
  C         2.4679960032       -0.0000069178        6.6187610974
  C         2.4679449263        1.4249116501        6.6187608814
  C         2.4679592106        0.0000223795        9.9281413157
  C         1.2340173481        0.7124398213        9.9281413453
  C         4.9359759159       -0.0000085055        6.6187605697
  C         4.9359267468        1.4249111646        6.6187607198
  C         4.9359416035        0.0000232906        9.9281414044
  C         3.7019994730        0.7124412755        9.9281415041
  C         7.4039581888       -0.0000089475        6.6187610250
  C         7.4039083232        1.4249108717        6.6187607990
  C         7.4039245639        0.0000222673        9.9281411902
  C         6.1699826318        0.7124410376        9.9281411533
  C         9.8719395061       -0.0000082887        6.6187608046
  C         9.8718905632        1.4249114882        6.6187609855
  C         9.8719031720        0.0000203733        9.9281412377
  C         8.6379642704        0.7124396597        9.9281411812
  C        -1.2339779791        2.1373292405        6.6187610145
  C        -1.2340255102        3.5622490046        6.6187609032
  C        -1.2340143336        2.1373573476        9.9281411709
  C        -2.4679570279        2.8497704443        9.9281410758
  C         1.2340055933        2.1373301549        6.6187608067
  C         1.2339540811        3.5622484010        6.6187609370
  C         1.2339684451        2.1373587944        9.9281412086
  C         0.0000266971        2.8497720798        9.9281411387
  C         3.7019851642        2.1373281500        6.6187608750
  C         3.7019359040        3.5622469577        6.6187608151
  C         3.7019506308        2.1373596738        9.9281412748
  C         2.4680093277        2.8497732547        9.9281414523
  C         6.1699672362        2.1373283520        6.6187608457
  C         6.1699178281        3.5622473271        6.6187609055
  C         6.1699321734        2.1373594382        9.9281413257
  C         4.9359908702        2.8497734808        9.9281412509
  C         8.6379492445        2.1373278452        6.6187607316
  C         8.6378999500        3.5622474094        6.6187607634
  C         8.6379118147        2.1373569458        9.9281411344
  C         7.4039725484        2.8497723122        9.9281409297
  C        -2.4679683676        4.2746644756        6.6187609928
  C        -2.4680165529        5.6995835973        6.6187609076
  C        -2.4680055544        4.2746890748        9.9281410854
  C        -3.7019490960        4.9871048638        9.9281410870
  C         0.0000147945        4.2746655117        6.6187609560
  C        -0.0000364220        5.6995823809        6.6187610268
  C        -0.0000222908        4.2746902759        9.9281414655
  C        -1.2339643642        4.9871056456        9.9281411078
  C         2.4679944163        4.2746641678        6.6187609863
  C         2.4679444037        5.6995822173        6.6187609698
  C         2.4679596897        4.2746910316        9.9281413060
  C         1.2340181546        4.9871068845        9.9281411476
  C         4.9359763271        4.2746633207        6.6187609116
  C         4.9359271363        5.6995815965        6.6187607913
  C         4.9359412844        4.2746905338        9.9281410503
  C         3.7019999860        4.9871068207        9.9281412098
  C         7.4039583416        4.2746637407        6.6187605358
  C         7.4039090244        5.6995821927        6.6187607436
  C         7.4039199476        4.2746890072        9.9281411854
  C         6.1699805860        4.9871062408        9.9281411228
  C        -3.7019599883        6.4119987023        6.6187610283
  C        -3.7020081289        7.8369168563        6.6187607630
  C        -3.7019974917        6.4120242069        9.9281410509
  C        -4.9359402761        7.1244403986        9.9281411328
  C        -1.2339757438        6.4119991077        6.6187607884
  C        -1.2340277597        7.8369158734        6.6187607986
  C        -1.2340134200        6.4120247291        9.9281412250
  C        -2.4679559035        7.1244410709        9.9281414019
  C         1.2340029525        6.4119982234        6.6187610496
  C         1.2339538027        7.8369152309        6.6187608356
  C         1.2339692235        6.4120251488        9.9281410617
  C         0.0000279436        7.1244414860        9.9281412275
  C         3.7019851033        6.4119977109        6.6187607236
  C         3.7019352559        7.8369156102        6.6187609916
  C         3.7019503610        6.4120248856        9.9281411617
  C         2.4680094067        7.1244418878        9.9281413798
  C         6.1699674141        6.4119982141        6.6187610442
  C         6.1699183097        7.8369162036        6.6187607959
  C         6.1699293833        6.4120232895        9.9281412986
  C         4.9359903991        7.1244410250        9.9281411997
  C        -4.9359504032        8.5493306469        6.6187607238
  C        -4.9359984649        9.9742496441        6.6187609249
  C        -4.9359878144        8.5493598077        9.9281412178
  C        -6.1699301717        9.2617761137        9.9281415761
  C        -2.4679675382        8.5493307973        6.6187610414
  C        -2.4680191129        9.9742486616        6.6187609751
  C        -2.4680054654        8.5493609519        9.9281413245
  C        -3.7019472460        9.2617777910        9.9281414753
  C         0.0000128996        8.5493293100        6.6187609442
  C        -0.0000368957        9.9742478551        6.6187608397
  C        -0.0000213553        8.5493609436        9.9281414802
  C        -1.2339636257        9.2617781728        9.9281413920
  C         2.4679939712        8.5493295134        6.6187607157
  C         2.4679447321        9.9742479942        6.6187609047
  C         2.4679604623        8.5493604745        9.9281413942
  C         1.2340196745        9.2617780068        9.9281410935
  C         4.9359763535        8.5493298876        6.6187612756
  C         4.9359266103        9.9742490050        6.6187607542
  C         4.9359389013        8.5493590553        9.9281408804
  C         3.7020001167        9.2617773132        9.9281410914
'''