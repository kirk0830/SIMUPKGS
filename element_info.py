# to find information of elements, this module should be import like:
# import element_info, but not from element_info import getElement,
# I will create a dictionary in this file.

mass_lib = {}
mass_lib['H'] = 1.008
mass_lib['He'] = 4.0026
mass_lib['Li'] = 6.94
mass_lib['Be'] = 9.0122
mass_lib['B'] = 10.81
mass_lib['C'] = 12.011
mass_lib['N'] = 14.007
mass_lib['O'] = 15.999
mass_lib['F'] = 18.998
mass_lib['Ne'] = 20.180
mass_lib['Na'] = 22.990
mass_lib['Mg'] = 24.305
mass_lib['Al'] = 26.982
mass_lib['Si'] = 28.085
mass_lib['P'] = 30.974
mass_lib['S'] = 32.06
mass_lib['Cl'] = 35.45
mass_lib['Ar'] = 39.948
mass_lib['K'] = 39.098
mass_lib['Ca'] = 40.078
mass_lib['Sc'] = 44.956
mass_lib['Ti'] = 47.867
mass_lib['V'] = 50.942
mass_lib['Cr'] = 51.996
mass_lib['Mn'] = 54.938
mass_lib['Fe'] = 55.845
mass_lib['Co'] = 58.933
mass_lib['Ni'] = 58.693
mass_lib['Cu'] = 64.546
mass_lib['Zn'] = 65.38
mass_lib['Ga'] = 69.723
mass_lib['Ge'] = 72.630
mass_lib['As'] = 74.922
mass_lib['Se'] = 78.971
mass_lib['Br'] = 79.904
mass_lib['Kr'] = 83.798

mass_lib['Ru'] = 101.07
mass_lib['Rh'] = 102.91
mass_lib['Pd'] = 106.42
mass_lib['Ag'] = 107.87
mass_lib['Ir'] = 192.22
mass_lib['Pt'] = 195.08
mass_lib['Au'] = 196.97

def getElement(symbol):

    try:
        return mass_lib[symbol]
    except KeyError:
        print('Element| ***error*** element information can not find! quit.')
        exit()