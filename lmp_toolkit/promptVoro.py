from lireVoronoi import lireVoronoi as lv
from lireVoronoi import VoroShape as vs

def allocateShape(freqList):
    
    VoroPrint = {}
    VoroPrint['PerfectIcosahedron'] = freqList[0]
    VoroPrint['ImperfectIcosahedron'] = freqList[1]
    VoroPrint['Tri-cappedTrigonalPrism'] = freqList[2]
    VoroPrint['Mono-cappedSquareArchimedeanAntiprism'] = freqList[3]
    return VoroPrint

header_list = [
    'Pd_300K',
    'Pd_650K',
    'PdO_300K',
    'PdO_650K'
]

for iheader in header_list:
    
    _, freq, _ = lv(
        header = iheader,
        nframe=40,
        atomtype = '2',
        motif = [
            vs('PerfectIcosahedron'), 
            vs('ImperfectIcosahedron'), 
            vs('Tri-cappedTrigonalPrism'), 
            vs('Mono-cappedSquareArchimedeanAntiprism')
            ],
        concat2file = False
    )

    print('Results: for {}, '.format(iheader))
    print(allocateShape(freq))
