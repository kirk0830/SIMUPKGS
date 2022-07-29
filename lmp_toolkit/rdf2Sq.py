# convert Radical Distribution Function (RDF) to Structural factor S(q), q denotes wavenumber, has unit of L^{-1}
# formula used is from wikipedia,
# S(q) = 1+4\pi\rho\frac{1}{q}\int{dr r sin(qr) [g(r)-1]}
# therefore q \in [0, 2\pi/r] ?

import numpy as np

def file2rdf(rdffile):
    
    r = []
    rdf = []
    with open(file = rdffile, mode = 'r') as rdff:
        line = 'start'
        while line:
            line = rdff.readline()[0:-1]
            if len(line) > 1:
                words = line.split(' ')
                r.append(float(words[0]))
                rdf.append(float(words[1]))
    r = np.array(r, dtype = float)
    rdf = np.array(rdf, dtype = float)
    return r, rdf

def rdf2Sq(r, rdf, rho = 1., qmax = 'default', nq = 200):
    
    '''
    rdf_2dlist (2d list, float, float): a two dimensional list where the pair distance is stored in the first column, value of rdf is stored in the second one.\n
    rho (float): density of present liquid, has unit of 1/Angstrom^3
    '''
    rmin = np.min(r)
    dr = r[1] - r[0]
    S = []
    
    def sinx_on_x(x):
        
        if x == 0.:
            return 1.
        else:
            return np.sin(x)/x
    
    def __f__(q, idx):
        
        return 4*np.pi*rho*(r[idx])**2*sinx_on_x(q*r[idx])*(rdf[idx]-1)
    
    # will use Simpson method to integrate, 1/2f(0)*dx + f(1)*dx + f(2)*dx + ... + f(N-1)*dx + 1/2f(N)*dx
    # if rmin == 0.:
    #     rmin = 0.001
    if qmax == 'default':
        qmax = 2*np.pi/rmin
    elif type(qmax) == str:
        raise TypeError
    elif type(qmax) == float:
        pass
    elif type(qmax) == int:
        qmax = float(qmax)
    else:
        raise TypeError
    qmin = 2*np.pi/np.max(r)

    q = np.arange(start = qmin, stop = qmax, step = (qmax - qmin)/nq, dtype = float)
    print('g(r)2S(q)| q the wavenumber, reciporal coordinate is in range: [{}, {}], delta q = {}'.format(qmin, qmax, (qmax - qmin)/nq))
    
    for iq in q:
        Sq = 0
        for idx_r in range(len(r)):
            #
            if idx_r != (len(r) - 1):
                Sq += 0.5*dr*__f__(iq, idx_r) + 0.5*dr*__f__(iq, idx_r+1)
                
        S.append(Sq+1)
        
    return q, S
