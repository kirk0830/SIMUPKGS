from rdf2Sq import *

r, g = file2rdf(rdffile = 'test.rdf')
q, S = rdf2Sq(r = r, rdf = g, rho = 0.08, qmax = 10, nq = 100)

drawgr = False
drawsq = True
from matplotlib import pyplot as plt

def drawline(x, y, pattern = {'marker': '', 'markerfacecolor': '', 'line': '-', 'color': 'blue'}, 
             title = '', xlabel = '', ylabel = '', fontsize = 15, size = '', comment = ''):
    if size != '':
        fg = plt.figure(figsize=size)
    else:
        fg = plt.figure()
    plt.title(title, fontsize = fontsize)
    plt.xlabel(xlabel, fontsize = fontsize)
    plt.ylabel(ylabel, fontsize = fontsize)
    plt.xticks(fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    if comment != '':
        plt.text(3, 9, comment, fontsize = fontsize)
    if pattern['marker'] == '':
        plt.plot(
            x, y, pattern['line'], 
            color=pattern['color']
            )
    elif pattern['markerfacecolor'] == '':
        plt.plot(
            x, y, pattern['line'], 
            color=pattern['color'], 
            marker=pattern['marker']
            )
    else:
        plt.plot(
            x, y, pattern['line'], 
            color=pattern['color'], 
            marker=pattern['marker'], 
            markerfacecolor=pattern['markerfacecolor']
            )
    plt.show()

pattern = {
    'line': '-',
    'color': 'blue',
    'marker': 'o',
    'markerfacecolor': 'none'
}
if drawgr:
    drawline(r, g, pattern = pattern,
            title = 'Radical distribution function', 
            xlabel = 'Pd-Pd pair distance, $r$ (Angstrom)',
            ylabel = '$g(r)$'
            )
if drawsq:
    drawline(q, S, pattern = pattern, 
            title = 'Structural factor', 
            xlabel = 'Reciporal distance/wavenumber $q$ (Angstrom$^{-1}$)',
            ylabel = '$S(q)$',
            size = (8,6),
            comment = '$S(q) = 1+4\pi\\rho\int{dr\ r^2\\frac{sin(qr)}{qr}[g(r)-1]}$'
            )
