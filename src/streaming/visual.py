'''
Created on Mar 24, 2012

@author: ckernan
'''

from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

def plot2DField(fields, show=True):
    def plotIndividual(field):
        p, v = field
        return plt.quiver(p[:, 0], p[:, 1], v[:, 0], v[:, 1])
    
    figure = plt.figure()
    if isinstance(fields, list):
        images = [(plotIndividual(f),) for f in fields]
        plot = animation.ArtistAnimation(figure, images, 
                                         interval=10, repeat=False)
    else:
        plot = plotIndividual(fields)
    
    if show: plt.show()
    else: return plot
    
def plotCorrelation(correlations, 
                    clear=True, 
                    points=True, 
                    show=True, 
                    line='b-', 
                    label='_nolegend_', 
                    **options):
    if clear: plt.clf()
    
    sums = np.zeros_like(correlations[0][1])
    
    for dist, corr in correlations:
        if points: plt.plot(dist, corr, 'k.', markersize=5)
        sums += corr
        
    plt.plot(correlations[0][0], sums/len(correlations), line, 
             linewidth=2, 
             label=label)
    _applyOptions(options)
    
    if show: plt.show()
    
def overlayMeanCorrelations(groupedCorrelations, 
                            show=True,
                            lines=['b-', 'r-', 'c-', 'm-', 'y-', 'k-',
                                   'b--', 'r--', 'c--', 'm--', 'y--', 'k--'],
                            **options):
    if len(groupedCorrelations) > len(lines):
        raise ValueError("Not enough line styles!")
    
    plt.clf()
    plt.hold(True)
    
    for line, (group, correlations) in zip(lines, groupedCorrelations.items()):
        first = True
        for c in correlations:
            label = group if first else "_nolegend_"
            first = False
            
            plotCorrelation(c, 
                            clear=False, 
                            show=False, 
                            points=False, 
                            line=line,
                            label=label,
                            **options)
    
    plt.legend()
    if show: plt.show()
    
def saveCurrentPlot(fileName):
    plt.savefig(fileName)
    
def _applyOptions(options):
    basics = dict(title=plt.title, xlabel=plt.xlabel, ylabel=plt.ylabel,
                  xlim=plt.xlim, ylim=plt.ylim)
    for key, func in basics.iteritems():
        if key in options: func(options[key])
    
