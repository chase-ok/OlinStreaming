'''
Created on Mar 24, 2012

@author: ckernan
'''

from matplotlib import pyplot as plt
from matplotlib import animation
from streaming import analysis
import numpy as np

def plot2DVectorFields(fields, show=True):
    def plotIndividual(field):
        p, v = field
        plt.quiver(p[:, 0], p[:, 1], v[:, 0], v[:, 1])
    return _plotAnimation(fields, plotIndividual, show=show)

def plot2DFields(fields, show=True):
    def plotIndividual(field):
        (x, y), values = field
        plt.pcolor(values, vmin=0.0, vmax=1.0)
    return _plotAnimation(fields, plotIndividual, show=show)

def _plotAnimation(collection, func, show=True):
    figure = plt.figure()
    plot = animation.FuncAnimation(figure, func, collection, interval=5, repeat=False)
    
    if show: plt.show()
    return plot

def plotCorrelation(correlation, 
                    clear=True, 
                    points=True, 
                    show=True, 
                    line='b-', 
                    label='_nolegend_',
                    alpha=1.0,
                    **options):
    if clear: plt.clf()
    
    if points:
        for x, y in correlation.points:
            plt.plot(x, y, '.k', markersize=5, label='_nolegend_')
        
    plt.plot(correlation.x, correlation.meanY, line, 
             linewidth=1, 
             label=label,
             alpha=alpha)
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
        for c in correlations:
            plotCorrelation(c, 
                            clear=False, 
                            show=False, 
                            points=False, 
                            line=line,
                            label="_nolegend_",
                            alpha=0.4,
                            **options)
        
        x, mean = analysis.meanOfCorrelations(correlations)
        plt.plot(x, mean, line, linewidth=3, label=group)
        _applyOptions(options)
    
    plt.legend()
    if show: plt.show()
    
def saveCurrentPlot(fileName):
    plt.savefig(fileName)
    
def _applyOptions(options):
    basics = dict(title=plt.title, xlabel=plt.xlabel, ylabel=plt.ylabel,
                  xlim=plt.xlim, ylim=plt.ylim)
    for key, func in basics.iteritems():
        if key in options: func(options[key])
    
