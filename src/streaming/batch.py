

from streaming import analysis, images, processing, visual, config
from matplotlib import pyplot as plt
import os
import pickle
from os import path

def fromDensities(**args):
    infos = []
    for density, folder in config.densities.items():
        infos.extend(images.loadSeqInfos(folder, group=density))
    return BatchRun(infos, **args)

class BatchRun(object):
    
    funcs = dict(particleDistance=analysis.particleDistance,
                 velocityCorrelation=analysis.velocityCorrelation,
                 directorCorrelation=analysis.directorCorrelation,
                 directorVelocityCorrelation=analysis.directorVelocityCorrelation,
                )
    
    def __init__(self, infos, outputFolder="../../output", minLength=50):
        self.infos = infos
        self.outputFolder = outputFolder
        self.minLength = 50
    
    def track(self, override=True):
        for info in self.infos:
            self._trackSeq(images.ImageSeq(info), override)
            
    def overlay(self, graph, **analysisOptions):
        results = self._performAnalysis(graph, **analysisOptions)
        results = self._organizeByGroup(results)
        
        if graph == "particleDistance":
            options = dict(title="Distance Correlation",
                           xlabel="Distance [microns]",
                           ylabel="# Particles [microns^-2]",
                           ylim=(0, 0.04))
        elif graph == "velocityCorrelation":
            options = dict(title="Velocity Correlation",
                           xlabel="Distance [microns]",
                           ylabel="Correlation [microns^-2]",
                           ylim=(-0.001, 0.02))
        elif graph == "directorCorrelation":
            options = dict(title="Director Correlation",
                           xlabel="Distance [microns]",
                           ylabel="Correlation [microns^-2]",
                           ylim=(0, 0.02))
        elif graph == "directorVelocityCorrelation":
            options = dict(title="Director-Velocity Correlation",
                           xlabel="Distance [microns]",
                           ylabel="Correlation [microns^-2]",
                           ylim=(0, 0.01))
        else:
            raise ValueError("Unknown graph: %s." % graph)
        
        visual.overlayMeanCorrelations(results, **options)        
    
    def _performAnalysis(self, analysisName, 
                         recalculate=False, 
                         skipMissing=True,
                         displayProgress=True):
        filePath = self._makeFilePath("analysis", analysisName + ".cache")
        func = self.funcs[analysisName]
        
        if recalculate or not path.exists(filePath):
            results = {}
        else:
            with open(filePath, 'rb') as f:
                results = pickle.load(f)
        
        for info in self.infos:
            if info in results: 
                continue
            
            if displayProgress:
                print "Analyzing " + info.uniqueName
            
            try:
                paths = self._loadPaths(info)
            except:
                if not skipMissing:
                    raise ValueError("Cannot find tracks for %s!"
                                     % info.uniqueName)
                
                if displayProgress:
                    print "No tracks found for %s, skipping." % info.uniqueName
                continue
            
            results[info] = func(paths)
        
        with open(filePath, 'wb') as f:
            pickle.dump(results, f)
        
        return results
    
    def _organizeByGroup(self, results):
        grouped = dict()
        for info, result in results.iteritems():
            grouped.setdefault(info.group, []).append(result)
        return grouped
    
    def _makeFilePath(self, subfolder, name):
        folder = path.join(self.outputFolder, subfolder)
        if not path.exists(folder): os.makedirs(folder)
        
        return path.join(folder, name)
        
    def _getPathsFile(self, info):
        return self._makeFilePath("paths", info.uniqueName + ".paths")
        
    def _trackSeq(self, seq, override):
        filePath = self._getPathsFile(seq.info)
        if override and path.exists(filePath): return
        
        paths = processing.trackBacteria(seq)
        with open(filePath, "wb") as f:
            pickle.dump(paths, f)
            
    def _loadPaths(self, info):
        with open(self._getPathsFile(info), "rb") as f:
            return pickle.load(f)
    
if __name__ == '__main__':
    batch = fromDensities()
    #batch.track(override=False)
    
    for func in batch.funcs:
        batch.overlay(func)
        visual.saveCurrentPlot(func + ".png")
    



