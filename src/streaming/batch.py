

from streaming import analysis, images, processing, visual, config
from matplotlib import pyplot as plt
import os
import numpy as np
import pickle
from os import path
from multiprocessing import Manager, Lock, Pool

def fromDensities(**args):
    infos = []
    for density, folder in config.densities.items():
        infos.extend(images.loadSeqInfos(folder, group=density))
    return BatchRun(infos, **args)

ANALYSIS_FUNCTIONS = \
    dict(particleDistance=analysis.particleDistance,
         velocityCorrelation=analysis.velocityCorrelation,
         directorCorrelation=analysis.directorCorrelation,
         directorVelocityCorrelation=analysis.directorVelocityCorrelation)

class BatchRun(object):
    
    funcs = ANALYSIS_FUNCTIONS
    
    def __init__(self, infos, 
                 outputFolder="../output", 
                 minLength=50, 
                 poolSize=4):
        self.infos = infos
        self.outputFolder = outputFolder
        self.minLength = 50
        self.poolSize = 4
    
    def track(self, recalculate=False):
        pool = Pool(self.poolSize)
        pool.map(_trackSeq, ((info, self._getPathsFile(info), recalculate) 
                             for info in self.infos))
            
    def overlay(self, graph, **analysisOptions):
        results = self._performAnalysis(graph, **analysisOptions)
        results = self._organizeByGroup(results)
        
        if graph == "particleDistance":
            options = dict(title="Distance Correlation",
                           xlabel="Distance [microns]",
                           ylabel="# Particles [microns^-2]",
                           ylim=(0, 0.15))
        elif graph == "velocityCorrelation":
            options = dict(title="Velocity Correlation",
                           xlabel="Distance [microns]",
                           ylabel="Correlation [microns^-2]",
                           ylim=(-0.001, 0.03))
        elif graph == "directorCorrelation":
            options = dict(title="Director Correlation",
                           xlabel="Distance [microns]",
                           ylabel="Correlation [microns^-2]",
                           ylim=(0, 0.03))
        elif graph == "directorVelocityCorrelation":
            options = dict(title="Director-Velocity Correlation",
                           xlabel="Distance [microns]",
                           ylabel="Correlation [microns^-2]",
                           ylim=(0, 0.03))
        else:
            raise ValueError("Unknown graph: %s." % graph)
        
        visual.overlayMeanCorrelations(results, show=False, **options)

    def dumpCorrelationToMatlab(self, correlation, fileName, **analysisOptions):
        results = self._performAnalysis(correlation, **analysisOptions)
        results = self._organizeByGroup(results)

        from mlabwrap import mlab

        groups = []
        for group, correlations in results.iteritems():

            wrapped = []
            for i, data in enumerate(correlations):
                info = _infoToMatlab(data.info, mlab)
                points = np.array(data.points)
                wrapped.append(mlab.struct("info", info, 
                                           "points", points, 
                                           "name", str(data.name)))

            groups.append(mlab.struct("group", group, 
                                      "correlations", _toStructArray(wrapped)))

        mlab._set(correlation, _toStructArray(groups))
        mlab.save(fileName, correlation)

    def dumpPathsToMatlab(self, fileName, **trackOptions):
        self.track(**trackOptions)

        from mlabwrap import mlab

        allPaths = []
        for info in self.infos:
            try:
                paths = self._loadPaths(info)
            except:
                continue

            print info.uniqueName

            frames = _makeCells([v.array for v in paths.velocities], mlab)
            allPaths.append(mlab.struct("info", _infoToMatlab(paths.info, mlab),
                                        "frames", _wrapCellValue(frames)))
        
        mlab._set("paths", _toStructArray(allPaths))
        mlab.save(fileName, "paths")
    
    def _performAnalysis(self, analysisName, 
                         recalculate=False, 
                         skipMissing=True,
                         radiuses=np.arange(0.5, 15, 0.5)):
        filePath = self._makeFilePath("analysis", analysisName + ".cache")
        
        pool = Pool(self.poolSize)
        manager = Manager()
        results = manager.dict()
        
        if not recalculate and path.exists(filePath):
            with open(filePath, 'rb') as f:
                existingResults = pickle.load(f)
                results.update(existingResults)
                
        def tasks():
            for info in self.infos:
                if info.uniqueName in results: continue
                print "Analyzing " + info.uniqueName
                
                pathsFile = self._getPathsFile(info)
                if not path.exists(pathsFile):
                    if not skipMissing:
                        raise ValueError("Cannot find tracks for %s!"
                                         % info.uniqueName)
                    
                    print "No tracks found for %s, skipping." % info.uniqueName
                    continue
                
                yield results, info, pathsFile, analysisName, radiuses
         
        pool.map(_doAnalysis, tasks())
        
        results = dict(results)
        with open(filePath, 'wb') as f:
            pickle.dump(results, f)
        
        return results
    
    def _organizeByGroup(self, results):
        grouped = dict()
        for info in self.infos:
            name = info.uniqueName
            if name in results:
                grouped.setdefault(info.group, []).append(results[name])
        return grouped
    
    def _makeFilePath(self, subfolder, name):
        folder = path.join(self.outputFolder, subfolder)
        if not path.exists(folder): os.makedirs(folder)
        
        return path.join(folder, name)
        
    def _getPathsFile(self, info):
        return self._makeFilePath("paths", info.uniqueName + ".paths")
            
    def _loadPaths(self, info):
        with open(self._getPathsFile(info), "rb") as f:
            return pickle.load(f)

def _trackSeq(arguments):
    info, filePath, recalculate = arguments
    
    if not recalculate and path.exists(filePath): return
    
    paths = processing.trackBacteria(images.ImageSeq(info))
    with open(filePath, "wb") as f:
        pickle.dump(paths, f)
        
def _doAnalysis(arguments):
    results, info, pathsFile, func, radiuses = arguments
    print pathsFile, func
    with open(pathsFile, 'rb') as f: paths = pickle.load(f)
    result = analysis.divideRadiuses(paths, ANALYSIS_FUNCTIONS[func], radiuses)
    results[info.uniqueName] = result

def _toStructArray(structs):
    array = None
    for i, struct in enumerate(structs):
        if array is None:
            array = struct
        else:
            array[i] = struct
    return array

def _infoToMatlab(info, connection):
    pairs = []
    for name, value in info.attributes.iteritems():
        if isinstance(value, unicode):
            value = str(value)
        pairs.extend([name, value])
    
    return connection.struct(*pairs)

def _makeCells(items, connection, blockSize=64):
    connection._set("makeCells", connection._do("@(varargin) varargin", nout=1))
    connection._set("joinCells", connection._do("@(c1, c2) [c1, c2]", nout=1))

    cells = connection.makeCells(*items[0:blockSize], nout=1)
    for offset in range(blockSize, len(items), blockSize):
        nextBlock = connection.makeCells(*items[offset:offset+blockSize], nout=1)
        cells = connection.joinCells(cells, nextBlock, nout=1)

    return cells

def _wrapCellValue(value, connection):
    connection._set("wrapValue", connection._do("@(v) {v}", nout=1))
    return connection.wrapValue(value, nout=1)

if __name__ == '__main__':
    batch = fromDensities()
    batch.track(override=True)
    for func in batch.funcs:
        batch.overlay(func)
        visual.saveCurrentPlot(func + ".png")
    



