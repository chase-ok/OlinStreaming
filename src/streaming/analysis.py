'''
Created on Mar 24, 2012

@author: ckernan
'''

import numpy as np
import math
from processing import PathData, VelocityMatrix
from utils import MatrixWrapper

EPSILON = 0.00001

class GriddedVelocityMatrix(VelocityMatrix):
    
    def __init__(self, binSize, cellMap, velocities):
        self.binSize = binSize
        self.cellMap = cellMap
        
        def index(i): return max([int(i), 0])
        
        cells = np.empty((velocities.shape[0], 1))
        for row in range(velocities.shape[0]):
            i, j = velocities[row, ("x", "y")] // binSize
            cells[row] = int(cellMap[index(i), index(j)])
        
        self.time = velocities.time # skipping VelocityMatrix init
        columns = velocities.columns + ["cell"]
        newArray = np.append(velocities.array, cells, 1)
        MatrixWrapper.__init__(self, columns, newArray)

class GriddedPathData(PathData):
    
    def __init__(self, paths, numBins=(10, 10)):
        self.binSize = paths.info.imageSize/numBins
        
        self.cellMap = {}
        cellCenters = []
        for i in range(numBins[0]):
            for j in range(numBins[1]):
                self.cellMap[i, j] = i*numBins[1] + j
                cellCenters.append([i, j]*self.binSize + self.binSize/2)
        self.cellCenters = np.array(cellCenters)
        
        velocities = [GriddedVelocityMatrix(self.binSize, self.cellMap, v)
                      for v in paths.velocities]
        PathData.__init__(self, paths.info, velocities)
        
    def makeCellData(self, func):
        return [func() for _ in range(len(self.cellMap))]

def _ensureGridded(paths):
    if isinstance(paths, GriddedPathData):
        return paths
    else:
        return GriddedPathData(paths)
    

def velocityField(paths):
    paths = _ensureGridded(paths)
    
    def calculateFrame(matrix):
        velocities = paths.makeCellData(lambda: np.zeros(2))
        counts = paths.makeCellData(lambda: 0.00001)
        
        for row in range(matrix.shape[0]):
            cell = int(matrix[row, "cell"])
            velocities[cell] += matrix[row, ("vx", "vy")]
            counts[cell] += 1
        
        meanVelocities = np.array([v/c for v, c in zip(velocities, counts)])
        return paths.cellCenters, meanVelocities
            
    return _calculateByFrame(paths, calculateFrame)

def particleDistance(paths):
    def calculateRadius(m, ref, inRange):
        return inRange.sum()
    
    return _byRadius(paths, calculateRadius)

def velocityCorrelation(paths, **radiusArgs):
    def calculateRadius(m, ref, inRange):
        velocities = m[inRange, ("vx", "vy")]
        numRows = velocities.shape[0]
        if numRows == 0: return None
        
        refVelocity = _toUnit(m[ref, ("vx", "vy")])
        corr = 0.0
        
        for row in range(numRows):
            corr += np.dot(refVelocity, _toUnit(velocities[row, :]))
        
        return corr/numRows
    
    return _byRadius(paths, calculateRadius, **radiusArgs)

def directorCorrelation(paths, **radiusArgs):
    def calculateRadius(m, ref, inRange):
        angles = m[inRange, "angle"]
        if angles.size == 0: return None
        
        refDirector = _toDirector(m[ref, "angle"])
        corr = 0.0
        
        for angle in angles:
            # abs because we can't tell if pointing forwards or backwards
            corr += abs(np.dot(refDirector, _toDirector(angle)))
        
        return corr/angles.size
    
    return _byRadius(paths, calculateRadius, **radiusArgs)

def directorVelocityCorrelation(paths, **radiusArgs):
    def calculateRadius(m, ref, inRange):
        velocities = m[inRange, ("vx", "vy")]
        numRows = velocities.shape[0]
        if numRows == 0: return None
        
        refDirector = _toDirector(m[ref, "angle"])
        corr = 0.0
        
        for row in range(numRows):
            corr += abs(np.dot(refDirector, _toUnit(velocities[row, :])))
        
        return corr/numRows
    
    return _byRadius(paths, calculateRadius, **radiusArgs)

def _toUnit(vector):
    mag = np.linalg.norm(vector)
    return vector/mag if mag > EPSILON else 0*vector

def _toDirector(angle):
    return np.array([np.cos(angle), np.sin(angle)])
    
def _makeBins(paths, cellData, shape=(10, 10)):
    shape = np.array(shape)
    binSize = paths.info.imageSize/shape
    
    positions = {}
    for row in range(shape[0]):
        for col in range(shape[1]):
            positions[row, col] = [row, col]*binSize + binSize/2
    
    return binSize, positions
    
def _calculateByFrame(paths, func):
    if isinstance(paths, PathData):
        return map(func, paths.velocities)
    else:
        return func(paths)

def _byRadius(paths, func, 
              radiuses=np.arange(1.0, 15, 1),
              divideArea=True):
    radiusesSq = radiuses**2
    
    def getAreaConstant(x, y, rSq1, rSq2):
        if divideArea:
            return math.pi*(rSq2 - rSq1)
        else:
            return 1.0
    
    def calculateFrame(matrix):
        data = np.zeros(len(radiuses))
        
        for row in range(matrix.shape[0]):
            x, y = matrix[row, ("x", "y")]
            distancesSq = (matrix[:, "x"] - x)**2 + (matrix[:, "y"] - y)**2
            
            for i in range(len(radiusesSq)):
                rSq1 = EPSILON if i == 0 else radiusesSq[i - 1]
                rSq2 = radiusesSq[i]
                
                particlesInRange = (distancesSq >= rSq1) & (distancesSq <= rSq2)
                result = func(matrix, row, particlesInRange)
                if result:
                    data[i] += result/getAreaConstant(x, y, rSq1, rSq2)
        
        data /= matrix.shape[0]
            
        return radiuses, data
    
    return _calculateByFrame(paths, calculateFrame)
