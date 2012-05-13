'''
Created on Mar 24, 2012

@author: ckernan
'''

import numpy as np
import math
import copy
from operator import itemgetter
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
        self.numBins = numBins
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
    
    def reshapeColumnAsGrid(self, col):
        return np.array(col).reshape(self.numBins)
    
    @property
    def cellCentersAsGrids(self):
        return self.reshapeColumnAsGrid(self.cellCenters[:, 0]), \
               self.reshapeColumnAsGrid(self.cellCenters[:, 1])
    

def _ensureGridded(paths, numBins=(12, 12), **others):
    if isinstance(paths, GriddedPathData):
        return paths
    else:
        return GriddedPathData(paths, numBins=numBins)
    

def velocityField(paths, **binOptions):
    paths = _ensureGridded(paths, **binOptions)
    
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

def localVelocityCorrelation(paths, **binOptions):
    paths = _ensureGridded(paths, **binOptions)
    
    def calculateFrame(matrix):
        velocities = paths.makeCellData(lambda: [])
        
        for row in range(matrix.shape[0]):
            cell = int(matrix[row, "cell"])
            velocities[cell].append(matrix[row, ("vx", "vy")])
        
        corr = map(_vectorCorrelation, velocities)
        return paths.cellCentersAsGrids, paths.reshapeColumnAsGrid(corr)
    
    return _calculateByFrame(paths, calculateFrame)

def particleDistance(paths, **radiusArgs):
    def calculateRadius(m, ref, inRange):
        return inRange.sum()
    
    return Correlation("Particle Distance", paths.info,
                       _byRadius(paths, calculateRadius, **radiusArgs))

def _vectorCorrelation(vectors):
    if len(vectors) == 0: return 0.0
    if len(vectors) == 1: return 1.0
    
    ref = _toUnit(vectors[0])
    corr = 0.0
    for vector in vectors[1:]:
        corr += np.dot(ref, _toUnit(vector))
    
    return corr/len(vectors)

#we're not gonna use the _vectorCorrelation function here because the matrix
#operations are much faster
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
    
    return Correlation("Velocity", paths.info,
                       _byRadius(paths, calculateRadius, **radiusArgs))

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
    
    return Correlation("Director", paths.info,
                       _byRadius(paths, calculateRadius, **radiusArgs))

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
    
    return Correlation("Director-Velocity", paths.info,
                       _byRadius(paths, calculateRadius, **radiusArgs))

class Correlation(object):
    
    def __init__(self, name, info, frames):
        self.info = info
        self._frames = frames
        self.name = name
        
    @property
    def x(self):
        return self._frames[0][0]
    
    @property
    def meanY(self):
        sums = np.zeros_like(self._frames[0][1])
        for _, ys in self._frames: sums += ys
        return sums/len(self._frames)
    
    @property
    def points(self):
        points = []
        for xs, ys in self._frames:
            points.extend(zip(xs, ys))
        return points
    
    @property
    def numFrames(self):
        return len(self._frames)
    
    def copy(self):
        return Correlation(self.name, self.info, copy.deepcopy(self._frames))
    
    def merge(self, other):
        assert other.numFrames == self.numFrames
        
        newFrames = []
        for selfFrame, otherFrame in zip(self._frames, other._frames):
            selfXY = list(enumerate(zip(*selfFrame)))
            otherXY = list(enumerate(zip(*otherFrame)))
            
            allXY = sorted(selfXY + otherXY, key=itemgetter(1))
            newFrames.append(zip(*(point for _, point in allXY)))
        
        return Correlation(self.name, self.info, newFrames)

def mergeCorrelations(correlations):
    base = correlations[0].copy()
    for c in correlations[1:]: base.merge(c)
    return base

def meanOfCorrelations(correlations):
    sums = np.zeros_like(correlations[0].x)
    for c in correlations:
        sums += c.meanY
    return correlations[0].x, sums/len(correlations)

def divideRadiuses(paths, correlFunc, radiuses, n=2):
    radiusSets = [radiuses[i:len(radiuses)-i:n] for i in range(n)]
    correlations = [correlFunc(paths, radiuses=rs) for rs in radiusSets]
    return mergeCorrelations(correlations)
    
def meanDensity(paths):
    area = paths.info.imageSize.prod()
    numFrames = float(len(paths.velocities))
    return (sum(m.shape[0] for m in paths.velocities)/area)/numFrames
    
def meanSpeed(paths):
    speed = 0.0
    for vel in paths.velocities:
        # axis=1 -> sum along rows
        speeds = np.sqrt((vel[:, ("vx", "vy")]**2).sum(axis=1))
        speed += speeds.sum()/len(speeds)
    return speed

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
    return map(func, paths.velocities)

def _byRadius(paths, func, 
              radiuses=np.arange(1.0, 15, 1),
              divideArea=True):
    if divideArea:
        circleArea = _createCircleAreaFunc(*paths.info.imageSize)
        def getAreaConstant(x, y, r1, r2):
            return circleArea(x, y, r2) - circleArea(x, y, r1)
    else:
        def getAreaConstant(x, y, r1, r2):
            return 1
    
    def calculateFrame(matrix):
        data = np.zeros(len(radiuses))
        
        for row in range(matrix.shape[0]):
            x, y = matrix[row, ("x", "y")]
            distancesSq = (matrix[:, "x"] - x)**2 + (matrix[:, "y"] - y)**2
            
            for i in range(len(radiuses)):
                r1 = EPSILON if i == 0 else radiuses[i - 1]
                r2 = radiuses[i]
                
                particlesInRange = (distancesSq >= r1**2) & (distancesSq <= r2**2)
                result = func(matrix, row, particlesInRange)
                if result:
                    data[i] += result/getAreaConstant(x, y, r1, r2)
        
        data /= matrix.shape[0]
            
        return radiuses, data
    
    return _calculateByFrame(paths, calculateFrame)

def _createCircleAreaFunc(width, height):
    table = _createClippingFuncTable(width, height)
    cache = {}
    
    def calculateArea(x, y, r):
        roundedX, roundedY = int(x), int(y)
        key = roundedX, roundedY, r
        try:
            return cache[key]
        except KeyError:
            assert r < width/2 and r < height/2
            area = table[roundedX - r >= 0, 
                         roundedX + r <= width, 
                         roundedY - r >= 0, 
                         roundedY + r <= height](roundedX, roundedY, r)
            cache[key] = area
            return area
    
    return calculateArea


def _createClippingFuncTable(width, height):
    def noClipping(x, y, r): 
        return math.pi*r**2
    
    def singleSideClipping(transform):
        def func(x, y, r):
            x = abs(transform(x, y))
            return x*math.sqrt(r**2 - x**2) + r**2*(math.pi/2 + math.asin(x/r))
        return func
    
    def cornerClipping(transform):
        def func(x, y, r):
            x, y = transform(x, y)
            x2, y2, r2 = x**2, y**2, r**2
            
            if x2 + y2 < r2:
                xIntercept = x + math.sqrt(r2 - y2)
                yIntercept = y + math.sqrt(r2 - x2)
                angle = 2*math.asin(math.sqrt(xIntercept*yIntercept)/(2*r))
                triangle = 0.5*xIntercept*yIntercept
                chord = 0.5*r2*(angle - math.sin(angle))
                return triangle + chord
            else:
                return x*math.sqrt(r2 - x2) + y*math.sqrt(r2 - y2) + \
                       r2*(math.asin(x/r) + math.asin(y/r))
        
        return func
    
    # (left, right, top, bottom)
    b = lambda s: tuple(c == "1" for c in s)
    return {b("1111"): noClipping,
            b("0111"): singleSideClipping(lambda x, y: x),
            b("1011"): singleSideClipping(lambda x, y: width - x),
            b("1101"): singleSideClipping(lambda x, y: y),
            b("1110"): singleSideClipping(lambda x, y: height - y),
            b("0101"): cornerClipping(lambda x, y: (x, y)),
            b("0110"): cornerClipping(lambda x, y: (x, height - y)),
            b("1001"): cornerClipping(lambda x, y: (width - x, y)),
            b("1010"): cornerClipping(lambda x, y: (width - x, height - y))}


    
