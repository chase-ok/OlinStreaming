
import numpy as np
cimport numpy as np

import math
import copy
from operator import itemgetter
from processing import PathData, VelocityMatrix
from utils import MatrixWrapper

cdef extern from "math.h":
    float fabs(float x)
    float cosf(float theta)
    float sinf(float theta)
    float acosf(float theta)

Num = np.float
ctypedef np.float_t Num_t

cdef Num_t EPSILON = 0.00001

cdef int maxInt(int x, int y):
    return x if x > y else y

cdef int index(Num_t i): 
    return maxInt(int(i), 0)

class GriddedVelocityMatrix(VelocityMatrix):
    
    def __init__(self, binSize, cellMap, velocities):
        self.binSize = binSize
        self.cellMap = cellMap

        cdef np.ndarray[np.int_t, ndim=2] cellMap_ = cellMap
        cdef np.ndarray[Num_t, ndim=1] cells = np.empty(velocities.shape[0], dtype=Num)
        cdef int row
        cdef Num_t i, j
        for row in range(velocities.shape[0]):
            i, j = velocities[row, ("x", "y")] // binSize
            cells[row] = cellMap_[index(i), index(j)]
        
        self.time = velocities.time # skipping VelocityMatrix init
        columns = velocities.columns + ["cell"]
        newArray = np.append(velocities.array, cells, 1)
        MatrixWrapper.__init__(self, columns, newArray)

class GriddedPathData(PathData):
    
    def __init__(self, paths, numBins=(10, 10)):
        self.numBins = numBins
        self.binSize = paths.info.imageSize/numBins

        cdef int numBinsX = numBins[0]
        cdef int numBinsY = numBins[1]
        cdef Num_t binSizeX = self.binSize[0]
        cdef Num_t binSizeY = self.binSize[1]
        cdef Num_t binSizeXHalf = binSizeX/2.0
        cdef Num_t binSizeYHalf = binSizeY/2.0
        
        cdef np.ndarray[np.int_t, ndim=2] cellMap = np.empty((numBinsX, numBinsY), dtype=np.int)
        cdef np.ndarray[Num_t, ndim=2] cellCenters = np.empty((numBinsX*numBinsY, 2), dtype=Num)
        
        cdef int i, j
        cdef int c = 0
        for i in range(numBinsX):
            for j in range(numBinsY):
                cellMap[i, j] = i*numBinsY + j
                cellCenters[c, 0] = i*binSizeX + binSizeXHalf
                cellCenters[c, 1] = j*binSizeY + binSizeYHalf
                c += 1

        self.cellMap = cellMap
        self.cellCenters = cellCenters
        self.cellCentersAsGrids = self.reshapeColumnAsGrid(self.cellCenters[:, 0]), \
                                  self.reshapeColumnAsGrid(self.cellCenters[:, 1])
        
        velocities = [GriddedVelocityMatrix(self.binSize, self.cellMap, v)
                      for v in paths.velocities]
        PathData.__init__(self, paths.info, velocities)
        
    def makeCellData(self, func):
        return [func() for _ in range(self.numCells)]
    
    def reshapeColumnAsGrid(self, col):
        return np.array(col).reshape(self.numBins)

    @property
    def numCells(self):
        return len(self.cellMap)
    

def _ensureGridded(paths, numBins=(12, 12), **others):
    if isinstance(paths, GriddedPathData):
        return paths
    else:
        return GriddedPathData(paths, numBins=numBins)
    

def velocityField(paths, **binOptions):
    paths = _ensureGridded(paths, **binOptions)    
    return _calculateByFrame(paths, lambda m: velocityField_calculateFrame(paths, m))

cdef object velocityField_calculateFrame(object paths, object matrix):
    cdef int n = paths.numCells
    cdef np.ndarray[Num_t, ndim=2] velocities = np.zeros((n, 2), dtype=Num)
    cdef np.ndarray[Num_t, ndim=1] counts = 0.000001*np.ones(n, dtype=Num)

    cdef object columns = ("vx", "vy")
    cdef int row, cell
    for row in range(matrix.shape[0]):
        cell = <int>(matrix[row, "cell"])
        velocities[cell, :] += matrix[row, columns]
        counts[cell] += 1.0
    
    cdef np.ndarray[Num_t, ndim=2] meanVelocities = np.empty(velocities.size, dtype=Num)
    for row in range(n):
        meanVelocities[row, 0] = velocities[row, 0]/counts[row]
        meanVelocities[row, 1] = velocities[row, 1]/counts[row]

    return paths.cellCenters, meanVelocities

def localVelocityCorrelation(paths, **binOptions):
    paths = _ensureGridded(paths, **binOptions)
    return _calculateByFrame(paths, lambda m: localVelocityCorrelation_calculateFrame(paths, m))
    
cdef object localVelocityCorrelation_calculateFrame(object paths, object matrix):
    cdef int n = paths.numCells
    cdef list velocities[n]
    cdef int i
    for i in range(n): velocities[i] = []
    
    cdef object columns = ("vx", "vy")
    cdef int cell
    for i in range(matrix.shape[0]):
        cell = <int>(matrix[row, "cell"])
        velocities[cell].append(matrix[row, columns])
    
    cdef np.ndarray[Num_t, ndim=1] corr = np.empty(n, dtype=Num)
    for i in range(n):
        corr[i] = _vectorCorrelation(velocities[i])

    return paths.cellCentersAsGrids, paths.reshapeColumnAsGrid(corr)

def particleDistance(paths, **radiusArgs):
    def calculateRadius(m, ref, inRange):
        return inRange.sum()
    
    return Correlation("Particle Distance", paths.info,
                       _byRadius(paths, calculateRadius, **radiusArgs))

cpdef Num_t _vectorCorrelation(list vectors):
    cdef int n = len(vectors)
    if n == 0: return 0.0
    if n == 1: return 1.0
    
    cdef np.ndarray[Num_t, ndim=1] ref = _toUnit(vectors[0])
    cdef Num_t corr = 0.0
    cdef int i
    for i in range(1, n):
        corr += np.dot(ref, _toUnit(vector))
    
    return corr/n

#we're not gonna use the _vectorCorrelation function here because the matrix
#operations are much faster
def velocityCorrelation(paths, **radiusArgs):
    return Correlation("Velocity", paths.info,
                       _byRadius(paths, velocityCorrelation_calculateRadius, **radiusArgs))

cdef object velocityCorrelation_calculateRadius(object m, int ref, np.ndarray[np.uint8_t, ndim=1] inRange):
    cdef np.ndarray[Num_t, ndim=2] velocities = m[inRange, ("vx", "vy")]
    cdef int numRows = velocities.shape[0]
    if numRows == 0: return None
    
    cdef np.ndarray[Num_t, ndim=1] refVelocity = _toUnit(m[ref, ("vx", "vy")])
    cdef Num_t corr = 0.0
    cdef int row
    for row in range(numRows):
        corr += np.dot(refVelocity, _toUnit(velocities[row, :]))
    
    return corr/numRows

def directorCorrelation(paths, **radiusArgs):
    return Correlation("Director", paths.info,
                       _byRadius(paths, calculateRadius, **radiusArgs))

cdef object directorCorrelation_calculateRadius(object m, int ref, np.ndarray[np.uint8_t, ndim=1] inRange):
    cdef np.ndarray[Num_t, ndim=1] angles = m[inRange, "angle"]
    if angles.size == 0: return None
    
    cdef Num_t refDirector = _toDirector(m[ref, "angle"])
    cdef Num_t corr = 0.0
    cdef int i
    for i in range(angles.size):
        # abs because we can't tell if pointing forwards or backwards
        corr += fabs(np.dot(refDirector, _toDirector(angle)))
    
    return corr/angles.size

def directorVelocityCorrelation(paths, **radiusArgs):
    return Correlation("Director-Velocity", paths.info,
                       _byRadius(paths, directorVelocityCorrelation_calculateRadius, **radiusArgs))

cdef object directorVelocityCorrelation_calculateRadius(object m, int ref, np.ndarray[np.uint8_t, ndim=1] inRange):
    cdef np.ndarray[Num_t, ndim=2] velocities = m[inRange, ("vx", "vy")]
    cdef int numRows = velocities.shape[0]
    if numRows == 0: return None
    
    cdef Num_t refDirector = _toDirector(m[ref, "angle"])
    cdef Num_t corr = 0.0
    cdef int row
    for row in range(numRows):
        corr += fabs(np.dot(refDirector, _toUnit(velocities[row, :])))
    
    return corr/numRows

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

cpdef np.ndarray[Num_t, ndim=1] _toUnit(np.ndarray[Num_t, ndim=1] vector):
    cdef Num_t mag = sqrt(vector[0]*vector[0] + vector[1]*vector[1])
    cdef np.ndarray[Num_t, ndim=1] unit = np.zeros(2, dtype=Num)
    if mag > EPSILON:
        unit[0] = vector[0]/mag
        unit[1] = vector[1]/mag
    return unit

cpdef void _modifyToUnit(np.ndarray[Num_t, ndim=1] vector):
    cdef Num_t mag = sqrt(vector[0]*vector[0] + vector[1]*vector[1])
    if mag > EPSILON:
        vector[0] = vector[0]/mag
        vector[1] = vector[1]/mag
    else:
        vector[0] = 0
        vector[1] = 0

cpdef np.ndarray[Num_t, ndim=1] _toDirector(Num_t angle):
    cdef np.ndarray[Num_t, ndim=1] v = np.empty(2, dtype=Num)
    v[0] = cos(angle)
    v[1] = sin(angle)
    return v
    
def _calculateByFrame(paths, func):
    return map(func, paths.velocities)

def _byRadius(object paths, object func, 
              np.ndarray[Num_t, ndim=1] radiuses=np.arange(1.0, 15, 1, dtype=Num),
              bool divideArea=True):
    cdef object getAreaConstant
    cdef object circleArea
    if divideArea:
        circleArea = _createCircleAreaFunc(*paths.info.imageSize)
        getAreaConstant = lambda x, y, r1, r2: circleArea(x, y, r2) - circleArea(x, y, r1)
    else:
        getAreaConstant = _byRadius_constantArea
    
    return _calculateByFrame(paths, lambda m: _byRadius_calculateFrame(m, radiuses, getAreaConstant))

cdef Num_t _byRadius_constantArea(Num_t x, Num_t y, Num_t r1, Num_t r2):
    return 1.0

cdef object _byRadius_calculateFrame(object matrix, np.ndarray[Num_t, ndim=1] radiuses, object getAreaConstant):
    cdef np.ndarray[Num_t, ndim=1] data = np.zeros(radiuses.size)
    
    cdef int row, i
    cdef Num_t x, y, r1, r2
    cdef np.ndarray[Num_t, ndim=1] distancesSq
    cdef np.ndarray[np.uint8_t, ndim=1] particlesInRange

    for row in range(matrix.shape[0]):
        x = matrix[row, "x"]
        y = matrix[row, "y"]

        distancesSq = (matrix[:, "x"] - x)**2 + (matrix[:, "y"] - y)**2
        
        for i in range(radiuses.size):
            r1 = EPSILON if i == 0 else radiuses[i - 1]
            r2 = radiuses[i]
            
            particlesInRange = (distancesSq >= r1**2) & (distancesSq <= r2**2)
            result = func(matrix, row, particlesInRange)
            if result:
                data[i] += result/getAreaConstant(x, y, r1, r2)
    
    data /= matrix.shape[0]
        
    return radiuses, data
    
DEF Pi = 3.14159265

def _createCircleAreaFunc(width, height):
    cdef dict table = _createClippingFuncTable(width, height)
    cdef dict cache = {}
    return lambda x, y, r: _createCircleAreaFunc_calculateArea(table, cache, width, height, x, y, r)

cdef Num_t _createCircleAreaFunc_calculateArea(dict table, 
                                               dict cache, 
                                               Num_t width, 
                                               Num_t height, 
                                               Num_t x, 
                                               Num_t y, 
                                               Num_t r):
    cdef int roundedX = <int>x
    cdef int roundedY = <int>y
    cdef object key = roundedX, roundedY, r
    cdef int tableIndex
    cdef Num_t area

    if key in cache:
        return cache[key]
    else:
        tableIndex = 0
        if roundedX - r >= 0:      tableIndex += 0b1000
        if roundedX + r <= width:  tableIndex += 0b0100
        if roundedY - r >= 0:      tableIndex += 0b0010
        if roundedY + r <= height: tableIndex += 0b0001
        area = table[tableIndex](roundedX, roundedY, r)
        cache[key] = area
        return area

def _createClippingFuncTable(width, height):
    def singleSideClipping(transform):
        return lambda x, y, r: _singleSideClipping(transform, x, y, r)
    
    def cornerClipping(transform):
        return lambda x, y, r: _cornerClipping(transform, x, y, r)
    
    # (left, right, top, bottom)
    b = lambda s: tuple(c == "1" for c in s)
    cdef dict table
    table[0b1111] = _noClipping
    table[0b0111] = singleSideClipping(lambda x, y: x)
    table[0b1011] = singleSideClipping(lambda x, y: width - x)
    table[0b1101] = singleSideClipping(lambda x, y: y)
    table[0b1110] = singleSideClipping(lambda x, y: height - y)
    table[0b0101] = cornerClipping(lambda x, y: (x, y))
    table[0b0110] = cornerClipping(lambda x, y: (x, height - y))
    table[0b1001] = cornerClipping(lambda x, y: (width - x, y))
    table[0b1010] = cornerClipping(lambda x, y: (width - x, height - y))

cdef Num_t _noClipping(Num_t x, Num_t y, Num_t r): 
    return Pi*r**2

cdef Num_t _singleSideClipping(object transform, Num_t x, Num_t y, Num_t r):
    x = fabs(transform(x, y))
    return x*sqrt(r**2 - x**2) + r**2*(Pi/2 + asin(x/r))

cdef Num_t _cornerClipping(object transform, Num_t x, Num_t y, Num_t r):
    x, y = transform(x, y)
    cdef Num_t x2 = x**2
    cdef Num_t y2 = y**2
    cdef Num_t r2 = r**2

    cdef Num_t xIntercept, yIntercept, angle, triangle, chord
    
    if x2 + y2 < r2:
        xIntercept = x + sqrt(r2 - y2)
        yIntercept = y + sqrt(r2 - x2)
        angle = 2*asin(sqrt(xIntercept*yIntercept)/(2*r))
        triangle = 0.5*xIntercept*yIntercept
        chord = 0.5*r2*(angle - sin(angle))
        return triangle + chord
    else:
        return x*sqrt(r2 - x2) + y*sqrt(r2 - y2) + \
               r2*(asin(x/r) + asin(y/r))
        