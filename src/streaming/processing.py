'''
Created on Mar 18, 2012

@author: ckernan
'''

import cv
import numpy as np
import math
from utils import MatrixWrapper
from streaming.images import *
from SimpleCV import ImageClass
from SimpleCV.Features import BlobMaker

def computeForegroundMasks(images, lowThresh=-1.8, highThresh=1.3):
    """ Works with grayscale only.
    """
    size = cv.GetSize(images[0])
    def createImage(): 
        return cv.CreateImage(size, cv.IPL_DEPTH_32F, 1)
    
    scratch = createImage()
    
    average = createImage()
    diff = createImage()
    averageDiff = createImage()
    previous = createImage()
    
    first = True
    for image in images:
        cv.ConvertScale(image, scratch) # to float
        
        if first:
            first = False
        else:
            cv.Acc(scratch, average)
            cv.AbsDiff(scratch, previous, diff)
            cv.Acc(diff, averageDiff)

        cv.Copy(scratch, previous)
    
    scale = 1.0/(len(images) - 1)
    cv.ConvertScale(average, average, scale)
    cv.ConvertScale(averageDiff, averageDiff, scale)
    
    # make sure diff isn't 0
    cv.AddS(averageDiff, cv.Scalar(1.0), averageDiff)
    
    # make thresholds
    def makeThreshold(scale):
        threshold = createImage()
        cv.ScaleAdd(averageDiff, scale, average, threshold)
        return threshold
    
    highThresh = makeThreshold(highThresh)
    lowThresh = makeThreshold(lowThresh)
    
    masks = []
    for image in images:
        mask = cv.CreateImage(size, cv.IPL_DEPTH_8U, 1)
        cv.ConvertScale(image, scratch) # to float
        cv.InRange(scratch, lowThresh, highThresh, mask)
        cv.SubRS(mask, 255, mask) # reverse
        masks.append(mask)
    return ImageSeq(images.info, masks)

def applyMasks(images, masks):
    copy = copyBlank(images[0])
    for image, mask in zip(images, masks):
        cv.Zero(copy)
        cv.Copy(image, copy, mask)
        cv.Copy(copy, image)

def eliminatePixelNoise(masks, numIters=1):
    element = cv.CreateStructuringElementEx(3, 3, 1, 1, cv.CV_SHAPE_CROSS)
    for mask in masks:
        cv.MorphologyEx(mask, mask, None, element, cv.CV_MOP_OPEN, numIters)
        cv.MorphologyEx(mask, mask, None, element, cv.CV_MOP_CLOSE, numIters)
        #cv.Smooth(mask, mask, cv.CV_MEDIAN, 3)

class EllipseMatrix(MatrixWrapper):
    
    def __init__(self):
        columns = ["time", "x", "y", "angle", "area"]
        MatrixWrapper.__init__(self, columns)
        
        self._prevTime = -1.0e99
        self.frameStarts = []
        
    def append(self, time, x, y, angle, area):
        if time > self._prevTime:
            self._prevTime = time
            self.frameStarts.append(self._count)
        
        MatrixWrapper.append(self, time, x, y, angle, area)

def findBlobs(images, masks, minArea=1, maxArea=-1):
    """ images and masks should contain SimpleCV images!
    """
    maker = BlobMaker()
    # TODO: Automatically update the recursion limit!
    return [maker.extractFromBinary(m, i, minArea, maxArea) 
            for m, i in zip(masks, images)]
    
def debugBlobs(images, blobs, videoFile=None):
    for blob in blobs: blob.draw()
    if videoFile:
        images.outputVideo(videoFile)
    else:
        images.displayMovie()

def findEllipsesWithSimpleCV(blobs):
    ellipses = EllipseMatrix()
    
    for i, blobSet in enumerate(blobs):
        time = i #i*masks.info.dt
        for blob in blobSet:
            try:
                x, y = blob.centroid()
            except ZeroDivisionError:
                x, y = blob.center()
            ellipses.append(time, x, y, TO_RADIANS*blob.angle(), blob.area())
    
    ellipses.compact()
    return ellipses

def findEllipses(masks, minArea=1.0, maxArea=200.0):
    # time, x, y, angle, area
    ellipses = EllipseMatrix()
    memory = cv.CreateMemStorage()
    
    for i, mask in enumerate(masks):
        time = i #i*masks.info.dt
        
        contours = cv.FindContours(mask, memory)
        for contour in _contourIterator(contours):
            if not contour or len(contour) == 0: continue
            
            area = cv.ContourArea(contour)
            if area < minArea or area > maxArea: continue
            
            minRect = cv.MinAreaRect2(contour)
            angle = TO_RADIANS*minRect[2]
            
            try:
                moments = cv.Moments(contour)
                x = moments.m10/moments.m00
                y = moments.m01/moments.m00
            except:
                boundingBox = cv.BoundingRect(contour)
                x = boundingBox[0] + boundingBox[2]/2
                y = boundingBox[1] + boundingBox[3]/2
            
            ellipses.append(time, x, y, angle, area)
    
    ellipses.compact()
    return ellipses

#def findEllipses(masks, minArea=1.0, maxArea=200.0, drawContours=False):
#    # time, x, y, angle, area
#    ellipses = EllipseMatrix()
#    
#    for i, mask in enumerate(masks):
#        time = i #i*masks.info.dt
#        
#        contours = cv.FindContours(mask, cv.CreateMemStorage())
#        for contour in _contourIterator(contours):
#            if len(contour) < 6: continue
#            
#            points = cv.CreateMat(1, len(contour), cv.CV_32FC2)
#            for row, point in enumerate(contour): points[0, row] = point
#            
#            area = cv.ContourArea(points)
#            if area < minArea or area > maxArea: continue
#            
#            (x, y), _, angle = cv.FitEllipse2(points)
#            ellipses.append(time, x, y, angle, area)
#    
#    ellipses.compact()
#    return ellipses

def _contourIterator(contour):
    while contour:
        yield contour
        contour = contour.h_next()
        
def drawContours(images, masks):
    for image, mask in zip(images, masks):
        contours = cv.FindContours(mask, cv.CreateMemStorage())
        lineColor = cv.Scalar(255, 255, 255, 255)
        cv.DrawContours(image, contours, lineColor, lineColor, 100)

class TracksMatrix(MatrixWrapper):
    
    def __init__(self, tracks):
        columns = ["x", "y", "angle", "area", "time", "id"]
        MatrixWrapper.__init__(self, columns, tracks)

class VelocityMatrix(MatrixWrapper):
    
    def __init__(self, time):
        self.time = time
        columns = ["x", "y", "vx", "vy", "angle"]
        MatrixWrapper.__init__(self, columns)

class PathData(object):
    
    def __init__(self, info, velocities):
        self.info = info
        self.velocities = velocities
    
TO_RADIANS = math.pi/180 

def trackBacteria(images, debug=False):
    from mlabwrap import mlab
    
    masks = computeForegroundMasks(images)
    
    #simpleImages = SimpleImageSeq(images)
    #blobs = findBlobs(simpleImages, SimpleImageSeq(masks))
    #if debug:
    #    debugBlobs(simpleImages, blobs)
    #ellipses = findEllipses(blobs)
    
    ellipses = findEllipses(masks)
    
    params = mlab.struct("mem", 1, "dim", 2, "good", 2, "quiet", 0)
    matrix = ellipses.selectColumns("x", "y", "angle", "area", "time")
    maxDist = 5
    while maxDist > 0.5:
        try:
            tracks = TracksMatrix(mlab.track(matrix, maxDist, params))
            break
        except:
            maxDist /= 2.0
    
    timeUnits = images.info.dt
    lengthUnits = images.info.pixel
    velocities = [VelocityMatrix(t*timeUnits) for t in range(1, len(images))]
    
    oldId = tracks[0, "id"]
    oldPosition = tracks[0, ("x", "y")]*lengthUnits
    oldTime = 0
    
    for row in range(1, tracks.shape[0]):
        position = tracks[row, ("x", "y")]*lengthUnits
        time = tracks[row, "time"]
        angle = tracks[row, "angle"]
        id_ = tracks[row, "id"]
        
        if oldId != id_:
            oldId = id_
        else:
            dt = (time - oldTime)*timeUnits
            velocity = (position - oldPosition)/dt
            combined = position.tolist() + velocity.tolist()
            velocities[int(time) - 1].append(*(combined + [angle]))
        
        oldPosition = position  
        oldTime = time
    
    for matrix in velocities: matrix.compact()
    return PathData(images.info, velocities)

def debugForegroundMaskSettings(images, lowThresh=-1.8, highThresh=1.4):
    simpleImages = SimpleImageSeq(images)
    masks = computeForegroundMasks(images, lowThresh, highThresh)
    blobs = findBlobs(simpleImages, SimpleImageSeq(masks))
    debugBlobs(simpleImages, blobs)

if __name__ == '__main__':
    from streaming import config
    infos = loadSeqInfos(config.densities["max"])
    seq = ImageSeq(infos[0])
    debugForegroundMaskSettings(seq)

    
    

        