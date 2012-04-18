'''
Created on Apr 10, 2012

@author: ckernan
'''
import unittest
import analysis
import processing
import images
import random
import numpy as np
import math

class AnalysisTest(unittest.TestCase):

    def setUp(self):
        self.perfect = _constructPerfectCorrelation(100.0)
        self.radiuses = np.arange(15, 80, 10)

    def testVelocityCorrelation(self):
        rc = analysis.velocityCorrelation(self.perfect, 
                                          radiuses=self.radiuses,
                                          divideArea=False)
        for _, correlation in rc:
            for c in correlation:
                self.assertAlmostEqual(c, 1.0)

    def testDirectorCorrelation(self):
        rc = analysis.directorCorrelation(self.perfect, 
                                          radiuses=self.radiuses,
                                          divideArea=False)
        for _, correlation in rc:
            for c in correlation:
                self.assertAlmostEqual(c, 1.0)
                
    def testDirectorVelocityCorrelation(self):
        rc = analysis.directorVelocityCorrelation(self.perfect, 
                                          radiuses=self.radiuses,
                                          divideArea=False)
        for _, correlation in rc:
            for c in correlation:
                self.assertAlmostEqual(c, 1.0)

def _constructPerfectCorrelation(frameSize):
    numFrames = 10
    numParticlesPerSide = 10
    images = _makeMockImages(numFrames, (frameSize, frameSize))
    
    spacing = frameSize/numParticlesPerSide
    def toPos(index): return index*spacing + spacing/2 
    
    velocities = []
    for t in range(numFrames):
        v = processing.VelocityMatrix(t)
        angle = random.random()*2*math.pi
        x = math.cos(angle)
        y = math.sin(angle)
        
        for i in range(numParticlesPerSide):
            for j in range(numParticlesPerSide):
                r = random.random()*10
                v.append(toPos(i), toPos(j), x*r, y*r, angle)
        
        v.compact()
        velocities.append(v)
    
    return processing.PathData(images, velocities)
        
def _makeMockImages(n, size=(100, 100)):
    info = images.ImageSeqInfo(name="MockImages", seriesNum=0, 
                               length=n, channel=0, imageSize=np.array(size))
    return images.ImageSeq(info, [None for _ in range(n)])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()