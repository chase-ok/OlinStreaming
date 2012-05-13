'''
Created on Mar 12, 2012

@author: ckernan
'''

import cv
import time
import os
from xml.dom import minidom
import re
import numpy as np
from SimpleCV import ImageClass, Stream

def loadImage(path, grayscale=True):
    return cv.LoadImage(path, not grayscale)

def copyImage(image):
    copy = copyBlank(image)
    cv.Copy(image, copy)
    return copy

def copyBlank(image):
    return cv.CreateImage(cv.GetSize(image), image.depth, image.channels)

def displayImage(image, title="Image"):
    cv.NamedWindow(title)
    cv.ShowImage(title, image)
    cv.WaitKey()
    cv.DestroyWindow(title)
    
_imageSeriesName = re.compile('Series([0-9]+)')
_propertiesFile = re.compile(r"Series[0-9]+_Properties\.xml")
    
def loadSeqInfos(folder, **args):
    return [parseConfigFile(os.path.join(folder, f), **args)
            for f in os.listdir(folder)
            if _propertiesFile.match(f)]
    
class ImageSeq(object):
    
    def __init__(self, info, images=None):
        self.info = info
        self._images = images or self._loadImages()
        
    def __len__(self): 
        return len(self._images)
    
    def __getitem__(self, i):
        return self._images.__getitem__(i)
    
    def __getslice__(self, i, j):
        return self._images.__getslice__(i, j)
    
    def __iter__(self):
        return iter(self._images)
    
    def __str__(self):
        return str(self.info) + "; in memory"
    
    def __repr__(self):
        return "ImageSeq(" + repr(self.info) + ")"
    
    def _loadImages(self):
        return [loadImage(path) for path in self.info.imagePaths]
        
    def copy(self):
        def copyImage(image):
            size = cv.GetSize(image)
            newImage = cv.CreateImage(size, image.depth, image.channels)
            cv.Copy(image, newImage)
            return newImage
        
        return ImageSeq(self.info, [copyImage(im) for im in self._images])
    
    def displayMovie(self, title="Movie", delay=0.01):
        cv.NamedWindow(title)
        for image in self._images:
            cv.ShowImage(title, image)
            cv.WaitKey() # TODO: FIXME
            time.sleep(delay)
        cv.DestroyWindow(title)

        
class SimpleImageSeq(ImageSeq):
    
    def __init__(self, sourceSeq):
        self.info = sourceSeq.info
        self._images = map(ImageClass.Image, sourceSeq)
    
    def copy(self):
        assert False #IMPLEMENT ME!
        
    def displayMovie(self, title="Movie", delay=0.01):
        for image in self._images:
            image.show()
            
    def outputVideo(self, path, fps=24):
        vs = Stream.VideoStream(path, fps=fps, framefill=False)
        for image in self._images: image.save(vs)

class ImageSeqInfo(object):
    
    def __init__(self, **options):
        required = ['name', 'seriesNum', 'length', 'channel']
        for req in required: 
            assert req in options
        
        for option, value in options.iteritems(): 
            setattr(self, option, value)
    
    @property
    def imagePaths(self):
        def imagePath(imageNum):
            relative = '{0}_t{1}_ch{2:02d}.tif'\
                       .format(self.name, imageNum, self.channel)
            return os.path.join(self.folder, relative)
        
        numDigits = len(str(self.length - 1))
        imageNums = (str(i).zfill(numDigits) for i in range(self.length))
        return [imagePath(num) for num in imageNums]

    @property
    def attributes(self):
        return self.__dict__
    
    def __str__(self):
        return "{0} in {1}: {2} images (channel {3})"\
               .format(self.name, os.path.basename(self.folder), self.length, 
                       self.channel)
        
    def __repr__(self):
        return "ImageSeqInfo(" + repr(self.attributes) + ")"
    
    @property
    def uniqueName(self):
        if hasattr(self, "group"):
            return self.group + "-" + self.name
        else:
            return self.startTime + "-" + self.name
    
def parseConfigFile(configFile, **args):
    xml = minidom.parse(configFile)
    dims = xml.getElementsByTagName('DimensionDescription')
    
    name = xml.getElementsByTagName('Name')[0].firstChild.data
    seriesNum = int(_imageSeriesName.match(name).group(1))
    length = int(dims[2].attributes['NumberOfElements'].value)
    pixel = np.array([float(dims[i].attributes['Voxel'].value) for i in [0, 1]])
    shape = np.array([int(dims[i].attributes['NumberOfElements'].value) for i in [0, 1]])
    duration = float(dims[2].attributes['Length'].value)
    timeString = xml.getElementsByTagName('StartTime')[0].firstChild.data
    startTime = timeString.split(' ')[1]

    return ImageSeqInfo(folder=os.path.dirname(configFile), 
                        configFile=configFile,
                        name=name, 
                        seriesNum=seriesNum, 
                        length=length, 
                        pixel=pixel,
                        dt=duration/length,
                        startTime=startTime,
                        imageSize=shape*pixel,
                        channel=0,
                        **args)

if __name__ == "__main__":
    import config
    infos = loadSeqInfos(config.densities['max'])
    test = ImageClass.Image(infos[0].imagePaths[0])
    test.show()
    
    
    