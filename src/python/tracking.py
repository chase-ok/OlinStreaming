""" 
Used as an easier bridge for matlab to read the images produced by the confocal.

It takes 3 arguments, the base directory of the data files, the image series 
number and the channel number to use.

It will print the following to stdout (assuming correct inputs), one line at a 
time:
    Start time (formatted as HH:MM:SS.FFF)
    dt
    Pixel size in x, y
    Image files (in order)
    
Usage:
    python tracking.py BASE_DIR SERIES_NUM CHANNEL_NUM
"""


import sys
import os.path
import shutil
from xml.dom import minidom
from optparse import OptionParser
    
def readProps(baseDir, seriesNum, channelNum):
    """
    Returns a dictionary of the relevant properties for a series from the
    Properties.xml file.
    
    Included properties:
        - seriesNum
        - baseDir
        - channelNum
        - numFrames
        - pixel (tuple of x, y)
        - startTime
        - dt
    
    NOTE: all in micrometers
    """
    props = { 'baseDir': baseDir, 'seriesNum': seriesNum, \
              'channelNum': channelNum }
    
    xml = minidom.parse(getPropsPath(baseDir, seriesNum))
    
    dimensions = xml.getElementsByTagName('DimensionDescription')
    
    props['numFrames'] = int(dimensions[2].attributes['NumberOfElements'].value)
    
    def parse_pixel(index): 
        return float(dimensions[index].attributes['Voxel'].value)
    props['pixel'] = (parse_pixel(0), parse_pixel(1))
    
    duration = float(dimensions[2].attributes['Length'].value)
    props['dt'] =  duration/props['numFrames']
    
    timeString = xml.getElementsByTagName('StartTime')[0].firstChild.data
    props['startTime'] = timeString.split(' ')[1]
    
    return props

def getPropsPath(baseDir, seriesNum):
    """ Returns the path to the props file based on the series number. """
    fileName = 'Series{0:03d}_Properties.xml'.format(seriesNum)
    path = os.path.join(baseDir, fileName)
    assert os.path.isfile(path), 'Properties file does not exist!'
    return path

def getTifPath(props, imageNum):
    """ Returns the path to the .tif file based on the image num. """
    if props['numFrames'] < 100:
        numDigits = 2
    else:
        numDigits = 3
    imageStr = str(imageNum).zfill(numDigits)
    fileName = 'Series{0:03d}_t{1}_ch{2:02d}.tif'.format(props['seriesNum'], \
                imageStr, props['channelNum'])
    return os.path.join(props['baseDir'], fileName)

def printImageNames(props):
    """ Prints all of the image paths to stdout. """
    for n in range(0, props['numFrames']):
        print getTifPath(props, n)

def printImportantProps(props):
    """ Prints the props needed by MATLAB. """
    print props['startTime']
    print props['dt']
    print props['pixel'][0]
    print props['pixel'][1]

def main():
    usage = 'usage: %prog <baseDir> <seriesNum> <channelNum>'
    
    parser = OptionParser(usage)
    (options, args) = parser.parse_args()
    
    if len(args) != 3:
        parser.print_usage()
        sys.exit(0)
    
    baseDir = args[0]
    seriesNum = int(args[1])
    channelNum = int(args[2])
    
    props = readProps(baseDir, seriesNum, channelNum)
    
    printImportantProps(props)
    printImageNames(props)
    
if __name__ == '__main__': main()
