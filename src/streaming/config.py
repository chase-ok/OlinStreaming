'''
Created on Mar 20, 2012

@author: ckernan
'''

import os

defaultRoot = "/home/ckernan/data/streaming"

densities = { "max": "Density 1", 
              "high": "High density", 
              "two-thirds": "Twothirds High density",
              "one-third": "Onethird high density",
              "one-fourth": "Onequarter high density",
              "one-eighth": "One eighth high density",
            }
for density, folder in densities.iteritems():
    densities[density] = os.path.join(defaultRoot, folder)