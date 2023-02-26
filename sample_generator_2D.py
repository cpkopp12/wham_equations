#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 21:12:20 2023

@author: cameronkopp
"""
from numpy import *
from numpy import random as r
from math import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import colors

# =============================================================================
# SAMPLE CONSTRUCTOR CLASS 
# =============================================================================
class SampleConstructor2D:
    """
        inputs: xmin, xmax (such that the sample will span range xmin-xmax 
                            along the x reaction coordinate) 
        ymin, ymax (such that the sample will span range ymin-ymax along the 
                    y reaction coordinate) 
        simNumberX (number of simulations to run over x coordinate range, 
                    will be equally spaced), 
        simNumberY (number of simulations to run over y coordinate range, 
                    will be equally spaced), 
        binSizeX (distance spanned by each bin over x coordinate, smaller 
                 binSize means more precision and a larger file),
        binSizeY (distance spanned by each bin over y coordinate, smaller 
                 binSize means more precision and a larger file),
        springConstant (wieghted MD simulations commonly refer to the biasing
                        potential coefficient as the spring constant,
                        here we use it to calculate the standard deviation of
                        our gaussian) *NOTE* to simplify the nature of the 
                        Gaussian probabilty function, the same spring constant
                        applies to both reaction coordinates,
        sampleNumber (number of samples taken from each simulation)
    """
    
    def __init__(self, xmin, xmax, ymin, ymax, simNumberX, simNumberY, binSizeX,
                 binSizeY, springConstant, sampleNumber):
        
        """
        Class Constructor: initalizes paramaters and sets up the necessary
        grid for the histogram.
        """
        #Read in parameters
        self.xmn = xmin
        self.xmx = xmax
        self.ymn = ymin
        self.ymx = ymax
        self.sNumX = simNumberX
        self.sNumY = simNumberY
        self.bsX = binSizeX
        self.bsY = binSizeY
        self.K = springConstant
        self.N = sampleNumber
        
        #constuct coordinate system, don't need to worry about combining
        #reaction coordinates until setting up binned histogram, will always
        #have indepent x,y coordinates + indecies
        #bin coordinates MUST be in center of bin
        
        #x reaction coordinate
        self.bcX = arange(self.xmn, self.xmx, self.bsX)+(self.bsX/2)
        self.biX = arange(0, ((self.xmx-self.xmn)/self.bsY), dtype=int)
        
        #xy reaction coordinate
        self.bcY = arange(self.ymn, self.ymx, self.bsY)+(self.bsY/2)
        self.biY = arange(0, ((self.ymx-self.ymn)/self.bsY), dtype=int)
        
        #construct the coordinates for the center of the biasing potential
        #for each simulation
        #for now, going to treat the coordinates of the baising potential
        #seperately along each reaction coordinate, though in the WHAM 
        #convergence we will need to combine thes in a [1X(simNumX x simNumY)]
        #array
        simstepX = (self.xmx-self.xmn)/self.sNumX
        simstepY = (self.ymx-self.ymn)/self.sNumY
        #coordinate of each baising potential
        self.muX = arange(self.bcX[0],self.bcX[self.bcX.size-1],simstepX)
        self.muY = arange(self.bcY[0],self.bcY[self.bcY.size-1],simstepY)
        #index for each sim along each coordinate axis
        self.simiX = arange(0, self.muX.size, dtype=int)
        self.simiY = arange(0, self.muY.size, dtype=int)
        
        #initalize histogram as zero (int) size of [biX.size X biY.size]
        histShape = (self.biX.size, self.biY.size)
        self.hist = zeros(histShape, dtype=int)
        
        
#%%
# =============================================================================
# TEST Cell
# =============================================================================


x0 = 0
xN = 3
y0 = 0
yN = 3
snX = 5
snY = 5
bsX =.5
bsY =.5
k = 1
samps =9

t1 =   SampleConstructor2D(x0, xN, y0, yN,snX,snY,bsX,bsY,k,samps) 
print(t1.bcY) 

print(t1.muX)    
print(t1.hist)
        
        
        
        
        
        
        
        
        
        
        
        
        