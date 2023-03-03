#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 05:10:32 2023

@author: cameronkopp
"""

from numpy import *
from numpy import random as r
from math import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import colors

# =============================================================================
# 2d WHAM CONVERGENCE CLASS
# =============================================================================

class WhamConvergence2D:
    """
    Contains a set of methods designed for optimizing convergence of the 
    whieghted histogram analysis method equations, which has proven to be
    a bottle neck in MD simulations. Hopefully you find these useful!
    """
    
    def __init__(self, xmin, xmax, ymin, ymax, simNumberX, simNumberY, 
                 binSizeX, binSizeY, springConstant, sampleNumber,
                 txtfile):
        """
        

        Parameters
        ----------
        xmin : float
            min of x reaction coordinate
        xmax : float
            max of x reaction coordinate
        ymin : float
            min of y reaction coordinate
        ymax : float
            max of y reaction coordinate
        simNumberX : int
            number of simulations to run over x coordinate range, 
            will be equally spaced 
        simNumberY : int
           number of simulations to run over y coordinate range, 
           will be equally spaced), 
        binSizeX: (distance spanned by each bin over x coordinate, smaller 
                 binSize means more precision and a larger file),
        binSizeY: (distance spanned by each bin over y coordinate, smaller 
                 binSize means more precision and a larger file),
        springConstant: (wieghted MD simulations commonly refer to the biasing
                        potential coefficient as the spring constant,
                        here we use it to calculate the standard deviation of
                        our gaussian) *NOTE* to simplify the nature of the 
                        Gaussian probabilty function, the same spring constant
                        applies to both reaction coordinates,
        sampleNumber: (number of samples taken from each simulation, i)
        
        txtfile: string of file from current directory

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
        self.Ni = sampleNumber
        
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
        
        #also need to initalize the baising potential matrix
        #and vector for set fi
        
        self.Cikl = zeros((self.simiX.size*self.simiY.size, self.biX.size, 
                           self.biY.size),dtype=float)
        self.fi = zeros((self.simiY.size*self.simiX.size),dtype=float)
        
        #load in histogram data from txt file
        with open(txtfile) as myFile:
            line = myFile.readline()
            # print(line.split())
            ix = 0
            while line:
                line = myFile.readline()
                row = line.split()
                # print(row)
                if (len(row) != 0):
                    maxyindex = len(row)
                    biy = arange(0,maxyindex)
                    for iy in biy:
                        self.hist[ix,iy] = int(row[iy])
                ix = ix + 1
                
        plt.contourf(self.biX, self.biY, self.hist)
         
        plt.show()
        
        return
    
    
    def setBiasingMatrix(self):
        """
        Stores the value of the biasing potential at each bin for each 
        simulation in self.Cikl
        """
        #shorten k
        k = self.K/2
        #cannot figure out a better way of dealing with the one d sim index i
        #the baising matrix is going to be a huge array
        i = 0
        for m1 in self.muX:
            for m2 in self.muY:
                for k in self.biX:
                    for l in self.biY:
                        self.Cikl[i,k,l] = exp(-k(((self.bcX[k]-m1)**2)+
                                                  ((self.bcY[l]-m2)**2)))
                i = i + 1
                
        return
    
    
        
#%%
# =============================================================================
# TEST CELL
# =============================================================================

filename ='2d-data-files/xmn0_xmx3_ymn0_sNumX3_sNumY30_bsX30_bsY0.01_k0.01_Ni4_5000.txt'
x0 = 0
xN = 3
y0 = 0
yN = 3
snX = 30
snY = 30
bsX =1/100
bsY =1/100
k = 4
samps = 5000

ct1 = WhamConvergence2D(x0, xN, y0, yN, snX, snY, bsX, bsY, k, samps, filename)
                
        
        