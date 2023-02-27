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
from mpl_toolkits import mplot3d

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
        sampleNumber (number of samples taken from each simulation, i)
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
        
        return
        
        
        
        
    def biasPotential(self, xt, yt, mx, my):
        """
        Parameters
        ----------
        xt : float
            value of x reaction cordinate randomly generated from 
            gaussian dist
        yt : float
            value of y reaction cordinate randomly generated from 
            gaussian dist
        mx: float
            value of the center of the potential along x reaction 
            coordinate
        mx: float
            value of the center of the potential along y reaction 
            coordinate

        Returns
        -------
        returns the probability distribution value associated with the
        harmonic biasing potential for the given test points

        """
        #constants: c1= standard deviation in terms of spring constant k
        #nc = normalization constant in two dimensions of equal 
        #   standard deviations (factor of 2 carried in with c1)
        c1 = self.K/2
        nc = c1/math.pi
        
        #contribution to wheighted prob/ dist/ function from biasing
        baisDist = nc *exp(-c1*(((xt-mx)**2)+(yt-my)**2))
        
        return baisDist
    
    def smoothTestFunction(self, xt, yt):
        """
        Will serve as our known unbiased probability distribution
        
        Need to make sure this function is always less than to use 
            accept-reject method

        Parameters
        ----------
        xt : float
            value of x reaction cordinate randomly generated from 
            gaussian dist
        yt : float
            value of y reaction cordinate randomly generated from 
            gaussian dist

        Returns
        -------
        The value of a messy and arbitrary function selected for its smoothness
        and that it is bounded by (0,1)

        """
        #NOTE THAT NORMALIZATION CONSTAN IS GENERATED NUMERICALLY
        #   BY INTEGRATING OF MIN AND MAX OF REACTION COORDS
        nc = 1/2.51388
        
        unbiasedDist = ((1/(1+(1-xt)**2)) * (((math.sin(3*xt)**2)+2)/3) * 
                        (1/(1+(2-yt)**2)) * (((math.sin(3*yt)**2)+2)/3))
        
        return unbiasedDist
    
    
    def dateGen2D(self):
        """
        Relies on a accept-reject algorithm, applied a few techniques to try 
        to speed it up but for large systems it will take some time.

        Returns
        -------
        A histogram corresponding to the binning of the samples generated 
        through this mock biased MD simulation along two reaction coordinates,
        it will be a [[biX.size X biY.size] 2d array of integers, and will be
        the resulting sumation over all of the simulations

        """
        
        #func for random gaussian points wants K -> sigma
        stdev = 1/sqrt(self.K)
        #covariance matrix
        cvmatrix =[[stdev**2,0],[0,stdev**2]];
        
        #each iteration through double for loop coresponds to the simulation
        #with the biasing potential centered at (mX,mY)
        
        for mX in self.muX:
            
            for mY in self.muY:
                
                print('(muX,muY) is ({},{})'.format(mX,mY))
                
                #reset accepted sample counter
                i = 0
                
                #loop until i reaches self.Ni = sample num per sim
                while(i < self.Ni):
                    
                    #generate test points for sample
                    xT, yT = r.multivariate_normal((mX, mY), cvmatrix)
                    
                    tval = (self.biasPotential(xT, yT, mX, mY) * 
                            self.smoothTestFunction(xT, yT))
                    #We know that tval will always be bounded by
                    #the value of the biasPotential dist alone,
                    #so the accept-reject method supports what follows
                    uniformT = r.uniform(0, self.biasPotential(xT, yT, mX, mY))
                    
                    #accept-reject condition
                    if (tval > uniformT):
                        #also check to make sure points land in histogram
                        if (xT < self.xmx) and (xT > self.xmn):
                            if (yT < self.ymx) and (yT > self.ymn):
                                #now we accept the point as valid sample
                                i = i + 1
                                #need to match xT, yT floats to a bin index
                                histBinFloatX = (xT/self.bsX)
                                histBinX = floor(histBinFloatX)
                                indexX = int(histBinX)
                                
                                histBinFloatY = (yT/self.bsY)
                                histBinY = floor(histBinFloatY)
                                indexY = int(histBinY)
                                
                                #update histogram
                                self.hist[indexX, indexY] = (
                                    self.hist[indexX, indexY] + 1)
                                
                                #print a message to track progress
                                if (i == floor(self.Ni/2)):
                                    print(
                                        'half done on ({},{})'.format(mX,mY))
                    
        plt.contourf(self.biX, self.biY, self.hist)
        
        plt.show()
        
        self.writeToFile()
        
        return self.hist
    
    
    def writeToFile(self):
        filename = ('2d-data-files/xmn{}_xmx{}_ymn{}_sNumX{}_sNumY{}_bsX{}_bsY{}_k{}_Ni{}_{}.txt'
                    .format(self.xmn, self.xmx, self.ymn, self.ymx, self.sNumX,
                            self.sNumY, self.bsX, self.bsY, self.K, self.Ni))
        with open(filename,'w') as myfile:
            for x in self.biX:
                myfile.write('\n');
                for y in self.biY:
                    myfile.write('{} '.format(self.hist[x,y]))
                    
        myfile.close()
        
        return
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
#%%
# =============================================================================
# TEST Cell
# =============================================================================


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

t1 =   SampleConstructor2D(x0, xN, y0, yN,snX,snY,bsX,bsY,k,samps) 
t1.dateGen2D()
        
        
        
        
        
        
        
        
        
        
        
        