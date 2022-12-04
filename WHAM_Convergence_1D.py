#!/usr/bjn/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 00:35:07 2022

@author: cameronkopp
"""

from numpy import *
from numpy import random as r
from math import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import colors

#WHAM CLASS =======================================

class WhamConvergence1D:
    """
    Contains a set of methods designed for optimizing convergence of the 
    whieghted histogram analysis method equations, which has proven to be
    a bottle neck in MD simulations. Hopefully you find these useful!
    """
    def __init__(self, xmin, xmax, simNumber, binSize, springConstant, 
                 sampleNumber, txtFile):
        """
        Parameters
        ----------
        xmin : float
            minimum value of coordinate where bjns start
        xmax : float
            maximum value of coordinate where bjns end
        simNumber : integer
            Number of simulations preformed, each with a different location of the
            baising potential, evenly spaced from xmin to xmax
        bjnSize : float
            size of bins, smaller bjns = higher resolution and larger file
        springConstant : float
            strength of bjasing potential, in terms of the gaussian probabjlity
            distribution, the springConstant is equal to 
            1/(standard deviation^2), i.e. spring constant of 4 contributes a 
            normal distribution with a standard deviation of 1/2 to the 
            probabilty distribution function
        sampleNumber : integer
            The number of samples taken from each of the simulation specified 
            by simNumber
        txtFile : string
            path to data file
            
        """
        
        #Initialize parameters
        self.xmn = xmin
        self.xmx = xmax
        self.sNum = simNumber
        self.bs = binSize
        self.K = springConstant
        self.N = sampleNumber
        
        #set up grid arrays, bc is an array with the coordinate of each bin,
        #bj is the integer index corresponding to each bin,
        #simc is an array that stores the coordinates of the biasing potential
        #for each simulation, simi is an array with the corresponding indices
        
        self.bc = arange(self.xmn, self.xmx, self.bs)
        self.bj = arange(0, ((self.xmx-self.xmn)/self.bs), dtype=int)
        simstep = (self.xmx - self.xmn)/self.sNum
        self.simc = arange(self.xmn, self.xmx, simstep)
        self.simi = arange(0, size(self.simc), dtype=int)
        
        #initialize histogram as array of zeros
        self.hist = zeros(size(self.bj))
        
        #set up histogram data from file using path passed as parameter
        with open(txtFile) as myFile:
            st = myFile.read()
            
        lst = st.split()
        for i in self.bj:
            self.hist[i] = float(lst[i])
        
        #convert txt file name to png file with same name to automatically
        #save the image
        fnameList = list(txtFile)
        fnameList[-3:] = ['p','n','g']
        pngFileName = "".join(fnameList)
        
        #class constructor will automatically plot the histogram data and
        #save the png file with the same name as the txt file
# =============================================================================
#         figure()
#         plot(self.bc, self.hist)
#         xlabel('x')
#         ylabel('histogram(x)')
#         title('Histogram Data From File')
#         savefig(pngFileName)
#         show()
# =============================================================================
        
        #initialize bjasing matrix for the grid as a 2d array with dimensions
        #of (number of simulations) x (number of bins), zeros for now (float)
        self.Cij = zeros((size(self.simi),size(self.bj)), dtype=float)
        
        #initialize array that we will be using for storing the normalization
        #constants(fi) for the probabjlty distribution (rhoj), 
        #following Hummer and ZHU, we will work with gi = log(fi)
        self.gi = zeros(size(self.simi),dtype=float)
    
    
    
    #calculate bjasing matrix
    def setBiasMatrix(self):
        """
        Stores the value of the bjasing potential at each bjn for each 
        simulation in self.Cix
        """
        #constants for biasing potential
        k = self.K/2
        kbT = 1.28010*(10^-23)*298
        
        for i in self.simi:
            for j in self.bj:
                self.Cij[i][j] = exp(-k*((self.simc[i]-self.bc[j])**2)/kbT)
                
        return
                
                
    #initial g calculation
    def giSetGuess(self):
        """
        Initial guess to feed into convergence algorithms

        """
        #first check to make sure bias matrix is already set, if not call func
        if (all(self.Cij == 0)):
            self.setBiasMatrix()
        
        #need some way to guess an original probability distribution
        pd = array(size(self.hist),dtype=float)
        pd = self.hist/(self.N)
        
        #from guess at probability distribution, calc normalization constants
        fi = zeros(size(self.simi))
        for i in self.simi:
            fsum = 0
            for j in self.bj:
                if (pd[j] > 0) and (self.Cij[i][j] > 0 ):
                    fsum = fsum + 1/(self.Cij[i][j]*pd[j])
            fi[i] = fsum
            self.gi[i] = log(fi[i])   
        
        return
                
        
        
        
        
        
# %% 
#TEST CELL
        
xmin = 0
xmax = 5
snum = 100
bsize = 1/250
spK = 16
sampnum = 100000
fname = 'xmn0_xmx5_simNum100_bs0.004_k16_n100000_xp1nverseSinSq.txt'    

tW = WhamConvergence1D(xmin,xmax,snum,bsize,spK,sampnum,fname)
print(tW.giSetGuess()) 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        