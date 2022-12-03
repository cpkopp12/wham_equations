#!/usr/bin/env python3
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
            minimum value of coordinate where bins start
        xmax : float
            maximum value of coordinate where bins end
        simNumber : integer
            Number of simulations preformed, each with a different location of the
            baising potential, evenly spaced from xmin to xmax
        binSize : float
            size of bins, smaller bins = higher resolution and larger file
        springConstant : float
            strength of biasing potential, in terms of the gaussian probability
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
        #bi is the integer index corresponding to each bin,
        #simc is an array that stores the coordinates of the biasing potential
        #for each simulation, simi is an array with the corresponding indices
        
        self.bc = arange(self.xmn, self.xmx, self.bs)
        self.bi = arange(0, ((self.xmx-self.xmn)/self.bs), dtype=int)
        simstep = (self.xmx - self.xmn)/self.sNum
        self.simc = arange(self.xmn, self.xmx, simstep)
        self.simi = arange(0, size(self.simc), dtype=int)
        
        #initialize histogram as array of zeros
        self.hist = zeros(size(self.bi))
        
        #set up histogram data from file using path passed as parameter
        with open(txtFile) as myFile:
            st = myFile.read()
            
        lst = st.split()
        for i in self.bi:
            self.hist[i] = float(lst[i])
            
        figure()
        plot(self.bc, self.hist)
        xlabel('x')
        ylabel('histogram(x)')
        title('Histogram Data From File')
        show()
        
        
        
# %% 
#TEST CELL
        
xmin = 0
xmax = 5
snum = 100
bsize = 1/200
spK = 16
sampnum = 100000
fname = 'xmn0_xmx5_simNum100_bs0.005_k16_n100000_xp1nverseSinSq'    

testWham = WhamConvergence1D(xmin,xmax,snum,bsize,spK,sampnum,fname)  
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        