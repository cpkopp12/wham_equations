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
        figure()
        plot(self.bc, self.hist)
        xlabel('x')
        ylabel('histogram(x)')
        title('Histogram Data From File')
        savefig(pngFileName)
        show()
        
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
        
        for i in self.simi:
            for j in self.bj:
                self.Cij[i][j] = exp(-k*((self.simc[i]-self.bc[j])**2))
                
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
    
    #WHAM optimization function calc
    def optFuncCalc(self, g):
        """
        calculates the optimization function for the WHAM equations, eq(19)
        of Hummer and Zhu,
        parameter g comes from the fact that in the line search we will be 
        passing test values for g, and we can always pass in self.gi as this
        parameter,
        returns a float
        """
        
        sum1 = 0
        sum2 = 0
        
        for i in self.simi:
            sum1 = sum1 + (self.N*g[i])
            
        for j in self.bj:
            if (self.hist[j] != 0):
                sum3 = 0
                for i in self.simi:
                    sum3 = sum3 + (self.N*self.Cij[i][j]*exp(g[i]))
                if (sum3 != 0):
                    sum2 = sum2 + (self.hist[j]*log(self.hist[j]/sum3))
        
        #Result of the calculation is the negative sum of the two summations
        optFunc = -sum1-sum2
        
        return optFunc   
    
    #Line search requires the derivative of the opt fuct with respect to gi
    def dgiOptFuncCalc(self, g):
        """
        calculates the derivative of the optimization function with respect to
        the elements of g, so this function returns an array floats of size
        self.gi
        """
        #initialize return array
        dOptFunc = zeros(size(self.gi))
        
        for i in self.simi:
            sum1 = 0 
            for j in self.bj:
                if (self.hist[j] != 0):
                    sum2 = 0
                    for i2 in self.simi:
                        sum2 = sum2 + (self.N*exp(g[i2]*self.Cij[i2][j]))
                    if (sum2 != 0):
                        if (j == 5) and (i == 1):
                            print('sum2')
                            print(sum2)
                        sum1 = sum1 + ((self.hist[j]*self.Cij[i][j])/sum2)
            #calculate the ith element of the derivaitive of the 
            # optimization function            
            dOptFunc[i] = (self.N * (exp(g[i]) * sum1)) - self.N
                           
        
        return dOptFunc
    
    
    
    #Backtracking Armijo Line Search
    def lineSearch(self, alpha0, tao, beta, hessian, dgiOptFunc,
                   bfgsLoopIteration, iterationLimit):
        """
        alpha0 : float
            Initial length of line search along search direction
        tao : float in (0,1)
            Length of line search decreases by this factor on each iteration
        beta : float in (0,1)
            coeficient in front of the Armijo condition, smaller factor means
            the function decreases less at each step, but a larger factor may
            make the number of iterations required at each step impractical
        hessian : matrix(size(self.gi) X size(self.gi)))
            used in BFGS algorithm to calculate search direction
        dgiOptFunc : array(size(self.gi))
            derivative of the optimization constant, calculated with
            self.dgiOptFuncCalc() method
        bfgsLoopIteration: it is useful to know what iteration the outer
            BFGS loop is up to, incase we want to periodically check contents
            of this function
        iterationLimit : int
            maximum number of backtracking iterations before the function
            returns a value which does not satisfy the armijo condition

        Returns
        -------
        self.gi which is closer to minimizing the optimation function than 
        before, because we're using the self.gi we don't actually return
        it from the function,
        

        """
        #initialize parameters
        a0 = alpha0
        t = tao
        b = beta
        H = hessian
        dgiA = dgiOptFunc
        bfgsi = bfgsLoopIteration
        iLim = iterationLimit
        gi0 = self.gi
        
        #pk is the line search vector
        pk = -1*dot(H,dgiA)
        print('pk')
        print(pk)
        
        #tolerance for armijo condition: f(x+(alpha * pk)) < f(x) + tol
        tol = b * dot(dgiA.T,pk)
        A0 = self.optFuncCalc(gi0)
        ATol = A0 + (a0 * tol)
        
        #calculate test gi by adding (alpha*pk)
        giTest = gi0 + (a0*pk)
        print('gi0: ',gi0)
        print('giTest: ', giTest)
        ATest = self.optFuncCalc(giTest)
        print('ATol0: ', ATol)
        print('ATest: ', ATest)
        #loop index
        l = 0
        al = a0
        while(ATol < ATest) and (l < iLim):
            l = l + 1
            alp1 = t * al
            giTest = gi0 + (alp1 * pk)
            ATest = self.optFuncCalc(giTest)
            ATol = A0 + (alp1 * tol)
            al = alp1
            print('l: ', l)
            print('ATest: ', ATest)
            print('ATol: ', ATol)
        
        print('gi0: ', gi0)
        print('giTest: ', giTest)
        print('diff: ',gi0-giTest)
        
        self.gi = giTest
        
        
        
        return
        
                                    
                
        
        
        
        
        
# %% 
#TEST CELL
        
# =============================================================================
# xmn = 0
# xmx = 5
# simnum = 10
# binsize = 1/20
# spK = 4
# sampnum = 1000
# fname = 'xmn0_xmx5_simNum10_bs0.05_k4_n1000_xp1nverseSinSq.txt'  
# =============================================================================

xmn = 0
xmx = 5
simnum = 200
binsize = 1/500
spK = 16
sampnum = 100000  
fname = "xmn0_xmx5_simNum200_bs0.002_k16_n100000_xp1nverseSinSq.txt"

# =============================================================================
# xmn = 0
# xmx = 5
# simnum = 50
# binsize = 1/125
# spK = 9
# sampnum = 10000
# fname = 'xmn0_xmx5_simNum50_bs0.008_k9_n10000_xp1nverseSinSq.txt'
# =============================================================================

tW = WhamConvergence1D(xmn,xmx,simnum,binsize,spK,sampnum,fname)
tW.giSetGuess()
optCalcTest = tW.optFuncCalc(tW.gi)
print(optCalcTest)
dOptTest = tW.dgiOptFuncCalc(tW.gi)
print(dOptTest)
print(tW.gi)

hess0 = identity(size(tW.gi))
a0c = 2
betac = 0.000001
taoc = 0.5
il = 15000

tW.lineSearch(a0c, taoc, betac, hess0, dOptTest, 0, il)

   
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        