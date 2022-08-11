#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:34:36 2022

@author: cameronkopp
"""

from numpy import *
from numpy import random as r
from math import *
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import colors



#==============================================================================
# SAMPLE CONSTRUCTION CLASS
#==============================================================================

class SampleConstructor:
    """
    inputs: xmin, xmax, simstep, binsize
    p(x) must be defined as < 1 for all x in xmin to x max
    """
    
    def __init__(self,xmn,xmx,sn,bs,k,fname):
        """
        constructor for class, set self. xmin,xmax,simstep,binsize,k(SPRING
        CONSTANT,USED FOR sigma IN GAUSSIAN)
        """
        # READ ARGS AND SET AS INStANCE VARS
        self.xmin = xmn
        self.xmax = xmx
        self.simnum = sn
        self.binsize = bs
        self.springk = k
        self.filename = fname
        
        # CALCULATE BIN COORDINATES
        self.binx = arange(self.xmin,self.xmax,self.binsize)
        self.binindex = arange(0,((self.xmax-self.xmin)*(1/self.binsize)),dtype=int)
        
        # Calculate bias coorndinates for simstep
        self.simstep = (self.xmax - self.xmin)/(self.simnum)
        self.simx = arange(self.xmin,self.xmax,self.simstep)
        self.simindex = arange(0,size(self.simx),dtype = int)
        
        # INITIALIZE histogram AS ZEROS AND SUM j OF histrogram: m(x)
        self.hist = zeros((size(self.simx),size(self.binindex)))
        self.m = zeros(size(self.binindex))
        
#==============================================================================
#     def sconx2(self,n):
#         """
#         return histogram of sample for p(x) = Cx**2, C = normalization constant
#         ars: ni is the number of samples
#         """
#         # INITIALIZE ni AS INSTANCE VAR
#         histfix = self.xmin / self.binsize
#         self.ni = n
#         self.nj = zeros(size(self.simx))
#         third = int(self.ni/3)
#         u = int(2*self.ni/2)
#         ag = (1/(sqrt(2*math.pi)*sqrt(2*self.springk)))
#         ap = (self.xmax**3)-(self.xmin**3)
#         #sample values
#         self.xlist = zeros((size(self.simindex),self.ni))
#         #START LOOP, OUTER DEFINES SIM NUMBER J
#         for j in self.simindex:
#             print(j)
#             # reset sample counter for each sim
#             i = 0
#             z = 0 
#             # center gaussian at simx[j]
#             mu = self.simx[j]
#             #WHILE LOOP GENERATES SAMPLES NI FOR SIM J
#             while(i < self.ni): #and (z < 100000):
# #                z = z + 1
#                 # u distributed test sample point
#                 xu = r.uniform(0,1)
#                 xt = (ap*xu)**(1/3)
#                 # MAKE SURE TEST SAMPLE POINT IN RANGE XMIN,XMAX
#                 pxt = (3/ap)*xt**2
#                 #value of p(x)*bias
#                 pgxt = (3/ap)*ag * exp(-self.springk*((xt-mu)**2))*(xt**2)
# #                         g = exp(-self.springk*((xt-mu))**2)
# #                         print(g)
#                 # uniform random 0 to gxt
#                 y = r.uniform(0,pxt)
#                 # test if pbxt < y, if it is save sample, if not back to 
#                 # top of the loop
#                 t = y - pgxt
#                 if t < 0:
#                     if (xt<self.xmax) and (xt>self.xmin):
#                         self.xlist[j][i] = xt
#                         histval = (xt * (1/self.binsize))-histfix
#                         h = floor(histval)
#                         ih = int(h)
#                         self.hist[j][ih] = self.hist[j][ih] + 1
#                         self.nj[j] = self.nj[j] + 1
#                     i = i + 1
# #                             print(i)
#                     if(i == third):
#                         print('1/3')
#                     if(i == u):
#                         print('2/3')
#         for j in self.simindex:
#             figure()
#             plot(self.binx,self.hist[j])
#             xlabel('x')
#             ylabel('y')
#             title('m')
#             show()
#             
#         for x in self.binindex:
#             for j in self.simindex:
#                 self.m[x] = self.m[x] + self.hist[j][x]
#         figure()
#         plot(self.binx,self.m)
#         xlabel('x')
#         ylabel('y')
#         title('m')
#         show()
#         return self.xlist,self.hist
#         
#     def sconx3(self,n):
#         """
#         return histogram of sample for p(x) = Cx**2, C = normalization constant
#         ars: ni is the number of samples
#         """
#         # INITIALIZE ni AS INSTANCE VAR
#         histfix = self.xmin / self.binsize
#         self.ni = n
#         self.nj = zeros(size(self.simx))
#         third = int(self.ni/3)
#         u = int(2*self.ni/2)
#         ag = (1/(sqrt(2*math.pi)*sqrt(2*self.springk)))
#         ap = (self.xmax**4)-(self.xmin**4)
#         #sample values
#         self.xlist = zeros((size(self.simindex),self.ni))
#         #START LOOP, OUTER DEFINES SIM NUMBER J
#         for j in self.simindex:
#             print(j)
#             # reset sample counter for each sim
#             i = 0
#             z = 0 
#             # center gaussian at simx[j]
#             mu = self.simx[j]
#             #WHILE LOOP GENERATES SAMPLES NI FOR SIM J
#             while(i < self.ni): 
#                 
#                 xu = r.uniform(0,1)
#                 xt = (ap*xu)**(1/4)
#                 # MAKE SURE TEST SAMPLE POINT IN RANGE XMIN,XMAX
#                 pxt = (4/ap)*xt**4
#                 #value of p(x)*bias
#                 pgxt = (1/ap)*ag * exp(-self.springk*((xt-mu)**2))*(xt**3)
# #                         g = exp(-self.springk*((xt-mu))**2)
# #                         print(g)
#                 # uniform random 0 to gxt
#                 y = r.uniform(0,pxt)
#                 # test if pbxt < y, if it is save sample, if not back to 
#                 # top of the loop
#                 t = y - pgxt
#                 if t < 0:
#                     if (xt<self.xmax) and (xt>self.xmin):
#                         self.xlist[j][i] = xt
#                         histval = (xt * (1/self.binsize))-histfix
#                         h = floor(histval)
#                         ih = int(h)
#                         self.hist[j][ih] = self.hist[j][ih] + 1
#                         self.nj[j] = self.nj[j] + 1
#                     i = i + 1
# #                             print(i)
#                     if(i == third):
#                         print('1/3')
#                     if(i == u):
#                         print('2/3')
#         for j in self.simindex:
#             figure()
#             plot(self.binx,self.hist[j])
#             xlabel('x')
#             ylabel('y')
#             title('m')
#             show()
#             
#         for x in self.binindex:
#             for j in self.simindex:
#                 self.m[x] = self.m[x] + self.hist[j][x]
#         figure()
#         plot(self.binx,self.m)
#         xlabel('')
#         ylabel('y')
#         title('m')
#         show()
#         return self.xlist,self.hist
#==============================================================================
        
        
    def sconXSquaredSinSquaredX(self,n):
        """
        return histogram of sample for p(x) = Cx**2, C = normalization constant
        ars: ni is the number of samples
        """
        # INITIALIZE ni AS INSTANCE VAR
        histfix = self.xmin / self.binsize
        self.ni = n
        self.nj = zeros(size(self.simx))
        third = int(self.ni/3)
        u = int(2*self.ni/2)
        ag = (sqrt(2*self.springk)/sqrt(2*math.pi))*(1/(self.xmax**2))
        #sample values
        k = self.springk
        self.xlist = zeros((size(self.simindex),self.ni))
        #START LOOP, OUTER DEFINES SIM NUMBER J
        for j in self.simindex:
#            if j == 25:
#                print(j)
#            if j == 50:
#                print(j)
#            if j == 75:
#                print(j)
            print(j)
            # reset sample counter for each sim
            i = 0
            # center gaussian at simx[j]
            mu = self.simx[j]
            #WHILE LOOP GENERATES SAMPLES NI FOR SIM J
            
            while(i < self.ni):
                #generate test point, xt
                xt = r.uniform(self.xmin,self.xmax)
                fmax = 1
                #value of p(x)*bias MUST BE DEFINED AS < 1 on xmin to xmax
                pjxt = ag*(xt**2)*exp(-(k/2)*((xt-mu)**2))*(math.sin(xt)**2)
                # uniform random 0 to fmax
                y = r.uniform(0,fmax)
                #test acceptance, if  y (uniform(0,1)) < pj(xt), accept xt as
                #valid sample point, if not reject and start again
                t = y - pjxt
                if t < 0:
                    if (xt<self.xmax) and (xt>self.xmin):
                        #modify list
                        self.xlist[j][i] = xt
                        #modify histogram
                        histval = (xt * (1/self.binsize))-histfix
                        h = floor(histval)
                        ih = int(h)
                        self.hist[j][ih] = self.hist[j][ih] + 1
                        #total sample points for sim j, important: some
                        #accepted xts dont fall in xmin, xmax, need this for 
                        #wham to work 
                        self.nj[j] = self.nj[j] + 1
                    i = i + 1
#                    print(i)
                    if(i == third):
                        print('1/3')
                    if(i == u):
                        print('2/3')
                        
                        
#        for j in self.simindex:
#            figure()
#            plot(self.binx,self.hist[j])
#            xlabel('x')
#            ylabel('y')
#            title('m')
#            show()
            
        for x in self.binindex:
            for j in self.simindex:
                self.m[x] = self.m[x] + self.hist[j][x]
        figure()
        plot(self.binx,self.m)
        xlabel('x')
        ylabel('m(x)')
        title(' ')
        show()
        return self.hist
        
    def sconXCubedSinSquared2PiX(self,n):
        """
        return histogram of sample for p(x) = Cx**2, C = normalization constant
        ars: ni is the number of samples
        """
        # INITIALIZE ni AS INSTANCE VAR
        pi = math.pi
        histfix = self.xmin / self.binsize
        self.ni = n
        self.nj = zeros(size(self.simx))
        third = int(self.ni/3)
        u = int(2*self.ni/2)
        ag = (sqrt(2*self.springk)/sqrt(2*math.pi))*(1/((self.xmax+2)**3))
        #sample values
        k = self.springk
        self.xlist = zeros((size(self.simindex),self.ni))
        #START LOOP, OUTER DEFINES SIM NUMBER J
        for j in self.simindex:
            print(j)
            # reset sample counter for each sim
            i = 0
            z = 0 
            # center gaussian at simx[j]
            mu = self.simx[j]
            #WHILE LOOP GENERATES SAMPLES NI FOR SIM J
            
            while(i < self.ni):
                #generate test point, xt
                xt = r.uniform(self.xmin,self.xmax)
                fmax = 1
                #value of p(x)*bias MUST BE DEFINED AS < 1 on xmin to xmax
                pjxt = ag*((xt+2)**3)*exp(-(k/2)*((xt-mu)**2))*(math.sin(2*pi*xt)**2)
                # uniform random 0 to fmax
                y = r.uniform(0,fmax)
                #test acceptance, if  y (uniform(0,1)) < pj(xt), accept xt as
                #valid sample point, if not reject and start again
                t = y - pjxt
                if t < 0:
                    if (xt<self.xmax) and (xt>self.xmin):
                        #modify list
                        self.xlist[j][i] = xt
                        #modify histogram
                        histval = (xt * (1/self.binsize))-histfix
                        h = floor(histval)
                        ih = int(h)
                        self.hist[j][ih] = self.hist[j][ih] + 1
                        #total sample points for sim j, important: some
                        #accepted xts dont fall in xmin, xmax, need this for 
                        #wham to work 
                        self.nj[j] = self.nj[j] + 1
                    i = i + 1
#                             print(i)
                    if(i == third):
                        print('1/3')
                    if(i == u):
                        print('2/3')
                        
                        
#        for j in self.simindex:
#            figure()
#            plot(self.binx,self.hist[j])
#            xlabel('x')
#            ylabel('y')
#            title('m')
#            show()
            
        for x in self.binindex:
            for j in self.simindex:
                self.m[x] = self.m[x] + self.hist[j][x]
        figure()
        plot(self.binx,self.m)
        xlabel('x')
        ylabel('m(x)')
        title(' ')
        show()
        
        with open(self.filename,'w') as myfile:
            for x in self.binindex:
                myfile.write('{} '.format(self.m[x]))
            
        
        return self.hist
        
#==============================================================================
#     def sconsinxfourth(self,n):
#         """
#         return histogram of sample for p(x) = Cx**2, C = normalization constant
#         ars: ni is the number of samples
#         """
#         # INITIALIZE ni AS INSTANCE VAR
#         histfix = self.xmin / self.binsize
#         self.ni = n
#         self.nj = zeros(size(self.simx))
#         third = int(self.ni/3)
#         u = int(2*self.ni/2)
#         ag = (1/(sqrt(2*math.pi)*sqrt(2*self.springk)))
#         ap = (self.xmax**4)-(self.xmin**4)
#         #sample values
#         self.xlist = zeros((size(self.simindex),self.ni))
#         #START LOOP, OUTER DEFINES SIM NUMBER J
#         for j in self.simindex:
#             print(j)
#             # reset sample counter for each sim
#             i = 0
#             z = 0 
#             # center gaussian at simx[j]
#             mu = self.simx[j]
#             #WHILE LOOP GENERATES SAMPLES NI FOR SIM J
#             while(i < self.ni): #and (z < 100000):
# #                z = z + 1
#                 # u distributed test sample point
#                 xu = r.uniform(self.xmin,self.xmax)
#                 # MAKE SURE TEST SAMPLE POINT IN RANGE XMIN,XMAX
#                 pxt = 2
#                 #value of p(x)*bias
#                 pgxt = ag * exp(-self.springk*((xu-mu)**2))*(math.sin(4*xu)**4)
# #                         g = exp(-self.springk*((xt-mu))**2)
# #                         print(g)
#                 # uniform random 0 to gxt
#                 y = r.uniform(0,pxt)
#                 # test if pbxt < y, if it is save sample, if not back to 
#                 # top of the loop
#                 t = y - pgxt
#                 if t < 0:
#                     if (xu<self.xmax) and (xu>self.xmin):
#                         self.xlist[j][i] = xu
#                         histval = (xu * (1/self.binsize))-histfix
#                         h = floor(histval)
#                         ih = int(h)
#                         self.hist[j][ih] = self.hist[j][ih] + 1
#                         self.nj[j] = self.nj[j] + 1
#                     i = i + 1
# #                             print(i)
#                     if(i == third):
#                         print('1/3')
#                     if(i == u):
#                         print('2/3')
#         for j in self.simindex:
#             figure()
#             plot(self.binx,self.hist[j])
#             xlabel('x')
#             ylabel('y')
#             title('m')
#             show()
#             
#         for x in self.binindex:
#             for j in self.simindex:
#                 self.m[x] = self.m[x] + self.hist[j][x]
#         figure()
#         plot(self.binx,self.m)
#         xlabel('x')
#         ylabel('y')
#         title('m')
#         show()
#         return self.xlist,self.hist
#==============================================================================
                    

#SET UP SAMPLE DATA  
xmin = 0
xmax = 5
simnum = 100
binsize = 1/100
k = 1
n = 5000
filename = '{}_{}_{}_{}_{}_{}_xcubedsinsquared.txt'.format(xmin,xmax,simnum,
                                                           binsize,k,n)
exampleSample = SampleConstructor(xmin,xmax,simnum,binsize,k,filename)
exampleHistogram = exampleSample.sconXCubedSinSquared2PiX(n)



# =============================================================================
# ##SET UP FOR AND CALL BFGS METHOD
# a0c = 2
# betac = 0.01
# taoc = 0.5
# il = 150
# lsl = 200
# etol = 0
# plotcheck = [2,5,10,25,50,100,150]
# #CALL initialGCALC FIRST, the gjBFGSCalc()
# g0 = exampleSample.initialGCalc()
# exampleG = exampleSample.gjBFGSCalc(il,etol,a0c,taoc,betac,lsl,plotcheck)
# =============================================================================

figure()
plot(exampleSample.binx,exampleSample.m)
xlabel('x')
ylabel('total bin count over all simulations')
title('m(x) for exampleSample')
show() 
