# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 15:08:55 2019

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
    
    def __init__(self,xmn,xmx,sn,bs,k):
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
                    



#==============================================================================
# WHAMEQS CLASS         
#==============================================================================

class WhamEqs1D(SampleConstructor):
    """
    stores variables and functions related to solving the wham eqs
    """
    def __init__(self,xmn,xmx,sn,bs,k):
        """
        constructer method: calcsulates/sets instance vars brange1, brange2,
        sframe_num
        """
        super().__init__(xmn,xmx,sn,bs,k)
        #CAN CALCULATE BAISING MATRIX, G0 HERE
        self.baisingMatrixCalc()
    
    #Biasing Matrix Calc
    def baisingMatrixCalc(self):
        """
        calc biasing matrix, initializes and calcs instance var cjx
        """
        self.cjx = zeros((size(self.simindex),size(self.binindex)),dtype=float)
        for i in self.simindex:
            for m in self.binindex:
                k = self.springk
                self.cjx[i][m] = exp(-k*(self.simx[i]-self.binx[m])**2)
        
    #INITIAL g CALC
    def initialGCalc(self):
        """
        initial self.p0(j X x),(1/fj)s = self.f0, and self.g0
        return: self.g0
        """
        self.p0 = self.hist/self.ni
        self.f0 = sum(self.p0*self.cjx,axis = 1)
        self.g0 = log(1/self.f0)
#        print(self.g0)
        return self.g0
        
    #WHAM CALC
    def whamCalc(self,g):
        """
        caclulates the wham equation in terms of gi, (Hummer, Zhu; 2012) eq 19
        args: g
        """
        sum1 = 0
        sum2 = 0
        for j in self.simindex:
            sum1 = sum1 + (self.nj[j]*g[j])
        for i in self.binindex:
            if(self.m[i] != 0):
                sum3 = 0
                for z in self.simindex:
                    sum3 = sum3 + (self.nj[z]*exp(g[z])*self.cjx[z][i])
                if(sum3 != 0):
                    sum2 = sum2 + (self.m[i]*log(self.m[i]/sum3))
        whamk = -1*(sum1 + sum2)
#        print(whamk)
        return whamk
        
    def dwhamCalc(self,g):
        """
        dwham in terms of g:(Hummer, Zhu; 2012) eq 20
        args: g
        output: dw with shape(g)
        """
        dw = zeros(size(self.simindex))
        for j in self.simindex:
            sum1 = 0
            for i in self.binindex:
                if(self.m[i] != 0):
                    sum2 = 0 
                    for z in self.simindex:
                        sum2 = sum2 + (self.nj[z]*exp(g[z])*self.cjx[z][i])
                    sum1 = sum1 +((self.m[i]*self.cjx[j][i])/sum2)
            dw[j] = ((exp(g[j])*sum1)-1)*self.nj[j]
#        print(dw)
        return dw
    
    
    #lsearchCalc
    def lSearch(self,a,t,b,i,iterationlimit):
        """
        lsearch function, needed within BFGS algorithm, a is alpha0, t is tao,
        b is beta, all parameters for running the search, i is the number of 
        iterations within BFGS, if calling outside of bfgs, i = 0, need to run
        all of the other methods first
        """
        alpha0 = a
        tao = t
        beta = b
        # self.gk and self.hk are calculated later in the BFGS loop, if this
        # is the first iteration we can set it here, this way can call it on
        # its own without BFGS, first work w/ gi, set self.gk after loop
        if (i == 0):
            gi = self.g0
            self.gk = self.g0
            self.hk = identity(size(self.simindex))
            self.dgk = self.dwhamCalc(self.g0)
        if(i != 0):
            gi = self.gk
        
        #set pk, tol, whamgi, whamtol, gt, whamgt, l=0 before loop
        pk = -1*dot(self.hk,self.dgk)
        tol = beta * dot(self.dgk.T,pk)
        whamgi = self.whamCalc(gi)
        whamtol = whamgi + (alpha0 * tol)
        gt = gi + (alpha0 * pk)
        whamgt = self.whamCalc(gt)
        l = 0
        alphal = alpha0
        while(whamgt > whamtol) and (l < iterationlimit):
            l = l + 1
            alphal1 = tao * alphal
            gt = gi + (alphal1 * pk)
            whamgt = self.whamCalc(gt)
            whamtol = whamgi + (alphal1 * tol)
            alphal = alphal1
        self.gk1 = gt
#        print(self.gk1)
        print('line search iterations: ',l)
        return gt
        
    #hk1 calc
    def hk1Calc(self):
        """
        hessian calc for each loop through bfgs, self.hk is already defined
        within line search, make gk1, delg, dgk, dgk1, deldgk also self vars
        update through each bfgs loop
        """
        ident = identity(size(self.simindex))
        
        # delgk, dgk1, deldgk
        self.delgk = self.gk1 - self.gk
        self.dgk1 = self.dwhamCalc(self.gk1)
        self.deldwhamk = self.dgk1 - self.dgk #name has to be more different
        #hk1 = hk + hcalc6 - hcalc5,
        #hcalc6 = (1 + hcalc1)*hcalc2
        #hcalc5 = (hcalc3 + hcalc4) / (deldw.T * delg)
        #hcalc1 = hc1n/hc1d = (deldwT*Hk*deldw)/(deldwT*delg)
        #hc1 is a working variable,
        hc1 = dot(self.deldwhamk.T,self.hk)
        hc1n = dot(hc1,self.deldwhamk)
        hc1d = dot(self.deldwhamk.T,self.delgk)
        hcalc1 = hc1n / hc1d
        
        #hcalc2 = hc2n/hc2d = (delg*delgT)/(delgT*deldw
        hc2n = dot(self.delgk,self.delgk.T)
        hc2d = dot(self.delgk.T,self.deldwhamk)
        hcalc2 = hc2n / hc2d
        
        #hcalc3 = (Hk*deldw*delgT)
        hc3 = dot(self.deldwhamk,self.delgk.T)
        hcalc3 = self.hk * hc3
        
        #hcalc4 = (hk*deldw*delgT)T
        hc4 = dot(self.deldwhamk,self.delgk.T)
        hcalc4t = self.hk * hc4
        hcalc4 = hcalc4t.T
        
        #hcalc5 (hcalc3 + hcalc4) / (deldw.T * delg)
        hc5n = hcalc3 + hcalc4
        hc5d = dot(self.deldwhamk.T,self.delgk)
        hcalc5 = hc5n / hc5d
        #hcalc6 = (1 + hcalc1)*hcalc2
        hcalc6 = (ident + hcalc1) * hcalc2
        
        #hk1 = hk + hcalc6 - hcacl5
        h1 = self.hk + hcalc6 - hcalc5
        self.hk1 = h1
        
        return h1 # ********GOING T0 RETURN ALL P1 ARRAYS HERE **************
        
    # BFGS method
    def gjBFGSCalc(self,llim,ltol,a0,t,b,iterationlimit,ploti):
        """
        BFGS algorithm applied to wham equations, needs all of the above 
        methods to work, need to run everything up to g0Calc before running
        this method
        ARGS: llim limits iterations, ltol is the error tolerance, (a0,t,b,
        iterationlimit) are the args for lsearch()
        """
        alpha0 = a0
        tao = t
        beta = b
        qlim = llim
        #constants for loop parameters
        lcon = llim
        tol = ltol
        q = 0 
        er = 10
        rg = zeros(self.simnum)
        plotcheck = ploti
        while (q < qlim): #and (er > tol):
            self.lSearch(alpha0,tao,beta,q,iterationlimit)
            self.hk1Calc()
#            print(self.delgk)
#            print(self.deldwhamk)
#            print(self.dgk)
#            print(self.dgk1)
#            w = sum(self.delgk)
#            er = abs(w)
            self.gk = self.gk1
            self.dgk = self.dgk1
            self.hk = self.hk1
            q = q + 1
            b = False
            for pc in plotcheck:
                if(pc == q):
                    ag = (sqrt(2*self.springk)/sqrt(2*math.pi))*(1/(self.xmax**2))
                    f0 = exp(self.gk)
                    self.rho = zeros(size(self.binindex))
                    for ai in self.binindex:
                        sum1 = 0
                        for ji in self.simindex:
                            sum1 = sum1 + self.nj[ji]*self.cjx[ji][ai]*f0[ji]
                        self.rho[ai] = (self.m[ai]/ag)/sum1
                        
                    figure()
                    plot(self.binx,self.rho)
                    xlabel('x')
                    ylabel('p(x)')
                    title('p(x), BFGS iteration = {}'.format(q))
                    show()  
                    w = self.whamCalc(self.gk)
                    print('BFGS iteration ',q,' wham ',w)
                    b = True
            
            rg = self.gk
            if b == False:
                w = self.whamCalc(rg)
                print('BFGS iteration ',q,' wham ',w)
#            print(self.delgk)
            
        
        return rg
        
    def iterativeSolve(self,ilimit):
        i = 0
        self.p0 = self.hist/self.ni
        f0 = zeros(size(self.simindex))
        f0 = 1/sum(self.p0*self.cjx,axis = 1)
        f1 = zeros(size(self.simindex))
        while (i < ilimit):
            i = i + 1
            for j in self.simindex:
                sum1= 0    
                for bi in self.binindex:
                    sum2 = 0
                    for j2 in self.simindex:
                        sum2 = sum2 + self.nj[j2]*f0[j2]*self.cjx[j2][bi]
                    sum1 = sum1 + self.m[bi]*self.cjx[j][bi]/sum2
                f1[j] = 1/sum1
#            print('f1 ', f1)
            g = log(f1)
            w1 = self.whamCalc(g)
            print('i ', i,' wham ',w1)
            f0 = f1
        ag = (sqrt(2*self.springk)/sqrt(2*math.pi))*(1/(self.xmax**2))
        self.rho = zeros(size(self.binindex))
        for ai in self.binindex:
            sum1 = 0
            for ji in self.simindex:
                sum1 = sum1 + self.nj[ji]*self.cjx[ji][ai]*f0[ji]
            self.rho[ai] = (self.m[ai]/ag)/sum1
            
        figure()
        plot(self.binx,self.rho)
        xlabel('x')
        ylabel('p(x)')
        title('p(x), BFGS iteration = {}'.format(i))
        show() 
            
            




#==============================================================================
# TESTS OF COMBINED CLASSES
#==============================================================================

#SET UP SAMPLE DATA  
xmin = 0
xmax = 5
simnum = 100
binsize = 1/100
k = 1
exampleSample = WhamEqs1D(xmin,xmax,simnum,binsize,k)
n = 5000
exampleHistogram = exampleSample.sconXCubedSinSquared2PiX(n)



##SET UP FOR AND CALL BFGS METHOD
a0c = 2
betac = 0.01
taoc = 0.5
il = 150
lsl = 200
etol = 0
plotcheck = [2,5,10,25,50,100,150]
#CALL initialGCALC FIRST, the gjBFGSCalc()
g0 = exampleSample.initialGCalc()
exampleG = exampleSample.gjBFGSCalc(il,etol,a0c,taoc,betac,lsl,plotcheck)

# =============================================================================
# figure()
# plot(exampleSample.binx,exampleSample.m)
# xlabel('x')
# ylabel('total bin count over all simulations')
# title('m(x) for exampleSample')
# show() 
# =============================================================================






