# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:11:00 2016

@author: Evyatar
"""


from math import *
import numpy as np
import scipy
from Instance import * 
from Linear import *

class InstanceLinear(Instance):
    
    def __init__(self,fr,fb,c):
        self.fr = fr
        self.fb = fb
        self.c=np.float64(c)

    def intersection(self):
        #t*ar+br=t*ab+bb
        return (self.fb.offset-self.fr.offset)/(self.fr.slope-self.fb.slope)

            

    def cr_z(self,z):
        if z==0:
            return self.fb.offset/self.fr.offset
        elif z==np.inf:
            return self.fr.slope/self.fb.slope
        else:
            if z<1:
                #interval y<z,z<=y<1,y>=1
                #case y<z
                cr_c1 = np.float64(1.0)
                #case z<=y<1,cr=f_{rb}(z)+c+f_b(y)/f_r(y)
                cr_c2 = maxLin(self.fb.slope,self.fb.offset+self.f_rb(z)+self.c,self.fr.slope,self.fr.offset,z,1.0)[0]
                #case y>1,cr=f_{rb}(z)+c+f_b(y)/f_b(y)
                cr_c3 = maxLin(self.fb.slope,self.fb.offset+self.f_rb(z)+self.c,self.fb.slope,self.fb.offset,1.0,np.inf)[0]
            else:
                #interval y<=1,1<y<z,y>=z
                #case y<=1
                cr_c1 = np.float64(1.0)
                #case 1<y<z,cr=f_r(y)/f_b(y)
                cr_c2 = maxLin(self.fr,self.fb,1.0,z)[0]
                #case y>=z,cr=f_{rb}(z)+c+f_b(y)/f_b(y)
                cr_c3 = maxLin(self.fb.slope,self.fb.offset+self.f_rb(z)+self.c,self.fb.slope,self.fb.offset,z,np.inf)[0]
            return max(cr_c1,cr_c2,cr_c3)
                

    #todo - g in linear
    def g(self,z):
        pass
    #todo - intersection linear
    def intersection_time(self):
        pass
    #todo right middle
    def middle_time(self):
        return 1.0
    def rho_inf(self):
        return maxLin(self.fr,self.fb,np.float64(1.0),np.inf)[0]
    
    def calc_t1(self,eps):
        #t1
        #f_{r}(z)+c =(1+eps)(f_{rb}(0) +c+f_b(z) 
        #frb(t1)-eps*f_b(t1)=(1+eps)f_{rb}(0)+eps*c
        
        #t1 = min(np.float64(1.0),)  
        lin = Linear(self.fr.offset-(1.0+eps)*self.fb.offset,self.fr.slope-(1.0+eps)*self.fb.slope)
        return min(np.float64(1.0),lin.inverse((1+eps)*self.f_rb(0.0)+eps*self.c))
    #todo
    def calc_t_N_minus_1(self,eps):
        #t_{N-1}\gets \min\{f_{r}^{-1}(c/ \epsilon) ,\min\{x:x\geq 1,\forall y\geq x,\frac{f_r(x)}{f_b(x)}\geq (1+\epsilon)\frac{f_r(y)}{f_b(y)} \}\} $
        val1 = self.fr.inverse(self.c/eps)
        
        maxRatio = self.rho_inf()
        if maxRatio==np.inf:
            val2=np.inf
        else:
            #fr = maxRatio/(1+eps) *fb
            lin = Linear(self.fr.offset-self.fb.offset*maxRatio/(1+eps) ,self.fr.slope-self.fb.slope*maxRatio/(1+eps) ) 
            val2 = lin.inverse(0.0)
        return min(val1,val2)


                   
#check
d1=Linear(0.0,1.0)     
d2=Linear(1.0,0.0)    
csr=InstanceLinear(d1,d2,1.0)
csr.approxOptStr(0.5)

