# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 22:26:08 2016

@author: Evyatar
"""

import numpy as np

def div(a,b):
    if b==0 and a>0:
        return np.inf
    if a==0 and b==0:
        return np.nan
    return a/b

class Linear(object):
    def __init__(self,slope,offset):
        self.offset=np.float64(offset)
        self.slope=np.float64(slope)
    def set_params(self,newOffset,newSlope):
        self.offset=np.float64(newOffset)
        self.slope=np.float64(newSlope)
    def __getitem__(self,t):
        if isinstance(t,slice):
            if t.step is None:
                return Linear(self.slope,self.offset)
        if t==np.inf and self.slope==0.0:
            return self.offset
        return np.float64(self.slope*t+self.offset)
    def __call__(self,t):
        return self.__getitem__(t)
        
    def value(self,t):
        return self.__getitem__(t)
    
    
    
    def inverse(self,x,flag="Min"):
        if x==np.inf:
            raise Exception("Error in Linear.inverse()")
        
        if x<self.offset:
            print "TODO minus"
            return np.float64(-1.0)
        elif self.slope==0 and self.offset<=x:
            return np.inf
        elif self.offset==x:
            return np.float64(0.0)  
        else:
            return (x-self.offset)/self.slope
    def __str__(self):
        return "(Slope={0},Offset={1})".format(self.slope,self.offset)
    def __add__(self,rhsLin):
        if isinstance(rhsLin,Linear):
            return Linear(self.slope+rhsLin.slope,self.offset+rhsLin.offset)
        else:
             return Linear(self.slope,self.offset+rhsLin)
    def __sub__(self,rhsLin):
        if isinstance(rhsLin,Linear):
             return Linear(self.slope-rhsLin.slope,self.offset-rhsLin.offset)
        else:
             return Linear(self.slope,self.offset-rhsLin)
    def __mul__(self,factor):
        return Linear(self.slope*factor,self.offset*factor)
    def __div__(self,factor):
        return Linear(self.slope/factor,self.offset/factor)



#return (maximum,maxima)
def maxLinOrg(a1,b1,a2,b2,A,B):
    #return max((a1*t+b1)/(a2*t+b2) for t\in [A,B]
    #derive = (a1b2-a2b1)/(a2*x+b2)^2
    if (a1*b2-a2*b1)==0.0:
        if a1==0.0 and a2==0.0:
            return (b1/b2,A)
        else:
            return (div(a1,a2),A)
    elif (a1*b2-a2*b1)>0.0:
        if B==np.inf:
            return (div(a1,a2),B)
        return ((a1*B+b1)/(a2*B+b2),B)
    else:
        return ((a1*A+b1)/(a2*A+b2),A)


def maxLin(a1,b1,a2,b2,A,B):
    return maxLinOrg(a1,b1,a2,b2,A,B)
def maxLin(l1,l2,A,B):
    return maxLinOrg(l1.slope,l1.offset,l2.slope,l2.offset,A,B)     



#TODO
def ratio_intersection(l1,l2,l3,l4):
    #print "ratio_intersection"
    a1,b1,a2,b2,a3,b3,a4,b4 = l1.slope,l1.offset,l2.slope,l2.offset,l3.slope,l3.offset,l4.slope,l4.offset
    A = a1*a4 - a2*a3
    B = a1*b4 + a4*b1 - a2*b3 - a3*b2
    C = b1*b4-b2*b3
    if A==0:
        if B==0:
            if C==0:
                print "Warn: A,B,C are zeros"
            return np.array([])
        return np.array([-C/B])
    delta = (B**2)-4*A*C
    if delta<0:
        print "No sulotion to ratio_intersection",A,B,C
        return np.array([])
    else:
        sol1 = (-B-np.sqrt(delta))/(2*A)
        sol2 = (-B+np.sqrt(delta))/(2*A)
        if ((sol1>=sup) and (sol1<=inf)):
            return np.array([sol1,sol2])
        if ((sol2>=sup) and (sol2<=inf)):
            return np.array([sol1,sol2])
        print "Sulotions to ratio_intersection not in range",A,B,C,sup,inf,sol1,sol2,delta,np.sqrt(delta)
        return np.array([sol1,sol2]) 
        
    
    





#check
d=Linear(0,1)