# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 22:26:08 2016

@author: Evyatar
"""
from Linear import *
import numpy as np
#this class is build from list of lins which start from point at time
class PieceWiseLinear(object):
    
    def __init__(self,lins,times):
        self.lins=lins
        self.times=times
        
    def __getitem__(self,t):
        if t>=self.times[-1]:
            return self.lins[-1][t]
        for ind in range(len(self.times)):
            if t>=self.times[ind] and t<self.times[ind+1]:
                return self.lins[ind][t]
        print "no item"
    def __call__(self,t):
        return self.__getitem__(t)
    def __str__(self):
        s=""
        for t,l in zip(self.times,self.lins):
            s+= "time:{0} - ".format(t)+str(l)+"\n"
        return s
    def __len__(self):
        return len(self.lins)
        
    def value(self,t):
        return self.__getitem__(t)
    def inverse(self,x):
        if x==np.inf:
            raise Exception("Error in PieceWiseLinear.inverse()")
        if x<self.lins[0].offset:
            return np.float64(-1.0)
        for ind in range(len(self.times)):
            inv_i = self.lins[ind].inverse(x)
            next_t = self.times[ind+1] if ind<len(self.times)-1 else np.inf
            #print "inv",ind,inv_i,next_t
            if inv_i>-1.0 and inv_i <= next_t:
                #print "return",inv_i
                return inv_i
        return np.inf
        
    

    
