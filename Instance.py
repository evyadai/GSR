# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:11:00 2016

@author: Evyatar
"""

import time
from math import *
import numpy as np
import scipy
import scipy.optimize

def div(a,b):
    if b==0 and a>0:
        return np.inf
    if a==0 and b==0:
        return np.nan
    return a/b



class Instance(object):
    def __init__(self):
        self.fr = fr
        self.fb = fb
        self.c=np.float64(c)
        self.norm =False
    def f_rb(self,x):
        #print "f_rb",self.fr[x]-self.fb[x],self.fr[x],self.fb[x]
        #t=input("")
        return self.fr[x]-self.fb[x]
    def opt(self,y):
        return min(self.fr[y],self.fb[y])
    
    def intersection(self):
        def diffrence(x):
            return self.fr(x)-self.fb(x)
        sol = optimize.root(diffrence, [0.0], method='hybr')
        return float(sol.x)

    def cost(self,z,y):
        if z==0:
            return self.fb(y)
        if y<z:
            return self.fr(y)
        if z==np.inf and y==np.inf:
            return self.fr(y)
        return (self.f_rb(z)+self.c+self.fb(y))
        
            
    def cr(self,z,y):
        if y==0 and self.opt(y)==0:
            if z>0:
                return np.float64(1.0)
            else:
                return np.inf
        if y==np.inf and self.fb[y]==np.inf:
            if z<np.inf:
                return 1.0
            else:
                return self.cr_inf() 
        return div(self.cost(z,y),self.opt(y))
            
    """   
    def cr_z(self,z):
        print "Unaccurate method due to local minima points"
        def crZ(y):
            return -self.cr(z,y)
        sol = optimize.minimize(crZ,[2.0],bounds=((0.0,None),) )
        #print(sol)
        return -float(sol.fun)
    def bfDet(self):
        pass        
    def cr_0(self):
        def neg_ratio(x):
            return -(self.fb(x)/self.fr(x))
        sol = optimize.minimize(neg_ratio,[0.0],bounds=((0,1),) )   
        return -float(sol.fun)
    def cr_inf(self):
        def neg_ratio(x):
            return -(self.fr(x)/self.fb(x))
        sol = optimize.minimize(neg_ratio,[0.0],bounds=((1.0,None),) )   
        #print(sol)        
        return -float(sol.fun)
    
        
    def g(self,z):
        def neg_ratio(x):
            return -(self.fr(x)/self.fb(x))
        sol = optimize.minimize(neg_ratio,[1.0],bounds=((1.0,z),) )   
        return -float(sol.fun)
    """
    def h(self,z):
        return div(self.fr[z]+self.c,self.fb[z])
    """
    def intersection_time(self):
        def diff(z):
            return self.g(z)-self.h(z)
        sol = optimize.root(diff,[1.0])
        if not sol.success:
            return float("inf")
        return float(sol.x)

    def middle_time(self):
        sol = optimize.minimize(self.h,[1.0],bounds=((1.0,self.intersection_time()),) )          
        return float(sol.x)
    """

    #TODO -det
    #cr_0
    #cr_inf
    #cr_z_m
    #cr_z?
    def cr_0(self):
        pass
    def cr_inf(self):
        pass
    def cr_zm(self):
        pass
    def cr_inf_overall(self):
        pass
    
    
    
    def optDetStr(self):
        if self.norm:
            normed=self
        else:
            normed = self.normalize()
        cr0 = normed.cr_0()
        crInf = normed.cr_inf_overall()
        zm,cr_zm = normed.middle_time_cr()
        if ((cr0 <=cr_zm) and (cr0 <= crInf)):
            res = (0,cr0)
        elif (cr_zm<=crInf):
            res =  (zm,cr_zm)
        else:
            res =(float("inf"),crInf)
        if self.norm:
            return res
        else:
            return (res[0]*normed.norm_time,res[1])

    
        
    def buildT(self,eps):
        T=[]
        T.append(0.0)
        ti=self.calc_t1(eps)
        
        while ti<1:
            T.append(ti)
            ti_new=min(np.float64(1.0),self.fr.inverse((1+eps)*self.fr[ti]),self.fb.inverse((1+eps)*self.fb[ti]))
            if ti==ti_new:
                raise Exception("Stuck at buildT.t_i="+str(ti))
            ti=ti_new
        
        #t_{N-1}
        t_N_minus_1=self.calc_t_N_minus_1(eps)
        #iter to t_{N-1}
        while ti<t_N_minus_1:
            #print ti
            T.append(ti)
            tr,tb=self.fr.inverse((1+eps)*self.fr[ti]),self.fb.inverse((1+eps)*self.fb[ti])
            tr,tb=tr if tr>-1  else np.inf,tb if tb>-1  else np.inf
            if ti>=tr:
                ti_new=tb
            elif ti>=tb:
                ti_new=tr
            else:
                ti_new=min(tr,tb)
            if ti>=ti_new:
                raise Exception("Stuck at buildT.t_i="+str(ti),tr,tb,(1+eps)*self.fb[ti])
            ti=ti_new
        T.append(t_N_minus_1)
        T.append(np.inf)
        return T
        
    def approxOptStr(self,eps):
        global C,f
        #delta  = eps**(2.0/3.0)
        T=self.buildT(eps)
        
        T_tag = T[:]
        if self.opt(0.0)==0.0:
            T_tag.remove(0.0)
        if self.cr_inf()==np.inf:
            T_tag.remove(np.inf)
        C=np.zeros((len(T),len(T_tag)+1))
        for ind_y,y in enumerate(T):
            for ind_z,z in enumerate(T_tag):
                if y==np.inf and self.opt(np.inf)==np.inf:
                    if z==np.inf :
                        C[ind_y,ind_z]=self.cr_inf()
                    else:
                        C[ind_y,ind_z]=np.float64(1.0)
                else:
                    C[ind_y,ind_z]=self.cr(z,y)
        
        C[:,-1]=-1.0
        f=np.zeros((len(T_tag)+1,))
        f[-1]=1.0
        A_eq=np.ones((1,len(T_tag)+1))
        A_eq[0,-1]=0.0
        b_ub=np.zeros_like(C[:,-1])
        #the problem are C*q-rho <= 0 
        b_eq=np.ones((1,))
        res = scipy.optimize.linprog(c=f,A_ub=C,b_ub=b_ub, A_eq=A_eq,b_eq=b_eq,options={"maxiter":10000}  )
        #tight startegy is the problem of C*q = rho'
        #con:  A_ub*x < b_ub  => constrain that rho is the maximal competitve ratio
        #C*q-rho' <0
        
        #con_eq: A_eq*x == b_eq => constainat that the sum of q is 1
        #the bound is q>=0,rho>0
        
        #target: min c^T *x => rho is minimal
        
        if (len(T_tag)!=len(T)):
            #print "tight - special with zeros"
            cr_tight=1.
        else:
            pass
            q_tight = np.dot( np.linalg.inv(C[:,:-1]),np.ones((C.shape[1]-1,)))
            
            cr_tight = 1/np.sum(q_tight)
            #print "cr_tight",cr_tight,np.sum(q_tight)
            #print np.linalg.inv(C[:,:-1])
            #print C[:,:-1],C[:,:-1].shape,len(T_tag),len(T)
            q_tight=q_tight*cr_tight
            if ((np.min(q_tight) <0.0) or (np.max(q_tight))):
                pass
               #print "tight is invalid",np.min(q_tight),np.max(q_tight)
           
        
        if res.status>0:
            print res.fun,res.status,res.nit
            raise Exception(res.message)
        q=list(res.x)[:-1]
        if self.opt(0.0)==0.0:
            q.insert(0,0.0)
        if self.cr_inf()==np.inf:
            q.append(0.0) 
        cr=res.fun
        #todo find tightness
        cr_lp = np.dot(C,res.x)
        #print cr,np.min(cr_lp),np.max(cr_lp)    ,res.x[-1]   
        #print cr_lp
        if (np.max(np.abs(cr_lp[1:])) > 1e-8)  :
            pass
            #print "untight (not 0)",np.max(np.abs(cr_lp[1:])),cr,self.cr_0()
        elif (np.max(np.abs(cr_lp)) > 1e-8)  :
            pass
            #print "untight (0)",self.cr_0(),np.max(np.abs(cr_lp[1:])),cr
        if (cr < cr_tight):
            pass
            #print "strong untight",cr,cr_tight,cr_tight-cr
        else:
            pass
            #print np.max(C),np.sum(np.sum(C[:,:-1],axis=0)),np.sum(np.sum(C[:,:-1],axis=1))
            #x_tight = np.hstack((q_tight,cr_tight))
            #print x_tight.shape
            #cr_lp_tight = np.dot(C,x_tight)
            #print np.max(cr_lp_tight),np.min(cr_lp_tight),np.sum(q_tight)
            #print np.min(res.x)
        #now we need to find tight strategy (using ODE - use the algorithm in tight)        
        
        debug_tup = C,f,A_eq,b_ub,b_eq
        return (cr,q,T,cr_lp,debug_tup)
        
        
        
        
        
        
        
        

