# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 19:11:00 2016

@author: Evyatar
"""


from math import *
import numpy as np
np.random.seed(1)

import scipy
from Instance import * 
from Linear import *
from PieceWiseLinear import *

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
#insert "%matplotlib tk"
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle


#TODO
#debug - check TODO minus+slope==0



class InstancePieceWiseLinear(Instance):
    
    def __init__(self,fr,fb,c):
        self.fr = fr
        self.fb = fb
        self.c=np.float64(c)
        self.times=fr.times
        self.norm_cost,self.norm_time=1.0,1.0
        if self.fr[1]==self.fb[1]==1:
            self.norm=True
        else:
            self.norm = False
    def __str__(self):
        s=""
        for t,lr,lb in zip(self.times,self.fr.lins,self.fb.lins):
            s=s+"Time>= "+str(t)+":  Rent: "+str(lr)+" Buy: "+str(lb)+"\n"
        return s
    def intersection(self):
        eps_num=1e-10
        for i,t in enumerate(self.times):
            int_t = (self.fb.lins[i].offset-self.fr.lins[i].offset)/(self.fr.lins[i].slope-self.fb.lins[i].slope)
            if i<(len(self.times)-1):            
                if int_t>=t and int_t<=self.times[i+1]:
                    return int_t
            
            elif int_t+eps_num>=t:
                    return int_t
                    
    def normalize(self):
        #we need to init new instance and return this        
        intersect=self.intersection()
        val_inter=self.fr[intersect]
        
        #print intersect,val_inter
        #all slopes need to be updated (*intersections/val_intersect)
        #all times need to be updated (/intersection)
        #all offsets need to be updated (/val_inter)
        
        normTimes =  [t/intersect for t in self.times]
        diagLinn = Linear(1.0,0.0)
        normFr = [ Linear(l.slope*intersect/val_inter,l.offset/val_inter)   for l in self.fr.lins]      
        normFb = [ Linear(l.slope*intersect/val_inter,l.offset/val_inter)   for l in self.fb.lins] 
        for ind in range(len(normTimes)):
            next_t = normTimes[ind+1] if ind<len(normTimes)-1 else np.inf
            if 1.>normTimes[ind] and 1.<=next_t:
                ind_intersect=ind+1

        
        if 1.0 not in self.times:
            normTimes.insert(ind_intersect,1.0)
            normFr.insert(ind_intersect,normFr[ind_intersect-1]) 
            normFb.insert(ind_intersect,normFb[ind_intersect-1]) 
        
    
        
        
        normed = InstancePieceWiseLinear( PieceWiseLinear(normFr,normTimes),PieceWiseLinear(normFb,normTimes),self.c/val_inter)
        normed.norm = True
        normed.norm_time=intersect
        normed.norm_cost = val_inter
        return normed
        
        
    def plot(self,num_fig=1) :
        t_max=self.times[-1]*2.0
        tt=self.times[:]
        tt.append(t_max)
        cost_r=[self.fr[t] for t in tt]
        cost_b=[self.fb[t] for t in tt]
        xx=[t for t in tt]   

        plt.figure(num_fig)
        plt.clf()
        plt.hold(True)
        plt.plot(xx,cost_r,'r')
        plt.plot(xx,cost_b,'b')
        plt.axis([0.0,xx[-1]*1.0,0.0,1.1*max(np.max(cost_r),np.max(cost_b))])
        plt.show()
        fs=20
        plt.xlabel("Time",fontsize=fs)
        plt.ylabel("Cost",fontsize=fs)
        #print cost_b,xx
        return plt
        
    def plot_q(self,eps,norm=True,num_fig=1):
        cr,q,T=self.approxOptStr(eps)
        Q=np.cumsum(q)
        if norm:
            xx=[t*self.norm_time for t in T]  
            Q=Q*self.norm_cost
        else:
            xx=T[:]
        plt.figure(num_fig)
        plt.hold(True)
        plt.title("Instance")
        plt.plot(xx,Q,'g')
        plt.legend(["Rent","buy","Probabilty"],loc=0)
        plt.xlabel("Time")   
        plt.ylabel("Cost (Probabilty)") 
        plt.text(xx[-2],self.norm_cost,"Probabilty 1.0")
        plt.text(xx[-2],0.5*self.norm_cost,"Probabilty 0.5")
        plt.show()
        
        print "plot_q",xx,Q
        
        
    def cr_z(self,z):
        if z==0:
            return self.fb.lins[0].offset/self.fr.lins[0].offset
        elif z==np.inf:
            if self.fr.lins[-1].slope==self.fb.lins[-1].slope:
                return self.fr[np.inf]/self.fb[np.inf]
            return self.fr.lins[-1].slope/self.fb.lins[-1].slope
        else:
            raise Exception("Error in InstnacePieceWiseLinear.cr_z")

    def cr_0(self):
        return div(self.fb.lins[0].offset,self.fr.lins[0].offset)
    def cr_inf(self):
        return div(self.fr.lins[-1].slope,self.fb.lins[-1].slope) if self.fr.lins[-1].slope>0 else cons.fr(self.times[-1]+10)/cons.fb(self.times[-1]+10)

    def zm(self):
        pass
    
    def cr_zm(self):
        return h(self.zm())



    #todo - g in linear g(z)=max fr/fb(y) (y>z)
    def g(self,z):
        pass
    #todo -h(z)= (fr(z)+c)/fb(z)
    def h(self,z):
        pass
    

    #we iterate on all the lines and find if there is inersection between h and g
    #if there is we find the exact intersection
    #otherwise we return inf    
    def intersection_time(self):      
        #print "inter"
        for ind,t in enumerate(self.times):
            if t>=1:
                #at point 1 g<h and then g always up
                #we need to divide to 2 parts:
                #1.h decrese or not up - then we need only to look on the start and end
                #2.h increase - then we need to solve 
                #TODO - try to checjk sopes or use ratio
                g1=self.fr[t]/self.fb[t]
                h1 = (self.fr[t]+self.c)/self.fb[t]
                if ind==len(self.times)-1:
                    next_t=np.inf
                    g2_n,g2_d=self.fr[next_t],self.fb[next_t]
                    if g2_d==np.inf:
                        g2=div(self.fr.lins[ind].slope,self.fb.lins[ind].slope)
                    else:
                        g2=div(g2_n,g2_d)
                    
                    h2_n,h2_d = self.fr[next_t]+self.c,self.fb[next_t]
                    if h2_d==np.inf:
                        h2=div(self.fr.lins[ind].slope,self.fb.lins[ind].slope)
                    else:
                        h2=div(h2_n,h2_d)
                else:
                    next_t=self.times[ind+1]
                    g2=self.fr[next_t]/self.fb[next_t]
                    h2 = (self.fr[next_t]+self.c)/self.fb[next_t]
                #print "inter",ind,g1,g2,h1,h2
                #TODO - find intersection better (maybe the intersection function will do this)
                #if max(g1,g2) >=min(h1,h2):
                if True:
                    #we have here intersection - what is the solve?
                    #if g1<=g2 then the function is monotone - simple solve
                    if g1<g2:
                        #we need to find zs:   
                        roots_intersection = ratio_intersection(self.fr.lins[ind],self.fb.lins[ind],self.fr.lins[ind]+self.c,self.fb.lins[ind])
                    if g1>=g2:
                    #then the function is const
                        roots_intersection = ratio_intersection(Linear(0,g1),Linear(0,1.0),self.fr.lins[ind]+self.c,self.fb.lins[ind])
                    inds = np.where( (roots_intersection>=t) & (roots_intersection<=next_t))[0]
                    if inds.size>0:
                        return roots_intersection[inds[0]]
                else:
                    pass
                   

        return float("inf")
        

    def middle_time_cr(self):
        z_s = self.intersection_time()
        #print "z_s at z_m",z_s
        #the competitve is always h(z) - we need to find the maximum
        cr_m = (self.fr[1]+self.c)/self.fb[1]
        z_m=1.0
        #because at each interval h always goes         
        for ind,t in enumerate(self.times):
            if t>=1 and t<=z_s:
                ht = (self.fr(t)+self.c)/self.fb(t)
                if ht<cr_m:
                    z_m=t
                    cr_m = ht
                #print "iter at z_m",ind,t,ht,cr_m
                    
        #at the end we need to check z_s or np.inf
        if t<z_s:
            t=z_s
            ht=self.cr_inf() if z_s==np.inf else (self.fr(t)+self.c)/self.fb(t)
            if ht<cr_m:
                z_m=t
                cr_m = ht
            #print "last iter at z_m",ht,cr_m,self.cr_inf() if z_s==np.inf else (self.fr(t)+self.c)/self.fb(t),z_s==np.inf,self.fr(t),self.c,self.fb(t),t
        #print "middle_time_cr",(z_m,cr_m)
        return (z_m,cr_m)
                
                
        """
                if ind==len(self.times)-1:
                    next_t=min(np.inf,z_s)
                    h2_n,h2_d = self.fr[next_t]+self.c,self.fb[next_t]
                    if h2_d==np.inf:
                        h2=div(self.fr.lins[ind].slope,self.fb.lins[ind].slope)
                    else:
                        h2=div(h2_n,h2_d)
                else:
                    next_t=min(self.times[ind+1],z_s)
                    h2 = (self.fr[next_t]+self.c)/self.fb[next_t]
                if cr_m > min(h1,h2):
                    cr_m = min(h1,h2)
                    if h1<=h2:
                        z_m=self.times[ind]
                    else:
                        z_m=next_t
        return (z_m,cr_m)
        
        """
        
    def middle_time(self):
        return self.middle_time_cr()[0]
        
    #todo 
    def between(self,t,ind):
        if ind==len(self.times)-1:
            if t>=self.times[ind]:
                return True
        elif (t>=self.times[ind]) and (t<=self.times[ind+1]):
            return True
        return False
    def calc_t1(self,eps):
        #t1
        #f_{r}(z)+c =(1+eps)(f_{rb}(0) +c+f_b(z) 
        #frb(t1)-eps*f_b(t1)=(1+eps)f_{rb}(0)+eps*c
        
        #t1 = min(np.float64(1.0),)
        
        for ind in range(len(self.times)):
            if self.times[ind]>=1.0:
                return 1.0
            else:
                lin = self.fr.lins[ind]-self.fb.lins[ind]*(1.0+eps)
                t1_candidate=lin.inverse((1+eps)*self.f_rb(0.0)+eps*self.c)
                if (t1_candidate>=self.times[ind]) and (t1_candidate<= self.times[ind+1]):
                    return min(np.float64(1.0),t1_candidate)
    #todo
    #fix situation when the ratio raise
    def calc_t_N_minus_1(self,eps):
        global lin
        #t_{N-1}\gets \min\{f_{r}^{-1}(c/ \epsilon) ,\min\{x:x\geq 1,\forall y\geq x,\frac{f_r(x)}{f_b(x)}\geq (1+\epsilon)\frac{f_r(y)}{f_b(y)} \}\} $

        if (self.c<=self.fr[1.0]*eps):
            return 1.0
        val1 = max(1.0,self.fr.inverse(self.c/eps))
        maxRatio = self.cr_inf()
        #print self.fr.lins[-1],val1,self.c/eps
        #print maxRatio/(1+eps),self.fr.inverse(self.c/eps)
        #print "val1",val1
        
        if maxRatio<np.inf:
            maxRatioLin = [maxLin(self.fr.lins[ind-1],self.fb.lins[ind-1],self.times[ind-1],self.times[ind])[0] for ind in range(1,len(self.times))]
            maxRatioLin.append(maxLin(self.fr.lins[-1],self.fb.lins[-1],self.times[-1],np.inf)[0])
            #print maxRatioLin,maxLin(self.fr.lins[-1],self.fb.lins[-1],self.times[-1],np.inf)
            #fr = maxRatio/(1+eps) *fb
            #print "max Ratios",maxRatioLin
            for ind in range(len(self.times)):
                if self.times[ind]<1.0:
                    continue
                #TODO - if ratio is down we need to consider from (ind+1) else ind
                maxRatio = max(maxRatioLin[ind:])
                #print ind,maxRatio,self.times[ind],eps
                
                lin_cond = self.fr.lins[ind]-self.fb.lins[ind]*(maxRatio/(1+eps) ) 
                val2_candidate = lin_cond.inverse(0.0) #this is the point where the ratio is fine
                if (val2_candidate<self.times[ind]) and (self.fr[self.times[ind]]/self.fb[self.times[ind]] >=maxRatio/(1+eps)):
                    #this case when the candidate is out of range (and we take the start of the range)                        
                    #print "oor",val2_candidate,self.times[ind],maxRatio,self.fr.lins[ind],self.fb.lins[ind]
                    val2_candidate = self.times[ind]
                    #print self.between(val2_candidate,ind)
                if self.between(val2_candidate,ind): #    (val2_candidate>=self.times[ind]) and (val2_candidate<= self.times[ind+1]) :
                    val2_cand_ratio = self.fr[val2_candidate]/self.fb[val2_candidate]
                    if  (val2_cand_ratio+1e-5>=(maxRatio/(1+eps))):
                        #print "return",ind,min(val1,val2_candidate),val2_cand_ratio,maxRatio/(1+eps),val1,val2_candidate
                        #print self.fr[50]/self.fb[50]
                        return min(val1,val2_candidate)
             
        return val1
        
    def valid(self):
        if min([self.times[1:][j]-self.times[:-1][j]for j in range(len(self.times)-1)]) < 0.:
            return False
        if min([min(self.fr.lins[i].slope,self.fb.lins[i].slope) for i in range(len(self.times))]) < 0:
            return False
        return True
                         
       


                   
#check
"""
pr=PieceWiseLinear([Linear(1.0/3,2.0/3.0),Linear(1.0/3,2.0/3.0),Linear(4.0/3,0)],[0.0,1.0,1.5])
pb=PieceWiseLinear([Linear(1.0,0.0),Linear(1.0,0.0),Linear(1.0,0.0)],[0.0,1.0,4.0/3])
exmp2=InstancePieceWiseLinear(pr,pb,1.0)



r0,b0=35.0,45.0
t1,t2=6.0,10.0
ar1,ab1=6.0,6.0
ending=True
if ending:
    endT=6.0
    pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1),Linear(r0+endT*ar1,0.0)],[0.0,t1,t2,t1+endT])
    pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1),Linear(b0-t2*ab1+ab1*(t1+endT),0.0)],[0.0,t1,t2,t1+endT])
else:
    pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1)],[0.0,t1,t2])
    pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1)],[0.0,t1,t2])
cons=InstancePieceWiseLinear(pr,pb,(b0-r0)*1.35)
cons.plot(False,2)
cons.normalize()
print cons.calc_t_N_minus_1(0.01)
cons.plot(False,3)
#plt.figure(3)
#plt.plot(np.arange(1.0,2.5,0.01),[cons.fr[i]/cons.fb[i] for i in np.arange(1.0,2.5,0.01)],'g-')
cons.plot_q(0.01,True,2)



%[b,a1,break,a2]
%HOT - 8G=40sh +1G=6sh
%12G - 50sh + 1G=6sh
%first we deal only what happen after the first 8G
%50s->1s , 1.5G-1t
% fr=[0.8,6*1.5/50,4/1.5,6*1.5/50];
% fb=[1,0,4/1.5,6*1.5/50];
% c=1.0;
%result : Never switch


% fr=[40,6,4,6];
% fb=[50,0,4,6];
% c=50;

%4.2 - we have no tight
% fr=[0.25,0.75,5/3,0.2];
% fb=[1,0,5/3,0.1];
% c=1;

%4.3 - we have unoptimal tight
%Intersection (t=8/9,cost=1)
fr=[1/3,0.75,4/3,0];
fb=[1,0,4/3,0];
c=1;





"""


#TODOS
#simulate instance of normalized linear ski rental
#check intersection
#check randomized competeive ratio
#check det cometitve ratio

def UT_linear():
    global i,xx,gg,hh
    print "Unit test linear ski rental problems"
    eps=0.1
    eps_num = 1e-10
    rr = np.arange(0.0,1.02,0.1)
    num_iter=100
    for eps in [0.1]:
        print "epsilon",eps
        for i in range(num_iter):
            a1 = max(np.random.choice(rr),0.05)
            a2=2
            while a2>=a1:
                a2 = np.random.choice(rr)
            #a1,a2=1.0,0.2
            b1,b2 = 1-a1,1-a2
            c = (b2-b1)*np.random.uniform(1,5)
            #a1,b1,a2,b2,c=0.5,0.5,0.0,1.0,0.684677189538
            
            print ",".join(str(u) for u in [a1,b1,a2,b2,c])
            pr=PieceWiseLinear([Linear(a1,b1),Linear(a1,b1)],[0,1])
            pb=PieceWiseLinear([Linear(a2,b2),Linear(a2,b2)],[0,1])
            cons=InstancePieceWiseLinear(pr,pb,c)
            #test 1
            
            
            xx=np.arange(1.0,20.0,0.01)
            gg= np.array([cons.fr[x]/cons.fb[x] for x in xx])
            
            gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
            hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
            gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
            hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
            #plt = cons.plot()
            #plt.hold(True)
            #plt.plot(xx,gg,"g")
            #plt.plot(xx,hh,"y")
            assert(np.abs(cons.intersection()-1.0)<eps_num)
            #print cons.calc_t1(eps)
            
            t1 = cons.calc_t1(eps)
            #print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
            assert (-eps_num <t1<1.0+eps_num)
            assert((cons.fr[t1]+c) <= eps_num+(1+eps)*(cons.f_rb(0)+c+cons.fb[t1]))
            assert((np.abs(t1-1.0)<eps_num) or np.abs(cons.fr[t1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t1]))<eps_num)
    
            alp = cons.calc_t_N_minus_1(eps)
            #print alp
            assert (alp>= 1.0)
            #print alp,cons.rho_inf(),cons.fr(alp)
            flag1 =(cons.fr[alp]/cons.fb[alp]) <= eps_num+ (1+eps)*cons.cr_inf()
            #print alp,c,cons.fr[alp]
            flag2 = (c<=eps_num+eps*cons.fr[alp])
            assert(flag1 or flag2)
            
            beta = (b2-b1)/c
            rho_prob_anylsis = np.exp(beta)/(np.exp(beta)-b2+b1)
            
            #print cons.approxOptStr(eps)[0],rho_prob_anylsis
            ratio = cons.approxOptStr(eps)[0]/rho_prob_anylsis
            assert ((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
    
            #check det
            rho_det_anylsis = min(div(b2,b1),c+1,div(a1,a2))
            t_det_anylsis = [0,1,np.inf][np.argmin([div(b2,b1),c+1,div(a1,a2)])]
            #print cons.intersection_time()
            #print cons.middle_time()
            #print cons.middle_time_cr()
            #print cons.cr_0(),cons.cr_inf()
            t_det_gen,rho_det_gen = cons.optDetStr()
            ratio = rho_det_gen/rho_det_anylsis
            print ratio,rho_det_gen,rho_det_anylsis
            assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
            
            if t_det_gen==np.inf and t_det_anylsis==np.inf:
                ratio=0.0
            else:
                ratio = t_det_gen-t_det_anylsis
                
            assert((ratio<=eps_num) and (ratio>= -eps_num))
    
            #we need to check cr_0,cr_inf
            
            
    
            #cons.plot(num_fig=i)
            #check norm - return new instance
            #check plot
            #print ""
            
            
            
    return True
 
 
def create_instance(br0,bb0,ar_s,ab_s,ts,c):
    lins_r,lins_b=[Linear(ar_s[0],br0)],[Linear(ab_s[0],bb0)]
    for ind,t in enumerate(ts[1:]):
        yr,yb=lins_r[-1][t],lins_b[-1][t]
        lins_r.append(Linear(ar_s[ind+1],-ar_s[ind+1]*t+yr))
        lins_b.append(Linear(ab_s[ind+1],-ab_s[ind+1]*t+yb))
    pr=PieceWiseLinear(lins_r,ts)
    pb=PieceWiseLinear(lins_b,ts)
    return InstancePieceWiseLinear(pr,pb,c)
     
 
 
def UT_PWL_5params():
    global cons
    print "Unit test piece-wise linear(5 parameters) ski rental problems"
    eps=0.05
    eps_num = 1e-10
    rr = np.arange(0.0,1.02,0.1)
    num_iter=100
    
    for i in range(num_iter):
        #print i
        ar0 = np.round(np.random.uniform(0.06,1),1)
        ab0 = np.round(np.random.uniform(0.0,ar0-0.07),1)
        ar1 = np.round(np.random.uniform(0.2,5.0),1)
        ab1=np.round(np.random.uniform(0.2,ar1-0.1),1)
        c=max(ar0-ab0,np.round(10**np.random.uniform(-0.3,1.0),1))
        

        #ar0,ab0,ar1,ab1,c =0.3,0.2,2.4,0.3,2.3
        #print ",".join([str(u) for u in [ar0,ab0,ar1,ab1,c]])
    
        cons = create_instance(1-ar0,1-ab0,[ar0,ar1],[ab0,ab1],[0.,1.],c)
        plt=cons.plot(num_fig=1)
        plt.hold(True)
        #print cons.optDetStr()
                
        xx=np.arange(1.0,20.0,0.1)
        gg= np.array([cons.fr[x]/cons.fb[x] for x in xx])
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        
        #plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        plt.plot(xx,gg,"g",xx,hh,"k")
        intr = cons.intersection()
        #print intr
        assert(np.abs(intr-1.0)<eps_num)
        
        t_1 = cons.calc_t1(eps)
        #print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num <t_1<1.0+eps_num)
        assert((cons.fr[t_1]+c) <= eps_num+(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))
        assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))<eps_num)


        alp = cons.calc_t_N_minus_1(eps)

        assert (alp>= 1.0)
        #print alp,cons.rho_inf(),cons.fr(alp)
        flag1 =(cons.fr[alp]/cons.fb[alp]) >= -eps_num+ cons.cr_inf()/(1+eps)
        #print alp,c,cons.fr[alp]
        flag2 = (c<=eps_num+eps*cons.fr[alp])
        #print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c
        
        assert(flag1 or flag2)
        
        

        cr0 = cons.cr_0()
        cr0_thr = div(1-ab0,1-ar0)
        assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf ))
        
        cri = cons.cr_inf()
        
        cri_thr = div(ar1,ab1) if ar1>0 else cons.fr[2]/cons.fb[2]
        assert(np.abs(cri_thr-cri)<eps_num or (min(cri,cri_thr)==np.inf ))
        #TODO - check det
        #we have cr_0 and cr_inf
         
        #fr = ar1(x-1)+1=ar1*x-ar1+1,fb=ab1(x-1)+1=ab1*x-ab1+1
        #(a1b2-a2b1)=ar1*(-ab1+1) - ab1*(-ar1+1)=ar1-ab1>0
        #g(z)= max(fr/fb,y<=z) = fr(z)/fb(z)
        #z_s = np.inf
        #min(h)
        #(a1b2-a2b1)=ar1*(-ab1+1) - ab1*(-ar1+1+c)=ar1-ab1-ab1*c
        hd =ar1-ab1-ab1*c
        cr_1= (cons.fr(1)+c)/cons.fb(1)
        if hd>=0:
            zm=1
            cr_zm=cr_1
        else:
            zm=None
            cr_zm = cri+eps #because h is decreasing and cr_inf <= h(z) z->inf
        
        
        
        rho_det_anylsis = min(cr0,cr_zm,cri)
        t_det_anylsis = [0,zm,np.inf][np.argmin(np.array([cr0,cr_zm,cri]))]
        t_det_gen,rho_det_gen = cons.optDetStr()
        ratio = rho_det_gen/rho_det_anylsis
        assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
        
        if t_det_gen==np.inf and t_det_anylsis==np.inf:
            ratio=0.0
        else:
            ratio = t_det_gen-t_det_anylsis
            
        assert((ratio<=eps_num) and (ratio>= -eps_num))  
        
        
        tup=cons.approxOptStr(eps)        
        ratio =tup[0]/rho_det_anylsis
        print "ratios 5params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio<= eps_num+(1+eps))
        #((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
    return True     
 
 
 
 
 
def UT_PWL_6params():
    global cons,gg,hh,xx
    print "Unit test piece-wise linear(6 parameters) ski rental problems"
    eps=0.05
    eps_num = 1e-10
    rr = np.arange(0.0,1.02,0.1)
    num_iter=100
    
    for i in range(num_iter):
        #print i
        ar0 = np.round(np.random.uniform(0.06,1),1)
        ab0 = np.round(np.random.uniform(0.0,ar0-0.07),1)
        ar1 = np.round(np.random.uniform(0.2,5.0),1)
        ab1=np.round(np.random.uniform(0.2,ar1-0.1),1)
        c=max(ar0-ab0,np.round(10**np.random.uniform(-0.3,1.0),1))
        t1=np.round(np.random.uniform(1.2,5.0),1)

        #ar0,ab0,ar1,ab1,c,t1 = 0.5,0.3,0.2,0.2,0.8,1.6
        #print ",".join([str(u) for u in [ar0,ab0,ar1,ab1,c,t1]])
    
        cons = create_instance(1-ar0,1-ab0,[ar0,ar0,ar1],[ab0,ab0,ab1],[0.,1.,t1],c)
        cons=cons.normalize()
        plt=cons.plot()
        plt.hold(True)
    
        assert (np.max(np.abs(np.array([cons.fr.lins[i](cons.times[i+1])-cons.fr.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps)
        assert (np.max(np.abs(np.array([cons.fb.lins[i](cons.times[i+1])-cons.fb.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps)
               
        xx=np.arange(1.0,50.0,0.01)
        gg= np.array([cons.fr[x]/cons.fb[x] for x in xx])
        
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        #plt.plot(xx,gg,"g",xx,hh,"k")
        #plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        plt.plot(xx,gg,"g")
        plt.plot(xx,hh,"y")
        
        assert(np.abs(cons.intersection()-1.0)<eps_num)
        
        t_1 = cons.calc_t1(eps)
        #print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num <t_1<1.0+eps_num)
        assert((cons.fr[t_1]+c) <= eps_num+(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))
        assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))<eps_num)


        alp = cons.calc_t_N_minus_1(eps)

        assert (alp>= 1.0)
        #print alp,cons.rho_inf(),cons.fr(alp)
        flag1 =(cons.fr[alp]/cons.fb[alp]) >= -eps_num+ cons.cr_inf()/(1+eps)
        #print alp,c,cons.fr[alp]
        flag2 = (c<=eps_num+eps*cons.fr[alp])
        #print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c
        
        assert(flag1 or flag2)
        
        

        cr0 = cons.cr_0()
        cr0_thr = div(1-ab0,1-ar0)
        assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf ))
        
        cri = cons.cr_inf()
        
        cri_thr = div(ar1,ab1) if ar1>0 else cons.fr[t1+10]/cons.fb[t1+10]
        #print cri_thr,cri
        assert(np.abs(cri_thr-cri)<eps_num or (min(cri,cri_thr)==np.inf ))
        #TODO - check det
        #we have cr_0 and cr_inf
         
        #fr = ar0(x-1)+1=ar0*x-ar0+1,fb=ab0(x-1)+1=ab1*x-ab1+1 
        #(a1b2-a2b1)=ar0*(-ab0+1) - ab0*(-ar0+1)=ar0-ab0>0
        #g(z)= max(fr/fb,y<=z) = fr(z)/fb(z)
        #1<t<t1 h>g   
         
        #fr = ar1(x-t1)+ar0(t1-1)+1=ar1*x+(ar0-ar1)*t1+1-ar0,fb=ab1*x+(ab0-ab1)*t1+1-ab0 
        #(a1b2-a2b1)=ar1*((ab0-ab1)*t1+1-ab0) -ab1*((ar0-ar1)*t1+1-ar0)=ar1*((ab0)*t1+1-ab0) -ab1*((ar0)*t1+1-ar0)
        #g(z)= max(fr/fb,y<=z) = fr(z)/fb(z)
        #t1<t           
        #if gd(t1)<0 =>g(>t1)=g(t1), else zs=inf 
         
        #first we will check the slopes         
        
        #we have 2 sections
        #1<t<t1 g=fr/fb if gd down then g=g(1) else g=g
        #h = (fr+c)/fb
        """
        h=lambda t:((cons.fr(t)+c)/cons.fb(t))
        hd1 = (cons.fr.lins[1].slope*cons.fb.lins[1].offset)-(cons.fb.lins[1].slope*(cons.fr.lins[1].offset+c))
        
        gd1 = (cons.fr.lins[1].slope*cons.fb.lins[1].offset)-(cons.fb.lins[1].slope*cons.fr.lins[1].offset)
        if gd1<=0:
            g=lambda t:(cons.fr(1)/cons.fb(1))
        else:
            h=lambda t:(cons.fr(t)/cons.fb(t))
            
        if g(t1)>=max(h(1),h(t1)):
            pass
            #we have intersection here
        else:
            #the same for the next interval
         
         
         
        gd = ar1*((ab0)*t1+1-ab0) -ab1*((ar0)*t1+1-ar0)
        if gd>0:
            pass
            #does h monotone inreases ? 
            #hd = rd*b-(r+c)*bd=
            
        
        
        
        #z_s = np.inf
        #min(h)
        #(a1b2-a2b1)=ar1*(-ab1+1) - ab1*(-ar1+1+c)=ar1-ab1-ab1*c
        hd =ar1-ab1-ab1*c
        cr_1= (cons.fr(1)+c)/cons.fb(1)
        if hd>=0:
            zm=1
            cr_zm=cr_1
        else:
            zm=None
            cr_zm = cri+eps #because h is decreasing and cr_inf <= h(z) z->inf
        """
        
        #we will find cr_zm numeric
        zs_numer =  xx[gg>=hh][0] if np.sum(gg>=hh)>0 else np.inf
        zm_numer = xx[np.argmin(hh[xx<=zs_numer])] if ((cri>np.min(hh[xx<=zs_numer])or zs_numer<np.inf)) else np.inf
        #cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm,cr_zm = cons.middle_time_cr()
        print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
        rho_det_anylsis = min(cr0,cr_zm,cri)
        t_det_anylsis = [0,zm,np.inf][np.argmin(np.array([cr0,cr_zm,cri]))]
        #print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen,rho_det_gen = cons.optDetStr()
        ratio = rho_det_gen/rho_det_anylsis
        assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
        
        if t_det_gen==np.inf and t_det_anylsis==np.inf:
            ratio=0.0
        else:
            ratio = t_det_gen-t_det_anylsis
            
        assert((ratio<=eps_num) and (ratio>= -eps_num))        
        #cons.approxOptStr(eps)
        tup=cons.approxOptStr(eps)        
        ratio =tup[0]/rho_det_anylsis
        #print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        #assert (ratio<= eps_num+(1+eps))
        #print "no assert of probalistic startegy at 6params"
    
    return True     
 
 
 
  
 
def exmp_det():
    global cons,gg,hh,xx,xx_e,cr_xx_e,cons_org,c,td,br0,bb0,ar1,dict_res
    print "Examples of instances for learning about determistic strategy"
    eps=0.05
    eps_num = 1e-10
    dict_res={}
    bb0=50.0
    ar1 = 6.0
    genr = ((br0,c,td) for br0 in np.arange(10.0,48.0,2) for c in np.arange((bb0-br0)+0.2,bb0-0.1,2) for td in np.arange((bb0-br0)/ar1+0.1,10.0,0.5))
    for i,(br0,c,td) in enumerate(genr):
        """
        #hot
        br0=40
        bb0=50
        ts=[0,12,18]
        ar_s = [0,6,0]
        ab_s = [0,0,0]
        c = 10
        """
        """
        #P_mobile
        br0=35
        bb0=40
        ts=[0,15,20,25,30]
        ar_s = [0,5,5,0,0]
        ab_s = [0,0,5,5,0]
        c =40
        """

        assert ((td-eps_num) > ((bb0-br0)/ar1)),"Bad TD "+",".join(str(u) for u in [c,br0,td])
        #print ",".join(str(si) for si in [c,br0,td])
        assert (c>(bb0-br0-eps_num)) ,"Bad C "+",".join(str(u) for u in [c,br0,td])       
        
        
        
        ts=[0,12,12+td]
        ar_s = [0,ar1,0]
        ab_s = [0,0,0]        
        
        cons_org = create_instance(br0,bb0,ar_s,ab_s,ts,c)
        cons=cons_org.normalize()
        #plt=cons.plot()
        #plt.hold(True)

        assert (np.max(np.abs(np.array([cons.fr.lins[i](cons.times[i+1])-cons.fr.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps)
        assert (np.max(np.abs(np.array([cons.fb.lins[i](cons.times[i+1])-cons.fb.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps)
               
        xx=np.arange(1.0,20.0,0.01)
        gg= np.array([cons.fr[x]/cons.fb[x] for x in xx])
        
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        #plt.plot(xx,gg,"g",xx,hh,"k")
        #plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        #plt.plot(xx,gg,"g")
        #plt.plot(xx,hh,"y")
        
        assert(np.abs(cons.intersection()-1.0)<eps_num)
        
        t_1 = cons.calc_t1(eps)
        #print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num <t_1<1.0+eps_num)
        assert((cons.fr[t_1]+c) <= eps_num+(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))
        #assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))<eps_num)
        #print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp>= 1.0)
        #print alp,cons.rho_inf(),cons.fr(alp)
        flag1 =(cons.fr[alp]/cons.fb[alp]) >= -eps_num+ cons.cr_inf()/(1+eps)
        #print alp,c,cons.fr[alp]
        flag2 = (c<=eps_num+eps*cons.fr[alp])
        #print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c
        
        #assert(flag1 or flag2)
        
        

        cr0 = cons.cr_0()
        cr0_thr = div(bb0,br0)
        #assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf ))
        #print "no assert of cr0"
        cri = cons.cr_inf()
        
        cri_thr = div(ar_s[-1],ab_s[-1]) if ar_s[-1]>0 else cons.fr[20]/cons.fb[20]
        #print cri_thr,cri
        assert(np.abs(cri_thr-cri)<eps_num or (min(cri,cri_thr)==np.inf ))
        
        zs_numer =  xx[gg>=hh][0] if np.sum(gg>=hh)>0 else np.inf
        zm_numer = xx[np.argmin(hh[xx<=zs_numer])] if ((cri>np.min(hh[xx<=zs_numer])or zs_numer<np.inf)) else np.inf
        #cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm,cr_zm = cons.middle_time_cr()
        #print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
        rho_det_anylsis = min(cr0,cr_zm,cri)
        t_det_anylsis = [0,zm,np.inf][np.argmin(np.array([cr0,cr_zm,cri]))]
        #print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen,rho_det_gen = cons.optDetStr()
        print "det str",br0,c,td,rho_det_gen,t_det_gen
        
        #we will build based on XX the competitve ratio
        xx_e = np.arange(0.0,20.0,0.05)
        xx_e[-1]=np.inf
        cr_xx_e =[ np.max([cons.cr(z,y) for y in xx_e]) for z in xx_e]
        #print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
        diff = np.min(cr_xx_e)-np.max([cons.cr(t_det_gen,y) for y in xx_e])
        assert((diff<=eps_num) and (diff>= -eps_num))  
        ratio = rho_det_gen/rho_det_anylsis
        assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
        
        if t_det_gen==np.inf and t_det_anylsis==np.inf:
            ratio=0.0
        else:
            ratio = t_det_gen-t_det_anylsis
            
        assert((ratio<=eps_num) and (ratio>= -eps_num))        
        #cons.approxOptStr(eps)
        tup=cons.approxOptStr(eps)        
        ratio =tup[0]/rho_det_anylsis
        #print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio<= eps_num+(1+eps)),"prob less good than det"
        #print "no assert of probalistic startegy at 6params"
        markersize=20   
        if False:
            if t_det_gen==0:
                plt.plot(c,td,"b*",markersize=markersize)
            elif t_det_gen<np.inf:
                plt.plot(c,td,"gx",markersize=markersize)
            elif t_det_gen==np.inf:
                plt.plot(c,td,"ro",markersize=markersize)
            else:
                print "no det choice"
        """
        simArtists = [plt.Line2D((0,1),(0,0), color=c, marker=mark, linestyle='-') for c,mark in zip(["b","g","r"],["*","x","*"])]
        plt.legend(simArtists,["Strategy starts at buy","Strategy switch at intersection","Strategy never switch"])
        plt.ylabel("Slope of the rent option",fontsize=20)
        plt.xlabel("Switch cost",fontsize=20)
        plt.title("Best strategy as function of switch cost and slope",fontsize=30)
        """
        dict_res[c,td,br0]=(rho_det_gen,t_det_gen)
        if i%100==-1%100:
            print "saved",i
            pickle.dump(dict_res,open("dict_res_det_c_td_br0.pkl","w+"))   
    
    
    pickle.dump(dict_res,open("dict_res_det_c_td_br0.pkl","w+"))   
    
    return True     
 


def exmp_prob():
    global cons,gg,hh,xx,cons_org,c,td,br0,bb0,ar1,dict_res,C,f,A_eq,b_ub,b_eq,q,T,t_det_anylsis,rho_det_anylsis,t_det_gen,rho_det_gen
    print "Examples of instances for learning about probabilties strategy"
    eps=0.05
    eps_num = 1e-10
    num_iter=10
    dict_res={}
    bb0=50.0
    ar1 = 6.0
    genr = ((br0,c,td) for br0 in np.arange(10.0,48.0,2) for c in np.arange((bb0-br0)+0.2,bb0-0.1,2) for td in np.arange((bb0-br0)/ar1+0.1,10.0,0.5))
    for i,(br0,c,td) in enumerate(genr):
        #br0,c,td=10.0,40.2,6.76666666667
        assert (td-eps_num) > ((bb0-br0)/ar1),"BAD TD"+",".join(str(si) for si in [br0,c,td])

        assert (c>(bb0-br0-eps_num)),"BAD C"+",".join(str(si) for si in [br0,c,td])      
        
        
        ts=[0,12,12+td]
        ar_s = [0,ar1,0]
        ab_s = [0,0,0]        
        
        cons_org = create_instance(br0,bb0,ar_s,ab_s,ts,c)
        cons=cons_org.normalize()
        #plt=cons.plot()
        #plt.hold(True)
    
        assert (np.max(np.abs(np.array([cons.fr.lins[i](cons.times[i+1])-cons.fr.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps),"BAD fr"+",".join(str(si) for si in [br0,c,td])
        assert (np.max(np.abs(np.array([cons.fb.lins[i](cons.times[i+1])-cons.fb.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps),"BAD fb"+",".join(str(si) for si in [br0,c,td])
               
        xx=np.arange(1.0,20.0,0.01)
        gg= np.array([cons.fr[x]/cons.fb[x] for x in xx])
        
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        #plt.plot(xx,gg,"g",xx,hh,"k")
        #plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        #plt.plot(xx,gg,"g")
        #plt.plot(xx,hh,"y")
        
        assert(np.abs(cons.intersection()-1.0)<eps_num),"BAD intersection"+",".join(str(si) for si in [br0,c,td])
        
        t_1 = cons.calc_t1(eps)
        #print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num <t_1<1.0+eps_num),"BAD t1_1 "+",".join(str(si) for si in [br0,c,td])
        assert((cons.fr[t_1]+cons.c) <= eps_num+(1+eps)*(cons.f_rb(0)+cons.c+cons.fb[t_1])),"BAD t1_2 "+",".join(str(si) for si in [br0,c,td])
        tol=0.5
        assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+cons.c-(1+eps)*(cons.f_rb(0)+cons.c+cons.fb[t_1]))<eps_num),"BAD t1_3 "+",".join(str(si) for si in [br0,c,td])+"|" +str(t_1)+"|"+str(np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1])))+"|"+str(cons.fr[1.0]+c)+"|"+str((1+eps)*(cons.f_rb(0)+c+cons.fb[1.0]))
        #print "no assert of t_1"

        alp = cons.calc_t_N_minus_1(eps)

        assert (alp>= 1.0),"BAD alp"+",".join(str(si) for si in [br0,c,td])
        #print alp,cons.rho_inf(),cons.fr(alp)
        flag1 =(cons.fr[alp]/cons.fb[alp]) >= -eps_num+ cons.cr_inf()/(1+eps)
        #print alp,c,cons.fr[alp]
        flag2 = (cons.c<=eps_num+eps*cons.fr[alp])
        #print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c
        
        assert(flag1 or flag2),"BAD flags"+",".join(str(si) for si in [br0,c,td])
        
        

        cr0 = cons.cr_0()
        cr0_thr = div(bb0,br0)
        assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf )),"BAD cr0"+",".join(str(si) for si in [br0,c,td])
        
        cri = cons.cr_inf()
        
        cri_thr = div(ar_s[-1],ab_s[-1]) if ar_s[-1]>0 else cons.fr[20]/cons.fb[20]
        #print cri_thr,cri
        assert(np.abs(cri_thr-cri)<eps_num or (min(cri,cri_thr)==np.inf )),"BAD cri"+",".join(str(si) for si in [br0,c,td])
        
        zs_numer =  xx[gg>=hh][0] if np.sum(gg>=hh)>0 else np.inf
        zm_numer = xx[np.argmin(hh[xx<=zs_numer])] if ((cri>np.min(hh[xx<=zs_numer])or zs_numer<np.inf)) else np.inf
        #cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
        zs = cons.intersection_time()
        zm,cr_zm = cons.middle_time_cr()
        #print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
        rho_det_anylsis = min(cr0,cr_zm,cri)
        t_det_anylsis = [0,zm,np.inf][np.argmin(np.array([cr0,cr_zm,cri]))]
        #print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
        t_det_gen,rho_det_gen = cons.optDetStr()
        print "prob str",c,td,br0
        
        #we will build based on XX the competitve ratio
        xx_e = np.arange(0.0,20.0,0.05)
        xx_e[-1]=np.inf
        cr_xx_e =[ np.max([cons.cr(z,y) for y in xx_e]) for z in xx_e]
        #print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
        diff = np.min(cr_xx_e)-np.max([cons.cr(t_det_gen,y) for y in xx_e])
        eps_diff=0.02        
        assert (np.abs(diff)<=eps_diff)   ,"bad opt det"+",".join(str(si) for si in [br0,c,td])+"|diff= "+str(diff)
        ratio = rho_det_gen/rho_det_anylsis
        assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps))) ,"bad opt det"+",".join(str(si) for si in [br0,c,td])
        
        if t_det_gen==np.inf and t_det_anylsis==np.inf:
            ratio=0.0
        else:
            ratio = t_det_gen-t_det_anylsis
            
        assert((ratio<=eps_num) and (ratio>= -eps_num))  ,"bad opt det ratio"+",".join(str(si) for si in [br0,c,td])       
        #cons.approxOptStr(eps)
        tup=cons.approxOptStr(eps)  
        q=tup[1]
        T=tup[2]
        C,f,A_eq,b_ub,b_eq=tup[4]
        ratio =tup[0]/rho_det_anylsis
        #print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
        assert (ratio<= eps_num+(1+eps)),"prob less good then det|ratio={}|det_cr={}|prob_cr={}".format(ratio,rho_det_anylsis,tup[0])+"|"+",".join(str(si) for si in [br0,c,td])
        #TODO - make some random probs for seeing if the LP works good
        markersize=20   
        if False:
            if t_det_gen==0:
                plt.plot(c,td,"b*",markersize=markersize)
            elif t_det_gen<np.inf:
                plt.plot(c,td,"gx",markersize=markersize)
            elif t_det_gen==np.inf:
                plt.plot(c,td,"ro",markersize=markersize)
            else:
                print "no det choice"
        #dict_res[c,td,ar1,br0,bb0]
        dict_res[c,td,br0]=(tup[0],q[0],q[-1],q,len(T))
        if i==-1:
            print "saved",i
            pickle.dump(dict_res,open("dict_res_prob_c_td_ar1_br0.pkl","w+"))   
    
    
    
    return True     


def exmp_p6():
    #TODO
    #add switch cost to text
    #add middle time to text (which is optimal)
    global cons,gg,hh,xx,cons_org,c,td,br0,bb0,ar1,dict_res,cr_xx_e,xx_e
    print "Example of instance for determinstic strategy "
    eps=0.05
    eps_num = 1e-10
    num_iter=1
    
    #$a_r^1,a_b^1,a_r^2,a_b^2,a_r^3 = 0.5,0,2,1,1$

    #$t_1,t_2,c=2 ,5 ,1$}

    br0,bb0 = 0.5,1
    ts=[0,2,5]
    ar_s = [0.5,2,2]
    ab_s = [0,1,1] 
    #if ar_s[-1] = 2 then optimal is 2.0 but z_m=inf else z_m =6.15 but it is not optimal
    c=1       
    
    cons_org = create_instance(br0,bb0,ar_s,ab_s,ts,c)
    cons=cons_org.normalize()
    plt=cons.plot()
    plt.hold(True)

    assert (np.max(np.abs(np.array([cons.fr.lins[i](cons.times[i+1])-cons.fr.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps),"BAD fr"+",".join(str(si) for si in [br0,c,td])
    assert (np.max(np.abs(np.array([cons.fb.lins[i](cons.times[i+1])-cons.fb.lins[i+1](cons.times[i+1]) for i in range(len(cons.times)-1)])))<eps),"BAD fb"+",".join(str(si) for si in [br0,c,td])
           
    xx=np.arange(1.0,20.0,0.01)
    gg= np.array([cons.fr[x]/cons.fb[x] for x in xx])
    
    gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
    hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
    #plt.plot(xx,gg,"g",xx,hh,"k")
    #plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
    gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
    hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
    plt.plot(xx,gg,"g")
    plt.plot(xx,hh,"y")
    
    plt.legend(["Rent","Buy","g(z)","h(z)"],loc=0)
    plt.text(0.5,1.33,"Switch cost:  1")
    plt.text(1.5,2.7,"middle time:  2")
    plt.plot(2,2.5,"y*")
    assert(np.abs(cons.intersection()-1.0)<eps_num),"BAD intersection"+",".join(str(si) for si in [br0,c,td])
    
    t_1 = cons.calc_t1(eps)
    #print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
    assert (-eps_num <t_1<1.0+eps_num),"BAD t1_1 "+",".join(str(si) for si in [br0,c,td])
    assert((cons.fr[t_1]+cons.c) <= eps_num+(1+eps)*(cons.f_rb(0)+cons.c+cons.fb[t_1])),"BAD t1_2 "+",".join(str(si) for si in [br0,c,td])
    tol=0.5
    assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+cons.c-(1+eps)*(cons.f_rb(0)+cons.c+cons.fb[t_1]))<eps_num),"BAD t1_3 "+",".join(str(si) for si in [br0,c,td])+"|" +str(t_1)+"|"+str(np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1])))+"|"+str(cons.fr[1.0]+c)+"|"+str((1+eps)*(cons.f_rb(0)+c+cons.fb[1.0]))
    #print "no assert of t_1"

    alp = cons.calc_t_N_minus_1(eps)

    assert (alp>= 1.0),"BAD alp"+",".join(str(si) for si in [br0,c,td])
    #print alp,cons.rho_inf(),cons.fr(alp)
    flag1 =(cons.fr[alp]/cons.fb[alp]) >= -eps_num+ cons.cr_inf()/(1+eps)
    #print alp,c,cons.fr[alp]
    flag2 = (cons.c<=eps_num+eps*cons.fr[alp])
    #print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c
    
    assert(flag1 or flag2),"BAD flags"+",".join(str(si) for si in [br0,c,td])
    
    

    cr0 = cons.cr_0()
    cr0_thr = div(bb0,br0)
    assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf )),"BAD cr0"+",".join(str(si) for si in [br0,c,td])
    
    cri = cons.cr_inf()
    
    cri_thr = div(ar_s[-1],ab_s[-1]) if ar_s[-1]>0 else cons.fr[20]/cons.fb[20]
    #print cri_thr,cri
    assert(np.abs(cri_thr-cri)<eps_num or (min(cri,cri_thr)==np.inf )),"BAD cri"+",".join(str(si) for si in [br0,c,td])
    
    zs_numer =  xx[gg>=hh][0] if np.sum(gg>=hh)>0 else np.inf
    g_zs_numer = gg[gg>=hh][0] if np.sum(gg>=hh)>0 else np.inf
    zm_numer = xx[np.argmin(hh[xx<=zs_numer])] if ((cri>np.min(hh[xx<=zs_numer])or zs_numer<np.inf)) else np.inf
    #cr_zm_numer = np.max(hh[xx<=zs_numer]) if ((cri>np.max(hh)or zs_numer<np.inf)) else cri
    zs = cons.intersection_time()
    zm,cr_zm = cons.middle_time_cr()
    print "zs",zs,zs_numer,zm,zm_numer,cri,cons.norm_time
    plt.plot(zs,g_zs_numer,"*")
    print zs,g_zs_numer
    rho_det_anylsis = min(cr0,cr_zm,cri)
    t_det_anylsis = [0,zm,np.inf][np.argmin(np.array([cr0,cr_zm,cri]))]
    #print "det inline",rho_det_anylsis,t_det_anylsis,cr_zm,cr0,cri
    t_det_gen,rho_det_gen = cons.optDetStr()
    print "det str",rho_det_gen,t_det_gen
    
    #we will build based on XX the competitve ratio
    xx_e = np.arange(0.0,20.0,0.05)
    xx_e[-1]=np.inf
    cr_xx_e =[ np.max([cons.cr(z,y) for y in xx_e]) for z in xx_e]
    #print "best cr and time",c,np.min(cr_xx_e),xx[np.argmin(cr_xx_e)],np.max([cons.cr(t_det_gen,y) for y in xx_e])
    diff = np.min(cr_xx_e)-np.max([cons.cr(t_det_gen,y) for y in xx_e])
    diff = 0 if np.min(cr_xx_e)==np.inf and     np.max([cons.cr(t_det_gen,y) for y in xx_e])==np.inf else diff
    assert((diff<=eps_num) and (diff>= -eps_num))  ,str(diff)+"|"+str(np.min(cr_xx_e))+"|"+str(np.max([cons.cr(t_det_gen,y) for y in xx_e]))
    ratio = rho_det_gen/rho_det_anylsis
    assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps))) ,"bad opt det"+",".join(str(si) for si in [br0,c,td])
    
    if t_det_gen==np.inf and t_det_anylsis==np.inf:
        ratio=0.0
    else:
        ratio = t_det_gen-t_det_anylsis
        
    assert((ratio<=eps_num) and (ratio>= -eps_num))  ,"bad opt det ratio"+",".join(str(si) for si in [br0,c,td])       
    #cons.approxOptStr(eps)
    tup=cons.approxOptStr(eps)        
    ratio =tup[0]/rho_det_anylsis
    #print "ratios 6params",ratio,rho_det_anylsis,tup[0],t_det_anylsis
    #assert (ratio<= eps_num+(1+eps)),"prob less good then det"
    print "assert  PROB!!!"

    
    
    return True     
































#ar1,ab1 - control the parameters before 1
#ar2,ab2 - control the parametres after 1 unti t1
#ar3 - the slope of rent after t2 - can be 0 or ab2
#t1 - the time when rent decrease his slope to ab2
#t2 - the time when buy decreses his slope to 0
#c - swith cost

def UT_PWL_8params():
    global cons
    print "Unit test linear ski rental problems"
    eps=0.05
    eps_num = 1e-10
    num_iter=100
    
    for i in range(num_iter):
        print i
        ar1 = np.round(np.random.uniform(0.06,1),1)
        ab1 = np.round(np.random.uniform(0.0,ar1-0.07),1)
        ar2 = np.round(np.random.uniform(0.2,5.0),1)
        ab2=np.round(np.random.uniform(0.2,ar2-0.1),1)
        c=min(ar1-ab1,np.round(10**np.random.uniform(-0.3,1.0),1))
        ar3=np.random.choice([0.0,ab2])
        t1=np.round(np.random.uniform(1.2,5.0),1)
        #TODO - make sure t2 is before intersection
        t2 = np.round(np.random.uniform(t1+0.1,5.0),1)
        ab3 = ab2 #0.0
        #ar1,ab1,ar2,ab2,ar3,t1,t2,c = 0.5, 0.3 ,0.2 ,0.2 ,0.0, 2.7 ,3.7 ,0.2
        #print ar1,ab1,ar2,ab2,ar3,t1,t2,c
    
        pr=PieceWiseLinear([Linear(ar1,1-ar1),Linear(ar2,1-ar2),Linear(ab2,-ab2*t1+ar2*t1+1-ar2),Linear(ar3,-ar3*t2+ab2*t2-ab2*t1+ar2*t1+1-ar2)],[0.0,1.,t1,t2])
        pb=PieceWiseLinear([Linear(ab1,1-ab1),Linear(ab2,1-ab2),Linear(ab2,1-ab2),Linear(ab3,-ab3*t2+ab2*t2+1-ab2)],[0.0,1.,t1,t2])
        cons=InstancePieceWiseLinear(pr,pb,c)
        
        #cons.plot(num_fig=i)
        #print cons.optDetStr()
                
        xx=np.arange(1.0,20.0,0.1)
        gg= np.array([cons.fr[x]/cons.fb[x] for x in xx])
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        #plt.plot(xx,gg,"g",xx,hh,"k")
        #plt.plot(xx,np.array([cons.fr[x]/cons.fb[x] for x in xx]),"g")
        gg= np.array([np.max(gg[:j+1]) for j in range(gg.size)])
        hh =  np.array([(cons.fr[x]+c)/cons.fb[x] for x in xx])
        
        
        assert(np.abs(cons.intersection()-1.0)<eps_num)
        
        t_1 = cons.calc_t1(eps)
        #print t1,(cons.fr[t1]+c) -(1+eps)*(cons.f_rb(0)+c+cons.fb[t1])
        assert (-eps_num <t_1<1.0+eps_num)
        assert((cons.fr[t_1]+c) <= eps_num+(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))
        assert((np.abs(t_1-1.0)<eps_num) or np.abs(cons.fr[t_1]+c-(1+eps)*(cons.f_rb(0)+c+cons.fb[t_1]))<eps_num)


        alp = cons.calc_t_N_minus_1(eps)

        assert (alp>= 1.0)
        #print alp,cons.rho_inf(),cons.fr(alp)
        flag1 =(cons.fr[alp]/cons.fb[alp]) >= eps_num+ cons.cr_inf()/(1+eps)
        #print alp,c,cons.fr[alp]
        flag2 = (c<=eps_num+eps*cons.fr[alp])
        #print flag1,flag2,alp,cons.fr[alp],cons.fb[alp],cons.cr_inf(),cons.c
        
        #assert(flag1 or flag2)
        
        #beta = (b2-b1)/c
        #rho_prob_anylsis = np.exp(beta)/(np.exp(beta)-b2+b1)
        
        #print cons.approxOptStr(eps)[0],rho_prob_anylsis
        #ratio = cons.approxOptStr(eps)[0]/rho_prob_anylsis
        #assert ((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))

        #check det
        #rho_det_anylsis = min(div(b2,b1),c+1,div(a1,a2))
        #t_det_anylsis = [0,1,np.inf][np.argmin([div(b2,b1),c+1,div(a1,a2)])]
   
        #t_det_gen,rho_det_gen = cons.optDetStr()
        #ratio = rho_det_gen/rho_det_anylsis
        #assert((ratio<=eps_num+(1+eps)) and (ratio>= -eps_num+1/(1+eps)))
        
        #if t_det_gen==np.inf and t_det_anylsis==np.inf:
        #    ratio=0.0
        #else:
        #    ratio = t_det_gen-t_det_anylsis
            
        #assert((ratio<=eps_num) and (ratio>= -eps_num))        
        
        
        cr0 = cons.cr_0()
        cr0_thr = div(1-ab1,1-ar1)
        assert(np.abs(cr0_thr-cr0)<eps_num or (min(cr0,cr0_thr)==np.inf ))
        
        cri = cons.cr_inf()
        
        cri_thr = np.inf if ar3>0 else cons.fr[50]/cons.fb[50]
        #assert(np.abs(cri_thr-cri)<eps_num or (min(cri,cri_thr)==np.inf ))
        #TODO - check det
        
        cons.approxOptStr(eps)
    
    return True    

#we may found some example a=5,b=1,t1=3,t2=5,c=2.5 but the optDetStr is a liitle bit strange
#we should debug this instance and add it to the UT
"""
    a,b,t1,c = 2,1,2,1.
    t2=5.
    ar = 0.5
    ar1,ab1,ar2,ab2,ar3 = 0.5,0.,2,1,1
    t1,t2,c=2,5,1.
"""







#TODO
#find examples - linear and cellular
#selection on examples
#$a_r^1,a_b^1,a_r^2,a_b^2,a_r^3 = 0.5,0,2,1,1$
#t_1,t_2,c=2 ,5 ,1$}

def exmp1_cellular():
    #first example of determinstic solution
    #HOT - 8G=40sh +1G=6sh
    #12G - 50sh + 1G=6sh
    num_exmp = 10
    i=0
    while i < num_exmp:
        r0,b0=np.random.uniform(35.0,40.0),np.random.uniform(45.0,70.)
        t1,t2=np.random.uniform(2.0,8.0),np.random.uniform(9.0,14.)
        ar1,ab1=6.0,6.0
        factor_c = np.random.uniform(1.0,2.0)
        ending=True
        if ending:
            endT=6.0
            pr=PieceWiseLinear([Linear(0.0,r0),Linear(ar1,r0-t1*ar1),Linear(ar1,r0-t1*ar1),Linear(0.0,r0+endT*ar1)],[0.0,t1,t2,t1+endT])
            pb=PieceWiseLinear([Linear(0.0,b0),Linear(0.0,b0),Linear(ab1,b0-t2*ab1),Linear(0.0,b0-t2*ab1+ab1*(t1+endT))],[0.0,t1,t2,t1+endT])
        else:
            pr=PieceWiseLinear([Linear(0.0,r0),Linear(ar1,r0-t1*ar1),Linear(ar1,r0-t1*ar1)],[0.0,t1,t2])
            pb=PieceWiseLinear([Linear(0.0,b0),Linear(0.0,b0),Linear(ab1,b0-t2*ab1)],[0.0,t1,t2])        
        
        cons=InstancePieceWiseLinear(pr,pb,(b0-r0)*factor_c)
        
        intr=cons.intersection()
        if intr==np.inf:
            continue
        if not cons.valid():
            continue
        
        cons.plot(num_fig=i+1)
        
        
        #cons_norm = cons.normalize()
        #cons_norm.plot(num_fig=2)
        
        z_opt,_=cons.optDetStr()
        z_plot = np.minimum(z_opt,plt.xticks()[0][-2])
        z_plot=np.maximum(z_plot,1)
        plt.plot(i+1);plt.plot([z_plot,z_plot],[0,plt.yticks()[0][-1]],"-k")
        if z_opt >intr and z_opt<np.inf:
            print i,True
        else:
            pass
            print False,intr,z_opt,z_plot
        i=i+1
        #print [str(lin) for lin in cons.fb.lins]
        #print cons.times
        #print cons_norm.optDetStr()
    #print cons_norm.middle_time(),cons_norm.intersection_time(),cons_norm.optDetStr(),cons_norm.c
    #TODO - plot the results of the DET on the graph
       
def exmp2_cellular():
    #second  example of DET - here we will show intersing properities of 
    #instances and their strategy - consult Boaz
    #HOT - 8G=40sh +1G=6sh
    #12G - 50sh + 1G=6sh
    #try only 2 slopes it will suppose to work
    r0,b0=35.0,45.0
    t1,t2=6.0,10.0
    ar1,ab1=6.0,6.0
    factor_c = 1.0
    ending=True
    if ending:
        endT=6.0
        pr=PieceWiseLinear([Linear(0.0,r0),Linear(ar1,r0-t1*ar1),Linear(ar1,r0-t1*ar1),Linear(0.0,r0+endT*ar1)],[0.0,t1,t2,t1+endT])
        pb=PieceWiseLinear([Linear(0.0,b0),Linear(0.0,b0),Linear(ab1,b0-t2*ab1),Linear(0.0,b0-t2*ab1+ab1*(t1+endT))],[0.0,t1,t2,t1+endT])
    else:
        pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1)],[0.0,t1,t2])
        pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1)],[0.0,t1,t2])
    cons=InstancePieceWiseLinear(pr,pb,(b0-r0)*factor_c)
    cons.plot(num_fig=1)

    cons_norm = cons.normalize()
    cons_norm.plot(num_fig=2)
           
    
def exmp3_cellular():
    #third  example of TIGHT (4.1) - we need to find example of instance without tight startegy
    #better with f_r unbounded (ending is false)
    
    #CODE function which find TIGHT and reports when they could not find one
    #FIND what is the optimal PROB
    #ANYLSIS why we could never find such strategy - HARD
    #we might use some points as milestone and show the unextinence of tight during bounds on CR
    #we can also show imposiblity for smooth startegy
    
    
    r0,b0=3 
    t1,t2=6.0,10.0
    ar1,ab1=6.0,6.0
    factor_c = 1.0
    ending=False
    if ending:
        endT=6.0
        pr=PieceWiseLinear([Linear(0.0,r0),Linear(ar1,r0-t1*ar1),Linear(ar1,r0-t1*ar1),Linear(0.0,r0+endT*ar1)],[0.0,t1,t2,t1+endT])
        pb=PieceWiseLinear([Linear(0.0,b0),Linear(0.0,b0),Linear(ab1,b0-t2*ab1),Linear(0.0,b0-t2*ab1+ab1*(t1+endT))],[0.0,t1,t2,t1+endT])
    else:
        #pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1)],[0.0,t1,t2])
        #pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1)],[0.0,t1,t2])
        
        pr=PieceWiseLinear([Linear(0.0,r0),Linear(ar1,r0-t1*ar1),Linear(ar1,r0-t1*ar1)],[0.0,t1,t2])
        pb=PieceWiseLinear([Linear(0.0,b0),Linear(0.0,b0),Linear(ab1,b0-t2*ab1)],[0.0,t1,t2])        
        
    cons=InstancePieceWiseLinear(pr,pb,(b0-r0)*factor_c)
    cons.plot(num_fig=1)

    cons_norm = cons.normalize()
    cons_norm.plot(num_fig=2)    
    
def exmp4_cellular():
    global cons,cr,cr_lp
    #forth  example of TIGHT (4.2) - show instace which is bound and has tight strategy
    #which are not optimal (HARD - maybe impossible)
    #we can maybe approximate the optimal PROB with EXP curve and then show tight is not optimally
    #Q - how could we be sure this is th eonly TIGHT ?
    eps=0.01
    print "exmp4 - tight which not optimal (bound case)"
    eps_num = 1e-10
    num_iter=1
    
    for i in range(num_iter):
        print i
        ar1 = np.round(np.random.uniform(0.06,1),1)
        ab1 = np.round(np.random.uniform(0.0,ar1-0.07),1)
        ar2 = np.round(np.random.uniform(0.3,5.0),1)
        ab2=np.round(np.random.uniform(0.1,ar2-0.1),1)
        c=max(ar1-ab1,np.round(10**np.random.uniform(-0.3,1.0),1))
        ar3=0.0 #ab2 #np.random.choice([0.0,ab2])
        t1=np.round(np.random.uniform(1.2,5.0),1)
        #TODO - make sure t2 is before intersection
        t2 = np.round(np.random.uniform(t1+0.1,5.0),1)
        ab3 =0.0# ab2 #0.0
        ar1,ab1,ar2,ab2,ar3,ab3,t1,t2,c = 0.9, 0.0 ,0.9 ,0.0 ,0.0, 0.0,1.5 ,2.0 ,2.0
        ar1,ab1,ar2,ab2,ar3,ab3,ar4,t1,t2,c = 0.9, 0.0 ,0.9 ,0.0 ,0.0, 0.0,5.0,1.1 ,3.0 ,1.0
        #print ar1,ab1,ar2,ab2,ar3,ab3,t1,t2,c
    
        pr=PieceWiseLinear([Linear(ar1,1-ar1),Linear(ar2,1-ar2),Linear(ar3,-ar3*t1+ar2*t1+1-ar2),Linear(ar4,-ar4*t2+ar3*t2-ar3*t1+ar2*t1+1-ar2)],[0.0,1.,t1,t2])
        pb=PieceWiseLinear([Linear(ab1,1-ab1),Linear(ab2,1-ab2),Linear(ab2,1-ab2),Linear(ab3,-ab3*t2+ab2*t2+1-ab2)],[0.0,1.,t1,t2])
        
        #example
        ar1,ab1,ar2,ab2,ar3,ab3,t1,t2 = 0.9,0.0,4.0,0.0,0.0,0.0,0.4,0.55
        #m(x-x0)+y0=n => n = -mx0+y0
        pr=PieceWiseLinear([Linear(ar1,1-ar1),Linear(ar2,-ar2*t1+ar1*t1+1-ar1),Linear(ar3,-ar3*t2+ar2*t2-ar2*t1+ar1*t1+1-ar1)],[0.0,t1,t2])
        pb=PieceWiseLinear([Linear(ab1,1-ab1),Linear(ab2,-ab2*t1+ab1*t1+1-ab1),Linear(ab3,-ab3*t2+ab2*t2-ab2*t1+ab1*t1+1-ab1)],[0.0,t1,t2])
        
        
        cons=InstancePieceWiseLinear(pr,pb,c)
        cons = create_instance(0.1,1.0,[0.9,0.9,0.0,0.3,0],[0.0,0.0,0.0,0.2],[0.0,1.0,1.3,2.2],1.0)
        cons.plot(num_fig=i)
        cons_norm = cons.normalize()
        cons_norm.plot(num_fig=i)
        cr,q,T,cr_lp = cons_norm.approxOptStr(eps)
        
        plt.plot(T,cr+cr_lp,"k*-")
        print len(T),np.cumsum(q).shape
        plt.plot(T[:-1],np.cumsum(q)[:-1],"g*-")
        print "CR",cr
        #print T
    return True    

  
  
  
def exmp5_cellular():
    #fifth  example of appendix - show parameters of instance for APPROX
    r0,b0=35.0,45.0
    t1,t2=6.0,10.0
    ar1,ab1=6.0,6.0
    factor_c = 1.0
    ending=False
    if ending:
        endT=6.0
        pr=PieceWiseLinear([Linear(0.0,r0),Linear(ar1,r0-t1*ar1),Linear(ar1,r0-t1*ar1),Linear(0.0,r0+endT*ar1)],[0.0,t1,t2,t1+endT])
        pb=PieceWiseLinear([Linear(0.0,b0),Linear(0.0,b0),Linear(ab1,b0-t2*ab1),Linear(0.0,b0-t2*ab1+ab1*(t1+endT))],[0.0,t1,t2,t1+endT])
    else:
        #pr=PieceWiseLinear([Linear(r0,0.0),Linear(r0-t1*ar1,ar1),Linear(r0-t1*ar1,ar1)],[0.0,t1,t2])
        #pb=PieceWiseLinear([Linear(b0,0.0),Linear(b0,0.0),Linear(b0-t2*ab1,ab1)],[0.0,t1,t2])
        
        pr=PieceWiseLinear([Linear(0.0,r0),Linear(ar1,r0-t1*ar1),Linear(ar1,r0-t1*ar1)],[0.0,t1,t2])
        pb=PieceWiseLinear([Linear(0.0,b0),Linear(0.0,b0),Linear(ab1,b0-t2*ab1)],[0.0,t1,t2])        
        
    cons=InstancePieceWiseLinear(pr,pb,(b0-r0)*factor_c)
    cons.plot(num_fig=1)

    cons_norm = cons.normalize()
    cons_norm.plot(num_fig=2)     
flag,flag2=None,None
#flag = UT_linear()
if not flag:
    print "Linear tests  has not been done"
    
#flag2 = UT_PWL_5params()
if not flag2:
    print "5 params tests  has not been done"

flag=flag and flag2 
if flag:
    print "Good Test"
print ""



#UT_PWL_6params()
#exmp_det()
#exmp_prob()
#exmp_p6()
#exmp1_cellular()

#exmp2_cellular()

#exmp3_cellular()

#exmp4_cellular()

#exmp5_cellular()

import pickle
def plot_results_det():
    global dict_res
    file_name="dict_res_det_c_td_br0.pkl"
    dict_res=pickle.load(open(file_name,"r"))
    markersize=18   
    fontsize=20
    
    
    dict_marker = {np.inf:"o",1.0:"d",0.0:"*"}
    fig = plt.figure(2)
    fig.clf()
    ax = Axes3D(fig)
    sc=ax.scatter(xs=[k[0] for k in dict_res.keys() if dict_res[k][1]==np.inf],ys=[k[1] for k in dict_res.keys()if dict_res[k][1]==np.inf], zs=[k[2] for k in dict_res.keys() if dict_res[k][1]==np.inf],marker="o",c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1]==np.inf],cmap='plasma',s=60)
    sc=ax.scatter(xs=[k[0] for k in dict_res.keys() if dict_res[k][1]==1.0],ys=[k[1] for k in dict_res.keys()if dict_res[k][1]==1.0], zs=[k[2] for k in dict_res.keys() if dict_res[k][1]==1.0],marker="d",c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1]==1.0],cmap='plasma',s=60)
    sc=ax.scatter(xs=[k[0] for k in dict_res.keys() if dict_res[k][1]==0.0],ys=[k[1] for k in dict_res.keys()if dict_res[k][1]==0.0], zs=[k[2] for k in dict_res.keys() if dict_res[k][1]==0.0],marker="*",c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1]==0.0],cmap='plasma',s=60)
    
    
    
    ax.set_xlabel("Switch cost",fontsize=fontsize)
    ax.set_ylabel("Time diffrenece",fontsize=fontsize)
    ax.set_zlabel("Initial rent price",fontsize=fontsize)
    ax.set_title("Determinstic strategies properties",fontsize=25)
    simArtists = [plt.Line2D((4,2,40),(4,2.001,40), color="k", marker=m, linestyle='-',markersize=15) for m in ["*","d","o"]]  
    
    ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
    #ax.legend(sc)   
    fig.colorbar(sc)
    return
    for k in dict_res.keys():
        (c,td) = k
        rho_det_gen,t_det_gen = dict_res[k]
        if t_det_gen==0:
            plt.plot(c,td,"b*",markersize=markersize)
        elif t_det_gen<np.inf:
            plt.plot(c,td,"gx",markersize=markersize)
        elif t_det_gen==np.inf:
            plt.plot(c,td,"ro",markersize=markersize)
        else:
            print "no det choice"
    plt.figure(2)
    plt.clf()
    plt.hold(True)
    plt.xlabel("Switch cost",fontsize=fontsize)
    plt.ylabel("Time diffrenece",fontsize=fontsize)
    colormap = plt.cm.gist_ncar
    #you can add cmap
    
    sc=plt.scatter([k[0] for k in dict_res.keys()],[k[1] for k in dict_res.keys()],s=100,c=[(dict_res[k][0]) for k in dict_res.keys()],cmap='plasma')
    plt.colorbar()
    
    simArtists = [ax.Line3D((4,2,40),(4,2.001,40), color=c, marker=m, linestyle='-') for m in ["*","d","0"]]  
    
    #ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
    ax.legend(sc)   
    
    
    
    plt.show()
    #3D - https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
    #Axes3D.scatter(xs, ys, zs=0, zdir='z', s=20, c=None, depthshade=True, *args, **kwargs)
    
    return


def plot_results_prob():
    global dict_res
    file_name="dict_res_prob_c_td_ar1_br0.pkl"
    dict_res=pickle.load(open(file_name,"r"))
    markersize=18   
    fontsize=20
    
    
    ind_feature = 1 #0-cr,1-q0,2-q_inf
    strs = ["Competetive ratio","Start at buy probability","Never switch probability"]
    for ind_feature in range(3):
        fig = plt.figure(3+ind_feature)
        fig.clf()
        ax = Axes3D(fig)
        sc=ax.scatter([k[0] for k in dict_res.keys() ],[k[1] for k in dict_res.keys()], zs=[k[3] for k in dict_res.keys() ],marker="o",c=[(dict_res[k][ind_feature]) for k in dict_res.keys() ],cmap='plasma',s=60)
    
        
        
        ax.set_xlabel("Switch cost",fontsize=fontsize)
        ax.set_ylabel("Time diffrenece",fontsize=fontsize)
        ax.set_zlabel("Initial rent price",fontsize=fontsize)
        ax.set_title("Probabilities strategies properties ("+strs[ind_feature]+")",fontsize=30)
        #simArtists = [plt.Line2D((4,2,40),(4,2.001,40), color="k", marker=m, linestyle='-',markersize=15) for m in ["*","d","o"]]  
        
        #ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
        #ax.legend(sc)   
        fig.colorbar(sc)
    return
    


def plot_results_det_2D():
    global dict_res_slice
    #first we need to slice
    file_name="dict_res_det_c_td_br0.pkl"
    dict_res=pickle.load(open(file_name,"r"))
    print dict_res.keys()[:2]
    target_c = np.arange(6.2,49,8)
    target_td = np.arange(1.1,10,1)
    target_br = np.arange(10,47,6)
    for ind,target in enumerate(target_br):
        dict_res_slice = {(k[0],k[1]):dict_res[k] for k in dict_res.keys() if abs(k[2]-target)<0.005}
        markersize=18   
        fontsize=20
        
        
        dict_marker = {np.inf:"o",1.0:"d",0.0:"*"}
        plt.figure(10+ind)
        plt.clf()
        for km in dict_marker.keys():
            plt.scatter(x=[k[0] for k in dict_res_slice.keys() if dict_res_slice[k][1]==km],y=[k[1] for k in dict_res_slice.keys()if dict_res_slice[k][1]==km],marker=dict_marker[km],c=[(dict_res_slice[k][0]) for k in dict_res_slice.keys() if dict_res_slice[k][1]==km],cmap='plasma',s=60)
        #plt.scatter(xs=[k[1] for k in dict_res.keys() if dict_res[k][1]==1.0],ys=[k[2] for k in dict_res.keys()if dict_res[k][1]==1.0],marker="d",c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1]==1.0],cmap='plasma',s=60)
        #plt.scatter(xs=[k[1] for k in dict_res.keys() if dict_res[k][1]==0.0],ys=[k[2] for k in dict_res.keys()if dict_res[k][1]==0.0], marker="*",c=[(dict_res[k][0]) for k in dict_res.keys() if dict_res[k][1]==0.0],cmap='plasma',s=60)
        plt.show()
    return
    
    
    ax.set_xlabel("Switch cost",fontsize=fontsize)
    ax.set_ylabel("Time diffrenece",fontsize=fontsize)
    ax.set_zlabel("Initial rent price",fontsize=fontsize)
    ax.set_title("Determinstic strategies properties",fontsize=25)
    simArtists = [plt.Line2D((4,2,40),(4,2.001,40), color="k", marker=m, linestyle='-',markersize=15) for m in ["*","d","o"]]  
    
    ax.legend(simArtists,["best startegy start with buy","best startegy switch at middle time","best strategy never buy"])
    #ax.legend(sc)   
    fig.colorbar(sc)
    return

exmp_det()
exmp_prob()
#plot_results_det_2D()    
#plot_results_prob()

