# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 14:32:18 2016

@author: Evyatar
"""
import numpy as np
from math import *
import numpy.linalg
from scipy import optimize
#q_0 = (r-1)/(bb/br-1)

#Q^i = k^i *exp(y*exp(ar-ab)/c) - opt'/(ar-ab)*r+ar/(ar-ab)


#Instance 4.2
#Intersection (t=1,cost=1)
#t<=5/3
#fr = 0.25+0.75t
#fb = 1
#t>=5/3
#fr=1.5+0.2(t-5/3)
#fb=1+0.1(t-5/3)

#VARS : q0,r,k0,k1,    *1*

"""
#q0 = (r-1)/3
eq1 = [1,-1.0/3.0,0,0,-1.0/3.0]
#Q0 (0<t<1)= k0*exp(0.75t)- r+1,Q(0)=q0
eq2 = [-1,-1,1,0,-1]

#Q1 (1<t<5/3)= k1*exp(0.75t)+1,Q0(1)=Q1(1)
eq3 = [0,-1,exp(0.75),-exp(0.75),0]

#Q2(5/3<t) = -r+2 , Q1(5/3) = Q2(5/3)
eq4 = [0.0,1.0,0.0,exp(1.25),1.0]

u=np.array([eq1,eq2,eq3,eq4])
import numpy.linalg

x1=numpy.linalg.solve(u[:,:4],u[:,4])
"""
#q_0 = (r-1)/(bb/br-1)

#Q^i = k^i *exp(y*exp(ar-ab)/c) - opt'/(ar-ab)*r+ar/(ar-ab)


#Instance 4.3
#Intersection (t=8/9,cost=1)
#t<=4/3
#fr = 1/3+0.75t
#fb = 1
#t>=4/3
#fr=4/3
#fb=1

#VARS : q0,r,k0,k1,k2   *1*


#q0 = (r-1)/2
eq1 = [1,-1.0/2.0,0,0,0,-1.0/2.0]
#Q0 (0<t<8/9)= k0*exp(0.75t)- r+1,Q(0)=q0
eq2 = [-1,-1,1,0,0,-1]

#Q1 (8/9<t<4/3)= k1*exp(0.75t)+1,Q0(8/9)=Q1(8/9)
eq3 = [0,-1,exp(2.0/3.0),-exp(2.0/3.0),0,0]

#Q2(4/3<t) = k2+1 , Q1(4/3) = Q2(4/3)
eq4 = [0.0,0.0,0.0,exp(1.0),-1.0,0.0]
#cost(inf)=1*r = q_inf*(4/3)+(1-q_inf)*cost(4/3)=q_inf*(4/3)+(1-q_inf)*r
#two possiblites :q_inf=0=>k2=-1 (impossible) or r=4/3
eq5 = [0.0,1.0,0.0,0.0,0.0,4.0/3.0]
u=np.array([eq1,eq2,eq3,eq4,eq5])


x2=numpy.linalg.solve(u[:,:5],u[:,5])


#our soultion is q0=0.2611,r=1.522,k0=q0+r-1~=0.7831,k1=0,k2=0
eq1 = [1,-1.0/2.0,0,-1.0/2.0]
#Q0 (0<t<8/9)= k0*exp(0.75t)- r+1,Q(0)=q0
eq2 = [-1,-1,1,-1]

#Q1 (8/9<t<4/3)= k1*exp(0.75t)+1,Q0(8/9)=Q1(8/9)
eq3 = [0,-1,exp(2.0/3.0),0]

#Q2(4/3<t) = k2+1 , Q1(4/3) = Q2(4/3)
eq4 = [0.0,0.0,0.0,0.0]
#cost(inf)=1*r = q_inf*(4/3)+(1-q_inf)*cost(4/3)=q_inf*(4/3)+(1-q_inf)*r
#two possiblites :q_inf=0=>k2=-1 (impossible) or r=4/3
eq5 = [0.0,1.0,0.0,4.0/3.0]
u=np.array([eq1,eq2,eq3])


x3=numpy.linalg.solve(u[:,:3],u[:,3])

def fun(x):
    return x[0]**2-1 if x[0]<0.5 else (x[0]-1)**2-1

sol = optimize.root(fun, [1.2], method='hybr')
#print(sol.x)




